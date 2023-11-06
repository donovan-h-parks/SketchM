use std::cmp::{max, min, Ordering};
use std::iter;
use std::sync::atomic::{self, AtomicU64};
use std::sync::Mutex;

use anyhow::{bail, Context, Result};
use log::info;
use num_format::ToFormattedString;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::prelude::ParallelBridge;
use serde::{Deserialize, Serialize};

use crate::config::LOCALE;
use crate::hashing::ItemHash;
use crate::index_sketches::{read_index, Index};
use crate::maybe_gzip_io::{maybe_gzip_csv_writer, maybe_gzip_reader};
use crate::progress::progress_bar;
use crate::sketch::{read_sketch, Sketch, SketchHeader, SketchType, WeightedSketch};

/// ANI and AF between two unweighted or weighted sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct AniAf<'a> {
    pub query_id: &'a str,
    pub ref_id: &'a str,
    pub ani: f32,
    pub af: f32,
}

/// Unweighted distance statistics between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct UnweightedDistances<'a> {
    pub query_id: &'a str,
    pub ref_id: &'a str,

    pub ani_jaccard: f32,
    pub ani_query: f32,
    pub ani_ref: f32,

    pub containment_query: f32,
    pub containment_ref: f32,
}

/// Weighted distance statistics between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct WeightedDistances<'a> {
    pub query_id: &'a str,
    pub ref_id: &'a str,

    pub ani_jaccard: f32,
    pub ani_query: f32,
    pub ani_ref: f32,

    pub ani_jaccard_weighted: f32,
    pub ani_query_weighted: f32,
    pub ani_ref_weighted: f32,

    pub containment_query: f32,
    pub containment_ref: f32,
    pub containment_query_weighted: f32,
    pub containment_ref_weighted: f32,
}

/// Calculate unweighted distances between query sketches and a reference index.
pub fn calc_sketch_distances_to_index(
    query_sketch_files: &[String],
    reference_index: &str,
    min_ani: f64,
    additional_stats: bool,
    single_genome_set: bool,
    output_file: &Option<String>,
    threads: usize,
) -> Result<()> {
    // load reference index
    info!("Loading reference index into memory:");
    let (index_header, index) = read_index(reference_index)?;
    let num_ref_genomes = index.genomes.len() as u32;
    info!(
        " - index contains {} genomes with {} unique k-mers",
        num_ref_genomes.to_formatted_string(&LOCALE),
        index.kmer_index.len().to_formatted_string(&LOCALE)
    );

    // Need to compress with multiple threads or this becomes the rate limiting
    // step. Setting this to the maximum number of threads specified gives good
    // performance even if these theads need to compete with processing of sketches.
    // Perhaps each thread should write to a seperate file and then these should be
    // merge at the end. Unclear if this would help though if the ultimate issue is
    // saturating the systems IO.
    let writer = maybe_gzip_csv_writer(output_file, threads)?;

    // calculate statistics between query and reference index
    info!("Calculating ANI between query sketches and reference index:");
    let safe_writer = Mutex::new(writer);
    let num_pairs = AtomicU64::new(0);
    for query_sketch_file in query_sketch_files {
        let mut reader = maybe_gzip_reader(query_sketch_file)
            .context(format!("Unable to read sketch file: {}", query_sketch_file))?;

        let sketch_header: SketchHeader = bincode::deserialize_from(&mut reader)?;
        let k = sketch_header.params.k();

        index_header
            .params
            .check_compatibility(&sketch_header.params)?;

        let progress_bar = progress_bar(sketch_header.num_sketches as u64);
        iter::from_fn(move || bincode::deserialize_from(&mut reader).ok())
            .enumerate()
            .par_bridge()
            .try_for_each(|(idx, sketch): (usize, SketchType)| {
                let sketch = sketch.to_sketch();

                let mut genome_stride_idx = num_ref_genomes;
                if single_genome_set {
                    // optimize calculations by reporting only lower triangle
                    genome_stride_idx = idx as u32;

                    // sanity check the query and reference sets are identical
                    if sketch.name != index.genomes[idx].0 {
                        bail!("Genome sets are not identical as required by `single-genome-set`.");
                    }
                }

                if additional_stats {
                    let dist_stats =
                        unweighted_dists_to_index(&sketch, &index, k, min_ani, genome_stride_idx);

                    num_pairs.fetch_add(dist_stats.len() as u64, atomic::Ordering::Acquire);

                    let mut locked_writer = safe_writer.lock().expect("Failed to obtain writer");
                    for ds in dist_stats {
                        locked_writer.serialize(ds)?;
                    }
                } else {
                    let ani_afs =
                        unweighted_ani_to_index(&sketch, &index, k, min_ani, genome_stride_idx);

                    num_pairs.fetch_add(ani_afs.len() as u64, atomic::Ordering::Acquire);

                    let mut locked_writer = safe_writer.lock().expect("Failed to obtain writer");
                    for ani_af in ani_afs {
                        locked_writer.serialize(ani_af)?;
                    }
                }

                progress_bar.inc(1);

                Ok(())
            })?;

        progress_bar.finish();
    }

    info!(
        "Found {} genome pairs meeting minimum ANI criterion.",
        num_pairs
            .load(atomic::Ordering::Relaxed)
            .to_formatted_string(&LOCALE)
    );

    Ok(())
}

/// Calculate common hashes between query sketch and reference index.
// In general, only the lower triangle of results is required when doing pairwise comparisons.
// This is archieved by incrementing the genome_stride_idx. If all pairwise results are desired,
// or when the query and reference sets are independent, the genome_stride_idx should be set to
// index.genomes.len().
fn common_hashes_to_index(sketch: &Sketch, index: &Index, genome_stride_idx: u32) -> Vec<u32> {
    let mut common_hash_counts: Vec<u32> = vec![0; genome_stride_idx as usize];
    for hash in &sketch.hashes {
        if let Some(genome_idx) = index.kmer_index.get(hash) {
            for idx in genome_idx
                .iter()
                .copied()
                .filter(|&idx| idx < genome_stride_idx)
            {
                common_hash_counts[idx as usize] += 1;
            }
        }
    }

    common_hash_counts
}

/// Calculate unweighted distance statistics between sketch and reference index.
fn unweighted_dists_to_index<'a>(
    sketch: &'a Sketch,
    index: &'a Index,
    k: u8,
    min_ani: f64,
    genome_stride_idx: u32,
) -> Vec<UnweightedDistances<'a>> {
    let common_hash_counts = common_hashes_to_index(sketch, index, genome_stride_idx);

    let mut dists = Vec::new();
    for (genome_idx, num_common_hashes) in common_hash_counts.into_iter().enumerate() {
        if num_common_hashes == 0 {
            continue;
        }

        let (gid, num_ref_hashes) = &index.genomes[genome_idx];
        let num_ref_hashes = *num_ref_hashes;
        let num_query_hashes = sketch.hash_count() as u32;

        // calculate Jaccard index based ANI estimate
        let total_hashes = num_query_hashes + num_ref_hashes - num_common_hashes;
        let jaccard = num_common_hashes as f32 / total_hashes as f32;
        let mash_dist = mash_distance(jaccard, k);
        let ani_jaccard = 100.0 * (1.0 - mash_dist);

        // calculate containment of query relative to reference, and reference relative to query
        let containment_query = num_common_hashes as f32 / num_query_hashes as f32;
        let containment_ref = num_common_hashes as f32 / num_ref_hashes as f32;

        // determine ANI from containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        let ani_query = 100.0 * containment_query.powf(exp);
        let ani_ref = 100.0 * containment_ref.powf(exp);

        if f32::max(ani_query, ani_ref) >= min_ani as f32 {
            let d = UnweightedDistances {
                query_id: &sketch.name,
                ref_id: gid,
                ani_jaccard,
                ani_query,
                ani_ref,
                containment_query,
                containment_ref,
            };

            dists.push(d);
        }
    }

    dists
}

/// Calculate unweighted ANI and AF between sketch and reference index.
fn unweighted_ani_to_index<'a>(
    sketch: &'a Sketch,
    index: &'a Index,
    k: u8,
    min_ani: f64,
    genome_stride_idx: u32,
) -> Vec<AniAf<'a>> {
    let common_hash_counts = common_hashes_to_index(sketch, index, genome_stride_idx);

    let mut dists = Vec::new();
    for (genome_idx, num_common_hashes) in common_hash_counts.into_iter().enumerate() {
        if num_common_hashes == 0 {
            continue;
        }

        let (gid, num_ref_hashes) = &index.genomes[genome_idx];
        let num_ref_hashes = *num_ref_hashes;
        let num_query_hashes = sketch.hash_count() as u32;

        // determine maximum containment between query and reference
        let max_containment = if num_query_hashes < num_ref_hashes {
            num_common_hashes as f32 / num_query_hashes as f32
        } else {
            num_common_hashes as f32 / num_ref_hashes as f32
        };

        // determine ANI for maximum containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        let ani = 100.0 * max_containment.powf(exp);

        if ani >= min_ani as f32 {
            let d = AniAf {
                query_id: &sketch.name,
                ref_id: gid,
                ani,
                af: max_containment,
            };

            dists.push(d);
        }
    }

    dists
}

/// Calculate unweighted distance statistics between query and reference sketches.
pub fn calc_sketch_distances(
    query_sketch_files: &[String],
    ref_sketch_files: &[String],
    min_ani: f64,
    additional_stats: bool,
    output_file: &Option<String>,
    threads: usize,
) -> Result<()> {
    // load all query sketches into memory
    info!("Loading query sketches into memory:");
    let mut prev_sketch_header: Option<SketchHeader> = None;
    let mut query_sketches = Vec::new();
    for query_sketch_file in query_sketch_files {
        let (sketch_header, sketches) = read_sketch(query_sketch_file)?;
        match prev_sketch_header {
            Some(ref header) => {
                header.params.check_compatibility(&sketch_header.params)?;
            }
            None => prev_sketch_header = Some(sketch_header),
        }

        for sketch in sketches {
            query_sketches.push(sketch.to_sketch());
        }
    }

    // load all reference sketches into memory
    info!("Loading reference sketches into memory:");
    let mut ref_sketches = Vec::new();
    for ref_sketch_file in ref_sketch_files {
        let (sketch_header, sketches) = read_sketch(ref_sketch_file)?;
        if let Some(ref header) = prev_sketch_header {
            header.params.check_compatibility(&sketch_header.params)?;
        }

        for sketch in sketches {
            ref_sketches.push(sketch.to_sketch());
        }
    }

    // create CSV writer using multiple threads to compressing output
    // if file name has a `gz` extension
    let writer = maybe_gzip_csv_writer(output_file, threads)?;

    // calculate statistics between query and reference sketches
    info!("Calculating ANI between query and reference sketches:");
    let progress_bar = progress_bar(query_sketches.len() as u64 * ref_sketches.len() as u64);
    let safe_writer = Mutex::new(writer);

    let k = prev_sketch_header
        .expect("Invalid sketch header")
        .params
        .k();
    for query_sketch in &query_sketches {
        ref_sketches.par_iter().for_each(|ref_sketch| {
            if additional_stats {
                let dist_stats = unweighted_distances(query_sketch, ref_sketch, k);
                if f32::max(dist_stats.ani_query, dist_stats.ani_ref) >= min_ani as f32 {
                    safe_writer
                        .lock()
                        .expect("Failed to obtain writer")
                        .serialize(dist_stats)
                        .expect("Failed to write distance stats to file");
                }
            } else {
                let ani_af = unweighted_ani(query_sketch, ref_sketch, k);
                if ani_af.ani >= min_ani as f32 {
                    safe_writer
                        .lock()
                        .expect("Failed to obtain writer")
                        .serialize(ani_af)
                        .expect("Failed to write ANI/AF to file");
                }
            }
            progress_bar.inc(1);
        });
    }

    progress_bar.finish();

    Ok(())
}

/// Calculate weighted distance statistics between query and reference sketches.
// This function is nearly identical to calc_sketch_distances() which is a poor
// replication of code, but the explicit use of Sketch and WeightedSketch (i.e.
// the inner types of SketchType) increase performance by ~20%.
pub fn calc_weighted_sketch_distances(
    query_sketch_files: &[String],
    ref_sketch_files: &[String],
    min_ani: f64,
    additional_stats: bool,
    output_file: &Option<String>,
    threads: usize,
) -> Result<()> {
    // load all query sketches into memory
    info!("Loading query sketches into memory:");
    let mut prev_sketch_header: Option<SketchHeader> = None;
    let mut query_sketches = Vec::new();
    for query_sketch_file in query_sketch_files {
        let (sketch_header, sketches) = read_sketch(query_sketch_file)?;
        match prev_sketch_header {
            Some(ref header) => {
                header.params.check_compatibility(&sketch_header.params)?;
            }
            None => prev_sketch_header = Some(sketch_header),
        }

        for sketch in sketches {
            query_sketches.push(sketch.to_weighted_sketch());
        }
    }

    // load all reference sketches into memory
    info!("Loading reference sketches into memory:");
    let mut ref_sketches = Vec::new();
    for ref_sketch_file in ref_sketch_files {
        let (sketch_header, sketches) = read_sketch(ref_sketch_file)?;
        if let Some(ref header) = prev_sketch_header {
            header.params.check_compatibility(&sketch_header.params)?;
        }

        for sketch in sketches {
            ref_sketches.push(sketch.to_weighted_sketch());
        }
    }

    // create CSV writer using multiple threads to compressing output
    // if file name has a `gz` extension
    let writer = maybe_gzip_csv_writer(output_file, threads)?;

    // calculate statistics between query and reference sketches
    info!("Calculating ANI between query and reference sketches:");
    let progress_bar = progress_bar(query_sketches.len() as u64 * ref_sketches.len() as u64);
    let safe_writer = Mutex::new(writer);

    let k = prev_sketch_header
        .expect("Invalid sketch header")
        .params
        .k();
    for query_sketch in &query_sketches {
        ref_sketches.par_iter().for_each(|ref_sketch| {
            if additional_stats {
                let dist_stats = weighted_distances(query_sketch, ref_sketch, k);
                if f32::max(dist_stats.ani_query, dist_stats.ani_ref) >= min_ani as f32 {
                    safe_writer
                        .lock()
                        .expect("Failed to obtain writer")
                        .serialize(dist_stats)
                        .expect("Failed to write distance stats to file")
                }
            } else {
                let ani_af = weighted_ani(query_sketch, ref_sketch, k);
                if ani_af.ani >= min_ani as f32 {
                    safe_writer
                        .lock()
                        .expect("Failed to obtain writer")
                        .serialize(ani_af)
                        .expect("Failed to write ANI/AF to file")
                }
            }
            progress_bar.inc(1);
        });
    }

    progress_bar.finish();

    Ok(())
}

/// Calculate the Mash distance from the Jaccard index.
// see: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x
fn mash_distance(jaccard: f32, k: u8) -> f32 {
    -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / k as f32
}

/// Determine number of common hashes between two sketches.
fn common_hashes(query_hashes: &[ItemHash], ref_hashes: &[ItemHash]) -> u64 {
    let mut num_common_hashes: u64 = 0;
    let mut q_iter = query_hashes.iter().peekable();
    let mut r_iter = ref_hashes.iter().peekable();
    while let (Some(q_hash), Some(r_hash)) = (q_iter.peek(), r_iter.peek()) {
        match q_hash.cmp(r_hash) {
            Ordering::Less => {
                q_iter.next();
            }
            Ordering::Greater => {
                r_iter.next();
            }
            Ordering::Equal => {
                num_common_hashes += 1;
                q_iter.next();
                r_iter.next();
            }
        }
    }

    num_common_hashes
}

/// Calculate unweighted distance statistics between two sketches.
pub fn unweighted_distances<'a>(
    query_sketch: &'a Sketch,
    ref_sketch: &'a Sketch,
    k: u8,
) -> UnweightedDistances<'a> {
    raw_unweighted_distance(
        &query_sketch.hashes,
        &ref_sketch.hashes,
        &query_sketch.name,
        &ref_sketch.name,
        k,
    )
}

/// Estimates statistics based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn raw_unweighted_distance<'a>(
    query_hashes: &[ItemHash],
    ref_hashes: &[ItemHash],
    query_name: &'a str,
    ref_name: &'a str,
    k: u8,
) -> UnweightedDistances<'a> {
    let num_common_hashes = common_hashes(query_hashes, ref_hashes);

    let mut ani_jaccard = 0.0f32;
    let mut containment_query = 0.0f32;
    let mut containment_ref = 0.0f32;
    let mut ani_query = 0.0f32;
    let mut ani_ref = 0.0f32;

    if num_common_hashes > 0 {
        // calculate Jaccard index based ANI estimate
        let total_hashes = query_hashes.len() as u64 + ref_hashes.len() as u64 - num_common_hashes;
        let jaccard = num_common_hashes as f32 / total_hashes as f32;
        let mash_dist = mash_distance(jaccard, k);
        ani_jaccard = 100.0 * (1.0 - mash_dist);

        // calculate containment of query relative to reference, and reference relative to query
        containment_query = num_common_hashes as f32 / query_hashes.len() as f32;
        containment_ref = num_common_hashes as f32 / ref_hashes.len() as f32;

        // determine ANI from containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani_query = 100.0 * containment_query.powf(exp);
        ani_ref = 100.0 * containment_ref.powf(exp);
    }

    UnweightedDistances {
        query_id: query_name,
        ref_id: ref_name,
        ani_jaccard,
        ani_query,
        ani_ref,
        containment_query,
        containment_ref,
    }
}

/// Calculate unweighted ANI and AF between two sketches.
pub fn unweighted_ani<'a>(query_sketch: &'a Sketch, ref_sketch: &'a Sketch, k: u8) -> AniAf<'a> {
    raw_unweighted_ani(
        &query_sketch.hashes,
        &ref_sketch.hashes,
        &query_sketch.name,
        &ref_sketch.name,
        k,
    )
}

/// Estimates ANI and AF based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn raw_unweighted_ani<'a>(
    query_hashes: &[ItemHash],
    ref_hashes: &[ItemHash],
    query_name: &'a str,
    ref_name: &'a str,
    k: u8,
) -> AniAf<'a> {
    let num_common_hashes = common_hashes(query_hashes, ref_hashes);

    let mut ani = 0.0f32;
    let mut max_containment = 0.0f32;

    if num_common_hashes > 0 {
        let num_query_hashes = query_hashes.len() as u64;
        let num_ref_hashes = ref_hashes.len() as u64;

        // determine maximum containment between query and reference
        max_containment = if num_query_hashes < num_ref_hashes {
            num_common_hashes as f32 / num_query_hashes as f32
        } else {
            num_common_hashes as f32 / num_ref_hashes as f32
        };

        // determine ANI for maximum containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani = 100.0 * max_containment.powf(exp);
    }

    AniAf {
        query_id: query_name,
        ref_id: ref_name,
        ani,
        af: max_containment,
    }
}

/// Weighted k-mer counts between sketches.
#[derive(Default)]
struct WeightedKmerCount {
    pub num_common_hashes: u64,
    pub min_hashes_weighted: u64,
    pub max_hashes_weighted: u64,
    pub weighted_query: u64,
    pub weighted_ref: u64,
}

/// Calculate weighted k-mer count statistics between sketches.
fn weighted_kmer_counts(
    query_sketch: &WeightedSketch,
    ref_sketch: &WeightedSketch,
) -> WeightedKmerCount {
    let mut wkc = WeightedKmerCount::default();

    // process k-mers (hashes) in common between sketches
    let mut q_iter = query_sketch.hashes.iter().peekable();
    let mut r_iter = ref_sketch.hashes.iter().peekable();
    while let (Some(q_item), Some(r_item)) = (q_iter.peek(), r_iter.peek()) {
        let (q_hash, q_count) = *q_item;
        let (r_hash, r_count) = *r_item;

        match q_hash.cmp(r_hash) {
            Ordering::Less => {
                wkc.weighted_query += *q_count as u64;
                q_iter.next();
            }
            Ordering::Greater => {
                wkc.weighted_ref += *r_count as u64;
                r_iter.next();
            }
            Ordering::Equal => {
                wkc.num_common_hashes += 1;
                wkc.min_hashes_weighted += min(*q_count as u64, *r_count as u64);
                wkc.max_hashes_weighted += max(*q_count as u64, *r_count as u64);

                wkc.weighted_query += *q_count as u64;
                wkc.weighted_ref += *r_count as u64;

                q_iter.next();
                r_iter.next();
            }
        }
    }

    // read rest of hashes since one sketch will run out of hashes first
    for (_hash, count) in q_iter {
        wkc.weighted_query += *count as u64;
        wkc.max_hashes_weighted += *count as u64;
    }

    for (_hash, count) in r_iter {
        wkc.weighted_ref += *count as u64;
        wkc.max_hashes_weighted += *count as u64;
    }

    wkc
}

/// Estimates weighted statistics based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn weighted_distances<'a>(
    query_sketch: &'a WeightedSketch,
    ref_sketch: &'a WeightedSketch,
    k: u8,
) -> WeightedDistances<'a> {
    let wkc = weighted_kmer_counts(query_sketch, ref_sketch);

    // calculate distance statistics
    let mut ani_jaccard = 0.0f32;
    let mut ani_query = 0.0f32;
    let mut ani_ref = 0.0f32;

    let mut ani_jaccard_weighted = 0.0f32;
    let mut ani_query_weighted = 0.0f32;
    let mut ani_ref_weighted = 0.0f32;

    let mut containment_query = 0.0f32;
    let mut containment_ref = 0.0f32;
    let mut containment_query_weighted = 0.0f32;
    let mut containment_ref_weighted = 0.0f32;

    if wkc.num_common_hashes > 0 {
        // calculate unweighted and weighted Jaccard index based ANI estimate
        let total_hashes = query_sketch.hashes.len() as u64 + ref_sketch.hashes.len() as u64
            - wkc.num_common_hashes;
        let jaccard = wkc.num_common_hashes as f32 / total_hashes as f32;
        let mash_dist = mash_distance(jaccard, k);
        ani_jaccard = 100.0 * (1.0 - mash_dist);

        let jaccard_w = wkc.min_hashes_weighted as f32 / wkc.max_hashes_weighted as f32;
        let mash_dist_w = mash_distance(jaccard_w, k);
        ani_jaccard_weighted = 100.0 * (1.0 - mash_dist_w);

        // calculate containment of query relative to reference, and reference relative to query
        containment_query = wkc.num_common_hashes as f32 / query_sketch.hashes.len() as f32;
        containment_ref = wkc.num_common_hashes as f32 / ref_sketch.hashes.len() as f32;

        containment_query_weighted = wkc.min_hashes_weighted as f32 / wkc.weighted_query as f32;
        containment_ref_weighted = wkc.min_hashes_weighted as f32 / wkc.weighted_ref as f32;

        // determine ANI from containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani_query = 100.0 * containment_query.powf(exp);
        ani_ref = 100.0 * containment_ref.powf(exp);

        ani_query_weighted = 100.0 * containment_query_weighted.powf(exp);
        ani_ref_weighted = 100.0 * containment_ref_weighted.powf(exp);
    }

    WeightedDistances {
        query_id: &query_sketch.name,
        ref_id: &ref_sketch.name,
        ani_jaccard,
        ani_query,
        ani_ref,
        ani_jaccard_weighted,
        ani_query_weighted,
        ani_ref_weighted,
        containment_query,
        containment_ref,
        containment_query_weighted,
        containment_ref_weighted,
    }
}

/// Estimates weighted ANI and AF based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn weighted_ani<'a>(
    query_sketch: &'a WeightedSketch,
    ref_sketch: &'a WeightedSketch,
    k: u8,
) -> AniAf<'a> {
    let wkc = weighted_kmer_counts(query_sketch, ref_sketch);

    let mut ani = 0.0f32;
    let mut max_containment = 0.0f32;

    if wkc.num_common_hashes > 0 {
        // determine maximum containment between query and reference
        max_containment = if wkc.weighted_query < wkc.weighted_ref {
            wkc.min_hashes_weighted as f32 / wkc.weighted_query as f32
        } else {
            wkc.min_hashes_weighted as f32 / wkc.weighted_ref as f32
        };

        // determine ANI for maximum containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani = 100.0 * max_containment.powf(exp);
    }

    AniAf {
        query_id: &query_sketch.name,
        ref_id: &ref_sketch.name,
        ani,
        af: max_containment,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_raw_unweighted_distance() {
        let d = raw_unweighted_distance(&[0, 1], &[2, 3], "q", "r", 3);
        assert_eq!(d.query_id, "q");
        assert_eq!(d.ref_id, "r");
        assert_eq!(d.ani_jaccard, 0.0);
        assert_eq!(d.ani_query, 0.0);
        assert_eq!(d.ani_ref, 0.0);
        assert_eq!(d.containment_query, 0.0);
        assert_eq!(d.containment_ref, 0.0);

        let d = raw_unweighted_distance(&[0, 1, 2], &[0, 1, 2], "q", "r", 3);
        assert_eq!(d.query_id, "q");
        assert_eq!(d.ref_id, "r");
        assert_eq!(d.ani_jaccard, 100.0);
        assert_eq!(d.ani_query, 100.0);
        assert_eq!(d.ani_ref, 100.0);
        assert_eq!(d.containment_query, 1.0);
        assert_eq!(d.containment_ref, 1.0);

        let d = raw_unweighted_distance(&[0, 1, 2], &[1, 2], "q", "r", 3);
        assert_eq!(d.query_id, "q");
        assert_eq!(d.ref_id, "r");
        assert_eq!(d.ani_jaccard, 100.0 * (1.0 - mash_distance(2. / 3., 3)));
        assert_eq!(d.ani_query, 100.0 * (2. / 3f32).powf(1. / 3f32));
        assert_eq!(d.ani_ref, 100.0 * (2. / 2f32).powf(1. / 3f32));
        assert_eq!(d.containment_query, 2. / 3.);
        assert_eq!(d.containment_ref, 2. / 2.);

        let d = raw_unweighted_distance(&[0, 2], &[1, 2], "q", "r", 3);
        assert_eq!(d.query_id, "q");
        assert_eq!(d.ref_id, "r");
        assert_eq!(d.ani_jaccard, 100.0 * (1.0 - mash_distance(1. / 3., 3)));
        assert_eq!(d.ani_query, 100.0 * (1. / 2f32).powf(1. / 3f32));
        assert_eq!(d.ani_ref, 100.0 * (1. / 2f32).powf(1. / 3f32));
        assert_eq!(d.containment_query, 1. / 2.);
        assert_eq!(d.containment_ref, 1. / 2.);
    }
}
