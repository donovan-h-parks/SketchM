use std::cmp::Ordering;
use std::cmp::{max, min};
use std::fs::File;
use std::io::stdout;
use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

use anyhow::Result;
use log::info;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::hashing::ItemHash;
use crate::progress::progress_bar;
use crate::sketch::{read_sketch, Sketch, SketchHeader, WeightedSketch};

/// Statistics for unweighted distance between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance<'a> {
    pub query_id: &'a str,
    pub ref_id: &'a str,

    pub ani_jaccard: f32,
    pub ani_query: f32,
    pub ani_ref: f32,

    pub containment_query: f32,
    pub containment_ref: f32,
}

/// Statistics for weighted distance between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct WeightedSketchDistance<'a> {
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

/// Calculate unweighted distance statistics between query and reference sketches.
pub fn calc_sketch_distances(
    query_sketch_files: &[String],
    ref_sketch_files: &[String],
    min_ani: f64,
    output_file: &Option<String>,
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

    // create writer
    let writer: Box<dyn Write + Send> = match output_file {
        Some(output_file) => {
            let out_path = Path::new(output_file);
            Box::new(File::create(out_path)?)
        }
        None => Box::new(stdout()),
    };

    let writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

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
            let dist_stats = unweighted_distances(query_sketch, ref_sketch, k);
            if f32::max(dist_stats.ani_query, dist_stats.ani_ref) >= min_ani as f32 {
                safe_writer
                    .lock()
                    .expect("Failed to obtain writer")
                    .serialize(dist_stats)
                    .expect("Failed to write sketch to file")
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
    output_file: &Option<String>,
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

    // create writer
    let writer: Box<dyn Write + Send> = match output_file {
        Some(output_file) => {
            let out_path = Path::new(output_file);
            Box::new(File::create(out_path)?)
        }
        None => Box::new(stdout()),
    };

    let writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

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
            let dist_stats = weighted_distances(query_sketch, ref_sketch, k);
            if f32::max(dist_stats.ani_query, dist_stats.ani_ref) >= min_ani as f32 {
                safe_writer
                    .lock()
                    .expect("Failed to obtain writer")
                    .serialize(dist_stats)
                    .expect("Failed to write sketch to file")
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

/// Calculate unweighted distance statistics between two sketches.
pub fn unweighted_distances<'a>(
    query_sketch: &'a Sketch,
    ref_sketch: &'a Sketch,
    k: u8,
) -> SketchDistance<'a> {
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
) -> SketchDistance<'a> {
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

    SketchDistance {
        query_id: query_name,
        ref_id: ref_name,
        ani_jaccard,
        ani_query,
        ani_ref,
        containment_query,
        containment_ref,
    }
}

/// Estimates weighted statistics based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn weighted_distances<'a>(
    query_sketch: &'a WeightedSketch,
    ref_sketch: &'a WeightedSketch,
    k: u8,
) -> WeightedSketchDistance<'a> {
    let mut num_common_hashes: u64 = 0;
    let mut min_hashes_weighted: u64 = 0;
    let mut max_hashes_weighted: u64 = 0;
    let mut weighted_query: u64 = 0;
    let mut weighted_ref: u64 = 0;

    // process k-mers (hashes) in common between sketches
    let mut q_iter = query_sketch.hashes.iter().peekable();
    let mut r_iter = ref_sketch.hashes.iter().peekable();
    while let (Some(q_item), Some(r_item)) = (q_iter.peek(), r_iter.peek()) {
        let (q_hash, q_count) = *q_item;
        let (r_hash, r_count) = *r_item;

        match q_hash.cmp(r_hash) {
            Ordering::Less => {
                weighted_query += *q_count as u64;
                q_iter.next();
            }
            Ordering::Greater => {
                weighted_ref += *r_count as u64;
                r_iter.next();
            }
            Ordering::Equal => {
                num_common_hashes += 1;
                min_hashes_weighted += min(*q_count as u64, *r_count as u64);
                max_hashes_weighted += max(*q_count as u64, *r_count as u64);

                weighted_query += *q_count as u64;
                weighted_ref += *r_count as u64;

                q_iter.next();
                r_iter.next();
            }
        }
    }

    // read rest of hashes since one sketch will run out of hashes first
    for (_hash, count) in q_iter {
        weighted_query += *count as u64;
        max_hashes_weighted += *count as u64;
    }

    for (_hash, count) in r_iter {
        weighted_ref += *count as u64;
        max_hashes_weighted += *count as u64;
    }

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

    if num_common_hashes > 0 {
        // calculate unweighted and weighted Jaccard index based ANI estimate
        let total_hashes =
            query_sketch.hashes.len() as u64 + ref_sketch.hashes.len() as u64 - num_common_hashes;
        let jaccard = num_common_hashes as f32 / total_hashes as f32;
        let mash_dist = mash_distance(jaccard, k);
        ani_jaccard = 100.0 * (1.0 - mash_dist);

        let jaccard_w = min_hashes_weighted as f32 / max_hashes_weighted as f32;
        let mash_dist_w = mash_distance(jaccard_w, k);
        ani_jaccard_weighted = 100.0 * (1.0 - mash_dist_w);

        // calculate containment of query relative to reference, and reference relative to query
        containment_query = num_common_hashes as f32 / query_sketch.hashes.len() as f32;
        containment_ref = num_common_hashes as f32 / ref_sketch.hashes.len() as f32;

        containment_query_weighted = min_hashes_weighted as f32 / weighted_query as f32;
        containment_ref_weighted = min_hashes_weighted as f32 / weighted_ref as f32;

        // determine ANI from containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani_query = 100.0 * containment_query.powf(exp);
        ani_ref = 100.0 * containment_ref.powf(exp);

        ani_query_weighted = 100.0 * containment_query_weighted.powf(exp);
        ani_ref_weighted = 100.0 * containment_ref_weighted.powf(exp);
    }

    WeightedSketchDistance {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_raw_distance() {
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
