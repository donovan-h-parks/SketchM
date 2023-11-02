use std::cmp::Ordering;
use std::fs::File;
use std::io::stdout;
use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

use anyhow::Result;
use log::info;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::progress::progress_bar;
use crate::sketch::{read_sketch, SketchHeader, SketchType};

/// Statistics for the distance between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance {
    pub query_id: String,
    pub ref_id: String,
    pub ani_jaccard: f32,
    pub ani_query: f32,
    pub ani_ref: f32,
    pub containment_query: f32,
    pub containment_ref: f32,
}

/// Calculate distance statistics between query and reference sketches.
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

        query_sketches.extend(sketches);
    }

    // load all reference sketches into memory
    info!("Loading reference sketches into memory:");
    let mut ref_sketches = Vec::new();
    for ref_sketch_file in ref_sketch_files {
        let (sketch_header, sketches) = read_sketch(ref_sketch_file)?;
        if let Some(ref header) = prev_sketch_header {
            header.params.check_compatibility(&sketch_header.params)?;
        }

        ref_sketches.extend(sketches);
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
            let dist_stats = distance(query_sketch, ref_sketch, k);
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

/// Calculate distance statistics between two sketches.
pub fn distance(query_sketch: &SketchType, ref_sketch: &SketchType, k: u8) -> SketchDistance {
    raw_distance(
        query_sketch,
        ref_sketch,
        query_sketch.name(),
        ref_sketch.name(),
        k,
    )
}

/// Estimates statistics based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn raw_distance(
    query_hashes: &SketchType,
    ref_hashes: &SketchType,
    query_name: &str,
    ref_name: &str,
    k: u8,
) -> SketchDistance {
    let mut num_common_hashes: u64 = 0;
    let mut q_iter = query_hashes.hashes().peekable();
    let mut r_iter = ref_hashes.hashes().peekable();
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
        let total_hashes =
            query_hashes.unique_hash_count() + ref_hashes.unique_hash_count() - num_common_hashes;
        let jaccard = num_common_hashes as f32 / total_hashes as f32;
        let mash_distance = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / k as f32;
        ani_jaccard = 100.0 * (1.0 - mash_distance);

        // calculate containment of query relative to reference, and reference relative to query
        containment_query = num_common_hashes as f32 / query_hashes.unique_hash_count() as f32;
        containment_ref = num_common_hashes as f32 / ref_hashes.unique_hash_count() as f32;

        // determine ANI from containment
        // see Hera et al., 2023, Genome Research
        let exp = 1.0 / k as f32;
        ani_query = 100.0 * containment_query.powf(exp);
        ani_ref = 100.0 * containment_ref.powf(exp);
    }

    SketchDistance {
        query_id: query_name.to_string(),
        ref_id: ref_name.to_string(),
        ani_jaccard,
        ani_query,
        ani_ref,
        containment_query,
        containment_ref,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kc(hashes: &[ItemHash]) -> KmerCounter {
        let mut kc = KmerCounter::new();
        for hash in hashes {
            kc.insert(*hash, 1);
        }

        kc
    }

    #[test]
    fn test_raw_distance() {
        let d = raw_distance(&kc(&[0, 1, 2]), &kc(&[1, 2]), "q", "r", 3);
        assert_eq!(d.containment, 2. / 2.);
        assert_eq!(d.abundance, 100.0 * 2. / 3.);
        assert_eq!(d.jaccard, 2. / 3.);
        assert_eq!(d.num_query_hashes, 3);
        assert_eq!(d.num_reference_hashes, 2);
        assert_eq!(d.num_common_hashes, 2);
        assert_eq!(d.total_hashes, 3);
        assert_eq!(d.weighted_query, 3);
        assert_eq!(d.weighted_ref, 2);
        assert_eq!(d.query, "q");
        assert_eq!(d.reference, "r");

        let d = raw_distance(&kc(&[0, 2]), &kc(&[1, 2]), "q", "r", 3);
        assert_eq!(d.containment, 1. / 2.);
        assert_eq!(d.abundance, 100.0 * 1. / 2.);
        assert_eq!(d.num_query_hashes, 2);
        assert_eq!(d.num_reference_hashes, 2);
        assert_eq!(d.jaccard, 1. / 3.);
        assert_eq!(d.num_common_hashes, 1);
        assert_eq!(d.total_hashes, 3);
        assert_eq!(d.weighted_query, 2);
        assert_eq!(d.weighted_ref, 2);

        let d = raw_distance(&kc(&[0, 1]), &kc(&[2, 3]), "q", "r", 3);
        assert_eq!(d.containment, 0. / 2.);
        assert_eq!(d.abundance, 100.0 * 0. / 2.);
        assert_eq!(d.jaccard, 0. / 2.);
        assert_eq!(d.num_query_hashes, 2);
        assert_eq!(d.num_reference_hashes, 2);
        assert_eq!(d.num_common_hashes, 0);
        assert_eq!(d.total_hashes, 4);
        assert_eq!(d.weighted_query, 2);
        assert_eq!(d.weighted_ref, 2);
    }
}
