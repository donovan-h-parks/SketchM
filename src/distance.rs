use std::cmp::Ordering;
use std::io::stdout;
use std::fs::File;
use std::path::Path;
use std::io::Write;

use anyhow::Result;
use serde::{Deserialize, Serialize};
use log::info;
use num_format::ToFormattedString;

use crate::hashing::ItemHash;
use crate::sketch::{read_sketch, Sketch, SketchHeader};
use crate::config::LOCALE;

/// Statistics for the distance between two sketches.
#[derive(Debug, Serialize, Deserialize)]
pub struct SketchDistance {
    pub query_id: String,
    pub ref_id: String,
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
    output_file: &Option<String>
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
            None => { prev_sketch_header = Some(sketch_header) },
        }

        query_sketches.extend(sketches);
    }
    info!(" - read {} query sketches", query_sketches.len().to_formatted_string(&LOCALE));

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
    info!(" - read {} reference sketches", ref_sketches.len().to_formatted_string(&LOCALE));

    // create writer
    let writer: Box<dyn Write> = match output_file {
        Some(output_file) => {
            let out_path = Path::new(output_file);
            Box::new(File::create(out_path)?)
        },
        None => {
            Box::new(stdout())        
        }
    };

    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(writer);
    
    // calculate statistics between query and reference sketches
    let k = prev_sketch_header.expect("Invalid sketch header").params.k();
    for query_sketch in &query_sketches {
        for ref_sketch in &ref_sketches {
            let dist_stats = distance(query_sketch, ref_sketch, k);
            if f32::max(dist_stats.ani_query, dist_stats.ani_ref) >= min_ani as f32 {
                writer.serialize(dist_stats)?;
            }
        }
    }

    Ok(())
}

/// Calculate distance statistics between two sketches.
pub fn distance(query_sketch: &Sketch, ref_sketch: &Sketch, k: u8) -> SketchDistance {
    raw_distance(
        &query_sketch.hashes,
        &ref_sketch.hashes,
        &query_sketch.name,
        &ref_sketch.name,
        k
    )
}

/// Estimates statistics based on the k-mer hash sets of two sketches.
///
/// If the hash values are not sorted the calculated distance statistics are 
/// invalid. Currently, this constraint is enforced when creating the sketch.
pub fn raw_distance(
    query_hashes: &[ItemHash],
    ref_hashes: &[ItemHash],
    query_name: &str,
    ref_name: &str,
    k: u8,
) -> SketchDistance {
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

    // calculate distance statistics
    let containment_query = num_common_hashes as f32 / query_hashes.len() as f32;
    let containment_ref = num_common_hashes as f32 / ref_hashes.len() as f32;

    // determine ANI from containment 
    // see Hera et al., 2023, Genome Research
    let exp = 1.0 / k as f32;
    let ani_query = 100.0 * containment_query.powf(exp);
    let ani_ref = 100.0 * containment_ref.powf(exp);

    SketchDistance {
        query_id: query_name.to_string(),
        ref_id: ref_name.to_string(),
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
