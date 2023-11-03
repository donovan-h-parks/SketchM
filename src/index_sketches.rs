use std::fs::File;
use std::io::BufWriter;

use anyhow::Result;
use log::info;
use nohash_hasher::IntMap;
use serde::{Deserialize, Serialize};

use crate::hashing::ItemHash;
use crate::progress::progress_bar;
use crate::sketch::read_sketch;

pub const INDEX_EXT: &str = ".idx";

#[derive(Default, Serialize, Deserialize)]
struct Index {
    pub genome_map: Vec<String>,
    pub kmer_index: IntMap<ItemHash, Vec<u32>>,
}

/// Create index indicating all genomes with a given k-mer.
// Genomes are stored in a vector where their index in the vector is the
// identifier used in the k-mer index. The k-mer index indicates all genomes
// that contain a given k-mer (i.e. <k-mer> -> [0, 5, 9])
pub fn index_sketches(sketch_files: &[String], output_file: &str) -> Result<()> {
    let mut index = Index::default();

    // add all sketches from all sketch files to the index
    for sketch_file in sketch_files {
        info!("Reading {sketch_file}:");
        let (sketch_header, sketches) = read_sketch(sketch_file)?;

        info!("Adding sketches to k-mer index:");
        let progress_bar = progress_bar(sketch_header.num_sketches as u64);
        for sketch in sketches {
            let genome_id = index.genome_map.len() as u32;
            index.genome_map.push(sketch.name().to_string());
            for hash in sketch.hashes() {
                let genomes = index.kmer_index.entry(hash).or_insert(Vec::new());
                genomes.push(genome_id);
            }

            progress_bar.inc(1);
        }
        progress_bar.finish();
    }

    // write out index
    info!("Writing index to file.");
    let mut out_file = output_file.to_string();
    if !out_file.ends_with(INDEX_EXT) {
        out_file += INDEX_EXT;
    };

    let fp = File::create(out_file)?;
    let writer = BufWriter::new(fp);
    bincode::serialize_into(writer, &index)?;

    Ok(())
}
