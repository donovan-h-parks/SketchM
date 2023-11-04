use std::fs::File;
use std::io::BufWriter;

use ahash::HashMap;
use anyhow::{Context, Result};
use log::info;
use num_format::ToFormattedString;
use serde::{Deserialize, Serialize};
use tinyvec::TinyVec;

use crate::config::LOCALE;
use crate::hashing::ItemHash;
use crate::maybe_gzip_io::maybe_gzip_reader;
use crate::progress::progress_bar;
use crate::sketch::{SketchHeader, SketchType};

pub const INDEX_EXT: &str = ".idx";

/// Structure for determining all genomes with a given k-mer.
// Genomes are stored as a vector where their index in the vector is the
// identifier used in the k-mer index. The name and number of hashes are
// stored for each genome as this information is needed for reporting and
// calculating distance statistics. The k-mer map indicates all genomes
// that contain a given k-mer as specified by integers indicating the index
// of the genome in the genomes vector (i.e. <k-mer> -> [0, 5, 9]).
#[derive(Default, Serialize, Deserialize)]
pub struct Index {
    pub genomes: Vec<(String, u32)>,

    // There is no reason to re-hash the k-mers here, but this is
    // currently necessary to avoid quadratic time when deseralizing
    // this HashMap. This is a known issue impacting any non-random hasher:
    // https://morestina.net/blog/1843/the-stable-hashmap-trap. TinyVec
    // can decreases computation time by ~25% and peak memory usage by ~15%
    // compared to Vec, though this depends on the reference set. A size of 3 is
    // used since ~70% of k-mers appear once, ~15% appear twice, and ~5% appear
    // 3 times in GTDB R214 when k=31. There is no optimal setting here though since
    // SketchM is general use and the input genomes could be all closely related strains
    // (e.g. dereplicating E. coli strains) or all highly diverse (e.g. GTDB R214
    // species representatives).
    pub kmer_index: HashMap<ItemHash, TinyVec<[u32; 3]>>,
}

/// Create index indicating all genomes with a given k-mer.
pub fn index_sketches(sketch_files: &[String], output_file: &str) -> Result<()> {
    let mut index = Index::default();

    // add all sketches from all sketch files to the index
    let mut total_kmers: u64 = 0;
    for sketch_file in sketch_files {
        info!("Adding sketches to k-mer index:");

        let mut reader = maybe_gzip_reader(sketch_file)
            .context(format!("Unable to read sketch file: {}", sketch_file))?;

        let sketch_header: SketchHeader = bincode::deserialize_from(&mut reader)?;
        let progress_bar = progress_bar(sketch_header.num_sketches as u64);

        for _ in 0..sketch_header.num_sketches {
            let sketch: SketchType = bincode::deserialize_from(&mut reader)?;

            let genome_id = index.genomes.len() as u32;

            index
                .genomes
                .push((sketch.name().to_string(), sketch.unique_hash_count() as u32));

            for hash in sketch.hashes() {
                let genomes = index.kmer_index.entry(hash).or_default();
                genomes.push(genome_id);
                total_kmers += 1;
            }

            progress_bar.inc(1);
        }
        progress_bar.finish();
    }

    info!(
        "Index contains {} genomes with {} unique and {} total k-mers.",
        index.genomes.len().to_formatted_string(&LOCALE),
        index.kmer_index.len().to_formatted_string(&LOCALE),
        total_kmers.to_formatted_string(&LOCALE)
    );

    // Histogram
    let mut counts: Vec<u64> = vec![0; 10];
    let mut other = 0;
    for gids in index.kmer_index.values() {
        if gids.len() <= 10 {
            counts[gids.len() - 1] += 1;
        } else {
            other += 1;
        }
    }

    info!("Histogram (No. genomes: occurences)");
    for (idx, count) in counts.into_iter().enumerate() {
        info!(
            " {}: {} ({:.2}%)",
            idx + 1,
            count.to_formatted_string(&LOCALE),
            100.0 * count as f32 / index.kmer_index.len() as f32
        );
    }
    info!(
        ">10: {} ({:.2}%)",
        other.to_formatted_string(&LOCALE),
        100.0 * other as f32 / index.kmer_index.len() as f32
    );

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

/// Read index.
pub fn read_index(index_file: &str) -> Result<Index> {
    let mut reader = maybe_gzip_reader(index_file)
        .context(format!("Unable to read index file: {}", index_file))?;

    let index: Index = bincode::deserialize_from(&mut reader)?;

    Ok(index)
}
