use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use anyhow::{Context, Result};
use bincode;
use itertools::Either;
use log::info;
use needletail::parse_fastx_reader;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::hashing::ItemHash;
use crate::maybe_gzip_io::maybe_gzip_reader;
use crate::progress::progress_bar;
use crate::sketch_params::SketchParams;

pub const SKETCH_VERSION: &str = "1";
pub const SKETCH_EXT: &str = ".sk";

#[derive(Debug)]
pub struct SeqFile {
    pub id: String,
    pub file: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SketchHeader {
    pub program: String,
    pub version: String,
    pub date_created: String,
    pub sketch_version: String,
    pub num_sketches: u32,
    pub params: SketchParams,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sketch {
    pub name: String,
    pub bp_count: u64,
    pub kmer_total_count: u64,
    pub hashes: Vec<ItemHash>,
}

impl Sketch {
    pub fn hash_count(&self) -> u64 {
        self.hashes.len() as u64
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WeightedSketch {
    pub name: String,
    pub bp_count: u64,
    pub kmer_total_count: u64,
    pub hashes: BTreeMap<ItemHash, u32>,
}

impl WeightedSketch {
    pub fn unique_hash_count(&self) -> u64 {
        self.hashes.len() as u64
    }

    pub fn weighted_hash_count(&self) -> u64 {
        self.hashes.values().map(|v| *v as u64).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }
}

#[derive(Serialize, Deserialize)]
pub enum SketchType {
    Unweighted(Sketch),
    Weighted(WeightedSketch),
}

impl SketchType {
    pub fn name(&self) -> &str {
        match self {
            SketchType::Unweighted(sketch) => &sketch.name,
            SketchType::Weighted(sketch) => &sketch.name,
        }
    }

    pub fn unique_hash_count(&self) -> u64 {
        match self {
            SketchType::Unweighted(sketch) => sketch.hash_count(),
            SketchType::Weighted(sketch) => sketch.unique_hash_count(),
        }
    }

    pub fn weighted_hash_count(&self) -> u64 {
        match self {
            SketchType::Unweighted(sketch) => sketch.hash_count(),
            SketchType::Weighted(sketch) => sketch.weighted_hash_count(),
        }
    }

    pub fn bp_count(&self) -> u64 {
        match self {
            SketchType::Unweighted(sketch) => sketch.bp_count,
            SketchType::Weighted(sketch) => sketch.bp_count,
        }
    }

    pub fn kmer_total_count(&self) -> u64 {
        match self {
            SketchType::Unweighted(sketch) => sketch.kmer_total_count,
            SketchType::Weighted(sketch) => sketch.kmer_total_count,
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            SketchType::Unweighted(sketch) => sketch.is_empty(),
            SketchType::Weighted(sketch) => sketch.is_empty(),
        }
    }

    /// An iterator through hash values and their associated weights.
    // see: https://stackoverflow.com/questions/77405985
    // Note: this is ~20% slower than calculating distance using
    //       the underlying concrete type.
    pub fn hashes(&self) -> impl Iterator<Item = u64> + '_ {
        match *self {
            SketchType::Unweighted(ref sketch) => Either::Left(sketch.hashes.iter().copied()),
            SketchType::Weighted(ref sketch) => Either::Right(sketch.hashes.keys().copied()),
        }
    }

    pub fn to_sketch(self) -> Sketch {
        match self {
            SketchType::Unweighted(sketch) => sketch,
            SketchType::Weighted(_) => panic!("Not an unweighted sketch"),
        }
    }

    pub fn to_weighted_sketch(self) -> WeightedSketch {
        match self {
            SketchType::Unweighted(_) => panic!("Not a weighted sketch"),
            SketchType::Weighted(sketch) => sketch,
        }
    }
}

/// Write sketch header information.
fn write_sketch_header<W>(
    writer: &mut W,
    num_sketches: u32,
    sketch_params: &SketchParams,
) -> Result<()>
where
    W: Write,
{
    let sketch_header = SketchHeader {
        program: env!("CARGO_PKG_NAME").to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        date_created: chrono::offset::Utc::now().to_string(),
        sketch_version: SKETCH_VERSION.to_string(),
        num_sketches,
        params: sketch_params.clone(),
    };

    let sketch_header_encoded = bincode::serialize(&sketch_header)?;
    writer.write_all(&sketch_header_encoded)?;

    Ok(())
}

pub fn sketch(
    genome_files: &Vec<SeqFile>,
    sketch_params: &SketchParams,
    output_file: &str,
) -> Result<()> {
    let fp = File::create(output_file)?;
    let mut writer = BufWriter::new(fp);

    // write out sketch header information
    write_sketch_header(&mut writer, genome_files.len() as u32, sketch_params)?;

    info!(
        "Calculating sketches with {} CPUs:",
        rayon::current_num_threads()
    );

    let progress_bar = progress_bar(genome_files.len() as u64);
    let safe_writer = Mutex::new(writer);

    genome_files.par_iter().for_each(|genome_file| {
        let sketch = sketch_file(genome_file, sketch_params).expect("Failed to create sketch");
        let sketch_encoded = bincode::serialize(&sketch).expect("Failed to convert sketch to JSON");
        safe_writer
            .lock()
            .expect("Failed to obtain writer")
            .write_all(&sketch_encoded)
            .expect("Failed to write sketch to file");

        progress_bar.inc(1);
    });

    progress_bar.finish();

    safe_writer
        .lock()
        .expect("Failed to obtain writer")
        .flush()?;

    Ok(())
}

/// Create sketch from sequence file.
pub fn sketch_file(seq_file: &SeqFile, sketch_params: &SketchParams) -> Result<SketchType> {
    let mut sketcher = sketch_params.create_sketcher();
    let reader = File::open(Path::new(&seq_file.file))
        .context(format!("Failed to open {}", seq_file.file))?;

    let mut fastx_reader = parse_fastx_reader(reader)?;
    while let Some(rec) = fastx_reader.next() {
        let record = rec?;
        sketcher.process_seq(&record);
    }

    if sketch_params.weighted() {
        return Ok(SketchType::Weighted(
            sketcher.to_weighted_sketch(&seq_file.id),
        ));
    }

    Ok(SketchType::Unweighted(sketcher.to_sketch(&seq_file.id)))
}

/// Read sketches.
pub fn read_sketch(sketch_file: &str) -> Result<(SketchHeader, Vec<SketchType>)> {
    let mut reader = maybe_gzip_reader(sketch_file)
        .context(format!("Unable to read sketch file: {}", sketch_file))?;

    let sketch_header: SketchHeader = bincode::deserialize_from(&mut reader)?;

    let progress_bar = progress_bar(sketch_header.num_sketches as u64);
    let mut sketches: Vec<SketchType> = Vec::new();
    for _ in 0..sketch_header.num_sketches {
        let sketch: SketchType = bincode::deserialize_from(&mut reader)?;
        sketches.push(sketch);
        progress_bar.inc(1);
    }

    progress_bar.finish();

    Ok((sketch_header, sketches))
}
