use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use anyhow::{Context, Result};
use bincode;
use log::info;
use needletail::parse_fastx_reader;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::hashing::ItemHash;
use crate::progress::progress_bar;
use crate::maybe_gzip_io::maybe_gzip_reader;
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
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }
}

/// Write sketch header information.
fn write_sketch_header<W>(writer: &mut W, num_sketches: u32, sketch_params: &SketchParams) -> Result<()>
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
pub fn sketch_file(seq_file: &SeqFile, sketch_params: &SketchParams) -> Result<Sketch> {
    let mut sketcher = sketch_params.create_sketcher();
    let reader = File::open(Path::new(&seq_file.file))
        .context(format!("Failed to open {}", seq_file.file))?;

    let mut fastx_reader = parse_fastx_reader(reader)?;
    while let Some(rec) = fastx_reader.next() {
        let record = rec?;
        sketcher.process_seq(&record);
    }

    Ok(sketcher.to_sketch(&seq_file.id))
}

/// Read sketches.
pub fn read_sketch(sketch_file: &str) -> Result<(SketchHeader, Vec<Sketch>)> {
    let mut reader = maybe_gzip_reader(sketch_file).context(format!(
        "Unable to read sketch file: {}",
        sketch_file
    ))?;

    let sketch_header: SketchHeader = bincode::deserialize_from(&mut reader)?;

    let progress_bar = progress_bar(sketch_header.num_sketches as u64);
    let mut sketches: Vec<Sketch> = Vec::new();
    for _ in 0..sketch_header.num_sketches {
        let sketch: Sketch = bincode::deserialize_from(&mut reader)?;
        sketches.push(sketch);
        progress_bar.inc(1);
    }

    progress_bar.finish();

    Ok((sketch_header, sketches))
}
