use std::env;
use std::fs::File;
use std::io::stdout;
use std::path::Path;
use std::time::Instant;
use std::io::Write;

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use log::{info};
use num_format::ToFormattedString;
use serde::Serialize;

use crate::cli::{Cli, Commands};
use crate::config::LOCALE;
use crate::genome::read_genome_path_file;
use crate::logging::setup_logger;
use crate::sketch::{sketch, read_sketch, SeqFile, SKETCH_EXT};
use crate::sketch_params::SketchParams;
use crate::distance::calc_sketch_distances;

mod cli;
pub mod config;
pub mod distance;
pub mod frac_min_hash;
pub mod genome;
pub mod hashing;
pub mod logging;
pub mod maybe_gzip_io;
pub mod progress;
pub mod sketch;
pub mod sketch_params;


/// Common initialization required by all commands.
fn init(threads: usize) -> Result<()> {
    const VERSION: &str = env!("CARGO_PKG_VERSION");
    info!("{} v{}", env!("CARGO_PKG_NAME"), VERSION);
    info!("{}", env::args().collect::<Vec<String>>().join(" "));

    info!("Using {} threads.", threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    Ok(())
}

/// Extract genome ID from input FASTA/Q file name.
pub fn genome_id_from_filename(seq_file: &str) -> String {
    let mut genome_id = Path::new(seq_file)
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();

    if genome_id.ends_with(".gz") {
        genome_id = genome_id.replace(".gz", "");
    }

    if genome_id.ends_with(".fq") {
        genome_id = genome_id.replace(".fq", "");
    } else if genome_id.ends_with(".fna") {
        genome_id = genome_id.replace(".fna", "");
    } else if genome_id.ends_with(".fa") {
        genome_id = genome_id.replace(".fa", "");
    } else if genome_id.ends_with(".fasta") {
        genome_id = genome_id.replace(".fasta", "");
    } else if genome_id.ends_with(".fastq") {
        genome_id = genome_id.replace(".fastq", "");
    }

    genome_id
}

/// Get writer to stdout or to specified file with format indicate by file extension.
fn output_results<T>(data: &[T], output_file: &Option<String>) -> Result<()>
where
    T: Serialize,
{
    if let Some(output_file) = output_file {
        let path = Path::new(output_file);
        let writer = File::create(path)
            .context(format!("Unable to create: {output_file}"))
            .unwrap();

        match path.extension() {
            None => {
                serde_json::to_writer_pretty(stdout(), &data)
                    .map_err(|_| anyhow!("Could not serialize JSON to file"))?;
            }
            Some(ext) => match ext.to_str() {
                Some("tsv") | Some("csv") => {
                    let mut wrt = if ext.to_str() == Some("tsv") {
                        csv::WriterBuilder::new()
                            .delimiter(b'\t')
                            .from_writer(writer)
                    } else {
                        csv::WriterBuilder::new().from_writer(writer)
                    };

                    for record in data {
                        wrt.serialize(record)?;
                    }
                    wrt.flush()?;
                }
                _ => {
                    // default to JSON output
                    serde_json::to_writer_pretty(writer, &data)
                        .map_err(|_| anyhow!("Could not serialize JSON to file"))?;
                }
            },
        }
    } else {
        // default to JSON when writing to stdout
        serde_json::to_writer_pretty(stdout(), &data)
            .map_err(|_| anyhow!("Could not serialize JSON to file"))?;
        writeln!(stdout())?; // write newline after JSON
        stdout().flush()?;
    }

    Ok(())
}

/// Run the sketch command.
fn run_sketch(args: &cli::SketchArgs) -> Result<()> {
    init(args.threads)?;

    // get genome files to sketch
    let mut input_genome_files = Vec::new();
    if let Some(genome_path_file) = &args.genome_path_file {
        info!("Reading path to genomic FASTA files:");

        let genome_files = read_genome_path_file(genome_path_file)?;

        for (gid, genome_file) in genome_files.into_iter() {
            let seq_file = SeqFile {
                id: gid,
                file: genome_file,
            };
            input_genome_files.push(seq_file);
        }

        info!(
            " - identified path to {} genomes",
            input_genome_files.len().to_formatted_string(&LOCALE)
        );
    }

    if let Some(genome_files) = &args.genome_files {
        for genome_file in genome_files {
            let seq_id = genome_id_from_filename(genome_file);

            let seq_file = SeqFile {
                id: seq_id,
                file: genome_file.clone(),
            };

            input_genome_files.push(seq_file);
        }
    }

    let mut out_file = args.output_file.clone();
    if !out_file.ends_with(SKETCH_EXT) {
        out_file += SKETCH_EXT;
    };

    if let Some(out_path) = Path::new(&out_file).parent() {
        std::fs::create_dir_all(out_path)?;
    }

    let sketch_params = SketchParams::new(args.kmer_length, args.scale, args.weighted);

    sketch(&input_genome_files, &sketch_params, &out_file)?;

    Ok(())
}

/// Run the distance command.
fn run_dist(args: &cli::DistArgs) -> Result<()> {
    init(args.threads)?;

    calc_sketch_distances(&args.query_sketches, &args.reference_sketches, args.min_ani, &args.output_file)?;

    Ok(())
}

/// Run the sketch info command.
fn run_info(args: &cli::InfoArgs) -> Result<()> {
    init(1)?;

    let (sketch_header, sketches) = read_sketch(&args.sketch_file)?;

    # [derive(Serialize)]
    struct SketchInfo {
        name: String,
        bp_count: u64,
        kmer_total_count: u64,
        sketch_hash_count: u64,
        k: u8,
        scale: u64,
    }

    let mut sketch_info = Vec::new();

    let mut genome_size_mean = 0.0f32;
    let mut hash_count_mean = 0.0f32;
    let scale_factor = 1.0 / sketches.len() as f32;
    for sketch in sketches {
        sketch_info.push(SketchInfo {
            name: sketch.name().to_string(),
            bp_count: sketch.bp_count(),
            kmer_total_count: sketch.kmer_total_count(),
            sketch_hash_count: sketch.unique_hash_count(),
            k: sketch_header.params.k(),
            scale: sketch_header.params.scale(),
        });

        genome_size_mean += scale_factor * sketch.bp_count() as f32;
        hash_count_mean += scale_factor * sketch.unique_hash_count() as f32;
    }
    info!("Average genome size: {:.1}", genome_size_mean);
    info!("Average hash count: {:.1}", hash_count_mean);

    output_results(&sketch_info, &args.output_file)?;

    Ok(())
}

fn main() -> Result<()> {
    let start = Instant::now();

    setup_logger();

    let cli = Cli::parse();
    match &cli.command {
        Commands::Sketch(args) => run_sketch(args)?,
        Commands::Dist(args) => run_dist(args)?,
        Commands::Info(args) => run_info(args)?,
    }

    info!("Elapsed time (sec): {:.2}", start.elapsed().as_secs_f32());
    info!("Done.");

    Ok(())
}
