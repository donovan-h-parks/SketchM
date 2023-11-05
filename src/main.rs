use std::env;
use std::fs::File;
use std::io::stdout;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use log::info;
use maybe_gzip_io::maybe_gzip_reader;
use num_format::ToFormattedString;
use serde::Serialize;

use crate::cli::{Cli, Commands};
use crate::config::LOCALE;
use crate::distance::{
    calc_sketch_distances, calc_sketch_distances_to_index, calc_weighted_sketch_distances,
};
use crate::genome::read_genome_path_file;
use crate::index_sketches::index_sketches;
use crate::logging::setup_logger;
use crate::sketch::{read_sketch, sketch, SeqFile, SketchHeader, SKETCH_EXT};
use crate::sketch_params::SketchParams;

mod cli;
pub mod config;
pub mod distance;
pub mod frac_min_hash;
pub mod genome;
pub mod hashing;
pub mod index_sketches;
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
                serde_json::to_writer_pretty(writer, &data)
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

    // calculate distances between reference sketches or index
    if let Some(ref_sketches) = &args.reference_sketches {
        // read sketch header to determine if we should
        // calculate unweighted or weighted distances
        let mut reader = maybe_gzip_reader(&args.query_sketches[0]).context(format!(
            "Unable to read sketch file: {}",
            &args.query_sketches[0]
        ))?;

        let sketch_header: SketchHeader = bincode::deserialize_from(&mut reader)?;

        if sketch_header.params.weighted() {
            info!("Calculating weighted distances between sketches.");
            calc_weighted_sketch_distances(
                &args.query_sketches,
                ref_sketches,
                args.min_ani,
                &args.output_file,
                args.threads,
            )?;
        } else {
            info!("Calculating unweighted distances between sketches.");
            calc_sketch_distances(
                &args.query_sketches,
                ref_sketches,
                args.min_ani,
                &args.output_file,
                args.threads,
            )?;
        }
    } else if let Some(ref_index) = &args.reference_index {
        calc_sketch_distances_to_index(
            &args.query_sketches,
            ref_index,
            args.min_ani,
            &args.output_file,
            args.threads,
        )?;
    }

    Ok(())
}

/// Run the index command.
fn run_index(args: &cli::IndexArgs) -> Result<()> {
    init(1)?;

    index_sketches(&args.sketches, &args.output_file)?;

    Ok(())
}

/// Run the sketch info command.
fn run_info(args: &cli::InfoArgs) -> Result<()> {
    init(1)?;

    let (sketch_header, sketches) = read_sketch(&args.sketch_file)?;

    #[derive(Serialize)]
    struct SketchInfo {
        name: String,
        bp_count: u64,
        kmer_total_count: u64,
        unique_hash_count: u64,
        weighted_hash_count: u64,
        k: u8,
        scale: u64,
        weighted: bool,
    }

    let mut sketch_info = Vec::new();

    let mut genome_size_mean = 0.0f32;
    let mut unique_hash_count_mean = 0.0f32;
    let mut weighted_hash_count_mean = 0.0f32;
    let scale_factor = 1.0 / sketches.len() as f32;
    for sketch in sketches.into_iter() {
        sketch_info.push(SketchInfo {
            name: sketch.name().to_string(),
            bp_count: sketch.bp_count(),
            kmer_total_count: sketch.kmer_total_count(),
            unique_hash_count: sketch.unique_hash_count(),
            weighted_hash_count: sketch.weighted_hash_count(),
            k: sketch_header.params.k(),
            scale: sketch_header.params.scale(),
            weighted: sketch_header.params.weighted(),
        });

        genome_size_mean += scale_factor * sketch.bp_count() as f32;
        unique_hash_count_mean += scale_factor * sketch.unique_hash_count() as f32;
        weighted_hash_count_mean += scale_factor * sketch.weighted_hash_count() as f32;
    }

    info!("K: {}", sketch_header.params.k());
    info!("Scale: {}", sketch_header.params.scale());
    info!("Weighted: {}", sketch_header.params.weighted());

    info!("Average genome size: {:.1}", genome_size_mean);
    info!("Average unique hash count: {:.1}", unique_hash_count_mean);
    info!(
        "Average weighted hash count: {:.1}",
        weighted_hash_count_mean
    );

    output_results(&sketch_info, &args.output_file)?;

    Ok(())
}

fn main() -> Result<()> {
    let start = Instant::now();

    setup_logger();

    let cli = Cli::parse();

    if cli.markdown_help {
        clap_markdown::print_help_markdown::<Cli>();
    } else {
        match &cli.command {
            Some(Commands::Sketch(args)) => run_sketch(args)?,
            Some(Commands::Dist(args)) => run_dist(args)?,
            Some(Commands::Index(args)) => run_index(args)?,
            Some(Commands::Info(args)) => run_info(args)?,
            None => {}
        }

        info!("Elapsed time (sec): {:.2}", start.elapsed().as_secs_f32());
        info!("Done.");
    }

    Ok(())
}
