use std::{
    ffi::OsStr,
    path::{Path, PathBuf},
};

use clap::{Parser, Subcommand};

use crate::sketch::{SeqFile, SKETCH_EXT};

const DEFAULT_K: u8 = 31;
const DEFAULT_SCALE: u64 = 1000;
const DEFAULT_ANI_THRESHOLD: f64 = 0.01;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(styles=get_styles())]
#[command(disable_help_subcommand = true)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,

    // hidden command used to create CLI markdown
    #[arg(long, hide = true)]
    pub markdown_help: bool,
}

#[derive(Subcommand)]
#[command(disable_help_subcommand = true)]
pub enum Commands {
    /// Create k-mer sketches for genomes
    #[command(arg_required_else_help = true)]
    Sketch(SketchArgs),

    /// Create k-mer index from sketches
    #[command(arg_required_else_help = true)]
    Index(IndexArgs),

    /// Compute distances between genome sketches
    #[command(arg_required_else_help = true)]
    Dist(DistArgs),

    /// Display information about sketch file
    #[command(arg_required_else_help = true)]
    Info(InfoArgs),

    /// Generate shell completion
    #[command(arg_required_else_help = true, hide = true)]
    ShellCompletion {
        /// Shell to generate the completions for
        #[arg(long, value_enum)]
        shell: clap_complete_command::Shell,
    },
}

#[derive(Parser)]
pub struct SketchArgs {
    /// Genome FASTA/Q file(s) to sketch
    #[arg(short = 'f', long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., group= "input")]
    pub genome_files: Option<Vec<PathBuf>>,

    /// File indicating path to genome files to sketch (TSV: genome ID followed by path to FASTA file)
    #[arg(short = 'p', long, help_heading = "Inputs", group = "input")]
    pub genome_path_file: Option<PathBuf>,

    /// Output sketch file
    #[arg(short, long, help_heading = "Output", requires = "input")]
    pub output_file: PathBuf,

    /// Length of k-mers to use
    #[arg(short, long, help_heading = "Sketching parameters", default_value_t = DEFAULT_K, value_parser = validate_kmer_length)]
    pub kmer_length: u8,

    /// Sketch scaling factor
    #[arg(short, long, help_heading = "Sketching parameters", default_value_t = DEFAULT_SCALE, value_parser = clap::value_parser!(u64).range(1..))]
    pub scale: u64,

    /// Generated sketch indicating number of times each k-mer occurs
    #[arg(short, long, help_heading = "Sketching parameters")]
    pub weighted: bool,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct IndexArgs {
    /// Genome sketches to index
    #[arg(short, long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., required=true)]
    pub sketches: Vec<PathBuf>,

    /// Output index file
    #[arg(short, long, help_heading = "Output")]
    pub output_file: PathBuf,
}

#[derive(Parser)]
pub struct DistArgs {
    /// Query genome sketches or genome FASTA/Q file(s)
    #[arg(short, long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., required=true, requires = "reference")]
    pub query_files: Vec<PathBuf>,

    /// Reference genome sketches or genome FASTA/Q file(s)
    #[arg(short, long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., group= "reference", conflicts_with = "reference_index", conflicts_with = "single_genome_set")]
    pub reference_files: Option<Vec<PathBuf>>,

    /// Reference k-mer index
    #[arg(short = 'i', long, help_heading = "Inputs", group = "reference")]
    pub reference_index: Option<PathBuf>,

    /// Output file [default: stdout]
    #[arg(short, long, help_heading = "Output")]
    pub output_file: Option<PathBuf>,

    /// Only report ANI values greater than or equal to this threshold [0, 100]
    #[arg(long, help_heading = "Output", default_value_t = DEFAULT_ANI_THRESHOLD, value_parser = validate_ani)]
    pub min_ani: f64,

    /// Report additional distance statistics
    #[arg(long, help_heading = "Output")]
    pub additional_stats: bool,

    /// Indicates query and reference genome sets are identical, and limits comparisons to lower triangle
    #[arg(long, help_heading = "Output")]
    pub single_genome_set: bool,

    /// Length of k-mers to use
    #[arg(short, long, help_heading = "Sketching parameters", default_value_t = DEFAULT_K, value_parser = validate_kmer_length)]
    pub kmer_length: u8,

    /// Sketch scaling factor
    #[arg(short, long, help_heading = "Sketching parameters", default_value_t = DEFAULT_SCALE, value_parser = clap::value_parser!(u64).range(1..))]
    pub scale: u64,

    /// Generated sketch indicating number of times each k-mer occurs
    #[arg(short, long, help_heading = "Sketching parameters")]
    pub weighted: bool,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct InfoArgs {
    /// Sketch file to query for information
    #[arg(short, long, help_heading = "Inputs")]
    pub sketch_file: PathBuf,

    /// Output file [default: stdout]
    #[arg(short, long, help_heading = "Output")]
    pub output_file: Option<PathBuf>,
}

/// Extract genome ID from input FASTA/Q file name.
pub fn genome_id_from_filename(seq_file: &Path) -> String {
    let mut genome_id = seq_file.file_name().unwrap().to_string_lossy().to_string();

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

/// Sort input files into sketches and sequence files
pub fn sort_input_files(input_files: &Vec<PathBuf>) -> (Vec<PathBuf>, Vec<SeqFile>) {
    let mut sketch_files = Vec::new();
    let mut seq_files = Vec::new();

    for input_file in input_files {
        if input_file.extension() == Some(OsStr::new(SKETCH_EXT)) {
            sketch_files.push(input_file.clone());
        } else {
            let seq_id = genome_id_from_filename(input_file);

            let seq_file = SeqFile {
                id: seq_id,
                file: input_file.clone(),
            };

            seq_files.push(seq_file);
        }
    }

    (sketch_files, seq_files)
}

fn validate_kmer_length(k: &str) -> Result<u8, String> {
    let k: u8 = k
        .parse()
        .map_err(|_| format!("`{k}` isn't a valid k-mer length"))?;

    if !(1..=32).contains(&k) {
        return Err("k-mer length must be in the range [1, 32]".to_string());
    }

    Ok(k)
}

fn validate_threads(threads: &str) -> Result<usize, String> {
    let threads: usize = threads
        .parse()
        .map_err(|_| format!("`{threads}` isn't a valid value"))?;

    if !(1..=1024).contains(&threads) {
        return Err("Threads  must be in the range [1, 1024]".to_string());
    }

    Ok(threads)
}

fn validate_ani(v: &str) -> Result<f64, String> {
    let v: f64 = v
        .parse()
        .map_err(|_| format!("`{v}` isn't a valid ANI value"))?;

    if !(0.0..=100.0).contains(&v) {
        return Err("ANI values must be in the range [0.0, 100.0]".to_string());
    }

    Ok(v)
}

fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

#[test]
fn test_verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
