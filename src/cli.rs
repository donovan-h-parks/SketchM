use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(disable_help_subcommand = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Create k-mer sketches for genomes
    Sketch(SketchArgs),

    /// Compute distances between genome sketches
    Dist(DistArgs),

    /// Create k-mer index from sketches
    Index(IndexArgs),

    /// Compute distances using k-mer index
    DistByIndex(DistByIndexArgs),

    /// Display information about sketch file
    Info(InfoArgs),
}

#[derive(Parser)]
pub struct SketchArgs {
    /// Genome FASTA/Q file(s) to sketch
    #[arg(short = 'f', long, value_delimiter = ' ', num_args = 1..)]
    pub genome_files: Option<Vec<String>>,

    /// File indicating path to genome files to sketch (TSV: genome ID followed by path to FASTA file)
    #[arg(short = 'p', long)]
    pub genome_path_file: Option<String>,

    /// Output sketch file
    #[arg(short, long, required = true)]
    pub output_file: String,

    /// Generated sketch indicating number of times each k-mer occurs
    #[arg(short, long)]
    pub weighted: bool,

    /// Length of k-mers to use
    #[arg(short, long, default_value_t = 31, value_parser = validate_kmer_length)]
    pub kmer_length: u8,

    /// Sketch scaling factor
    #[arg(short, long, default_value_t = 1000)]
    pub scale: u64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct DistArgs {
    /// Query genome sketches
    #[arg(short, long, value_delimiter = ' ', num_args = 1..)]
    pub query_sketches: Vec<String>,

    /// Reference genome sketches
    #[arg(short, long, value_delimiter = ' ', num_args = 1..)]
    pub reference_sketches: Vec<String>,

    /// Output file [default: stdout]
    #[arg(short, long)]
    pub output_file: Option<String>,

    /// Only report ANI values above this threshold [0, 100]
    #[arg(long, default_value_t = 0.01, value_parser = validate_ani)]
    pub min_ani: f64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct IndexArgs {
    /// Genome sketches to index
    #[arg(short, long, value_delimiter = ' ', num_args = 1..)]
    pub sketches: Vec<String>,

    /// Output index file
    #[arg(short, long)]
    pub output_file: String,
}

#[derive(Parser)]
pub struct DistByIndexArgs {
    /// Query genome sketches
    #[arg(short, long, value_delimiter = ' ', num_args = 1..)]
    pub query_sketches: Vec<String>,

    /// Reference k-mer index
    #[arg(short, long, value_delimiter = ' ', num_args = 1..)]
    pub reference_index: String,

    /// Output file [default: stdout]
    #[arg(short, long)]
    pub output_file: Option<String>,

    /// Only report ANI values above this threshold [0, 100]
    #[arg(long, default_value_t = 0.01, value_parser = validate_ani)]
    pub min_ani: f64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct InfoArgs {
    /// Sketch file to query for information
    #[arg(short, long)]
    pub sketch_file: String,

    /// Output file [default: stdout]
    #[arg(short, long)]
    pub output_file: Option<String>,
}

fn validate_kmer_length(k: &str) -> Result<u8, String> {
    let k: u8 = k
        .parse()
        .map_err(|_| format!("`{k}` isn't a valid k-mer length"))?;

    if !(0..=32).contains(&k) {
        return Err("K-mer length must be in the range [1, 32]".to_string());
    }

    Ok(k)
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

#[test]
fn test_verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
