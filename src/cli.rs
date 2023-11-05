use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
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
pub enum Commands {
    /// Create k-mer sketches for genomes.
    #[command(arg_required_else_help = true)]
    Sketch(SketchArgs),

    /// Create k-mer index from sketches.
    #[command(arg_required_else_help = true)]
    Index(IndexArgs),

    /// Compute distances between genome sketches.
    #[command(arg_required_else_help = true)]
    Dist(DistArgs),

    /// Display information about sketch file.
    #[command(arg_required_else_help = true)]
    Info(InfoArgs),
}

#[derive(Parser)]
pub struct SketchArgs {
    /// Genome FASTA/Q file(s) to sketch
    #[arg(short = 'f', long, value_delimiter = ' ', num_args = 1.., group= "input")]
    pub genome_files: Option<Vec<String>>,

    /// File indicating path to genome files to sketch (TSV: genome ID followed by path to FASTA file)
    #[arg(short = 'p', long, group = "input")]
    pub genome_path_file: Option<String>,

    /// Output sketch file
    #[arg(short, long, requires = "input")]
    pub output_file: String,

    /// Generated sketch indicating number of times each k-mer occurs
    #[arg(short, long)]
    pub weighted: bool,

    /// Length of k-mers to use
    #[arg(short, long, default_value_t = 31, value_parser = validate_kmer_length)]
    pub kmer_length: u8,

    /// Sketch scaling factor
    #[arg(short, long, default_value_t = 1000, value_parser = clap::value_parser!(u64).range(1..))]
    pub scale: u64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct IndexArgs {
    /// Genome sketches to index
    #[arg(short, long, value_delimiter = ' ', num_args = 1.., required=true)]
    pub sketches: Vec<String>,

    /// Output index file
    #[arg(short, long)]
    pub output_file: String,
}

#[derive(Parser)]
pub struct DistArgs {
    /// Query genome sketches
    #[arg(short, long, value_delimiter = ' ', num_args = 1.., required=true, requires = "reference")]
    pub query_sketches: Vec<String>,

    /// Reference genome sketches
    #[arg(short, long, value_delimiter = ' ', num_args = 1.., group= "reference", conflicts_with = "reference_index")]
    pub reference_sketches: Option<Vec<String>>,

    /// Reference k-mer index
    #[arg(short = 'i', long, group = "reference")]
    pub reference_index: Option<String>,

    /// Output file [default: stdout]
    #[arg(short, long)]
    pub output_file: Option<String>,

    /// Only report ANI values above this threshold [0, 100]
    #[arg(long, default_value_t = 0.01, value_parser = validate_ani)]
    pub min_ani: f64,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 1, value_parser = validate_threads)]
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

#[test]
fn test_verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
