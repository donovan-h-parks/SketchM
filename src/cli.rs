use clap::{Parser, Subcommand};

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
    pub genome_files: Option<Vec<String>>,

    /// File indicating path to genome files to sketch (TSV: genome ID followed by path to FASTA file)
    #[arg(short = 'p', long, help_heading = "Inputs", group = "input")]
    pub genome_path_file: Option<String>,

    /// Output sketch file
    #[arg(short, long, help_heading = "Output", requires = "input")]
    pub output_file: String,

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
    pub sketches: Vec<String>,

    /// Output index file
    #[arg(short, long, help_heading = "Output")]
    pub output_file: String,
}

#[derive(Parser)]
pub struct DistArgs {
    /// Query genome sketches or genome FASTA/Q file(s)
    #[arg(short, long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., required=true, requires = "reference")]
    pub query_files: Vec<String>,

    /// Reference genome sketches or genome FASTA/Q file(s)
    #[arg(short, long, help_heading = "Inputs", value_delimiter = ' ', num_args = 1.., group= "reference", conflicts_with = "reference_index", conflicts_with = "single_genome_set")]
    pub reference_files: Option<Vec<String>>,

    /// Reference k-mer index
    #[arg(short = 'i', long, help_heading = "Inputs", group = "reference")]
    pub reference_index: Option<String>,

    /// Output file [default: stdout]
    #[arg(short, long, help_heading = "Output")]
    pub output_file: Option<String>,

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
    pub sketch_file: String,

    /// Output file [default: stdout]
    #[arg(short, long, help_heading = "Output")]
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
