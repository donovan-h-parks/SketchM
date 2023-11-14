use std::{
    path::{Path, PathBuf},
    str,
};

use ahash::{HashMap, HashMapExt};
use anyhow::{Context, Result};
use serde::Deserialize;

/// Parse TSV file indicating genome IDs and path to genomic FASTA/Q file for genome.
pub fn read_genome_path_file(input_path: &Path) -> Result<HashMap<String, PathBuf>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(input_path)
        .context(format!("Unable to open file: {}", input_path.display()))?;

    #[derive(Deserialize)]
    pub struct GenomeFile {
        pub id: String,
        pub file: PathBuf,
    }

    let mut genome_paths = HashMap::new();
    for record in reader.deserialize() {
        let genome_path: GenomeFile = record?;
        genome_paths.insert(genome_path.id, genome_path.file);
    }

    Ok(genome_paths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_genome_path_file() {
        let genome_path_file = Path::new("tests/data/genome_paths.tsv");
        let genome_paths =
            read_genome_path_file(&genome_path_file).expect("Invalid genome path file");

        assert_eq!(genome_paths.len(), 2);
        assert_eq!(
            genome_paths.get("Genome1"),
            Some(&PathBuf::from("./tests/data/genome1.fna.gz"))
        );
        assert_eq!(
            genome_paths.get("Genome2"),
            Some(&PathBuf::from("./tests/data/genome2.fna.gz"))
        );
    }
}
