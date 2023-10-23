use std::collections::BTreeSet;

use needletail::parser::SequenceRecord;

use crate::hashing::{dna_hashes, ItemHash};
use crate::sketch::Sketch;
use crate::sketch_params::SketchParams;

#[derive(Clone, Debug)]
pub struct FracMinHash {
    hashes: BTreeSet<ItemHash>,
    kmer_total_count: u64,
    bp_count: u64,
    kmer_length: u8,
    scale: u64,
    max_hash: u64,
}

impl FracMinHash {
    pub fn new(kmer_length: u8, scale: u64) -> Self {
        FracMinHash {
            hashes: BTreeSet::new(),
            kmer_length,
            kmer_total_count: 0,
            bp_count: 0,
            scale,
            max_hash: ItemHash::max_value() / scale,
        }
    }

    pub fn process_seq(&mut self, seq: &SequenceRecord) {
        self.bp_count += seq.num_bases() as u64;
        self.kmer_total_count += seq.num_bases() as u64 - self.kmer_length as u64 + 1;

        dna_hashes(&seq.seq(), &mut self.hashes, self.max_hash, self.kmer_length);
    }

    pub fn hashes_in_sketch(&self) -> u64 {
        self.hashes.len() as u64
    }

    pub fn kmer_total_count(&self) -> u64 {
        self.kmer_total_count
    }

    pub fn bp_count(&self) -> u64 {
        self.bp_count
    }

    pub fn parameters(&self) -> SketchParams {
        SketchParams::new(self.kmer_length, self.scale)
    }

    pub fn to_sketch(self, name: &str) -> Sketch {
        Sketch {
            name: name.to_string(),
            bp_count: self.bp_count(),
            kmer_total_count: self.kmer_total_count,
            hashes: self.hashes.into_iter().collect(),
        }
    }
}
