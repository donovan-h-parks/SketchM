use std::collections::BTreeMap;

use needletail::parser::SequenceRecord;

use crate::hashing::{dna_hashes, ItemHash};
use crate::sketch::{Sketch, WeightedSketch};

#[derive(Clone, Debug)]
pub struct FracMinHash {
    hashes: BTreeMap<ItemHash, u32>,
    kmer_length: u8,
    max_hash: u64,
    kmer_total_count: u64,
    bp_count: u64,
}

impl FracMinHash {
    pub fn new(kmer_length: u8, scale: u64) -> Self {
        FracMinHash {
            hashes: BTreeMap::new(),
            kmer_length,
            max_hash: ItemHash::max_value() / scale,
            kmer_total_count: 0,
            bp_count: 0,
        }
    }

    pub fn process_seq(&mut self, seq: &SequenceRecord) {
        self.bp_count += seq.num_bases() as u64;
        self.kmer_total_count += seq.num_bases() as u64 - self.kmer_length as u64 + 1;

        dna_hashes(
            &seq.seq(),
            &mut self.hashes,
            self.max_hash,
            self.kmer_length,
        );
    }

    pub fn unique_hash_count(&self) -> u64 {
        self.hashes.len() as u64
    }

    pub fn weighted_hash_count(&self) -> u64 {
        self.hashes.values().map(|v| *v as u64).sum()
    }

    pub fn kmer_total_count(&self) -> u64 {
        self.kmer_total_count
    }

    pub fn bp_count(&self) -> u64 {
        self.bp_count
    }

    pub fn to_sketch(self, name: &str) -> Sketch {
        Sketch {
            name: name.to_string(),
            bp_count: self.bp_count(),
            kmer_total_count: self.kmer_total_count,
            hashes: self.hashes.into_keys().collect(),
        }
    }

    pub fn to_weighted_sketch(self, name: &str) -> WeightedSketch {
        WeightedSketch {
            name: name.to_string(),
            bp_count: self.bp_count(),
            kmer_total_count: self.kmer_total_count,
            hashes: self.hashes,
        }
    }
}
