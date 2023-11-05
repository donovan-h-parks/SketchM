use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

use crate::frac_min_hash::FracMinHash;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SketchParams {
    kmer_length: u8,
    scale: u64,
    weighted: bool,
}

impl Default for SketchParams {
    fn default() -> Self {
        SketchParams {
            kmer_length: 31,
            scale: 1000,
            weighted: false,
        }
    }
}

impl SketchParams {
    pub fn new(kmer_length: u8, scale: u64, weighted: bool) -> Self {
        SketchParams {
            kmer_length,
            scale,
            weighted,
        }
    }

    pub fn create_sketcher(&self) -> FracMinHash {
        FracMinHash::new(self.kmer_length, self.scale)
    }

    pub fn k(&self) -> u8 {
        self.kmer_length
    }

    pub fn scale(&self) -> u64 {
        self.scale
    }

    pub fn weighted(&self) -> bool {
        self.weighted
    }

    /// Return true if sketch parameters are identical.
    pub fn check_compatibility(&self, other: &SketchParams) -> Result<bool> {
        if self.k() != other.k() {
            bail!(
                "Sketch has k = {}, but other sketch has k = {}",
                self.k(),
                other.k()
            );
        }

        // It is possible to define different sketch statistics so they can
        // be calculated across sketches with different scales. This generally
        // means just using the smaller of the two scales. However, this is
        // rarely the desire intent so at least for now this is being explicitly
        // checked and different scales reported as an error.
        if self.scale() != other.scale() {
            bail!(
                "Sketch has scale = {}, but other sketch has scale = {}",
                self.scale(),
                other.scale()
            );
        }

        // It is possible to convert a weighted sketch into an unweighted sketch,
        // but this conversion is not currently supported.
        if self.weighted() != other.weighted() {
            bail!(
                "Sketch has weighted = {}, but other sketch has weighted = {}",
                self.weighted(),
                other.weighted()
            );
        }

        Ok(true)
    }
}
