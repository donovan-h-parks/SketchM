// Methods in this file were ported from minimap2 by Heng Li or modified from skani by Jim Shaw.

// minimap2 MIT License
//
// Copyright (c) 2018-     Dana-Farber Cancer Institute
//               2017-2018 Broad Institute, Inc.
//
// skani MIT License
//
// Copyright (c) 2022 Jim Shaw
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

use std::collections::BTreeMap;

pub type ItemHash = u64;

const NT_TO_BYTE: [u8; 256] = {
    let mut table = [0; 256];

    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;

    table
};

/// Thomas Wang's integer hash function.
// Ported from minimap2 and following Rust implementation by Anicet Ebou.
// https://gist.github.com/lh3/974ced188be2f90422cc#file-inthash-c
// https://aebou.rbind.io/post/a-rust-glimpse-at-thomas-wang-integer-hash-function
// Further reading: https://gist.github.com/badboy/6267743
#[inline]
pub fn tw_hash64(kmer: ItemHash) -> ItemHash {
    let mut hash = kmer;

    hash = (!hash).wrapping_add(hash << 21); // key = (key << 21) - key - 1
    hash = hash ^ (hash >> 24);

    hash = hash.wrapping_add(hash << 3).wrapping_add(hash << 8); // key * 265
    hash = hash ^ (hash >> 14);

    hash = hash.wrapping_add(hash << 2).wrapping_add(hash << 4); // key * 21
    hash = hash ^ (hash >> 28);

    hash = hash.wrapping_add(hash << 31);

    hash
}

/// Determine hashes in sequence satisfying maximum k-mer hash criterion.
// Modified from the fmh_seeds method by Jim Shaw in skani.
pub fn dna_hashes(seq: &[u8], hashes: &mut BTreeMap<ItemHash, u32>, max_hash: ItemHash, k: u8) {
    let k = k as usize;

    if seq.len() < k {
        return;
    }

    let mut fwd_kmer: ItemHash = 0;
    let mut rev_kmer: ItemHash = 0;

    let rev_shift_dist = 2 * (k - 1);
    let fwd_mask = ItemHash::MAX >> (std::mem::size_of::<ItemHash>() * 8 - 2 * k);

    for i in 0..k - 1 {
        let nuc_f = NT_TO_BYTE[seq[i] as usize] as ItemHash;
        fwd_kmer <<= 2;
        fwd_kmer |= nuc_f;

        let nuc_r = 3 - nuc_f;
        rev_kmer >>= 2;
        rev_kmer |= nuc_r << rev_shift_dist;
    }

    for i in k - 1..seq.len() {
        let nuc_f = NT_TO_BYTE[seq[i] as usize] as ItemHash;
        fwd_kmer <<= 2;
        fwd_kmer |= nuc_f;
        fwd_kmer &= fwd_mask;

        let nuc_r = 3 - nuc_f;
        rev_kmer >>= 2;
        rev_kmer |= nuc_r << rev_shift_dist;

        let canonical_kmer_marker = if fwd_kmer < rev_kmer {
            fwd_kmer
        } else {
            rev_kmer
        };

        let hash = tw_hash64(canonical_kmer_marker);
        if hash < max_hash {
            let count = hashes.entry(hash).or_insert(0);
            *count = count.saturating_add(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Inverse of Thomas Wang's integer hash function.
    // https://aebou.rbind.io/post/a-rust-glimpse-at-thomas-wang-integer-hash-function
    fn tw_hash64i(hashed_key: u64) -> u64 {
        let mut key = hashed_key;

        // Invert h_key = h_key.wrapping_add(h_key << 31)
        let mut tmp: u64 = key.wrapping_sub(key << 31);
        key = key.wrapping_sub(tmp << 31);

        // Invert h_key = h_key ^ h_key >> 28;
        tmp = key ^ key >> 28;
        key ^= tmp >> 28;

        // Invert h_key = h_key.wrapping_add(h_key << 2).wrapping_add(h_key << 4)
        key = key.wrapping_mul(14933078535860113213u64);

        // Invert h_key = h_key ^ h_key >> 14;
        tmp = key ^ key >> 14;
        tmp = key ^ tmp >> 14;
        tmp = key ^ tmp >> 14;
        key ^= tmp >> 14;

        // Invert h_key = h_key.wrapping_add(h_key << 3).wrapping_add(h_key << 8)
        key = key.wrapping_mul(15244667743933553977u64);

        // Invert h_key = h_key ^ h_key >> 24
        tmp = key ^ key >> 24;
        key ^= tmp >> 24;

        // Invert h_key = (!h_key).wrapping_add(h_key << 21)
        tmp = !key;
        tmp = !(key.wrapping_sub(tmp << 21));
        tmp = !(key.wrapping_sub(tmp << 21));
        key = !(key.wrapping_sub(tmp << 21));

        key
    }

    #[test]
    fn test_hashing() {
        assert_eq!(tw_hash64i(tw_hash64(0)), 0);
        assert_eq!(tw_hash64i(tw_hash64(u64::MAX)), u64::MAX);
        assert_eq!(tw_hash64i(tw_hash64(27)), 27);
        assert_eq!(tw_hash64i(tw_hash64(108)), 108);
        assert_eq!(tw_hash64i(tw_hash64(177)), 177);
    }

    #[test]
    fn test_bit_kmer_value() {
        let mut hashes = BTreeMap::new();

        dna_hashes(b"AAAA", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(tw_hash64(0), *hashes.iter().next().unwrap().0); // AAAA = 00000000b = 0

        hashes.clear();
        dna_hashes(b"TTTT", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(tw_hash64(0), *hashes.iter().next().unwrap().0); // AAAA = 0 < TTTT

        hashes.clear();
        dna_hashes(b"CCCC", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(tw_hash64(85), *hashes.iter().next().unwrap().0); // CCCC = 01010101b = 85

        hashes.clear();
        dna_hashes(b"GGGG", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(tw_hash64(85), *hashes.iter().next().unwrap().0); // CCCC = 85 < GGGG
    }

    #[test]
    fn test_canonical_kmer() {
        let mut hashes = BTreeMap::new();
        dna_hashes(b"AAAAAAAA", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 5);
        assert_eq!(hashes.get(&tw_hash64(0)), Some(&5));

        let mut hashes = BTreeMap::new();
        dna_hashes(b"TTTTTTTT", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 5);
        assert_eq!(hashes.get(&tw_hash64(0)), Some(&5));

        let mut hashes = BTreeMap::new();
        dna_hashes(b"CCCCCCCC", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 5);
        assert_eq!(hashes.get(&tw_hash64(85)), Some(&5));

        let mut hashes = BTreeMap::new();
        dna_hashes(b"GGGGGGGG", &mut hashes, u64::MAX, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 5);
        assert_eq!(hashes.get(&tw_hash64(85)), Some(&5));
    }

    #[test]
    fn test_simple_seq() {
        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, u64::MAX, 4);

        // kmer | rev  | smallest | binary    | decimal
        // ACGT | ACGT | ACGT     | 00011011b | 27
        // CGTA | TACG | CGTA     | 01101100b | 108
        // GTAC | GTAC | GTAC     | 10110001b | 177
        // TACG | CGTA | CGTA     | 01101100b | 108
        // ACGT | ACGT | ACGT     | 00011011b | 27

        assert_eq!(hashes.len(), 3);
        assert_eq!(hashes.values().sum::<u32>(), 5);

        assert_eq!(hashes.get(&tw_hash64(27)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(108)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(177)), Some(&1));
    }

    #[test]
    fn test_mixed_case_seq() {
        // should produce same restul as test_simple_seq()

        let mut hashes = BTreeMap::new();
        dna_hashes(b"AcgTaCGt", &mut hashes, u64::MAX, 4);

        assert_eq!(hashes.len(), 3);
        assert_eq!(hashes.values().sum::<u32>(), 5);

        assert_eq!(hashes.get(&tw_hash64(27)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(108)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(177)), Some(&1));
    }

    #[test]
    fn test_ambiguous_bases() {
        // Ambiguous bases are treated as an A; a common convention
        // in high performance bioinformatic software. Should produce
        // same results as test_simple_seq().

        let mut hashes = BTreeMap::new();
        dna_hashes(b"NCGTnCGT", &mut hashes, u64::MAX, 4);

        assert_eq!(hashes.len(), 3);
        assert_eq!(hashes.values().sum::<u32>(), 5);

        assert_eq!(hashes.get(&tw_hash64(27)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(108)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(177)), Some(&1));
    }

    #[test]
    fn test_filter_hashes() {
        // kmer | rev  | smallest | binary    | decimal | hash
        // ACGT | ACGT | ACGT     | 00011011b | 27      | 12564563040126408309
        // CGTA | TACG | CGTA     | 01101100b | 108     | 13364770925836396135
        // GTAC | GTAC | GTAC     | 10110001b | 177     | 8958356766268387398
        // TACG | CGTA | CGTA     | 01101100b | 108     | 13364770925836396135
        // ACGT | ACGT | ACGT     | 00011011b | 27      | 12564563040126408309

        assert_eq!(tw_hash64(27), 12564563040126408309);
        assert_eq!(tw_hash64(108), 13364770925836396135);
        assert_eq!(tw_hash64(177), 8958356766268387398);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 8958356766268387398, 4);
        assert_eq!(hashes.len(), 0);
        assert_eq!(hashes.values().sum::<u32>(), 0);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 8958356766268387398 + 1, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 1);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 12564563040126408309, 4);
        assert_eq!(hashes.len(), 1);
        assert_eq!(hashes.values().sum::<u32>(), 1);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 12564563040126408309 + 1, 4);
        assert_eq!(hashes.len(), 2);
        assert_eq!(hashes.values().sum::<u32>(), 3);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 13364770925836396135, 4);
        assert_eq!(hashes.len(), 2);
        assert_eq!(hashes.values().sum::<u32>(), 3);

        let mut hashes = BTreeMap::new();
        dna_hashes(b"ACGTACGT", &mut hashes, 13364770925836396135 + 1, 4);
        assert_eq!(hashes.len(), 3);
        assert_eq!(hashes.values().sum::<u32>(), 5);
    }

    #[test]
    fn test_rev_comp_hashes() {
        // test that reverse complement hashes are being properly calculated
        // and selected for odd k
        //
        // kmer | rev | smallest | binary  | decimal | hash
        // ACG  | CGT | ACG      | 000110b | 6       | 12564563040126408309
        // CGT  | ACG | ACG      | 000110b | 6       | 12564563040126408309
        // GTT  | AAC | AAC      | 000001b | 1       | 8958356766268387398
        // TTT  | AAA | AAA      | 000000b | 0       | 13364770925836396135

        let mut hashes = BTreeMap::new();

        dna_hashes(b"ACGTTT", &mut hashes, u64::MAX, 3);
        assert_eq!(hashes.len(), 3);
        assert_eq!(hashes.values().sum::<u32>(), 4);

        assert_eq!(hashes.get(&tw_hash64(6)), Some(&2));
        assert_eq!(hashes.get(&tw_hash64(1)), Some(&1));
        assert_eq!(hashes.get(&tw_hash64(0)), Some(&1));
    }
}
