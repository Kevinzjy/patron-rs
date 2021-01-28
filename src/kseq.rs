use log::{info, error};

use debruijn::Kmer;
use debruijn::dna_string::*;
use debruijn::kmer::{Kmer10, Kmer15};
use debruijn::Vmer;

pub fn split_reads(seq: DnaString) -> Vec<Kmer10> {
    let mask = (1 << 2 * 10) - 1;
    let mut markers: Vec<Kmer10> = Vec::new();

    for offset in 0..seq.len()-15+1 {
        let kmer: Kmer15 = seq.get_kmer(offset);
        let minimizer = get_minimizer(kmer, &mask);
        markers.push(minimizer);
    }
    markers
}

pub fn get_minimizer(kmer: Kmer15, mask: &u64) -> Kmer10 {
    let kseq = kmer.to_string();
    let kstr = DnaString::from_dna_string(&kseq);

    let mut minimizer: Kmer10 = kstr.get_kmer(0);
    let mut minhash: u64 = hash64(minimizer.to_u64(), mask);
    for i in 1..15-10+1 {
        let cand: Kmer10 = kstr.get_kmer(i);
        let cand_hash = hash64(cand.to_u64(), mask) << 8 | 10;
        if cand_hash < minhash {
            minhash = cand_hash;
            minimizer = cand;
        }
    }
    minimizer
}

/// Hash function from minimap2
/// https://github.com/lh3/minimap2/blob/master/sketch.c
fn hash64(kmer: u64, mask: &u64) -> u64 {
    let mut key = kmer;
    key = (!key).wrapping_add(key<<21) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key<<3)) + (key<<8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key<<2)) + (key<<4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key<<31)) & mask;
    key
}

#[test]
fn test_hash64() {
    let mask = (1 << 2 * 15) - 1;
    let test1 = DnaString::from_dna_string("AAAAAAACCCTTTTT");
    let test2 = DnaString::from_dna_string("AAAAAGGGTTTTTTT");
    let kmer1: Kmer15 = test1.get_kmer(0);
    let kmer2: Kmer15 = test1.get_kmer(0);
    let hash1 = hash64(kmer1.to_u64(), &mask);
    let hash2 = hash64(kmer1.to_u64(), &mask);
    assert_eq! (hash1, hash2);
}
