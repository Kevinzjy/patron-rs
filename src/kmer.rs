use std::str;
use std::mem;

use byteorder::{ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};

enum Base {
    A,
    T,
    C,
    G,
}

pub fn get_minimizers(seq: &str) {
    let k = 5;
    let w = 15;
    let shift1 = 2 * (k - 1);
    let mask = (1 << 2 * k) - 1;

    let mut kmer = [0; 2];
    let mut minimizers: Vec<(u64, usize)> = Vec::with_capacity(seq.len());
    let mut kmers: Vec<(u64, usize)> = Vec::with_capacity(seq.len());

    let bseq = seq.as_bytes();
    let mut l = 0;
    let mut tmp_min = (mask, 0);
    for i in 0..bseq.len() {
        let c = match bseq[i] {
            65 => 0, // A
            67 => 1, // C
            71 => 2, // G
            84 => 3, // T
            _ => continue,
        };
        
        // Strand
        kmer[0] = (kmer[0] << 2 | c) & mask;
        kmer[1] = (kmer[1] >> 2) | (3^c) << shift1;
        if kmer[0] == kmer[1] {
            continue; // Skip symmetry kmers
        }
        let z = (kmer[0] < kmer[1]) as usize;
        let info = (hash64(kmer[z], &mask), i);

        // Skip first kmers
        if info.1 < k {
            continue;
        }

        if info.1 < w - k  {
            if info.0 <= tmp_min.0 {
                tmp_min = info;
            }
        } else {
            if info.1 > tmp_min.1 + w - k {
                let st = tmp_min.1 + 1 - w;
                let en = i - w;
                println! ("{} {:?} {:?}", i, tmp_min, info, );
                tmp_min = kmers[tmp_min.1 + 1 - w];
                for j in st..en {
                    if kmers[j].0 <= tmp_min.0 {
                        tmp_min = kmers[j];
                    }
                }
                if info.0 <= tmp_min.0 {
                    minimizers.push(info);
                    tmp_min = info;
                } else {
                    minimizers.push(tmp_min);
                }
            } else {
                if info.0 <= tmp_min.0 {
                    minimizers.push(info);
                    tmp_min = info;
                }
            }
        }
        kmers.push(info);
    }

    println!("{:?}", kmers);
    println!("{:?}", minimizers);
}

pub fn kmer_to_uint(kmer: &str) {
    let mut res: i32 = 0;
    // let mut bs = [0u8; mem::size_of::<i32>()];

    for base in kmer.as_bytes() {
        res <<= 2;
        match base {
            65 => res += 0, // A
            67 => res += 1, // C
            71 => res += 2, // G
            84 => res += 3, // T
            _ => eprint!("Illegal base: {:?}", base),
        }
        println!("{:020b}", res);
    }
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
fn test_kmer() {
    let seq = "AAGTCGATCGAAGCTGATCGATCGATCGTGCTACGTGATGATGCTAGCCTGACTGATCGTAGCAGC";
    get_minimizers(&seq);
}