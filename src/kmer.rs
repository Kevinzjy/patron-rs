use std::str;
use std::mem;
use log::{info, error};

use byteorder::{ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};

enum Base {
    A,
    T,
    C,
    G,
}

pub struct Khmer {
    k: f64,
    seq: String,
}

// pub fn convert_actg(base: String) {
//     match base ==
// }

pub fn split_reads() {
    let seq = String::from("get_kmer");
    let mut kmers: Vec<String> = Vec::new();
    for offset in 0..seq.len() - 15 + 1{
        let kmer = String::from(&seq);
    }
}

pub fn kmer_to_uint(kmer: String) {
    let mut res: i32 = 0;
    // let mut bs = [0u8; mem::size_of::<i32>()];

    for base in kmer.as_bytes() {
        res <<= 2;
        match base {
            65 => res += 0, // A
            67 => res += 1, // C
            71 => res += 2, // G
            84 => res += 3, // T
            _ => error! ("Illegal base: {:?}", base),
        }
        println!("{:020b}", res);
    }
}

#[test]
fn test_kmer() {
    let seq = String::from("AGTCGATCGA");
    kmer_to_uint(seq);
}