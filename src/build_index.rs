// Build pseudo index
use std::path::Path;
use log::{info, error};

use bio::io::{fasta, fastq};

use crate::utils;
use crate::kseq;

use debruijn::Kmer;
use debruijn::dna_string::*;
use debruijn::kmer::Kmer15;
use debruijn::Vmer;


pub fn read_transcripts (file_name: &Path) {
    let buf_reader = utils::read_fastq(file_name).unwrap();
    let reader = fastq::Reader::new(buf_reader);

    let mut n = 1;
    for result in reader.records() {
        let record = result.unwrap();
        // println! {"Record {}", n};
        // println! {"{} {}", record.id(), record.seq().len()};

        let seq = DnaString::from_acgt_bytes(record.seq());
        let slice1 = seq.slice(10, 40);

        let first_kmer: Kmer15 = slice1.get_kmer(0);
        // println! {"{:?}", first_kmer};

        n += 1;
        if n >= 10 {
            break;
        }
    }

    kseq::hash64();
}