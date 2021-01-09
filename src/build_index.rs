// Build pseudo index
use std::path::Path;
// use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use log::{info, error};

use bio::io::{fasta, fastq};

use crate::utils;
use crate::kseq;

use hashbrown::HashMap;
use rustc_hash::FxHashMap;
use debruijn::Kmer;
use debruijn::dna_string::*;
use debruijn::kmer::Kmer10;
use debruijn::Vmer;

pub fn read_transcripts (file_name: &Path) {
    let buf_reader = utils::read_fastq(file_name).unwrap();
    let reader = fastq::Reader::new(buf_reader);

    // let mut kmer_index = HashMap::<u64, Vec<String>>::default();
    let mut kmer_index: FxHashMap<u64, Vec<String>> = FxHashMap::default();

    let mut n = 1;
    for result in reader.records() {
        let record = result.unwrap();
        let read_id = String::from(record.id());
        if record.seq().len() <= 100 {
            continue;
        }

        let seq = DnaString::from_acgt_bytes(&record.seq());
        let minimizers = kseq::split_reads(&seq);

        for x in minimizers {
            kmer_index.entry(x.to_u64())
                .or_default()
                .push(read_id.clone());
        }

        n += 1;
        if n >= 1_000_000 {
            break;
        }
    }
    info! ("Loaded {} reads", n);
    let test_seq = DnaString::from_dna_string("AACCTCTCCA");
    let test_kmer: Kmer10 = test_seq.get_kmer(0);
    let test_key = &test_kmer.to_u64();
    let read_ids = kmer_index.get(&test_key).unwrap();
    println! ("Kmer: {}, host reads: {:?}", test_kmer.to_string(), read_ids.len());
}
