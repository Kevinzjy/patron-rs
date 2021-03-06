// Build pseudo index
use std::path::Path;
use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::cmp::Reverse;
use std::io::{Write, BufWriter};

use bio::io::fastq;
use debruijn::dna_string::*;
use debruijn::Kmer;

use crate::utils;
use crate::kseq;


pub fn align_reads(file_name: &Path, index: HashMap::<u64, HashSet<String>>) {
    let buf_reader = utils::read_fastx(file_name).unwrap();
    let reader = fastq::Reader::new(buf_reader);

    let mut output = BufWriter::new(File::create("/tmp/pseudo_alignment.txt").expect("Unable to create file"));
    let mut n_aligned = 0;
    for result in reader.records() {
        let record = result.unwrap();
        let read_id = String::from(record.id());

        let seq = DnaString::from_acgt_bytes(record.seq());
        let minimizers = kseq::split_reads(seq);

        let mut aligned_contigs = HashMap::<&String, u64>::new();
        let mut n = 0;
        for x in minimizers {
            n += 1;
            let key = x.to_u64();
            let contig = match index.get(&key) {
                Some(val) => val,
                None => continue,
            };

            for i in contig {
                *aligned_contigs.entry(i).or_insert(0) += 1;
            }

        }
        let mut count_vec: Vec<_> = aligned_contigs.iter().collect();
        &count_vec.sort_by_key(|k| Reverse(k.1));

        let mut j = 0;
        for gene_id in count_vec.iter() {
            write!(output, "{}\t{}\t{}\t{}\n", read_id, gene_id.0, gene_id.1, n).expect("Write Error");
            j += 1;
            if j >= 10 {
                break;
            }
        }

        n_aligned += 1;
        if n_aligned % 10_000 == 0 {
            utils::info(format!("Aligned {} reads", n_aligned));
            break;
        }
    }
}
