// Build pseudo index
use std::path::Path;
use std::collections::{HashMap, HashSet};

use bio::io::{fasta, fastq};

use debruijn::Kmer;
use debruijn::dna_string::*;
use debruijn::kmer::Kmer5;
use debruijn::Vmer;

use crate::utils;
use crate::kseq;


pub fn read_transcripts (file_name: &Path) -> HashMap::<u64, HashSet<String>> {
    let buf_reader = utils::read_fastx(file_name).unwrap();
    let reader = fasta::Reader::new(buf_reader);

    let mut kmer_index = HashMap::<u64, HashSet<String>>::default();

    let mut n = 1;
    for result in reader.records() {
        let record = result.unwrap();
        let read_id = String::from(record.id());
        let id_array: Vec<&str> = read_id.split('|').collect();
        let tscp_id = String::from(id_array[1]);

        if record.seq().len() <= 100 {
            continue;
        }

        let seq = DnaString::from_acgt_bytes(record.seq());
        let minimizers = kseq::split_reads(seq);

        for x in minimizers {
            kmer_index.entry(x.to_u64())
                .or_default()
                .insert(read_id.clone());
        }

        n += 1;
        if n % 10_000 == 0 {
            utils::info(format!("Loaded {} reads", n));
        }
    }
    utils::info(format!("Loaded {} reads", n));

    utils::info("Save index to file");
    let db_name = Path::new("/tmp/test.db");
    utils::save_index(&db_name, &kmer_index);
    
    kmer_index
}


// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_index() {
//         let mut test_idx = HashMap::<u64, HashSet<String>>::default();
//         let seq = DnaString::from_dna_string("ACAGCAGCAGCACGTATGACAGATAGTGACAGCAGTTTGTGACCGCAAGAGCAGTAATATGATG");
//         let minimizers = kseq::split_reads(seq);

//         for x in minimizers {
//             test_idx.entry(x.to_u64())
//                 .or_default()
//                 .insert(String::from(x.to_string()));
//         }
//         let seq2 = DnaString::from_dna_string("ACAGCAGCAGCACGTATGACAG");
//         let kmer: Kmer5 = seq2.get_kmer(0);
//         let tags: Vec<String> = test_idx.get(&kmer.to_u64()).unwrap().iter().cloned().collect();
//         assert_eq! (&tags[0], "ACAGC");
//     }
// }
