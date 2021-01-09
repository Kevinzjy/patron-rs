// Utilities
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use log::{info, error};


use flate2::read::MultiGzDecoder;
use failure::Error;


pub fn read_fastq (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    if file_name.extension().unwrap() == "gz" {
        _read_fastq_gz(file_name)
    } else {
        _read_fastq(file_name)
    }
}

fn _read_fastq_gz (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    let fastq_fn = File::open(file_name).unwrap();
    let gz = MultiGzDecoder::new(fastq_fn);
    let buf_reader = BufReader::with_capacity(32*1024, gz);
    Ok(Box::new(buf_reader))
}


fn _read_fastq (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    let fastq_fn = File::open(file_name).unwrap();
    let buf_reader = BufReader::with_capacity(32*1024, fastq_fn);
    Ok(Box::new(buf_reader))
}


fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
