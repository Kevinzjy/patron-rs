// Utilities
use std::fs::File;
use std::fmt::Display;
use std::io::{BufRead, BufWriter, BufReader};
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::DefaultHasher;
use std::process;

use bincode;
use chrono::Local;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::{ZlibDecoder, MultiGzDecoder};
use failure::Error;

// Reading Fastq
pub fn read_fastx (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    if file_name.extension().unwrap() == "gz" {
        _read_fastx_gz(file_name)
    } else {
        _read_fastx(file_name)
    }
}

fn _read_fastx_gz (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    let fastq_fn = File::open(file_name).unwrap();
    let gz = MultiGzDecoder::new(fastq_fn);
    let buf_reader = BufReader::with_capacity(32*1024, gz);
    Ok(Box::new(buf_reader))
}

fn _read_fastx (file_name: &Path) -> Result<Box<dyn BufRead>, Error> {
    let fastq_fn = File::open(file_name).unwrap();
    let buf_reader = BufReader::with_capacity(32*1024, fastq_fn);
    Ok(Box::new(buf_reader))
}

// Calculate Hash values
fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

// Serialize alignment index
pub fn save_index(file_name: &Path, data: &HashMap::<u64, HashSet<String>>) {
    let writer = BufWriter::new(File::create(file_name).unwrap());
    let mut encoder = ZlibEncoder::new(writer, Compression::default());
    bincode::serialize_into(&mut encoder, &data).unwrap();
}

pub fn load_index(file_name: &Path) -> HashMap::<u64, HashSet<String>> {
    let reader = BufReader::new(File::open(file_name).unwrap());
    let mut decoder = ZlibDecoder::new(reader);
    let decoded = bincode::deserialize_from(&mut decoder).unwrap();
    decoded
}

// A simple logging function
pub fn debug<T: Display> (info: T) {
    _logger("DEBUG", info);
}

pub fn info<T: Display> (info: T) {
    _logger("INFO", info);
}

pub fn warn<T: Display> (info: T) {
    _logger("WARN", info);
}

pub fn error<T: Display> (info: T) {
    _logger("ERROR", info);
    panic! ("Error occured!");
}

fn _logger<T: Display> (level: &str, info: T) {
    let date = Local::now();
    eprintln!("[{}] |{:5}| - {}", date.format("%a %Y-%m-%d %H:%M:%S"), level, info);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logger() {
        debug(String::from("Debugging"));
        info("Information");
        warn("Warning");
    }

    #[test]
    #[should_panic]
    fn test_error() {
        error("Test error");
    }
}
