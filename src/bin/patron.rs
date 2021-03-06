use std::path::Path;
use std::collections::HashMap;
use serde::Deserialize;

// extern crate bio;

use docopt::Docopt;
use failure::Error;

use patron::{build_index, align, utils};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
Pseudo-Alignment of Transcriptome Reads from Oxford Nanopore Sequencing

Usage:
  patron-rs [--num-threads=<n>] -r FASTA <reads>
  patron-rs --help | --version

Options:
    -n --num-threads N  Number of worker theads [default: 2]
    -r FASTA            Reference fasta
    -h --help           Show this screen
    -v --version        Show version
";

#[derive(Clone, Debug, Deserialize)]
struct Args {
    arg_reads: String,

    flag_r: String,
    flag_num_threads: usize,

    flag_version: bool,
    flag_help: bool,
}


fn main() -> Result<(), Error> {
    // Parse parameters
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    if args.flag_version {
        print! ("{} {}", PKG_NAME, PKG_VERSION);
        return Ok(());
    }

    // Init logger
    utils::info("Start running PATRON");

    // Generate index from reference fasta
    let pseudo_index;
    let idx_file = Path::new("/tmp/test.db");
    if idx_file.exists() {
        utils::info("Loading pre-built index");
        pseudo_index = utils::load_index(idx_file);
    } else {
        let reference = Path::new(&args.flag_r);
        match reference.exists() {
            true => utils::info("Building index"),
            false => utils::error(format!("Can not find file: {}", &args.flag_r)),
        }
        pseudo_index = build_index::read_transcripts(&reference);
    }
    
    // Align reads
    let reads = Path::new(&args.arg_reads);
    match reads.exists() {
        true => utils::info("Mapping reads"),
        false => utils::error(format!("Can not find file: {}", &args.arg_reads)),
    }
    align::align_reads(&reads, pseudo_index);

    // GC
    utils::info("All Finished!");
    Ok(())
}
