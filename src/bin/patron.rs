use std::path::Path;
use log::{info, error};
use serde::Deserialize;

// extern crate bio;

use docopt::Docopt;
use failure::Error;

use patron::{build_index, utils, kseq};

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
        println! {"{} {}", PKG_NAME, PKG_VERSION};
        return Ok(());
    }

    // Init logger
    pretty_env_logger::init_timed();

    info!("Start running PATRON-RS");

    // Generate index from reference fasta
    let reference = Path::new(&args.flag_r);
    match reference.exists() {
        true => info! ("Loading fastq reads"),
        false => error! ("Can not find file: {}", &args.flag_r),
    }

    build_index::read_transcripts(&reference);

    Ok(())
}
