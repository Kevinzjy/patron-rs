#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
RUST_BACKTRACE=1 ./target/release/patron -n 4 -r /home/zhangjy/data/database/genome/hg38/gencode.v36.transcripts.fa /home/zhangjy/data/01.ONT_methods/full_spectrum_Sequencing/FullSpeed/Normal_12h.fastq.gz
