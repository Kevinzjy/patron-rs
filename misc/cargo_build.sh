#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
./target/release/patron -n 4 -r /home/zhangjy/data/01.ONT_methods/full_spectrum_Sequencing/FullSpeed/Normal_12h.fastq.gz /home/zhangjy/data/01.ONT_methods/full_spectrum_Sequencing/FullSpeed/Normal_12h.fastq.gz
