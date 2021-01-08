use debruijn::Kmer;
use std::num::Wrapping;

use debruijn::kmer::Kmer10;
use debruijn::kmer::Kmer15;


pub fn hash64(){
    let mask = 1 << 2 * 10 - 1;
    println! {"Mask: {}", mask};

    let kmer = Kmer10::from_u64(2283182);
    println! {"Kmer: {}", kmer.to_string()};

    let mut key = kmer.to_u64();
    println! {"Hash: {}", &key};

    println! {"Test: {}", (!key + (key << 21)) & mask};

    key = (!key + (key << 21)) & mask;
    println! {"Hash: {}", &key};
    key = key ^ key >> 24;
    println! {"Hash: {}", &key};
    key = ((key + (key << 3)) + (key << 8)) & mask;
    println! {"Hash: {}", &key};
    key = key ^ key >> 14;
    println! {"Hash: {}", &key};
    key = ((key + (key << 2)) + (key << 4)) & mask;
    println! {"Hash: {}", &key};
    key = key ^ key >> 28;
    println! {"Hash: {}", &key};
    key = (key + (key << 31)) & mask;
    println! {"Hash: {}", key};
}

// def hash64(key, mask):
//     """
//     https://github.com/lh3/minimap2/blob/c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd/sketch.c
//     """
//     key = ~key + (key << 21) & mask
//     key = key ^ key >> 24
//     key = ((key + (key << 3)) + (key << 8)) & mask
//     key = key ^ key >> 14
//     key = ((key + (key << 2)) + (key << 4)) & mask
//     key = key ^ key >> 28
//     key = (key + (key << 31)) & mask
//     return key


// def get_minimizer(pos, bases, w):
//     mask = (1 << 2*w) - 1
//     minimizer = bases[0:w]
//     min_hash = hash64(kmer_to_uint(minimizer), mask) << 8 | pos[w] - pos[0]
//     for i in range(1, len(bases) - w):
//         candidate = bases[i:i+w]
//         candidate_span = pos[i+w] - pos[i]
//         candidate_hash = hash64(kmer_to_uint(candidate), mask) << 8 | candidate_span

//         rev = revcomp(''.join(candidate))
//         rev_hash = hash64(kmer_to_uint(rev), mask) << 8 | candidate_span

//         if rev_hash < candidate_hash:
//             candidate_hash = rev_hash
//             candidate = rev

//         if candidate_hash < min_hash:
//             minimizer = candidate
//             min_hash = candidate_hash
//     return kmer_to_uint(minimizer)
