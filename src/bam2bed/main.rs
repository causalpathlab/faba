use anyhow::{self, Context, Error, Result};

// use bio;
// use bio_types::annot::spliced::SeqSplicedStranded;
// use bio_types::sequence;

// use bio::alphabets::dna::revcomp;
use clap::{Args, Parser, Subcommand};

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read};

use std::path::Path;

use std::cmp::{max, min};

use std::sync::{Arc, Mutex};
use std::{str, thread};

// use env_logger;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct EpiArgs {
    /// foreground BAM file
    #[arg(short, long)]
    fg_bam: Box<str>,

    /// background BAM file
    #[arg(short, long)]
    bg_bam: Box<str>,
}

#[derive(Debug)]
struct DnaFreq {
    a: usize,   // number of A's
    t: usize,   // number of T's
    g: usize,   // number of G's
    c: usize,   // number of C's
    tot: usize, // total
    gpos: i64,  // genomic position
}

struct DnaFreqVecs {
    forward: Vec<DnaFreq>,
    reverse: Vec<DnaFreq>,
}

fn get_dna_freq(
    arc_bam: &Arc<Mutex<bam::IndexedReader>>,
    region: (&str, i64, i64),
) -> Result<DnaFreqVecs> {
    let (_, lb, ub) = region;

    let mut bam_reader = arc_bam.lock().expect("unable to lock the reader");

    bam_reader
        .fetch(region)
        .expect("unable to fetch the region");

    if lb >= ub {
        return Err(anyhow::anyhow!("lb >= ub"));
    }

    let nn = max(ub - lb, 0i64) as usize;
    let mut reverse_freq = Vec::with_capacity(nn);
    let mut forward_freq = Vec::with_capacity(nn);

    for g in lb..ub {
        forward_freq.push(DnaFreq {
            a: 0,
            t: 0,
            g: 0,
            c: 0,
            tot: 0,
            gpos: g,
        });
        reverse_freq.push(DnaFreq {
            a: 0,
            t: 0,
            g: 0,
            c: 0,
            tot: 0,
            gpos: g,
        });
    }

    // Iter aligned read and reference positions on a basepair level
    // https://docs.rs/rust-htslib/latest/src/rust_htslib/bam/ext.rs.html#135
    // [read_pos, genome_pos]

    for rr in bam_reader.rc_records() {
        match rr {
            Ok(rec) => {
                if rec.is_duplicate() {
                    continue;
                }

                let seq = rec.seq().as_bytes();

                for [rpos, gpos] in rec.aligned_pairs() {
                    let (r, g, v) = (rpos as usize, gpos as usize, gpos - lb);

                    if g < (lb as usize) || g >= (ub as usize) || v < 0 {
                        continue;
                    }

                    let bp = seq[r];

                    let freq = match rec.is_reverse() {
                        true => &mut reverse_freq[v as usize],
                        _ => &mut forward_freq[v as usize],
                    };

                    debug_assert_eq!(freq.gpos, gpos);
                    freq.tot += 1;
                    match bp {
                        b'A' | b'a' => freq.a += 1,
                        b'T' | b't' => freq.t += 1,
                        b'G' | b'g' => freq.g += 1,
                        b'C' | b'c' => freq.c += 1,
                        _ => (),
                    }
                }
            }
            _ => {
                // report error message?
            }
        }
    }

    // This doesn't take into account strand information and duplicates
    //
    // for pp in bam_reader.pileup() {
    //     let pileup = pp.expect("failed to instantiate pileup");
    //     let gpos = pileup.pos() as i64; // this is u32, so might not be ideal
    //     if gpos < lb || gpos >= ub {
    //         continue;
    //     }
    //     let v = gpos - lb;
    //     let freq = &mut ret[v as usize];
    //     debug_assert_eq!(freq.gpos, gpos);

    //     for aa in pileup.alignments() {
    //         if aa.is_del() || aa.is_refskip() {
    //             continue;
    //         }

    //         let r = aa.qpos().unwrap();
    //         let bp = aa.record().seq()[r];

    //         // if aa.record().is_reverse() {
    //         //     dbg!(aa.record());
    //         // }

    //         freq.tot += 1;
    //         match bp {
    //             b'A' | b'a' => freq.a += 1,
    //             b'T' | b't' => freq.t += 1,
    //             b'G' | b'g' => freq.g += 1,
    //             b'C' | b'c' => freq.c += 1,
    //             _ => (),
    //         }
    //     }
    // }

    Ok(DnaFreqVecs {
        forward: forward_freq,
        reverse: reverse_freq,
    })
}

fn main() -> anyhow::Result<()> {
    /////////////////////////
    // set up event logger //
    /////////////////////////
    // env_logger::Builder::new()
    //     .target(env_logger::Target::Stderr)
    //     .init();

    let args = EpiArgs::parse();

    let bam_file_bg = args.bg_bam.as_ref();
    let _ = check_bam_index(bam_file_bg, None);

    let bam_file_fg = args.fg_bam.as_ref();
    let _ = check_bam_index(bam_file_fg, None);

    dbg!(&bam_file_fg);
    dbg!(&bam_file_bg);

    // shared index reader
    let arc_bam_bg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_bg)?));

    let arc_bam_fg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_fg)?));

    // need to figure out chromosome names and boundaries
    // let br = bam::Reader::from_path(bam_file_bg)?;
    // let hdr = br.header();
    // let mut chr2tid: HashMap<Box<str>, usize> = HashMap::new();
    // for (tid, tgt) in hdr.target_names().iter().enumerate() {
    //     let chr_name = str::from_utf8(tgt).unwrap_or(".");
    //     chr2tid.insert(chr_name.into(), tid);
    // }

    // chr18:34220983-34318581
    let chr_name = "chr18";
    let (lb, ub) = (34304689 as i64, 34304694 as i64);

    // let mut br = bam::IndexedReader::from_path(bam_file_name)?;
    // br.fetch((chr_name, lb, ub))?;
    // let _ = get_dna_freq(&mut br, lb, ub);

    // thread::spawn(move || {
    let region = (chr_name, lb, ub);
    let count_bg = get_dna_freq(&arc_bam_bg, region).unwrap();
    let count_fg = get_dna_freq(&arc_bam_fg, region).unwrap();

    //
    for (r, g) in (lb..ub).enumerate() {
        let bg = &count_bg.forward[r];
        let fg = &count_fg.forward[r];
        debug_assert_eq!(g, bg.gpos);
        debug_assert_eq!(g, fg.gpos);

        dbg!(bg);
        dbg!(fg);
    }
    // });

    Ok(())
}

fn check_bam_index(bam_file_name: &str, idx_file_name: Option<&str>) -> anyhow::Result<Box<str>> {
    // log::info!("Checking BAM index");

    let idx_file = match idx_file_name {
        Some(x) => String::from(x),
        None => format!("{}.bai", bam_file_name),
    };

    if Path::new(&idx_file).exists() {
        return Ok(idx_file.into_boxed_str());
    }

    let ncore = thread::available_parallelism()
        .expect("failed to figure out number of cores")
        .get();

    // log::info!(
    //     "Creating a new index file {} using {} cores",
    //     &idx_file,
    //     &ncore
    // );

    // need to build an index for this bam file
    bam::index::build(
        bam_file_name,
        Some(&idx_file),
        bam::index::Type::Bai,
        ncore as u32,
    )?;

    Ok(idx_file.into_boxed_str())
}

// FAVA
//
// Filtering
// Altered
//
//
// parser from
// BAM file
// Expression
// A
// N
//

// // let mut jobs = vec![];
// // let block_size = 1024 as usize;

// fn make_blocks(max_size: i64, block_size: i64) -> Vec<(i64, i64)> {
//     let mut jobs = vec![];
//     for lb in (0..max_size).step_by(block_size as usize) {
//         let ub = min(max_size, lb + block_size);
//         jobs.push((lb, ub));
//     }
//     return jobs;
// }

// let mut seq_blocks = HashMap::<i32, Vec<(i64, i64)>>::new();
// let block_size = 10000 as i64;

// let hdr = br.header();

// for (tid, k) in hdr.target_names().iter().enumerate() {
//     let name = str::from_utf8(k).unwrap_or(".");
//     let max_size = hdr.target_len(tid as u32).unwrap() as i64;
//     seq_blocks.insert(tid as i32, make_blocks(max_size, block_size));
// }

// fn count_c2u(rec: &bam::Record) -> usize {
//     for cigar in rec.cigar().iter() {
//         //
//     }

//     0
// }

// for (tid, v) in seq_blocks.iter() {
//     dbg!(tid);
//     dbg!(v.len());

//     let mut bam_idx =
//         bam::IndexedReader::from_path(bam_file_name).expect("failed to open BAM IndexedReader");

//     for (lb, ub) in v {
//         bam_idx
//             .fetch((*tid, *lb, *ub))
//             .expect("failed to fetch BAM");

//         for rr in bam_idx.records() {
//             let rec = rr.expect("failed to fetch the alignment");

//             assert_eq!(rec.tid(), *tid);

//             // extract 10x cell barcode
//             // if let Ok(x) = rec.aux(b"CB") {
//             //     dbg!(x);
//             // }

//             // extract 10x UMI barcode
//             // if let Ok(x) = rec.aux(b"UB") {
//             //     dbg!(x);
//             // }

//             // let seq = rec.seq().as_bytes();
//             // dbg!(str::from_utf8(&seq));
//         }
//     }
// }

// some helper portions
// let pos = rec.pos(); // Position of the read

// rec.is_reverse();
// rec.seq_len();
// rec.seq();

// // for c in seq {
// //     let cc = c as char;
// // }

// // rec.reference_start();

// // rec.base(0);

// let seq = &rec.seq().as_bytes();
// let seq_chars = str::from_utf8(&seq)?.chars();

// for (i, c) in seq_chars.enumerate() {
//     //
// }

// use bio::alphabets::rna;

// let rev_seq = rna::revcomp(seq);

// dbg!(rev_seq);

// for b in 0..rec.len() {
//     //
//     let c = rec.base(b) as char;
//     //
// }

// // for (i, c) in y.chars().enumerate() {
// //     //
// // }

// // for x in rec.seq().i
// //     //
// // }

// // bam::Record::is_reverse(&self)

// // if let Ok(mods) = rec.basemods_iter() {
// //     //
// // }

// match rec.strand() {
//     _ => {}
// }

// fg.fetch((tid, regions.0, regions.1))
//     .expect("failed to fetch the region of interest");

//
// for rr_fg in fg.records() {
//     let rec = rr_fg?;
//     // rec_fg.
// }

// fn proc(args: &EpiArgs, tid: u32, regions: (i64, i64)) -> Result<()> {
//     use bam::ext::BamRecordExtensions;
//     use bio_types::sequence::SequenceRead;
//     //
//     let mut fg = bam::IndexedReader::from_path(args.fg_bam.as_ref())?;
//     let mut bg = bam::IndexedReader::from_path(args.bg_bam.as_ref())?;

//     bg.fetch((tid, regions.0, regions.1))
//         .expect("failed to fetch the region of interest");

//     fg.fetch((tid, regions.0, regions.1))
//         .expect("failed to fetch the region of interest");

//     let ref_base = b'C';
//     let alt_base: u8 = b'U';

//     // 1. construct background frequency map
//     let mut ref_count = HashMap::<i64, usize>::new();
//     let mut alt_count = HashMap::<i64, usize>::new();

//     // build a count matrix to keep track of mutations
//     // C -> U

//     for rr in bg.rc_records() {
//         let rec = rr?;
//         if rec.is_duplicate() {
//             continue;
//         };

//         // https://docs.rs/rust-htslib/latest/src/rust_htslib/bam/ext.rs.html#342
//         // fn aligned_pairs(&self) -> IterAlignedPairs {
//         //     IterAlignedPairs {
//         //         genome_pos: self.pos(),
//         //         read_pos: 0,
//         //         cigar: self.cigar().take().0,
//         //         remaining_match_bp: 0,
//         //         cigar_index: 0,
//         //     }
//         // }
//         //

//         let seq = rec.seq().as_bytes();

//         // Iter aligned read and reference positions on a basepair level
//         // https://docs.rs/rust-htslib/latest/src/rust_htslib/bam/ext.rs.html#135
//         // [read_pos, genome_pos]
//         for [read_pos, genome_pos] in rec.aligned_pairs() {
//             let (r, g) = (read_pos as usize, genome_pos as usize);
//             assert_eq!(rec.base(r), seq[r]);
//         }

//         // Iter over <([read_start, read_stop], [genome_start, genome_stop])
//         // blocks of continously aligned reads.
//         // for ([seq_lb, seq_ub], [gen_lb, gen_ub]) in rec.aligned_block_pairs() {
//         //     //
//         // }
//     }

//     Ok(())
// }
