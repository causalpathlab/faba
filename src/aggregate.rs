use super::AggregateArgs;
use anyhow;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read};
use std::cmp::{max, min};
use std::sync::{Arc, Mutex};
use std::{str, thread};

use crate::util::check_bam_index;

pub fn run_aggregate(args: &AggregateArgs) -> anyhow::Result<()> {
    // Visit all the alignments and figure out

    let nthread_max = thread::available_parallelism()
        .expect("failed to figure out number of cores")
        .get();

    let nthread = match args.threads {
        Some(x) => min(nthread_max, x),
        None => nthread_max,
    };

    let (bam_file_bg, bam_file_fg) = (args.bg_bam.as_ref(), args.fg_bam.as_ref());

    check_bam_index(bam_file_bg, args.bg_bai.as_deref())
        .expect("check index for the background BAM");
    check_bam_index(bam_file_fg, args.fg_bai.as_deref())
        .expect("check index for the foreground BAM");

    let block_size = match args.bsize {
        Some(bs) => bs,
        _ => 10_000,
    } as i64;

    let mut jobs = vec![];

    let br = bam::Reader::from_path(bam_file_fg)?;
    let hdr = br.header();

    for (tid, name) in hdr.target_names().iter().enumerate() {
        let max_size = hdr.target_len(tid as u32).unwrap() as i64;
        let chr_name = Box::new(str::from_utf8(name).unwrap());
        jobs.push((chr_name, make_blocks(max_size, block_size)));
    }

    // shared index reader
    let arc_bam_bg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_bg)?));

    let arc_bam_fg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_fg)?));

    rayon::ThreadPoolBuilder::new()
        .num_threads(nthread as usize)
        .build_global()
        .unwrap();

    for (chr, blocks) in jobs.iter() {
        //
        let chr_name = *(chr.as_ref());

        dbg!(chr_name);

        let _ = blocks.par_iter().map(|(lb, ub)| {
            let region = (chr_name, *lb, *ub);
            let bg = get_dna_freq(&arc_bam_bg, region);
            let fg = get_dna_freq(&arc_bam_fg, region);
            match (bg, fg) {
                (Ok(freq_bg), Ok(freq_fg)) => {
                    freq_bg.forward;
                    freq_bg.reverse;
                    freq_fg.forward;
                    freq_fg.reverse;
                    //
                }
                _ => {
                    // do nothing --> ignore errors
                }
            }
        });

        // blocks.par_bridge();
    }

    Ok(())
}

fn make_blocks(max_size: i64, block_size: i64) -> Vec<(i64, i64)> {
    let mut jobs = vec![];
    for lb in (0..max_size).step_by(block_size as usize) {
        let ub = min(max_size, lb + block_size);
        jobs.push((lb, ub));
    }
    return jobs;
}

///////////////////////////
// DNA frequency vectors //
///////////////////////////

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

////////////////////////////////////////////
// Extract DNA base pair frequency tables //
////////////////////////////////////////////

fn get_dna_freq(
    arc_bam: &Arc<Mutex<bam::IndexedReader>>,
    region: (&str, i64, i64),
) -> anyhow::Result<DnaFreqVecs> {
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

                // TODO: cell barcode umi
                // extract 10x cell barcode
                // if let Ok(cb) = rec.aux(b"CB") {
                //     dbg!(x);
                // }

                // extract 10x UMI barcode
                // if let Ok(umi) = rec.aux(b"UB") {
                //     dbg!(x);
                // }

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

    Ok(DnaFreqVecs {
        forward: forward_freq,
        reverse: reverse_freq,
    })
}
