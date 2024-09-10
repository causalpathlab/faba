use anyhow;
use clap::Args;

#[derive(Args)]
pub struct CaseControlArgs {
    /// foreground BAM file
    #[arg(short, long)]
    fg_bam: Box<str>,

    /// background BAM file
    #[arg(short, long)]
    bg_bam: Box<str>,

    /// foreground BAI file (default: <FG_BAM>.bai)
    #[arg(long)]
    fg_bai: Option<Box<str>>,

    /// background BAI file (default: <BG_BAM>.bai)
    #[arg(long)]
    bg_bai: Option<Box<str>>,

    /// number of threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// block size (default: 10000)
    #[arg(long)]
    bsize: Option<usize>,

    /// output file header
    #[arg(short, long)]
    output: Option<Box<str>>,
}

use crate::dnafreq::*;
use crate::util::check_bam_index;
use fastapprox::faster as fa;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::cmp::min;
use std::sync::{Arc, Mutex};
use std::{str, thread};

pub fn run_case_control(args: &CaseControlArgs) -> anyhow::Result<()> {
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
        jobs.push((chr_name, make_intervals(max_size, block_size)));
    }

    // shared index reader
    let arc_bam_bg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_bg)?));

    let arc_bam_fg = Arc::new(Mutex::new(bam::IndexedReader::from_path(bam_file_fg)?));

    rayon::ThreadPoolBuilder::new()
        .num_threads(nthread as usize)
        .build_global()
        .unwrap();

    // TODO: scan to figure out potential variant sites
    // MAF in one side != MAF in the other side

    for (chr, blocks) in jobs {
        let chr_name = *(chr.as_ref());

        blocks.iter().par_bridge().for_each(|(lb, ub)| {
            let region = (chr_name, *lb, *ub);

            if let (Ok(fg), Ok(bg)) = (
                get_dna_freq(&arc_bam_fg, region),
                get_dna_freq(&arc_bam_bg, region),
            ) {
                for s in fg.samples() {
                    if !bg.has_sample(s) {
                        continue;
                    }

                    if let (Some(stat1), Some(stat0)) = (fg.get_forward(s), bg.get_forward(s)) {
                        debug_assert_eq!(stat1.len(), stat0.len());
                    }

                    if let (Some(stat1), Some(stat0)) = (fg.get_reverse(s), bg.get_reverse(s)) {
                        debug_assert_eq!(stat1.len(), stat0.len());

                        let nn = stat1.len();
                        for j in 0..nn {
                            let fg_bp = &stat1[j];
                            let bg_bp = &stat0[j];

                            if let (Some(maf_fg), Some(maf_bg)) = (
                                fg_bp.major_allele_frequency(),
                                bg_bp.major_allele_frequency(),
                            ) {
                                if (maf_fg.0 == maf_bg.0)
                                    && (maf_fg.1 > 0.9f32 && maf_bg.1 > 0.9f32)
                                {
                                    //
                                } else {
				    dbg!(fg_bp);
				    dbg!(bg_bp);
				}
                            }
                        }
                    }

                    //
                }
                // For each sample, we collect statistics
                // for s in fg.forward.keys() {
                //     if let (Some(stat1), Some(stat0)) = (fg.get_forward(s), bg.get_forward(s)) {
                //         debug_assert_eq!(stat1.len(), stat0.len());
                //         let nn = stat1.len();
                //         for j in 0..nn {
                //             let fg_bp = &stat1[j];
                //             let bg_bp = &stat0[j];

                //             if fg_bp.tot + bg_bp.tot < 2f32 {
                //                 continue;
                //             }

                //             // todo: major allele frequency difference

                //             let score = local_bayes_factor(FreqStat {
                //                 fg: fg_bp,
                //                 bg: bg_bp,
                //             });

                //             if score > 3f32 {
                //                 dbg!(s);
                //                 dbg!(score);
                //                 dbg!(&stat1[j]);
                //                 dbg!(&stat0[j]);
                //             }
                //         }
                //     }
                // }

                // // ignore empty block
                // for j in 0..bg.forward.len() {
                //     if bg.forward[j].tot > 0 && fg.forward[j].tot > 0 {
                //         //
                //     }
                // }
            }
        });
        dbg!((&chr_name, blocks.len()));
    }

    Ok(())
}

/// Other utilities
/// make a vector of intervals
fn make_intervals(max_size: i64, block_size: i64) -> Vec<(i64, i64)> {
    let mut jobs = vec![];
    for lb in (0..max_size).step_by(block_size as usize) {
        let ub = min(max_size, lb + block_size);
        jobs.push((lb, ub));
    }
    return jobs;
}

/// a wrapper for input frequency stat argument
pub struct FreqStat<'a> {
    fg: &'a DnaFreq,
    bg: &'a DnaFreq,
}

/// statistical assessment of variant calling
/// by computing vanilla Bayes factor
pub fn local_bayes_factor(stat: FreqStat) -> f32 {
    let fg = stat.fg;
    let bg = stat.bg;

    let a0 = 0.25f32;

    let lgamma_ratio = |a: f32, b: f32, pc: f32| -> f32 {
        fa::ln_gamma(a + pc) + fa::ln_gamma(b + pc) - fa::ln_gamma(a + b + pc + pc)
    };

    // lgamma_ratio(fg.a, bg.a, a0)
    //     + lgamma_ratio(fg.t, bg.t, a0)
    //     + lgamma_ratio(fg.g, bg.g, a0)
    //     + lgamma_ratio(fg.c, bg.c, a0)
    // - lgamma_ratio(fg.tot, bg.tot, a0 * 4.)
    0f32
}
