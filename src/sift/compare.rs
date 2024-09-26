// use crate::sift::rules;
// use crate::util::bam::*;
// use crate::util::dna::*;
// use crate::util::misc::make_intervals;

use anyhow;

// use rayon::prelude::*;
// use rust_htslib::bam::{self, Read};
// use std::cmp::min;
// use std::collections::HashSet;
// use std::sync::{Arc, Mutex};
// use std::{str, thread};

use crate::sift::sifter::*;

use super::CaseControlArgs as RunArgs;

/// Sift through BAM records to identify potential variant sites.  The
/// resulting sites may not be necessarily true hits, but they can
/// form a good starting candidate pool.
///
/// Later, routines in `aggregate` could revisit them all, regardless
/// of significance levels, and collect sufficient statistics for
/// further tests.
///
/// So, we just need to return a vector of (chr: Box<str>, lb: i64, ub: i64)
///
pub fn search_case_control(args: &RunArgs) -> anyhow::Result<()> {
    let block_size = args.block_size;
    let nthread_max = std::thread::available_parallelism()?.get();
    let nthread = match args.threads {
        Some(x) => nthread_max.min(x),
        None => nthread_max,
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(nthread as usize)
        .build_global()
        .unwrap();

    println!("Establishing BAM file sifters");

    let mut bam_fg = BamSifter::from_file(args.fg_bam.as_ref(), args.fg_bai.as_deref(), block_size);
    let mut bam_bg = BamSifter::from_file(args.bg_bam.as_ref(), args.bg_bai.as_deref(), block_size);

    println!("Searching for variable positions");

    bam_fg.sweep_variable_positions();
    bam_bg.sweep_variable_positions();

    // update variable positions to each other
    bam_bg.add_missing_positions(&bam_fg);
    bam_fg.add_missing_positions(&bam_bg);

    println!("Collecting sufficient statistics");

    // For each variable position

    // let  = bam_sifter_bg.get_forward_variable_positions();

    // Combine these positions

    // for (chr, blocks) in jobs {
    //     // Step 1. Make a list of variant sites: chr, lb, ub applying
    //     // a set of simple rules.
    //     let fg_var_positions = find_variable_positions(chr.as_ref(), &blocks, arc_bam_fg.clone());

    //     let bg_var_positions = find_variable_positions(chr.as_ref(), &blocks, arc_bam_bg.clone());

    //     // Step 2. Output BED format
    // }

    Ok(())
}
