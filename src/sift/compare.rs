use crate::sift::rules;
use crate::util::bam::*;
use crate::util::dna::*;
use crate::util::misc::make_intervals;

use anyhow;

use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::cmp::min;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};
use std::{str, thread};

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
pub fn search(args: &RunArgs) -> anyhow::Result<()> {
    // Visit all the alignments and figure out

    let nthread_max = thread::available_parallelism()?.get();

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

    // Step 1. Find variable base positions
    let find_variable_positions = |chr_name: &str,
                                   blocks: &Vec<(i64, i64)>,
                                   arc_bam: Arc<Mutex<bam::IndexedReader>>|
     -> HashSet<i64> {
        let arc_var_pos = Arc::new(Mutex::new(HashSet::new()));

        blocks.iter().par_bridge().for_each(|(lb, ub)| {
            let region = (chr_name, *lb, *ub);
            let base_filter = rules::BaseFilters::new();
            let mut temp = vec![];

            if let Ok(freq_map) = get_dna_base_freq(&arc_bam, region) {
                for samp in freq_map.samples() {
                    if let Some(stats) = freq_map.get_forward(samp) {
                        for bs in stats {
                            if base_filter.is_variable(bs) {
                                temp.push(bs.position());
                            }
                        }
                    }
                }
            }

            let mut ret = arc_var_pos.lock().unwrap();
            for j in temp {
                ret.insert(j);
            }
        });

        let ret = arc_var_pos.lock().expect("").clone();
        ret
    };

    for (chr, blocks) in jobs {
        // Step 1. Make a list of variant sites: chr, lb, ub applying
        // a set of simple rules.
        let fg_var_positions = find_variable_positions(chr.as_ref(), &blocks, arc_bam_fg.clone());

        let bg_var_positions = find_variable_positions(chr.as_ref(), &blocks, arc_bam_bg.clone());

        // Step 2. Output BED format
    }

    Ok(())
}
