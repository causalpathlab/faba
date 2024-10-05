use crate::util::misc::make_intervals;

use anyhow;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

use crate::sift::*;

struct DirectedPositions {
    forward_positions: Vec<i64>,
    reverse_positions: Vec<i64>,
}

struct DirectedStats {
    forward: HashMap<(BamSample, Box<str>), Vec<DnaBaseStat>>,
    reverse: HashMap<(BamSample, Box<str>), Vec<DnaBaseStat>>,
}

struct DirectedSets {
    forward_positions: HashSet<i64>,
    reverse_positions: HashSet<i64>,
}

pub struct BamSifter {
    bam_reader: bam::IndexedReader,
    jobs: Vec<(Box<str>, Vec<(i64, i64)>)>,
    forward_variable_map: HashMap<Box<str>, HashSet<i64>>,
    reverse_variable_map: HashMap<Box<str>, HashSet<i64>>,
    forward_stat: HashMap<(BamSample, Box<str>), Vec<DnaBaseStat>>,
    reverse_stat: HashMap<(BamSample, Box<str>), Vec<DnaBaseStat>>,
}

#[allow(dead_code)]
impl BamSifter {
    ///
    /// create a wrapper for BAM file sifting routines
    ///
    /// * `bam_file` - alignment file name
    /// * `bai_file` - index file name
    ///
    pub fn from_file(bam_file: &str, bai_file: Option<&str>, block_size: Option<usize>) -> Self {
        //
        let block_size = match block_size {
            Some(x) => x as i64,
            _ => 10_000i64,
        };

        //
        // read header information
        //
        let br = bam::Reader::from_path(bam_file)
            .expect(&format!("failed to initialize BAM file: {}", bam_file));

        let hdr = br.header();

        let mut chr_interval_jobs = vec![];

        for (tid, name) in hdr.target_names().iter().enumerate() {
            let max_size = hdr.target_len(tid as u32).unwrap() as i64;
            let name_ = String::from_utf8(name.to_vec()).unwrap();
            let chr_name = name_.into_boxed_str();
            chr_interval_jobs.push((chr_name, make_intervals(max_size, block_size)));
        }

        // check index for a fast look-up
        let index_file = check_bam_index(bam_file, bai_file)
            .expect(&format!("failed to generate index for: {}", bam_file));

        BamSifter {
            bam_reader: bam::IndexedReader::from_path_and_index(bam_file, &index_file)
                .expect("failed to create indexed reader"),
            jobs: chr_interval_jobs,
            forward_variable_map: HashMap::new(),
            reverse_variable_map: HashMap::new(),
            forward_stat: HashMap::new(),
            reverse_stat: HashMap::new(),
        }
    }

    /// Sweep all the blocks to identify variable positions. This will
    /// fill in the found variable positions in forward_variable_map
    /// and reverse_variable_map.
    ///
    pub fn sweep_variable_positions(&mut self) -> anyhow::Result<()> {
        for (chr, blocks) in self.jobs.iter() {
            let fvar_set = self
                .forward_variable_map
                .entry(chr.clone())
                .or_insert(HashSet::new());

            let rvar_set = self
                .reverse_variable_map
                .entry(chr.clone())
                .or_insert(HashSet::new());

            let forward_arc = Arc::new(Mutex::new(fvar_set));
            let reverse_arc = Arc::new(Mutex::new(rvar_set));

            let bam_arc = Arc::new(Mutex::new(&mut self.bam_reader));

            blocks
                .iter()
                .par_bridge()
                .try_for_each(|(lb, ub)| -> anyhow::Result<()> {
                    let region = (chr.as_ref(), *lb, *ub);
                    let base_filter = rules::BaseFilters::new();
                    let mut forward = vec![];
                    let mut reverse = vec![];

                    if let Ok(freq_map) = get_dna_base_freq(&bam_arc, region) {
                        for samp in freq_map.samples() {
                            // forward direction : 5 -> 3
                            if let Some(stats) = freq_map.get_forward(samp) {
                                for bs in stats {
                                    if base_filter.is_variable(bs) {
                                        forward.push(bs.position());
                                    }
                                }
                            }
                            // reverse direction : 3 -> 5
                            if let Some(stats) = freq_map.get_reverse(samp) {
                                for bs in stats {
                                    if base_filter.is_variable(bs) {
                                        reverse.push(bs.position());
                                    }
                                }
                            }
                        }
                    }

                    forward_arc
                        .lock()
                        .expect("failed to lock forward")
                        .extend(forward);
                    reverse_arc
                        .lock()
                        .expect("failed to lock reverse")
                        .extend(reverse);

                    Ok(())
                })?;
        }
        Ok(())
    }

    /// add (potentially) missed variable positions
    pub fn add_missed_positions(&mut self, other: &BamSifter) {
        self.add_forward_positions(other.get_forward_variable_positions());
        self.add_reverse_positions(other.get_reverse_variable_positions());
    }

    fn add_forward_positions(&mut self, new_pos_to_add: &HashMap<Box<str>, HashSet<i64>>) {
        for (k, new_pos_set) in new_pos_to_add {
            if let Some(pos_set) = self.forward_variable_map.get_mut(k) {
                pos_set.extend(new_pos_set);
            }
        }
    }

    fn add_reverse_positions(&mut self, new_pos_to_add: &HashMap<Box<str>, HashSet<i64>>) {
        for (k, new_pos_set) in new_pos_to_add {
            if let Some(pos_set) = self.reverse_variable_map.get_mut(k) {
                pos_set.extend(new_pos_set);
            }
        }
    }

    pub fn get_forward_variable_positions(&self) -> &HashMap<Box<str>, HashSet<i64>> {
        &self.forward_variable_map
    }

    pub fn get_reverse_variable_positions(&self) -> &HashMap<Box<str>, HashSet<i64>> {
        &self.reverse_variable_map
    }

    /// Populate statistics. This will accumulate sufficient
    /// statistics of the variable positions previously found by
    /// [`sweep_variable_positions`].
    ///
    pub fn populate_statistics(&mut self) {
        let fstat_arc = Arc::new(Mutex::new(&mut self.forward_stat));
        let rstat_arc = Arc::new(Mutex::new(&mut self.reverse_stat));
        let bam_arc = Arc::new(Mutex::new(&mut self.bam_reader));

        for (chr, positions) in self.forward_variable_map.iter() {
            positions.iter().par_bridge().for_each(|x| {
                let _chr = chr.as_ref();
                let _bp = *x;
                let region = (_chr, _bp, _bp + 1);

                let mut fstat = fstat_arc.lock().expect("unable to lock fstat");
                let mut rstat = rstat_arc.lock().expect("unable to lock rstat");

                if let Ok(freq_map) = get_dna_base_freq(&bam_arc, region) {
                    for samp in freq_map.samples() {
                        let fstat_vec = fstat.entry((samp.clone(), chr.clone())).or_insert(vec![]);

                        if let Some(statvec) = freq_map.get_forward(samp) {
                            for bs in statvec {
                                fstat_vec.push(bs.clone());
                            }
                        }

                        let rstat_vec = rstat.entry((samp.clone(), chr.clone())).or_insert(vec![]);

                        if let Some(statvec) = freq_map.get_reverse(samp) {
                            for bs in statvec {
                                rstat_vec.push(bs.clone());
                            }
                        }
                    }
                }
            });
        }
    }
}
