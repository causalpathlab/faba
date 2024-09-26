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

struct DirectedSets {
    forward_positions: HashSet<i64>,
    reverse_positions: HashSet<i64>,
}

pub struct BamSifter {
    arc_bam: Arc<Mutex<bam::IndexedReader>>,
    jobs: Vec<(Box<str>, Vec<(i64, i64)>)>,
    forward_variable_map: HashMap<Box<str>, HashSet<i64>>,
    reverse_variable_map: HashMap<Box<str>, HashSet<i64>>,
    forward_stat_map: HashMap<Box<str>, HashMap<BamSample, Vec<DnaBaseStat>>>,
    reverse_stat_map: HashMap<Box<str>, HashMap<BamSample, Vec<DnaBaseStat>>>,
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
            arc_bam: Arc::new(Mutex::new(
                bam::IndexedReader::from_path_and_index(bam_file, &index_file)
                    .expect("failed to create indexed reader"),
            )),
            jobs: chr_interval_jobs,
            forward_variable_map: HashMap::new(),
            reverse_variable_map: HashMap::new(),
            forward_stat_map: HashMap::new(),
            reverse_stat_map: HashMap::new(),
        }
    }

    pub fn sweep_variable_positions(&mut self) {
        for (chr, blocks) in self.jobs.iter() {
            if let Ok(var_positions) = self.find_variable_positions(chr, blocks) {
                self.forward_variable_map
                    .insert(chr.clone(), var_positions.forward_positions);
                self.reverse_variable_map
                    .insert(chr.clone(), var_positions.reverse_positions);
            }
        }
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

    pub fn add_missing_positions(&mut self, other: &BamSifter) {
        self.add_forward_positions(other.get_forward_variable_positions());
        self.add_reverse_positions(other.get_reverse_variable_positions());
    }

    pub fn get_forward_variable_positions(&self) -> &HashMap<Box<str>, HashSet<i64>> {
        &self.forward_variable_map
    }

    pub fn get_reverse_variable_positions(&self) -> &HashMap<Box<str>, HashSet<i64>> {
        &self.reverse_variable_map
    }

    /// Sift bases within the given region by checking `is_variable`
    ///
    /// * `chr_name` - chromosome/sequence name
    /// * `start` - lower bound
    /// * `end` - upper bound
    ///
    fn get_statistics(&self, chr_name: &str, start: i64, end: i64) -> anyhow::Result<()> {
        let region = (chr_name, start, end);

        // todo: collect statistics

        if let Ok(freq_map) = get_dna_base_freq(&self.arc_bam, region) {
            for samp in freq_map.samples() {
                if let Some(stat) = freq_map.get_forward(samp) {
                    for bs in stat {
                        // let x = b.bi_allelic_stat();
                        // x.a1;
                        // x.a2;
                        // x.n1;
                        // x.n2;
                        bs.position();
                        // b.get(b)
                    }
                }
            }
        }

        Ok(())
    }

    fn populate_statistics() {

        // let forward_var_pos = Arc::new(Mutex::new(HashSet::new()));
        // let reverse_var_pos = Arc::new(Mutex::new(HashSet::new()));
    }

    /// Sift bases within the given region by checking `is_variable`
    ///
    /// * `chr_name` - chromosome/sequence name
    /// * `start` - lower bound
    /// * `end` - upper bound
    ///
    fn variable_bases(
        &self,
        chr_name: &str,
        start: i64,
        end: i64,
    ) -> anyhow::Result<DirectedPositions> {
        let region = (chr_name, start, end);
        let base_filter = rules::BaseFilters::new();
        let mut forward = vec![];
        let mut reverse = vec![];

        if let Ok(freq_map) = get_dna_base_freq(&self.arc_bam, region) {
            for samp in freq_map.samples() {
                if let Some(stats) = freq_map.get_forward(samp) {
                    for bs in stats {
                        if base_filter.is_variable(bs) {
                            forward.push(bs.position());
                        }
                    }
                }

                if let Some(stats) = freq_map.get_reverse(samp) {
                    for bs in stats {
                        if base_filter.is_variable(bs) {
                            reverse.push(bs.position());
                        }
                    }
                }
            }
        } else {
            return Err(anyhow::anyhow!("empty stat map"));
        }

        Ok(DirectedPositions {
            forward_positions: forward,
            reverse_positions: reverse,
        })
    }

    /// Sift all blocks
    ///
    /// * `chr_name` - chromosome/sequence name
    /// * `blocks` - genomic intervals
    ///
    fn find_variable_positions(
        &self,
        chr_name: &str,
        blocks: &Vec<(i64, i64)>,
    ) -> anyhow::Result<DirectedSets> {
        //
        // keep track of genomic locations
        let forward_var_pos = Arc::new(Mutex::new(HashSet::new()));

        let reverse_var_pos = Arc::new(Mutex::new(HashSet::new()));

        blocks
            .iter()
            .par_bridge()
            .try_for_each(|(lb, ub)| -> anyhow::Result<()> {
                if let Ok(positions) = self.variable_bases(chr_name, *lb, *ub) {
                    match forward_var_pos.try_lock() {
                        Ok(mut ret) => {
                            for j in positions.forward_positions.iter() {
                                ret.insert(*j);
                            }
                        }
                        Err(_) => {
                            return Err(anyhow::anyhow!("failed to lock forward"));
                        }
                    };

                    match reverse_var_pos.try_lock() {
                        Ok(mut ret) => {
                            for j in positions.reverse_positions.iter() {
                                ret.insert(*j);
                            }
                        }
                        Err(_) => {
                            return Err(anyhow::anyhow!("failed to lock reverse"));
                        }
                    };
                }

                Ok(())
            })?;

        let forward = forward_var_pos
            .lock()
            .expect("failed to copy forward")
            .clone();

        let reverse = reverse_var_pos
            .lock()
            .expect("failed to copy reverse")
            .clone();

        Ok(DirectedSets {
            forward_positions: forward,
            reverse_positions: reverse,
        })
    }
}
