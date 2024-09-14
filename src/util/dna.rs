use crate::util::bam::*;

use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux, Read};
use std::cmp::max;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum Dna {
    A,
    T,
    G,
    C,
}

#[allow(dead_code)]
pub struct DnaBaseStat {
    data: [(Dna, f32); 4],
    gpos: i64,
}

#[allow(dead_code)]
impl DnaBaseStat {
    fn new(gpos: i64) -> Self {
        DnaBaseStat {
            data: [
                (Dna::A, 0f32),
                (Dna::T, 0f32),
                (Dna::G, 0f32),
                (Dna::C, 0f32),
            ],
            gpos,
        }
    }

    pub fn position(&self) -> i64 {
        self.gpos
    }

    pub fn set(&mut self, b: Dna, val: f32) {
        match b {
            Dna::A => self.data[0].1 = val,
            Dna::T => self.data[1].1 = val,
            Dna::G => self.data[2].1 = val,
            Dna::C => self.data[3].1 = val,
        }
    }

    pub fn add(&mut self, b: Dna, val: f32) {
        match b {
            Dna::A => self.data[0].1 += val,
            Dna::T => self.data[1].1 += val,
            Dna::G => self.data[2].1 += val,
            Dna::C => self.data[3].1 += val,
        }
    }

    pub fn get(&self, b: Dna) -> f32 {
        match b {
            Dna::A => self.data[0].1,
            Dna::T => self.data[1].1,
            Dna::G => self.data[2].1,
            Dna::C => self.data[3].1,
        }
    }

    pub fn most_frequent(&self) -> &(Dna, f32) {
        self.data
            .iter()
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            .unwrap()
    }

    pub fn second_most_frequent(&self) -> &(Dna, f32) {
        let mfa = self.most_frequent();
        self.data
            .iter()
            .filter(|s| s.0 != mfa.0)
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            .unwrap()
    }

    pub fn bi_allelic_stat(&self) -> BiAllele {
        let fst = self.most_frequent();
        let snd = self
            .data
            .iter()
            .filter(|s| s.0 != fst.0)
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            .unwrap();

        BiAllele {
            a1: fst.0.clone(),
            a2: snd.0.clone(),
            n1: fst.1,
            n2: snd.1,
        }
    }
}

#[allow(dead_code)]
pub struct BiAllele {
    pub a1: Dna,
    pub a2: Dna,
    pub n1: f32,
    pub n2: f32,
}

/// Extract DNA base pair frequency tables in multi-threaded visits
/// over BAM file reader. Here, we only go through aligned reads.
///
pub fn get_dna_base_freq(
    arc_bam: &Arc<Mutex<bam::IndexedReader>>,
    region: (&str, i64, i64),
) -> anyhow::Result<DnaStatMap> {
    let (_, lb, ub) = region;

    let mut bam_reader = arc_bam.lock().expect("unable to lock the reader");

    bam_reader
        .fetch(region)
        .expect("unable to fetch the region");

    if lb >= ub {
        return Err(anyhow::anyhow!("lb >= ub"));
    }

    let bam_records: Vec<bam::Record> = bam_reader
        .records()
        .into_iter()
        .filter_map(Result::ok)
        .filter(|rec| !rec.is_duplicate())
        .collect();

    if bam_records.len() < 1 {
        return Err(anyhow::anyhow!("Empty region"));
    }

    // map: sample -> forward/reverse frequency vectors
    let mut ret = DnaStatMap::new();
    ret.new_sample(&Sample::Combined, lb, ub);
    // dbg!("added combined");

    for rec in bam_records {
        let mut sample_id = Sample::Combined;

        // https://docs.rs/rust-htslib/0.47.0/rust_htslib/bam/record/enum.Aux.html
        // extract 10x cell barcode
        if let Ok(aux) = rec.aux(b"CB") {
            if let Aux::String(cb) = aux {
                sample_id = Sample::Barcode(cb.into());
                if !ret.has_sample(&sample_id) {
                    // dbg!("added new cell barcode");
                    ret.new_sample(&sample_id, lb, ub);
                }
            }
        }

        // extract 10x UMI barcode
        // if let Ok(umi) = rec.aux(b"UB") {
        //     dbg!(umi);
        // }

        let seq = rec.seq().as_bytes();

        //
        // Iter aligned read and reference positions on a basepair level
        // https://docs.rs/rust-htslib/latest/src/rust_htslib/bam/ext.rs.html#135
        // [read_pos, genome_pos]
        //
        for [rpos, gpos] in rec.aligned_pairs() {
            let (r, g, v) = (rpos as usize, gpos as usize, gpos - lb);

            if g < (lb as usize) || g >= (ub as usize) || v < 0 {
                continue;
            }

            let bp = seq[r];

            let freq = match rec.is_reverse() {
                true => ret.get_reverse_base_mut(&sample_id, v as usize),
                _ => ret.get_forward_base_mut(&sample_id, v as usize),
            };

            if let Some(freq) = freq {
                debug_assert_eq!(freq.gpos, gpos);
                match bp {
                    b'A' | b'a' => freq.add(Dna::A, 1.),
                    b'T' | b't' => freq.add(Dna::T, 1.),
                    b'G' | b'g' => freq.add(Dna::G, 1.),
                    b'C' | b'c' => freq.add(Dna::C, 1.),
                    _ => (),
                }
            }
        }
    }

    Ok(ret)
}

/// DNA frequency map from forward and reverse strands
#[allow(dead_code)]
pub struct DnaStatMap {
    forward: HashMap<usize, Vec<DnaBaseStat>>,
    reverse: HashMap<usize, Vec<DnaBaseStat>>,
    samp2id: HashMap<Sample, usize>,
    id2samp: Vec<Sample>,
}

#[allow(dead_code)]
impl DnaStatMap {
    /// Create Dna Stat Map containing both forward and reverse
    /// directions.
    ///
    /// map: sample_id -> dna stat vector
    fn new() -> Self {
        DnaStatMap {
            forward: HashMap::new(),
            reverse: HashMap::new(),
            samp2id: HashMap::new(),
            id2samp: vec![],
        }
    }

    pub fn has_sample(&self, key: &Sample) -> bool {
        self.samp2id.contains_key(key)
    }

    pub fn samples(&self) -> &Vec<Sample> {
        &self.id2samp
    }

    // pub fn nsamples(&self) -> usize {
    //     self.id2samp.len()
    // }

    pub fn new_sample(&mut self, key: &Sample, lb: i64, ub: i64) {
        if !self.has_sample(&key) {
            let id = self.id2samp.len();
            self.samp2id.insert(key.clone(), id); //
            self.id2samp.push(key.clone()); // check

            debug_assert_eq!(self.id2samp.len(), id + 1);
            debug_assert_eq!(self.samp2id.len(), id + 1);

            let nn = max(ub - lb, 0i64) as usize;

            let forward_freq = self
                .forward
                .entry(id)
                .or_insert_with(|| Vec::with_capacity(nn));

            let reverse_freq = self
                .reverse
                .entry(id)
                .or_insert_with(|| Vec::with_capacity(nn));

            for g in lb..ub {
                forward_freq.push(DnaBaseStat::new(g));
                reverse_freq.push(DnaBaseStat::new(g));
            }
            debug_assert_eq!(forward_freq.len(), nn);
            debug_assert_eq!(reverse_freq.len(), nn);
        }
    }

    pub fn get_forward(&self, key: &Sample) -> Option<&Vec<DnaBaseStat>> {
        self.samp2id.get(&key).and_then(|id| self.forward.get(id))
    }

    pub fn get_reverse(&self, key: &Sample) -> Option<&Vec<DnaBaseStat>> {
        self.samp2id.get(&key).and_then(|id| self.reverse.get(id))
    }

    pub fn get_forward_base_mut(&mut self, key: &Sample, at: usize) -> Option<&mut DnaBaseStat> {
        self.samp2id
            .get(&key)
            .and_then(|id| self.forward.get_mut(id))
            .and_then(|vv| vv.get_mut(at))
    }

    pub fn get_reverse_base_mut(&mut self, key: &Sample, at: usize) -> Option<&mut DnaBaseStat> {
        self.samp2id
            .get(&key)
            .and_then(|id| self.reverse.get_mut(id))
            .and_then(|vv| vv.get_mut(at))
    }
}

// Old implementation
// /// DNA base-level frequency
// #[derive(Debug)]
// pub struct DnaBaseFreq {
//     a: f32,
//     t: f32,
//     g: f32,
//     c: f32,
//     gpos: i64,
// }
//
// impl DnaBaseFreq {
//     /// Identify most frequent allele and report it with empirical
//     /// frequency. Caveat: We may have ambiguous results if two or
//     /// more alelles tie.
//     pub fn major_allele_frequency(&self) -> Option<(Dna, f32)> {
//         let tot = self.total();
//         if tot > 0f32 {
//             let max = self.a.max(self.t).max(self.g).max(self.c);
//             if max == self.a {
//                 Some((Dna::A, self.a / tot))
//             } else if max == self.t {
//                 Some((Dna::T, self.t / tot))
//             } else if max == self.g {
//                 Some((Dna::G, self.g / tot))
//             } else {
//                 Some((Dna::C, self.c / tot))
//             }
//         } else {
//             None
//         }
//     }
//     /// Total number of reads mapped on this BP
//     pub fn total(&self) -> f32 {
//         self.a + self.t + self.g + self.c
//     }
// }
