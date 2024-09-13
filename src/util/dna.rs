use crate::util::bam::*;
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Aux, Read};
use std::cmp::max;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

///
/// DNA base-level frequency
#[derive(Debug)]
pub struct DnaBaseFreq {
    a: f32,
    t: f32,
    g: f32,
    c: f32,
    gpos: i64,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Dna {
    A,
    T,
    G,
    C,
}

impl DnaBaseFreq {
    /// Identify most frequent allele and report it with empirical
    /// frequency. Caveat: We may have ambiguous results if two or
    /// more alelles tie.
    pub fn major_allele_frequency(&self) -> Option<(Dna, f32)> {
        let tot = self.total();
        if tot > 0f32 {
            let max = self.a.max(self.t).max(self.g).max(self.c);

            if max == self.a {
                Some((Dna::A, self.a / tot))
            } else if max == self.t {
                Some((Dna::T, self.t / tot))
            } else if max == self.g {
                Some((Dna::G, self.g / tot))
            } else {
                Some((Dna::C, self.c / tot))
            }
        } else {
            None
        }
    }

    /// Total number of reads mapped on this BP
    pub fn total(&self) -> f32 {
        self.a + self.t + self.g + self.c
    }
}

/// Extract DNA base pair frequency tables in multi-threaded visits
/// over BAM file reader. Here, we only go through aligned reads.
///
pub fn get_dna_base_freq(
    arc_bam: &Arc<Mutex<bam::IndexedReader>>,
    region: (&str, i64, i64),
) -> anyhow::Result<DnaFreqMap> {
    let (_, lb, ub) = region;

    let mut bam_reader = arc_bam.lock().expect("unable to lock the reader");

    bam_reader
        .fetch(region)
        .expect("unable to fetch the region");

    if lb >= ub {
        return Err(anyhow::anyhow!("lb >= ub"));
    }

    // (cb) -> forward/reverse frequency vectors
    let mut ret = DnaFreqMap::new();
    ret.new_sample(&Sample::Combined, lb, ub);
    // dbg!("added combined");

    // Iter aligned read and reference positions on a basepair level
    // https://docs.rs/rust-htslib/latest/src/rust_htslib/bam/ext.rs.html#135
    // [read_pos, genome_pos]

    for rr in bam_reader.rc_records() {
        match rr {
            Ok(rec) => {
                if rec.is_duplicate() {
                    continue;
                }

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

                for [rpos, gpos] in rec.aligned_pairs() {
                    let (r, g, v) = (rpos as usize, gpos as usize, gpos - lb);

                    if g < (lb as usize) || g >= (ub as usize) || v < 0 {
                        continue;
                    }

                    let bp = seq[r];

                    let freq = match rec.is_reverse() {
                        true => ret.get_reverse_bp_mut(&sample_id, v as usize),
                        _ => ret.get_forward_bp_mut(&sample_id, v as usize),
                    };

                    if let Some(freq) = freq {
                        debug_assert_eq!(freq.gpos, gpos);
                        match bp {
                            b'A' | b'a' => freq.a += 1f32,
                            b'T' | b't' => freq.t += 1f32,
                            b'G' | b'g' => freq.g += 1f32,
                            b'C' | b'c' => freq.c += 1f32,
                            _ => (),
                        }
                    }
                }
            }
            _ => {
                // report error message?
            }
        }
    }

    Ok(ret)
}

/// DNA frequency map from forward and reverse strands
pub struct DnaFreqMap {
    forward: HashMap<usize, Vec<DnaBaseFreq>>,
    reverse: HashMap<usize, Vec<DnaBaseFreq>>,
    samp2id: HashMap<Sample, usize>,
    id2samp: Vec<Sample>,
}

impl DnaFreqMap {
    fn new() -> Self {
        DnaFreqMap {
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
                forward_freq.push(DnaBaseFreq {
                    a: 0f32,
                    t: 0f32,
                    g: 0f32,
                    c: 0f32,
                    gpos: g,
                });
                reverse_freq.push(DnaBaseFreq {
                    a: 0f32,
                    t: 0f32,
                    g: 0f32,
                    c: 0f32,
                    gpos: g,
                });
            }
            debug_assert_eq!(forward_freq.len(), nn);
            debug_assert_eq!(reverse_freq.len(), nn);
        }
    }

    pub fn get_forward(&self, key: &Sample) -> Option<&Vec<DnaBaseFreq>> {
        self.samp2id.get(&key).and_then(|id| self.forward.get(id))
    }

    pub fn get_reverse(&self, key: &Sample) -> Option<&Vec<DnaBaseFreq>> {
        self.samp2id.get(&key).and_then(|id| self.reverse.get(id))
    }

    pub fn get_forward_bp_mut(&mut self, key: &Sample, at: usize) -> Option<&mut DnaBaseFreq> {
        self.samp2id
            .get(&key)
            .and_then(|id| self.forward.get_mut(id))
            .and_then(|vv| vv.get_mut(at))
    }

    pub fn get_reverse_bp_mut(&mut self, key: &Sample, at: usize) -> Option<&mut DnaBaseFreq> {
        self.samp2id
            .get(&key)
            .and_then(|id| self.reverse.get_mut(id))
            .and_then(|vv| vv.get_mut(at))
    }
}
