use crate::util::dna::*;
// use std::cmp::{max, min};
// use fastapprox::faster as fa;;

#[allow(dead_code)]
pub struct BaseFilters {
    max_major_allele_cutoff: f32,
    min_minor_allele_cutoff: f32,
}

#[allow(dead_code)]
impl BaseFilters {
    pub fn new() -> Self {
        BaseFilters {
            max_major_allele_cutoff: 1_f32 - 1e-4_f32,
            min_minor_allele_cutoff: 1e-4_f32,
        }
    }

    pub fn b_allele_frequency(&self, stat: &DnaBaseStat) -> f32 {
        let stat = stat.bi_allelic_stat();
        stat.n1 / (stat.n1 + stat.n2).max(1_f32)
    }

    pub fn is_variable(&self, stat: &DnaBaseStat) -> bool {
        let stat = stat.bi_allelic_stat();
        stat.n1 > 0_f32 && stat.n2 > 0_f32
    }

    pub fn is_near_zero_variance(&self, stat: &DnaBaseStat) -> bool {
        stat.most_frequent().1 > self.max_major_allele_cutoff
    }
}
