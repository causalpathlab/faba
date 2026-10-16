//! Posterior summaries of the chains, label-free and without a cells x cells
//! matrix.
//!
//! The posterior similarity `P_ij` (how often cells `i` and `j` share a
//! cluster) is never formed. The least-squares partition (Dahl 2006) is the
//! draw `pi` minimising `sum_ij (1[pi_i = pi_j] - P_ij)^2`, which up to a
//! constant is `sum_a n_a^2 - (2/T) sum_t sum_ab n_ab(pi, t)^2`, from the
//! contingency counts `n_ab` of `pi` against each draw `t`. A cell's
//! confidence, its mean co-clustering with its own cluster, is likewise
//! `(1/T) sum_t n_(pi_i, t_i) / n_(pi_i)`.

use super::crp::{Chain, Sample};
use legume_numeric::mcmc::engine::split_rhat_chains;
use rayon::prelude::*;

/// What the chains say, together.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Posterior {
    /// `K`, clusters of at least `min_cells` cells, and the share of draws with that many.
    pub min_cells: usize,
    pub k: Vec<(usize, f64)>,
    pub alpha_mean: f64,
    /// Split R-hat across chains of `K` and of the log posterior: infinite
    /// (written `"inf"`) when chains sit still at different values.
    #[serde(serialize_with = "finite_or_inf")]
    pub rhat_k: f64,
    #[serde(serialize_with = "finite_or_inf")]
    pub rhat_log_posterior: f64,
    pub chains: usize,
    pub draws: usize,
    /// Each chain's mean log posterior: one far below the others sits in a
    /// worse mode.
    pub chain_log_posterior: Vec<f64>,
    /// Adjusted Rand index between each pair of chains' last draws (1: the same partition).
    pub chain_agreement: Vec<f64>,
    pub accepted_splits: usize,
    pub accepted_merges: usize,
}

/// JSON has no infinity: write it as `"inf"`.
fn finite_or_inf<S: serde::Serializer>(v: &f64, s: S) -> Result<S::Ok, S::Error> {
    if v.is_finite() {
        s.serialize_f64(*v)
    } else {
        s.serialize_str("inf")
    }
}

/// The chosen partition and each cell's confidence in it.
pub struct Partition {
    /// Clusters numbered 0.. by size, largest first.
    pub z: Vec<u16>,
    pub k: usize,
    pub confidence: Vec<f64>,
}

/// Counts of `x` against `y` (labels below `kx`, `ky`), row-major.
fn contingency(x: &[u16], kx: usize, y: &[u16], ky: usize) -> Vec<u32> {
    let mut n = vec![0u32; kx * ky];
    for (a, b) in x.iter().zip(y) {
        n[*a as usize * ky + *b as usize] += 1;
    }
    n
}

fn sum_sq(n: &[u32]) -> f64 {
    n.iter().map(|&v| (v as f64) * (v as f64)).sum()
}

/// Clusters of at least `min_cells` cells in a draw.
pub fn substantial(s: &Sample, min_cells: usize) -> usize {
    let mut n = vec![0usize; s.k];
    s.z.iter().for_each(|&z| n[z as usize] += 1);
    n.iter().filter(|&&v| v >= min_cells).count()
}

pub fn posterior(chains: &[Chain], min_cells: usize) -> Posterior {
    let draws: Vec<&Sample> = chains.iter().flat_map(|c| &c.samples).collect();
    let mut hist: std::collections::BTreeMap<usize, usize> = Default::default();
    for s in &draws {
        *hist.entry(substantial(s, min_cells)).or_default() += 1;
    }
    let t = draws.len().max(1) as f64;
    let series = |f: &dyn Fn(&Sample) -> f64| -> Vec<Vec<f32>> {
        chains
            .iter()
            .map(|c| c.samples.iter().map(|s| f(s) as f32).collect())
            .collect()
    };
    Posterior {
        min_cells,
        k: hist.into_iter().map(|(k, n)| (k, n as f64 / t)).collect(),
        alpha_mean: draws.iter().map(|s| s.alpha).sum::<f64>() / t,
        rhat_k: split_rhat_chains(&series(&|s| substantial(s, min_cells) as f64)) as f64,
        rhat_log_posterior: split_rhat_chains(&series(&|s| s.log_post)) as f64,
        chains: chains.len(),
        draws: draws.len(),
        chain_log_posterior: chains
            .iter()
            .map(|c| {
                c.samples.iter().map(|s| s.log_post).sum::<f64>() / c.samples.len().max(1) as f64
            })
            .collect(),
        chain_agreement: {
            let last: Vec<&Sample> = chains.iter().filter_map(|c| c.samples.last()).collect();
            let mut v = Vec::new();
            for a in 0..last.len() {
                for b in a + 1..last.len() {
                    v.push(ari(last[a], last[b]));
                }
            }
            v
        },
        accepted_splits: chains.iter().map(|c| c.splits).sum(),
        accepted_merges: chains.iter().map(|c| c.merges).sum(),
    }
}

/// Adjusted Rand index of two draws' partitions.
pub fn ari(x: &Sample, y: &Sample) -> f64 {
    let n = contingency(&x.z, x.k, &y.z, y.k);
    let c2 = |v: f64| v * (v - 1.0) / 2.0;
    let index: f64 = n.iter().map(|&v| c2(v as f64)).sum();
    let rows: f64 = (0..x.k)
        .map(|a| c2(n[a * y.k..(a + 1) * y.k].iter().sum::<u32>() as f64))
        .sum();
    let cols: f64 = (0..y.k)
        .map(|b| c2((0..x.k).map(|a| n[a * y.k + b]).sum::<u32>() as f64))
        .sum();
    let expected = rows * cols / c2(x.z.len() as f64);
    let max = (rows + cols) / 2.0;
    if max == expected {
        1.0
    } else {
        (index - expected) / (max - expected)
    }
}

/// The least-squares partition among the draws, and each cell's confidence.
pub fn partition(chains: &[Chain]) -> Partition {
    let draws: Vec<&Sample> = chains.iter().flat_map(|c| &c.samples).collect();
    assert!(!draws.is_empty(), "no draws to summarise");
    let t = draws.len() as f64;
    // Cross term of each candidate against every draw.
    let loss: Vec<f64> = draws
        .par_iter()
        .map(|cand| {
            let own = sum_sq(&contingency(&cand.z, cand.k, &cand.z, cand.k));
            let cross: f64 = draws
                .iter()
                .map(|d| sum_sq(&contingency(&cand.z, cand.k, &d.z, d.k)))
                .sum();
            own - 2.0 * cross / t
        })
        .collect();
    let best = (0..draws.len())
        .min_by(|&a, &b| loss[a].total_cmp(&loss[b]))
        .unwrap();
    let chosen = draws[best];
    // Renumber by size, largest first.
    let mut sizes: Vec<(usize, usize)> = (0..chosen.k)
        .map(|k| (chosen.z.iter().filter(|&&z| z as usize == k).count(), k))
        .collect();
    sizes.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut rank = vec![0u16; chosen.k];
    for (r, &(_, k)) in sizes.iter().enumerate() {
        rank[k] = r as u16;
    }
    let z: Vec<u16> = chosen.z.iter().map(|&k| rank[k as usize]).collect();
    let k = chosen.k;
    let n_own: Vec<f64> = (0..k)
        .map(|a| z.iter().filter(|&&x| x as usize == a).count() as f64)
        .collect();
    let mut confidence = vec![0.0f64; z.len()];
    for d in &draws {
        let n = contingency(&z, k, &d.z, d.k);
        for (i, c) in confidence.iter_mut().enumerate() {
            let a = z[i] as usize;
            *c += n[a * d.k + d.z[i] as usize] as f64 / n_own[a];
        }
    }
    confidence.iter_mut().for_each(|c| *c /= t);
    Partition { z, k, confidence }
}
