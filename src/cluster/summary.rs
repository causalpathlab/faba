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
    /// Held-out reads' log-likelihood per read, over every kept draw.
    pub heldout_per_read: Option<f64>,
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

/// Counts of `x` against `y` (labels below `kx`, `ky`), row-major, into `n`.
fn contingency_into(n: &mut Vec<u32>, x: &[u16], kx: usize, y: &[u16], ky: usize) {
    n.clear();
    n.resize(kx * ky, 0);
    for (a, b) in x.iter().zip(y) {
        n[*a as usize * ky + *b as usize] += 1;
    }
}

fn contingency(x: &[u16], kx: usize, y: &[u16], ky: usize) -> Vec<u32> {
    let mut n = Vec::new();
    contingency_into(&mut n, x, kx, y, ky);
    n
}

/// Cells per cluster of labels below `k`.
pub(super) fn sizes(z: &[u16], k: usize) -> Vec<usize> {
    let mut n = vec![0usize; k];
    z.iter().for_each(|&x| n[x as usize] += 1);
    n
}

/// Labels below `k` renumbered by cluster size, largest first (ties by label).
pub(super) fn renumber_by_size(z: &mut [u16], k: usize) {
    let n = sizes(z, k);
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| n[b].cmp(&n[a]).then(a.cmp(&b)));
    let mut rank = vec![0u16; k];
    for (r, &a) in order.iter().enumerate() {
        rank[a] = r as u16;
    }
    z.iter_mut().for_each(|x| *x = rank[*x as usize]);
}

/// Least-squares candidates scored: draws evenly thinned to at most this
/// many, each scored against every draw.
const CANDIDATES: usize = 100;

fn sum_sq(n: &[u32]) -> f64 {
    n.iter().map(|&v| (v as f64) * (v as f64)).sum()
}

/// Clusters of at least `min_cells` cells in a draw.
pub fn substantial(s: &Sample, min_cells: usize) -> usize {
    sizes(&s.z, s.k).iter().filter(|&&v| v >= min_cells).count()
}

/// Split R-hat of a scalar's per-chain series, centred on its overall mean
/// before the f32 cast (which would round a log joint in the millions to
/// whole units).
pub fn centred_rhat(chains: &[Vec<f64>]) -> f64 {
    let n = chains.iter().map(Vec::len).sum::<usize>().max(1) as f64;
    let mean = chains.iter().flatten().sum::<f64>() / n;
    let series: Vec<Vec<f32>> = chains
        .iter()
        .map(|c| c.iter().map(|&x| (x - mean) as f32).collect())
        .collect();
    split_rhat_chains(&series) as f64
}

/// `f` of each pair of `xs`, in order.
pub fn pairwise<T>(xs: &[T], f: impl Fn(&T, &T) -> f64) -> Vec<f64> {
    let mut out = Vec::new();
    for a in 0..xs.len() {
        for b in a + 1..xs.len() {
            out.push(f(&xs[a], &xs[b]));
        }
    }
    out
}

/// Adjusted Rand index between each pair of labellings.
pub fn pairwise_ari(labels: &[Vec<u16>]) -> Vec<f64> {
    pairwise(labels, |a, b| ari(a, b))
}

/// Each chain's mean log joint, the split R-hat across chains, and the best
/// chain (highest mean), from each chain's log joint at its kept draws.
pub fn chain_scores(log_joint: &[Vec<f64>]) -> (Vec<f64>, f64, usize) {
    let means: Vec<f64> = log_joint
        .iter()
        .map(|c| c.iter().sum::<f64>() / c.len().max(1) as f64)
        .collect();
    let best = (0..means.len())
        .max_by(|&a, &b| means[a].total_cmp(&means[b]))
        .unwrap_or(0);
    (means, centred_rhat(log_joint), best)
}

/// Each row's most weighted column (`k` columns, row-major), ties to the
/// first.
pub fn main_labels(weights: &[f32], k: usize) -> Vec<u16> {
    weights
        .chunks(k)
        .map(|row| {
            (0..k)
                .max_by(|&a, &b| row[a].total_cmp(&row[b]).then(b.cmp(&a)))
                .unwrap_or(0) as u16
        })
        .collect()
}

/// Held-out log-likelihood per read from each draw's `(log-likelihood,
/// reads)`, none without held-out reads.
pub fn per_read(draws: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    let (ll, n) = draws.fold((0.0, 0.0), |(l, n), d| (l + d.0, n + d.1));
    (n > 0.0).then(|| ll / n)
}

pub fn posterior(chains: &[Chain], min_cells: usize) -> Posterior {
    let draws: Vec<&Sample> = chains.iter().flat_map(|c| &c.samples).collect();
    // Substantial clusters per draw, per chain, once for both uses.
    let ks: Vec<Vec<f32>> = chains
        .iter()
        .map(|c| {
            c.samples
                .iter()
                .map(|s| substantial(s, min_cells) as f32)
                .collect()
        })
        .collect();
    let mut hist: std::collections::BTreeMap<usize, usize> = Default::default();
    for &k in ks.iter().flatten() {
        *hist.entry(k as usize).or_default() += 1;
    }
    let t = draws.len().max(1) as f64;
    Posterior {
        min_cells,
        k: hist.into_iter().map(|(k, n)| (k, n as f64 / t)).collect(),
        alpha_mean: draws.iter().map(|s| s.alpha).sum::<f64>() / t,
        rhat_k: split_rhat_chains(&ks) as f64,
        // Centred before the f32 cast, which would round a log posterior in
        // the millions to whole units.
        rhat_log_posterior: centred_rhat(
            &chains
                .iter()
                .map(|c| c.samples.iter().map(|s| s.log_post).collect())
                .collect::<Vec<_>>(),
        ),
        chains: chains.len(),
        draws: draws.len(),
        chain_log_posterior: chains
            .iter()
            .map(|c| {
                c.samples.iter().map(|s| s.log_post).sum::<f64>() / c.samples.len().max(1) as f64
            })
            .collect(),
        chain_agreement: pairwise_ari(
            &chains
                .iter()
                .filter_map(|c| c.samples.last().map(|s| s.z.clone()))
                .collect::<Vec<_>>(),
        ),
        accepted_splits: chains.iter().map(|c| c.splits).sum(),
        accepted_merges: chains.iter().map(|c| c.merges).sum(),
        heldout_per_read: per_read(draws.iter().map(|d| d.heldout)),
    }
}

/// Adjusted Rand index of two labellings of the same cells.
pub fn ari(x: &[u16], y: &[u16]) -> f64 {
    let k = |z: &[u16]| z.iter().map(|&v| v as usize + 1).max().unwrap_or(0);
    let (kx, ky) = (k(x), k(y));
    let n = contingency(x, kx, y, ky);
    let c2 = |v: f64| v * (v - 1.0) / 2.0;
    let index: f64 = n.iter().map(|&v| c2(v as f64)).sum();
    let rows: f64 = (0..kx)
        .map(|a| c2(n[a * ky..(a + 1) * ky].iter().sum::<u32>() as f64))
        .sum();
    let cols: f64 = (0..ky)
        .map(|b| c2((0..kx).map(|a| n[a * ky + b]).sum::<u32>() as f64))
        .sum();
    let expected = rows * cols / c2(x.len() as f64);
    let max = (rows + cols) / 2.0;
    if max == expected {
        1.0
    } else {
        (index - expected) / (max - expected)
    }
}

/// The least-squares partition among the draws (scored over an even
/// thinning of them, see [`CANDIDATES`]), and each cell's confidence.
pub fn partition(chains: &[Chain]) -> Partition {
    let draws: Vec<&Sample> = chains.iter().flat_map(|c| &c.samples).collect();
    assert!(!draws.is_empty(), "no draws to summarise");
    let t = draws.len() as f64;
    let step = draws.len().div_ceil(CANDIDATES);
    // Cross term of each candidate against every draw.
    let loss: Vec<(usize, f64)> = (0..draws.len())
        .step_by(step)
        .collect::<Vec<_>>()
        .par_iter()
        .map(|&i| {
            let cand = draws[i];
            let mut n = Vec::new();
            let mut cross = 0.0;
            for d in &draws {
                contingency_into(&mut n, &cand.z, cand.k, &d.z, d.k);
                cross += sum_sq(&n);
            }
            let own: f64 = sizes(&cand.z, cand.k).iter().map(|&v| (v * v) as f64).sum();
            (i, own - 2.0 * cross / t)
        })
        .collect();
    let best = loss.iter().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap().0;
    let chosen = draws[best];
    let (mut z, k) = (chosen.z.clone(), chosen.k);
    renumber_by_size(&mut z, k);
    let n_own: Vec<f64> = sizes(&z, k).into_iter().map(|v| v as f64).collect();
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
