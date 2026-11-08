//! The cells' side of the blockmodel: each cell's links per cell role, and
//! the cell roles' global weights. This is all a link's cell role touches
//! outside its own reads, so it is what models of several modalities, each
//! with its own gene roles and reads, would share.

use super::Weights;
use crate::cluster::crp::ln_urn;
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::Distribution;

/// A chain's cell-side counts: links per cell and cell role (`cells x A`)
/// and the global weights over the cell roles.
#[derive(Debug, Clone)]
pub struct Urn {
    pub(super) n: Vec<u32>,
    pub(super) beta: Vec<f64>,
}

/// The cells' side of the model: how many cells and cell roles, and the
/// prior on each cell's role weights, `theta_c ~ Dirichlet(alpha beta)`.
#[derive(Debug, Clone, Copy)]
pub struct Cells {
    pub(super) cells: usize,
    pub(super) roles: usize,
    pub(super) alpha: f64,
    pub(super) weights: Weights,
}

impl Cells {
    /// No links counted, uniform global weights.
    pub(super) fn urn(&self) -> Urn {
        Urn {
            n: vec![0; self.cells * self.roles],
            beta: vec![1.0 / self.roles as f64; self.roles],
        }
    }

    /// A link of cell `c` in role `a` counted in (`add`) or out.
    pub(super) fn count(&self, urn: &mut Urn, c: usize, a: usize, add: bool) {
        let n = &mut urn.n[c * self.roles + a];
        if add {
            *n += 1;
        } else {
            *n -= 1;
        }
    }

    /// The log weight of role `a` for a link of cell `c` (that link
    /// counted out): `ln (n_ca + alpha beta_a)`.
    pub(super) fn ln_weight(&self, urn: &Urn, c: usize, a: usize) -> f64 {
        (urn.n[c * self.roles + a] as f64 + self.alpha * urn.beta[a]).ln()
    }

    /// Cell `c`'s urn: its links per role under its Dirichlet.
    pub(super) fn term(&self, urn: &Urn, c: usize) -> f64 {
        let row = &urn.n[c * self.roles..(c + 1) * self.roles];
        ln_urn(row, self.alpha, &urn.beta)
    }

    /// Every cell's urn.
    pub(super) fn ln_joint(&self, urn: &Urn) -> f64 {
        (0..self.cells).map(|c| self.term(urn, c)).sum()
    }

    /// Each cell's posterior-mean role weights, `cells x A`, row-major.
    pub(super) fn means(&self, urn: &Urn) -> Vec<f64> {
        role_weights(&urn.n, self.alpha, &urn.beta)
    }

    /// The global weights from their conditional (HDP only): each cell's
    /// table counts in each role (an Antoniak draw), then
    /// `beta ~ Dirichlet(gamma / A + tables)`.
    pub(super) fn resample(&self, urn: &mut Urn, rng: &mut SmallRng) {
        let Weights::Hdp(gamma) = self.weights else {
            return;
        };
        let aa = self.roles;
        let mut tables = vec![0.0f64; aa];
        for c in 0..self.cells {
            for (a, t) in tables.iter_mut().enumerate() {
                let ab = self.alpha * urn.beta[a];
                for j in 0..urn.n[c * aa + a] {
                    if rng.random::<f64>() < ab / (ab + j as f64) {
                        *t += 1.0;
                    }
                }
            }
        }
        let draws: Vec<f64> = tables
            .iter()
            .map(|&t| {
                rand_distr::Gamma::new(gamma / aa as f64 + t, 1.0)
                    .unwrap()
                    .sample(rng)
                    .max(1e-300)
            })
            .collect();
        let total: f64 = draws.iter().sum();
        urn.beta = draws.iter().map(|d| d / total).collect();
    }
}

/// Row-major role weights given the counts: each row's share of its links
/// in each role around `conc * base`.
pub(super) fn role_weights(counts: &[u32], conc: f64, base: &[f64]) -> Vec<f64> {
    counts
        .chunks(base.len())
        .flat_map(|row| {
            let total = row.iter().sum::<u32>() as f64 + conc;
            row.iter()
                .zip(base)
                .map(move |(&x, &b)| (x as f64 + conc * b) / total)
        })
        .collect()
}
