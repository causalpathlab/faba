//! A side of the blockmodel, the cells' or the genes': each unit's links
//! per role, and the roles' global weights. The cells' side is all a
//! link's cell role touches outside its own reads, so it is what models of
//! several modalities, each with its own gene roles and reads, would share.

use super::Weights;
use crate::cluster::fast;
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::Distribution;

/// A chain's counts on one side: links per unit and role (`units x roles`)
/// and the global weights over the roles.
#[derive(Debug, Clone)]
pub struct Urn {
    pub(super) n: Vec<u32>,
    pub(super) beta: Vec<f64>,
}

/// One side of the model: how many units (cells or genes) and roles, and
/// the prior on each unit's role weights, `theta_c ~ Dirichlet(alpha
/// beta)` (cells) or `psi_g ~ Dirichlet(gene_alpha beta_gene)` (genes),
/// the global weights `beta` drawn (HDP) or fixed at `1 / roles`.
#[derive(Debug, Clone, Copy)]
pub struct Side {
    pub(super) units: usize,
    pub(super) roles: usize,
    pub(super) alpha: f64,
    pub(super) weights: Weights,
}

impl Side {
    /// No links counted, uniform global weights.
    pub(super) fn urn(&self) -> Urn {
        Urn {
            n: vec![0; self.units * self.roles],
            beta: vec![1.0 / self.roles as f64; self.roles],
        }
    }

    /// A link of unit `c` in role `a` counted in (`add`) or out.
    pub(super) fn count(&self, urn: &mut Urn, c: usize, a: usize, add: bool) {
        let n = &mut urn.n[c * self.roles + a];
        if add {
            *n += 1;
        } else {
            *n -= 1;
        }
    }

    /// The log weight of role `a` for a link of unit `c` (that link
    /// counted out): the change in the unit's urn as role `a` gains it,
    /// `ln (n_ca + alpha beta_a)`.
    pub(super) fn ln_weight(&self, urn: &Urn, c: usize, a: usize) -> f64 {
        fast::urn_step(self.alpha * urn.beta[a], urn.n[c * self.roles + a])
    }

    /// Unit `c`'s urn: its links per role under its Dirichlet.
    fn term(&self, urn: &Urn, c: usize) -> f64 {
        fast::ln_urn(self.row(urn, c), self.alpha, &urn.beta)
    }

    /// Unit `c`'s links per role.
    fn row<'u>(&self, urn: &'u Urn, c: usize) -> &'u [u32] {
        &urn.n[c * self.roles..(c + 1) * self.roles]
    }

    /// The change in unit `c`'s urn as roles `r` and `r2` exchange their
    /// links.
    pub(super) fn exchange(&self, urn: &Urn, c: usize, pair: (usize, usize)) -> f64 {
        fast::urn_exchange(self.row(urn, c), self.alpha, &urn.beta, pair)
    }

    /// Every unit's urn.
    pub(super) fn ln_joint(&self, urn: &Urn) -> f64 {
        (0..self.units).map(|c| self.term(urn, c)).sum()
    }

    /// Each unit's posterior-mean role weights, `units x roles`, row-major.
    pub(super) fn means(&self, urn: &Urn) -> Vec<f64> {
        role_weights(&urn.n, self.alpha, &urn.beta)
    }

    /// The global weights from their conditional (HDP only): each unit's
    /// table counts in each role (an Antoniak draw), then
    /// `beta ~ Dirichlet(gamma / roles + tables)`. One role's weight is 1,
    /// and is not drawn.
    pub(super) fn resample(&self, urn: &mut Urn, rng: &mut SmallRng) {
        let (Weights::Hdp(gamma), 2..) = (self.weights, self.roles) else {
            return;
        };
        let roles = self.roles;
        // A unit's first link in a role opens a table for certain; the
        // `j`-th after it one with probability `ab / (ab + j)`.
        let ab: Vec<f64> = urn.beta.iter().map(|b| self.alpha * b).collect();
        let mut tables = vec![0.0f64; roles];
        for row in urn.n.chunks_exact(roles) {
            for ((t, &ab), &n) in tables.iter_mut().zip(&ab).zip(row) {
                if n > 0 {
                    *t += 1.0;
                }
                for j in 1..n {
                    if rng.random::<f64>() * (ab + j as f64) < ab {
                        *t += 1.0;
                    }
                }
            }
        }
        for (b, &t) in urn.beta.iter_mut().zip(&tables) {
            *b = rand_distr::Gamma::new(gamma / roles as f64 + t, 1.0)
                .unwrap()
                .sample(rng)
                .max(1e-300);
        }
        let total: f64 = urn.beta.iter().sum();
        urn.beta.iter_mut().for_each(|b| *b /= total);
    }
}

/// Row-major role weights given the counts: each row's share of its links
/// in each role around `conc * base`.
fn role_weights(counts: &[u32], conc: f64, base: &[f64]) -> Vec<f64> {
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
