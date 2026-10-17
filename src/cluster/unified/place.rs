//! Cells outside the fit, placed against a fitted chain.

use super::{ln_mean, slice, Links, Means, Model};
use crate::cluster::crp::{log_sum_exp, sample_log};
use rand::rngs::SmallRng;
use rand::RngExt;
use std::cell::RefCell;

/// Sweeps of a placed cell's own Gibbs, and those discarded before its
/// averages begin.
const PLACE_SWEEPS: usize = 60;
const PLACE_BURN_IN: usize = 20;

/// Placed cell `c`'s own random stream, from `seed`.
fn cell_rng(seed: u64, c: usize) -> SmallRng {
    <SmallRng as rand::SeedableRng>::seed_from_u64(
        seed ^ (c as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15),
    )
}

impl Model<'_> {
    /// Each site's log fractions `(ln p, ln (1 - p))` in each block
    /// (`sites x K`) at a chain's `means`: the profile's posterior mean
    /// given the block's mean reads there and its mean shift (finite
    /// strength only; at infinite strength they hang on the cell's level).
    fn plug_in(&self, means: &Means) -> Vec<(f64, f64)> {
        let kk = self.blocks();
        if !self.finite() {
            return Vec::new();
        }
        (0..self.mu.len() * kk)
            .map(|at| {
                let (s, k) = (at / kk, at % kk);
                ln_mean(
                    self.pseudo(s, means.shift[k]),
                    means.site_on[at],
                    means.site_off[at],
                )
            })
            .collect()
    }

    /// Cells outside the fit (`links`, numbered by the model's genes),
    /// placed against a fitted chain's `means` (with the sites'): each
    /// cell's links relabelled by collapsed Gibbs over its cell roles, the
    /// gene roles summed out at the genes' mean weights, alternating at
    /// infinite strength with its level by slice sampling, and its role
    /// weights and level averaged over the sweeps after a burn-in. The
    /// cell's own urn around the mean global weights and its level are the
    /// model's; the gene weights, profiles and shifts are plug-in
    /// estimates, not integrated over. Each cell's role weights, `cells x
    /// A`, row-major, and its level (none at finite strength).
    pub fn place(&self, links: &Links, means: &Means, seed: u64) -> (Vec<f32>, Vec<f64>) {
        use rayon::prelude::*;
        let bb = self.settings.gene_roles;
        let kk = self.blocks();
        let ln_p = self.plug_in(means);
        let link_gene = links.genes_of(self.gene_of);
        let ln_gene: Vec<f64> = means
            .gene
            .iter()
            .map(|&w| (w as f64).max(1e-300).ln())
            .collect();
        let prior: Vec<f64> = means.beta.iter().map(|b| self.settings.alpha * b).collect();
        // Link `l`'s reads' log-likelihood in cell role `a` at level `v`
        // (read at infinite strength only), its gene role summed out, each
        // gene role's term put in `terms` first.
        let like = |l: usize, a: usize, v: f64, terms: &mut [f64]| {
            let (s, m, u) = links.of(l);
            let g = link_gene[l] as usize;
            for (b, x) in terms.iter_mut().enumerate() {
                let k = a * bb + b;
                *x = ln_gene[g * bb + b]
                    + if self.finite() {
                        (0..s.len())
                            .map(|i| {
                                let (p, q) = ln_p[s[i] as usize * kk + k];
                                m[i] * p + u[i] * q
                            })
                            .sum::<f64>()
                    } else {
                        (0..s.len())
                            .map(|i| self.binomial(s[i] as usize, m[i], u[i], means.shift[k] + v))
                            .sum::<f64>()
                    };
            }
            log_sum_exp(terms.iter().copied())
        };
        let (weights, levels): (Vec<Vec<f32>>, Vec<f64>) = (0..links.cells)
            .into_par_iter()
            .map(|c| self.place_cell(links, c, &like, &prior, seed))
            .unzip();
        let levels = if self.finite() { Vec::new() } else { levels };
        (weights.concat(), levels)
    }

    /// Placed cell `c`'s role weights and level (see [`Self::place`]): its
    /// links' cell roles by Gibbs under `like` (link, cell role, level and
    /// a gene-role buffer) and its urn around `prior`, its level by slice
    /// sampling at infinite strength. At finite strength `like` does not
    /// read the level, so it is taken once per link and role.
    fn place_cell(
        &self,
        links: &Links,
        c: usize,
        like: &impl Fn(usize, usize, f64, &mut [f64]) -> f64,
        prior: &[f64],
        seed: u64,
    ) -> (Vec<f32>, f64) {
        let aa = self.settings.cell_roles;
        let mut rng = cell_rng(seed, c);
        let own: Vec<usize> = (links.cell_ptr[c]..links.cell_ptr[c + 1]).collect();
        let terms = RefCell::new(vec![0.0f64; self.settings.gene_roles]);
        let score = |l: usize, a: usize, v: f64| like(l, a, v, &mut terms.borrow_mut());
        let table: Vec<f64> = if self.finite() {
            own.iter()
                .flat_map(|&l| (0..aa).map(move |a| (l, a)))
                .map(|(l, a)| score(l, a, 0.0))
                .collect()
        } else {
            Vec::new()
        };
        // Link `own[i]` in cell role `a` at level `v`.
        let at = |i: usize, a: usize, v: f64| {
            if self.finite() {
                table[i * aa + a]
            } else {
                score(own[i], a, v)
            }
        };
        let mut row = vec![0u32; aa];
        let mut z: Vec<usize> = own.iter().map(|_| rng.random_range(0..aa)).collect();
        z.iter().for_each(|&a| row[a] += 1);
        let total = z.len() as f64 + self.settings.alpha;
        let var = self.settings.cell_sd * self.settings.cell_sd;
        let (mut v, mut level) = (0.0f64, 0.0f64);
        let mut weights = vec![0.0f64; aa];
        let mut w = vec![0.0f64; aa];
        for sweep in 0..PLACE_SWEEPS {
            for i in 0..own.len() {
                row[z[i]] -= 1;
                for (a, x) in w.iter_mut().enumerate() {
                    *x = (row[a] as f64 + prior[a]).ln() + at(i, a, v);
                }
                z[i] = sample_log(&mut w, &mut rng);
                row[z[i]] += 1;
            }
            if !self.finite() {
                let target = |x: f64| {
                    -x * x / (2.0 * var)
                        + z.iter().enumerate().map(|(i, &a)| at(i, a, x)).sum::<f64>()
                };
                v = slice(v, &target, 1.0, &mut rng);
            }
            if sweep >= PLACE_BURN_IN {
                for (a, x) in weights.iter_mut().enumerate() {
                    *x += (row[a] as f64 + prior[a]) / total;
                }
                level += v;
            }
        }
        let n = (PLACE_SWEEPS - PLACE_BURN_IN) as f64;
        (
            weights.into_iter().map(|x| (x / n) as f32).collect(),
            level / n,
        )
    }
}
