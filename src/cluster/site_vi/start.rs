//! The start of the site blockmodel's variational fit by k-means: the cells
//! clustered by their links' residuals (spherical k-means++, sparse over the
//! genes), the genes by their mean residual in each cell cluster, and `q`
//! read off the labels, softened (see [`super`] and the count fit's start,
//! whose k-means this shares).

use super::{omega, Fit, Globals, Locals};
use crate::cluster::count::vi::{soft, spherical, Rows};
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::kmeans::{kmeans_rows_seeded, KmeansMetric, KmeansRowsOpts};
use rand::rngs::SmallRng;
use rand::RngExt;
use rayon::prelude::*;

/// Lloyd iterations at most for the genes (the cells' are the count
/// start's).
const ITERS: usize = 30;

/// The start of a fit: its globals and locals.
pub struct Start {
    pub(super) globals: Globals,
    pub(super) locals: Locals,
}

/// What every restart's start reads alike: each link's residual (its
/// `sum (m - n p) / sqrt(sum n p (1 - p) + 1)` at its sites' baselines `p`),
/// and each cell's row of them over its genes, less its mean and of unit
/// length.
pub struct Profiles {
    residual: Vec<f64>,
    v: Vec<f32>,
}

impl Profiles {
    pub fn new(fit: &Fit) -> Self {
        let data = fit.data;
        let links = data.links;
        let residual: Vec<f64> = (0..links.len())
            .into_par_iter()
            .map(|l| {
                let (sites, m, u) = links.of(l);
                let (mut top, mut var) = (0.0, 0.0);
                for ((&s, &m), &u) in sites.iter().zip(m).zip(u) {
                    let p = 1.0 / (1.0 + (-data.mu[s as usize]).exp());
                    top += m - (m + u) * p;
                    var += (m + u) * p * (1.0 - p);
                }
                top / (var + 1.0).sqrt()
            })
            .collect();
        let v = (0..data.cells())
            .into_par_iter()
            .flat_map_iter(|c| {
                let r = &residual[data.links_of(c)];
                let mean = r.iter().sum::<f64>() / r.len().max(1) as f64;
                let norm = (r.iter().map(|x| (x - mean).powi(2)).sum::<f64>())
                    .sqrt()
                    .max(1e-12);
                r.iter().map(move |x| ((x - mean) / norm) as f32)
            })
            .collect();
        Self { residual, v }
    }
}

/// The k-means start: `A` clusters of cells by their links' residuals
/// (see [`Profiles`]); then `B` clusters of genes by their mean residual in
/// each cell cluster. `tick` told each Lloyd iteration of the cells (1..),
/// `stop` ends them early. Refused when there is no link.
pub fn kmeans_start(
    fit: &Fit,
    profiles: &Profiles,
    rng: &mut SmallRng,
    tick: &(dyn Fn(usize) + Sync),
    stop: &(dyn Fn() -> bool + Sync),
) -> anyhow::Result<Start> {
    let data = fit.data;
    let links = data.links;
    anyhow::ensure!(links.len() > 0, "no reads to start from");
    let (aa, bb) = (fit.prior.cell_roles, fit.prior.gene_roles);
    let (cc, gg) = (data.cells(), data.genes);
    let residual = &profiles.residual;
    let rows = Rows {
        ptr: &links.cell_ptr,
        index: &data.link_gene,
        v: &profiles.v,
        dim: gg,
    };
    let z = spherical(&rows, aa, rng, tick, stop);

    // Each gene's mean residual in each cell cluster; the genes clustered.
    let (mut h, mut n) = (vec![0f64; gg * aa], vec![0f64; gg * aa]);
    for (c, &zc) in z.iter().enumerate() {
        for l in data.links_of(c) {
            let at = data.link_gene[l] as usize * aa + zc;
            h[at] += residual[l];
            n[at] += 1.0;
        }
    }
    let h = DMatrix::from_row_iterator(
        gg,
        aa,
        h.iter().zip(&n).map(|(h, n)| (h / n.max(1.0)) as f32),
    );
    let opts = KmeansRowsOpts {
        k: bb,
        max_iter: ITERS,
        seed: rng.random(),
        metric: KmeansMetric::Cosine,
        min_changed_frac: 0.0,
        init_sample: 0,
    };
    let y = kmeans_rows_seeded(&h, &opts).labels;

    // q from the softened labels: each cell's and gene's links spread over
    // its roles; the shifts at their prior; each cell's level fitted alone.
    let pr = fit.prior;
    let (a0, b0) = (pr.alpha / aa as f64, pr.gene_alpha / bb as f64);
    let mut gene_links = vec![0f64; gg];
    data.link_gene
        .iter()
        .for_each(|&g| gene_links[g as usize] += 1.0);
    let p = (0..gg * bb)
        .map(|i| b0 + gene_links[i / bb] * soft(y[i / bb], i % bb, bb))
        .collect();
    let k = aa * bb;
    let shift_var = pr.shift_sd * pr.shift_sd;
    let globals = Globals {
        p,
        d: vec![0.0; k],
        v: vec![shift_var; k],
    };
    let cell_var = pr.cell_sd * pr.cell_sd;
    // Each cell's role weights, and its level fitted alone.
    type Cell = (Vec<f64>, (f64, f64));
    let cells: Vec<Cell> = (0..cc)
        .into_par_iter()
        .map(|c| {
            let own = data.links_of(c).len() as f64;
            let t = (0..aa).map(|a| a0 + own * soft(z[c], a, aa)).collect();
            let (mut prec, mut num) = (1.0 / cell_var, 0.0);
            for l in data.links_of(c) {
                let (sites, m, u) = links.of(l);
                for ((&site, &m), &u) in sites.iter().zip(m).zip(u) {
                    let (n, mu) = (m + u, data.mu[site as usize]);
                    let w = omega(n, (mu * mu + shift_var).sqrt());
                    prec += w;
                    num += m - n / 2.0 - w * mu;
                }
            }
            (t, (num / prec, 1.0 / prec))
        })
        .collect();
    let locals = Locals {
        t: cells.iter().flat_map(|c| c.0.iter().copied()).collect(),
        e: cells.iter().map(|c| c.1 .0).collect(),
        s: cells.iter().map(|c| c.1 .1).collect(),
    };
    Ok(Start { globals, locals })
}
