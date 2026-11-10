//! The start of a variational fit by k-means: the cells clustered by their
//! normalized profiles (spherical k-means, k-means++ seeds), the genes by
//! their profiles over those clusters, and the globals of `q` read off the
//! labels, softened (see [`super`]). Only the nonzero counts are visited;
//! the centroids are roles x genes, never cells x genes.

use super::{exposure, Fit, Globals};
use crate::cluster::count::{column_totals, sum_parts, Counts, Settings};
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::kmeans::{kmeans_rows_seeded, KmeansMetric, KmeansRowsOpts};
use rand::rngs::SmallRng;
use rand::RngExt;
use rayon::prelude::*;

/// The share of each label spread evenly over every role, so that no role
/// starts empty (a role at nothing stays there).
const SOFT: f64 = 0.5;

/// Lloyd iterations at most, for the cells and for the genes.
pub const ITERS: usize = 30;

/// Label `label`'s softened share of role `j`, of `k` roles.
fn soft(label: usize, j: usize, k: usize) -> f64 {
    SOFT / k as f64 + if j == label { 1.0 - SOFT } else { 0.0 }
}

/// What every restart's start reads alike: each cell's total and its
/// profile, `ln(1 + x / total * 1e4)` of unit length, one per nonzero.
pub struct Profiles {
    total: Vec<f64>,
    v: Vec<f32>,
}

impl Profiles {
    pub fn new(counts: &Counts) -> Self {
        let total: Vec<f64> = (0..counts.cells())
            .into_par_iter()
            .map(|c| counts.of(c).1.iter().map(|&x| x as f64).sum())
            .collect();
        let v = (0..counts.cells())
            .into_par_iter()
            .flat_map_iter(|c| {
                let per = 1e4 / total[c].max(1.0);
                let row: Vec<f64> = (counts.of(c).1.iter())
                    .map(|&x| (x as f64 * per).ln_1p())
                    .collect();
                let norm = row.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-12);
                row.into_iter().map(move |v| (v / norm) as f32)
            })
            .collect();
        Self { total, v }
    }
}

/// The start of a fit: the globals, and each cell's cluster, whose
/// softened shares give its first units by role.
pub struct Start {
    pub(super) globals: Globals,
    z: Vec<usize>,
}

impl Start {
    /// Cell `c`'s first units in role `a` of `aa`, of its `total`.
    pub(super) fn units(&self, c: usize, a: usize, aa: usize, total: f64) -> f64 {
        total * soft(self.z[c], a, aa)
    }
}

/// The k-means start: `settings.cell_roles` clusters of cells, then
/// `settings.gene_roles` of genes, each seeded by k-means++ from `rng`;
/// `tick` told each Lloyd iteration of the cells done (1..), `stop` ends
/// them early. Refused when there is no count.
pub fn kmeans_start(
    counts: &Counts,
    profiles: &Profiles,
    settings: &Settings,
    rng: &mut SmallRng,
    tick: &(dyn Fn(usize) + Sync),
    stop: &(dyn Fn() -> bool + Sync),
) -> anyhow::Result<Start> {
    anyhow::ensure!(!counts.x.is_empty(), "no counts to start from");
    let (gg, aa, bb) = (counts.genes, settings.cell_roles, settings.gene_roles);
    let rows = Rows {
        ptr: &counts.ptr,
        index: &counts.gene,
        v: &profiles.v,
        dim: gg,
    };
    let z = spherical(&rows, aa, rng, tick, stop);
    let total = &profiles.total;

    // Each gene's profile over the cell clusters: its mean share of a
    // cell's counts in each; then the genes clustered by it.
    let mut size = vec![0f64; aa];
    z.iter().for_each(|&a| size[a] += 1.0);
    let h = fold_cells(counts, gg * aa, |h, c| {
        let (g, x) = counts.of(c);
        let per = 1.0 / total[c].max(1.0) / size[z[c]];
        for (&g, &x) in g.iter().zip(x) {
            h[g as usize * aa + z[c]] += x as f64 * per;
        }
    });
    let h = DMatrix::from_row_iterator(gg, aa, h.iter().map(|&v| v as f32));
    let opts = KmeansRowsOpts {
        k: bb,
        max_iter: ITERS,
        seed: rng.random(),
        metric: KmeansMetric::Cosine,
        min_changed_frac: 0.0,
        init_sample: 0,
    };
    let y = kmeans_rows_seeded(&h, &opts).labels;

    // The labels' statistics, each label softened: every count to its
    // cell's and gene's roles by their shares. As the shares are a constant
    // plus a label's, the sums need only the hard labels': the units per
    // block `k`, per cell cluster `t` and per gene `x`.
    let hard = fold_cells(counts, aa * bb + aa + gg, |s, c| {
        let (g, x) = counts.of(c);
        for (&g, &x) in g.iter().zip(x) {
            let g = g as usize;
            s[z[c] * bb + y[g]] += x as f64;
            s[aa * bb + aa + g] += x as f64;
        }
        s[aa * bb + z[c]] += total[c];
    });
    let (k, rest) = hard.split_at(aa * bb);
    let (t, x) = rest.split_at(aa);
    let all: f64 = t.iter().sum();
    let m: Vec<f64> = (0..gg * bb)
        .map(|i| x[i / bb] * soft(y[i / bb], i % bb, bb))
        .collect();
    // Each cell cluster's units by gene role, then each block's.
    let by_role = |a: usize, b: usize| t[a] * SOFT / bb as f64 + (1.0 - SOFT) * k[a * bb + b];
    let every: Vec<f64> = (0..bb)
        .map(|b| (0..aa).map(|a| by_role(a, b)).sum())
        .collect();
    let blocks: Vec<f64> = (0..aa * bb)
        .map(|i| SOFT / aa as f64 * every[i % bb] + (1.0 - SOFT) * by_role(i / bb, i % bb))
        .collect();
    let theta_total: Vec<f64> = (0..aa)
        .map(|a| SOFT / aa as f64 * all + (1.0 - SOFT) * t[a])
        .collect();
    // psi: each gene role's units spread over the genes, summing to one.
    let per_role = column_totals(&m, bb);
    let psi: Vec<f64> = (0..gg * bb)
        .map(|i| m[i] / per_role[i % bb].max(1e-12))
        .collect();
    let fit = Fit::new(counts, settings);
    let globals = globals(&fit, &m, &blocks, &theta_total, &psi);
    Ok(Start { globals, z })
}

/// `f` over every cell of `counts` into a sum of `len` zeros, the same
/// whatever the threads (see [`sum_parts`]).
fn fold_cells(counts: &Counts, len: usize, f: impl Fn(&mut [f64], usize) + Sync) -> Vec<f64> {
    sum_parts(counts.cells(), len, PARTS, f)
}

/// The fixed parts the sums over cells are split into.
const PARTS: usize = 64;

/// Sparse rows of unit length: row `i`'s entries `v[ptr[i]..ptr[i + 1]]`
/// at columns `index`, `dim` columns in all.
struct Rows<'a> {
    ptr: &'a [usize],
    index: &'a [u32],
    v: &'a [f32],
    dim: usize,
}

impl Rows<'_> {
    fn len(&self) -> usize {
        self.ptr.len() - 1
    }

    /// Row `i`'s entries, by column.
    fn row(&self, i: usize) -> impl Iterator<Item = (usize, f32)> + '_ {
        let r = self.ptr[i]..self.ptr[i + 1];
        (self.index[r.clone()].iter().zip(&self.v[r])).map(|(&j, &v)| (j as usize, v))
    }

    /// Row `i`'s dot product with each of the `k` dense centroids `mu`.
    fn dots(&self, i: usize, mu: &[f32], k: usize, out: &mut [f32]) {
        out[..k].fill(0.0);
        for (j, v) in self.row(i) {
            for (o, m) in out[..k].iter_mut().zip(mu.chunks_exact(self.dim)) {
                *o += v * m[j];
            }
        }
    }
}

/// Spherical k-means of the sparse `rows` into `k` clusters: k-means++
/// seeds (by the cosine distance `1 - dot`), then at most [`ITERS`] Lloyd
/// iterations, each told to `tick`, an emptied cluster seeded again with
/// the row farthest from its own. Each row's cluster.
fn spherical(
    rows: &Rows,
    k: usize,
    rng: &mut SmallRng,
    tick: &(dyn Fn(usize) + Sync),
    stop: &(dyn Fn() -> bool + Sync),
) -> Vec<usize> {
    let (n, dim) = (rows.len(), rows.dim);
    let mut mu = vec![0f32; k * dim];
    let set = |mu: &mut [f32], j: usize, i: usize| {
        let row = &mut mu[j * dim..(j + 1) * dim];
        row.fill(0.0);
        rows.row(i).for_each(|(c, v)| row[c] = v);
    };
    // k-means++: each next seed drawn by its distance to the nearest seed.
    let mut near = vec![1f32; n];
    for j in 0..k {
        let total: f64 = near.iter().map(|&d| d.max(0.0) as f64).sum();
        let i = match total > 0.0 {
            true => {
                let mut u = rng.random::<f64>() * total;
                (near.iter().position(|&d| {
                    u -= d.max(0.0) as f64;
                    u <= 0.0
                }))
                .unwrap_or(n - 1)
            }
            false => rng.random_range(0..n),
        };
        set(&mut mu, j, i);
        let seed = &mu[j * dim..(j + 1) * dim];
        near.par_iter_mut().enumerate().for_each(|(r, d)| {
            let mut dot = [0f32];
            rows.dots(r, seed, 1, &mut dot);
            *d = d.min(1.0 - dot[0]);
        });
    }
    let mut z = vec![usize::MAX; n];
    for iter in 1..=ITERS {
        // Every row is assigned once, whatever is asked.
        if iter > 1 && stop() {
            break;
        }
        // Assign: the nearest centroid, and how near.
        let best: Vec<(usize, f32)> = (0..n)
            .into_par_iter()
            .map_init(
                || vec![0f32; k],
                |out, i| {
                    rows.dots(i, &mu, k, out);
                    (0..k).fold((0, f32::MIN), |b, j| match out[j] > b.1 {
                        true => (j, out[j]),
                        false => b,
                    })
                },
            )
            .collect();
        let moved = (best.iter().zip(&z)).filter(|(b, z)| b.0 != **z).count();
        best.iter().zip(&mut z).for_each(|(b, z)| *z = b.0);
        // Update: each centroid the normalized sum of its rows, a thread's
        // sums at a time.
        let sums = sum_parts(n, k * dim + k, PARTS, |s: &mut [f32], i| {
            let j = z[i];
            s[k * dim + j] += 1.0;
            rows.row(i).for_each(|(c, v)| s[j * dim + c] += v);
        });
        mu.copy_from_slice(&sums[..k * dim]);
        let empty: Vec<usize> = (0..k).filter(|&j| sums[k * dim + j] == 0.0).collect();
        // The rows farthest from their centroids, only when some cluster
        // emptied.
        let mut far: Vec<usize> = (0..n).collect();
        if !empty.is_empty() && empty.len() < n {
            far.select_nth_unstable_by(empty.len(), |&a, &b| best[a].1.total_cmp(&best[b].1));
        }
        for (&j, &i) in empty.iter().zip(&far) {
            set(&mut mu, j, i);
            z[i] = j;
        }
        for j in (0..k).filter(|j| !empty.contains(j)) {
            let row = &mut mu[j * dim..(j + 1) * dim];
            let norm = row.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-12);
            row.iter_mut().for_each(|v| *v /= norm);
        }
        tick(iter);
        if moved == 0 {
            break;
        }
    }
    z
}

/// The globals of `q` from a point of the model: units per gene and gene
/// role `m` (genes x B) and per block `blocks` (A x B), the cells' theta
/// summed per role, and psi (genes x B). Each factor's conditional given
/// those, its rates from `q`'s own moments in turn.
fn globals(fit: &Fit, m: &[f64], blocks: &[f64], theta_total: &[f64], psi: &[f64]) -> Globals {
    let pr = fit.prior;
    let (aa, bb) = (fit.aa, fit.bb);
    let gg = m.len() / bb;
    let psi_total = column_totals(psi, bb);
    let mut g = Globals {
        aa,
        bb,
        p: m.iter().map(|&m| pr.eta + m).collect(),
        pq: vec![0.0; gg * bb],
        qr: (psi.chunks_exact(bb))
            .map(|row| pr.f0 + row.iter().sum::<f64>())
            .collect(),
        l: blocks.iter().map(|&k| pr.kappa + k).collect(),
        lr: (0..aa * bb)
            .map(|k| pr.kappa + theta_total[k / bb] * psi_total[k % bb])
            .collect(),
    };
    let rho_shape = pr.rho_shape(bb);
    let f = exposure(&g, theta_total);
    for (row, qr) in g.pq.chunks_exact_mut(bb).zip(&g.qr) {
        let e_rho = rho_shape / qr;
        row.iter_mut().zip(&f).for_each(|(pq, f)| *pq = e_rho + f);
    }
    g
}
