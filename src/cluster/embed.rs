//! A small dense picture of each cell, built once, for where chains start.
//! Nothing here enters the likelihood, so the sampler stays exact whatever
//! the picture.

use super::read::Weighted;
use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::kmeans::{kmeans_rows_seeded, KmeansMetric, KmeansRowsOpts};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

/// Dimensions of the projection.
const DIM: usize = 32;

/// Fixed, so a run's starts do not depend on its `--seed`.
const SEED: u64 = 0x5eed_c0de;

/// Each cell's random projection of its per-site residuals,
/// `(on - p0 n) / sqrt(n p0 (1 - p0) + 1)` at the site's pooled fraction
/// `p0`, scaled to unit length: cells x `DIM`.
pub fn embed(reads: &Weighted, p0: &[f64]) -> DMatrix<f32> {
    let cells = reads.cells();
    let mut rng = SmallRng::seed_from_u64(SEED);
    let proj: Vec<f32> = (0..p0.len() * DIM)
        .map(|_| StandardNormal.sample(&mut rng))
        .collect();
    let mut z = DMatrix::<f32>::zeros(cells, DIM);
    for c in 0..cells {
        let mut row = [0.0f32; DIM];
        let (site, on, off) = reads.of(c);
        for i in 0..site.len() {
            let s = site[i] as usize;
            let (n, p) = (on[i] + off[i], p0[s]);
            let r = ((on[i] - p * n) / (n * p * (1.0 - p) + 1.0).sqrt()) as f32;
            for (d, x) in row.iter_mut().enumerate() {
                *x += r * proj[s * DIM + d];
            }
        }
        let norm = row.iter().map(|x| x * x).sum::<f32>().sqrt().max(1e-12);
        for (d, x) in row.iter().enumerate() {
            z[(c, d)] = x / norm;
        }
    }
    z
}

/// A k-means partition of the picture into `k` groups.
pub fn kmeans(z: &DMatrix<f32>, k: usize, seed: u64) -> Vec<u32> {
    let fit = kmeans_rows_seeded(
        z,
        &KmeansRowsOpts {
            k,
            max_iter: 20,
            seed,
            metric: KmeansMetric::Cosine,
            min_changed_frac: 0.001,
            init_sample: 0,
        },
    );
    fit.labels.into_iter().map(|l| l as u32).collect()
}
