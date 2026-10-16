//! A small dense picture of each cell, built once, for what the sampler only
//! needs to be roughly right: which cells are alike (split-merge anchors) and
//! where chains start. Nothing here enters the likelihood, so the sampler
//! stays exact whatever the picture.

use legume_numeric::matrix::dmatrix_io::DMatrix;
use legume_numeric::matrix::kmeans::{kmeans_rows_seeded, KmeansMetric, KmeansRowsOpts};
use legume_numeric::matrix::knn::knn_rows;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

/// Dimensions of the projection.
const DIM: usize = 32;

/// Fixed, so a run's anchors and starts do not depend on its `--seed`.
const SEED: u64 = 0x5eed_c0de;

/// Each cell's random projection of its per-site residuals,
/// `(on - p0 n) / sqrt(n p0 (1 - p0) + 1)` at the site's pooled fraction
/// `p0`, scaled to unit length: cells x `DIM`.
pub fn embed(ptr: &[usize], site: &[u32], on: &[f64], off: &[f64], p0: &[f64]) -> DMatrix<f32> {
    let cells = ptr.len() - 1;
    let mut rng = SmallRng::seed_from_u64(SEED);
    let proj: Vec<f32> = (0..p0.len() * DIM)
        .map(|_| StandardNormal.sample(&mut rng))
        .collect();
    let mut z = DMatrix::<f32>::zeros(cells, DIM);
    for c in 0..cells {
        let mut row = [0.0f32; DIM];
        for i in ptr[c]..ptr[c + 1] {
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

/// Each cell's `k` nearest cells in the picture.
pub fn neighbours(z: &DMatrix<f32>, k: usize) -> Vec<Vec<u32>> {
    let (idx, _) = knn_rows(z, k);
    idx.into_iter()
        .map(|n| n.into_iter().map(|j| j as u32).collect())
        .collect()
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
