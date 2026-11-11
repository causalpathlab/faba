//! Small functions the fits and their tables share.

/// A site's pooled converted fraction from `on` of `total` reads, with half
/// a read of pseudocount: the centre of every cluster's prior there.
pub(in crate::cluster) fn pooled_fraction(on: f64, total: f64) -> f64 {
    (on + 0.5) / (total + 1.0)
}

/// Each of `sites` sites' [`pooled_fraction`] over entries `(site, on,
/// off)`.
pub(in crate::cluster) fn pooled_fractions(
    site: &[u32],
    on: &[f64],
    off: &[f64],
    sites: usize,
) -> Vec<f64> {
    let (mut pon, mut tot) = (vec![0.0f64; sites], vec![0.0f64; sites]);
    for (i, &s) in site.iter().enumerate() {
        pon[s as usize] += on[i];
        tot[s as usize] += on[i] + off[i];
    }
    (0..sites)
        .map(|s| pooled_fraction(pon[s], tot[s]))
        .collect()
}

/// The logit `ln (p / (1 - p))`, the inverse of the sigmoid.
pub(in crate::cluster) fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

/// Each row's most weighted column (`k` columns, row-major), ties to the
/// first.
pub(in crate::cluster) fn main_labels(weights: &[f32], k: usize) -> Vec<u16> {
    weights
        .chunks(k)
        .map(|row| {
            (0..k)
                .max_by(|&a, &b| row[a].total_cmp(&row[b]).then(b.cmp(&a)))
                .unwrap_or(0) as u16
        })
        .collect()
}
