//! The heatmaps' numbers and orders: each row's z across the cell modules,
//! and a diagonal order of rows and columns.

use super::Count;

/// A row's methylated fractions standardized across its cell modules;
/// none for a cell module without reads.
pub fn standardized(reads: &[Count]) -> Vec<Option<f64>> {
    let f: Vec<f64> = reads.iter().filter_map(|r| r.fraction()).collect();
    let n = f.len().max(1) as f64;
    let mean = f.iter().sum::<f64>() / n;
    let sd = (f.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
    let sd = if sd > 0.0 { sd } else { 1.0 };
    reads
        .iter()
        .map(|r| r.fraction().map(|x| (x - mean) / sd))
        .collect()
}

/// A diagonal order of a matrix's rows and columns by reciprocal averaging
/// on its positive part: a row sits by the columns where it is high, a
/// column by the rows. Missing entries count as zero.
pub fn diagonal(z: &[Vec<Option<f64>>]) -> (Vec<usize>, Vec<usize>) {
    let cols = z.first().map_or(0, Vec::len);
    let p = |r: usize, c: usize| z[r][c].unwrap_or(0.0).max(0.0);
    let mut col_pos: Vec<usize> = (0..cols).collect();
    let mut row_pos: Vec<usize> = (0..z.len()).collect();
    for _ in 0..50 {
        row_pos = rank(
            (0..z.len())
                .map(|r| mean_at((0..cols).map(|c| (p(r, c), col_pos[c]))))
                .collect(),
        );
        let next = rank(
            (0..cols)
                .map(|c| mean_at((0..z.len()).map(|r| (p(r, c), row_pos[r]))))
                .collect(),
        );
        if next == col_pos {
            break;
        }
        col_pos = next;
    }
    (order(&row_pos), order(&col_pos))
}

/// Each score's rank, ties by index, so the order settles.
fn rank(score: Vec<f64>) -> Vec<usize> {
    let mut by: Vec<usize> = (0..score.len()).collect();
    by.sort_by(|&a, &b| score[a].total_cmp(&score[b]).then(a.cmp(&b)));
    order(&by)
}

/// The weighted mean of positions `(weight, position)`; last for a line
/// with no weight.
fn mean_at(weights: impl Iterator<Item = (f64, usize)>) -> f64 {
    let (w, s) = weights.fold((0.0, 0.0), |(w, s), (x, at)| (w + x, s + x * at as f64));
    if w > 0.0 {
        s / w
    } else {
        f64::MAX
    }
}

/// The inverse of a permutation: which item is at each place.
fn order(pos: &[usize]) -> Vec<usize> {
    let mut out = vec![0; pos.len()];
    for (i, &k) in pos.iter().enumerate() {
        out[k] = i;
    }
    out
}
