//! A chain's roles matched to a reference chain's, since each chain numbers
//! its roles arbitrarily: cell roles by their cells' memberships, then gene
//! roles by their genes', the block shifts breaking ties, each an optimal
//! one-to-one assignment (Hungarian algorithm).

use super::Means;

/// Most the block shifts' distance, scaled to at most 1, weighs against
/// the gene memberships' overlap: enough to tell apart gene roles whose
/// genes tie, too little to overrule overlaps that differ by more than
/// `B` times it in total.
const SHIFT_TIE: f64 = 1e-9;

/// A chain's roles in the reference's numbering: its cell role `a` is the
/// reference's `cell[a]`, its gene role `b` the reference's `gene[b]`;
/// with how far each reference role's members match its counterpart's
/// (weighted Jaccard of the memberships, `sum min / sum max` over the
/// cells or genes: 1 when the same), and whether it is in use (some cell's
/// or gene's main role, in either chain), by reference role.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Alignment {
    pub cell: Vec<usize>,
    pub gene: Vec<usize>,
    pub cell_quality: Vec<f64>,
    pub gene_quality: Vec<f64>,
    pub cell_used: Vec<bool>,
    pub gene_used: Vec<bool>,
}

impl Alignment {
    /// Chain `m`'s alignment to itself (`aa` cell roles, `bb` gene roles):
    /// every role its own, of quality 1, in use as in `m`.
    pub fn identity(m: &Means, aa: usize, bb: usize) -> Self {
        Alignment {
            cell: (0..aa).collect(),
            gene: (0..bb).collect(),
            cell_quality: vec![1.0; aa],
            gene_quality: vec![1.0; bb],
            cell_used: used(&m.cell, aa),
            gene_used: used(&m.gene, bb),
        }
    }

    /// The worst match of any role in use, cell or gene (1 with none): a
    /// role no cell or gene has as its main one holds only its weight's
    /// smoothing, which need not match from chain to chain.
    pub fn worst(&self) -> f64 {
        let cell = self.cell_quality.iter().zip(&self.cell_used);
        let gene = self.gene_quality.iter().zip(&self.gene_used);
        cell.chain(gene)
            .filter(|q| *q.1)
            .map(|q| *q.0)
            .fold(1.0, f64::min)
    }

    /// Block `k` (`a B + b`) in the reference's numbering.
    pub fn block(&self, k: usize) -> usize {
        let bb = self.gene.len();
        self.cell[k / bb] * bb + self.gene[k % bb]
    }

    /// Every block's number in the reference's, by block (see [`block`]).
    ///
    /// [`block`]: Alignment::block
    pub fn blocks(&self) -> Vec<usize> {
        (0..self.cell.len() * self.gene.len())
            .map(|k| self.block(k))
            .collect()
    }
}

/// Weighted Jaccard of column `x` of `p` and column `y` of `q` (`k`
/// columns each, row-major): 1 for two empty columns.
fn jaccard(p: &[f32], x: usize, q: &[f32], y: usize, k: usize) -> f64 {
    let (mut low, mut high) = (0.0f64, 0.0f64);
    for (r, s) in p.chunks(k).zip(q.chunks(k)) {
        let (u, v) = (f64::from(r[x]), f64::from(s[y]));
        low += u.min(v);
        high += u.max(v);
    }
    if high > 0.0 {
        low / high
    } else {
        1.0
    }
}

/// Whether each of `p`'s `k` columns (row-major) is some row's largest,
/// or ties it.
fn used(p: &[f32], k: usize) -> Vec<bool> {
    let mut out = vec![false; k];
    for row in p.chunks(k) {
        let top = row.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        row.iter()
            .zip(&mut out)
            .for_each(|(&w, u)| *u = *u || w == top);
    }
    out
}

/// Whether each reference role is in use in either chain: in `reference`
/// (`ours`), or the chain's role matched to it, by `to`, in the chain
/// (`theirs`).
fn either(ours: Vec<bool>, theirs: &[bool], to: &[usize]) -> Vec<bool> {
    let mut out = ours;
    to.iter()
        .zip(theirs)
        .for_each(|(&j, &u)| out[j] = out[j] || u);
    out
}

/// `chain`'s roles matched to `reference`'s (`aa` cell roles, `bb` gene
/// roles).
pub fn align(reference: &Means, chain: &Means, aa: usize, bb: usize) -> Alignment {
    // Chain role `i` against reference role `j` at `i n + j`.
    let cell_score: Vec<f64> = (0..aa * aa)
        .map(|at| jaccard(&chain.cell, at / aa, &reference.cell, at % aa, aa))
        .collect();
    let cell = assign(&cell_score, aa);
    let gene_overlap: Vec<f64> = (0..bb * bb)
        .map(|at| jaccard(&chain.gene, at / bb, &reference.gene, at % bb, bb))
        .collect();
    let distance: Vec<f64> = (0..bb * bb)
        .map(|at| {
            let (y, x) = (at / bb, at % bb);
            (0..aa)
                .map(|a| (chain.shift[a * bb + y] - reference.shift[cell[a] * bb + x]).abs())
                .sum()
        })
        .collect();
    let farthest = distance.iter().copied().fold(0.0, f64::max);
    let gene_score: Vec<f64> = gene_overlap
        .iter()
        .zip(&distance)
        .map(|(&o, &d)| {
            if farthest > 0.0 {
                o - SHIFT_TIE * d / farthest
            } else {
                o
            }
        })
        .collect();
    let gene = assign(&gene_score, bb);
    let quality = |score: &[f64], to: &[usize]| {
        let n = to.len();
        let mut q = vec![0.0; n];
        to.iter()
            .enumerate()
            .for_each(|(i, &j)| q[j] = score[i * n + j]);
        q
    };
    Alignment {
        cell_quality: quality(&cell_score, &cell),
        gene_quality: quality(&gene_overlap, &gene),
        cell_used: either(used(&reference.cell, aa), &used(&chain.cell, aa), &cell),
        gene_used: either(used(&reference.gene, bb), &used(&chain.gene, bb), &gene),
        cell,
        gene,
    }
}

/// The one-to-one assignment of rows to columns maximising the total
/// `score` (`n x n`, row-major): each row's column, a score not finite
/// counted below every finite one. The Hungarian algorithm with
/// potentials, `O(n^3)`.
pub fn assign(score: &[f64], n: usize) -> Vec<usize> {
    let finite = score.iter().copied().filter(|x| x.is_finite());
    let worst = finite.clone().fold(f64::INFINITY, f64::min);
    let spread = finite.fold(f64::NEG_INFINITY, f64::max) - worst;
    let floor = if worst.is_finite() {
        worst - 1.0 - spread * n as f64
    } else {
        0.0
    };
    // Minimise the cost `-score`, rows and columns numbered from 1, column
    // 0 a sentinel: `owner[j]` the row on column `j`, `via[j]` the column
    // before `j` on the augmenting path.
    let cost = |i: usize, j: usize| {
        let x = score[(i - 1) * n + j - 1];
        -(if x.is_finite() { x } else { floor })
    };
    let (mut u, mut v) = (vec![0.0f64; n + 1], vec![0.0f64; n + 1]);
    let (mut owner, mut via) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    for i in 1..=n {
        owner[0] = i;
        let mut j0 = 0;
        let mut slack = vec![f64::INFINITY; n + 1];
        let mut used = vec![false; n + 1];
        loop {
            used[j0] = true;
            let i0 = owner[j0];
            let (mut delta, mut j1) = (f64::INFINITY, 0);
            for j in 1..=n {
                if used[j] {
                    continue;
                }
                let reduced = cost(i0, j) - u[i0] - v[j];
                if reduced < slack[j] {
                    slack[j] = reduced;
                    via[j] = j0;
                }
                if slack[j] < delta {
                    delta = slack[j];
                    j1 = j;
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[owner[j]] += delta;
                    v[j] -= delta;
                } else {
                    slack[j] -= delta;
                }
            }
            j0 = j1;
            if owner[j0] == 0 {
                break;
            }
        }
        while j0 != 0 {
            let j1 = via[j0];
            owner[j0] = owner[j1];
            j0 = j1;
        }
    }
    let mut column = vec![0usize; n];
    for j in 1..=n {
        column[owner[j] - 1] = j - 1;
    }
    column
}
