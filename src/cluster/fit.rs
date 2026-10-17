//! Given the chosen partition: each cluster's reads per site, cells outside
//! the fit assigned to a cluster, and every cell's level offset.

use super::crp::{ln_sigmoid, pooled_fraction};
use super::read::{Reads, Weighted};
use super::summary::Partition;

/// Pseudo-reads pulling a batch's fraction at a site toward the site's pooled
/// one, so a thinly read site is barely reweighted.
const SHRINK: f64 = 20.0;

/// Largest density ratio a read's weight starts from (truncated importance
/// sampling); self-normalising then scales it, so a final weight can exceed
/// it (the largest is recorded).
const MAX_WEIGHT: f64 = 10.0;

/// Importance weights that take the batches out of the reads.
///
/// The target is reads without batch: at site `s` a read is converted with
/// the site's pooled fraction `p_s`; batch `b` delivers `p_sb` instead (its
/// own fraction there, shrunk toward `p_s`). Each read is weighted by the
/// ratio of the two, `p_s / p_sb` if converted and `(1 - p_s) / (1 - p_sb)`
/// if not, so a batch's weighted fraction at a site is the target's. Ratios
/// are capped at [`MAX_WEIGHT`], then scaled per batch and site so its
/// weighted reads add up to its raw reads (self-normalised): every batch
/// brings the evidence it has, no more and no less.
pub struct BatchWeights {
    /// Per entry of the reads, its converted and unconverted reads' weights.
    pub on: Vec<f64>,
    pub off: Vec<f64>,
    /// Per batch, Kish's effective share of its reads, `(sum w)^2 / (n sum w^2)`.
    pub effective_share: Vec<f64>,
    /// Share of reads whose ratio was capped, and the largest final weight.
    pub capped_share: f64,
    pub max_weight: f64,
}

/// [`BatchWeights`] for `reads`, `batch_of[c]` being cell `c`'s batch in `0..batches`.
pub fn batch_weights(
    reads: &Reads,
    batch_of: &[u16],
    batches: usize,
    sites: usize,
) -> BatchWeights {
    let (mut on, mut tot) = (vec![0.0f64; sites * batches], vec![0.0f64; sites * batches]);
    for (c, &b) in batch_of.iter().enumerate() {
        let (s, m, u) = reads.of(c);
        for i in 0..s.len() {
            let at = s[i] as usize * batches + b as usize;
            on[at] += m[i] as f64;
            tot[at] += (m[i] + u[i]) as f64;
        }
    }
    // Per site and batch: the two weights, capped, then self-normalised.
    let (mut w_on, mut w_off) = (vec![1.0f64; sites * batches], vec![1.0f64; sites * batches]);
    let (mut capped, mut reads_n) = (0.0f64, 0.0f64);
    for s in 0..sites {
        let r = s * batches..(s + 1) * batches;
        let p = pooled_fraction(on[r.clone()].iter().sum(), tot[r.clone()].iter().sum());
        for at in r {
            if tot[at] == 0.0 {
                continue;
            }
            let pb = (on[at] + SHRINK * p) / (tot[at] + SHRINK);
            let (raw_on, raw_off) = (p / pb, (1.0 - p) / (1.0 - pb));
            let (a, b) = (raw_on.min(MAX_WEIGHT), raw_off.min(MAX_WEIGHT));
            let (m, u) = (on[at], tot[at] - on[at]);
            capped += m * f64::from(u8::from(raw_on > MAX_WEIGHT))
                + u * f64::from(u8::from(raw_off > MAX_WEIGHT));
            reads_n += tot[at];
            let scale = tot[at] / (a * m + b * u);
            w_on[at] = a * scale;
            w_off[at] = b * scale;
        }
    }
    let mut out = BatchWeights {
        on: Vec::with_capacity(reads.site.len()),
        off: Vec::with_capacity(reads.site.len()),
        effective_share: vec![0.0; batches],
        capped_share: capped / reads_n.max(1.0),
        max_weight: 0.0,
    };
    let (mut sw, mut sw2, mut n) = (
        vec![0.0f64; batches],
        vec![0.0f64; batches],
        vec![0.0f64; batches],
    );
    for (c, &b) in batch_of.iter().enumerate() {
        let (s, m, u) = reads.of(c);
        let b = b as usize;
        for i in 0..s.len() {
            let at = s[i] as usize * batches + b;
            let (wo, wf) = (w_on[at], w_off[at]);
            out.on.push(wo);
            out.off.push(wf);
            let (m, u) = (m[i] as f64, u[i] as f64);
            if m > 0.0 {
                out.max_weight = out.max_weight.max(wo);
            }
            if u > 0.0 {
                out.max_weight = out.max_weight.max(wf);
            }
            sw[b] += wo * m + wf * u;
            sw2[b] += wo * wo * m + wf * wf * u;
            n[b] += m + u;
        }
    }
    out.effective_share = (0..batches)
        .map(|b| {
            if sw2[b] > 0.0 {
                sw[b] * sw[b] / (n[b] * sw2[b])
            } else {
                1.0
            }
        })
        .collect();
    out
}

/// The fit's final clusters: those of at least `min_cells` cells of the
/// chosen partition, cells of smaller ones moved to their most probable such
/// cluster (its probability their confidence), then numbered 0.. by size,
/// largest first. `raw` and `weighted` are the fitted cells' reads, as read
/// and as the fit saw them.
pub fn settle(
    raw: &Reads,
    weighted: &Weighted,
    part: &Partition,
    min_cells: usize,
    sites: usize,
    strength: f64,
) -> (Partition, Tiles) {
    let sizes = super::summary::sizes(&part.z, part.k);
    // `part` is numbered by size, so the substantial clusters come first.
    let k = sizes.iter().filter(|&&n| n >= min_cells).count().max(1);
    let mut z: Vec<u16> = part
        .z
        .iter()
        .map(|&z| if (z as usize) < k { z } else { OUTSIDE })
        .collect();
    let mut confidence = part.confidence.clone();
    if z.contains(&OUTSIDE) {
        let tiles = Tiles::new(raw, weighted, &z, k, sites, strength);
        for c in 0..z.len() {
            if z[c] == OUTSIDE {
                let (s, m, u) = weighted.of(c);
                let p = tiles.place(s, m, u, &sizes[..k]);
                z[c] = p.cluster;
                confidence[c] = p.probability;
            }
        }
        // Folding the small clusters in can reorder the sizes: renumber.
        super::summary::renumber_by_size(&mut z, k);
    }
    let tiles = Tiles::new(raw, weighted, &z, k, sites, strength);
    (Partition { z, k, confidence }, tiles)
}

/// A cell left out of [`Tiles::new`].
pub const OUTSIDE: u16 = u16::MAX;

/// Fitted cells' reads per site and cluster (sites x clusters, small: a
/// pair of counts per site for each cluster): raw, for the table, and
/// weighted, for the posterior.
pub struct Tiles {
    pub k: usize,
    on: Vec<u32>,
    off: Vec<u32>,
    /// Posterior mean and standard deviation of each site's fraction per
    /// cluster, and the mean's logit.
    mean: Vec<f64>,
    sd: Vec<f64>,
    logit: Vec<f64>,
}

/// Where a cell goes: its most probable cluster, that probability, and the
/// cell's level offset there.
pub struct Placement {
    pub cluster: u16,
    pub probability: f64,
    pub offset: f64,
}

impl Tiles {
    /// `z[c]` is fitted cell `c`'s cluster in `0..k`, or [`OUTSIDE`].
    pub fn new(
        raw: &Reads,
        weighted: &Weighted,
        z: &[u16],
        k: usize,
        sites: usize,
        strength: f64,
    ) -> Self {
        let (mut on, mut off) = (vec![0u32; sites * k], vec![0u32; sites * k]);
        let (mut won, mut woff) = (vec![0f64; sites * k], vec![0f64; sites * k]);
        for (c, &zc) in z.iter().enumerate().filter(|(_, &zc)| zc != OUTSIDE) {
            let (s, m, u) = raw.of(c);
            let (_, wm, wu) = weighted.of(c);
            for i in 0..s.len() {
                let at = s[i] as usize * k + zc as usize;
                on[at] += m[i];
                off[at] += u[i];
                won[at] += wm[i];
                woff[at] += wu[i];
            }
        }
        let (mut mean, mut sd, mut logit) = (Vec::new(), Vec::new(), Vec::new());
        for s in 0..sites {
            let r = s * k..(s + 1) * k;
            let m: f64 = won[r.clone()].iter().sum();
            let n: f64 = m + woff[r.clone()].iter().sum::<f64>();
            let p0 = pooled_fraction(m, n);
            for at in r {
                let x = strength * p0 + won[at];
                let y = strength * (1.0 - p0) + woff[at];
                let t = x + y;
                mean.push(x / t);
                sd.push((x * y / (t * t * (t + 1.0))).sqrt());
                logit.push((x / y).ln());
            }
        }
        Tiles {
            k,
            on,
            off,
            mean,
            sd,
            logit,
        }
    }

    pub fn on(&self, s: usize, k: usize) -> u32 {
        self.on[s * self.k + k]
    }

    pub fn total(&self, s: usize, k: usize) -> u32 {
        self.on[s * self.k + k] + self.off[s * self.k + k]
    }

    /// Posterior mean and standard deviation of site `s`'s fraction in
    /// cluster `k`, from the weighted reads.
    pub fn posterior(&self, s: usize, k: usize) -> (f64, f64) {
        (self.mean[s * self.k + k], self.sd[s * self.k + k])
    }

    /// Where a cell goes: each cluster scored at the cell's own level offset
    /// there (so a cell is placed by its pattern, not its overall level),
    /// clusters weighted by their sizes.
    pub fn place<T: Copy + Into<f64>>(
        &self,
        sites: &[u32],
        on: &[T],
        off: &[T],
        sizes: &[usize],
    ) -> Placement {
        let total: usize = sizes.iter().sum();
        let scored: Vec<(f64, f64)> = (0..self.k)
            .map(|k| {
                let b = self.offset(sites, on, off, k);
                let fit: f64 = (0..sites.len())
                    .map(|i| {
                        let eta = self.logit[sites[i] as usize * self.k + k] + b;
                        on[i].into() * ln_sigmoid(eta) + off[i].into() * ln_sigmoid(-eta)
                    })
                    .sum();
                let prior = (sizes[k].max(1) as f64 / total.max(1) as f64).ln();
                (fit - b * b / 2.0 + prior, b)
            })
            .collect();
        let mx = scored.iter().map(|s| s.0).fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = scored.iter().map(|s| (s.0 - mx).exp()).sum();
        let (k, &(score, offset)) = scored
            .iter()
            .enumerate()
            .max_by(|a, b| a.1 .0.total_cmp(&b.1 .0))
            .expect("a fit has a cluster");
        Placement {
            cluster: k as u16,
            probability: (score - mx).exp() / sum,
            offset,
        }
    }

    /// A cell's level offset in cluster `k`: the shift `b` of every site's
    /// logit fraction that best fits its reads, under a N(0, 1) prior.
    pub fn offset<T: Copy + Into<f64>>(&self, sites: &[u32], on: &[T], off: &[T], k: usize) -> f64 {
        let mut b = 0.0;
        for _ in 0..20 {
            let (mut g, mut h) = (-b, 1.0);
            for i in 0..sites.len() {
                let eta = self.logit[sites[i] as usize * self.k + k] + b;
                let p = 1.0 / (1.0 + (-eta).exp());
                let (m, u) = (on[i].into(), off[i].into());
                let n = m + u;
                g += m - n * p;
                h += n * p * (1.0 - p);
            }
            let step = (g / h).clamp(-1.0, 1.0);
            b += step;
            if step.abs() < 1e-6 {
                break;
            }
        }
        b
    }
}
