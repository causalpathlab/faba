//! Given the chosen partition: each cluster's reads per site, cells outside
//! the fit assigned to a cluster, and every cell's level offset.

use super::read::Reads;
use super::summary::Partition;

/// Pseudo-reads pulling a batch's fraction at a site toward the site's pooled
/// one, so a thinly read site is barely reweighted.
const SHRINK: f64 = 20.0;

/// Largest importance weight a read gets (truncated importance sampling).
const MAX_WEIGHT: f64 = 10.0;

/// Importance weights that take the batches out of the reads.
///
/// The target is reads without batch: at site `s` a read is converted with
/// the site's pooled fraction `p_s`; batch `b` delivers `p_sb` instead (its
/// own fraction there, shrunk toward `p_s`). Each read is weighted by the
/// ratio of the two, `p_s / p_sb` if converted and `(1 - p_s) / (1 - p_sb)`
/// if not, so a batch's weighted fraction at a site is the target's. Weights
/// are capped at [`MAX_WEIGHT`], then scaled per batch and site so its
/// weighted reads add up to its raw reads (self-normalised): every batch
/// brings the evidence it has, no more and no less.
pub struct BatchWeights {
    /// Per entry of the reads, its converted and unconverted reads' weights.
    pub on: Vec<f64>,
    pub off: Vec<f64>,
    /// Per batch, Kish's effective share of its reads, `(sum w)^2 / (n sum w^2)`.
    pub effective_share: Vec<f64>,
    /// Share of reads whose weight was capped.
    pub capped_share: f64,
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
        let p =
            (on[r.clone()].iter().sum::<f64>() + 0.5) / (tot[r.clone()].iter().sum::<f64>() + 1.0);
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
/// chosen partition, numbered 0.. by size; cells of smaller ones go to their
/// most probable such cluster, its probability their confidence.
pub fn settle(
    reads: &Reads,
    part: &Partition,
    min_cells: usize,
    sites: usize,
    strength: f64,
) -> (Partition, Tiles) {
    let mut sizes = vec![0usize; part.k];
    part.z.iter().for_each(|&k| sizes[k as usize] += 1);
    // `part` is numbered by size, so the substantial clusters come first.
    let k = sizes.iter().filter(|&&n| n >= min_cells).count().max(1);
    let kept: Vec<u16> = part.z.iter().map(|&z| z.min(k as u16 - 1)).collect();
    let main: Vec<usize> = (0..part.z.len())
        .filter(|&c| (part.z[c] as usize) < k)
        .collect();
    let main_reads = Reads::from_cells(
        main.iter()
            .map(|&c| {
                let (s, m, u) = reads.of(c);
                (0..s.len()).map(|i| (s[i], m[i], u[i])).collect()
            })
            .collect(),
    );
    let main_z: Vec<u16> = main.iter().map(|&c| part.z[c]).collect();
    let tiles = Tiles::new(&main_reads, &main_z, k, sites, strength);
    let (mut z, mut confidence) = (kept, part.confidence.clone());
    for c in (0..part.z.len()).filter(|&c| part.z[c] as usize >= k) {
        let (s, m, u) = reads.of(c);
        let p = tiles.assign(s, m, u, &sizes[..k]);
        let (best, prob) = p
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(j, &v)| (j as u16, v))
            .unwrap();
        z[c] = best;
        confidence[c] = prob;
    }
    let tiles = Tiles::new(reads, &z, k, sites, strength);
    (Partition { z, k, confidence }, tiles)
}

/// Fitted cells' reads per site and cluster (sites x clusters, small: a
/// pair of counts per site for each cluster).
pub struct Tiles {
    pub k: usize,
    on: Vec<u32>,
    off: Vec<u32>,
    /// The prior `Beta(a_s, b_s)` of every site's fraction.
    a: Vec<f64>,
    b: Vec<f64>,
}

impl Tiles {
    /// `z[c]` is fitted cell `c`'s cluster in `0..k`.
    pub fn new(reads: &Reads, z: &[u16], k: usize, sites: usize, strength: f64) -> Self {
        let (mut on, mut off) = (vec![0u32; sites * k], vec![0u32; sites * k]);
        for (c, &zc) in z.iter().enumerate() {
            let (s, m, u) = reads.of(c);
            for i in 0..s.len() {
                let at = s[i] as usize * k + zc as usize;
                on[at] += m[i];
                off[at] += u[i];
            }
        }
        let (a, b) = (0..sites)
            .map(|s| {
                let (m, n): (u64, u64) = (0..k).fold((0, 0), |(m, n), j| {
                    let at = s * k + j;
                    (m + on[at] as u64, n + (on[at] + off[at]) as u64)
                });
                let p0 = (m as f64 + 0.5) / (n as f64 + 1.0);
                (strength * p0, strength * (1.0 - p0))
            })
            .unzip();
        Tiles { k, on, off, a, b }
    }

    pub fn on(&self, s: usize, k: usize) -> u32 {
        self.on[s * self.k + k]
    }

    pub fn total(&self, s: usize, k: usize) -> u32 {
        self.on[s * self.k + k] + self.off[s * self.k + k]
    }

    /// Posterior mean and standard deviation of site `s`'s fraction in cluster `k`.
    pub fn posterior(&self, s: usize, k: usize) -> (f64, f64) {
        let x = self.a[s] + self.on(s, k) as f64;
        let y = self.b[s] + (self.total(s, k) - self.on(s, k)) as f64;
        let n = x + y;
        (x / n, (x * y / (n * n * (n + 1.0))).sqrt())
    }

    /// A cell's cluster probabilities, each cluster scored at the cell's own
    /// level offset there (so a cell is placed by its pattern, not its
    /// overall level), clusters weighted by their sizes.
    pub fn assign(&self, sites: &[u32], on: &[u32], off: &[u32], sizes: &[usize]) -> Vec<f64> {
        let total: usize = sizes.iter().sum();
        let ll: Vec<f64> = (0..self.k)
            .map(|k| {
                let b = self.offset(sites, on, off, k);
                let fit: f64 = (0..sites.len())
                    .map(|i| {
                        let p = self.posterior(sites[i] as usize, k).0;
                        let eta = (p / (1.0 - p)).ln() + b;
                        // ln sigmoid(eta) and ln sigmoid(-eta), stably.
                        let (lp, lq) = (-(-eta).exp().ln_1p(), -eta.exp().ln_1p());
                        on[i] as f64 * lp + off[i] as f64 * lq
                    })
                    .sum();
                fit - b * b / 2.0 + (sizes[k].max(1) as f64 / total.max(1) as f64).ln()
            })
            .collect();
        let mx = ll.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = ll.iter().map(|l| (l - mx).exp()).collect();
        let sum: f64 = w.iter().sum();
        w.into_iter().map(|x| x / sum).collect()
    }

    /// A cell's level offset in cluster `k`: the shift `b` of every site's
    /// logit fraction that best fits its reads, under a N(0, 1) prior.
    pub fn offset(&self, sites: &[u32], on: &[u32], off: &[u32], k: usize) -> f64 {
        let eta: Vec<f64> = sites
            .iter()
            .map(|&s| {
                let p = self.posterior(s as usize, k).0;
                (p / (1.0 - p)).ln()
            })
            .collect();
        let mut b = 0.0;
        for _ in 0..20 {
            let (mut g, mut h) = (-b, 1.0);
            for i in 0..sites.len() {
                let p = 1.0 / (1.0 + (-(eta[i] + b)).exp());
                let n = (on[i] + off[i]) as f64;
                g += on[i] as f64 - n * p;
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
