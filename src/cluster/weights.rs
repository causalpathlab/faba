//! Importance weights that take the batches out of the reads.

use super::math::pooled_fraction;
use super::read::Reads;

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
