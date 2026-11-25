//! One scale per site quantity, chosen once and used by every figure of
//! it: a histogram, a 2D panel's axis and its marginal, `faba qc`'s 2D
//! view and its page. A quantity is drawn on its shown value (a p-value as
//! `-log10 p`, the rest as stored), spread linearly or on `ln(1 + v)`, over
//! one data-beans [`Axis`] and one set of bins.

use super::args::SiteFilterArgs;
use super::pairs::Axis as Quantity;
use super::sites::Criterion;
use data_beans::figure::Axis;

/// How a quantity's shown value is spread along its axes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Scale {
    /// As it is.
    Linear,
    /// On `ln(1 + v)`: zero has a place, long tails fold in.
    Log,
}

impl Scale {
    fn parse(s: &str) -> anyhow::Result<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "linear" => Ok(Scale::Linear),
            "log" => Ok(Scale::Log),
            other => anyhow::bail!("scale {other}: linear or log"),
        }
    }
}

/// Each quantity's [`Scale`], by `Quantity as usize`, and whether
/// `--scale` asked for it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
pub struct Scales([Scale; 7], #[serde(skip)] [bool; 7]);

impl Default for Scales {
    fn default() -> Self {
        Self(Quantity::ALL.map(Quantity::default_scale), [false; 7])
    }
}

impl Scales {
    pub fn of(&self, q: Quantity) -> Scale {
        self.0[q as usize]
    }

    /// Whether `--scale` named `q`'s scale, not left it to its default.
    pub fn asked(&self, q: Quantity) -> bool {
        self.1[q as usize]
    }

    /// The defaults with `items` (`q=linear|log`, as `--scale` takes them)
    /// over them. Log is refused where it cannot apply (signed log odds).
    pub fn parse<S: AsRef<str>>(items: &[S]) -> anyhow::Result<Self> {
        let mut out = Self::default();
        for item in items {
            let item = item.as_ref();
            let Some((q, s)) = item.split_once('=') else {
                anyhow::bail!("--scale {item}: write it as QUANTITY=linear|log, e.g. pv=log");
            };
            let (q, s) = (Quantity::parse(q)?, Scale::parse(s)?);
            anyhow::ensure!(
                s == Scale::Linear || q.can_log(),
                "--scale {item}: {} is signed and has no log scale; use linear",
                q.key()
            );
            out.0[q as usize] = s;
            out.1[q as usize] = true;
        }
        Ok(out)
    }
}

/// Where `v` sits along a line of `scale`: `ln(1 + v)` on a log scale, `v`
/// on a linear one; the transform every [`Ruler`]'s axis reads (a log
/// ruler's is [`Axis::log1p`]), whatever its range. NaN stays NaN.
pub fn along(scale: Scale, v: f64) -> f64 {
    match scale {
        Scale::Log if !v.is_nan() => (v.max(0.0) + 1.0).ln(),
        _ => v,
    }
}

/// Where criterion `c`'s threshold at `raw` is drawn on quantity `q` on
/// `scale`: at its shown value; a count's on a linear scale, whose bars
/// are centred on counts, half a count out, between the last bar it
/// removes and the first it keeps.
pub fn edge(q: Quantity, scale: Scale, c: Criterion, raw: f64) -> f64 {
    let v = q.shown(raw);
    match (scale, q.integer()) {
        (Scale::Linear, true) if c.is_max() => v + 0.5,
        (Scale::Linear, true) => v - 0.5,
        _ => v,
    }
}

/// A quantity's axis and bins, as every figure of it reads them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ruler {
    pub axis: Axis,
    pub bins: usize,
    pub scale: Scale,
}

/// Sites at most the range is read on.
const SAMPLE: usize = 20_000;

/// The 0.1% and 99.9% quantiles of the finite values of `shown` (read on
/// a sample of at most [`SAMPLE`]), where a ruler over them spans; none
/// without a finite value.
pub fn spread(shown: &[f32]) -> Option<(f64, f64)> {
    let step = shown.len().div_ceil(SAMPLE).max(1);
    let v: Vec<f32> = (shown.iter().step_by(step))
        .copied()
        .filter(|x| x.is_finite())
        .collect();
    (!v.is_empty()).then(|| {
        let qs = legume_numeric::matrix::utils::quantiles(&v, &[0.001, 0.999]);
        (qs[0] as f64, qs[1] as f64)
    })
}

impl Ruler {
    /// For quantity `q` on `scale` over its shown values `shown`, reaching
    /// `edges` too, in about `bins` bins: from the values' 0.1% to 99.9%
    /// quantiles (a value past them, an infinity too, is drawn at its end;
    /// see [`spread`] and [`Ruler::spanning`]).
    #[cfg(test)]
    pub fn new(q: Quantity, scale: Scale, shown: &[f32], edges: &[f64], bins: usize) -> Self {
        Self::spanning(q, scale, spread(shown), edges, bins)
    }

    /// For quantity `q` on `scale` over values spread as `spread` (else 0
    /// to 1), reaching `edges` too, in about `bins` bins. An edit ratio
    /// spans 0 to 1, a quantity that cannot be negative starts at 0. A
    /// count on a linear scale gets bins of whole numbers, each centred on
    /// its count.
    pub fn spanning(
        q: Quantity,
        scale: Scale,
        spread: Option<(f64, f64)>,
        edges: &[f64],
        bins: usize,
    ) -> Self {
        let bins = bins.max(1);
        let (mut lo, mut hi) = spread.unwrap_or((0.0, 1.0));
        match q {
            Quantity::EditRatio => (lo, hi) = (0.0, 1.0),
            Quantity::LogOdds => {}
            _ => lo = lo.min(0.0),
        }
        let reach = (edges.iter()).fold(Axis::linear(lo, hi), |a, &e| a.including(e));
        (lo, hi) = (reach.lo, reach.hi);
        if hi <= lo {
            hi = lo + 1.0;
        }
        match (scale, q.integer()) {
            (Scale::Log, _) => Self {
                axis: Axis::log1p(lo.max(0.0), hi),
                bins,
                scale,
            },
            (Scale::Linear, false) => Self {
                axis: Axis::linear(lo, hi),
                bins,
                scale,
            },
            (Scale::Linear, true) => {
                let (lo, hi) = (lo.floor(), hi.ceil());
                let counts = hi - lo + 1.0;
                let width = (counts / bins as f64).ceil().max(1.0);
                let n = (counts / width).ceil() as usize;
                let start = lo - 0.5;
                Self {
                    axis: Axis::linear(start, start + n as f64 * width),
                    bins: n,
                    scale,
                }
            }
        }
    }

    /// [`Ruler::spanning`] `spread` (see [`spread`]) and the thresholds of
    /// `f` on `q` that `on` holds drawn, each at its [`edge`]; with those
    /// edges, ascending. Every figure of a quantity, a panel or `faba qc`'s
    /// 2D view, builds its ruler so.
    pub fn cut(
        (q, scale): (Quantity, Scale),
        spread: Option<(f64, f64)>,
        f: &SiteFilterArgs,
        on: impl Fn(Criterion) -> bool,
        bins: usize,
    ) -> (Self, Vec<f64>) {
        let mut edges: Vec<f64> = (q.criteria().iter().copied())
            .filter(|&c| on(c))
            .map(|c| edge(q, scale, c, c.get(f)))
            .collect();
        edges.sort_by(f64::total_cmp);
        edges.dedup();
        (Self::spanning(q, scale, spread, &edges, bins), edges)
    }

    /// [`Ruler::clamp`] on a copy of `shown`.
    #[cfg(test)]
    pub fn clamped(&self, shown: &[f32]) -> Vec<f32> {
        let mut v = shown.to_vec();
        self.clamp(&mut v);
        v
    }

    /// `shown` with infinities at the axis' ends, so they are counted in
    /// its end bins; NaN (no value) stays, and is left out.
    pub fn clamp(&self, shown: &mut [f32]) {
        let (lo, hi) = (self.axis.lo as f32, self.axis.hi as f32);
        for v in shown {
            match *v {
                f32::INFINITY => *v = hi,
                f32::NEG_INFINITY => *v = lo,
                _ => {}
            }
        }
    }

    /// Where `v` sits along this ruler's line (see [`along`]).
    pub fn along(&self, v: f64) -> f64 {
        along(self.scale, v)
    }

    /// The axis' ends along its line.
    pub fn span(&self) -> (f64, f64) {
        (self.along(self.axis.lo), self.along(self.axis.hi))
    }

    /// The axis' ticks along its line, each with its label.
    pub fn ticks(&self) -> Vec<(f64, String)> {
        (self.axis.ticks().into_iter())
            .map(|t| (self.along(t), data_beans::figure::short(t)))
            .collect()
    }
}

#[cfg(test)]
#[path = "tests/scale.rs"]
mod tests;
