//! The view style's two arms: the signal arm (WT) and the control arm
//! (MUT) of a site table, each site's own reads on each, so a histogram
//! or a 2D panel draws the control's sites behind the signal's (or in
//! front, `--front mut`). Only the reads and the edit ratio exist on each
//! arm; the p-value, the fold and the log odds compare the arms, and the
//! cells are counted on the signal's `_site` matrices only, so those are
//! drawn for WT alone. Nothing here decides a site: the thresholds and
//! their checks are [`super::sites::Criterion`]'s, read off the fail bits.

use super::pairs::{Axis, Sites};
use super::sites::SiteTable;

/// The control arm's bars and squares: light gray.
pub const MUT_FILL: &str = "#c9c9c9";
/// A signal site some threshold filters out: white, outlined.
pub const FILTERED_FILL: &str = "#ffffff";
/// The outline of what is drawn white.
pub const OUTLINE: &str = crate::figure::INK;
/// A threshold's line, in either palette.
pub const THRESHOLD: &str = "#e2461c";

/// Which arm is drawn in front (`--front`, `f` in the view).
#[derive(
    clap::ValueEnum,
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Front {
    /// The signal arm over the control's.
    #[default]
    Wt,
    /// The control arm over the signal's.
    Mut,
}

impl Front {
    /// The other arm in front.
    pub fn next(self) -> Self {
        match self {
            Front::Wt => Front::Mut,
            Front::Mut => Front::Wt,
        }
    }

    /// The arm's name.
    pub fn name(self) -> &'static str {
        match self {
            Front::Wt => "WT",
            Front::Mut => "MUT",
        }
    }
}

/// How the p-value's ticks read (`--p-axis`, `u` in the view): the bars
/// sit on `-log10 p` either way.
#[derive(
    clap::ValueEnum,
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum PAxis {
    /// The p-value itself: 1, 0.2, 1e-3.
    #[default]
    P,
    /// `-log10 p`: 0, 0.7, 3.
    Neglog10,
}

impl PAxis {
    /// The other reading.
    pub fn next(self) -> Self {
        match self {
            PAxis::P => PAxis::Neglog10,
            PAxis::Neglog10 => PAxis::P,
        }
    }

    /// What its ticks read, as an axis title.
    pub fn name(self) -> &'static str {
        match self {
            PAxis::P => "p-value",
            PAxis::Neglog10 => "-log10 p",
        }
    }

    /// The tick text at `v` along `-log10 p`.
    pub fn label(self, v: f64) -> String {
        match self {
            PAxis::P => fmt_p(10f64.powf(-v)),
            PAxis::Neglog10 => data_beans::figure::short(v),
        }
    }
}

/// How a count is spread along an axis or a shade (`--yscale`,
/// `--marginal-scale`, `--count-scale`): as it is, on its square root, or
/// on `ln(1 + n)`, so a tall background and a small series both read.
#[derive(
    clap::ValueEnum,
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum CountScale {
    Linear,
    #[default]
    Sqrt,
    Log1p,
}

impl CountScale {
    /// The next a key steps to: linear, sqrt, log1p, linear.
    pub fn next(self) -> Self {
        match self {
            CountScale::Linear => CountScale::Sqrt,
            CountScale::Sqrt => CountScale::Log1p,
            CountScale::Log1p => CountScale::Linear,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            CountScale::Linear => "linear",
            CountScale::Sqrt => "sqrt",
            CountScale::Log1p => "log1p",
        }
    }

    /// `n` along the scale.
    pub fn of(self, n: f64) -> f64 {
        let n = n.max(0.0);
        match self {
            CountScale::Linear => n,
            CountScale::Sqrt => n.sqrt(),
            CountScale::Log1p => n.ln_1p(),
        }
    }

    /// The count at `t` along the scale.
    pub fn count(self, t: f64) -> f64 {
        match self {
            CountScale::Linear => t,
            CountScale::Sqrt => t * t,
            CountScale::Log1p => t.exp_m1(),
        }
    }

    /// A few round counts to tick an axis of counts reaching `most` on: 0,
    /// the round count nearest its middle, and the least round count at or
    /// above `most`, where the axis ends; each where it sits
    /// ([`CountScale::of`]).
    pub fn ticks(self, most: f64) -> Vec<f64> {
        let mut out = vec![0.0];
        if most > 0.0 {
            let top = round_count(most, true);
            let mid = round_count(self.count(self.of(top) / 2.0), false);
            for v in [mid, top] {
                if v > out[out.len() - 1] {
                    out.push(v);
                }
            }
        }
        out
    }
}

/// `v` rounded to a round count (1, 2, 2.5 or 5 times a power of ten):
/// the nearest, or the least at or above `v` when `up`.
fn round_count(v: f64, up: bool) -> f64 {
    if v < 1.0 {
        return if up { 1.0 } else { v.round() };
    }
    let p = 10f64.powf(v.log10().floor());
    let steps = [1.0, 2.0, 2.5, 5.0, 10.0].map(|m| m * p);
    match up {
        true => steps.into_iter().find(|&s| s >= v).unwrap_or(10.0 * p),
        false => (steps.into_iter())
            .min_by(|a, b| (a - v).abs().total_cmp(&(b - v).abs()))
            .unwrap_or(p),
    }
}

/// Each quantity's histogram y scale (`--yscale`, `y` in the view), by
/// [`Axis`] `as usize`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct YScales([CountScale; 7]);

impl YScales {
    pub fn of(&self, q: Axis) -> CountScale {
        self.0[q as usize]
    }

    pub fn set(&mut self, q: Axis, s: CountScale) {
        self.0[q as usize] = s;
    }

    /// `--yscale`: `q=SCALE` for one quantity, or `SCALE` for every one,
    /// later over earlier, over the default (sqrt).
    pub fn parse(items: &[String]) -> anyhow::Result<Self> {
        use clap::ValueEnum;
        let scale = |s: &str| {
            CountScale::from_str(s.trim(), true)
                .map_err(|_| anyhow::anyhow!("--yscale {s}: linear, sqrt or log1p"))
        };
        let mut out = Self::default();
        for item in items {
            match item.split_once('=') {
                Some((q, s)) => {
                    let q = Axis::parse(q).map_err(|e| anyhow::anyhow!("--yscale {item}: {e}"))?;
                    out.set(q, scale(s)?);
                }
                None => out = Self([scale(item)?; 7]),
            }
        }
        Ok(out)
    }
}

/// Where the view style's key goes (`--key-at`).
#[derive(
    clap::ValueEnum,
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum KeyAt {
    /// Beside the plot, at its right: the plot narrows by the key's width;
    /// a key that would leave the plot under 60% of the panel's width is
    /// drawn smaller, else at the top.
    #[default]
    Right,
    /// Above the plot, its series in a row (two, or shorter names, where
    /// one does not fit): the plot keeps its width, losing that row's
    /// height.
    Top,
    /// Over the plot, in its top right corner: the plot keeps its width.
    Inside,
    /// In a figure of its own, `{stem}_key`: none in the panel.
    Alone,
}

/// How the view style draws the arms: which is in front, how the p-value
/// reads, whether its key stays with `--panel-only` (`--key`) and where,
/// and how each count is spread: a histogram's y, a marginal's depth, a
/// 2D square's shade.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArmLook {
    pub front: Front,
    pub p_axis: PAxis,
    pub key: bool,
    pub key_at: KeyAt,
    pub y: YScales,
    pub marginal: CountScale,
    pub count: CountScale,
}

impl Default for ArmLook {
    fn default() -> Self {
        Self {
            front: Front::default(),
            p_axis: PAxis::default(),
            key: false,
            key_at: KeyAt::default(),
            y: YScales::default(),
            // One busy bar (p near 1, an edit ratio near 0) leaves the
            // others a depth that reads.
            marginal: CountScale::Log1p,
            // A square's shade as it always was.
            count: CountScale::Log1p,
        }
    }
}

/// `p` in two significant figures, as a p-value reads: 1, 0.2, 0.043,
/// then 3.2e-4 below 0.01.
pub fn fmt_p(p: f64) -> String {
    let trim = |s: String| match s.contains('.') {
        true => s.trim_end_matches('0').trim_end_matches('.').to_string(),
        false => s,
    };
    if !p.is_finite() || p <= 0.0 {
        return "0".into();
    }
    // As read in two figures: a threshold of 0.01 kept in single precision
    // reads 0.01, not 1e-2.
    let p: f64 = format!("{p:.1e}").parse().unwrap_or(p);
    if p >= 0.01 {
        let decimals = (1 - p.log10().floor() as i32).max(0) as usize;
        return trim(format!("{p:.decimals$}"));
    }
    let s = format!("{p:.1e}");
    match s.split_once('e') {
        Some((m, e)) => format!("{}e{e}", trim(m.to_string())),
        None => s,
    }
}

/// The round `-log10 p` a p-value axis caps at: the greatest of 10, 20,
/// 30, 50, 100, 200, 300 at most `hi` (the largest shown), so its last bar
/// or square holds every p at or below it; raised past every threshold in
/// `cuts`. None when `hi` is under the least, the axis then its own.
pub fn p_cap(hi: f64, cuts: impl Iterator<Item = f64>) -> Option<f64> {
    const LADDER: [f64; 7] = [10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0];
    let mut cap = LADDER.into_iter().rev().find(|&k| k <= hi)?;
    for d in cuts {
        if d >= cap {
            cap = LADDER.into_iter().find(|&k| k > d)?;
        }
    }
    Some(cap)
}

/// The p-value's cap in every panel of `s` at the thresholds `f`, the
/// histogram's and the 2D panels' alike: [`p_cap`] of where the p-values
/// the cut reads spread (see [`super::scale::spread`]), past the p-value
/// threshold when on.
pub fn p_cap_of(s: &Sites, f: &super::args::SiteFilterArgs) -> Option<f64> {
    let hi = s.column(Axis::Pvalue).spread?.1;
    let on = (Axis::Pvalue.criteria().iter())
        .filter(|c| s.criteria.contains(c) && !c.is_off(f))
        .map(|c| Axis::Pvalue.shown(c.get(f)))
        .filter(|d| d.is_finite());
    p_cap(hi, on)
}

/// Whether `q` has a value on each arm of `t`: the reads and the edit
/// ratio, where the table has a control arm.
pub fn per_arm(q: Axis, t: &SiteTable) -> bool {
    t.has_control() && matches!(q, Axis::Coverage | Axis::Converted | Axis::EditRatio)
}

/// Site `i`'s value of `q` on the signal arm, as shown: its own reads (not
/// the control's added) where `q` has a value on each arm, else `q`'s.
pub fn wt(s: &Sites, q: Axis, i: usize) -> f64 {
    let t = s.table;
    match q {
        _ if !per_arm(q, t) => q.value(s, i),
        Axis::Coverage => t.coverage[i] as f64,
        Axis::Converted => t.converted[i] as f64,
        _ => t.edit_ratio[i] as f64,
    }
}

/// Site `i`'s value of `q` on the control arm, as shown: NaN where `q`
/// has none there (an edit ratio over no control reads, a quantity that
/// compares the arms).
pub fn mutant(s: &Sites, q: Axis, i: usize) -> f64 {
    let t = s.table;
    match q {
        _ if !per_arm(q, t) => f64::NAN,
        Axis::Coverage => t.control_coverage[i] as f64,
        Axis::Converted => t.control_converted[i] as f64,
        _ => match t.control_coverage[i] {
            0 => f64::NAN,
            n => (t.control_converted[i] as f32 / n as f32) as f64,
        },
    }
}

/// Whether a knob's histogram draws `q` on each arm of `t`: where
/// [`per_arm`] has it, and the p-value (the control arm tested alone, see
/// [`SiteTable::control_pvalues`]) and the cells (counted on the control
/// arm's own batches) where the table has them. The signal arm is drawn
/// on the values its thresholds read ([`wt`]). The 2D panels keep to
/// [`per_arm`].
pub fn per_arm_hist(q: Axis, t: &SiteTable) -> bool {
    per_arm(q, t)
        || match q {
            Axis::Pvalue => t.has_control(),
            Axis::Cells => t.control_cells.is_some(),
            _ => false,
        }
}

/// Site `i`'s value of `q` on the control arm as a knob's histogram shows
/// it: its own p-value against sequencing error, its own batches' cells,
/// else [`mutant`]'s.
pub fn mutant_hist(s: &Sites, q: Axis, i: usize) -> f64 {
    let t = s.table;
    match q {
        Axis::Pvalue if per_arm_hist(q, t) => shown_p(t.control_pvalues()[i]),
        Axis::Cells => t.control_cells.as_ref().map_or(f64::NAN, |c| c[i] as f64),
        _ => mutant(s, q, i),
    }
}

/// A p-value as shown: `-log10 p`, infinite at 0.
fn shown_p(p: f32) -> f64 {
    match p as f64 {
        p if p <= 0.0 => f64::INFINITY,
        p => (-p.log10()).max(0.0),
    }
}

/// One quantity's values on each arm, along its line (see
/// [`super::scale::along`]), for the 2D panels; and where they spread
/// together.
pub struct ArmColumn {
    pub wt: Vec<f32>,
    pub mutant: Vec<f32>,
    pub spread: Option<(f64, f64)>,
}

impl ArmColumn {
    /// Quantity `q` of `s` on each arm, on the scale the cache holds for
    /// it.
    pub fn of(s: &Sites, q: Axis) -> Self {
        use rayon::prelude::*;
        let scale = s.cache.scales.of(q);
        let read = |value: fn(&Sites, Axis, usize) -> f64| -> Vec<f32> {
            (0..s.fails.len())
                .into_par_iter()
                .map(|i| value(s, q, i) as f32)
                .collect()
        };
        let (wt_shown, mut_shown) = (read(wt), read(mutant));
        let both: Vec<f32> = wt_shown.iter().chain(&mut_shown).copied().collect();
        let along = |v: &[f32]| -> Vec<f32> {
            (v.par_iter())
                .map(|&x| super::scale::along(scale, x as f64) as f32)
                .collect()
        };
        Self {
            wt: along(&wt_shown),
            mutant: along(&mut_shown),
            spread: super::scale::spread(&both),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A p-value reads in two significant figures, small ones in
    /// scientific notation.
    #[test]
    fn a_p_value_reads_in_two_figures() {
        for (p, want) in [
            (1.0, "1"),
            (0.1995, "0.2"),
            (0.05, "0.05"),
            (0.0398, "0.04"),
            (0.0427, "0.043"),
            (1e-3, "1e-3"),
            (3.16e-4, "3.2e-4"),
            (1e-30, "1e-30"),
        ] {
            assert_eq!(fmt_p(p), want, "{p}");
        }
        assert_eq!(PAxis::P.label(2.0), "0.01");
        assert_eq!(PAxis::Neglog10.label(2.0), "2");
    }
}
