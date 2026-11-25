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
#[derive(clap::ValueEnum, Debug, Clone, Copy, Default, PartialEq, Eq, Hash, serde::Serialize)]
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
#[derive(clap::ValueEnum, Debug, Clone, Copy, Default, PartialEq, Eq, Hash, serde::Serialize)]
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

/// How the view style draws the arms: which is in front, how the p-value
/// reads, and whether its key stays with `--panel-only` (`--key`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ArmLook {
    pub front: Front,
    pub p_axis: PAxis,
    pub key: bool,
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
