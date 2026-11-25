//! A knob's histogram as the view draws it, on screen, in its figure and
//! as a `--style view` panel: the control arm's sites (MUT) light gray
//! behind, the signal arm's (WT) in front, white and outlined where some
//! threshold filters the site out and dark where every one keeps it (or
//! the control arm in front, `--front mut`). The bars are the knob's own
//! bins, each threshold's bin split at it, so its dashed line sits on a
//! cut between bars with the words for either side beside it. Every count
//! is read off the fail bits: nothing here decides a site.

use super::*;
use crate::figure::{set_width, text_width, FAINT};
use crate::qc::arms::{self, ArmLook, CountScale, Front, KeyAt, PAxis};
use crate::qc::site_plot::Titles;
use ratatui::buffer::Buffer;

/// The key's names of the series: the control arm's, the signal arm's
/// filtered out by any threshold, and its retained by every one.
const KEY_MUT: &str = "background (MUT)";
/// The control arm's name on the p-value, drawn on its own test against
/// sequencing error, not the threshold's WT-vs-MUT test.
const KEY_MUT_P: &str = "background (MUT, vs error)";
pub(super) const KEY_FILTERED: &str = "filtered out (WT)";
pub(super) const KEY_KEPT: &str = "retained (WT)";
/// The series' short names, for a key at the top too narrow for the long.
const SHORT_MUT: &str = "MUT";
const SHORT_FILTERED: &str = "WT filtered";
const SHORT_KEPT: &str = "WT kept";

/// Rows of words a threshold's line takes over the bars, in points of the
/// 8-point layout.
const WORDS_ROW: f64 = 9.0;
/// The words' size, in points of the 8-point layout.
const WORDS: f64 = 7.0;

/// One knob's histogram on both arms: its bars' edges as shown, ascending,
/// and what each bar holds.
pub(super) struct ArmHist {
    pub(super) c: Criterion,
    /// Bar `j` holds the values from `edges[j]` up to `edges[j + 1]`
    /// (that edge itself in the next bar).
    pub(super) edges: Vec<f64>,
    /// Each edge along the drawn line (the bins' scale), for placing it.
    pub(super) at: Vec<f64>,
    /// Whether an edge is a threshold's alone, named by its words, not a
    /// tick.
    split: Vec<bool>,
    /// The signal arm's sites each bar holds: kept by every threshold, and
    /// filtered out by any.
    pub(super) kept: Vec<usize>,
    pub(super) filtered: Vec<usize>,
    /// The control arm's, where the quantity has a value on it.
    pub(super) mutant: Option<Vec<usize>>,
    /// Each threshold drawn: its criterion, its value and the edge it is
    /// at.
    pub(super) cuts: Vec<(Criterion, f64, usize)>,
    /// Where infinities were drawn, as shown: at the high end, the low.
    pub(super) clamped: (Option<f64>, Option<f64>),
    /// Whether the table has a control arm at all.
    has_control: bool,
}

/// Position of `v` along a line of `scale`, as its bins lay it out.
fn line(scale: Scale, v: f64) -> f64 {
    match scale {
        Scale::Linear => v,
        Scale::Sqrt => v.max(0.0).sqrt(),
        Scale::Log => v.max(0.0).ln_1p(),
    }
}

/// Where bin `k` of `bins` starts: the least value it holds, found by
/// halving between values below and in it; the values span `lo..=hi`. On
/// a log or square-root line, which starts at 0, none below it.
fn real_edge(bins: &Binning, k: i32, (lo, hi): (f64, f64)) -> f64 {
    if bins.scale != Scale::Linear && bins.key(0.0) >= k {
        return 0.0;
    }
    let span = (hi - lo).abs().max(1.0);
    let (mut a, mut b, mut step) = (lo, hi, span);
    while bins.key(a) >= k {
        a -= step;
        step *= 2.0;
        if !a.is_finite() {
            return a;
        }
    }
    step = span;
    while bins.key(b) < k {
        b += step;
        step *= 2.0;
        if !b.is_finite() {
            return b;
        }
    }
    loop {
        let m = a + (b - a) / 2.0;
        if m <= a || m >= b {
            return b;
        }
        match bins.key(m) >= k {
            true => b = m,
            false => a = m,
        }
    }
}

/// The round `-log10 p` a p-value's ticks may sit at past 0.05.
const ROUND_NEG_LOG10_P: [f64; 15] = [
    2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 150.0, 200.0, 300.0,
];

/// Where threshold `c` at `d` (as shown) splits a bar: the values passing
/// it on one side, failing it on the other. An upper bound read as it is
/// (the edit ratio's) keeps the value at it, so splits just past it.
fn split_at(c: Criterion, d: f64) -> f64 {
    match c.is_max() && !c.neg_log() {
        true => d.next_up(),
        false => d,
    }
}

impl ArmHist {
    /// Knob `c` of `s` at the thresholds `f` on the bins `bins` (the knob's
    /// column's), reaching the values of both arms and the thresholds.
    pub(super) fn new(s: &Sites, c: Criterion, f: &SiteFilterArgs, bins: Binning) -> Self {
        use rayon::prelude::*;
        let q = Axis::of(c);
        let n = s.fails.len();
        let wt: Vec<f64> = (0..n).into_par_iter().map(|i| arms::wt(s, q, i)).collect();
        let mutant: Option<Vec<f64>> = (arms::per_arm_hist(q, s.table)).then(|| {
            (0..n)
                .into_par_iter()
                .map(|i| arms::mutant_hist(s, q, i))
                .collect()
        });
        let values = || wt.iter().chain(mutant.iter().flatten()).copied();
        let (mut lo, mut hi) = (values().filter(|v| v.is_finite()))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| {
                (l.min(v), h.max(v))
            });
        if lo > hi {
            (lo, hi) = (0.0, 0.0);
        }
        // The thresholds drawn: those of the quantity that apply and are on,
        // the axis reaching each.
        let on: Vec<(Criterion, f64)> = (q.criteria().iter().copied())
            .filter(|k| s.criteria.contains(k) && !k.is_off(f))
            .map(|k| (k, k.get(f)))
            .filter(|&(k, raw)| k.display(raw).is_finite())
            .collect();
        for &(k, raw) in &on {
            (lo, hi) = (lo.min(k.display(raw)), hi.max(k.display(raw)));
        }
        // A p-value's last bar holds every p at or below a round cap, the
        // 2D panels' too.
        let cap = (c == Criterion::MaxPv)
            .then(|| arms::p_cap_of(s, f))
            .flatten();
        if let Some(k) = cap {
            hi = k;
        }
        let clamped = (
            values().any(|v| v == f64::INFINITY || v > hi).then_some(hi),
            values().any(|v| v == f64::NEG_INFINITY).then_some(lo),
        );
        let kmin = bins.key(lo);
        let kmax = bins.key(hi).min(kmin + MAX_BINS - 1);
        let mut edges: Vec<f64> = (kmin..=kmax + 1)
            .map(|k| match c.is_integer() {
                true => bins.lower_edge(k) as f64,
                false => real_edge(&bins, k, (lo, hi)),
            })
            .collect();
        // Whole counts: a bin holding none is no bar.
        edges.dedup();
        // No edge past the cap: the last bar ends at it.
        if cap.is_some() {
            while edges.len() > 2 && edges[edges.len() - 2] >= hi {
                edges.pop();
            }
            let last = edges.len() - 1;
            edges[last] = hi;
        }
        let mut split = vec![false; edges.len()];
        let last = |e: &[f64]| e[e.len() - 1];
        for &(k, raw) in &on {
            let e = split_at(k, k.display(raw)).clamp(edges[0], last(&edges));
            if let Err(at) = edges.binary_search_by(|x| x.total_cmp(&e)) {
                edges.insert(at, e);
                split.insert(at, true);
            }
        }
        let cuts = (on.iter())
            .map(|&(k, raw)| {
                let e = split_at(k, k.display(raw)).clamp(edges[0], last(&edges));
                let at = edges.partition_point(|&x| x < e);
                (k, raw, at)
            })
            .collect();
        let bars = edges.len() - 1;
        let bar = |v: f64| -> Option<usize> {
            let v = match v {
                f64::INFINITY => hi,
                f64::NEG_INFINITY => lo,
                v if v.is_nan() => return None,
                v => v,
            };
            Some((edges.partition_point(|&e| e <= v).saturating_sub(1)).min(bars - 1))
        };
        let (mut kept, mut filtered) = (vec![0; bars], vec![0; bars]);
        for (i, &v) in wt.iter().enumerate() {
            if let Some(b) = bar(v) {
                match s.fails[i] == 0 {
                    true => kept[b] += 1,
                    false => filtered[b] += 1,
                }
            }
        }
        let mutant = mutant.map(|m| {
            let mut out = vec![0; bars];
            m.iter().filter_map(|&v| bar(v)).for_each(|b| out[b] += 1);
            out
        });
        let at = edges.iter().map(|&e| line(bins.scale, e)).collect();
        Self {
            c,
            edges,
            at,
            split,
            kept,
            filtered,
            mutant,
            cuts,
            clamped,
            has_control: s.table.has_control(),
        }
    }

    /// Bars it draws.
    pub(super) fn bars(&self) -> usize {
        self.kept.len()
    }

    /// The signal arm's sites in bar `j`.
    pub(super) fn total(&self, j: usize) -> usize {
        self.kept[j] + self.filtered[j]
    }

    /// The control arm's sites in bar `j`, none without it.
    fn mutant_at(&self, j: usize) -> usize {
        self.mutant.as_ref().map_or(0, |m| m[j])
    }

    /// `v` (as shown) as a tick reads it: a whole count, a p-value as
    /// `p_axis` says, or in three figures.
    pub(super) fn tick_text(&self, v: f64, p_axis: PAxis) -> String {
        match self.c {
            c if c.is_integer() => format!("{}", v.round() as i64),
            Criterion::MaxPv => p_axis.label(v),
            c => c.fmt_short(v),
        }
    }

    /// The ticks it could take, along its line: every edge but a
    /// threshold's own, then where infinities were drawn, said so ("<= p"
    /// past which a p-value reads 0); the second flag says it is wanted
    /// first.
    pub(super) fn tick_choices(&self, scale: Scale, p_axis: PAxis) -> Vec<(f64, String, bool)> {
        let mut out: Vec<(f64, String, bool)> = match self.c {
            Criterion::MaxPv => self.round_p_ticks(scale, p_axis),
            _ => (self.edges.iter().zip(&self.at))
                .zip(&self.split)
                .filter(|(_, &split)| !split)
                .map(|((&e, &at), _)| (at, self.tick_text(e, p_axis), false))
                .collect(),
        };
        if let Some(hi) = self.clamped.0 {
            let said = match (self.c, p_axis) {
                (Criterion::MaxPv, PAxis::P) => format!("<= {}", arms::fmt_p(10f64.powf(-hi))),
                _ => format!(">= {}", self.tick_text(hi, p_axis)),
            };
            out.push((line(scale, hi), said, true));
        }
        if let Some(lo) = self.clamped.1 {
            out.push((
                line(scale, lo),
                format!("<= {}", self.tick_text(lo, p_axis)),
                true,
            ));
        }
        out
    }

    /// A p-value's ticks at round p (1, 0.5, 0.2, 0.1, 0.05, 0.01, then
    /// round powers of ten: 1e-3, 1e-5, 1e-10, 1e-20 ...), each placed
    /// where that p sits along the line, within the drawn range; the cap,
    /// if any, is said apart.
    fn round_p_ticks(&self, scale: Scale, p_axis: PAxis) -> Vec<(f64, String, bool)> {
        let round = [0.0, 2f64.log10(), 5f64.log10(), 1.0, 20f64.log10()]
            .into_iter()
            .chain(ROUND_NEG_LOG10_P);
        let (a, b) = (self.edges[0], self.edges[self.edges.len() - 1]);
        let end = self.clamped.0;
        round
            .filter(|&d| d >= a && d <= b && Some(d) != end)
            .map(|d| (line(scale, d), self.tick_text(d, p_axis), false))
            .collect()
    }

    /// The words either side of threshold `c` at `raw`: the low end of
    /// the axis first.
    pub(super) fn words(&self, c: Criterion, raw: f64) -> (String, String) {
        let v = match c {
            Criterion::MaxPv => arms::fmt_p(raw),
            c => c.fmt_short(raw),
        };
        let per_arm = self.mutant.is_some();
        match c {
            Criterion::MaxPv => (format!("p > {v}"), format!("p <= {v}")),
            Criterion::MaxEditRatio => (format!("<= {v}"), format!("> {v}")),
            Criterion::MinCells => (format!("< {v} cells"), format!(">= {v} cells")),
            // Its threshold is on both arms' reads added.
            Criterion::MinCoverage if per_arm => {
                (format!("WT+MUT < {v}"), format!("WT+MUT >= {v}"))
            }
            Criterion::MinCoverage | Criterion::MinConverted => {
                (format!("< {v} reads"), format!(">= {v} reads"))
            }
            _ => (format!("< {v}"), format!(">= {v}")),
        }
    }

    /// The series its key lists, the control arm (the background behind)
    /// first, then the signal arm's filtered out and its retained; then a
    /// note for a quantity with no control arm drawn.
    fn legend(&self) -> (Vec<(Series, &'static str)>, Option<&'static str>) {
        let mut items: Vec<(Series, &str)> = Vec::new();
        if self.mutant.is_some() {
            items.push((Series::Mutant, self.mutant_name()));
        }
        items.extend([(Series::Filtered, KEY_FILTERED), (Series::Kept, KEY_KEPT)]);
        let note = match (self.mutant.is_some(), self.has_control) {
            (true, _) => None,
            (false, true) => Some("WT only"),
            (false, false) => Some("no MUT arm"),
        };
        (items, note)
    }

    /// A key at the top's rows within `room` across at `words` points of
    /// type, each item its series (none for the note) and name: the long
    /// names in one row, else the short, else the long in two rows, else
    /// the short in two, else the short as many rows as they take.
    fn top_key(&self, words: f64, room: f64) -> Vec<Vec<(Option<Series>, &'static str)>> {
        let (items, note) = self.legend();
        let short = |s: Series| match s {
            Series::Mutant => SHORT_MUT,
            Series::Filtered => SHORT_FILTERED,
            Series::Kept => SHORT_KEPT,
        };
        let named = |short_names: bool| -> Vec<(Option<Series>, &'static str)> {
            (items.iter())
                .map(|&(s, long)| (Some(s), if short_names { short(s) } else { long }))
                .chain(note.map(|n| (None, n)))
                .collect()
        };
        // Items in rows, each row taking what fits after the last.
        let wrap = |all: Vec<(Option<Series>, &'static str)>| {
            let mut rows: Vec<Vec<(Option<Series>, &'static str)>> = Vec::new();
            for item in all {
                match rows.last_mut() {
                    Some(row)
                        if top_row_width(row, words) + KEY_ITEM_GAP + item_width(item, words)
                            <= room =>
                    {
                        row.push(item)
                    }
                    _ => rows.push(vec![item]),
                }
            }
            rows
        };
        [(false, 1), (true, 1), (false, 2), (true, 2)]
            .into_iter()
            .map(|(short_names, most)| (wrap(named(short_names)), most))
            .find(|(rows, most)| rows.len() <= *most)
            .map_or_else(|| wrap(named(true)), |(rows, _)| rows)
    }

    /// The control arm's name in the key: on the p-value, the test it is
    /// drawn on.
    pub(super) fn mutant_name(&self) -> &'static str {
        match self.c {
            Criterion::MaxPv => KEY_MUT_P,
            _ => KEY_MUT,
        }
    }

    /// The series in the order drawn: the arm behind first.
    fn order(&self, front: Front) -> Vec<Series> {
        let wt = [Series::Filtered, Series::Kept];
        match (self.mutant.is_some(), front) {
            (false, _) => wt.to_vec(),
            (true, Front::Wt) => [Series::Mutant, wt[0], wt[1]].to_vec(),
            (true, Front::Mut) => [wt[0], wt[1], Series::Mutant].to_vec(),
        }
    }

    /// Bar `j`'s height in `series`: a filtered bar stands on the kept,
    /// its top the signal arm's sites.
    fn value(&self, series: Series, j: usize) -> usize {
        match series {
            Series::Mutant => self.mutant_at(j),
            Series::Filtered => self.total(j),
            Series::Kept => self.kept[j],
        }
    }

    /// The tallest bar of either arm.
    pub(super) fn tallest(&self) -> usize {
        (0..self.bars())
            .map(|j| self.total(j).max(self.mutant_at(j)))
            .max()
            .unwrap_or(0)
    }
}

/// One of a histogram's series, as drawn and keyed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Series {
    Mutant,
    Filtered,
    Kept,
}

impl Series {
    /// Drawn in the box `b`: the palette's kept bar colour `kept`.
    fn draw(self, c: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64), kept: &str) {
        match self {
            Series::Mutant => c.rect(x, y, w, h, arms::MUT_FILL),
            Series::Kept => c.rect(x, y, w, h, kept),
            Series::Filtered => c.outlined((x, y, w, h), arms::FILTERED_FILL, (arms::OUTLINE, 0.4)),
        }
    }
}

/// A threshold's words as placed: where they start, on which row over
/// the bars.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Placed {
    pub(super) text: String,
    pub(super) from: f64,
    pub(super) row: usize,
}

/// Room left between a threshold's line and its words, either side.
const WORDS_GAP: f64 = 3.0;

impl ArmHist {
    /// Each threshold's words placed over the bars, inside `(left, right)`:
    /// each side's beside the line on that side, on the first row where
    /// they are clear of words placed before; a side with no room for its
    /// words has none, never the other side's place. `x_of` places a
    /// position along the line, `width` measures a text.
    pub(super) fn place_words(
        &self,
        x_of: &dyn Fn(f64) -> f64,
        (left, right): (f64, f64),
        width: impl Fn(&str) -> f64,
    ) -> Vec<Placed> {
        let mut cuts = self.cuts.clone();
        cuts.sort_by_key(|cut| cut.2);
        let mut placed: Vec<Placed> = Vec::new();
        let mut put = |text: String, from: f64| {
            let to = from + width(&text);
            let clear = |row: usize| {
                (placed.iter().filter(|p| p.row == row))
                    .all(|p| to + 2.0 <= p.from || from >= p.from + width(&p.text) + 2.0)
            };
            let row = (0..).find(|&r| clear(r)).expect("some row is clear");
            placed.push(Placed { text, from, row });
        };
        for &(c, raw, e) in &cuts {
            let lx = x_of(self.at[e]);
            let (low, high) = self.words(c, raw);
            let (wl, wh) = (width(&low), width(&high));
            if lx - WORDS_GAP - wl >= left {
                put(low, lx - WORDS_GAP - wl);
            }
            if lx + WORDS_GAP + wh <= right {
                put(high, lx + WORDS_GAP);
            }
        }
        placed
    }
}

/// How a knob's histogram is drawn beyond its bars: its y scale, its arms,
/// its titles.
pub(super) struct HistLook {
    pub(super) y_scale: CountScale,
    /// The scale its bins are on.
    pub(super) scale: Scale,
    pub(super) look: ArmLook,
    pub(super) title: String,
    pub(super) x_title: String,
    pub(super) y_title: String,
}

/// A key item's width at `words` points of type: its swatch and the room
/// after it, then its name (a note has no swatch).
fn item_width((series, name): (Option<Series>, &str), words: f64) -> f64 {
    let swatch = if series.is_some() {
        KEY_SWATCH_ROOM
    } else {
        0.0
    };
    swatch + set_width(name, words)
}

/// A row of a key at the top: its items and the gaps between them.
fn top_row_width(row: &[(Option<Series>, &str)], words: f64) -> f64 {
    let items: f64 = row.iter().map(|&i| item_width(i, words)).sum();
    items + KEY_ITEM_GAP * row.len().saturating_sub(1) as f64
}

/// A key's swatch and the room after it, before its name.
const KEY_SWATCH_ROOM: f64 = 10.0;
/// Room between a key at the top's items in a row.
const KEY_ITEM_GAP: f64 = 6.0;
/// A key's row height.
const KEY_ROW: f64 = 10.0;

/// The least of a panel's width a key at the right leaves its plot.
const PLOT_LEAST: f64 = 0.6;
/// The sizes a key at the right may be drawn at, of its own, the largest
/// that leaves the plot [`PLOT_LEAST`] first.
const KEY_SHRINK: [f64; 2] = [1.0, 0.85];

/// The most spaces between the ticks a histogram spreads them over.
const TICKS: usize = 6;

/// Of `choices` (position, text, wanted first), those whose text fits
/// clear of the others, each placed by `x_of`, at `size`, inside
/// `bounds`: the wanted first, then the most of the rest that fit clear
/// of each other, each the nearest to one of a few evenly spaced places
/// along the line (as many as fit, measured by their texts' widths, up to
/// [`TICKS`] spaces).
fn pick_ticks(
    choices: &[(f64, String, bool)],
    x_of: impl Fn(f64) -> f64,
    size: f64,
    (left, right): (f64, f64),
) -> Vec<(f64, String)> {
    // A tick's place across and half its label's room.
    let place = |at: f64, text: &str| (x_of(at), text_width(text, size) / 2.0 + 1.5);
    let fits = |(x, half): (f64, f64)| x - half >= left && x + half <= right;
    let clear = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).abs() >= a.1 + b.1;
    let mut wanted: Vec<(f64, f64, &str)> = Vec::new();
    for (at, text, _) in choices.iter().filter(|c| c.2) {
        let p = place(*at, text);
        if fits(p) && wanted.iter().all(|w| clear(p, (w.0, w.1))) {
            wanted.push((p.0, p.1, text));
        }
    }
    let plain: Vec<(f64, f64, &str)> = (choices.iter().filter(|c| !c.2))
        .map(|(at, text, _)| {
            let p = place(*at, text);
            (p.0, p.1, text.as_str())
        })
        .filter(|&(x, half, _)| {
            fits((x, half)) && wanted.iter().all(|w| clear((x, half), (w.0, w.1)))
        })
        .collect();
    let (lo, hi) = (plain.iter()).fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), c| {
        (l.min(c.0), h.max(c.0))
    });
    // For each count of spaces, the nearest to each evenly spaced place;
    // the most that are all clear of each other.
    let mut best: Vec<(f64, f64, &str)> = Vec::new();
    for spaces in 1..=TICKS {
        let mut set: Vec<(f64, f64, &str)> = Vec::new();
        for k in 0..=spaces {
            let target = lo + (hi - lo) * k as f64 / spaces as f64;
            let nearest =
                (plain.iter()).min_by(|a, b| (a.0 - target).abs().total_cmp(&(b.0 - target).abs()));
            if let Some(&n) = nearest {
                if !set.iter().any(|s| s.0 == n.0) {
                    set.push(n);
                }
            }
        }
        let all_clear = (0..set.len())
            .all(|i| (i + 1..set.len()).all(|j| clear((set[i].0, set[i].1), (set[j].0, set[j].1))));
        if all_clear && set.len() > best.len() {
            best = set;
        }
    }
    let mut chosen: Vec<(f64, String)> = (wanted.into_iter().chain(best))
        .map(|(x, _, t)| (x, t.to_string()))
        .collect();
    chosen.sort_by(|a, b| a.0.total_cmp(&b.0));
    chosen
}

impl ArmHist {
    /// Drawn in the box `(x, y, w, h)` of `canvas` as `how` says: the
    /// title over it unless alone, the y counts left, the bars, each
    /// threshold's line with its words over the bars, the ticks at the
    /// bars' edges, and the key at the right (alone, only with `--key`).
    pub(super) fn draw(
        &self,
        canvas: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        how: &HistLook,
    ) {
        let style = *canvas.style();
        let alone = style.panel_only;
        let ys = how.y_scale;
        // The axis ends at its top tick, a round count over the tallest bar.
        let ticks = ys.ticks(self.tallest() as f64);
        let tallest = ys.of(ticks[ticks.len() - 1]);
        // Round counts where they sit, as a fraction of the axis.
        let labels: Vec<(f64, String)> = (ticks.into_iter())
            .map(|v| match tallest > 0.0 {
                true => (ys.of(v) / tallest, crate::tui::ui::compact(v)),
                false => (0.0, "0".to_string()),
            })
            .collect();
        let size8 = canvas.label_size(8.0);
        let widest = |texts: &mut dyn Iterator<Item = String>, size: f64| {
            texts.map(|t| text_width(&t, size)).fold(0.0, f64::max)
        };
        let left = 14.0 + widest(&mut labels.iter().map(|l| l.1.clone()), size8) + 5.0;
        let keyed = (!alone || how.look.key) && how.look.key_at != KeyAt::Alone;
        let words = canvas.label_size(WORDS);
        // At the right only where the plot keeps most of the panel: a key
        // too wide drawn smaller, else above the plot.
        let (key_at, key_f) = match how.look.key_at {
            KeyAt::Right => (KEY_SHRINK.into_iter())
                .find(|&f| w - left - 6.0 - self.key_size(words * f).0 - 2.0 >= PLOT_LEAST * w)
                .map_or((KeyAt::Top, 1.0), |f| (KeyAt::Right, f)),
            at => (at, 1.0),
        };
        // A key at the top: its rows across the panel, up to the plot's
        // right end, taking their height off the plot's.
        let top_rows = match keyed && key_at == KeyAt::Top {
            true => self.top_key(words, w - 8.0),
            false => Vec::new(),
        };
        let key_w = match keyed && key_at == KeyAt::Right {
            true => self.key_size(words * key_f).0 + 2.0,
            false => 0.0,
        };
        // The x title on a line, else broken after its "; " on two.
        let x_lines: Vec<&str> = match how.x_title.split_once("; ") {
            Some((a, b)) if text_width(&how.x_title, size8) > w - 4.0 => {
                vec![&how.x_title[..a.len() + 1], b]
            }
            _ => vec![how.x_title.as_str()],
        };
        let line_h = size8 * 1.2;
        let bottom = 28.0 + (x_lines.len() - 1) as f64 * line_h;
        let head_top = if alone { 6.0 } else { 20.0 };
        let (top, right) = (head_top + top_rows.len() as f64 * KEY_ROW, 6.0 + key_w);
        let (px, py) = (x + left, y + top);
        let (pw, ph) = ((w - left - right).max(1.0), (h - top - bottom).max(1.0));
        if !alone {
            let title = figure::fit_text(&how.title, w - left, canvas.label_size(11.0));
            canvas.bold(px, y + 13.0, &title, 11.0, Anchor::Start, INK);
        }
        let (a0, a1) = (self.at[0], self.at[self.at.len() - 1]);
        let span = (a1 - a0).max(1e-12);
        let x_of = |a: f64| px + (a - a0) / span * pw;
        // The bars stand under the rows of words the thresholds' lines take.
        // Over the bars the y labels' column is free: the words may reach
        // into it, not into the key, beside the plot or over its corner.
        let right_end = match keyed && key_at == KeyAt::Inside {
            true => px + pw - self.key_size(words).0 - 4.0,
            false => px + pw + 6.0,
        };
        let bounds = (x + 2.0, right_end);
        let placed = self.place_words(&x_of, bounds, |t: &str| set_width(t, words));
        // Room under the last row for the top count's label beside it.
        let head = match placed.iter().map(|p| p.row + 1).max() {
            None => 0.0,
            Some(rows) => rows as f64 * WORDS_ROW + 5.0,
        };
        let bh = (ph - head).max(1.0);
        let base = py + ph;
        for (f, label) in &labels {
            let gy = base - f * bh;
            if *f > 0.0 {
                canvas.line(px, gy, px + pw, gy, FAINT, 0.4);
            }
            canvas.text(px - 5.0, gy + 3.0, label, 8.0, Anchor::End, MUTED);
        }
        canvas.vtext(x + 8.0, py + ph / 2.0, &how.y_title, 8.0, MUTED);
        let kept_tone = style.palette.roles().bar;
        let height = |v: usize| match tallest > 0.0 {
            true => ys.of(v as f64) / tallest * bh,
            false => 0.0,
        };
        for series in self.order(how.look.front) {
            for j in 0..self.bars() {
                let bar_h = height(self.value(series, j));
                let (x0, x1) = (x_of(self.at[j]), x_of(self.at[j + 1]));
                let gap = if x1 - x0 > 3.0 { 0.15 * (x1 - x0) } else { 0.0 };
                if bar_h > 0.0 {
                    let b = (x0 + gap / 2.0, base - bar_h, x1 - x0 - gap, bar_h);
                    series.draw(canvas, b, kept_tone);
                }
            }
        }
        canvas.line(px, base, px + pw, base, INK, 0.6);
        // The y axis up to the bars' full scale: over it the words' rows.
        canvas.line(px, base - bh, px, base, INK, 0.6);
        for &(_, _, e) in &self.cuts {
            let lx = x_of(self.at[e]);
            canvas.dashed(lx, py, lx, base, arms::THRESHOLD, 0.9);
        }
        for p in &placed {
            let ty = py + 7.0 + p.row as f64 * WORDS_ROW;
            canvas.text(p.from, ty, &p.text, WORDS, Anchor::Start, INK);
        }
        let choices = self.tick_choices(how.scale, how.look.p_axis);
        for (tx, text) in pick_ticks(&choices, x_of, size8, (x + 1.0, x + w - 1.0)) {
            canvas.line(tx, base, tx, base + 3.0, INK, 0.6);
            canvas.text(tx, base + 11.0, &text, 8.0, Anchor::Middle, INK);
        }
        let below = y + h - 5.0 - (x_lines.len() - 1) as f64 * line_h;
        for (k, line) in x_lines.iter().enumerate() {
            let (at, said) = figure::place_title(canvas, line, px + pw / 2.0, (x, w), 8.0);
            let ty = below + k as f64 * line_h;
            canvas.text(at, ty, &said, 8.0, Anchor::Middle, MUTED);
        }
        match (keyed, key_at) {
            (true, KeyAt::Inside) => {
                let (kw, kh) = self.key_size(words);
                let (lx, ly) = (px + pw - kw, py);
                canvas.rect(lx - 2.0, ly, kw + 2.0, kh, "#ffffff");
                self.draw_key(canvas, (lx, ly), kept_tone, 1.0);
            }
            (true, KeyAt::Top) => {
                // The rows from one left edge, the widest ending at the
                // plot's right.
                let widest = (top_rows.iter())
                    .map(|row| top_row_width(row, words))
                    .fold(0.0, f64::max);
                let from = (px + pw - widest).max(x + 2.0);
                for (k, row) in top_rows.iter().enumerate() {
                    let ry = y + head_top - 4.0 + k as f64 * KEY_ROW;
                    self.draw_key_row(canvas, row, (from, ry), kept_tone, words);
                }
            }
            (true, _) => self.draw_key(canvas, (px + pw + 8.0, py), kept_tone, key_f),
            (false, _) => {}
        }
    }

    /// A key at the top's `row` from `(lx, ly)` across: each item's swatch
    /// in `kept_tone` for the retained, then its name, measured at `words`.
    fn draw_key_row(
        &self,
        canvas: &mut Canvas,
        row: &[(Option<Series>, &str)],
        (mut lx, ly): (f64, f64),
        kept_tone: &str,
        words: f64,
    ) {
        for &(series, said) in row {
            let y = ly + 2.0;
            if let Some(series) = series {
                series.draw(canvas, (lx, y, 7.0, 7.0), kept_tone);
                lx += KEY_SWATCH_ROOM;
            }
            canvas.text(lx, y + 6.3, said, WORDS, Anchor::Start, MUTED);
            lx += set_width(said, words) + KEY_ITEM_GAP;
        }
    }

    /// Its key's width and height at `words` points of type: a swatch and
    /// its series' name a row, then the note; all in proportion to the
    /// type.
    fn key_size(&self, words: f64) -> (f64, f64) {
        let f = words / WORDS;
        let (items, note) = self.legend();
        let names = items.iter().map(|i| text_width(i.1, words) + 10.0 * f);
        let w = names
            .chain(note.map(|n| text_width(n, words)))
            .fold(0.0, f64::max);
        let rows = items.len() + usize::from(note.is_some());
        (w + 2.0, 4.0 + rows as f64 * 10.0 * f)
    }

    /// Its key, from `(lx, ly)` down, `f` of its own size: each series'
    /// swatch in `kept_tone` for the retained, its name, then the note.
    fn draw_key(&self, canvas: &mut Canvas, (lx, ly): (f64, f64), kept_tone: &str, f: f64) {
        let (items, note) = self.legend();
        let (row_h, swatch, size) = (10.0 * f, 7.0 * f, WORDS * f);
        for (row, (series, said)) in items.iter().enumerate() {
            let y = ly + 2.0 + row as f64 * row_h;
            series.draw(canvas, (lx, y, swatch, swatch), kept_tone);
            canvas.text(lx + 10.0 * f, y + 6.3 * f, said, size, Anchor::Start, MUTED);
        }
        if let Some(note) = note {
            let y = ly + 2.0 + items.len() as f64 * row_h;
            canvas.text(lx, y + 6.3 * f, note, size, Anchor::Start, MUTED);
        }
    }

    /// Its key alone, as a figure in `style` just its size: `--key-at
    /// alone`'s, one key for a row of panels.
    pub(super) fn key_figure(&self, style: &figure::style::Style) -> String {
        let probe = Canvas::styled(1.0, 1.0, style);
        let (kw, kh) = self.key_size(probe.label_size(WORDS));
        let mut canvas = Canvas::styled(kw + 4.0, kh + 2.0, style);
        let kept_tone = style.palette.roles().bar;
        self.draw_key(&mut canvas, (2.0, 0.0), kept_tone, 1.0);
        canvas.finish()
    }

    /// Drawn into the terminal's own cells in `area`, without images, as
    /// [`ArmHist::draw`] draws it: the control arm dim, the signal arm's
    /// filtered shaded and its kept solid (the arm in front over the
    /// other), the threshold's line and words accented, the ticks below.
    pub(super) fn render_text(&self, buf: &mut Buffer, area: Rect, how: &HistLook) {
        let put = |buf: &mut Buffer, x: u16, y: u16, s: &str, style: Style| {
            buf[(x, y)]
                .set_symbol(s)
                .set_style(Style::reset().patch(style));
        };
        let [plot, axis, labels] = Layout::vertical([
            Constraint::Min(1),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(area);
        let [gutter, chart] =
            Layout::horizontal([Constraint::Length(GUTTER), Constraint::Min(1)]).areas(plot);
        let words_rows = u16::from(!self.cuts.is_empty());
        if chart.width < 2 || chart.height <= words_rows {
            return;
        }
        let bars_top = chart.top() + words_rows;
        let rows = (chart.bottom() - bars_top) as usize;
        let width = chart.width as usize;
        let (a0, a1) = (self.at[0], self.at[self.at.len() - 1]);
        let span = (a1 - a0).max(1e-12);
        let col_of = |a: f64| (((a - a0) / span * width as f64) as usize).min(width - 1) as u16;
        // The bar under each column's middle.
        let bar_of = |col: usize| {
            let a = a0 + (col as f64 + 0.5) / width as f64 * span;
            (self.at.partition_point(|&e| e <= a).saturating_sub(1)).min(self.bars() - 1)
        };
        let ys = how.y_scale;
        let tallest = ys.of(self.tallest() as f64);
        let cells = rows * 8;
        let units = |v: usize| match v > 0 && tallest > 0.0 {
            true => ((ys.of(v as f64) / tallest * cells as f64).round() as usize).clamp(1, cells),
            false => 0,
        };
        const EIGHTHS: [&str; 8] = ["▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"];
        let cut_cols: Vec<u16> = (self.cuts.iter())
            .map(|cut| col_of(self.at[cut.2]))
            .collect();
        for &cx in &cut_cols {
            for y in bars_top..chart.bottom() {
                put(buf, chart.x + cx, y, "┊", ACCENTED);
            }
        }
        for series in self.order(how.look.front) {
            for col in 0..width {
                let u = units(self.value(series, bar_of(col)));
                for (j, y) in (bars_top..chart.bottom()).rev().enumerate() {
                    let fill = u.saturating_sub(j * 8).min(8);
                    if fill == 0 {
                        break;
                    }
                    let (glyph, style) = match series {
                        Series::Mutant => (EIGHTHS[fill - 1], DIM),
                        Series::Filtered => ("░", PLAIN),
                        _ => (EIGHTHS[fill - 1], PLAIN),
                    };
                    put(buf, chart.x + col as u16, y, glyph, style);
                }
            }
        }
        // The words either side of each line, on the row over the bars,
        // placed as the figure places them, a column a character.
        let col_x = |a: f64| f64::from(chart.x + col_of(a));
        let bounds = (f64::from(chart.x), f64::from(chart.right()));
        let count = |t: &str| t.chars().count() as f64;
        for p in self
            .place_words(&col_x, bounds, count)
            .iter()
            .filter(|p| p.row == 0)
        {
            let from = (p.from.max(f64::from(chart.x))) as u16;
            let room = chart.right().saturating_sub(from) as usize;
            buf.set_stringn(from, chart.y, &p.text, room, ACCENTED);
        }
        // The y axis: the count at the top and halfway, on the y scale.
        let gx = gutter.right() - 1;
        for y in gutter.top()..gutter.bottom() {
            put(buf, gx, y, "│", DIM);
        }
        let mut y_label = |y: u16, f: f64| {
            let s = crate::tui::ui::compact(ys.count(f * tallest));
            let x = gx.saturating_sub(1 + s.len() as u16).max(gutter.x);
            buf.set_string(x, y, &s, DIM);
            put(buf, gx, y, "┤", DIM);
        };
        if tallest > 0.0 {
            y_label(bars_top, 1.0);
            if rows >= 6 {
                y_label(bars_top + rows as u16 / 2, 0.5);
            }
        }
        for x in axis.left()..axis.right() {
            let sym = match x.cmp(&gx) {
                std::cmp::Ordering::Less => " ",
                std::cmp::Ordering::Equal => "└",
                std::cmp::Ordering::Greater => "─",
            };
            put(buf, x, axis.y, sym, DIM);
        }
        // The ticks, each label starting at its edge, clear of the last.
        let choices = self.tick_choices(how.scale, how.look.p_axis);
        let x_of = |a: f64| f64::from(chart.x + col_of(a));
        let bounds = (f64::from(labels.x), f64::from(labels.right()));
        for (x, text) in pick_ticks(&choices, x_of, 1.0 / 0.55, bounds) {
            let x = x as u16;
            put(buf, x, axis.y, "┴", DIM);
            let half = text.chars().count() as u16 / 2;
            buf.set_string(x.saturating_sub(half).max(labels.x), labels.y, &text, DIM);
        }
    }
}

impl SitePicker<'_> {
    /// The focused knob's histogram on both arms, at the thresholds now.
    pub(super) fn binned_arms(&self) -> ArmHist {
        let c = self.criterion();
        ArmHist::new(&self.view().sites(), c, &self.filter, self.column.hist.bins)
    }
}

/// The knob's x axis' own title: what it shows, each arm's own reads where
/// it has a value on each (the coverage's threshold adds the two), a
/// p-value as `p_axis` reads it (each arm's own where it has one).
pub(super) fn arm_axis(
    c: Criterion,
    own: &'static str,
    per_arm: bool,
    p_axis: PAxis,
) -> &'static str {
    match c {
        // The signal arm on the threshold's own test, the control arm on
        // its own reads tested against sequencing error.
        Criterion::MaxPv if per_arm => match p_axis {
            PAxis::P => "p: WT vs MUT (Fisher); MUT vs error",
            PAxis::Neglog10 => "-log10 p: WT vs MUT (Fisher); MUT vs error",
        },
        Criterion::MaxPv => p_axis.name(),
        Criterion::MinCoverage if per_arm => "reads per arm",
        Criterion::MinConverted if per_arm => "converted reads per arm",
        _ => own,
    }
}

/// The knob's titles: its x axis' as asked (`--xlab`, `t` in the view)
/// or its own, and its counts' as asked or "sites".
pub(super) fn arm_titles(c: Criterion, own: String, titles: &Titles) -> (String, String) {
    (titles.x_or(Axis::of(c), &own), titles.sites())
}
