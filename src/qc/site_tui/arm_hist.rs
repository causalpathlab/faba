//! A knob's histogram as the view draws it, on screen, in its figure and
//! as a `--style view` panel: the control arm's sites (MUT) light gray
//! behind, the signal arm's (WT) in front, white and outlined where some
//! threshold filters the site out and dark where every one keeps it (or
//! the control arm in front, `--front mut`). The bars are the knob's own
//! bins, each threshold's bin split at it, so its dashed line sits on a
//! cut between bars with the words for either side beside it. Every count
//! is read off the fail bits: nothing here decides a site.

use super::*;
use crate::figure::{text_width, FAINT};
use crate::qc::arms::{self, ArmLook, Front, PAxis};
use crate::qc::site_plot::Titles;
use ratatui::buffer::Buffer;

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
        let mutant: Option<Vec<f64>> = (arms::per_arm(q, s.table)).then(|| {
            (0..n)
                .into_par_iter()
                .map(|i| arms::mutant(s, q, i))
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
        let clamped = (
            values().any(|v| v == f64::INFINITY).then_some(hi),
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
        let mut out: Vec<(f64, String, bool)> = (self.edges.iter().zip(&self.at))
            .zip(&self.split)
            .filter(|(_, &split)| !split)
            .map(|((&e, &at), _)| (at, self.tick_text(e, p_axis), false))
            .collect();
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

    /// The series its key lists, the one in front first, then a note for
    /// a quantity with no control arm drawn.
    fn legend(&self, front: Front) -> (Vec<(Series, &'static str)>, Option<&'static str>) {
        let wt = [(Series::Kept, "WT kept"), (Series::Filtered, "WT filtered")];
        let mut items: Vec<(Series, &str)> = Vec::new();
        match (self.mutant.is_some(), front) {
            (true, Front::Mut) => {
                items.push((Series::Mutant, "MUT"));
                items.extend(wt);
            }
            (true, Front::Wt) => {
                items.extend(wt);
                items.push((Series::Mutant, "MUT"));
            }
            (false, _) => items.extend(wt),
        }
        if !self.cuts.is_empty() {
            items.push((Series::Line, "threshold"));
        }
        let note = match (self.mutant.is_some(), self.has_control) {
            (true, _) => None,
            (false, true) => Some("WT only"),
            (false, false) => Some("no MUT arm"),
        };
        (items, note)
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
            Series::Line => 0,
        }
    }

    /// The tallest bar of either arm.
    fn tallest(&self) -> usize {
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
    Line,
}

impl Series {
    /// Drawn in the box `b`: the palette's kept bar colour `kept`.
    fn draw(self, c: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64), kept: &str) {
        match self {
            Series::Mutant => c.rect(x, y, w, h, arms::MUT_FILL),
            Series::Kept => c.rect(x, y, w, h, kept),
            Series::Filtered => c.outlined((x, y, w, h), arms::FILTERED_FILL, (arms::OUTLINE, 0.4)),
            Series::Line => c.dashed(x, y + h / 2.0, x + w, y + h / 2.0, arms::THRESHOLD, 0.9),
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
    /// words has them on the other, an arrow pointing back across the
    /// line. `x_of` places a position along the line, `width` measures a
    /// text.
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
            match lx - WORDS_GAP - wl >= left {
                true => put(low, lx - WORDS_GAP - wl),
                false => put(format!("← {low}"), lx + WORDS_GAP),
            }
            match lx + WORDS_GAP + wh <= right {
                true => put(high, lx + WORDS_GAP),
                false => {
                    let text = format!("{high} →");
                    let from = (lx - WORDS_GAP - width(&text)).max(left);
                    put(text, from);
                }
            }
        }
        placed
    }
}

/// How a knob's histogram is drawn beyond its bars: its y scale, its arms,
/// its titles.
pub(super) struct HistLook {
    pub(super) y_scale: Scale,
    /// The scale its bins are on.
    pub(super) scale: Scale,
    pub(super) look: ArmLook,
    pub(super) title: String,
    pub(super) x_title: String,
    pub(super) y_title: String,
}

/// Of `choices` (position, text, wanted first), those whose text fits
/// clear of the others, about six: the wanted first, then evenly spread;
/// each placed by `x_of`, at `size`, inside `bounds`.
fn pick_ticks(
    choices: &[(f64, String, bool)],
    x_of: impl Fn(f64) -> f64,
    size: f64,
    (left, right): (f64, f64),
) -> Vec<(f64, String)> {
    let mut chosen: Vec<(f64, f64, String)> = Vec::new();
    let mut take = |at: f64, text: &str| {
        let (x, half) = (x_of(at), text_width(text, size) / 2.0 + 1.5);
        let clear = (chosen.iter()).all(|&(cx, ch, _)| (x - cx).abs() >= half + ch);
        if clear && x - half >= left && x + half <= right {
            chosen.push((x, half, text.to_string()));
        }
    };
    for (at, text, _) in choices.iter().filter(|c| c.2) {
        take(*at, text);
    }
    let plain: Vec<&(f64, String, bool)> = choices.iter().filter(|c| !c.2).collect();
    let every = plain.len().div_ceil(6).max(1);
    for (at, text, _) in plain.iter().step_by(every) {
        take(*at, text);
    }
    chosen.sort_by(|a, b| a.0.total_cmp(&b.0));
    chosen.into_iter().map(|(x, _, t)| (x, t)).collect()
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
        let tallest = figure::scaled(ys, self.tallest() as f64);
        let y_label = |f: f64| match tallest > 0.0 {
            true => crate::tui::ui::compact(figure::unscaled(ys, f * tallest)),
            false => "0".to_string(),
        };
        let labels = [0.0, 0.5, 1.0].map(|f| (f, y_label(f)));
        let size8 = canvas.label_size(8.0);
        let widest = |texts: &mut dyn Iterator<Item = String>, size: f64| {
            texts.map(|t| text_width(&t, size)).fold(0.0, f64::max)
        };
        let left = 14.0 + widest(&mut labels.iter().map(|l| l.1.clone()), size8) + 5.0;
        let (items, note) = self.legend(how.look.front);
        let keyed = !alone || how.look.key;
        let words = canvas.label_size(WORDS);
        let key_w = match keyed {
            true => {
                let texts = items
                    .iter()
                    .map(|i| i.1.to_string())
                    .chain(note.map(String::from));
                12.0 + widest(&mut texts.into_iter(), words) + 2.0
            }
            false => 0.0,
        };
        let (top, right, bottom) = (if alone { 6.0 } else { 20.0 }, 6.0 + key_w, 28.0);
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
        // into it, not into the key.
        let bounds = (x + 2.0, px + pw + 6.0);
        let placed = self.place_words(&x_of, bounds, |t: &str| text_width(t, words));
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
            true => figure::scaled(ys, v as f64) / tallest * bh,
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
        let (at, x_title) = figure::place_title(canvas, &how.x_title, px + pw / 2.0, (x, w), 8.0);
        canvas.text(at, y + h - 5.0, &x_title, 8.0, Anchor::Middle, MUTED);
        if keyed {
            let lx = px + pw + 8.0;
            for (row, (series, said)) in items.iter().enumerate() {
                let ly = py + 2.0 + row as f64 * 10.0;
                series.draw(canvas, (lx, ly, 7.0, 7.0), kept_tone);
                canvas.text(lx + 10.0, ly + 6.3, said, WORDS, Anchor::Start, MUTED);
            }
            if let Some(note) = note {
                let ly = py + 2.0 + items.len() as f64 * 10.0;
                canvas.text(lx, ly + 6.3, note, WORDS, Anchor::Start, MUTED);
            }
        }
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
        let tallest = figure::scaled(ys, self.tallest() as f64);
        let cells = rows * 8;
        let units = |v: usize| match v > 0 && tallest > 0.0 {
            true => ((figure::scaled(ys, v as f64) / tallest * cells as f64).round() as usize)
                .clamp(1, cells),
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
            let s = crate::tui::ui::compact(figure::unscaled(ys, f * tallest));
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
/// p-value as `p_axis` reads it; then the bins' scale.
pub(super) fn arm_axis(
    c: Criterion,
    own: &'static str,
    per_arm: bool,
    p_axis: PAxis,
) -> &'static str {
    match c {
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
