//! The two grids, gene modules and sites, each row a gene module or a site
//! against every cell module: mirror bars with their methylated %, a %
//! heatmap or a z heatmap, the heatmaps' boxes sized by their reads.

use super::bars::{mirror, peak};
use super::{ModulesView, Spot};
use crate::cluster::browse::{standardized, Count, Modules};
use crate::site_analysis::pileup::{channel_names, query_symbol};
use crate::tui::hidden_marks;
use data_beans::interactive::ui::{compact, DIM, HIGHLIGHT, PLAIN};
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Style};
use std::ops::Range;

/// Columns per cell module in the grids: a bar and a 3-digit %, or a box.
pub(super) const CELL: u16 = 6;

/// Columns per cell module when not all fit at [`CELL`]: a bar or a box,
/// and a gap; no numbers.
pub(super) const COMPACT: u16 = 3;

/// Width of the module labels, `M12  1234 genes`.
pub(super) const MODULE_LABEL: u16 = 17;

/// Room after a module row's cells for its shift range, ` Δ1.23`.
pub(super) const MODULE_EXTRA: u16 = 8;

/// Widest site label, `• GENE1 123456789`.
pub(super) const SITE_LABEL: u16 = 24;

/// Width of a site's counts, `reads m u cells`.
pub(super) const COUNTS: u16 = 27;

/// How a grid draws its blocks; `v` cycles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Look {
    Bars,
    Percent,
    Z,
}

impl Look {
    pub(super) fn next(self) -> Self {
        match self {
            Look::Bars => Look::Percent,
            Look::Percent => Look::Z,
            Look::Z => Look::Bars,
        }
    }
}

/// The % heatmap's scale; `y` cycles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ramp {
    Sqrt,
    Linear,
    Log10,
}

impl Ramp {
    pub(super) fn next(self) -> Self {
        match self {
            Ramp::Sqrt => Ramp::Linear,
            Ramp::Linear => Ramp::Log10,
            Ramp::Log10 => Ramp::Sqrt,
        }
    }

    pub(super) fn name(self) -> &'static str {
        match self {
            Ramp::Sqrt => "sqrt",
            Ramp::Linear => "linear",
            Ramp::Log10 => "log10",
        }
    }

    /// A fraction on this scale; log10 runs from 0.1% to 100%.
    pub fn of(self, f: f64) -> f64 {
        match self {
            Ramp::Sqrt => f.max(0.0).sqrt(),
            Ramp::Linear => f,
            Ramp::Log10 => (f.max(1e-3).log10() + 3.0) / 3.0,
        }
    }
}

/// A grid's colour ranges over all its rows: the % heatmap's (transformed)
/// 2nd and 98th percentiles, the z heatmap's 98th percentile of |z|, and the
/// deepest block, for the boxes' sizes.
pub struct Stretch {
    pub lo: f64,
    pub hi: f64,
    pub z_top: f64,
    pub deepest: f64,
}

impl Stretch {
    pub fn of<'r>(rows: impl Iterator<Item = &'r [Count]>, ramp: Ramp) -> Self {
        let (mut v, mut z, mut deepest) = (Vec::new(), Vec::new(), 0.0f64);
        for reads in rows {
            v.extend(
                reads
                    .iter()
                    .filter_map(|r| r.fraction())
                    .map(|f| ramp.of(f)),
            );
            z.extend(standardized(reads).into_iter().flatten().map(f64::abs));
            deepest = reads.iter().fold(deepest, |m, r| m.max(r.total));
        }
        let at = |x: &mut Vec<f64>, q: f64| {
            x.sort_by(f64::total_cmp);
            let i = (q * x.len().saturating_sub(1) as f64).round() as usize;
            x.get(i).copied()
        };
        Self {
            lo: at(&mut v, 0.02).unwrap_or(0.0),
            hi: at(&mut v, 0.98).unwrap_or(1.0),
            z_top: at(&mut z, 0.98).filter(|&t| t > 0.0).unwrap_or(1.0),
            deepest,
        }
    }

    /// A % box's gray: black at `lo`, white at `hi`, clipped beyond.
    pub fn gray(&self, ramp: Ramp, f: f64) -> u8 {
        let span = (self.hi - self.lo).max(1e-12);
        (255.0 * ((ramp.of(f) - self.lo) / span).clamp(0.0, 1.0)).round() as u8
    }

    /// A z box's colour: white at 0, full red (high) or blue (low) at
    /// `z_top` and beyond; and how far toward it, 0 to 1.
    pub fn diverging(&self, z: f64) -> (Color, f64) {
        let t = (z.abs() / self.z_top).min(1.0);
        let fade = (255.0 * (1.0 - t)).round() as u8;
        let c = if z >= 0.0 {
            Color::Rgb(255, fade, fade)
        } else {
            Color::Rgb(fade, fade, 255)
        };
        (c, t)
    }
}

/// What a frame needs worked out over a whole module, kept with the key it
/// was worked out for, so a frame works it out only when the key changes:
/// by default a grid's colour ranges (for the scale, and for the sites grid
/// the gene role too).
pub(in crate::cluster::browse) struct Kept<K, V = Stretch>(Option<(K, V)>);

impl<K, V> Default for Kept<K, V> {
    fn default() -> Self {
        Self(None)
    }
}

impl<K: Copy + PartialEq, V> Kept<K, V> {
    /// The key it was worked out for; none before it was.
    pub(in crate::cluster::browse) fn key(&self) -> Option<K> {
        self.0.as_ref().map(|&(k, _)| k)
    }

    /// What `make` works out for `key`, anew only when it was worked out
    /// for another key.
    fn of(&mut self, key: K, make: impl FnOnce() -> V) -> &V {
        if self.key() != Some(key) {
            self.0 = Some((key, make()));
        }
        &self.0.as_ref().expect("just kept").1
    }
}

impl<K: Copy + PartialEq> Kept<K> {
    /// The ranges worked out by `make` for `key` when `look` is a heatmap
    /// and they were worked out for another key; bars need none.
    fn update(&mut self, look: Look, key: K, make: impl FnOnce() -> Stretch) {
        if look != Look::Bars && self.key() != Some(key) {
            self.0 = Some((key, make()));
        }
    }

    /// The ranges `look` draws with: none for bars.
    fn get(&self, look: Look) -> Option<&Stretch> {
        self.0
            .as_ref()
            .filter(|_| look != Look::Bars)
            .map(|(_, s)| s)
    }
}

/// The rows a grid shows, from row `first` of all of them, and how: the
/// label column's width and head, the row highlighted, the look, and the
/// colour ranges over all the rows (for a heatmap).
struct Placing<'s> {
    label_w: u16,
    /// What the labels count, over their last column.
    head: &'static str,
    /// Room kept right of the cells for each row's numbers.
    extra_w: u16,
    first: usize,
    selected: usize,
    look: Look,
    stretch: Option<&'s Stretch>,
}

/// A grid's corner `room` wide: the cell modules `hidden`, then `head`
/// ending a column short of the room, when both fit apart.
fn corner(hidden: &str, head: &str, room: usize) -> String {
    let lead = hidden.chars().count();
    let gap = room.saturating_sub(1 + lead + head.chars().count());
    if head.is_empty() || (gap == 0 && lead > 0) {
        hidden.to_string()
    } else {
        format!("{hidden}{}{head}", " ".repeat(gap))
    }
}

/// One row of a grid, its reads in the panes' column order.
struct Row<'r> {
    label: String,
    reads: &'r [Count],
    extra: String,
    spot: Spot,
}

impl ModulesView<'_> {
    /// The modules grid: each gene module's reads in each cell module.
    pub(super) fn draw_modules(&mut self, buf: &mut Buffer, area: Rect) {
        let m = self.m;
        let look = self.looks[0];
        // Colours over every module's blocks, worked out once per scale.
        let ramp = self.ramp;
        self.module_stretch.update(look, ramp, || {
            Stretch::of(m.order.iter().map(|&g| m.block(g)), ramp)
        });
        let order = self.module_order();
        let fit = (area.height.saturating_sub(1) / 2).max(1) as usize;
        let at = order.iter().position(|&g| g == self.module).unwrap_or(0);
        let first = (at + 1).saturating_sub(fit);
        let rows: Vec<Row> = order[first..order.len().min(first + fit)]
            .iter()
            .map(|&g| Row {
                label: format!("M{:<3}{:>5} genes", g + 1, m.gene_sizes[g]),
                reads: m.block(g),
                extra: format!("Δ{:.2}", m.shift_range(g)),
                spot: Spot::Module(g),
            })
            .collect();
        let placing = Placing {
            label_w: MODULE_LABEL,
            head: "main genes",
            extra_w: MODULE_EXTRA,
            first,
            selected: at,
            look,
            stretch: self.module_stretch.get(look),
        };
        self.draw_grid(buf, area, &rows, placing);
    }

    /// The sites grid: the module's sites in the list's order, one row
    /// each, scrolled to keep the one in focus in view, each with its
    /// counts over every cell read; a pinned gene's sites marked.
    pub(super) fn draw_sites(&mut self, buf: &mut Buffer, area: Rect) {
        let Some(g) = self.role() else { return };
        self.rows = (area.height.saturating_sub(1) / 2).max(1) as usize;
        self.select_site(self.site);
        let label_w = self.site_label_w(g);
        let members = &self.m.members[g];
        let list = &self.list;
        let first = self.from.min(list.len().saturating_sub(self.rows));
        let rows: Vec<Row> = list[first..(first + self.rows).min(list.len())]
            .iter()
            .map(|&i| {
                let s = &self.m.sites[members[i].site];
                let mark = if self.pinned(g, i) { PIN_MARK } else { "" };
                Row {
                    label: site_label(mark, &s.gene, s.pos),
                    reads: &members[i].reads,
                    extra: self.site_counts(members[i].site),
                    spot: Spot::Site(i),
                }
            })
            .collect();
        let look = self.looks[1];
        // Colours over every site of the module, not only those shown,
        // worked out once per module and scale.
        let ramp = self.ramp;
        self.site_stretch.update(look, (g, ramp), || {
            Stretch::of(members.iter().map(|s| s.reads.as_slice()), ramp)
        });
        let selected = self.listed_at().unwrap_or(0);
        let placing = Placing {
            label_w,
            head: "",
            extra_w: COUNTS + 1,
            first,
            selected,
            look,
            stretch: self.site_stretch.get(look),
        };
        let shown = self.draw_grid(buf, area, &rows, placing);
        let x = area.x + label_w + shown as u16 * self.cell_w;
        // The two channels by their first letters: m/u, e/u.
        let (on, off) = channel_names(&self.m.modality);
        let initial = |s: &str| s.chars().next().unwrap_or('?').to_string();
        let head = format!(
            "{:>6} {:>6} {:>6} {:>6}",
            "reads",
            initial(on),
            initial(off),
            "cells"
        );
        if x < area.right() {
            buf.set_stringn(x, area.y, head, (area.right() - x) as usize, DIM);
        }
        scrollbar(buf, area, first, self.rows, list.len());
    }

    /// The sites grid's label column for gene role `g`: its widest site
    /// label and a gap, over all its sites (so it holds as the list
    /// scrolls), at most [`SITE_LABEL`]; each gene's widest worked out
    /// once per role, its pin mark added as pinned.
    pub(in crate::cluster::browse) fn site_label_w(&mut self, g: usize) -> u16 {
        let m = self.m;
        let widths = self.site_labels.of(g, || label_widths(m, g));
        let pinned = self.pins.iter().filter(|p| p.role == g);
        let mark = PIN_MARK.chars().count();
        let marked = pinned.filter_map(|p| widths.get(p.gene)).map(|w| w + mark);
        let widest = widths.iter().copied().chain(marked).max().unwrap_or(0);
        (widest as u16).min(SITE_LABEL)
    }

    /// A site's raw reads over every cell read: all, converted, unconverted,
    /// and the cells with any.
    fn site_counts(&self, site: usize) -> String {
        match &self.m.stats {
            Some(s) => format!(
                "{:>6} {:>6} {:>6} {:>6}",
                compact((s.on[site] + s.off[site]) as f64),
                compact(s.on[site] as f64),
                compact(s.off[site] as f64),
                compact(s.cells[site] as f64)
            ),
            None => format!("{:>6} {:>6} {:>6} {:>6}", "-", "-", "-", "-"),
        }
    }

    /// The `rows` shown, placed as `at` says, two lines a row, under the
    /// names of the cell modules that fit from the first shown; those
    /// hidden either side counted in the corner. How many are shown.
    fn draw_grid(&self, buf: &mut Buffer, area: Rect, rows: &[Row], at: Placing) -> usize {
        let cell_w = self.cell_w;
        let cols = self.visible(cells_room(area.width, at.label_w, at.extra_w));
        for (k, &c) in self.m.columns[cols.clone()].iter().enumerate() {
            let x = area.x + at.label_w + k as u16 * cell_w;
            let name = format!("C{}", c + 1);
            // A compact head too long to keep a gap gives just the number.
            let name = if name.len() < cell_w as usize {
                name
            } else {
                (c + 1).to_string()
            };
            if x + 2 < area.right() {
                buf.set_stringn(x, area.y, name, cell_w as usize, DIM);
            }
        }
        let hidden = hidden_marks(cols.start, self.m.cell_roles - cols.end);
        let room = at.label_w.saturating_sub(1) as usize;
        buf.set_stringn(area.x, area.y, corner(&hidden, at.head, room), room, DIM);
        let fit = (area.height.saturating_sub(1) / 2) as usize;
        for (k, row) in rows.iter().take(fit).enumerate() {
            let selected = at.first + k == at.selected;
            let y = area.y + 1 + 2 * k as u16;
            self.hits.add(Rect::new(area.x, y, area.width, 2), row.spot);
            let style = if selected { HIGHLIGHT } else { PLAIN };
            let room = at.label_w.saturating_sub(1) as usize;
            buf.set_stringn(area.x, y, &row.label, room, style);
            let width = area.width.saturating_sub(at.label_w + 1);
            let grid = Rect::new(area.x + at.label_w, y, width, 2);
            let cells = (cols.clone(), cell_w);
            let end = match at.stretch {
                Some(s) => heat_row(
                    buf,
                    grid,
                    row.reads,
                    cells,
                    (at.look, self.ramp),
                    s,
                    selected,
                ),
                None => grid_row(buf, grid, row.reads, cells),
            };
            if end + 1 < area.right() {
                let room = (area.right() - end - 1) as usize;
                buf.set_stringn(end + 1, y, &row.extra, room, DIM);
            }
        }
        cols.len()
    }
}

/// The mark before a pinned gene's site labels.
const PIN_MARK: &str = "• ";

/// A site's label in the sites grid: `mark`, the gene's symbol and the
/// site's position. [`label_widths`] measures it.
fn site_label(mark: &str, gene: &str, pos: i64) -> String {
    format!("{mark}{} {pos}", query_symbol(gene))
}

/// Each gene of role `g`'s widest site label (see [`site_label`]),
/// unpinned, and a gap.
fn label_widths(m: &Modules, g: usize) -> Vec<usize> {
    let members = &m.members[g];
    let digits = |i: usize| crate::tui::digits(m.sites[members[i].site].pos);
    m.genes(g)
        .iter()
        .map(|span| {
            let symbol = query_symbol(m.gene_of(g, span.start)).chars().count();
            symbol + 1 + span.clone().map(digits).max().unwrap_or(0) + 1
        })
        .collect()
}

/// The columns a grid `width` wide inside its frame has for its cells, a
/// `label_w` label column and `extra_w` after them kept.
pub(super) fn cells_room(width: u16, label_w: u16, extra_w: u16) -> u16 {
    width.saturating_sub(label_w + extra_w)
}

/// One grid row of mirror bars for the cell modules `cols` of `reads`, in
/// cells `cell_w` wide, each with its methylated % to its right when there
/// is room, on the row's own scale over all of them; the column after the
/// last cell drawn.
fn grid_row(
    buf: &mut Buffer,
    area: Rect,
    reads: &[Count],
    (cols, cell_w): (Range<usize>, u16),
) -> u16 {
    let peak = peak(reads);
    let mut end = area.x;
    for (k, r) in reads[cols].iter().enumerate() {
        let x = area.x + k as u16 * cell_w;
        if x + cell_w > area.right() + 1 {
            break;
        }
        for dx in 0..2 {
            mirror(buf, x + dx, area.y, 1, *r, peak);
        }
        if let (Some(f), true) = (r.fraction(), cell_w == CELL) {
            let pct = format!("{:>3.0}", 100.0 * f);
            buf.set_string(x + 2, area.y, pct, DIM);
        }
        end = x + cell_w - 1;
    }
    end
}

/// One grid row of boxes two lines high for the cell modules `cols` of
/// `reads`, in cells `cell_w` wide, each sized by its reads against the
/// grid's deepest block and coloured by its % (on `ramp`) or its z (across
/// all of them); its number in it when `numbered` and there is room. The
/// column after the last cell drawn.
fn heat_row(
    buf: &mut Buffer,
    area: Rect,
    reads: &[Count],
    (cols, cell_w): (Range<usize>, u16),
    (look, ramp): (Look, Ramp),
    stretch: &Stretch,
    numbered: bool,
) -> u16 {
    let z = standardized(reads);
    let numbered = numbered && cell_w == CELL;
    let mut end = area.x;
    for (k, c) in cols.enumerate() {
        let r = &reads[c];
        let x = area.x + k as u16 * cell_w;
        if x + cell_w > area.right() + 1 {
            break;
        }
        end = x + cell_w - 1;
        let Some(f) = r.fraction() else { continue };
        let (fill, ink, text) = match look {
            Look::Z => {
                let z = z[c].unwrap_or(0.0);
                let (fill, t) = stretch.diverging(z);
                let ink = if t < 0.6 { Color::Black } else { Color::White };
                (fill, ink, format!("{z:+.1}"))
            }
            _ => {
                let g = stretch.gray(ramp, f);
                let ink = if g > 140 { Color::Black } else { Color::White };
                (Color::Rgb(g, g, g), ink, format!("{:.0}", 100.0 * f))
            }
        };
        let (w, h) = box_size(r.total, stretch.deepest, cell_w);
        let bx = x + (cell_w - 1 - w) / 2;
        for dy in 0..h {
            for dx in 0..w {
                buf.set_string(bx + dx, area.y + dy, " ", Style::new().bg(fill));
            }
        }
        if numbered {
            // Centred in the cell; the part over the box in its ink.
            let len = text.chars().count() as u16;
            let tx = x + (CELL - 1).saturating_sub(len) / 2;
            for (k, ch) in text.chars().enumerate() {
                let cx = tx + k as u16;
                let style = if (bx..bx + w).contains(&cx) {
                    Style::new().bg(fill).fg(ink)
                } else {
                    PLAIN
                };
                buf.set_string(cx, area.y, ch.to_string(), style);
            }
        }
    }
    end
}

/// A box's columns and lines: its area in proportion to `reads` against
/// `deepest`, at most `cell_w - 1` columns by two lines, at least one cell.
pub fn box_size(reads: f64, deepest: f64, cell_w: u16) -> (u16, u16) {
    let a = if deepest > 0.0 {
        (reads / deepest).clamp(0.0, 1.0).sqrt()
    } else {
        0.0
    };
    let w = ((a * (cell_w - 1) as f64).round() as u16).clamp(1, cell_w - 1);
    let h = if a >= 0.5 { 2 } else { 1 };
    (w, h)
}

/// A scrollbar down the right edge of `area` (below its head row) when
/// `len` rows do not fit in `shown`.
fn scrollbar(buf: &mut Buffer, area: Rect, from: usize, shown: usize, len: usize) {
    if len <= shown || area.height < 3 {
        return;
    }
    let x = area.right() - 1;
    let h = area.height - 1;
    let thumb = (h as usize * shown).div_ceil(len).max(1) as u16;
    let top = ((h - thumb) as usize * from / (len - shown)) as u16;
    for k in 0..h {
        let on = (top..top + thumb).contains(&k);
        buf.set_string(x, area.y + 1 + k, if on { "┃" } else { "│" }, DIM);
    }
}
