//! The modules grid: a box per block (or a Hinton square), the gene
//! modules' names and sizes on the left, the cell modules' on top, and
//! the colour scale on the right; with [`GridArgs`], the modules' sizes and
//! reads as bars on its margins (gene modules left, cell modules top or
//! bottom), each gene module's top GO terms on its right, a small key
//! under the plot alone, and the blocks' values at the font size.
//!
//! Every margin bar is `--marginal-size` of the grid's shorter side deep,
//! so the margins are laid out in cells, as the grid is: the cells and
//! the bars are as large as the page allows, the page kept as asked.

use super::{
    cells, fill, hinton_fill, values, width, GridArgs, Ink, Look, Refused, Side, BLOCK_EDGE,
    BLOCK_OFF, BLOCK_OFF_GRAY, BLOCK_UP, BLOCK_UP_GRAY, READ_OFF, READ_OFF_GRAY, READ_UP,
    READ_UP_GRAY, SIZE_BAR,
};
use crate::cluster::browse::go::{q_said, Table};
use crate::cluster::browse::view::{Ramp, Stretch};
use crate::cluster::browse::{Count, Modules};
use crate::figure::size::Page;
use crate::figure::{svg, Anchor, Canvas, FAINT, INK, MUTED};

/// The outline of each Hinton square: gray90, very thin.
const HINTON_EDGE: (&str, f64) = ("#e5e5e5", 0.1);

/// What a margin's bars show.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Strip {
    /// Each module's size: genes, or cells.
    Size,
    /// Each module's reads summed over its blocks, as mirror bars.
    Reads,
}

/// The strips on `side`, nearest the grid first: the sizes, then the
/// reads.
fn strips(grid: &GridArgs, side: Side) -> Vec<Strip> {
    let on = [
        (Strip::Size, &grid.size_bars),
        (Strip::Reads, &grid.read_bars),
    ];
    (on.into_iter())
        .filter(|(_, sides)| sides.contains(&side))
        .map(|(s, _)| s)
        .collect()
}

/// A count of reads in a few characters.
fn reads_said(n: f64) -> String {
    match n < 999.95 {
        true => format!("{n:.0}"),
        false => crate::tui::ui::compact(n),
    }
}

/// A margin's modules: each one's size, or its reads, in the grid's order.
struct Margin {
    sizes: Vec<f64>,
    reads: Vec<Count>,
    /// The largest size, as its label says it.
    size_label: String,
}

impl Margin {
    /// The taller side of the tallest of the reads, the top of their
    /// scale, as the browser's mirror bars take it.
    fn peak(&self) -> f64 {
        (self.reads.iter())
            .map(|r| r.up().max(r.off()))
            .fold(0.0, f64::max)
    }

    /// What `strip` is labelled with: its largest value.
    fn label(&self, strip: Strip, per_cell: bool) -> String {
        match (strip, per_cell) {
            (Strip::Size, _) => self.size_label.clone(),
            (Strip::Reads, false) => format!("{} reads", reads_said(self.peak())),
            (Strip::Reads, true) => format!("{:.2} UMIs/cell", self.peak()),
        }
    }
}

/// The gene modules' margin (the rows') and the cell modules' (the
/// columns'), the reads each module's blocks in the grid summed.
fn margins(m: &Modules, rows: &[usize]) -> (Margin, Margin) {
    let big = rows.iter().map(|&g| m.gene_sizes[g]).max().unwrap_or(0);
    let genes = Margin {
        sizes: rows.iter().map(|&g| m.gene_sizes[g] as f64).collect(),
        reads: (rows.iter())
            .map(|&g| m.block(g).iter().copied().sum())
            .collect(),
        size_label: format!("{big} genes"),
    };
    let sizes: Vec<f64> = m.columns.iter().map(|&c| m.cell_sizes[c]).collect();
    let most = sizes.iter().copied().fold(0.0, f64::max);
    let cells = Margin {
        sizes,
        reads: (0..m.columns.len())
            .map(|j| rows.iter().map(|&g| m.block(g)[j]).sum())
            .collect(),
        size_label: format!("{} cells", m.cells_label(most)),
    };
    (genes, cells)
}

/// Gene module `g`'s name on the grid: `M5 (420 genes)`, or with
/// `--short-module-names` `M5 (420)`.
pub(super) fn module_name(m: &Modules, g: usize, grid: &GridArgs) -> String {
    let unit = if grid.short_module_names {
        ""
    } else {
        " genes"
    };
    format!("M{} ({}{unit})", g + 1, m.gene_sizes[g])
}

/// `s`'s width in `f`-point GO text.
pub(in crate::cluster::browse) fn go_width(s: &str, f: f64) -> f64 {
    crate::figure::set_width(s, f)
}

/// `s` as it fits `room` points of `f`-point GO text: whole, or cut at the
/// end of a word (a space, `_` or `-` after it) with an ellipsis, cut
/// within a word only when no word's end keeps half of what fits.
fn shortened(s: &str, room: f64, f: f64) -> String {
    if go_width(s, f) <= room {
        return s.to_string();
    }
    let chars: Vec<char> = s.chars().collect();
    let mut used = go_width("…", f);
    let fits = (chars.iter())
        .take_while(|&&c| {
            used += crate::figure::glyph_width(c, f);
            used <= room
        })
        .count();
    let word = (fits / 2..=fits)
        .rev()
        .find(|&i| i > 0 && matches!(chars[i], ' ' | '_' | '-'));
    let kept: String = chars[..word.unwrap_or(fits)].iter().collect();
    kept.trim_end_matches([' ', '_', '-']).to_string() + "…"
}

/// The widest the GO terms beside the grid are kept, in characters, when
/// the page leaves less beside the grid.
const GO_CHARS: f64 = 28.0;
/// A GO term's line, in font sizes; and the gap kept between one row's
/// lines and the next's, at least.
const GO_LINE: f64 = 1.15;
const GO_GAP: f64 = 0.5;
/// A line's capitals' height, in font sizes.
const CAP: f64 = 0.75;

/// The GO terms on each row: each term's name and the q after it.
type Terms = Vec<Vec<(String, String)>>;

fn go_terms(table: &Table, rows: &[usize], n: usize) -> Terms {
    (rows.iter())
        .map(|&g| {
            let picked = table.of(g).map_or(Vec::new(), |t| t.representative(n));
            (picked.into_iter())
                .map(|h| (h.named().to_string(), format!(" (q {})", q_said(h.q))))
                .collect()
        })
        .collect()
}

/// How many of `n` lines a row `cell` points deep holds in `f`-point text,
/// and their spacing: each line no less than the font, and [`GO_GAP`]
/// between the row's lines and the next row's, so that rows never run
/// into each other; fewer lines where `n` do not fit; one at least.
fn go_lines(cell: f64, f: f64, n: usize) -> (usize, f64) {
    let room = cell - GO_GAP * f;
    let held = (room / f).floor().max(1.0) as usize;
    let shown = n.min(held).max(1);
    let spacing = match shown {
        1 => 0.0,
        k => (GO_LINE * f).min(room / k as f64),
    };
    (shown, spacing)
}

/// The reference |z| values the Hinton key's squares show: up to three
/// round values no more than `top`, the largest last.
fn references(top: f64) -> Vec<f64> {
    const NICE: [f64; 12] = [0.1, 0.2, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0];
    let under: Vec<f64> = NICE.into_iter().filter(|&v| v <= top).collect();
    match under.len() {
        0 if top > 0.0 => vec![top],
        n => under[n.saturating_sub(3)..].to_vec(),
    }
}

/// A block's value as `--values` prints it: z to a tenth, the shift to a
/// hundredth, the % (or UMIs per cell) as the grid writes it.
fn value_text(m: &Modules, look: Look, v: f64, g: usize, j: usize) -> String {
    match look {
        Look::Z | Look::Hinton => format!("{v:+.1}"),
        Look::Shift => format!("{v:+.2}"),
        Look::Percent | Look::Bars => m.block(g)[j].number().unwrap_or_default(),
    }
}

/// How far toward its full colour a block under the bars is drawn.
const UNDER_BARS: f64 = 0.4;

/// A block's bars: their width a share of the cell, their full height a
/// share of half of it, so a row's bars never touch the next row's.
const BAR_WIDTH: f64 = 0.4;
const BAR_REACH: f64 = 0.92;

/// The blocks' mirror bars, as the browser draws them: converted reads
/// rising above a zero line, unconverted hanging below, each row on its
/// own scale, as the browser's rows are: its `peaks` reads (the taller
/// side of its tallest) filling a side.
struct Bars {
    peaks: Vec<f64>,
    up: &'static str,
    off: &'static str,
}

impl Bars {
    /// The bars over `rows`, in grays or the accent.
    fn of(m: &Modules, rows: &[usize], gray: bool) -> Self {
        let peaks = (rows.iter())
            .map(|&g| super::super::view::peak(m.block(g)))
            .collect();
        let (up, off) = match gray {
            true => (BLOCK_UP_GRAY, BLOCK_OFF_GRAY),
            false => (BLOCK_UP, BLOCK_OFF),
        };
        Self { peaks, up, off }
    }

    /// The bar width in a `cell`.
    fn width(cell: f64) -> f64 {
        BAR_WIDTH * (cell - 2.0)
    }

    /// Block `r`'s bar in row `i`, in the cell at `(x, y)`, `cell` points
    /// a side, its left edge at `bx`, its zero line across the cell's
    /// middle.
    fn draw(&self, c: &mut Canvas, (i, r): (usize, Count), (bx, y): (f64, f64), cell: f64) {
        let (bw, zero) = (Self::width(cell), y + cell / 2.0);
        let side = BAR_REACH * (cell - 2.0) / 2.0 / self.peaks[i].max(1e-12);
        let (up, off) = (side * r.up(), side * r.off());
        c.rect(bx, zero - up, bw, up, self.up);
        c.outlined((bx, zero, bw, off), self.off, BLOCK_EDGE);
    }

    /// The key to the bars from `(x, y)`: each side's swatch and what it
    /// is, and their scale.
    /// Transposed, a gene module's bars are a column's.
    fn key(&self, c: &mut Canvas, (x, y): (f64, f64), f: f64, (per_cell, tr): (bool, bool)) {
        let s = 0.8 * f;
        let mut at = x;
        let sides: &[(&str, bool)] = match per_cell {
            true => &[("UMIs per cell", true)],
            false => &[("converted", true), ("unconverted", false)],
        };
        for &(said, up) in sides {
            match up {
                true => c.rect(at, y, s, s, self.up),
                false => c.outlined((at, y, s, s), self.off, BLOCK_EDGE),
            }
            at += s + 0.3 * f;
            c.text(at, y + s, said, f, Anchor::Start, MUTED);
            at += width(said, f) + 0.8 * f;
        }
        let scale = if tr {
            "one scale a column"
        } else {
            "one scale a row"
        };
        c.text(at, y + s, scale, f, Anchor::Start, MUTED);
    }
}

/// Where a drawn grid put what a screen places beside it: the left edge
/// of the GO terms, when they are drawn.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Placed {
    pub go_x: Option<f64>,
    /// The top of the grid's first row, and the rows' pitch (a cell).
    pub rows_top: f64,
    pub pitch: f64,
    /// How much of the page's width the grid takes, from its left.
    pub width: f64,
}

/// The figure of `m` under `look` with `grid`'s extras, `go` the terms for
/// `--go-labels`, titled `title`, laid out for `size`.
pub(super) fn heatmap(
    m: &Modules,
    look: Look,
    grid: &GridArgs,
    go: Option<&Table>,
    title: &str,
    size: Page,
    ink: Ink,
) -> Result<String, Refused> {
    placed(m, look, grid, go, title, size, ink).map(|(svg, _)| svg)
}

/// [`heatmap`], and where it put the GO terms.
pub(super) fn placed(
    m: &Modules,
    look: Look,
    grid: &GridArgs,
    go: Option<&Table>,
    title: &str,
    size: Page,
    ink: Ink,
) -> Result<(String, Placed), Refused> {
    let (rows, values) = values(m, look);
    let all: Vec<f64> = values.iter().flatten().flatten().copied().collect();
    let diverging = !matches!(look, Look::Percent | Look::Bars);
    // The blocks' reads as bars: the look's own, or over another's.
    let bars = (look == Look::Bars || grid.block_bars).then(|| Bars::of(m, &rows, ink.gray));
    let top = match diverging {
        true => all.iter().fold(0.0f64, |a, v| a.max(v.abs())),
        false => all.iter().fold(0.0f64, |a, &v| a.max(v)),
    };
    // The colours' top: the largest value, or with `--grid-cap` that
    // quantile of the blocks' values (|values| for a diverging look), the
    // blocks past it at full colour; sizes keep the largest.
    let shade_top = grid
        .cap
        .map_or(top, |q| quantile(&all, q, diverging).min(top));
    let capped = shade_top < top;
    // The Hinton squares coloured as the browser's are, full red or blue
    // at its top |z|; in grays as the other looks.
    let z_top = Stretch::of(rows.iter().map(|&g| m.block(g)), Ramp::Sqrt).z_top;
    let z_top = if capped { z_top.min(shade_top) } else { z_top };
    let square_fill = |v: f64| match ink.gray {
        true => fill(v, shade_top, diverging),
        false => hinton_fill(v, z_top),
    };
    let what = match (look, m.per_cell()) {
        (Look::Shift, true) => "log rate",
        (Look::Shift, false) => "shift (logit, over each site's baseline)",
        (Look::Percent, true) => "UMIs per cell",
        (Look::Percent, false) => "methylated %",
        (Look::Z, _) => "z across the cell modules",
        (Look::Hinton, _) => "z across the cell modules; square area by |z|",
        (Look::Bars, _) => "reads in each block",
    };
    let f = ink.f;
    // Transposed, the cell modules are the rows and the gene modules the
    // columns: every block, name, margin and label goes to the other axis.
    let tr = grid.transpose;
    let names: Vec<String> = rows.iter().map(|&g| module_name(m, g, grid)).collect();
    let cell_names: Vec<String> = (m.columns.iter()).map(|&c| format!("C{}", c + 1)).collect();
    let (cols, n) = (m.columns.len(), rows.len());
    // The grid's columns and rows as drawn.
    let (across, down) = match tr {
        false => (cols, n),
        true => (n, cols),
    };
    let (genes, cell_margin) = margins(m, &rows);
    let (left, over, under) = (
        strips(grid, Side::Left),
        strips(grid, Side::Top),
        strips(grid, Side::Bottom),
    );
    // The strips drawn by the grid's rows (left of it), over its columns
    // and under them, and the modules each strip is of. Transposed, the
    // gene modules' go over the columns and the cell modules' (top or
    // bottom) left of the rows.
    let cell_strips: Vec<Strip> = over.iter().chain(&under).copied().collect();
    let (side_strips, top_strips, bottom_strips): (&[Strip], &[Strip], &[Strip]) = match tr {
        false => (&left, &over, &under),
        true => (&cell_strips, &left, &[]),
    };
    let (side_margin, cross_margin) = match tr {
        false => (&genes, &cell_margin),
        true => (&cell_margin, &genes),
    };
    let per_cell = m.per_cell();
    let said: Vec<String> = (left.iter().map(|&s| genes.label(s, per_cell)))
        .chain(cell_strips.iter().map(|&s| cell_margin.label(s, per_cell)))
        .collect();
    let row_names = if tr { &cell_names } else { &names };
    let label = (row_names.iter().chain(&said))
        .map(|n| width(n, f))
        .fold(0.0, f64::max)
        + f;
    // A transposed grid's gene module names, turned, over its columns.
    let names_h = match tr {
        true => names.iter().map(|s| width(s, f)).fold(0.0, f64::max) + 0.8 * f,
        false => 0.0,
    };
    // The plot alone: no title over it, no scale beside it, no note under;
    // the bars alone have no colour scale.
    let (mut head, scale, mut foot) = (
        2.0 * f + ink.room(2.0 * f),
        if look == Look::Bars {
            0.0
        } else {
            ink.room(6.0 * f)
        },
        f + ink.room(2.0 * f),
    );
    // Each margin is `depth` cells deep, a gap from the grid and the next.
    let gap = 0.4 * f;
    let depth = grid.marginal_size * cols.min(n).max(1) as f64;
    let (nl, nt, nb) = (
        side_strips.len() as f64,
        top_strips.len() as f64,
        bottom_strips.len() as f64,
    );
    if tr {
        // The title's room, the names turned, then the strips over the
        // columns, a gap from the grid and each other.
        head = ink.room(2.0 * f) + names_h + nt * gap;
    } else if nt > 0.0 {
        // From half a line under the title down, then the columns' names.
        head = ink.room(2.0 * f) + 2.1 * f + (nt - 1.0) * gap;
    }
    // The left margins' labels go under the grid, or over it when only the
    // bottom has margins; a line each.
    let labels_over = !tr && nl > 0.0 && nt == 0.0 && nb > 0.0;
    if labels_over {
        head += ((nl - 1.0) * 1.1 * f - 0.6 * f).max(0.0);
    }
    foot += nb * gap;
    // The key under the plot alone: a colour scale, and for the Hinton
    // look squares at reference |z|, their sides in cells.
    let small_key = ink.panel_only && grid.key;
    let refs = match look {
        Look::Hinton if small_key => references(top),
        _ => Vec::new(),
    };
    let key_cells = refs.last().map_or(0.0, |v| (v / top.max(1e-12)).sqrt());
    if small_key && look != Look::Bars {
        foot += 2.3 * f;
    }
    // The bars' key, a line: under the plot alone with `--key`, else
    // over the note.
    let bars_key = bars.is_some() && (small_key || !ink.panel_only);
    if bars_key {
        foot += 1.4 * f;
    }
    let terms = grid.go_labels.zip(go).map(|(k, t)| go_terms(t, &rows, k));
    // A row with no term: blank, or saying so with `--go-none`.
    let none = (go.filter(|_| grid.go_none)).map(|t| format!("none at q <= {}", t.q));
    // The width the labels need, whole.
    let go_need = terms.as_ref().map(|terms| {
        (terms.iter().flatten())
            .map(|(name, q)| go_width(name, f) + go_width(q, f))
            .chain(none.iter().map(|s| go_width(s, f)))
            .fold(0.0, f64::max)
    });
    // Each gene module's shift range beside its row (under its column,
    // turned, when transposed), as the browser has it, when the bars are
    // drawn (or the screen asks for it).
    let ranges: Vec<String> = match bars.is_some() || grid.ranges {
        true => (rows.iter())
            .map(|&g| format!("Δ{:.2}", m.shift_range(g)))
            .collect(),
        false => Vec::new(),
    };
    let range_w = (ranges.iter())
        .map(|s| 0.6 * f + width(s, f))
        .fold(0.0, f64::max);
    // The grid laid out with the labels' width kept to `GO_CHARS` at most;
    // then the labels take what they need of the room the page leaves.
    let kept = go_need.map(|w| w.min(width("x", f) * GO_CHARS));
    let go_band = kept.map_or(0.0, |w| 0.6 * f + w);
    let go_top = tr && grid.go_side == Some(super::GoSide::Top);
    // Under a transposed grid: its ranges and GO terms under the columns,
    // beside the left strips' labels, a line each.
    let under_cols = range_w + if go_top { 0.0 } else { go_band };
    let lines_h = nl * 1.1 * f;
    match tr {
        true => {
            foot += under_cols.max(lines_h);
            if go_top {
                head += go_band;
            }
        }
        false if nl > 0.0 && !labels_over => foot += lines_h,
        false => {}
    }
    let beside = match tr {
        true => 0.0,
        false => range_w + go_band,
    };
    let fixed = (label + nl * gap + beside + scale, head + foot);
    let units = (
        across as f64 + nl * depth,
        down as f64 + (nt + nb) * depth + key_cells,
    );
    // Square boxes as large as the page allows (44 points without one).
    let (cell, (w, h)) = cells(size, 44.0, fixed, units)?;
    let left_over = match tr {
        false => w - fixed.0 - cell * units.0,
        true => h - fixed.1 - cell * units.1,
    };
    let go_w = (go_need.zip(kept)).map(|(need, kept)| need.min(kept + left_over.max(0.0)));
    let grown = go_w.zip(kept).map_or(0.0, |(w, k)| w - k);
    let beside = match tr {
        true => 0.0,
        false => range_w + go_w.map_or(0.0, |w| 0.6 * f + w),
    };
    let d = depth * cell;
    // The grid under its top margins (and a transposed grid's GO terms
    // over them, as wide as they grew).
    let head = head + nt * d + if go_top { grown } else { 0.0 };
    let font = (cell * 0.28).min(f);
    let gx = label + nl * (d + gap);
    let right = gx + cell * across as f64;
    let bottom = head + cell * down as f64;
    // Row `i`'s (gene module's) block in column `j` (cell module's), its
    // top left corner.
    let at = |i: usize, j: usize| match tr {
        false => (gx + cell * j as f64, head + cell * i as f64),
        true => (gx + cell * i as f64, head + cell * j as f64),
    };
    let (up_fill, off_fill) = match ink.gray {
        true => (READ_UP_GRAY, READ_OFF_GRAY),
        false => (READ_UP, READ_OFF),
    };
    let placed = Placed {
        go_x: go_w.filter(|_| !tr).map(|_| right + range_w + 0.6 * f),
        rows_top: head,
        pitch: cell,
        width: (fixed.0 + cell * units.0 + if tr { 0.0 } else { grown }).min(w),
    };
    let drawn = svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        // The cell modules' names over the columns, or left of the rows.
        for (j, said) in cell_names.iter().enumerate() {
            let (x, y) = at(0, j);
            match tr {
                false => {
                    let x = x + cell / 2.0;
                    c.text(x, head - f * 0.6, said, f, Anchor::Middle, INK);
                }
                true => {
                    let y = y + cell / 2.0 + f / 2.5;
                    c.text(label - f * 0.6, y, said, f, Anchor::End, INK);
                }
            }
        }
        // Over the top strips, the transposed grid's gene modules' names.
        let strips_top = head - nt * (d + gap);
        for (i, (row, name)) in values.iter().zip(&names).enumerate() {
            let (x, y) = at(i, 0);
            match tr {
                false => c.text(
                    label - f * 0.6,
                    y + cell / 2.0 + f / 2.5,
                    name,
                    f,
                    Anchor::End,
                    INK,
                ),
                true => {
                    let x = x + cell / 2.0 + f / 2.5;
                    c.rotated(x, strips_top - 0.4 * f, name, f, Anchor::Start, INK);
                }
            }
            for (j, v) in row.iter().enumerate() {
                let (x, y) = at(i, j);
                let centre = (x + cell / 2.0, y + cell / 2.0);
                // A value at the font size, where it fits in the cell.
                let asked = (*v)
                    .filter(|_| grid.values)
                    .map(|v| value_text(m, look, v, rows[i], j))
                    .filter(|t| width(t, f) <= cell - 3.0 && 1.2 * f <= cell - 2.0);
                if let Some(b) = &bars {
                    // Under the bars, the look's own block lightened; under
                    // the bars alone, their zero line.
                    match (look, v) {
                        (Look::Bars, _) => {
                            let (x0, x1) = (x + 1.0, x + cell - 1.0);
                            c.line(x0, centre.1, x1, centre.1, FAINT, 0.4);
                        }
                        (Look::Hinton, Some(v)) => {
                            let side = (cell - 2.0) * (v.abs() / top.max(1e-12)).min(1.0).sqrt();
                            let colour = square_fill(v * UNDER_BARS).0;
                            let (cx, cy) = (centre.0 - side / 2.0, centre.1 - side / 2.0);
                            c.rect(cx, cy, side, side, &colour);
                        }
                        (Look::Hinton, None) => {}
                        (_, v) => {
                            let colour = v.map_or("#ffffff".into(), |v| {
                                ink.fill(v * UNDER_BARS, shade_top, diverging).0
                            });
                            c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, &colour);
                        }
                    }
                    // The bar centred, or left of its value where that fits
                    // beside it.
                    let bw = Bars::width(cell);
                    let beside = asked.filter(|t| width(t, f) + bw + 2.5 <= cell - 2.0);
                    let bx = match beside {
                        Some(_) => x + 1.5,
                        None => centre.0 - bw / 2.0,
                    };
                    b.draw(c, (i, m.block(rows[i])[j]), (bx, y), cell);
                    if let Some(t) = beside {
                        let tx = x + cell - 1.5;
                        c.text(tx, centre.1 + f / 2.5, &t, f, Anchor::End, INK);
                    }
                    continue;
                }
                if look == Look::Hinton {
                    // A true square, centred, none for a block without reads.
                    let Some(v) = v else { continue };
                    let side = (cell - 2.0) * (v.abs() / top.max(1e-12)).min(1.0).sqrt();
                    let (colour, dark) = square_fill(*v);
                    let (cx, cy) = centre;
                    let at = (cx - side / 2.0, cy - side / 2.0, side, side);
                    c.outlined(at, &colour, HINTON_EDGE);
                    if let Some(t) = asked {
                        let under = dark && side >= width(&t, f) + 2.0 && side >= 1.2 * f;
                        let text_ink = if under { "#ffffff" } else { INK };
                        c.text(cx, cy + f / 2.5, &t, f, Anchor::Middle, text_ink);
                    }
                    continue;
                }
                let Some(v) = v else {
                    c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, "#ffffff");
                    continue;
                };
                let (colour, dark) = ink.fill(*v, shade_top, diverging);
                c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, &colour);
                let text_ink = if dark { "#ffffff" } else { INK };
                let (text, size) = match (grid.values, asked) {
                    (true, Some(t)) => (t, f),
                    (true, None) => continue,
                    (false, _) => (
                        match look {
                            Look::Shift | Look::Z | Look::Hinton => format!("{v:+.2}"),
                            // As the grid writes it: a whole %, or UMIs per cell.
                            Look::Percent | Look::Bars => {
                                m.block(rows[i])[j].number().unwrap_or_default()
                            }
                        },
                        font,
                    ),
                };
                c.text(
                    centre.0,
                    centre.1 + size / 2.5,
                    &text,
                    size,
                    Anchor::Middle,
                    text_ink,
                );
            }
        }
        // The margins by the rows (left of the grid), nearest the grid
        // first, each bar by its row: sizes growing out from the grid,
        // reads as mirror bars (converted out, unconverted in) on one scale.
        let mut lines_at = match labels_over {
            true => head - 0.6 * f,
            false => bottom + nb * (d + gap) + f,
        };
        for (k, &strip) in side_strips.iter().enumerate() {
            let x1 = gx - gap - k as f64 * (d + gap);
            for i in 0..down {
                let y = head + cell * i as f64 + 1.0;
                match strip {
                    Strip::Size => {
                        let most = side_margin.sizes.iter().copied().fold(0.0, f64::max);
                        let len = d * side_margin.sizes[i] / most.max(1e-12);
                        c.rect(x1 - len, y, len, cell - 2.0, SIZE_BAR);
                    }
                    Strip::Reads => {
                        let (zero, peak) = (x1 - d / 2.0, side_margin.peak().max(1e-12));
                        let r = side_margin.reads[i];
                        let (up, off) = (d / 2.0 * r.up() / peak, d / 2.0 * r.off() / peak);
                        c.rect(zero - up, y, up, cell - 2.0, up_fill);
                        c.rect(zero, y, off, cell - 2.0, off_fill);
                    }
                }
            }
            let said = side_margin.label(strip, per_cell);
            c.text(x1, lines_at, &said, f, Anchor::End, MUTED);
            lines_at += if labels_over { -1.1 * f } else { 1.1 * f };
        }
        // The margins by the columns, over or under the grid, nearest it
        // first, each bar by its column, labelled on the left.
        let sides = [(&top_strips, -1.0), (&bottom_strips, 1.0)];
        for (strips, dir) in sides {
            for (k, &strip) in strips.iter().enumerate() {
                // The strip's edge nearer the grid (past the columns' names
                // over an upright grid), and its far one.
                let near = match dir < 0.0 {
                    true if tr => head - gap - k as f64 * (d + gap),
                    true => head - 1.6 * f - k as f64 * (d + gap),
                    false => bottom + gap + k as f64 * (d + gap),
                };
                let (y0, y1) = match dir < 0.0 {
                    true => (near - d, near),
                    false => (near, near + d),
                };
                for j in 0..across {
                    let x = gx + cell * j as f64 + 1.0;
                    match strip {
                        Strip::Size => {
                            let most = cross_margin.sizes.iter().copied().fold(0.0, f64::max);
                            let len = d * cross_margin.sizes[j] / most.max(1e-12);
                            let y = if dir < 0.0 { y1 - len } else { y0 };
                            c.rect(x, y, cell - 2.0, len, SIZE_BAR);
                        }
                        Strip::Reads => {
                            let (zero, peak) = ((y0 + y1) / 2.0, cross_margin.peak().max(1e-12));
                            let r = cross_margin.reads[j];
                            let (up, off) = (d / 2.0 * r.up() / peak, d / 2.0 * r.off() / peak);
                            c.rect(x, zero - up, cell - 2.0, up, up_fill);
                            c.rect(x, zero, cell - 2.0, off, off_fill);
                        }
                    }
                }
                let said = cross_margin.label(strip, per_cell);
                let y = (y0 + y1) / 2.0 + 0.35 * f;
                c.text(gx - gap, y, &said, f, Anchor::End, MUTED);
            }
        }
        // Each gene module's shift range right of its row, or under its
        // column, turned.
        for (i, said) in ranges.iter().enumerate() {
            let (x, y) = at(i, 0);
            match tr {
                false => {
                    let y = y + cell / 2.0 + f / 2.5;
                    c.text(right + 0.6 * f, y, said, f, Anchor::Start, MUTED);
                }
                true => {
                    let x = x + cell / 2.0 + f / 2.5;
                    c.rotated(x, bottom + 0.6 * f, said, f, Anchor::End, MUTED);
                }
            }
        }
        // Each gene module's top GO terms by its row, right of the grid
        // and the ranges (under or over its column, turned, when
        // transposed), as many lines as the row holds.
        if let (Some(terms), Some(go_w)) = (&terms, go_w) {
            let x = right + range_w + 0.6 * f;
            if !ink.panel_only && !tr {
                let q = go.map_or(0.0, |t| t.q);
                c.text(
                    x,
                    head - f * 0.6,
                    &format!("GO (q <= {q})"),
                    f,
                    Anchor::Start,
                    MUTED,
                );
            }
            // A transposed grid's terms: from under the ranges down, or
            // from over the names up.
            let (go_y, go_anchor) = match go_top {
                false => (bottom + range_w + 0.6 * f, Anchor::End),
                true => (strips_top - names_h - 0.2 * f, Anchor::Start),
            };
            for (i, row) in terms.iter().enumerate() {
                let (cx, cy) = at(i, 0);
                let mid = match tr {
                    false => cy + cell / 2.0,
                    true => cx + cell / 2.0,
                };
                let line = |c: &mut Canvas, along: f64, said: &str, ink: &str| match tr {
                    false => c.text(x, along, said, f, Anchor::Start, ink),
                    true => c.rotated(along, go_y, said, f, go_anchor, ink),
                };
                if row.is_empty() {
                    if let Some(none) = &none {
                        let said = shortened(none, go_w, f);
                        line(c, mid + 0.35 * f, &said, MUTED);
                    }
                    continue;
                }
                // A line a term, the lines' capitals centred on the row.
                let (shown, spacing) = go_lines(cell, f, row.len());
                let block = (shown - 1) as f64 * spacing + CAP * f;
                let first = mid - block / 2.0 + CAP * f;
                for (k, (name, q)) in row.iter().take(shown).enumerate() {
                    let along = first + k as f64 * spacing;
                    let name = shortened(name, go_w - go_width(q, f), f);
                    line(c, along, &format!("{name}{q}"), INK);
                }
            }
        }
        let low = if diverging { -shade_top } else { 0.0 };
        let tick = |v: f64| match look {
            Look::Z | Look::Hinton => format!("{v:.1}"),
            _ => format!("{v:.2}"),
        };
        // A capped scale's ends say the blocks past them share their colour.
        let end = |v: f64, said: String| match (capped, v) {
            (true, v) if v > 0.0 && v >= shade_top => format!("≥{said}"),
            (true, v) if v < 0.0 && v <= -shade_top => format!("≤{said}"),
            _ => said,
        };
        let below = match tr {
            true => under_cols.max(lines_h) + if go_top { 0.0 } else { grown },
            false if nl > 0.0 && !labels_over => lines_h,
            false => 0.0,
        };
        let ky = bottom + nb * (d + gap) + below + 0.6 * f;
        let x0 = 0.5 * f;
        if let Some(b) = bars.as_ref().filter(|_| bars_key) {
            // Under the colour key, or over the note, or in its place.
            let y = match (ink.panel_only, look) {
                (false, _) => h - 3.0 * f,
                (true, Look::Bars) => ky,
                // Under the Hinton key's squares and their labels too.
                (true, _) => ky + (2.3 * f).max(key_cells * (cell - 2.0) + 1.4 * f),
            };
            b.key(c, (x0, y), f, (per_cell, tr));
        }
        if small_key && look != Look::Bars {
            // Under everything, from the page's left: the colour scale
            // (labelled at its ends and middle), then for the Hinton look
            // squares at reference |z|, each its area as in the grid.
            let ends = [
                end(low, tick(low)),
                tick((shade_top + low) / 2.0),
                end(shade_top, tick(shade_top)),
            ];
            let said = ends.iter().map(|t| width(t, f)).sum::<f64>() + 2.0 * f;
            let rw = said.max(6.0 * f).min(right - x0);
            // Lightened as the blocks are under the bars.
            let shade = if bars.is_some() { UNDER_BARS } else { 1.0 };
            let key_fill = |along: f64| {
                let v = low + (shade_top - low) * along;
                match look {
                    Look::Hinton => square_fill(v * shade).0,
                    _ => ink.fill(v * shade, shade_top, diverging).0,
                }
            };
            let ends = [ends[0].as_str(), ends[1].as_str(), ends[2].as_str()];
            super::key_strip(c, (x0, ky, rw), f, 24, key_fill, ends);
            if !refs.is_empty() {
                let sides: Vec<f64> = (refs.iter())
                    .map(|v| (cell - 2.0) * (v / top.max(1e-12)).sqrt())
                    .collect();
                let room = |v: &f64, s: &f64| s.max(width(&format!("{v}"), f)) + 0.8 * f;
                let start = x0 + rw + 1.5 * f;
                let lead = width("|z|", f) + 0.6 * f;
                // As many as fit, the largest kept.
                let mut from = 0;
                while from < refs.len() {
                    let need: f64 = (refs[from..].iter().zip(&sides[from..]))
                        .map(|(v, s)| room(v, s))
                        .sum();
                    if start + lead + need <= right {
                        break;
                    }
                    from += 1;
                }
                if from < refs.len() {
                    let tallest = sides.last().copied().unwrap_or(0.0);
                    let label_y = ky + tallest + 1.0 * f;
                    c.text(start, label_y, "|z|", f, Anchor::Start, MUTED);
                    let mut at = start + lead;
                    for (v, s) in refs[from..].iter().zip(&sides[from..]) {
                        let slot = room(v, s) - 0.8 * f;
                        let cx = at + slot / 2.0;
                        c.rect(cx - s / 2.0, ky, *s, *s, &square_fill(*v * shade).0);
                        c.text(cx, label_y, &format!("{v}"), f, Anchor::Middle, MUTED);
                        at += slot + 0.8 * f;
                    }
                }
            }
        }
        if ink.panel_only {
            return;
        }
        // A transposed grid's GO terms' q said in the note.
        let what = match (tr, &terms, go) {
            (true, Some(_), Some(t)) => format!("{what}; GO terms at q <= {}", t.q),
            _ => what.to_string(),
        };
        // A transposed grid is narrow: its note kept to the page.
        let what = match tr {
            true => crate::figure::fit_text(&what, w - label - 0.5 * f, f),
            false => what,
        };
        if look == Look::Bars {
            c.text(label, h - f, &what, f, Anchor::Start, MUTED);
            return;
        }
        // The scale: a column of steps from the top value down to its other
        // end, labelled at both ends and the middle, right of the GO terms.
        let scale_fill = |along: f64| {
            let v = shade_top - (shade_top - low) * along;
            match look {
                Look::Hinton => square_fill(v).0,
                _ => ink.fill(v, shade_top, diverging).0,
            }
        };
        let mid = (shade_top + low) / 2.0;
        let said = [
            end(shade_top, format!("{shade_top:.2}")),
            end(mid, format!("{mid:.2}")),
            end(low, format!("{low:.2}")),
        ];
        let at = (right + beside + f, head, cell * down as f64);
        let said = [said[0].as_str(), said[1].as_str(), said[2].as_str()];
        super::scale_column(c, at, f, scale_fill, said);
        c.text(label, h - f, &what, f, Anchor::Start, MUTED);
    });
    Ok((drawn, placed))
}

/// The `q` quantile (in (0, 1]) of `values`, of their sizes when
/// `diverging`: the least that `q` of them are no more than.
pub(super) fn quantile(values: &[f64], q: f64, diverging: bool) -> f64 {
    let mut v: Vec<f64> = (values.iter())
        .map(|&v| if diverging { v.abs() } else { v })
        .collect();
    if v.is_empty() {
        return 0.0;
    }
    let at = ((q * v.len() as f64).ceil() as usize).clamp(1, v.len()) - 1;
    *v.select_nth_unstable_by(at, f64::total_cmp).1
}
