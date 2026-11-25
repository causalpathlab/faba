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
    cells, fill, values, width, GridArgs, Ink, Look, Refused, Side, BLOCK_EDGE, BLOCK_OFF,
    BLOCK_OFF_GRAY, BLOCK_UP, BLOCK_UP_GRAY, READ_OFF, READ_OFF_GRAY, READ_UP, READ_UP_GRAY,
    SIZE_BAR,
};
use crate::cluster::browse::go::Table;
use crate::cluster::browse::{Count, Modules};
use crate::figure::size::Page;
use crate::figure::{svg, Anchor, Canvas, FAINT, INK, MUTED};

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

/// `s` shortened with an ellipsis to fit `room` points of `f`-point text,
/// a character taken as wide as a capital (a term is often in capitals).
fn fit(s: &str, room: f64, f: f64) -> String {
    let fits = (room / (0.66 * f)).floor().max(0.0) as usize;
    match s.chars().count() > fits {
        true => (s.chars().take(fits.saturating_sub(1)))
            .chain(['…'])
            .collect(),
        false => s.to_string(),
    }
}

/// The widest the GO terms beside the grid are, in characters.
const GO_CHARS: f64 = 28.0;
/// A GO term's line, in font sizes.
const GO_LINE: f64 = 1.15;

/// The GO terms on each row: each term's name and q, up to `n`.
type Terms = Vec<Vec<(String, String)>>;

fn go_terms(table: &Table, rows: &[usize], n: usize) -> Terms {
    (rows.iter())
        .map(|&g| {
            let hits = table.of(g).map_or(&[][..], |t| &t.hits);
            (hits.iter().take(n))
                .map(|h| (h.named().to_string(), format!("{:.1e}", h.q)))
                .collect()
        })
        .collect()
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
    fn key(&self, c: &mut Canvas, (x, y): (f64, f64), f: f64, per_cell: bool) {
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
        c.text(at, y + s, "one scale a row", f, Anchor::Start, MUTED);
    }
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
    let (rows, values) = values(m, look);
    let all: Vec<f64> = values.iter().flatten().flatten().copied().collect();
    let diverging = !matches!(look, Look::Percent | Look::Bars);
    // The blocks' reads as bars: the look's own, or over another's.
    let bars = (look == Look::Bars || grid.block_bars).then(|| Bars::of(m, &rows, ink.gray));
    let top = match diverging {
        true => all.iter().fold(0.0f64, |a, v| a.max(v.abs())),
        false => all.iter().fold(0.0f64, |a, &v| a.max(v)),
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
    let names: Vec<String> = (rows.iter())
        .map(|&g| format!("M{} ({} genes)", g + 1, m.gene_sizes[g]))
        .collect();
    let (cols, n) = (m.columns.len(), rows.len());
    let (genes, cell_margin) = margins(m, &rows);
    let (left, over, under) = (
        strips(grid, Side::Left),
        strips(grid, Side::Top),
        strips(grid, Side::Bottom),
    );
    let per_cell = m.per_cell();
    let said: Vec<String> = (left.iter().map(|&s| genes.label(s, per_cell)))
        .chain((over.iter().chain(&under)).map(|&s| cell_margin.label(s, per_cell)))
        .collect();
    let label = (names.iter().chain(&said))
        .map(|n| width(n, f))
        .fold(0.0, f64::max)
        + f;
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
    let (nl, nt, nb) = (left.len() as f64, over.len() as f64, under.len() as f64);
    if !over.is_empty() {
        // From half a line under the title down, then the columns' names.
        head = ink.room(2.0 * f) + 2.1 * f + (nt - 1.0) * gap;
    }
    // The left margins' labels go under the grid, or over it when only the
    // bottom has margins; a line each.
    let labels_over = !left.is_empty() && over.is_empty() && !under.is_empty();
    if labels_over {
        head += ((nl - 1.0) * 1.1 * f - 0.6 * f).max(0.0);
    }
    foot += nb * gap;
    if !left.is_empty() && !labels_over {
        foot += nl * 1.1 * f;
    }
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
    let none = go.map_or(String::new(), |t| format!("none at q <= {}", t.q));
    let go_w = terms.as_ref().map(|terms| {
        let widest = (terms.iter().flatten())
            .map(|(name, q)| width(name, f) + width(q, f) + 0.5 * f)
            .chain([width(&none, f)])
            .fold(0.0, f64::max);
        widest.min(width("x", f) * GO_CHARS)
    });
    // Each row's shift range right of the grid, as the browser has it,
    // when the bars are drawn.
    let ranges: Vec<String> = match bars {
        Some(_) => (rows.iter())
            .map(|&g| format!("Δ{:.2}", m.shift_range(g)))
            .collect(),
        None => Vec::new(),
    };
    let range_w = (ranges.iter())
        .map(|s| 0.6 * f + width(s, f))
        .fold(0.0, f64::max);
    let beside = range_w + go_w.map_or(0.0, |w| 0.6 * f + w);
    let fixed = (label + nl * gap + beside + scale, head + foot);
    let units = (
        cols as f64 + nl * depth,
        n as f64 + (nt + nb) * depth + key_cells,
    );
    // Square boxes as large as the page allows (44 points without one).
    let (cell, (w, h)) = cells(size, 44.0, fixed, units)?;
    let d = depth * cell;
    // The grid under its top margins.
    let head = head + nt * d;
    let font = (cell * 0.28).min(f);
    let gx = label + nl * (d + gap);
    let right = gx + cell * cols as f64;
    let bottom = head + cell * n as f64;
    let (up_fill, off_fill) = match ink.gray {
        true => (READ_UP_GRAY, READ_OFF_GRAY),
        false => (READ_UP, READ_OFF),
    };
    Ok(svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        for (j, &role) in m.columns.iter().enumerate() {
            let x = gx + cell * (j as f64 + 0.5);
            c.text(
                x,
                head - f * 0.6,
                &format!("C{}", role + 1),
                f,
                Anchor::Middle,
                INK,
            );
        }
        for (i, (row, name)) in values.iter().zip(&names).enumerate() {
            let y = head + cell * i as f64;
            c.text(
                label - f * 0.6,
                y + cell / 2.0 + f / 2.5,
                name,
                f,
                Anchor::End,
                INK,
            );
            for (j, v) in row.iter().enumerate() {
                let x = gx + cell * j as f64;
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
                            let colour = fill(v * UNDER_BARS, top, diverging).0;
                            let (cx, cy) = (centre.0 - side / 2.0, centre.1 - side / 2.0);
                            c.rect(cx, cy, side, side, &colour);
                        }
                        (Look::Hinton, None) => {}
                        (_, v) => {
                            let colour = v.map_or("#ffffff".into(), |v| {
                                ink.fill(v * UNDER_BARS, top, diverging).0
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
                    let (colour, dark) = fill(*v, top, diverging);
                    let (cx, cy) = centre;
                    c.rect(cx - side / 2.0, cy - side / 2.0, side, side, &colour);
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
                let (colour, dark) = ink.fill(*v, top, diverging);
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
        // The gene modules' margins, nearest the grid first, each bar by
        // its row: sizes growing out from the grid, reads as mirror bars
        // (converted out, unconverted in) on one scale.
        let mut lines_at = match labels_over {
            true => head - 0.6 * f,
            false => bottom + nb * (d + gap) + f,
        };
        for (k, &strip) in left.iter().enumerate() {
            let x1 = gx - gap - k as f64 * (d + gap);
            for i in 0..n {
                let y = head + cell * i as f64 + 1.0;
                match strip {
                    Strip::Size => {
                        let most = genes.sizes.iter().copied().fold(0.0, f64::max);
                        let len = d * genes.sizes[i] / most.max(1e-12);
                        c.rect(x1 - len, y, len, cell - 2.0, SIZE_BAR);
                    }
                    Strip::Reads => {
                        let (zero, peak) = (x1 - d / 2.0, genes.peak().max(1e-12));
                        let r = genes.reads[i];
                        let (up, off) = (d / 2.0 * r.up() / peak, d / 2.0 * r.off() / peak);
                        c.rect(zero - up, y, up, cell - 2.0, up_fill);
                        c.rect(zero, y, off, cell - 2.0, off_fill);
                    }
                }
            }
            let said = genes.label(strip, per_cell);
            c.text(x1, lines_at, &said, f, Anchor::End, MUTED);
            lines_at += if labels_over { -1.1 * f } else { 1.1 * f };
        }
        // The cell modules' margins, over or under the grid, nearest it
        // first, each bar by its column, labelled on the left.
        let sides = [(&over, -1.0), (&under, 1.0)];
        for (strips, dir) in sides {
            for (k, &strip) in strips.iter().enumerate() {
                // The strip's edge nearer the grid, and its far one.
                let near = match dir < 0.0 {
                    true => head - 1.6 * f - k as f64 * (d + gap),
                    false => bottom + gap + k as f64 * (d + gap),
                };
                let (y0, y1) = match dir < 0.0 {
                    true => (near - d, near),
                    false => (near, near + d),
                };
                for j in 0..cols {
                    let x = gx + cell * j as f64 + 1.0;
                    match strip {
                        Strip::Size => {
                            let most = cell_margin.sizes.iter().copied().fold(0.0, f64::max);
                            let len = d * cell_margin.sizes[j] / most.max(1e-12);
                            let y = if dir < 0.0 { y1 - len } else { y0 };
                            c.rect(x, y, cell - 2.0, len, SIZE_BAR);
                        }
                        Strip::Reads => {
                            let (zero, peak) = ((y0 + y1) / 2.0, cell_margin.peak().max(1e-12));
                            let r = cell_margin.reads[j];
                            let (up, off) = (d / 2.0 * r.up() / peak, d / 2.0 * r.off() / peak);
                            c.rect(x, zero - up, cell - 2.0, up, up_fill);
                            c.rect(x, zero, cell - 2.0, off, off_fill);
                        }
                    }
                }
                let said = cell_margin.label(strip, per_cell);
                let y = (y0 + y1) / 2.0 + 0.35 * f;
                c.text(gx - gap, y, &said, f, Anchor::End, MUTED);
            }
        }
        // Each row's shift range right of the grid.
        for (i, said) in ranges.iter().enumerate() {
            let y = head + cell * (i as f64 + 0.5) + f / 2.5;
            c.text(right + 0.6 * f, y, said, f, Anchor::Start, MUTED);
        }
        // Each gene module's top GO terms on its row, right of the grid
        // and the ranges, as many lines as the row holds.
        if let (Some(terms), Some(go_w)) = (&terms, go_w) {
            let x = right + range_w + 0.6 * f;
            if !ink.panel_only {
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
            for (i, row) in terms.iter().enumerate() {
                let yc = head + cell * (i as f64 + 0.5);
                if row.is_empty() {
                    let said = fit(&none, go_w, f);
                    c.text(x, yc + 0.35 * f, &said, f, Anchor::Start, MUTED);
                    continue;
                }
                let held = (((cell - 2.0) / (GO_LINE * f)).floor() as usize).max(1);
                let shown = row.len().min(held);
                let first = yc - (shown as f64 - 1.0) * GO_LINE * f / 2.0 + 0.35 * f;
                for (k, (name, q)) in row.iter().take(shown).enumerate() {
                    let y = first + k as f64 * GO_LINE * f;
                    let room = go_w - width(q, f) - 0.5 * f;
                    c.text(x, y, &fit(name, room, f), f, Anchor::Start, INK);
                    c.text(x + go_w, y, q, f, Anchor::End, MUTED);
                }
            }
        }
        let low = if diverging { -top } else { 0.0 };
        let tick = |v: f64| match look {
            Look::Z | Look::Hinton => format!("{v:.1}"),
            _ => format!("{v:.2}"),
        };
        let below = (!left.is_empty() && !labels_over).then_some(nl * 1.1 * f);
        let ky = bottom + nb * (d + gap) + below.unwrap_or(0.0) + 0.6 * f;
        let x0 = 0.5 * f;
        if let Some(b) = bars.as_ref().filter(|_| bars_key) {
            // Under the colour key, or over the note, or in its place.
            let y = match (ink.panel_only, look) {
                (false, _) => h - 3.0 * f,
                (true, Look::Bars) => ky,
                // Under the Hinton key's squares and their labels too.
                (true, _) => ky + (2.3 * f).max(key_cells * (cell - 2.0) + 1.4 * f),
            };
            b.key(c, (x0, y), f, per_cell);
        }
        if small_key && look != Look::Bars {
            // Under everything, from the page's left: the colour scale
            // (labelled at its ends and middle), then for the Hinton look
            // squares at reference |z|, each its area as in the grid.
            let ends = [tick(low), tick((top + low) / 2.0), tick(top)];
            let said = ends.iter().map(|t| width(t, f)).sum::<f64>() + 2.0 * f;
            let rw = said.max(6.0 * f).min(right - x0);
            // Lightened as the blocks are under the bars.
            let shade = if bars.is_some() { UNDER_BARS } else { 1.0 };
            let steps = 24;
            for k in 0..steps {
                let v = low + (top - low) * (k as f64 + 0.5) / steps as f64;
                let colour = match look {
                    Look::Hinton => fill(v * shade, top, diverging).0,
                    _ => ink.fill(v * shade, top, diverging).0,
                };
                let step = rw / steps as f64;
                c.rect(x0 + step * k as f64, ky, step + 0.3, 0.7 * f, &colour);
            }
            let ty = ky + 1.7 * f;
            c.text(x0, ty, &ends[0], f, Anchor::Start, MUTED);
            c.text(x0 + rw / 2.0, ty, &ends[1], f, Anchor::Middle, MUTED);
            c.text(x0 + rw, ty, &ends[2], f, Anchor::End, MUTED);
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
                        c.rect(
                            cx - s / 2.0,
                            ky,
                            *s,
                            *s,
                            &fill(*v * shade, top, diverging).0,
                        );
                        c.text(cx, label_y, &format!("{v}"), f, Anchor::Middle, MUTED);
                        at += slot + 0.8 * f;
                    }
                }
            }
        }
        if ink.panel_only {
            return;
        }
        if look == Look::Bars {
            c.text(label, h - f, what, f, Anchor::Start, MUTED);
            return;
        }
        // The scale: a column of steps from the top value down to its other
        // end, labelled at both ends and the middle, right of the GO terms.
        let (sx, sy, sh, steps) = (right + beside + f, head, cell * n as f64, 40);
        for k in 0..steps {
            let v = top - (top - low) * (k as f64 + 0.5) / steps as f64;
            let y = sy + sh * k as f64 / steps as f64;
            c.rect(
                sx,
                y,
                f * 1.2,
                sh / steps as f64 + 0.5,
                &ink.fill(v, top, diverging).0,
            );
        }
        let mid = (top + low) / 2.0;
        for (v, y) in [
            (top, sy + f * 0.4),
            (mid, sy + sh / 2.0 + f * 0.35),
            (low, sy + sh),
        ] {
            c.text(
                sx + f * 1.6,
                y,
                &format!("{v:.2}"),
                f * 0.9,
                Anchor::Start,
                MUTED,
            );
        }
        c.text(label, h - f, what, f, Anchor::Start, MUTED);
    }))
}
