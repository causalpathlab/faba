//! A fit's figures, as `--figure` writes them (a PDF and a PNG) and MCP's
//! `cluster_view` returns them: the modules grid (gene modules against cell
//! modules, each block coloured by its fitted shift, a count fit's log rate,
//! its methylated %, UMIs per cell, or that standardized across the cell
//! modules), the trace (each restart's or chain's ELBO or log joint by
//! step), and the cell modules against the batches. Each is laid out for a
//! size when one is given (inches), its text at a size that reads in print.

use super::order::standardized;
use super::Modules;
use crate::figure::size::{self, Page};
use crate::figure::{esc, svg, Anchor, Canvas, FAINT, INK, MUTED};

/// Which figure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Plot {
    /// Gene modules against cell modules.
    Modules,
    /// Each restart's (or chain's) ELBO (or log joint) by step.
    Trace,
    /// Cell modules against the batches: each one's share in each.
    Batches,
}

/// A figure's text size in points: 9 at its own size, less on a small
/// page (to 5, the least that reads in print), by the sides asked for.
fn font(page: Page) -> f64 {
    let side = [page.0, page.1]
        .into_iter()
        .flatten()
        .fold(f64::INFINITY, f64::min);
    match side.is_finite() {
        true => (side / 40.0).clamp(5.0, 9.0),
        false => 9.0,
    }
}

/// The least cell side a grid is drawn with, in points: its values need
/// the room.
const CELL_LEAST: f64 = 12.0;

/// Square cells, `cols` by `rows`, beside `fixed` points across and down,
/// `own` points each at the grid's own size: a side not asked for follows
/// the grid's own proportions, and the cells are as large as the page
/// allows. The cell's side and the figure's size; refused when the cells
/// would be too small to read (under [`CELL_LEAST`]).
fn cells(
    page: Page,
    own: f64,
    (fw, fh): (f64, f64),
    (cols, rows): (usize, usize),
) -> Result<(f64, (f64, f64)), String> {
    let (cols, rows) = (cols.max(1) as f64, rows.max(1) as f64);
    let (ow, oh) = (fw + own * cols, fh + own * rows);
    let (w, h) = match page {
        (None, None) => return Ok((own, (ow, oh))),
        (Some(w), None) => (w, w * oh / ow),
        (None, Some(h)) => (h * ow / oh, h),
        (Some(w), Some(h)) => (w, h),
    };
    let cell = ((w - fw) / cols).min((h - fh) / rows);
    if cell < CELL_LEAST {
        return Err(format!(
            "the grid needs at least {:.1} by {:.1} inches for its labels and cells",
            (fw + CELL_LEAST * cols) / size::PT,
            (fh + CELL_LEAST * rows) / size::PT
        ));
    }
    Ok((cell, (w, h)))
}

/// The width of `s` in `f`-point text, about.
fn width(s: &str, f: f64) -> f64 {
    s.chars().count() as f64 * f * 0.56
}

/// The figure's title at its top left, in `f + 1`-point text, shortened
/// from the front to fit a figure `w` points wide.
fn heading(c: &mut Canvas, title: &str, f: f64, w: f64) {
    let shown = title_text(title, f, w - f);
    c.text(f / 2.0, 1.6 * f, &esc(&shown), f + 1.0, Anchor::Start, INK);
}

/// The title in `f + 1`-point text, at most `room` points wide: shortened
/// from the front when it would not fit.
fn title_text(title: &str, f: f64, room: f64) -> String {
    let fits = (room / width("x", f + 1.0)) as usize;
    let n = title.chars().count();
    match n > fits && fits > 1 {
        true => std::iter::once('…')
            .chain(title.chars().skip(n + 1 - fits))
            .collect(),
        false => title.to_string(),
    }
}

/// The figure `plot` of `m` under `look`, titled `title`, laid out for
/// the page `size` (see [`crate::figure::SizeArgs`]); why there is none
/// (no trace, one batch).
pub fn draw(
    m: &Modules,
    plot: Plot,
    look: Look,
    title: &str,
    size: Page,
) -> Result<String, String> {
    match plot {
        Plot::Modules => heatmap(m, look, title, size),
        Plot::Trace => trace(m, title, size),
        Plot::Batches => batches(m, title, size),
    }
}

/// What colours a block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Look {
    /// The block's fitted shift (a count fit's log rate), on a diverging
    /// scale.
    Shift,
    /// Its methylated % (a count fit's UMIs per cell), light to dark.
    Percent,
    /// That standardized across the cell modules, within each gene module.
    Z,
}

/// The diverging pair and its neutral midpoint, and the sequential ramp's
/// ends.
const LOW: (f64, f64, f64) = (42.0, 120.0, 214.0);
const MID: (f64, f64, f64) = (240.0, 239.0, 236.0);
const HIGH: (f64, f64, f64) = (227.0, 73.0, 72.0);
const DARK: (f64, f64, f64) = (26.0, 72.0, 140.0);

/// `a` to `b` by `t` (0 to 1), as `#rrggbb`.
fn mix(a: (f64, f64, f64), b: (f64, f64, f64), t: f64) -> String {
    let t = t.clamp(0.0, 1.0);
    let c = |x: f64, y: f64| (x + (y - x) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", c(a.0, b.0), c(a.1, b.1), c(a.2, b.2))
}

/// A value's fill: diverging around 0 to `top` either way, or sequential
/// from 0 to `top`; and whether it is dark enough for light ink.
fn fill(v: f64, top: f64, diverging: bool) -> (String, bool) {
    let t = (v.abs() / top.max(1e-12)).min(1.0);
    let far = match (diverging, v >= 0.0) {
        (false, _) => DARK,
        (true, true) => HIGH,
        (true, false) => LOW,
    };
    (mix(MID, far, t), t > 0.6)
}

/// Each gene module's value in each cell module under `look`, rows the
/// gene modules with sites in the diagonal order, columns as the panes show
/// them; none for a block without reads.
fn values(m: &Modules, look: Look) -> (Vec<usize>, Vec<Vec<Option<f64>>>) {
    let rows = m.diagonal.clone();
    let values = rows
        .iter()
        .map(|&g| match look {
            Look::Shift => (m.columns.iter()).map(|&c| Some(m.shift(c, g))).collect(),
            Look::Percent => (m.block(g).iter())
                .map(|r| {
                    r.fraction()
                        .map(|f| if m.per_cell() { f } else { 100.0 * f })
                })
                .collect(),
            Look::Z => standardized(m.block(g)),
        })
        .collect();
    (rows, values)
}

/// The figure of `m` under `look`, titled `title`: a box per block with its
/// value, the gene modules' names and sizes on the left, the cell modules'
/// on top, and the colour scale on the right.
fn heatmap(m: &Modules, look: Look, title: &str, size: Page) -> Result<String, String> {
    let (rows, values) = values(m, look);
    let all: Vec<f64> = values.iter().flatten().flatten().copied().collect();
    let diverging = look != Look::Percent;
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
    };
    let f = font(size);
    let names: Vec<String> = (rows.iter())
        .map(|&g| format!("M{} ({} genes)", g + 1, m.gene_sizes[g]))
        .collect();
    let label = names.iter().map(|n| width(n, f)).fold(0.0, f64::max) + f;
    let (head, scale, foot) = (4.0 * f, 6.0 * f, 3.0 * f);
    let (cols, n) = (m.columns.len(), rows.len());
    // Square boxes as large as the page allows (44 points without one).
    let (cell, (w, h)) = cells(size, 44.0, (label + scale, head + foot), (cols, n))?;
    let font = (cell * 0.28).min(f);
    Ok(svg(w, h, |c: &mut Canvas| {
        heading(c, title, f, w);
        for (j, &role) in m.columns.iter().enumerate() {
            let x = label + cell * (j as f64 + 0.5);
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
                let x = label + cell * j as f64;
                let Some(v) = v else {
                    c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, "#ffffff");
                    continue;
                };
                let (colour, dark) = fill(*v, top, diverging);
                c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, &colour);
                let text = match look {
                    Look::Shift | Look::Z => format!("{v:+.2}"),
                    Look::Percent => format!("{v:.0}"),
                };
                let ink = if dark { "#ffffff" } else { INK };
                c.text(
                    x + cell / 2.0,
                    y + cell / 2.0 + font / 2.5,
                    &text,
                    font,
                    Anchor::Middle,
                    ink,
                );
            }
        }
        // The scale: a column of steps from the top value down to its other
        // end, labelled at both ends and the middle.
        let (sx, sy, sh, steps) = (label + cell * cols as f64 + f, head, cell * n as f64, 40);
        let low = if diverging { -top } else { 0.0 };
        for k in 0..steps {
            let v = top - (top - low) * (k as f64 + 0.5) / steps as f64;
            let y = sy + sh * k as f64 / steps as f64;
            c.rect(
                sx,
                y,
                f * 1.2,
                sh / steps as f64 + 0.5,
                &fill(v, top, diverging).0,
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

/// A plot's frame: its y range ticked at the ends and the middle on the
/// left of `(x, y, w, h)`, its x range under it, and their labels.
#[allow(clippy::too_many_arguments)]
fn axes(
    c: &mut Canvas,
    (x, y, w, h): (f64, f64, f64, f64),
    (x0, x1): (f64, f64),
    (y0, y1): (f64, f64),
    (xl, yl): (&str, &str),
    f: f64,
) {
    c.line(x, y + h, x + w, y + h, MUTED, 0.8);
    c.line(x, y, x, y + h, MUTED, 0.8);
    for (v, at) in [(y1, y), ((y0 + y1) / 2.0, y + h / 2.0), (y0, y + h)] {
        c.line(x - f * 0.3, at, x, at, MUTED, 0.8);
        c.text(
            x - f * 0.5,
            at + f * 0.35,
            &tick(v),
            f * 0.9,
            Anchor::End,
            MUTED,
        );
    }
    for k in 0..=4 {
        let (v, at) = (x0 + (x1 - x0) * k as f64 / 4.0, x + w * k as f64 / 4.0);
        c.line(at, y + h, at, y + h + f * 0.3, MUTED, 0.8);
        c.text(
            at,
            y + h + f * 1.3,
            &format!("{v:.0}"),
            f * 0.9,
            Anchor::Middle,
            MUTED,
        );
    }
    c.text(x + w / 2.0, y + h + f * 2.6, xl, f, Anchor::Middle, INK);
    c.rotated(x - TICK * f - f, y + h / 2.0, yl, f, Anchor::Middle, INK);
}

/// An axis value in a few characters: thousands and millions shortened.
fn tick(v: f64) -> String {
    match v.abs() {
        a if a >= 1e6 => format!("{:.3}M", v / 1e6),
        a if a >= 1e4 => format!("{:.1}k", v / 1e3),
        _ => format!("{v:.2}"),
    }
}

/// A tick label's width at most, in font sizes.
const TICK: f64 = 4.5;

/// The trace: each restart's or chain's line over its steps, the start
/// (k-means) or burn-in's end marked, in a colour of its own.
fn trace(m: &Modules, title: &str, size: Page) -> Result<String, String> {
    let t = m.fit.trace.as_ref().map_err(Clone::clone)?;
    let points: Vec<(f64, f64)> = t
        .chains
        .iter()
        .flatten()
        .map(|&(x, v)| (x as f64, v))
        .collect();
    if points.is_empty() {
        return Err("the trace holds no points".into());
    }
    let vi = t.quantity == crate::cluster::progress::ELBO;
    let (x0, x1) = points
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
    let (y0, y1) = points
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
    let f = font(size);
    let (w, h) = size::fill(size, || (480.0, 300.0));
    let (left, top, right, bottom) = (
        TICK * f + 2.5 * f,
        3.0 * f,
        width("restart 0", f) + f,
        3.5 * f,
    );
    let (pw, ph) = (w - left - right, h - top - bottom);
    let sx = |x: f64| left + pw * (x - x0) / (x1 - x0).max(1e-12);
    let sy = |v: f64| top + ph * (1.0 - (v - y0) / (y1 - y0).max(1e-12));
    let colours = [
        "#2a78d6", "#e34948", "#1f9e89", "#d89b23", "#7c5cc4", "#c2549d", "#5c7d8a", "#8a5c3b",
    ];
    Ok(svg(w, h, |c: &mut Canvas| {
        heading(c, title, f, w);
        let names = match vi {
            true => ("k-means iteration, then epoch", "ELBO"),
            false => ("sweep", "log joint"),
        };
        axes(c, (left, top, pw, ph), (x0, x1), (y0, y1), names, f);
        if t.burn_in as f64 > x0 && (t.burn_in as f64) < x1 {
            let x = sx(t.burn_in as f64);
            c.dashed(x, top, x, top + ph, FAINT, 0.8);
        }
        for (i, chain) in t.chains.iter().enumerate() {
            let line: Vec<(f64, f64)> = chain.iter().map(|&(x, v)| (sx(x as f64), sy(v))).collect();
            let colour = colours[i % colours.len()];
            c.polyline(&line, colour, 1.5);
            let name = format!("{} {i}", if vi { "restart" } else { "chain" });
            c.text(
                left + pw + f,
                top + 1.3 * f * (i as f64 + 1.0),
                &name,
                f,
                Anchor::Start,
                colour,
            );
        }
    }))
}

/// The cell modules against the batches, by condition: each module's share
/// in each batch (its cells by batch size), its evenness and label.
fn batches(m: &Modules, title: &str, size: Page) -> Result<String, String> {
    let b = m.batches.as_ref().ok_or("the cells are of one batch")?;
    let d = &b.design;
    let order = d.order();
    let f = font(size);
    let said: Vec<String> = (m.columns.iter())
        .map(|&r| {
            let s = &b.spread[r];
            let even = s.evenness.map_or("-".into(), |e| format!("{e:.2}"));
            format!("{even}  {}", s.label.word())
        })
        .collect();
    let longest =
        |names: &mut dyn Iterator<Item = &String>| names.map(|s| width(s, f)).fold(0.0, f64::max);
    // Under the title, the batches' names upright, read from below.
    let (label, head) = (3.5 * f, 2.6 * f + longest(&mut d.batches.iter()) * 1.1 + f);
    let (right, foot) = (longest(&mut said.iter()) + 2.0 * f, 3.0 * f);
    let (cols, n) = (order.len(), m.columns.len());
    let (cell, (w, h)) = cells(size, 40.0, (label + right, head + foot), (cols, n))?;
    let font = (cell * 0.28).min(f);
    Ok(svg(w, h, |c: &mut Canvas| {
        heading(c, title, f, w);
        for (j, &bt) in order.iter().enumerate() {
            let x = label + cell * (j as f64 + 0.5);
            c.rotated(
                x + f * 0.35,
                head - f * 0.6,
                &esc(&d.batches[bt]),
                f,
                Anchor::Start,
                INK,
            );
        }
        let right_x = label + cell * cols as f64 + f;
        c.text(
            right_x,
            head - f * 0.6,
            "even  label",
            f,
            Anchor::Start,
            MUTED,
        );
        for (i, (&role, said)) in m.columns.iter().zip(&said).enumerate() {
            let s = &b.spread[role];
            let y = head + cell * i as f64;
            c.text(
                label - f * 0.6,
                y + cell / 2.0 + f / 2.5,
                &format!("C{}", role + 1),
                f,
                Anchor::End,
                INK,
            );
            let top = s.per_batch.iter().copied().fold(0.0, f64::max);
            for (j, &bt) in order.iter().enumerate() {
                let p = s.per_batch[bt];
                let (colour, dark) = fill(p, top, false);
                let x = label + cell * j as f64;
                c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, &colour);
                let ink = if dark { "#ffffff" } else { INK };
                c.text(
                    x + cell / 2.0,
                    y + cell / 2.0 + font / 2.5,
                    &format!("{:.0}%", 100.0 * p),
                    font,
                    Anchor::Middle,
                    ink,
                );
            }
            let ink = if s.label.is_batch() {
                crate::figure::ACCENT
            } else {
                INK
            };
            c.text(
                right_x,
                y + cell / 2.0 + f / 2.5,
                &esc(said),
                f,
                Anchor::Start,
                ink,
            );
        }
        let listing = d.listing().join(" · ");
        c.text(label, h - f, &esc(&listing), f * 0.9, Anchor::Start, MUTED);
    }))
}
