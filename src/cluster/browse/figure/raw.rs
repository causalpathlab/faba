//! The raw view: the modules grid's layout, each block the heatmap of its
//! cell module's cells (rows when transposed, the strongest member first)
//! against its gene module's genes or sites, rasterised (a PNG a block,
//! each pixel pooling the cells and features it covers), under the
//! modules' names, their lines and the colour key, drawn as vectors.

use super::grid::{module_name, quantile, Placed};
use super::{key_strip, scale_column, width, GridArgs, Ink, Look, Refused};
use crate::cluster::browse::raw::{Raw, RawArgs, RawBands, RawValue};
use crate::cluster::browse::Modules;
use crate::figure::size::Page;
use crate::figure::{hex, mixed, svg, Anchor, Canvas, INK, MUTED};

/// Pixels a point a figure's tiles are drawn with at most, about 290
/// dpi; a view on screen asks for its own.
pub const PRINT: f64 = 4.0;
/// A block's side at the figure's own size, as the grid's.
const OWN: f64 = 44.0;
/// The least the plot is drawn in, a side, in points.
const LEAST: f64 = 24.0;

/// The ramp's ends, light to dark, in colour and in grays: lighter than the
/// grid's light end is, so a pixel with too few reads, white, stands apart
/// from one at 0.
const LOW: (f64, f64, f64) = (214.0, 228.0, 242.0);
const HIGH: (f64, f64, f64) = (26.0, 72.0, 140.0);
const GRAY_LOW: (f64, f64, f64) = (222.0, 222.0, 222.0);
const GRAY_HIGH: (f64, f64, f64) = (30.0, 30.0, 30.0);
const BLANK: (u8, u8, u8) = (255, 255, 255);
/// The modules' lines over the tiles.
const EDGE: (&str, f64) = ("#ffffff", 0.6);

/// What the raw figure is drawn from, and how: the reads, the options
/// asked for, the modules' bands (the screen's gene tracks need equal
/// ones), the tiles' pixels a point, and whether only its layout is wanted
/// (no tiles drawn).
#[derive(Clone, Copy)]
pub struct RawDraw<'a> {
    pub raw: &'a Raw,
    pub shown: &'a RawArgs,
    pub bands: RawBands,
    pub per_point: f64,
    pub measure: bool,
}

impl<'a> RawDraw<'a> {
    /// `raw` drawn as `shown` asks, for print.
    pub fn of(raw: &'a Raw, shown: &'a RawArgs) -> Self {
        Self {
            raw,
            shown,
            bands: shown.bands,
            per_point: PRINT,
            measure: false,
        }
    }
}

/// `v`'s colour on a ramp to `top`.
fn colour(v: f64, top: f64, gray: bool) -> (u8, u8, u8) {
    let t = v / top.max(1e-12);
    match gray {
        true => mixed(GRAY_LOW, GRAY_HIGH, t),
        false => mixed(LOW, HIGH, t),
    }
}

/// What the colour is, said.
pub(super) fn said(value: RawValue, unit: crate::cluster::browse::raw::RawUnit) -> String {
    let unit = unit.noun();
    match value {
        RawValue::Fraction => format!("methylated fraction, each cell at each {unit}"),
        RawValue::Methylated => format!("converted reads per cell and {unit}"),
        RawValue::Log1p => format!("log(1 + reads) per cell and {unit}"),
    }
}

/// Bands over `length` points, one for each of `sizes`: as wide as each
/// other, or by size; each one's start and width.
fn bands(sizes: &[usize], length: f64, by: RawBands) -> Vec<(f64, f64)> {
    let total = sizes.iter().sum::<usize>().max(1) as f64;
    let mut at = 0.0;
    (sizes.iter())
        .map(|&n| {
            let w = match by {
                RawBands::Equal => length / sizes.len().max(1) as f64,
                RawBands::Size => length * n as f64 / total,
            };
            at += w;
            (at - w, w)
        })
        .collect()
}

/// `pixels` (cells by features, row by row, `rows` by `cols`) as a PNG,
/// features down and cells across unless `cells_down`; written fast, as
/// its blocks of colour compress well anyway.
fn png(pixels: &[(u8, u8, u8)], (rows, cols): (usize, usize), cells_down: bool) -> Vec<u8> {
    use image::codecs::png::{CompressionType, FilterType, PngEncoder};
    let (w, h) = match cells_down {
        true => (cols, rows),
        false => (rows, cols),
    };
    let img = image::RgbImage::from_fn(w as u32, h as u32, |x, y| {
        let (r, q) = match cells_down {
            true => (y as usize, x as usize),
            false => (x as usize, y as usize),
        };
        let (a, b, c) = pixels[r * cols + q];
        image::Rgb([a, b, c])
    });
    let mut out = Vec::new();
    let encoder =
        PngEncoder::new_with_quality(&mut out, CompressionType::Fast, FilterType::NoFilter);
    img.write_with_encoder(encoder)
        .expect("a PNG is written to memory");
    out
}

/// The raw view of `m` as `how` asks, its rows the gene modules in the
/// grid's order under `look` (the cell modules, transposed by `grid`),
/// titled `title`, laid out for `size`; and where it put its rows, for the
/// screen's gene tracks.
pub(super) fn placed(
    m: &Modules,
    how: RawDraw,
    grid: &GridArgs,
    look: Look,
    title: &str,
    (size, ink): (Page, Ink),
) -> Result<(String, Placed), Refused> {
    let RawDraw { raw, shown, .. } = how;
    let f = ink.f;
    let tr = grid.transpose;
    let genes: Vec<usize> = super::values(m, look).0;
    let cells: &[usize] = &m.columns;
    let gene_names: Vec<String> = genes.iter().map(|&g| module_name(m, g, grid)).collect();
    let cell_names: Vec<String> = (cells.iter())
        .map(|&c| format!("C{} ({})", c + 1, raw.cells(c)))
        .collect();
    let gene_sizes: Vec<usize> = genes.iter().map(|&g| raw.features(g)).collect();
    let cell_sizes: Vec<usize> = cells.iter().map(|&c| raw.cells(c)).collect();
    // The rows' modules and the columns'.
    let ((row_names, row_sizes), (col_names, col_sizes)) = match tr {
        false => ((&gene_names, &gene_sizes), (&cell_names, &cell_sizes)),
        true => ((&cell_names, &cell_sizes), (&gene_names, &gene_sizes)),
    };
    let label = row_names.iter().map(|n| width(n, f)).fold(0.0, f64::max) + f;
    let names_w = col_names.iter().map(|n| width(n, f)).fold(0.0, f64::max);
    // The key: a column right of the plot, or under the plot alone with
    // `--key`, a line.
    let small_key = ink.panel_only && grid.key;
    let scale = ink.room(6.0 * f);
    let foot = ink.room(2.0 * f) + if small_key { 2.6 * f } else { 0.0 } + 0.5 * f;
    let (across, down) = (col_names.len().max(1) as f64, row_names.len().max(1) as f64);
    let page = |names_h: f64| {
        let fixed = (label + scale + 0.5 * f, ink.room(2.0 * f) + names_h + foot);
        let own = (fixed.0 + OWN * across, fixed.1 + OWN * down);
        let (w, h) = match size {
            (None, None) => own,
            (Some(w), None) => (w, w * own.1 / own.0),
            (None, Some(h)) => (h * own.0 / own.1, h),
            (Some(w), Some(h)) => (w, h),
        };
        (fixed, (w, h))
    };
    // The columns' names across when they fit their bands, else turned.
    let flat = |names_h: f64| {
        let (fixed, (w, _)) = page(names_h);
        let bands = bands(col_sizes, w - fixed.0, how.bands);
        let least = bands.iter().map(|b| b.1).fold(f64::INFINITY, f64::min);
        names_w + 0.4 * f <= least
    };
    let turned = !flat(1.6 * f);
    let names_h = if turned { names_w + 0.8 * f } else { 1.6 * f };
    let (fixed, (w, h)) = page(names_h);
    let (pw, ph) = (w - fixed.0, h - fixed.1);
    if pw < LEAST || ph < LEAST {
        let need = |side: Option<f64>, fixed: f64| side.map(|s| s.max(fixed + LEAST) + 1.0);
        return Err(Refused::Small((
            need(size.0, fixed.0),
            need(size.1, fixed.1),
        )));
    }
    let (x0, y0) = (label, ink.room(2.0 * f) + names_h);
    let placed = Placed {
        go_x: None,
        rows_top: y0,
        pitch: ph / down,
        width: (x0 + pw + scale).min(w),
    };
    if how.measure {
        return Ok((String::new(), placed));
    }
    let rows = bands(row_sizes, ph, how.bands);
    let cols = bands(col_sizes, pw, how.bands);
    // Each block's tile, cells by features, its pixels' values.
    let mut tiles = Vec::new();
    let mut all = Vec::new();
    for (i, &(y, bh)) in rows.iter().enumerate() {
        for (j, &(x, bw)) in cols.iter().enumerate() {
            let (c, g, cell_pts, gene_pts) = match tr {
                false => (cells[j], genes[i], bw, bh),
                true => (cells[i], genes[j], bh, bw),
            };
            let px = |pts: f64| (pts * how.per_point).ceil().max(1.0) as usize;
            let tile = raw.tile(c, g, (px(cell_pts), px(gene_pts)));
            let values: Vec<Option<f64>> = (tile.pixels.iter())
                .map(|p| p.value(shown.value, shown.min_reads))
                .collect();
            all.extend(values.iter().flatten());
            tiles.push(((x0 + x, y0 + y, bw, bh), (tile.rows, tile.cols), values));
        }
    }
    let top = all.iter().copied().fold(0.0, f64::max);
    // The colours' top: the 99th percentile, or `--grid-cap`'s.
    let shade_top = quantile(&all, grid.cap.unwrap_or(0.99), false).min(top);
    let capped = shade_top < top;
    let tick = |v: f64| match shown.value {
        RawValue::Fraction => format!("{:.0}%", 100.0 * v),
        _ => format!("{v:.2}"),
    };
    let end = |v: f64| match capped && v >= shade_top {
        true => format!("≥{}", tick(v)),
        false => tick(v),
    };
    let what = said(shown.value, raw.unit);
    let blank = format!("white: under {} reads", shown.min_reads);
    let ends = [tick(0.0), tick(shade_top / 2.0), end(shade_top)];
    let drawn = svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        for (rect, dims, values) in &tiles {
            let pixels: Vec<(u8, u8, u8)> = (values.iter())
                .map(|v| v.map_or(BLANK, |v| colour(v, shade_top, ink.gray)))
                .collect();
            c.image(*rect, &png(&pixels, *dims, tr));
        }
        // The modules' lines between the bands.
        for &(y, _) in rows.iter().skip(1) {
            c.line(x0, y0 + y, x0 + pw, y0 + y, EDGE.0, EDGE.1);
        }
        for &(x, _) in cols.iter().skip(1) {
            c.line(x0 + x, y0, x0 + x, y0 + ph, EDGE.0, EDGE.1);
        }
        for (name, &(y, bh)) in row_names.iter().zip(&rows) {
            let y = y0 + y + bh / 2.0 + f / 2.5;
            c.text(x0 - 0.6 * f, y, name, f, Anchor::End, INK);
        }
        for (name, &(x, bw)) in col_names.iter().zip(&cols) {
            let x = x0 + x + bw / 2.0;
            match turned {
                true => c.rotated(x + f / 2.5, y0 - 0.4 * f, name, f, Anchor::Start, INK),
                false => c.text(x, y0 - 0.6 * f, name, f, Anchor::Middle, INK),
            }
        }
        if small_key {
            // Under the plot: the ramp, its ends and middle, and what
            // white is.
            let (ky, kw) = (y0 + ph + 0.6 * f, (pw / 2.0).max(6.0 * f).min(pw));
            let fill = |along: f64| hex(colour(shade_top * along, shade_top, ink.gray));
            let said = [ends[0].as_str(), ends[1].as_str(), ends[2].as_str()];
            key_strip(c, (x0, ky, kw), f, 40, fill, said);
            c.text(x0 + kw + f, ky + 0.7 * f, &blank, f, Anchor::Start, MUTED);
        }
        if ink.panel_only {
            return;
        }
        // The scale right of the plot, top value first.
        let fill = |along: f64| hex(colour(shade_top * (1.0 - along), shade_top, ink.gray));
        let said = [ends[2].as_str(), ends[1].as_str(), ends[0].as_str()];
        scale_column(c, (x0 + pw + f, y0, ph), f, fill, said);
        let n: usize = cell_sizes.iter().sum();
        let k: usize = gene_sizes.iter().sum();
        let unit = raw.unit.noun();
        let note = format!("{what}; {n} cells × {k} {unit}s; {blank}");
        let note = crate::figure::fit_text(&note, w - 1.0 * f, f);
        c.text(0.5 * f, h - f, &note, f, Anchor::Start, MUTED);
    });
    Ok((drawn, placed))
}
