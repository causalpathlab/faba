//! The raw view: the modules grid's layout, each block the heatmap of its
//! cell module's cells (rows when transposed, the strongest member first)
//! against its gene module's genes or sites, rasterised (a PNG a block,
//! each pixel pooling the cells and features it covers), under the
//! modules' names, their lines and the colour key, drawn as vectors.

use super::grid::{module_name, quantile, Placed};
use super::{key_strip, scale_column, width, GridArgs, Ink, Look, Refused};
use crate::cluster::browse::raw::{Raw, RawArgs, RawBands, RawSort, RawValue};
use crate::cluster::browse::Modules;
use crate::figure::size::Page;
use crate::figure::{hex, mixed, svg, Anchor, Canvas, INK, MUTED};
use crate::qc::arms::Rgb;

/// A block's side at the figure's own size, as the grid's.
const OWN: f64 = 44.0;
/// Pixels a point a shaded (`--raw-ramp`) tile is drawn at, for print.
const PRINT: f64 = 4.0;
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
/// ones), the most dots a point (a screen's pixels; none for print, where
/// `--raw-dot` alone says), and whether only its layout is wanted (no
/// tiles drawn).
#[derive(Clone, Copy)]
pub struct RawDraw<'a> {
    pub raw: &'a Raw,
    pub shown: &'a RawArgs,
    pub bands: RawBands,
    pub per_point: Option<f64>,
    pub measure: bool,
}

impl<'a> RawDraw<'a> {
    /// `raw` drawn as `shown` asks, for print.
    pub fn of(raw: &'a Raw, shown: &'a RawArgs) -> Self {
        Self {
            raw,
            shown,
            bands: shown.bands,
            per_point: None,
            measure: false,
        }
    }
}

/// The ramp's ends, light to dark: to `--raw-colour`'s (a dot's) from a
/// tint of it, else in colour or in grays.
fn ramp_ends(asked: Option<Rgb>, gray: bool) -> Ends {
    match (asked, gray) {
        (Some((r, g, b)), _) => {
            let high = (f64::from(r), f64::from(g), f64::from(b));
            let low = mixed((255.0, 255.0, 255.0), high, 0.15);
            ((f64::from(low.0), f64::from(low.1), f64::from(low.2)), high)
        }
        (None, true) => (GRAY_LOW, GRAY_HIGH),
        (None, false) => (LOW, HIGH),
    }
}

/// A ramp's ends, light to dark.
type Ends = ((f64, f64, f64), (f64, f64, f64));

/// `v`'s colour on a ramp to `top`.
fn colour(v: f64, top: f64, (low, high): Ends) -> (u8, u8, u8) {
    mixed(low, high, v / top.max(1e-12))
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
/// other, or by size, each within `cap` (the least and the most) of an
/// equal band; each one's start and width.
fn bands(sizes: &[usize], length: f64, by: RawBands, cap: (f64, f64)) -> Vec<(f64, f64)> {
    let n = sizes.len().max(1) as f64;
    let mut share: Vec<f64> = match by {
        RawBands::Equal => vec![1.0 / n; sizes.len()],
        RawBands::Size => {
            let total = sizes.iter().sum::<usize>().max(1) as f64;
            sizes.iter().map(|&k| k as f64 / total).collect()
        }
    };
    // Each share its size's times `t`, clamped to the cap: `t` found by
    // halving so that they fill the length (they can, the cap taking one
    // equal band in).
    let (lo, hi) = (cap.0 / n, cap.1 / n);
    let scaled = |s: f64, t: f64| (s * t).clamp(lo, hi);
    let filled = |t: f64| share.iter().map(|&s| scaled(s, t)).sum::<f64>();
    let (mut low, mut high) = (0.0, 1.0);
    while filled(high) < 1.0 && high < 1e12 {
        high *= 2.0;
    }
    for _ in 0..100 {
        let mid = (low + high) / 2.0;
        match filled(mid) < 1.0 {
            true => low = mid,
            false => high = mid,
        }
    }
    share.iter_mut().for_each(|s| *s = scaled(*s, high));
    let total: f64 = share.iter().sum::<f64>().max(1e-12);
    let mut at = 0.0;
    (share.iter())
        .map(|s| {
            let w = length * s / total;
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
        let bands = bands(col_sizes, w - fixed.0, how.bands, shown.band_cap());
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
    let rows = bands(row_sizes, ph, how.bands, shown.band_cap());
    let cols = bands(col_sizes, pw, how.bands, shown.band_cap());
    // Each block's tile, cells by features, its pixels' values.
    // Square dots, as many across as down a point: a dot `--raw-dot`
    // points a side (no finer than the screen's pixels); a shade at the
    // screen's pixels, or print's.
    let per_point = match shown.ramp {
        true => how.per_point.unwrap_or(PRINT),
        false => (1.0 / shown.dot).min(how.per_point.unwrap_or(f64::INFINITY)),
    };
    let mut tiles = Vec::new();
    let mut all = Vec::new();
    for (i, &(y, bh)) in rows.iter().enumerate() {
        for (j, &(x, bw)) in cols.iter().enumerate() {
            let (c, g, cell_pts, gene_pts) = match tr {
                false => (cells[j], genes[i], bw, bh),
                true => (cells[i], genes[j], bh, bw),
            };
            let px = |pts: f64| (pts * per_point).round().max(1.0) as usize;
            let tile = raw.tile(c, g, (px(cell_pts), px(gene_pts)), shown);
            // A shade's value; or a dot (1), any pair under it a dot.
            let values: Vec<Option<f64>> = (tile.pixels.iter())
                .map(|p| match shown.ramp {
                    true => p.value(shown.value, shown.min_reads),
                    false => p.dot().then_some(1.0),
                })
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
    // The threshold is a beta whatever the shades are.
    let beta = format!("β≥{:.0}%", 100.0 * shown.threshold);
    let end = |v: f64| match capped && v >= shade_top {
        true => format!("≥{}", tick(v)),
        false => tick(v),
    };
    let what = said(shown.value, raw.unit);
    // What a dot is, and what white is.
    let unit = raw.unit.noun();
    let (what, blank) = match shown.ramp {
        true => (what, format!("white: under {} reads", shown.min_reads)),
        false => (
            format!(
                "dots where a cell-{unit} pair has {beta}, ≥{} of ≥{} reads converted",
                shown.min_methylated, shown.min_reads
            ),
            "white: otherwise".to_string(),
        ),
    };
    let dots = format!("{beta}, ≥{} converted; {blank}", shown.min_methylated);
    let ends = [tick(0.0), tick(shade_top / 2.0), end(shade_top)];
    let drawn = svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        // A dot at the threshold or above, or a shade on the ramp.
        let ramp = ramp_ends(shown.colour, ink.gray);
        let dot = colour(1.0, 1.0, ramp);
        // A dot's value is 1, the top of its ramp.
        let pixel = |v: Option<f64>| v.map_or(BLANK, |v| colour(v, shade_top, ramp));
        for (rect, dims, values) in &tiles {
            let pixels: Vec<(u8, u8, u8)> = values.iter().map(|&v| pixel(v)).collect();
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
        if small_key && !shown.ramp {
            // Under the plot: a dot, and what it and white are.
            let ky = y0 + ph + 0.6 * f;
            c.rect(x0, ky, 0.7 * f, 0.7 * f, &hex(dot));
            c.text(x0 + 1.1 * f, ky + 0.7 * f, &dots, f, Anchor::Start, MUTED);
        } else if small_key {
            // Under the plot: the ramp, its ends and middle, and what
            // white is.
            let (ky, kw) = (y0 + ph + 0.6 * f, (pw / 2.0).max(6.0 * f).min(pw));
            let fill = |along: f64| hex(colour(shade_top * along, shade_top, ramp));
            let said = [ends[0].as_str(), ends[1].as_str(), ends[2].as_str()];
            key_strip(c, (x0, ky, kw), f, 40, fill, said);
            c.text(x0 + kw + f, ky + 0.7 * f, &blank, f, Anchor::Start, MUTED);
        }
        if ink.panel_only {
            return;
        }
        // The scale right of the plot, top value first; or a dot, and
        // what it is.
        match shown.ramp {
            true => {
                let fill = |along: f64| hex(colour(shade_top * (1.0 - along), shade_top, ramp));
                let said = [ends[2].as_str(), ends[1].as_str(), ends[0].as_str()];
                scale_column(c, (x0 + pw + f, y0, ph), f, fill, said);
            }
            false => {
                let (sx, said) = (x0 + pw + f, &beta);
                c.rect(sx, y0, 1.2 * f, 1.2 * f, &hex(dot));
                c.text(sx, y0 + 2.4 * f, said, 0.9 * f, Anchor::Start, MUTED);
            }
        }
        let n: usize = cell_sizes.iter().sum();
        let k: usize = gene_sizes.iter().sum();
        let sorted = match shown.sort {
            RawSort::Value => "; sorted by dots in each block",
            RawSort::Reads => "; sorted by reads in each block",
            RawSort::None => "",
        };
        let white = if shown.ramp {
            format!("; {blank}")
        } else {
            String::new()
        };
        let note = format!("{what}; {n} cells × {k} {unit}s{sorted}{white}");
        let note = crate::figure::fit_text(&note, w - 1.0 * f, f);
        c.text(0.5 * f, h - f, &note, f, Anchor::Start, MUTED);
    });
    Ok((drawn, placed))
}

#[cfg(test)]
#[path = "../tests/raw_bands.rs"]
mod tests;
