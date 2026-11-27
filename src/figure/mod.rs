//! Figures behind the full-screen views: each view draws itself as an SVG
//! here, and [`save_with_png`] writes that one drawing as a vector PDF and a PNG, so
//! the exported figure is the same whatever the terminal could show.

use std::fmt::Write as _;
use std::sync::{Arc, OnceLock};

/// Width of `HistPlot`'s y-axis gutter, in terminal columns.
pub const GUTTER: u16 = 6;

pub mod gallery;
pub mod logo;
pub mod size;
pub mod style;
pub mod term;

pub use size::SizeArgs;

use crate::tui::ui::{compact, help_line, input_line, Scale, DIM, HIGHLIGHT, PLAIN};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use ratatui_image::picker::Picker;
use ratatui_image::protocol::Protocol;
use ratatui_image::Image;

/// Text and axes.
pub const INK: &str = "#1a1a1a";
/// Bars.
pub const BAR: &str = "#4d4d4d";
/// Secondary text and axis labels.
pub const MUTED: &str = "#6b6b6b";
/// Bars behind a front subset, gridlines.
pub const FAINT: &str = "#cfcfcf";
/// The one accent, as on screen: what a threshold drops, the cursor.
pub const ACCENT: &str = "#d97757";

/// `#rrggbb` as its parts, each 0 to 255.
pub fn rgb(hex: &str) -> (f64, f64, f64) {
    let p = |i: usize| {
        hex.get(i..i + 2)
            .and_then(|h| u8::from_str_radix(h, 16).ok())
    };
    let p = |i: usize| p(i).unwrap_or(0) as f64;
    (p(1), p(3), p(5))
}

/// `a` to `b` at `t` (0 `a`, 1 `b`), as `#rrggbb`.
pub fn mix(a: (f64, f64, f64), b: (f64, f64, f64), t: f64) -> String {
    hex(mixed(a, b, t))
}

/// A colour's parts as `#rrggbb`.
pub fn hex((r, g, b): (u8, u8, u8)) -> String {
    format!("#{r:02x}{g:02x}{b:02x}")
}

/// [`mix`], as its parts.
pub fn mixed(a: (f64, f64, f64), b: (f64, f64, f64), t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    let c = |x: f64, y: f64| (x + (y - x) * t).round() as u8;
    (c(a.0, b.0), c(a.1, b.1), c(a.2, b.2))
}
const FONT: &str = "Helvetica, Arial, 'DejaVu Sans', sans-serif";

/// Escape text for SVG content and attributes.
pub(crate) fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub enum Anchor {
    Start,
    Middle,
    End,
}

impl Anchor {
    fn svg(self) -> &'static str {
        match self {
            Anchor::Start => "start",
            Anchor::Middle => "middle",
            Anchor::End => "end",
        }
    }
}

/// An SVG drawing in points, white background, in a [`style::Style`]:
/// drawn at its 8-point layout on a page shrunk or grown by the style's
/// font ([`Canvas::page`]), inside one scale to the size asked, so its
/// text, margins and strokes scale together.
pub struct Canvas {
    width: f64,
    height: f64,
    /// The white page under the drawing, at the size asked.
    ground: String,
    body: String,
    style: style::Style,
}

impl Canvas {
    pub fn new(width: f64, height: f64) -> Self {
        Self::styled(width, height, &style::Style::default())
    }

    /// A drawing `width` by `height` points in `style`: laid out for
    /// 8-point text on [`Canvas::page`], its labels printed scaled to the
    /// font, or every one at the font when asked so.
    pub fn styled(width: f64, height: f64, style: &style::Style) -> Self {
        let mut c = Self {
            width,
            height,
            ground: String::new(),
            body: String::new(),
            style: *style,
        };
        c.rect(0.0, 0.0, width, height, "#ffffff");
        c.ground = std::mem::take(&mut c.body);
        c
    }

    /// The page it is laid out on, for 8-point text: the size asked over
    /// the font's scale ([`style::Style::em`]).
    pub fn page(&self) -> (f64, f64) {
        let em = self.style.em();
        (self.width / em, self.height / em)
    }

    /// The look it is drawn in.
    pub fn style(&self) -> &style::Style {
        &self.style
    }

    /// A label's size as written, for one laid out at `size` points: as
    /// given, or the layout's own 8 points when every label is set at the
    /// font, which the scale prints at that font.
    fn size(&self, size: f64) -> String {
        self.label_size(size).to_string()
    }

    /// The size, in the layout's points, a label laid out at `size` is
    /// written at (see [`Canvas::size`]): what its width is reckoned by.
    pub fn label_size(&self, size: f64) -> f64 {
        match self.style.uniform {
            true => style::OWN_FONT,
            false => size,
        }
    }

    /// A filled rectangle; empty ones (zero or negative size) are skipped,
    /// as SVG renderers reject them.
    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        // Written to two decimals: what would read as zero is skipped too.
        if !(w >= 0.005 && h >= 0.005 && w.is_finite() && h.is_finite()) {
            return;
        }
        let _ = writeln!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{:.2}" height="{:.2}" fill="{fill}"/>"#,
            w.max(0.0),
            h.max(0.0)
        );
    }

    /// A filled rectangle `(x, y, w, h)` outlined `width` points in
    /// `stroke`; empty ones skipped, as [`Self::rect`] skips them.
    pub fn outlined(
        &mut self,
        (x, y, w, h): (f64, f64, f64, f64),
        fill: &str,
        (stroke, width): (&str, f64),
    ) {
        if !(w >= 0.005 && h >= 0.005 && w.is_finite() && h.is_finite()) {
            return;
        }
        let _ = writeln!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}" stroke="{stroke}" stroke-width="{width:.2}"/>"#
        );
    }

    /// A raster image, `png` stretched over `(x, y, w, h)`, each pixel a
    /// sharp-edged block (no smoothing between pixels); empty ones skipped,
    /// as [`Self::rect`] skips them.
    pub fn image(&mut self, (x, y, w, h): (f64, f64, f64, f64), png: &[u8]) {
        use base64::Engine;
        if !(w >= 0.005 && h >= 0.005 && w.is_finite() && h.is_finite()) {
            return;
        }
        let data = base64::engine::general_purpose::STANDARD.encode(png);
        let _ = writeln!(
            self.body,
            r#"<image x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" preserveAspectRatio="none" image-rendering="optimizeSpeed" style="image-rendering:pixelated" href="data:image/png;base64,{data}"/>"#
        );
    }

    pub fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
        let _ = writeln!(
            self.body,
            r#"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="{width}"/>"#
        );
    }

    pub fn dashed(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
        let _ = writeln!(
            self.body,
            r#"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="{width}" stroke-dasharray="3 2"/>"#
        );
    }

    pub fn text(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        let size = self.size(size);
        self.text_at(x, y, s, &size, anchor, fill);
    }

    /// [`Canvas::text`] at `size` as it is, a uniform font or not: for a
    /// label shrunk to fit a small room (see [`Canvas::bold_exact`]).
    pub fn text_exact(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        self.text_at(x, y, s, &exact(size), anchor, fill);
    }

    fn text_at(&mut self, x: f64, y: f64, s: &str, size: &str, anchor: Anchor, fill: &str) {
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-family="{FONT}" font-size="{size}" text-anchor="{}" fill="{fill}">{}</text>"#,
            anchor.svg(),
            esc(s)
        );
    }

    /// [`Canvas::text`] over a white halo, to read where it lies over bars.
    pub fn text_halo(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        let size = self.size(size);
        let _ = writeln!(
            self.body,
            r##"<text x="{x:.2}" y="{y:.2}" font-family="{FONT}" font-size="{size}" text-anchor="{}" fill="{fill}" stroke="#ffffff" stroke-width="2" stroke-linejoin="round" paint-order="stroke">{}</text>"##,
            anchor.svg(),
            esc(s)
        );
    }

    /// Text turned to read bottom to top, centred on `(x, y)`.
    pub fn vtext(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str) {
        let size = self.size(size);
        let _ = writeln!(
            self.body,
            r#"<text transform="translate({x:.2} {y:.2}) rotate(-90)" font-family="{FONT}" font-size="{size}" text-anchor="middle" fill="{fill}">{}</text>"#,
            esc(s)
        );
    }

    /// Text turned 90° counter-clockwise, reading bottom to top, anchored
    /// at `(x, y)` as `anchor` says along its length.
    pub fn rotated(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        let size = self.size(size);
        let _ = writeln!(
            self.body,
            r#"<text transform="translate({x:.2} {y:.2}) rotate(-90)" font-family="{FONT}" font-size="{size}" text-anchor="{}" fill="{fill}">{}</text>"#,
            anchor.svg(),
            esc(s)
        );
    }

    pub fn bold(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        let size = self.size(size);
        self.bold_at(x, y, s, &size, anchor, fill);
    }

    /// [`Canvas::bold`] at `size` as it is, a uniform font or not: for a
    /// title shrunk to fit its room (see [`place_bold_title`]).
    pub fn bold_exact(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
        self.bold_at(x, y, s, &exact(size), anchor, fill);
    }

    fn bold_at(&mut self, x: f64, y: f64, s: &str, size: &str, anchor: Anchor, fill: &str) {
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-family="{FONT}" font-size="{size}" font-weight="bold" text-anchor="{}" fill="{fill}">{}</text>"#,
            anchor.svg(),
            esc(s)
        );
    }

    /// One letter stretched to fill the box `(x, y, w, h)`, as sequence
    /// logos draw them. Metrics are Helvetica Bold's, in em: advance, and the
    /// ink's bottom and top (round letters overshoot baseline and cap).
    pub fn glyph(&mut self, x: f64, y: f64, w: f64, h: f64, ch: char, fill: &str) {
        if w <= 0.0 || h <= 0.0 {
            return;
        }
        const SIZE: f64 = 100.0;
        let (advance, bottom, top) = match ch {
            'C' | 'G' | 'O' | 'U' => (0.74, -0.015, 0.735),
            'T' => (0.61, 0.0, 0.72),
            _ => (0.72, 0.0, 0.72),
        };
        let (sx, sy) = (w / (advance * SIZE), h / ((top - bottom) * SIZE));
        let _ = writeln!(
            self.body,
            r#"<text transform="translate({:.2} {:.2}) scale({sx:.4} {sy:.4})" font-family="{FONT}" font-size="{SIZE}" font-weight="bold" text-anchor="middle" fill="{fill}">{}</text>"#,
            x + w / 2.0,
            y + h + bottom * SIZE * sy,
            esc(&ch.to_string())
        );
    }

    /// An open line through `points`.
    pub fn polyline(&mut self, points: &[(f64, f64)], stroke: &str, width: f64) {
        let pts: Vec<String> = points
            .iter()
            .map(|(x, y)| format!("{x:.2},{y:.2}"))
            .collect();
        let _ = writeln!(
            self.body,
            r#"<polyline points="{}" fill="none" stroke="{stroke}" stroke-width="{width}"/>"#,
            pts.join(" ")
        );
    }

    /// A canvas with no size or background, for shapes to place inside
    /// another drawing with [`Canvas::into_body`].
    pub fn layer() -> Self {
        Self {
            width: 0.0,
            height: 0.0,
            ground: String::new(),
            body: String::new(),
            style: style::Style::default(),
        }
    }

    /// The shapes drawn so far, without the `<svg>` wrapper.
    pub fn into_body(self) -> String {
        self.body
    }

    /// The drawing as SVG: the page, then what was drawn, scaled from its
    /// layout to the size asked (no scale at 8 points).
    pub fn finish(self) -> String {
        let em = self.style.em();
        let body = match em == 1.0 {
            true => self.body,
            false => format!("<g transform=\"scale({em})\">\n{}</g>\n", self.body),
        };
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">
{ground}{body}</svg>
"#,
            w = self.width,
            h = self.height,
            ground = self.ground,
        )
    }
}

/// A size as written exactly: to a hundredth, no trailing zeros.
fn exact(size: f64) -> String {
    let size = format!("{size:.2}");
    size.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Points per terminal cell in a figure saved from a view: a cell is about
/// twice as tall as wide, and the figures' 8- to 12-point text about a row
/// high, as the terminal's is.
pub const CELL: (f64, f64) = (6.0, 12.0);

/// The size in points of a figure saved from a view that took `area` on
/// screen, cells `cell` points each, so it keeps the view's shape; `None`
/// before it was drawn. Both sides grow alike until neither is smaller
/// than `min` (a view states its own), lest a narrow terminal squeeze the
/// text out.
pub fn page(area: Rect, min: (f64, f64), (cw, ch): (f64, f64)) -> Option<(f64, f64)> {
    (area.width > 0 && area.height > 0).then(|| {
        let (w, h) = (f64::from(area.width) * cw, f64::from(area.height) * ch);
        let grow = (min.0 / w).max(min.1 / h).max(1.0);
        // `w * (min / w)` can round a hair under `min`: never below it.
        ((w * grow).max(min.0), (h * grow).max(min.1))
    })
}

/// A figure's panels top to bottom: room that never stretches (title,
/// gaps, fixed rows), and panels each at least its `least` height, which
/// share what a taller page adds in proportion to those.
pub struct Panels<'a> {
    pub fixed: f64,
    pub least: &'a [f64],
}

impl Panels<'_> {
    /// The least page height: the fixed room and every panel at its least,
    /// for [`page`].
    pub fn height(&self) -> f64 {
        self.fixed + self.least.iter().sum::<f64>()
    }

    /// Each panel's height on a page `height` tall.
    pub fn heights(&self, height: f64) -> Vec<f64> {
        let grow = (height - self.fixed) / self.least.iter().sum::<f64>();
        self.least.iter().map(|m| m * grow).collect()
    }
}

/// Points per cell for a terminal font `font` pixels a cell: as high as
/// [`CELL`], as wide for that as the font; [`CELL`] for a font of no size.
pub fn font_cell(font: ratatui_image::FontSize) -> (f64, f64) {
    if font.width == 0 || font.height == 0 {
        return CELL;
    }
    (
        CELL.1 * f64::from(font.width) / f64::from(font.height),
        CELL.1,
    )
}

/// The `(width, height)` an SVG from [`Canvas::finish`] declares; zero
/// where one is missing.
pub(crate) fn svg_size(svg: &str) -> (f64, f64) {
    let get = |key: &str| -> f64 {
        let at = svg.find(&format!(" {key}=\"")).map(|i| i + key.len() + 3);
        at.and_then(|i| svg[i..].split('"').next()?.parse().ok())
            .unwrap_or(0.0)
    };
    (get("width"), get("height"))
}

/// Figures from [`Canvas::finish`] tiled into one, `cols` to a row, in
/// order: one column stacks them top to bottom, as many columns as figures
/// lays them side by side. Each cell is as wide as the widest figure and
/// each row as tall as its tallest.
pub fn tile_svgs(svgs: &[String], cols: usize) -> String {
    let cols = cols.max(1);
    let cell = svgs.iter().map(|s| svg_size(s).0).fold(0.0, f64::max);
    let mut y = 0.0;
    let mut body = String::new();
    for row in svgs.chunks(cols) {
        for (k, svg) in row.iter().enumerate() {
            // Each nested at its offset, as it was drawn.
            let x = cell * k as f64;
            body.push_str(&svg.replacen("<svg ", &format!("<svg x=\"{x}\" y=\"{y}\" "), 1));
        }
        y += row.iter().map(|s| svg_size(s).1).fold(0.0, f64::max);
    }
    let width = cell * cols.min(svgs.len()).max(1) as f64;
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{y}" viewBox="0 0 {width} {y}">
{body}</svg>
"#
    )
}

/// Figure `b` from [`Canvas::finish`] to the right of figure `a`, their
/// tops aligned, each its own width: a key beside a row of panels.
pub fn beside_svgs(a: &str, b: &str) -> String {
    let ((aw, ah), (bw, bh)) = (svg_size(a), svg_size(b));
    let (width, height) = (aw + bw, ah.max(bh));
    let b = b.replacen("<svg ", &format!("<svg x=\"{aw}\" y=\"0\" "), 1);
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
{a}{b}</svg>
"#
    )
}

/// A `w` x `h` drawing made by `draw`, as SVG.
pub fn svg(w: f64, h: f64, draw: impl FnOnce(&mut Canvas)) -> String {
    let mut c = Canvas::new(w, h);
    draw(&mut c);
    c.finish()
}

/// `v` on a view's y scale, as the terminal histogram draws it.
pub(crate) fn scaled(scale: Scale, v: f64) -> f64 {
    match scale {
        Scale::Log => (v.max(0.0) + 1.0).log10(),
        Scale::Sqrt => v.max(0.0).sqrt(),
        Scale::Linear => v.max(0.0),
    }
}

/// The round values an axis on `scale` up to `max` (scaled) may label,
/// each with its share of the axis, from the least: powers of ten on a
/// log axis (10, 100, 1k, ...), evenly spread counts on a square-root one,
/// `spread` of the axis apart at least (see [`even_counts`]); none on a
/// linear one.
pub(crate) fn round_y_ticks(scale: Scale, max: f64, spread: f64) -> Vec<(f64, String)> {
    let top = unscaled(scale, max);
    let values: Vec<f64> = match scale {
        Scale::Linear => Vec::new(),
        Scale::Log => (1..16).map(|k| 10f64.powi(k)).collect(),
        Scale::Sqrt => even_counts(top, spread),
    };
    (values.into_iter())
        .filter(|&v| v <= top * (1.0 + 1e-9) && max > 0.0)
        .map(|v| (scaled(scale, v) / max, round_label(v)))
        .collect()
}

/// The least share of a square-root axis between its top two labels.
pub(crate) const SQRT_SPREAD: f64 = 0.17;

/// The counts a square-root axis up to `top` labels: the multiples of a
/// round step (1, 1.5, 2, 2.5 or 5 times a power of ten, whole), the least
/// whose four or fewer up to the top stand `spread` of the axis apart at
/// the top, where they crowd most.
fn even_counts(top: f64, spread: f64) -> Vec<f64> {
    let steps = (0..16).flat_map(|k| [1.0, 1.5, 2.0, 2.5, 5.0].map(|m| m * 10f64.powi(k)));
    for step in steps.filter(|s: &f64| s.fract() == 0.0) {
        let n = (top / step * (1.0 + 1e-9)).floor() as usize;
        if n == 0 {
            break;
        }
        let share = |k: usize| (k as f64 * step / top).sqrt();
        if n <= 4 && share(1) >= spread && (n == 1 || share(n) - share(n - 1) >= spread) {
            return (1..=n).map(|k| k as f64 * step).collect();
        }
    }
    Vec::new()
}

/// The round values a y axis on `scale` up to `max` (scaled) labels on a
/// stretch `depth` deep (a plot, or a mirror's arm), each with its share of
/// it: a label's height and a little apart, from the top down and none at
/// zero's (a square-root axis's chosen so, to stay evenly spread).
fn kept_round(c: &Canvas, scale: Scale, max: f64, depth: f64) -> Vec<(f64, String)> {
    let apart = 1.3 * c.label_size(8.0) / depth;
    let spread = match scale {
        Scale::Sqrt => apart.max(SQRT_SPREAD),
        _ => SQRT_SPREAD,
    };
    let mut kept: Vec<(f64, String)> = Vec::new();
    for (share, label) in round_y_ticks(scale, max, spread).into_iter().rev() {
        let last = kept.last().map_or(f64::INFINITY, |k| k.0);
        if share <= last - apart && share >= apart {
            kept.push((share, label));
        }
    }
    kept
}

/// A round count as an axis says it: `20`, `500`, `1k`, `2.5k`, `1M`.
fn round_label(v: f64) -> String {
    match v {
        v if v >= 1e6 => format!("{}M", v / 1e6),
        v if v >= 1e3 => format!("{}k", v / 1e3),
        v => format!("{v}"),
    }
}

pub(crate) fn unscaled(scale: Scale, t: f64) -> f64 {
    match scale {
        Scale::Log => 10f64.powf(t) - 1.0,
        Scale::Sqrt => t * t,
        Scale::Linear => t,
    }
}

/// About how wide `s` sets at `size` points.
pub fn text_width(s: &str, size: f64) -> f64 {
    data_beans::figure::text_width(s, size)
}

/// How wide a position label `s` (digits and a point) sets at `size`
/// points: [`set_width`], a point or comma half a digit's width.
pub fn tick_width(s: &str, size: f64) -> f64 {
    set_width(s, size)
}

/// A character's width at `size` points, Helvetica's advance rounded up
/// (a name in capitals sets wider than one in lower case, so one width
/// for all would leave the one short and the other running over).
pub fn glyph_width(c: char, size: f64) -> f64 {
    let em = match c {
        'i' | 'j' | 'l' | 'I' | '.' | ',' | ':' | ';' | '\'' | '!' | '|' | ' ' => 0.28,
        'f' | 't' | 'r' | '(' | ')' | '[' | ']' | '-' | '/' => 0.34,
        'c' | 'k' | 's' | 'v' | 'x' | 'y' | 'z' | 'J' => 0.5,
        'm' | 'M' => 0.84,
        'w' => 0.73,
        'W' => 0.95,
        'L' => 0.56,
        'F' | 'T' | 'Z' => 0.62,
        'C' | 'D' | 'H' | 'N' | 'R' | 'U' => 0.73,
        'G' | 'O' | 'Q' => 0.78,
        '…' => 1.0,
        c if c.is_uppercase() => 0.67,
        _ => 0.56,
    };
    em * size
}

/// How wide `s` sets at `size` points, character by character (see
/// [`glyph_width`]): closer than [`text_width`] where texts sit side by
/// side.
pub fn set_width(s: &str, size: f64) -> f64 {
    s.chars().map(|c| glyph_width(c, size)).sum()
}

/// `s` as it fits `room` points across at `size`: as it is, else without
/// a parenthetical at its end, else without what follows a colon too,
/// else cut short with an ellipsis.
pub fn fit_text(s: &str, room: f64, size: f64) -> String {
    if text_width(s, size) <= room {
        return s.to_string();
    }
    let head = (s.strip_suffix(')').and_then(|t| t.rsplit_once(" ("))).map_or(s, |(head, _)| head);
    // Then what it is, without what is said of it after a colon.
    let what = head.split_once(": ").map_or(head, |(what, _)| what);
    for shorter in [head, what] {
        if !shorter.is_empty() && text_width(shorter, size) <= room {
            return shorter.to_string();
        }
    }
    let fits = ((room / (size * 0.55)).floor() as usize).saturating_sub(1);
    let cut: String = s.chars().take(fits).collect();
    format!("{}…", cut.trim_end())
}

/// How wide `s` sets in bold at `size`, roughly: Helvetica Bold's
/// capitals and digits, which a symbol is made of, set wider than
/// [`text_width`]'s average letter.
pub fn bold_width(s: &str, size: f64) -> f64 {
    let em = |ch: char| match ch {
        'I' | 'i' | 'j' | 'l' | '.' | ',' => 0.28,
        'M' | 'm' => 0.86,
        'W' => 0.94,
        'w' => 0.78,
        'A'..='Z' => 0.72,
        '0'..='9' => 0.556,
        'a'..='z' => 0.6,
        _ => 0.34,
    };
    s.chars().map(em).sum::<f64>() * size
}

/// Where a bold title `s` wanted centred at `at` goes in the box from `x`
/// to `x + w`, and the size it is written at: centred where it fits so,
/// else moved over as little as keeps it inside, else shrunk to the box.
/// Laid out at `size` on `c`.
pub fn place_bold_title(c: &Canvas, s: &str, at: f64, (x, w): (f64, f64), size: f64) -> (f64, f64) {
    let size = c.label_size(size);
    let (lo, hi) = (x + 1.0, x + w - 1.0);
    let wide = bold_width(s, size);
    let size = match wide > hi - lo {
        true => size * (hi - lo) / wide,
        false => size,
    };
    let half = bold_width(s, size) / 2.0;
    (at.clamp(lo + half, (hi - half).max(lo + half)), size)
}

/// Where a title `s` wanted centred at `at` goes in the box from `x` to
/// `x + w`, and what of it: centred where it fits so; else moved over as
/// little as keeps it inside; else shortened to the box (see
/// [`fit_text`]). Laid out at `size` on `c`.
pub fn place_title(c: &Canvas, s: &str, at: f64, (x, w): (f64, f64), size: f64) -> (f64, String) {
    let size = c.label_size(size);
    let (lo, hi) = (x + 2.0, x + w - 2.0);
    let text = fit_text(s, hi - lo, size);
    let half = text_width(&text, size) / 2.0;
    (at.clamp(lo + half, (hi - half).max(lo + half)), text)
}

/// A bar panel matching a terminal histogram: one bar per value, an optional
/// front subset over faint full bars, accent bars, x tick labels by bar, a
/// dashed pointer, site marks under the axis, and dashed dividers between
/// spans of bars.
pub struct Bars<'a> {
    pub values: &'a [f64],
    pub front: Option<&'a [f64]>,
    pub accent: &'a dyn Fn(usize) -> bool,
    /// Colour of a bar that is not accented, by bar; [`BAR`] for one tone.
    pub colour: &'a dyn Fn(usize) -> &'static str,
    pub y_scale: Scale,
    /// Top of the y axis; `None` scales to the tallest bar.
    pub y_max: Option<f64>,
    pub ticks: Vec<(usize, String)>,
    pub pointer: Option<usize>,
    pub marks: Vec<usize>,
    /// Full-height dashed lines on the left edge of these bars. With any,
    /// `ticks` name the spans between them and get no tick stub, which
    /// would read as one more boundary.
    pub dividers: Vec<usize>,
    pub title: String,
    pub x_title: String,
    pub y_title: String,
    /// Tick labels turned 90° counter-clockwise, for long ones (positions).
    pub rotated: bool,
}

/// The room around a bar panel, in points of its 8-point layout: left,
/// right, top and bottom, the top without a title when drawn alone.
/// A [`style::Style::tight`] panel takes only the room its labels need:
/// no title, the left as its widest y label (`widest`, laid out) and its
/// y title when `titled` (see [`tight_left`]), the bottom as its ticks and
/// the x title, little on the right.
fn bar_margins(
    c: &Canvas,
    (ticks, rotated): (&[(usize, String)], bool),
    widest: f64,
    titled: bool,
) -> [f64; 4] {
    let top = if c.style().panel_only || c.style().tight {
        8.0
    } else {
        20.0
    };
    let room = tick_room(ticks, rotated).max(c.style().tick_room);
    if c.style().tight {
        let bottom = match rotated {
            true => 25.0 + room,
            false => 28.0_f64.max(25.0 + c.style().tick_room),
        };
        return [tight_left(widest, titled), TIGHT_RIGHT, top, bottom];
    }
    [52.0, 12.0, top, 38.0 + room]
}

/// The room right of a tight panel's plot, laid out.
pub(crate) const TIGHT_RIGHT: f64 = 4.0;

/// A tight panel's left margin, laid out, for y labels `widest` across:
/// room for them and, when `titled`, the y title beside them; a panel
/// whose y axis goes untitled (one not first in a row of panels) gives
/// that room to its plot.
fn tight_left(widest: f64, titled: bool) -> f64 {
    match titled {
        true => (widest + 22.0).min(52.0),
        false => (widest + 8.0).min(38.0),
    }
}

/// The widest of `labels` as `c` writes 8-point ones, laid out.
fn widest(c: &Canvas, labels: &[&str]) -> f64 {
    let size = c.label_size(8.0);
    labels
        .iter()
        .map(|l| text_width(l, size))
        .fold(0.0, f64::max)
}

/// A bar panel's frame: the plot's box `(x, y, w, h)` inside the panel's
/// margins, and the palette's colours for what a threshold drops and for
/// its line.
struct BarFrame {
    plot: (f64, f64, f64, f64),
    dropped: &'static str,
    line: &'static str,
}

impl BarFrame {
    /// The frame of a bar panel in `(x, y, w, h)` with `ticks` (turned when
    /// `rotated`), its `title` drawn over it unless it is drawn alone.
    fn draw(
        c: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        ticks: (&[(usize, String)], bool),
        (title, y_title): (&str, &str),
        y_labels: &[&str],
    ) -> Self {
        let widest = widest(c, y_labels);
        let [left, right, top, bottom] = bar_margins(c, ticks, widest, !y_title.is_empty());
        // A tight page's rows go untitled: the figure's title and the
        // arms' names say what they are.
        if !c.style().panel_only && !c.style().tight {
            c.bold(x + left, y + 13.0, title, 11.0, Anchor::Start, INK);
        }
        let roles = c.style().palette.roles();
        Self {
            plot: (x + left, y + top, w - left - right, h - top - bottom),
            dropped: roles.dropped,
            line: roles.line,
        }
    }
}

/// The left margin a tight panel of `y_labels` takes on `c`, laid out,
/// with its y title when `titled`: where its plot starts.
fn left_of(c: &Canvas, y_labels: &[&str], titled: bool) -> f64 {
    tight_left(widest(c, y_labels), titled)
}

impl Bars<'_> {
    /// The y axis's top (scaled) and its labels at 0, half and full scale.
    fn y_axis(&self) -> (f64, [String; 3]) {
        let tallest = self
            .values
            .iter()
            .map(|&v| scaled(self.y_scale, v))
            .fold(0.0, f64::max);
        let max = self
            .y_max
            .map_or(tallest, |m| scaled(self.y_scale, m).max(tallest));
        let labels = [0.0, 0.5, 1.0].map(|f| {
            if max > 0.0 {
                compact(unscaled(self.y_scale, f * max))
            } else {
                "0".to_string()
            }
        });
        (max, labels)
    }

    /// Where its plot starts in a tight panel `h` deep on `c`, laid out:
    /// past the labels it writes at that depth (see [`Bars::written`]).
    pub fn left(&self, c: &Canvas, h: f64) -> f64 {
        let written = self.written(c, h);
        let labels: Vec<&str> = written.iter().map(|w| w.1.as_str()).collect();
        left_of(c, &labels, !self.y_title.is_empty())
    }

    /// The y labels it writes in a panel `h` deep, each with its share of
    /// the plot's depth: on a square-root or log axis the round values that
    /// stand apart there (see [`kept_round`]) over zero; else zero, half
    /// and full scale.
    fn written(&self, c: &Canvas, h: f64) -> Vec<(f64, String)> {
        let (max, labels) = self.y_axis();
        let titled = !self.y_title.is_empty();
        let [_, _, top, bottom] = bar_margins(c, (&self.ticks, self.rotated), 0.0, titled);
        let kept = match self.y_scale {
            Scale::Linear => Vec::new(),
            scale => kept_round(c, scale, max, h - top - bottom),
        };
        match kept.is_empty() {
            true => [0.0, 0.5, 1.0].into_iter().zip(labels).collect(),
            false => std::iter::once((0.0, "0".to_string()))
                .chain(kept)
                .collect(),
        }
    }

    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) -> (f64, f64, f64, f64) {
        let (max, _) = self.y_axis();
        let written = self.written(c, h);
        let labels: Vec<&str> = written.iter().map(|w| w.1.as_str()).collect();
        let BarFrame {
            plot: (px, py, pw, ph),
            dropped,
            line,
        } = BarFrame::draw(
            c,
            (x, y, w, h),
            (&self.ticks, self.rotated),
            (&self.title, &self.y_title),
            &labels,
        );
        let n = self.values.len().max(1);
        let bw = pw / n as f64;
        let height = |v: f64| {
            if max <= 0.0 {
                0.0
            } else {
                scaled(self.y_scale, v) / max * ph
            }
        };

        // y gridlines and labels where they were chosen.
        for (f, label) in &written {
            let f = *f;
            let gy = py + ph - f * ph;
            if f > 0.0 {
                c.line(px, gy, px + pw, gy, FAINT, 0.4);
            }
            c.text(px - 5.0, gy + 3.0, label, 8.0, Anchor::End, MUTED);
        }
        y_title(c, &self.y_title, x, (py, ph));

        let gap = if bw > 3.0 { 0.15 * bw } else { 0.0 };
        let mut bars = |values: &[f64], colour: &dyn Fn(usize) -> &'static str| {
            for (i, &v) in values.iter().enumerate() {
                let bh = height(v);
                if bh > 0.0 {
                    let bh = bh.max(MIN_BAR);
                    let (bx, w) = bar_span(px, bw, gap, i);
                    c.rect(bx, py + ph - bh, w, bh, colour(i));
                }
            }
        };
        let tone = |i: usize| {
            if (self.accent)(i) {
                dropped
            } else {
                (self.colour)(i)
            }
        };
        match self.front {
            Some(front) => {
                bars(self.values, &|_| FAINT);
                bars(front, &tone);
            }
            None => bars(self.values, &tone),
        }

        if let Some(p) = self.pointer {
            let cx = px + (p as f64 + 0.5) * bw;
            c.dashed(cx, py, cx, py + ph, line, 0.8);
        }
        c.line(px, py + ph, px + pw, py + ph, INK, 0.6);
        c.line(px, py, px, py + ph, INK, 0.6);
        for &m in &self.marks {
            let mx = px + (m as f64 + 0.5) * bw;
            c.line(mx, py + ph, mx, py + ph + 3.0, MUTED, 0.5);
        }
        for &d in &self.dividers {
            let dx = px + d as f64 * bw;
            c.dashed(dx, py, dx, py + ph + 4.0, MUTED, 0.8);
        }
        let stub = self.dividers.is_empty();
        let at = (px, bw, py + ph, (x, x + w));
        x_ticks(c, &self.ticks, at, stub, self.rotated);
        let (at, x_title) = place_title(c, &self.x_title, px + pw / 2.0, (x, w), 8.0);
        c.text(at, y + h - 6.0, &x_title, 8.0, Anchor::Middle, MUTED);
        (px, py, pw, ph)
    }
}

/// A panel's y title `s` at `x`, along its plot from `py`, `ph` tall: on a
/// tight page shortened where it would run past; none when empty.
fn y_title(c: &mut Canvas, s: &str, x: f64, (py, ph): (f64, f64)) {
    if s.is_empty() {
        return;
    }
    let s = match c.style().tight {
        true => fit_text(s, ph, c.label_size(8.0)),
        false => s.to_string(),
    };
    c.vtext(x + 12.0, py + ph / 2.0, &s, 8.0, MUTED);
}

/// A bar panel ([`Bars`], [`Diverging`], [`Mirror`]) as a figure lays it
/// out: where its plot starts on a tight page, the position ticks it is
/// given for that, and its drawing.
pub trait Panel {
    /// Where its plot starts in a tight panel `h` deep on `c`, laid out.
    fn left(&self, c: &Canvas, h: f64) -> f64;
    /// Its position ticks, turned or flat.
    fn set_ticks(&mut self, ticks: Vec<(usize, String)>, rotated: bool);
    /// Drawn in the box `(x, y, w, h)`; the box its bars take.
    fn draw_in(&self, c: &mut Canvas, bbox: (f64, f64, f64, f64)) -> (f64, f64, f64, f64);
}

macro_rules! panel {
    ($($t:ident),*) => {$(
        impl Panel for $t<'_> {
            fn left(&self, c: &Canvas, h: f64) -> f64 {
                $t::left(self, c, h)
            }
            fn set_ticks(&mut self, ticks: Vec<(usize, String)>, rotated: bool) {
                self.ticks = ticks;
                self.rotated = rotated;
            }
            fn draw_in(&self, c: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64)) -> (f64, f64, f64, f64) {
                self.draw(c, x, y, w, h)
            }
        }
    )*};
}
panel!(Bars, Diverging, Mirror);

/// The least a bar is drawn across, and one of any value up, laid out: a
/// hairline, so a bar on a small page, or of a few reads beside many,
/// still shows, as a terminal's block does.
pub const MIN_BAR: f64 = 0.75;

/// Where bar `i` of bars `bw` apart from `px` is drawn, `gap` between
/// them, and how wide: at least [`MIN_BAR`], centred on its place.
fn bar_span(px: f64, bw: f64, gap: f64, i: usize) -> (f64, f64) {
    let w = (bw - gap).max(MIN_BAR);
    (px + (i as f64 + 0.5) * bw - w / 2.0, w)
}

/// Signed bars around a zero line: positive up in the accent, negative
/// down in grey; `None` draws no bar. Ticks and pointer as [`Bars`].
pub struct Diverging<'a> {
    pub values: &'a [Option<f64>],
    pub ticks: Vec<(usize, String)>,
    pub pointer: Option<usize>,
    pub title: String,
    pub x_title: String,
    pub y_title: String,
    /// Label for a value on the y axis.
    pub label: &'a dyn Fn(f64) -> String,
    /// Tick labels turned, as [`Bars::rotated`].
    pub rotated: bool,
}

impl Diverging<'_> {
    /// Where its plot starts in a tight panel `h` deep on `c`, laid out.
    pub fn left(&self, c: &Canvas, h: f64) -> f64 {
        let (up, down) = split_signed(self.values);
        self.mirror(&up, &down).left(c, h)
    }

    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) -> (f64, f64, f64, f64) {
        let (up, down) = split_signed(self.values);
        self.mirror(&up, &down).draw(c, x, y, w, h)
    }

    /// The mirror it draws: `up` in the accent over `down` in grey.
    fn mirror<'b>(&'b self, up: &'b [f64], down: &'b [f64]) -> Mirror<'b> {
        let max = up.iter().chain(down).fold(0.0, |m: f64, &v| m.max(v));
        Mirror {
            up: Half::plain(up, ACCENT),
            down: Half::plain(down, BAR),
            y_scale: Scale::Linear,
            y_max: Some(max),
            y_labels: [(self.label)(max), (self.label)(0.0), (self.label)(-max)],
            ticks: self.ticks.clone(),
            pointer: self.pointer,
            title: self.title.clone(),
            x_title: self.x_title.clone(),
            y_title: self.y_title.clone(),
            rotated: self.rotated,
        }
    }
}

/// One arm of a [`Mirror`], for placing its name: the upper or the lower,
/// the y labels' lines across the plot, its bars (each `bw` across, `n` of
/// them) and how far bar `i` reaches from the zero line, laid out.
struct Arm<'a> {
    upper: bool,
    lines: &'a [f64],
    bars: (f64, usize),
    reach: &'a dyn Fn(usize) -> f64,
}

impl Arm<'_> {
    /// Where the arm's `name` is written on the plot `(px, py, pw, ph)`:
    /// its anchor point and anchor, and whether it is clear of the bars.
    /// Always at the left, so every panel's names are in one place: at the
    /// arm's outer corner, clear of the y labels' lines and of the bars,
    /// else nudged toward the zero line; where nothing there is clear of
    /// both, clear of the lines where it overlaps the bars least (drawn
    /// over a halo); the corner where nothing is clear of the lines.
    fn name_at(&self, c: &Canvas, name: &str, (px, py, ph): (f64, f64, f64)) -> (f64, bool) {
        let size = c.label_size(8.0);
        // Measured glyph by glyph: an arm's name is mostly capitals, which
        // set wider than the average letter `text_width` reckons, and must
        // not touch a bar it is said to be clear of.
        let wide = set_width(name, size);
        let (above, below) = (0.75 * size, 0.25 * size);
        let zero = py + ph / 2.0;
        let corner = if self.upper { py + 9.0 } else { py + ph - 3.0 };
        // Baselines from the corner toward the zero line, a point apart.
        let room = if self.upper {
            zero - 2.0 - below - corner
        } else {
            corner - above - zero - 2.0
        };
        let baselines = (0..=room.max(0.0).floor() as usize).map(|k| match self.upper {
            true => corner + k as f64,
            false => corner - k as f64,
        });
        let left = px + NAME_INSET;
        let off_lines = |t: f64, b: f64| self.lines.iter().all(|&g| g <= t - 1.0 || g >= b + 1.0);
        // How far into the bars a name from `l` to `l + wide`, `t` to `b`,
        // runs: the most any under it, or one either side (a hairline bar
        // is wider than its place), reaches past it.
        let over = |l: f64, t: f64, b: f64| {
            let (bw, n) = self.bars;
            if bw <= 0.0 || n == 0 {
                return 0.0;
            }
            let first = (((l - px) / bw).floor() as i64 - 1).max(0) as usize;
            let last = (((l + wide - px) / bw).floor() as i64 + 1).clamp(0, n as i64 - 1) as usize;
            (first..=last)
                .map(|i| match (self.reach)(i) {
                    r if r <= 0.0 => 0.0,
                    r if self.upper => b + 1.0 - (zero - r),
                    r => zero + r + 1.0 - t,
                })
                .fold(0.0, f64::max)
        };
        let mut least: Option<(f64, f64)> = None;
        for y in baselines {
            let (t, b) = (y - above, y + below);
            if !off_lines(t, b) {
                continue;
            }
            let into = over(left, t, b);
            if into <= 0.0 {
                return (y, true);
            }
            if least.is_none_or(|m| into < m.0) {
                least = Some((into, y));
            }
        }
        (least.map_or(corner, |m| m.1), false)
    }
}

/// How far in from the plot's left an arm's name starts, laid out.
const NAME_INSET: f64 = 4.0;

/// Height of an 8-point character's advance, roughly, for laying out
/// turned labels.
const CHAR_W: f64 = 4.6;

/// Room below the axis that turned tick labels take beyond flat ones.
pub(crate) fn tick_room(ticks: &[(usize, String)], rotated: bool) -> f64 {
    if !rotated {
        return 0.0;
    }
    let longest = ticks
        .iter()
        .map(|(_, l)| l.chars().count())
        .max()
        .unwrap_or(0);
    (longest as f64 * CHAR_W - 6.0).max(0.0)
}

/// Where a flat tick label `half` across either side of its middle is
/// centred, for a tick at `tx`, inside `(lo, hi)`: on the tick, else moved
/// in as little as keeps it inside while the tick stays under it; none
/// where it cannot.
pub(crate) fn flat_label(tx: f64, half: f64, (lo, hi): (f64, f64)) -> Option<f64> {
    if hi - lo < 2.0 * half {
        return None;
    }
    let at = tx.clamp(lo + half, hi - half);
    ((at - tx).abs() <= 0.75 * half).then_some(at)
}

/// The x ticks under the axis at `base`, bar `i` at `px + (i + ½)·bw`:
/// a stub (with `stub`) and the label, flat and centred, or turned 90°
/// counter-clockwise to end just under the stub. A flat label at an end
/// moves in to stay inside the panel's `span` across (see [`flat_label`]).
/// Thinned to fit: a label that would run into the one before it, or out
/// of the span, is left out, stub and all.
fn x_ticks(
    c: &mut Canvas,
    ticks: &[(usize, String)],
    (px, bw, base, span): (f64, f64, f64, (f64, f64)),
    stub: bool,
    rotated: bool,
) {
    let size = c.label_size(8.0);
    let mut last = f64::NEG_INFINITY;
    for (i, label) in ticks {
        let tx = px + (*i as f64 + 0.5) * bw;
        let (at, (lo, hi)) = if rotated {
            (tx, (tx + 3.0 - 0.75 * size, tx + 3.0 + 0.25 * size))
        } else {
            let half = tick_width(label, size) / 2.0;
            let Some(at) = flat_label(tx, half, span) else {
                continue;
            };
            (at, (at - half, at + half))
        };
        if lo < last + 0.25 * size || lo < span.0 || hi > span.1 {
            continue;
        }
        last = hi;
        if stub {
            c.line(tx, base, tx, base + 4.0, INK, 0.6);
        }
        if rotated {
            c.rotated(tx + 3.0, base + 6.0, label, 8.0, Anchor::End, INK);
        } else {
            c.text(at, base + 13.0, label, 8.0, Anchor::Middle, INK);
        }
    }
}

/// Signed values split into their positive and negative parts, as sizes;
/// `None` is zero on both.
pub fn split_signed(values: &[Option<f64>]) -> (Vec<f64>, Vec<f64>) {
    let part = |sign: f64| {
        let v = values.iter();
        v.map(|v| v.map_or(0.0, |v| (sign * v).max(0.0))).collect()
    };
    (part(1.0), part(-1.0))
}

/// One side of a [`Mirror`]: bars, optionally with a part drawn in front
/// (in the accent, the bars behind then faint), and a name in the corner.
pub struct Half<'a> {
    pub values: &'a [f64],
    pub front: Option<&'a [f64]>,
    pub colour: &'static str,
    pub name: String,
}

impl<'a> Half<'a> {
    pub fn plain(values: &'a [f64], colour: &'static str) -> Self {
        Half {
            values,
            front: None,
            colour,
            name: String::new(),
        }
    }
}

/// Two bar series on one scale around a zero line, one growing up and the
/// other down, as a Miami plot. Ticks and pointer as [`Bars`].
pub struct Mirror<'a> {
    pub up: Half<'a>,
    pub down: Half<'a>,
    pub y_scale: Scale,
    /// Top of either side; `None` scales to the tallest bar.
    pub y_max: Option<f64>,
    /// Axis labels at the top, the zero line and the bottom; empty ones are
    /// the scale's own.
    pub y_labels: [String; 3],
    pub ticks: Vec<(usize, String)>,
    pub pointer: Option<usize>,
    pub title: String,
    pub x_title: String,
    pub y_title: String,
    /// Tick labels turned, as [`Bars::rotated`].
    pub rotated: bool,
}

impl Mirror<'_> {
    /// The y axis: its top (scaled) and its labels at the top, the zero
    /// line and the bottom.
    fn y_axis(&self) -> (f64, [String; 3]) {
        let all = self.up.values.iter().chain(self.down.values);
        let tallest = all.map(|&v| scaled(self.y_scale, v)).fold(0.0, f64::max);
        let max = self
            .y_max
            .map_or(tallest, |m| scaled(self.y_scale, m).max(tallest));
        let own = compact(unscaled(self.y_scale, max));
        let labels = [0, 1, 2].map(|k| match self.y_labels[k].as_str() {
            "" if k == 1 => "0".to_string(),
            "" => own.clone(),
            l => l.to_string(),
        });
        (max, labels)
    }

    /// Where its plot starts in a tight panel `h` deep on `c`, laid out:
    /// past the labels it writes at that depth (see [`Mirror::written`]).
    pub fn left(&self, c: &Canvas, h: f64) -> f64 {
        let (written, _) = self.written(c, h);
        let written: Vec<&str> = written.iter().map(String::as_str).collect();
        left_of(c, &written, !self.y_title.is_empty())
    }

    /// The y labels it writes in a panel `h` deep, and the round values
    /// among them (see [`Mirror::kept_round`]): its plot's depth does not
    /// hang on their width, so they are chosen at that depth (the round
    /// values each arm holds, else its own top, zero and bottom).
    fn written(&self, c: &Canvas, h: f64) -> (Vec<String>, Vec<(f64, String)>) {
        let titled = !self.y_title.is_empty();
        let [_, _, top, bottom] = bar_margins(c, (&self.ticks, self.rotated), 0.0, titled);
        let kept = self.kept_round(c, h - top - bottom);
        let written = match kept.is_empty() {
            true => self.y_axis().1.to_vec(),
            false => std::iter::once("0".to_string())
                .chain(kept.iter().map(|k| k.1.clone()))
                .collect(),
        };
        (written, kept)
    }

    /// The round values its y axis labels up each arm on a plot `ph`
    /// deep, each a label's height and a little apart, from the top down
    /// and none past the zero's (a square-root axis's chosen so, to stay
    /// evenly spread).
    fn kept_round(&self, c: &Canvas, ph: f64) -> Vec<(f64, String)> {
        // A square-root or log axis of its own labels; else its top and
        // zero.
        let own_labels = self.y_labels.iter().all(String::is_empty);
        match self.y_scale {
            Scale::Linear => Vec::new(),
            _ if !own_labels => Vec::new(),
            scale => kept_round(c, scale, self.y_axis().0, ph / 2.0),
        }
    }

    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) -> (f64, f64, f64, f64) {
        let (max, labels) = self.y_axis();
        // The left margin made for the labels written at the plot's depth
        // (the same as the frame's: its margins' top and bottom do not hang
        // on the labels), so none runs past the panel's edge.
        let (written, kept) = self.written(c, h);
        let written: Vec<&str> = written.iter().map(String::as_str).collect();
        let BarFrame {
            plot: (px, py, pw, ph),
            line,
            ..
        } = BarFrame::draw(
            c,
            (x, y, w, h),
            (&self.ticks, self.rotated),
            (&self.title, &self.y_title),
            &written,
        );
        let front_tone = c.style().palette.roles().front;
        let n = self.up.values.len().max(self.down.values.len()).max(1);
        let bw = pw / n as f64;
        let zero = py + ph / 2.0;
        let height = |v: f64| {
            if max <= 0.0 {
                0.0
            } else {
                scaled(self.y_scale, v) / max * ph / 2.0
            }
        };
        // Where the y labels' lines run, the zero line's among them.
        let mut lines = vec![zero];
        match kept.is_empty() {
            true => {
                for (f, label) in [1.0, 0.0, -1.0].into_iter().zip(&labels) {
                    let gy = zero - f * ph / 2.0;
                    if f != 0.0 {
                        c.line(px, gy, px + pw, gy, FAINT, 0.4);
                        lines.push(gy);
                    }
                    c.text(px - 5.0, gy + 3.0, label, 8.0, Anchor::End, MUTED);
                }
            }
            false => {
                c.text(px - 5.0, zero + 3.0, "0", 8.0, Anchor::End, MUTED);
                for (share, label) in &kept {
                    for sign in [-1.0, 1.0] {
                        let gy = zero + sign * share * ph / 2.0;
                        c.line(px, gy, px + pw, gy, FAINT, 0.4);
                        lines.push(gy);
                        c.text(px - 5.0, gy + 3.0, label, 8.0, Anchor::End, MUTED);
                    }
                }
            }
        }
        y_title(c, &self.y_title, x, (py, ph));
        let gap = if bw > 3.0 { 0.15 * bw } else { 0.0 };
        for (half, sign) in [(&self.up, -1.0), (&self.down, 1.0)] {
            let mut bars = |values: &[f64], colour: &'static str| {
                for (i, &v) in values.iter().enumerate() {
                    // A bar of anything shows, as the browser shows it.
                    let bh = match height(v) {
                        h if h > 0.0 => h.max(MIN_BAR),
                        _ => 0.0,
                    };
                    let (bx, w) = bar_span(px, bw, gap, i);
                    let top = if sign < 0.0 { zero - bh } else { zero };
                    c.rect(bx, top, w, bh, colour);
                }
            };
            match half.front {
                Some(front) => {
                    bars(half.values, FAINT);
                    bars(front, front_tone);
                }
                None => bars(half.values, half.colour),
            }
        }
        let plot = (px, py, pw, ph);
        for (half, upper) in [(&self.up, true), (&self.down, false)] {
            // How far each bar of the arm reaches from the zero line.
            let reach = |i: usize| {
                let v = half.values.get(i).copied().unwrap_or(0.0);
                let f = half.front.and_then(|f| f.get(i)).copied().unwrap_or(0.0);
                match height(v.max(f)) {
                    h if h > 0.0 => h.max(MIN_BAR),
                    _ => 0.0,
                }
            };
            let arm = Arm {
                upper,
                lines: &lines,
                bars: (bw, n),
                reach: &reach,
            };
            let nx = plot.0 + NAME_INSET;
            match arm.name_at(c, &half.name, (plot.0, plot.1, plot.3)) {
                (ny, true) => c.text(nx, ny, &half.name, 8.0, Anchor::Start, MUTED),
                (ny, false) => c.text_halo(nx, ny, &half.name, 8.0, Anchor::Start, MUTED),
            }
        }
        if let Some(p) = self.pointer {
            let cx = px + (p as f64 + 0.5) * bw;
            c.dashed(cx, py, cx, py + ph, line, 0.8);
        }
        c.line(px, zero, px + pw, zero, INK, 0.6);
        c.line(px, py, px, py + ph, INK, 0.6);
        let at = (px, bw, py + ph, (x, x + w));
        x_ticks(c, &self.ticks, at, true, self.rotated);
        let (at, x_title) = place_title(c, &self.x_title, px + pw / 2.0, (x, w), 8.0);
        c.text(at, y + h - 6.0, &x_title, 8.0, Anchor::Middle, MUTED);
        (px, py, pw, ph)
    }
}

fn fontdb() -> Arc<usvg::fontdb::Database> {
    static DB: OnceLock<Arc<usvg::fontdb::Database>> = OnceLock::new();
    DB.get_or_init(|| {
        let mut db = usvg::fontdb::Database::new();
        db.load_system_fonts();
        Arc::new(db)
    })
    .clone()
}

/// `prefix` without a trailing `.pdf` or `.png`, so either spelling works.
fn strip_extension(prefix: &str) -> &str {
    prefix
        .strip_suffix(".pdf")
        .or_else(|| prefix.strip_suffix(".png"))
        .unwrap_or(prefix)
}

/// The files a save under `prefix` writes: the PDF, then the PNG.
pub(crate) fn saved_paths(prefix: &str) -> [String; 2] {
    let prefix = strip_extension(prefix);
    [format!("{prefix}.pdf"), format!("{prefix}.png")]
}

/// The PDF a save under `prefix` writes.
fn pdf_path(prefix: &str) -> String {
    format!("{}.pdf", strip_extension(prefix))
}

/// Whether a save under `prefix` would replace a file already there.
fn taken(prefix: &str) -> bool {
    saved_paths(prefix)
        .iter()
        .any(|p| std::path::Path::new(p).exists())
}

/// `prefix`, or failing that the first of `{prefix}_2`, `{prefix}_3`, …
/// that no save would replace a file under.
fn free_name(prefix: &str) -> String {
    let prefix = strip_extension(prefix);
    (1..)
        .map(|n| match n {
            1 => prefix.to_string(),
            n => format!("{prefix}_{n}"),
        })
        .find(|p| !taken(p))
        .expect("some suffix is free")
}

/// SVG parsing options, with the system fonts loaded once.
fn options() -> usvg::Options<'static> {
    usvg::Options {
        fontdb: fontdb(),
        ..usvg::Options::default()
    }
}

/// Write `svg` as `{prefix}.pdf` (vector) alone; returns the path written.
pub fn save(svg: &str, prefix: &str) -> anyhow::Result<Vec<String>> {
    let tree = usvg::Tree::from_str(svg, &options())?;
    Ok(vec![write_pdf(&tree, prefix)?])
}

/// `svg` as PNG bytes, `scale` pixels per point, on white.
pub fn png(svg: &str, scale: f32) -> anyhow::Result<Vec<u8>> {
    png_of(&usvg::Tree::from_str(svg, &options())?, scale)
}

/// `tree` as PNG bytes, `scale` pixels per point, on white.
fn png_of(tree: &usvg::Tree, scale: f32) -> anyhow::Result<Vec<u8>> {
    let pixmap = term::pixmap(tree, scale)?;
    let (w, h) = (pixmap.width(), pixmap.height());
    let picture = image::RgbaImage::from_raw(w, h, pixmap.take())
        .ok_or_else(|| anyhow::anyhow!("pixel buffer size mismatch"))?;
    let mut out = std::io::Cursor::new(Vec::new());
    picture.write_to(&mut out, image::ImageFormat::Png)?;
    Ok(out.into_inner())
}

/// `svg` as `{prefix}.pdf` (vector, its labels as text) and `{prefix}.png`
/// at two pixels per point beside it, from one parse: what `--figure`
/// writes without a view; the paths written.
pub fn save_with_png(svg: &str, prefix: &str) -> anyhow::Result<Vec<String>> {
    let tree = usvg::Tree::from_str(svg, &options())?;
    let [_, png_path] = saved_paths(prefix);
    let pdf = write_pdf(&tree, prefix)?;
    std::fs::write(&png_path, png_of(&tree, 2.0)?)?;
    Ok(vec![pdf, png_path])
}

/// A file a figure is written as (`--format`; MCP's `formats`).
#[derive(
    clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    /// The drawing itself.
    Svg,
    /// Vector, its labels as text: what `s` saves.
    Pdf,
    /// Two pixels per point, on white.
    Png,
}

impl Format {
    /// What its files end in.
    pub fn extension(self) -> &'static str {
        match self {
            Format::Svg => "svg",
            Format::Pdf => "pdf",
            Format::Png => "png",
        }
    }
}

/// `svg` as `{prefix}.{format}` for each of `formats`, the PDF as `s`
/// saves it; refused, nothing written, when one of them is there. The
/// paths written, in `formats`' order.
pub fn save_as(svg: &str, prefix: &str, formats: &[Format]) -> anyhow::Result<Vec<String>> {
    let prefix = strip_extension(prefix);
    let paths: Vec<String> = (formats.iter())
        .map(|f| format!("{prefix}.{}", f.extension()))
        .collect();
    if let Some(p) = paths.iter().find(|p| std::path::Path::new(p).exists()) {
        anyhow::bail!("{p} exists; not replacing it");
    }
    let tree = usvg::Tree::from_str(svg, &options())?;
    for (f, path) in formats.iter().zip(&paths) {
        match f {
            Format::Svg => std::fs::write(path, svg)?,
            Format::Pdf => {
                write_pdf(&tree, prefix)?;
            }
            Format::Png => std::fs::write(path, png_of(&tree, 2.0)?)?,
        }
    }
    Ok(paths)
}

/// The PDF of `tree` under `prefix`.
fn write_pdf(tree: &usvg::Tree, prefix: &str) -> anyhow::Result<String> {
    let path = pdf_path(prefix);
    let pdf = svg2pdf::to_pdf(tree, Default::default(), Default::default())
        .map_err(|e| anyhow::anyhow!("PDF conversion failed: {e:?}"))?;
    std::fs::write(&path, pdf)?;
    Ok(path)
}

/// `svg` as `{prefix}.pdf`, and a picture of it for the gallery's thumbnail, drawn in
/// memory only.
fn save_picture(svg: &str, prefix: &str) -> anyhow::Result<(Vec<String>, image::RgbaImage)> {
    let tree = usvg::Tree::from_str(svg, &options())?;
    let pdf = write_pdf(&tree, prefix)?;
    let pixmap = term::pixmap(&tree, 1.0)?;
    let (w, h) = (pixmap.width(), pixmap.height());
    let picture = image::RgbaImage::from_raw(w, h, pixmap.take())
        .ok_or_else(|| anyhow::anyhow!("pixel buffer size mismatch"))?;
    Ok((vec![pdf], picture))
}

/// What a key did to a [`LineInput`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Edit {
    /// The line is still open (or the key was not for it).
    Typing,
    Cancelled,
    /// Enter, with the trimmed text.
    Submitted(String),
}

/// A one-line text input: typing, Backspace, Esc to cancel, Enter to submit.
pub struct LineInput {
    buf: Option<String>,
    max: usize,
}

impl LineInput {
    pub fn new(max: usize) -> Self {
        Self { buf: None, max }
    }

    pub fn open(&mut self, initial: &str) {
        self.buf = Some(initial.to_string());
    }

    pub fn text(&self) -> Option<&str> {
        self.buf.as_deref()
    }

    pub fn active(&self) -> bool {
        self.buf.is_some()
    }

    pub fn handle(&mut self, key: KeyEvent) -> Edit {
        let Some(buf) = self.buf.as_mut() else {
            return Edit::Typing;
        };
        // The apply key submits too.
        let code = match key.code {
            _ if crate::tui::is_apply(&key) => KeyCode::Enter,
            code => code,
        };
        match code {
            KeyCode::Char(ch) if buf.len() < self.max => buf.push(ch),
            KeyCode::Backspace => {
                buf.pop();
            }
            KeyCode::Esc => {
                self.buf = None;
                return Edit::Cancelled;
            }
            KeyCode::Enter => {
                let text = self.buf.take().unwrap_or_default();
                return Edit::Submitted(text.trim().to_string());
            }
            _ => {}
        }
        Edit::Typing
    }
}

/// The save-as prompt every view shares: it opens with a suggested name
/// no save has used here yet, Enter confirms, Esc closes; the outcome shows
/// until the next key. It never replaces a file: a name already taken
/// keeps the prompt open, offering the next free one.
pub struct SavePrompt {
    input: LineInput,
    status: Option<Result<String, String>>,
    default: String,
    /// Why the prompt is still open: the name tried was taken.
    taken: Option<String>,
}

impl SavePrompt {
    pub fn new(default: &str) -> Self {
        Self {
            input: LineInput::new(512),
            status: None,
            default: default.to_string(),
            taken: None,
        }
    }

    pub fn active(&self) -> bool {
        self.input.active()
    }

    pub fn open(&mut self) {
        self.input.open(&free_name(&self.default));
        self.status = None;
        self.taken = None;
    }

    /// Handle a key while open: `Some(prefix)` once the user confirms.
    pub fn handle(&mut self, key: KeyEvent) -> Option<String> {
        match self.input.handle(key) {
            Edit::Submitted(name) if taken(&name) => {
                self.input.open(&free_name(&name));
                self.taken = Some(format!("{} is already there", strip_extension(&name)));
                None
            }
            Edit::Submitted(name) if !name.is_empty() => {
                self.default = name.clone();
                self.taken = None;
                Some(name)
            }
            _ => None,
        }
    }

    /// Record how a save went, for the footer.
    pub fn report(&mut self, result: anyhow::Result<Vec<String>>) {
        self.status = Some(match result {
            Ok(paths) => Ok(format!("saved {}", paths.join(", "))),
            Err(e) => Err(format!("save failed: {e}")),
        });
    }

    /// Clear a shown outcome (on the next key).
    pub fn dismiss(&mut self) {
        self.status = None;
    }

    /// The name being typed, while the prompt is open.
    pub fn text(&self) -> Option<&str> {
        self.input.text()
    }

    /// Why the prompt offers another name, after a taken one was tried.
    pub fn taken(&self) -> Option<&str> {
        self.taken.as_deref()
    }

    /// The footer while the prompt or an outcome is showing.
    pub fn footer(&self) -> Option<Line<'static>> {
        if let Some(buf) = self.input.text() {
            return Some(input_line(
                "save as (.pdf): ",
                buf,
                &[("Enter", "save"), ("Esc", "back")],
            ));
        }
        self.status.as_ref().map(|s| match s {
            Ok(msg) => Line::from(Span::styled(format!(" {msg}"), DIM)),
            Err(msg) => status_line(msg),
        })
    }
}

/// A footer message in the accent colour.
pub fn status_line(msg: &str) -> Line<'static> {
    Line::from(Span::styled(format!(" {msg}"), HIGHLIGHT))
}

/// What [`Controls::key`] made of a key.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Key {
    /// Not a controls key: the view handles it.
    Pass,
    /// Used by the controls; nothing for the view to do.
    Used,
    /// A save was confirmed under this name: draw the figure and call
    /// [`Controls::save`].
    Save(String),
}

/// What every view's figure needs: the save prompt (`s`), which shows the
/// figures saved here beside the name, and drawing plots as images where the
/// terminal can (`i` switches to text and back).
pub struct Controls {
    save: SavePrompt,
    picker: Option<Picker>,
    /// Images wanted (`i`), where the terminal has them.
    images: bool,
    /// Ask the attached terminal for images ([`Controls::detect`]).
    detecting: bool,
    /// The attach ([`legume_tui::attaches`]) the picker was asked at;
    /// none before the first ask.
    picked_at: Option<u64>,
    /// Which view saves, for the gallery: the default name's first word.
    what: String,
    /// The saves here, read when the prompt first opens.
    gallery: Option<gallery::Gallery>,
    /// Thumbnails ready to draw, by thumbnail file; `None` where one failed.
    thumbs: std::collections::HashMap<std::path::PathBuf, Option<Protocol>>,
    /// Where a save goes, as the prompt says, taken when it opens.
    here: String,
}

/// Rows a thumbnail takes in the save pop-up.
const THUMB_ROWS: u16 = 5;

impl Controls {
    /// Text plots only; [`Controls::detect`] asks the terminal for images.
    pub fn new(default_name: &str) -> Self {
        Self {
            save: SavePrompt::new(default_name),
            picker: None,
            images: true,
            detecting: false,
            picked_at: None,
            what: default_name
                .split_once('_')
                .map_or(default_name, |(view, _)| view)
                .into(),
            gallery: None,
            thumbs: Default::default(),
            here: String::new(),
        }
    }

    /// Use images if the terminal has them: asked of the terminal attached
    /// at the first draw ([`Controls::refresh`]), and again after each
    /// attach.
    pub fn detect(mut self) -> Self {
        self.detecting = true;
        self
    }

    #[cfg(test)]
    pub fn with_picker(mut self, picker: Picker) -> Self {
        self.picker = Some(picker);
        self.images = true;
        self
    }

    /// Ask the terminal attached now for images, once per attach: call it
    /// first in a view's `render`, with the screen up. A new terminal (or
    /// the same one back: kitty's images were deleted as it attached) has
    /// the thumbnails made again; plots made with the old picker are made
    /// again by [`term::PlotImage`] itself.
    pub fn refresh(&mut self) {
        if !self.detecting {
            return;
        }
        let attach = legume_tui::attaches();
        if self.picked_at == Some(attach) {
            return;
        }
        self.picked_at = Some(attach);
        self.picker = legume_tui::current().and_then(term::picker);
        self.thumbs.clear();
    }

    /// The size of a figure saved from a view that took `area`, as
    /// [`page`] says, its cells shaped as the terminal's font where the
    /// terminal has said, else as [`CELL`].
    pub fn page(&self, area: Rect, min: (f64, f64)) -> Option<(f64, f64)> {
        let cell = self
            .picker
            .as_ref()
            .map_or(CELL, |p| font_cell(p.font_size()));
        page(area, min, cell)
    }

    /// Handle `key` if it is for the save prompt, `s` or `i`.
    pub fn key(&mut self, key: KeyEvent) -> Key {
        if self.save.active() {
            return match self.save.handle(key) {
                Some(prefix) => Key::Save(prefix),
                None => Key::Used,
            };
        }
        self.save.dismiss();
        match key.code {
            KeyCode::Char('s') => {
                self.save.open();
                self.gallery.get_or_insert_with(gallery::Gallery::here);
                self.here = std::env::current_dir()
                    .map(|d| crate::tui::tilde(&d))
                    .unwrap_or_default();
            }
            KeyCode::Char('i') if self.picker.is_some() => self.images = !self.images_on(),
            _ => return Key::Pass,
        }
        Key::Used
    }

    /// Save `svg` under the confirmed `prefix`, list it with its thumbnail,
    /// and show how it went.
    pub fn save(&mut self, svg: &str, prefix: &str) {
        let saved = save_picture(svg, prefix).map(|(paths, picture)| {
            let gallery = self.gallery.get_or_insert_with(gallery::Gallery::here);
            if let Err(e) = gallery.add(std::path::Path::new(&paths[0]), &self.what, &picture) {
                log::debug!("not listed in the gallery: {e}");
            }
            paths
        });
        self.save.report(saved);
    }

    /// Save each `(suffix, svg)` of `files` under the confirmed `prefix`
    /// and the suffix, list the first with its thumbnail, and show how it
    /// went: every file, or the first that failed.
    pub fn save_many(&mut self, files: &[(String, String)], prefix: &str) {
        let prefix = strip_extension(prefix);
        let mut written = Vec::new();
        for (k, (suffix, svg)) in files.iter().enumerate() {
            match save_picture(svg, &format!("{prefix}{suffix}")) {
                Ok((paths, picture)) => {
                    if k == 0 {
                        let gallery = self.gallery.get_or_insert_with(gallery::Gallery::here);
                        let first = std::path::Path::new(&paths[0]);
                        if let Err(e) = gallery.add(first, &self.what, &picture) {
                            log::debug!("not listed in the gallery: {e}");
                        }
                    }
                    written.extend(paths);
                }
                Err(e) => return self.save.report(Err(e)),
            }
        }
        self.save.report(Ok(written));
    }

    /// While the save prompt is open, a pop-up over `area`: the figures
    /// saved here on the left, newest first, each with its thumbnail where
    /// the terminal draws images; the name being typed on the right. Call
    /// it last in a view's `render`.
    pub fn render_save(&mut self, frame: &mut Frame, area: Rect) {
        let Some(name) = self.save.text() else {
            return;
        };
        let w = area.width.saturating_sub(4).clamp(40, 100);
        let h = area.height.saturating_sub(2).clamp(8, 28);
        let inner = crate::tui::popup_frame(frame, area, w, h, " save the figure ".into());
        let entries = self.gallery.as_ref().map_or(&[][..], |g| g.entries());
        let strip = if entries.is_empty() {
            0
        } else {
            (inner.width / 2).min(36)
        };
        let [left, right] =
            Layout::horizontal([Constraint::Length(strip), Constraint::Fill(1)]).areas(inner);
        let lines = vec![
            Line::from(Span::styled(" save as (.pdf)", DIM)),
            Line::from(""),
            Line::from(Span::styled(format!(" {name}▏"), HIGHLIGHT)),
            match self.save.taken() {
                Some(why) => status_line(&format!("{why}; not replacing it")),
                None => Line::from(""),
            },
            Line::from(Span::styled(format!(" in {}", self.here), DIM)),
            Line::from(""),
            help_line(&[("Enter", "save"), ("Esc", "back")]),
        ];
        frame.render_widget(Paragraph::new(lines), right);
        let Some(gallery) = self.gallery.as_ref().filter(|_| strip > 0) else {
            return;
        };

        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs());
        let rows = if self.picker.is_some() { THUMB_ROWS } else { 0 };
        let each = rows + 3;
        let title = format!(" saved here ({})", entries.len());
        frame.render_widget(
            Paragraph::new(Line::from(Span::styled(title, DIM))),
            Rect::new(left.x, left.y, left.width, 1),
        );
        let fits = usize::from(left.height.saturating_sub(1) / each);
        for (k, e) in entries.iter().take(fits).enumerate() {
            let y = left.y + 1 + k as u16 * each;
            if let (Some(picker), true) = (&self.picker, rows > 0) {
                let pic = Rect::new(left.x + 1, y, left.width.saturating_sub(2), rows);
                let made = self
                    .thumbs
                    .entry(gallery.thumb(e))
                    .or_insert_with_key(|path| {
                        term::protocol(picker, image::open(path).ok()?, pic)
                    });
                if let Some(p) = made {
                    frame.render_widget(Image::new(p), pic);
                }
            }
            let text = vec![
                Line::from(Span::styled(format!(" {}", e.name()), PLAIN)),
                Line::from(Span::styled(
                    format!(" {} · {}", gallery::ago(e.when, now), e.what),
                    DIM,
                )),
            ];
            frame.render_widget(
                Paragraph::new(text),
                Rect::new(left.x, y + rows, left.width, 2),
            );
        }
    }

    /// The picker to draw images with, when images are on.
    pub fn images(&self) -> Option<&Picker> {
        self.picker.as_ref().filter(|_| self.images)
    }

    /// Header suffix naming the plot mode, when there is a choice.
    pub fn tag(&self) -> &'static str {
        match (&self.picker, self.images_on()) {
            (None, _) => "",
            (Some(_), true) => " · image",
            (Some(_), false) => " · text",
        }
    }

    /// The `i` and `s` help keys.
    pub fn help_keys(&self, keys: &mut Vec<(&'static str, &'static str)>) {
        if self.picker.is_some() {
            keys.push(("i", "image/text"));
        }
        keys.push(("s", "save"));
    }

    /// The footer while the save prompt or its outcome is showing.
    pub fn footer(&self) -> Option<Line<'static>> {
        self.save.footer()
    }

    /// Whether plots draw as images (`i`), where the terminal can.
    pub fn images_on(&self) -> bool {
        self.images && self.picker.is_some()
    }

    /// Draw plots as images or not, as a view kept from before; text
    /// only where the terminal has no images.
    pub fn set_images(&mut self, on: bool) {
        self.images = on;
    }

    /// Whether the save prompt is open, where the apply key submits.
    pub fn prompting(&self) -> bool {
        self.save.active()
    }
}

/// Each label of a drawing as its font size, text and baseline (none
/// when turned), as printed: scaled as the drawing's scale says, the size
/// to two decimals. Sequence-logo letters (scaled to their boxes) are left
/// out.
#[cfg(test)]
pub(crate) fn labels(svg: &str) -> Vec<(String, String, Option<f64>)> {
    let attr = |l: &str, key: &str| -> Option<String> {
        let rest = l.split(&format!(" {key}=\"")).nth(1)?;
        Some(rest.split('"').next()?.to_string())
    };
    let mut scale = 1.0;
    let mut out = Vec::new();
    for l in svg.lines() {
        if let Some(s) = l.strip_prefix(r#"<g transform="scale("#) {
            scale = s
                .split(')')
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1.0);
        } else if l == "</g>" {
            scale = 1.0;
        }
        if !l.starts_with("<text") || l.contains(" scale(") {
            continue;
        }
        let text = l.split('>').nth(1).and_then(|r| r.split('<').next());
        let y = attr(l, "y").and_then(|y| y.parse::<f64>().ok());
        let size = attr(l, "font-size").unwrap_or_default();
        let size = match size.parse::<f64>() {
            Ok(pt) => ((pt * scale * 100.0).round() / 100.0).to_string(),
            Err(_) => size,
        };
        out.push((size, text.unwrap_or("").into(), y.map(|y| y * scale)));
    }
    out
}

/// Each `<text>` of `svg` as the box it takes, left, top, right, bottom,
/// with what it says (unescaped): about 0.55 of its size a character
/// across (a point or comma half that, as [`tick_width`]), 0.75 above its
/// baseline and 0.25 below; a turned one centred on its anchor. Sequence-logo letters (scaled to their boxes) are left out.
#[cfg(test)]
pub(crate) fn text_boxes(svg: &str) -> Vec<((f64, f64, f64, f64), String)> {
    let attr = |l: &str, name: &str| -> Option<String> {
        let key = format!(" {name}=\"");
        let start = l.find(&key)? + key.len();
        Some(l[start..start + l[start..].find('"')?].to_string())
    };
    let num = |l: &str, name: &str| attr(l, name).unwrap().parse::<f64>().unwrap();
    (svg.lines()
        .filter(|l| l.starts_with("<text") && !l.contains(" scale(")))
    .map(|l| {
        let s = l[l.find('>').unwrap() + 1..l.rfind("</text>").unwrap()]
            .replace("&lt;", "<")
            .replace("&gt;", ">")
            .replace("&quot;", "\"")
            .replace("&apos;", "'")
            .replace("&amp;", "&");
        let size = num(l, "font-size");
        let w = tick_width(&s, size);
        let (up, down) = (0.75 * size, 0.25 * size);
        let b = match attr(l, "transform") {
            Some(t) => {
                let inner = &t["translate(".len()..t.find(')').unwrap()];
                let mut it = inner.split([' ', ',']).map(|v| v.parse::<f64>().unwrap());
                let (x, y) = (it.next().unwrap(), it.next().unwrap());
                (x - up, y - w / 2.0, x + down, y + w / 2.0)
            }
            None => {
                let (x, y) = (num(l, "x"), num(l, "y"));
                let left = match attr(l, "text-anchor").as_deref() {
                    Some("middle") => x - w / 2.0,
                    Some("end") => x - w,
                    _ => x,
                };
                (left, y - up, left + w, y + down)
            }
        };
        (b, s)
    })
    .collect()
}

/// A short fingerprint of a drawing (FNV-1a), for tests that pin a
/// figure's output.
#[cfg(test)]
pub(crate) fn fingerprint(svg: &str) -> u64 {
    svg.bytes().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

#[cfg(test)]
#[path = "tests/figure.rs"]
mod tests;
