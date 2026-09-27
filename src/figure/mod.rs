//! Figures behind the full-screen views: each view draws itself as an SVG
//! here, and [`save`] writes that one drawing as a vector PDF and a PNG, so
//! the exported figure is the same whatever the terminal could show.

use std::fmt::Write as _;
use std::sync::{Arc, OnceLock};

pub mod term;

use data_beans::interactive::ui::{compact, input_line, Scale, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui_image::picker::Picker;

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
const FONT: &str = "Helvetica, Arial, 'DejaVu Sans', sans-serif";

/// PNG pixels per SVG unit (the SVG is laid out in points: 3x is 216 dpi).
const PNG_SCALE: f32 = 3.0;

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

/// An SVG drawing in points, white background.
pub struct Canvas {
    width: f64,
    height: f64,
    body: String,
}

impl Canvas {
    pub fn new(width: f64, height: f64) -> Self {
        let mut c = Self {
            width,
            height,
            body: String::new(),
        };
        c.rect(0.0, 0.0, width, height, "#ffffff");
        c
    }

    /// A filled rectangle; empty ones (zero or negative size) are skipped,
    /// as SVG renderers reject them.
    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        if !(w > 0.0 && h > 0.0) {
            return;
        }
        let _ = writeln!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{:.2}" height="{:.2}" fill="{fill}"/>"#,
            w.max(0.0),
            h.max(0.0)
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
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-family="{FONT}" font-size="{size}" text-anchor="{}" fill="{fill}">{}</text>"#,
            anchor.svg(),
            esc(s)
        );
    }

    /// Text turned to read bottom to top, centred on `(x, y)`.
    pub fn vtext(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str) {
        let _ = writeln!(
            self.body,
            r#"<text transform="translate({x:.2} {y:.2}) rotate(-90)" font-family="{FONT}" font-size="{size}" text-anchor="middle" fill="{fill}">{}</text>"#,
            esc(s)
        );
    }

    pub fn bold(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: Anchor, fill: &str) {
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
            body: String::new(),
        }
    }

    /// The shapes drawn so far, without the `<svg>` wrapper.
    pub fn into_body(self) -> String {
        self.body
    }

    pub fn finish(self) -> String {
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">
{body}</svg>
"#,
            w = self.width,
            h = self.height,
            body = self.body
        )
    }
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

fn unscaled(scale: Scale, t: f64) -> f64 {
    match scale {
        Scale::Log => 10f64.powf(t) - 1.0,
        Scale::Sqrt => t * t,
        Scale::Linear => t,
    }
}

/// A bar panel matching a terminal histogram: one bar per value, an optional
/// front subset over faint full bars, accent bars, x tick labels by bar, a
/// dashed pointer, and site marks under the axis.
pub struct Bars<'a> {
    pub values: &'a [f64],
    pub front: Option<&'a [f64]>,
    pub accent: &'a dyn Fn(usize) -> bool,
    pub y_scale: Scale,
    /// Top of the y axis; `None` scales to the tallest bar.
    pub y_max: Option<f64>,
    pub ticks: Vec<(usize, String)>,
    pub pointer: Option<usize>,
    pub marks: Vec<usize>,
    pub title: String,
    pub x_title: String,
    pub y_title: String,
}

impl Bars<'_> {
    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) {
        let (left, right, top, bottom) = (52.0, 12.0, 20.0, 38.0);
        let (px, py, pw, ph) = (x + left, y + top, w - left - right, h - top - bottom);
        c.bold(x + left, y + 13.0, &self.title, 11.0, Anchor::Start, INK);
        let n = self.values.len().max(1);
        let bw = pw / n as f64;
        let tallest = self
            .values
            .iter()
            .map(|&v| scaled(self.y_scale, v))
            .fold(0.0, f64::max);
        let max = self
            .y_max
            .map_or(tallest, |m| scaled(self.y_scale, m).max(tallest));
        let height = |v: f64| {
            if max <= 0.0 {
                0.0
            } else {
                scaled(self.y_scale, v) / max * ph
            }
        };

        // y gridlines and labels at 0, half and full scale.
        for f in [0.0, 0.5, 1.0] {
            let gy = py + ph - f * ph;
            if f > 0.0 {
                c.line(px, gy, px + pw, gy, FAINT, 0.4);
            }
            let label = if max > 0.0 {
                compact(unscaled(self.y_scale, f * max))
            } else {
                "0".into()
            };
            c.text(px - 5.0, gy + 3.0, &label, 8.0, Anchor::End, MUTED);
        }
        c.vtext(x + 12.0, py + ph / 2.0, &self.y_title, 8.0, MUTED);

        let gap = if bw > 3.0 { 0.15 * bw } else { 0.0 };
        let mut bars = |values: &[f64], colour: &dyn Fn(usize) -> &'static str| {
            for (i, &v) in values.iter().enumerate() {
                let bh = height(v);
                if bh > 0.0 {
                    let bx = px + i as f64 * bw + gap / 2.0;
                    c.rect(bx, py + ph - bh, bw - gap, bh, colour(i));
                }
            }
        };
        let tone = |i: usize| if (self.accent)(i) { ACCENT } else { BAR };
        match self.front {
            Some(front) => {
                bars(self.values, &|_| FAINT);
                bars(front, &tone);
            }
            None => bars(self.values, &tone),
        }

        if let Some(p) = self.pointer {
            let cx = px + (p as f64 + 0.5) * bw;
            c.dashed(cx, py, cx, py + ph, ACCENT, 0.8);
        }
        c.line(px, py + ph, px + pw, py + ph, INK, 0.6);
        c.line(px, py, px, py + ph, INK, 0.6);
        for &m in &self.marks {
            let mx = px + (m as f64 + 0.5) * bw;
            c.line(mx, py + ph, mx, py + ph + 3.0, MUTED, 0.5);
        }
        for (i, label) in &self.ticks {
            let tx = px + (*i as f64 + 0.5) * bw;
            c.line(tx, py + ph, tx, py + ph + 4.0, INK, 0.6);
            c.text(tx, py + ph + 13.0, label, 8.0, Anchor::Middle, INK);
        }
        c.text(
            px + pw / 2.0,
            y + h - 6.0,
            &self.x_title,
            8.0,
            Anchor::Middle,
            MUTED,
        );
    }
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
}

impl Diverging<'_> {
    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) {
        let (up, down) = split_signed(self.values);
        let max = up.iter().chain(&down).fold(0.0, |m: f64, &v| m.max(v));
        Mirror {
            up: Half::plain(&up, ACCENT),
            down: Half::plain(&down, BAR),
            y_scale: Scale::Linear,
            y_max: Some(max),
            y_labels: [(self.label)(max), (self.label)(0.0), (self.label)(-max)],
            ticks: self.ticks.clone(),
            pointer: self.pointer,
            title: self.title.clone(),
            x_title: self.x_title.clone(),
            y_title: self.y_title.clone(),
        }
        .draw(c, x, y, w, h);
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
}

impl Mirror<'_> {
    pub fn draw(&self, c: &mut Canvas, x: f64, y: f64, w: f64, h: f64) {
        let (left, right, top, bottom) = (52.0, 12.0, 20.0, 38.0);
        let (px, py, pw, ph) = (x + left, y + top, w - left - right, h - top - bottom);
        c.bold(x + left, y + 13.0, &self.title, 11.0, Anchor::Start, INK);
        let n = self.up.values.len().max(self.down.values.len()).max(1);
        let bw = pw / n as f64;
        let all = self.up.values.iter().chain(self.down.values);
        let tallest = all.map(|&v| scaled(self.y_scale, v)).fold(0.0, f64::max);
        let max = self
            .y_max
            .map_or(tallest, |m| scaled(self.y_scale, m).max(tallest));
        let zero = py + ph / 2.0;
        let height = |v: f64| {
            if max <= 0.0 {
                0.0
            } else {
                scaled(self.y_scale, v) / max * ph / 2.0
            }
        };
        let own = compact(unscaled(self.y_scale, max));
        for (f, label) in [(1.0, 0), (0.0, 1), (-1.0, 2)] {
            let gy = zero - f * ph / 2.0;
            if f != 0.0 {
                c.line(px, gy, px + pw, gy, FAINT, 0.4);
            }
            let label = match self.y_labels[label].as_str() {
                "" if f == 0.0 => "0",
                "" => &own,
                l => l,
            };
            c.text(px - 5.0, gy + 3.0, label, 8.0, Anchor::End, MUTED);
        }
        c.vtext(x + 12.0, py + ph / 2.0, &self.y_title, 8.0, MUTED);
        let gap = if bw > 3.0 { 0.15 * bw } else { 0.0 };
        for (half, sign) in [(&self.up, -1.0), (&self.down, 1.0)] {
            let mut bars = |values: &[f64], colour: &'static str| {
                for (i, &v) in values.iter().enumerate() {
                    let bh = height(v);
                    let bx = px + i as f64 * bw + gap / 2.0;
                    let top = if sign < 0.0 { zero - bh } else { zero };
                    c.rect(bx, top, bw - gap, bh, colour);
                }
            };
            match half.front {
                Some(front) => {
                    bars(half.values, FAINT);
                    bars(front, ACCENT);
                }
                None => bars(half.values, half.colour),
            }
        }
        c.text(px + 4.0, py + 9.0, &self.up.name, 8.0, Anchor::Start, MUTED);
        c.text(
            px + 4.0,
            py + ph - 3.0,
            &self.down.name,
            8.0,
            Anchor::Start,
            MUTED,
        );
        if let Some(p) = self.pointer {
            let cx = px + (p as f64 + 0.5) * bw;
            c.dashed(cx, py, cx, py + ph, ACCENT, 0.8);
        }
        c.line(px, zero, px + pw, zero, INK, 0.6);
        c.line(px, py, px, py + ph, INK, 0.6);
        let base = py + ph;
        for (i, label) in &self.ticks {
            let tx = px + (*i as f64 + 0.5) * bw;
            c.line(tx, base, tx, base + 4.0, INK, 0.6);
            c.text(tx, base + 13.0, label, 8.0, Anchor::Middle, INK);
        }
        c.text(
            px + pw / 2.0,
            y + h - 6.0,
            &self.x_title,
            8.0,
            Anchor::Middle,
            MUTED,
        );
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

/// SVG parsing options, with the system fonts loaded once.
fn options() -> usvg::Options<'static> {
    usvg::Options {
        fontdb: fontdb(),
        ..usvg::Options::default()
    }
}

/// Write `svg` as `{prefix}.pdf` (vector) and `{prefix}.png`; returns the
/// paths written.
pub fn save(svg: &str, prefix: &str) -> anyhow::Result<Vec<String>> {
    let prefix = strip_extension(prefix);
    let tree = usvg::Tree::from_str(svg, &options())?;
    let pdf = svg2pdf::to_pdf(&tree, Default::default(), Default::default())
        .map_err(|e| anyhow::anyhow!("PDF conversion failed: {e:?}"))?;
    let pdf_path = format!("{prefix}.pdf");
    std::fs::write(&pdf_path, pdf)?;
    let png_path = format!("{prefix}.png");
    term::pixmap(&tree, PNG_SCALE)?.save_png(&png_path)?;
    Ok(vec![pdf_path, png_path])
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
        match key.code {
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

/// The save-as prompt every view shares: it opens with a suggested name,
/// Enter confirms, Esc closes; the outcome shows until the next key.
pub struct SavePrompt {
    input: LineInput,
    status: Option<Result<String, String>>,
    default: String,
}

impl SavePrompt {
    pub fn new(default: &str) -> Self {
        Self {
            input: LineInput::new(512),
            status: None,
            default: default.to_string(),
        }
    }

    pub fn active(&self) -> bool {
        self.input.active()
    }

    pub fn open(&mut self) {
        self.input.open(&self.default);
        self.status = None;
    }

    /// Handle a key while open: `Some(prefix)` once the user confirms.
    pub fn handle(&mut self, key: KeyEvent) -> Option<String> {
        match self.input.handle(key) {
            Edit::Submitted(name) if !name.is_empty() => {
                self.default = name.clone();
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

    /// The footer while the prompt or an outcome is showing.
    pub fn footer(&self) -> Option<Line<'static>> {
        if let Some(buf) = self.input.text() {
            return Some(input_line(
                "save as (.pdf + .png): ",
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

/// What every view's figure needs: the save prompt (`s`), and drawing plots
/// as images where the terminal can (`i` switches to text and back).
pub struct Controls {
    save: SavePrompt,
    picker: Option<Picker>,
    images: bool,
}

impl Controls {
    /// Text plots only; [`Controls::detect`] asks the terminal for images.
    pub fn new(default_name: &str) -> Self {
        Self {
            save: SavePrompt::new(default_name),
            picker: None,
            images: false,
        }
    }

    /// Use images if the terminal has them. Call before the view opens.
    pub fn detect(mut self) -> Self {
        self.picker = term::picker();
        self.images = self.picker.is_some();
        self
    }

    #[cfg(test)]
    pub fn with_picker(mut self, picker: Picker) -> Self {
        self.picker = Some(picker);
        self.images = true;
        self
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
            KeyCode::Char('s') => self.save.open(),
            KeyCode::Char('i') if self.picker.is_some() => self.images ^= true,
            _ => return Key::Pass,
        }
        Key::Used
    }

    /// Save `svg` under the confirmed `prefix` and show how it went.
    pub fn save(&mut self, svg: &str, prefix: &str) {
        self.save.report(save(svg, prefix));
    }

    /// The picker to draw images with, when images are on.
    pub fn images(&self) -> Option<&Picker> {
        self.picker.as_ref().filter(|_| self.images)
    }

    /// Header suffix naming the plot mode, when there is a choice.
    pub fn tag(&self) -> &'static str {
        match (&self.picker, self.images) {
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
}

#[cfg(test)]
#[path = "tests/figure.rs"]
mod tests;
