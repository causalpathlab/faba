//! Figures behind the full-screen views: each view draws itself as an SVG
//! here, and [`save`] writes that one drawing as a vector PDF and a PNG, so
//! the exported figure is the same whatever the terminal could show.

use std::fmt::Write as _;
use std::sync::{Arc, OnceLock};

use data_beans::interactive::ui::{compact, input_line, Scale, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};

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

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
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

    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
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

/// `v` on a view's y scale, as the terminal histogram draws it.
fn scaled(scale: Scale, v: f64) -> f64 {
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
        let max = self
            .values
            .iter()
            .map(|&v| scaled(self.y_scale, v))
            .fold(0.0, f64::max);
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
pub fn strip_extension(prefix: &str) -> &str {
    prefix
        .strip_suffix(".pdf")
        .or_else(|| prefix.strip_suffix(".png"))
        .unwrap_or(prefix)
}

/// Write `svg` as `{prefix}.pdf` (vector) and `{prefix}.png`; returns the
/// paths written.
pub fn save(svg: &str, prefix: &str) -> anyhow::Result<Vec<String>> {
    let prefix = strip_extension(prefix);
    let opt = usvg::Options {
        fontdb: fontdb(),
        ..usvg::Options::default()
    };
    let tree = usvg::Tree::from_str(svg, &opt)?;

    let pdf = svg2pdf::to_pdf(&tree, Default::default(), Default::default())
        .map_err(|e| anyhow::anyhow!("PDF conversion failed: {e:?}"))?;
    let pdf_path = format!("{prefix}.pdf");
    std::fs::write(&pdf_path, pdf)?;

    let size = tree
        .size()
        .to_int_size()
        .scale_by(PNG_SCALE)
        .ok_or_else(|| anyhow::anyhow!("figure has no size"))?;
    let mut pixmap = resvg::tiny_skia::Pixmap::new(size.width(), size.height())
        .ok_or_else(|| anyhow::anyhow!("figure too large to rasterise"))?;
    pixmap.fill(resvg::tiny_skia::Color::WHITE);
    resvg::render(
        &tree,
        resvg::tiny_skia::Transform::from_scale(PNG_SCALE, PNG_SCALE),
        &mut pixmap.as_mut(),
    );
    let png_path = format!("{prefix}.png");
    pixmap.save_png(&png_path)?;
    Ok(vec![pdf_path, png_path])
}

/// The save-as prompt every view shares: `s` opens it with a suggested
/// name, Enter confirms, Esc closes; the outcome shows until the next key.
pub struct SavePrompt {
    input: Option<String>,
    status: Option<Result<String, String>>,
    default: String,
}

impl SavePrompt {
    pub fn new(default: &str) -> Self {
        Self {
            input: None,
            status: None,
            default: default.to_string(),
        }
    }

    pub fn active(&self) -> bool {
        self.input.is_some()
    }

    pub fn open(&mut self) {
        self.input = Some(self.default.clone());
        self.status = None;
    }

    /// Handle a key while open: `Some(prefix)` once the user confirms.
    pub fn handle(&mut self, key: KeyEvent) -> Option<String> {
        let buf = self.input.as_mut()?;
        match key.code {
            KeyCode::Char(ch) if buf.len() < 512 => buf.push(ch),
            KeyCode::Backspace => {
                buf.pop();
            }
            KeyCode::Esc => self.input = None,
            KeyCode::Enter => {
                let name = self.input.take().unwrap_or_default();
                let name = name.trim();
                if !name.is_empty() {
                    self.default = name.to_string();
                    return Some(name.to_string());
                }
            }
            _ => {}
        }
        None
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
        if let Some(buf) = &self.input {
            return Some(input_line(
                "save as (.pdf + .png): ",
                buf,
                &[("Enter", "save"), ("Esc", "back")],
            ));
        }
        self.status.as_ref().map(|s| match s {
            Ok(msg) => Line::from(Span::styled(format!(" {msg}"), DIM)),
            Err(msg) => Line::from(Span::styled(format!(" {msg}"), HIGHLIGHT)),
        })
    }
}

#[cfg(test)]
#[path = "tests/figure.rs"]
mod tests;
