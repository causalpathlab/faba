//! Full-screen sequence logo for `faba pwm --interactive`.
//!
//! One stack of letters per position around the sites, tallest on top, in
//! information content (bits) or frequency, in the usual logo colours; the
//! site itself (position 0) is marked in the accent colour. A cursor reads
//! out a position's counts.

use data_beans::interactive::ui::{header, help_line, panel, Screen, DIM, HIGHLIGHT, PLAIN};
use ratatui::buffer::Buffer;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::data::dna::DnaBaseCount;
use crate::figure::term::PlotImage;
use crate::figure::{self, Anchor, Canvas, Controls, Key, ACCENT, INK, MUTED};

/// Logo order: alphabetical, as sequence logos are usually keyed.
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// The usual sequence-logo colours, in [`BASES`] order: A green, C blue,
/// G amber, T red.
const BASE_RGB: [(u8, u8, u8); 4] = [
    (0x10, 0x96, 0x48),
    (0x25, 0x5c, 0x99),
    (0xf7, 0xb3, 0x2b),
    (0xd6, 0x28, 0x39),
];

fn base_hex(b: usize) -> String {
    let (r, g, bl) = BASE_RGB[b];
    format!("#{r:02x}{g:02x}{bl:02x}")
}

fn base_style(b: usize) -> Style {
    let (r, g, bl) = BASE_RGB[b];
    Style::new()
        .fg(Color::Rgb(r, g, bl))
        .add_modifier(Modifier::BOLD)
}

use crate::figure::GUTTER;

/// Widest column a position gets.
const MAX_COL: u16 = 5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Bits,
    Freq,
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Bits => "bits",
            Mode::Freq => "frequency",
        }
    }

    /// Height of the full y axis.
    fn top(self) -> f64 {
        match self {
            Mode::Bits => 2.0,
            Mode::Freq => 1.0,
        }
    }
}

/// One position's base counts, in [`BASES`] order.
struct Column {
    counts: [u64; 4],
}

impl Column {
    fn total(&self) -> u64 {
        self.counts.iter().sum()
    }

    fn freqs(&self) -> [f64; 4] {
        let n = self.total().max(1) as f64;
        self.counts.map(|c| c as f64 / n)
    }

    /// Information content in bits: 2 minus the Shannon entropy.
    fn bits(&self) -> f64 {
        if self.total() == 0 {
            return 0.0;
        }
        let h: f64 = self
            .freqs()
            .iter()
            .filter(|&&p| p > 0.0)
            .map(|&p| -p * p.log2())
            .sum();
        (2.0 - h).max(0.0)
    }

    /// Letter heights on the `mode` axis, in [`BASES`] order.
    fn heights(&self, mode: Mode) -> [f64; 4] {
        let scale = match mode {
            Mode::Bits => self.bits(),
            Mode::Freq => 1.0,
        };
        self.freqs().map(|p| p * scale)
    }
}

/// 5x7 bitmaps of the bases, in [`BASES`] order, for letters big enough to
/// draw in blocks.
const FONT: [[&str; 7]; 4] = [
    [
        ".###.", "#...#", "#...#", "#####", "#...#", "#...#", "#...#",
    ],
    [
        ".###.", "#...#", "#....", "#....", "#....", "#...#", ".###.",
    ],
    [
        ".###.", "#...#", "#....", "#.###", "#...#", "#...#", ".###.",
    ],
    [
        "#####", "..#..", "..#..", "..#..", "..#..", "..#..", "..#..",
    ],
];

/// Whether base `b`'s bitmap is inked at `row` (of 7) and `col` (of 5).
fn bitmap(b: usize, row: usize, col: usize) -> bool {
    FONT[b][row.min(6)].as_bytes()[col.min(4)] == b'#'
}

/// Bases from the shortest letter to the tallest: bottom to top in a stack.
fn stack_order(heights: &[f64; 4]) -> [usize; 4] {
    let mut order = [0, 1, 2, 3];
    order.sort_by(|&a, &b| heights[a].total_cmp(&heights[b]));
    order
}

/// Split `rows` whole rows between letters in proportion to `heights`
/// (largest remainder), so a stack's rows add up to its rounded total.
fn allot(heights: &[f64; 4], rows: f64) -> [usize; 4] {
    let exact = heights.map(|h| h * rows);
    let total = exact.iter().sum::<f64>().round() as usize;
    let mut out = exact.map(|x| x.floor() as usize);
    let mut order: Vec<usize> = (0..4).collect();
    order.sort_by(|&a, &b| {
        let ra = exact[a] - exact[a].floor();
        let rb = exact[b] - exact[b].floor();
        rb.total_cmp(&ra)
    });
    let short = total.saturating_sub(out.iter().sum());
    for &i in order.iter().take(short) {
        out[i] += 1;
    }
    out
}

/// State of the logo view, independent of the terminal so it can be tested.
pub struct PwmView {
    title: String,
    /// Relative position of the first column.
    first: i64,
    columns: Vec<Column>,
    mode: Mode,
    cursor: usize,
    /// First column on screen, when they do not all fit.
    offset: usize,
    controls: Controls,
    plot: PlotImage,
    done: bool,
}

impl PwmView {
    pub fn new(title: &str, pwm: &[DnaBaseCount], window: i64) -> Self {
        let columns = pwm
            .iter()
            .map(|c| Column {
                counts: [c.count_a(), c.count_c(), c.count_g(), c.count_t()].map(|x| x as u64),
            })
            .collect();
        Self {
            title: title.to_string(),
            first: -window,
            columns,
            mode: Mode::Bits,
            cursor: window.max(0) as usize,
            offset: 0,
            controls: Controls::new("pwm_logo"),
            plot: PlotImage::default(),
            done: false,
        }
    }

    fn rel(&self, j: usize) -> i64 {
        self.first + j as i64
    }

    fn move_cursor(&mut self, delta: isize) {
        let last = self.columns.len().saturating_sub(1) as isize;
        self.cursor = (self.cursor as isize + delta).clamp(0, last) as usize;
    }

    fn stats_line(&self) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let Some(col) = self.columns.get(self.cursor) else {
            return Line::from(dim(" no positions".into()));
        };
        let mut spans = vec![
            dim(" position ".into()),
            Span::styled(format!("{:+}", self.rel(self.cursor)), HIGHLIGHT),
            dim("   reads ".into()),
            Span::raw(col.total().to_string()),
            dim("   ".into()),
        ];
        for (b, p) in BASES.iter().zip(col.freqs()) {
            spans.push(Span::styled(
                format!("{b} "),
                PLAIN.add_modifier(Modifier::BOLD),
            ));
            spans.push(Span::raw(format!("{p:.3}  ")));
        }
        spans.push(dim(" IC ".into()));
        spans.push(Span::raw(format!("{:.3} bits", col.bits())));
        Line::from(spans)
    }

    /// The logo as a figure: letters stretched to their heights.
    fn figure(&self) -> String {
        let n = self.columns.len().max(1) as f64;
        let width = 56.0 + (560.0 / n).clamp(8.0, 40.0) * n + 20.0;
        let mut c = Canvas::new(width, 330.0);
        c.bold(16.0, 22.0, &self.title, 12.0, Anchor::Start, INK);
        self.draw_logo(&mut c, (0.0, 30.0, width, 300.0), None);
        c.finish()
    }

    /// The logo in the box `(x, y, w, h)`, letters stretched to their
    /// heights; `cursor` marks a position under the axis.
    fn draw_logo(
        &self,
        c: &mut Canvas,
        (x0, y0, w, h): (f64, f64, f64, f64),
        cursor: Option<usize>,
    ) {
        let n = self.columns.len().max(1);
        let (left, top, bottom) = (x0 + 56.0, y0 + 10.0, 50.0);
        let cw = ((w - 76.0) / n as f64).max(1.0);
        let ph = (h - 10.0 - bottom).max(1.0);
        let top_value = self.mode.top();
        for f in [0.0, 0.5, 1.0] {
            let y = top + ph - f * ph;
            c.line(left - 4.0, y, left, y, INK, 0.6);
            let label = format!("{}", f * top_value);
            c.text(left - 7.0, y + 3.0, &label, 8.0, Anchor::End, MUTED);
        }
        c.line(left, top, left, top + ph, INK, 0.6);
        c.line(left, top + ph, left + cw * n as f64, top + ph, INK, 0.6);
        c.vtext(x0 + 18.0, top + ph / 2.0, self.mode.name(), 9.0, MUTED);
        for (j, col) in self.columns.iter().enumerate() {
            let x = left + j as f64 * cw;
            let heights = col.heights(self.mode);
            let order = stack_order(&heights);
            let mut y = top + ph;
            for &b in &order {
                let lh = heights[b] / top_value * ph;
                y -= lh;
                c.glyph(x + 0.5, y, cw - 1.0, lh, BASES[b], &base_hex(b));
            }
            let rel = self.rel(j);
            if rel == 0 {
                // The site: an accent tick under the axis.
                c.rect(x + 0.5, top + ph + 1.0, cw - 1.0, 2.5, ACCENT);
            }
            if rel % 5 == 0 || n <= 25 {
                let colour = if rel == 0 { ACCENT } else { INK };
                let label = rel.to_string();
                c.text(
                    x + cw / 2.0,
                    top + ph + 13.0,
                    &label,
                    8.0,
                    Anchor::Middle,
                    colour,
                );
            }
            if cursor == Some(j) {
                c.rect(x + 0.5, top + ph + 17.0, cw - 1.0, 3.0, INK);
            }
        }
        let mid = left + cw * n as f64 / 2.0;
        c.text(
            mid,
            top + ph + 34.0,
            "position relative to the site",
            8.0,
            Anchor::Middle,
            MUTED,
        );
    }

    /// Draw the stacks into `area`: letters over all rows but the last two
    /// (axis and labels), y gutter on the left.
    fn render_logo(&mut self, buf: &mut Buffer, area: Rect) {
        let [plot, axis, labels] = Layout::vertical([
            Constraint::Min(1),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(area);
        let [gutter, chart] =
            Layout::horizontal([Constraint::Length(GUTTER), Constraint::Min(1)]).areas(plot);
        let n = self.columns.len();
        if chart.width == 0 || chart.height == 0 || n == 0 {
            return;
        }
        let cw = (chart.width / n as u16).clamp(1, MAX_COL);
        let fit = (chart.width / cw) as usize;
        // Keep the cursor on screen.
        if self.cursor < self.offset {
            self.offset = self.cursor;
        } else if self.cursor >= self.offset + fit {
            self.offset = self.cursor + 1 - fit;
        }

        let rows = chart.height as f64 / self.mode.top();
        for (slot, j) in (self.offset..n.min(self.offset + fit)).enumerate() {
            let x0 = chart.x + slot as u16 * cw;
            let heights = self.columns[j].heights(self.mode);
            let alloted = allot(&heights, rows);
            let order = stack_order(&heights);
            // Smallest at the bottom, tallest on top. A letter with room
            // for it is drawn from its bitmap; a thin one repeats its glyph.
            let mut bottom = chart.bottom();
            for &b in &order {
                let rows = (alloted[b] as u16).min(bottom - chart.top());
                let big = cw >= 3 && rows >= 3;
                for k in 0..rows {
                    let y = bottom - rows + k;
                    for dx in 0..cw.min(chart.right().saturating_sub(x0)) {
                        let symbol = if !big {
                            BASES[b].to_string()
                        } else if bitmap(
                            b,
                            k as usize * 7 / rows as usize,
                            dx as usize * 5 / cw as usize,
                        ) {
                            "█".to_string()
                        } else {
                            continue;
                        };
                        buf[(x0 + dx, y)]
                            .set_symbol(&symbol)
                            .set_style(Style::reset().patch(base_style(b)));
                    }
                }
                bottom -= rows;
            }
        }

        // y axis: the top of the scale and its middle.
        let gx = gutter.right() - 1;
        for y in gutter.top()..gutter.bottom() {
            buf[(gx, y)].set_symbol("│").set_style(DIM);
        }
        let top = self.mode.top();
        let ylabels = [
            (gutter.top(), top),
            (gutter.top() + gutter.height / 2, top / 2.0),
        ];
        for (y, v) in ylabels {
            let s = format!("{v}");
            buf.set_string(gx.saturating_sub(1 + s.len() as u16), y, &s, DIM);
            buf[(gx, y)].set_symbol("┤").set_style(DIM);
        }

        // x axis: baseline, then every fifth position and the site labelled.
        for x in axis.left()..axis.right() {
            let sym = if x < gx {
                " "
            } else if x == gx {
                "└"
            } else {
                "─"
            };
            buf[(x, axis.y)].set_symbol(sym).set_style(DIM);
        }
        let mut next_free = labels.x;
        for (slot, j) in (self.offset..n.min(self.offset + fit)).enumerate() {
            let x = chart.x + slot as u16 * cw + (cw - 1) / 2;
            let rel = self.rel(j);
            if j == self.cursor {
                buf[(x, axis.y)].set_symbol("▲").set_style(HIGHLIGHT);
            }
            if rel % 5 == 0 {
                let s = rel.to_string();
                if x >= next_free && x + (s.len() as u16) <= labels.right() {
                    let style = if rel == 0 { HIGHLIGHT } else { DIM };
                    buf.set_string(x, labels.y, &s, style);
                    next_free = x + s.len() as u16 + 1;
                }
            }
        }
    }
}

impl crate::tui::View for PwmView {}

impl Screen for PwmView {
    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn handle_key(&mut self, key: KeyEvent) {
        match self.controls.key(key) {
            Key::Pass => self.plot.invalidate(),
            Key::Used => return,
            Key::Save(prefix) => return self.controls.save(&self.figure(), &prefix),
        }
        match key.code {
            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1),
            KeyCode::Home => self.cursor = 0,
            KeyCode::End => self.cursor = self.columns.len().saturating_sub(1),
            KeyCode::Char('0') => self.cursor = (-self.first).max(0) as usize,
            KeyCode::Char('m') => {
                self.mode = match self.mode {
                    Mode::Bits => Mode::Freq,
                    Mode::Freq => Mode::Bits,
                }
            }
            KeyCode::Char('q') | KeyCode::Esc | KeyCode::Enter => self.done = true,
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let [top, body, stats, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(frame.area());
        let n = self.columns.len();
        frame.render_widget(
            header(
                "pwm",
                &self.title,
                &format!(
                    "{} positions · {}{}",
                    n,
                    self.mode.name(),
                    self.controls.tag()
                ),
            ),
            top,
        );
        let block = panel(format!(" {} ", self.mode.name()), true);
        let inner = block.inner(body);
        frame.render_widget(block, body);
        let mut plot = std::mem::take(&mut self.plot);
        let drawn = self.controls.images().is_some_and(|picker| {
            plot.render(frame, inner, picker, |w, h| {
                figure::svg(w, h, |c| {
                    self.draw_logo(c, (0.0, 0.0, w, h), Some(self.cursor))
                })
            })
        });
        self.plot = plot;
        if !drawn {
            self.render_logo(frame.buffer_mut(), inner);
        }
        frame.render_widget(Paragraph::new(self.stats_line()), stats);
        let help = self.controls.footer().unwrap_or_else(|| {
            let mut keys = vec![("←/→", "position"), ("0", "site"), ("m", "bits/frequency")];
            self.controls.help_keys(&mut keys);
            keys.push(("q", "quit"));
            help_line(&keys)
        });
        frame.render_widget(help, footer);
        self.controls.render_save(frame, frame.area());
    }
}

/// Show the logo full screen until the user quits.
pub fn show_pwm(title: &str, pwm: &[DnaBaseCount], window: i64) -> anyhow::Result<()> {
    let mut view = PwmView::new(title, pwm, window);
    view.controls = Controls::new("pwm_logo").detect();
    crate::tui::run_view(&mut view)
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
