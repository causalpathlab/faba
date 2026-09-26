//! Full-screen sequence logo for `faba pwm --interactive`.
//!
//! One stack of letters per position around the sites, tallest on top, in
//! information content (bits) or frequency. Letters are told apart by glyph,
//! not colour; the site itself (position 0) is in the accent colour. A cursor
//! reads out a position's counts.

use data_beans::interactive::ui::{
    header, help_line, panel, run_screen, Screen, ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use ratatui::buffer::Buffer;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::data::dna::DnaBaseCount;

/// Logo order: alphabetical, as sequence logos are usually keyed.
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Width of the y-axis gutter.
const GUTTER: u16 = 6;

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
            let mut order: Vec<usize> = (0..4).collect();
            order.sort_by(|&a, &b| heights[a].total_cmp(&heights[b]));
            let base = if self.rel(j) == 0 { ACCENTED } else { PLAIN };
            let style = base.add_modifier(Modifier::BOLD);
            // Smallest at the bottom, tallest on top.
            let mut y = chart.bottom();
            for &b in &order {
                for _ in 0..alloted[b] {
                    if y == chart.top() {
                        break;
                    }
                    y -= 1;
                    for x in x0..(x0 + cw).min(chart.right()) {
                        buf[(x, y)]
                            .set_symbol(&BASES[b].to_string())
                            .set_style(Style::reset().patch(style));
                    }
                }
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

impl Screen for PwmView {
    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn handle_key(&mut self, key: KeyEvent) {
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
                &format!("{} positions · {}", n, self.mode.name()),
            ),
            top,
        );
        let block = panel(format!(" {} ", self.mode.name()), true);
        let inner = block.inner(body);
        frame.render_widget(block, body);
        self.render_logo(frame.buffer_mut(), inner);
        frame.render_widget(Paragraph::new(self.stats_line()), stats);
        frame.render_widget(
            help_line(&[
                ("←/→", "position"),
                ("0", "site"),
                ("m", "bits/frequency"),
                ("q", "quit"),
            ]),
            footer,
        );
    }
}

/// Show the logo full screen until the user quits.
pub fn show_pwm(title: &str, pwm: &[DnaBaseCount], window: i64) -> anyhow::Result<()> {
    run_screen(&mut PwmView::new(title, pwm, window))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
