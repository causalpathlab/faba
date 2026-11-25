//! The gene modules' GO terms (see [`super::super::go`]), read on a thread
//! of their own the first time they are asked for, never before: a column
//! right of the modules grid on `a` (annotate), each gene module's
//! strongest terms on its row, `a` again to hide or show it; and full pane
//! on `G`, one gene module's terms at a time, ↑/↓ to the next, as the
//! modules list steps them.

use super::super::go::{Hit, Table, DEFAULT_Q};
use super::grid::{COMPACT, COUNTS, MODULE_EXTRA, MODULE_LABEL, SITE_LABEL};
use super::{ModulesView, Pane};
use crate::tui::ui::{help_line, panel, DIM, HIGHLIGHT, PLAIN};
use ratatui::buffer::Buffer;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::Rect;
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Wrap};
use ratatui::Frame;
use std::sync::{Arc, OnceLock};

/// The terms as they arrive: none until read, or why there are none.
pub(super) type Pending = Arc<OnceLock<Result<Table, String>>>;

/// The GO column's width at most, and the least it is shown at.
const GO_W: u16 = 44;
const GO_LEAST: u16 = 14;

/// What the footer says when the GO column has no room.
const NO_ROOM: &str = "no room for the GO column: a wider screen shows it; G opens GO";

/// The least the modules and sites panes need beside the GO column: their
/// labels and numbers, one compact cell module each, and their frames.
const PANES_LEAST: u16 =
    MODULE_LABEL + COMPACT + MODULE_EXTRA + 2 + SITE_LABEL + COMPACT + COUNTS + 4;

impl ModulesView<'_> {
    /// `a`: show the GO column, the terms read on a thread the first time;
    /// again, hide or show it, nothing read anew. After a failed read, it
    /// reads again with the column on, to show what it reads; with no room
    /// for the column, it still turns it on and off, saying which.
    pub(super) fn toggle_go_column(&mut self) {
        let retry = matches!(self.go_table, Some(Err(_))) && self.go.is_none();
        if !self.read_go() {
            return;
        }
        self.go_column = retry || !self.go_column;
        // With no room the column cannot show either way: say which it is.
        if self.go_no_room {
            let said = if self.go_column {
                NO_ROOM
            } else {
                "GO column off"
            };
            self.status = Some(said.into());
        }
    }

    /// `G`: open the GO pane on the gene module in focus, the terms read
    /// on a thread the first time.
    pub(super) fn open_go(&mut self) {
        if self.read_go() {
            self.full = Some(super::Full::Go);
        }
    }

    /// Start reading the terms unless they are read or being read: whether
    /// this view has GO to read (if not, the footer says so).
    fn read_go(&mut self) -> bool {
        let Some(ask) = self.go_ask.clone() else {
            self.status = Some("no GO for this view".into());
            return false;
        };
        // Read once; a read that failed is tried again (its files fetched
        // since, say).
        if self.go.is_some() || matches!(self.go_table, Some(Ok(_))) {
            return true;
        }
        self.go_table = None;
        let pending: Pending = Arc::default();
        self.go = Some(pending.clone());
        let source = ask.source(self.m);
        let (genes, roles) = (self.m.gene_main.clone(), self.m.gene_roles);
        std::thread::spawn(move || {
            let table = source.and_then(|s| {
                super::super::go::terms(&genes, roles, &s, ask.min_weight, ask.q)
                    .map_err(|e| format!("{e:#}"))
            });
            let _ = pending.set(table);
        });
        true
    }

    /// The GO column's width in a screen `width` wide: none when it is not
    /// asked for or there is no room for it, which the footer says once.
    pub(super) fn fit_go_column(&mut self, width: u16) -> u16 {
        let room = width.saturating_sub(PANES_LEAST).min(GO_W);
        let fits = self.go_column && room >= GO_LEAST;
        if self.go_column && !fits && !self.go_no_room {
            self.status = Some(NO_ROOM.into());
        }
        self.go_no_room = self.go_column && !fits;
        self.go_w = if fits { room } else { 0 };
        self.go_w
    }

    /// The GO column in `area`, its head on the grid's head line, then
    /// gene modules `roles` on their rows, two lines each, as the grid
    /// draws them.
    pub(super) fn draw_go_column(&self, buf: &mut Buffer, area: Rect, roles: &[usize]) {
        // A column's gap from the grid.
        let (x, w) = (area.x + 1, area.width.saturating_sub(1) as usize);
        // The cutoff asked for (`--go-q`), the one the terms were held to.
        let q = self.go_ask.as_ref().map_or(DEFAULT_Q, |a| a.q);
        buf.set_stringn(x, area.y, format!("GO (q <= {q})"), w, DIM);
        let row_y = |k: usize| area.y + 1 + 2 * k as u16;
        let table = match &self.go_table {
            Some(Ok(table)) => table,
            Some(Err(why)) => {
                buf.set_stringn(x, row_y(0), why, w, HIGHLIGHT);
                return;
            }
            None => {
                buf.set_stringn(x, row_y(0), "reading GO...", w, DIM);
                return;
            }
        };
        let fit = (area.height.saturating_sub(1) / 2) as usize;
        for (k, &g) in roles.iter().take(fit).enumerate() {
            let hits = table.of(g).map_or(&[][..], |m| &m.hits);
            if hits.is_empty() {
                buf.set_stringn(x, row_y(k), format!("none at q <= {q}"), w, DIM);
                continue;
            }
            for (dy, line) in packed(hits, w).into_iter().enumerate() {
                let mut at = x;
                for (text, style) in line {
                    at = buf.set_stringn(at, row_y(k) + dy as u16, text, w, style).0;
                }
            }
        }
    }

    /// A key while the GO pane is open: `G` or Esc closes it, ↑/↓ (j/k)
    /// steps the gene module as the modules list does (and with it the
    /// focus there), `q` quits.
    pub(super) fn go_key(&mut self, key: KeyEvent) {
        let by = match key.code {
            KeyCode::Char('G') | KeyCode::Esc => return self.full = None,
            KeyCode::Char('q') => return self.done = true,
            KeyCode::Down | KeyCode::Char('j') => 1,
            KeyCode::Up | KeyCode::Char('k') => -1,
            _ => return,
        };
        let pane = std::mem::replace(&mut self.pane, Pane::Modules);
        self.step(by);
        self.pane = pane;
    }

    /// The terms, once read: kept, and whether they arrived just now.
    pub(super) fn take_go(&mut self) -> bool {
        let Some(table) = self.go.as_ref().and_then(|p| p.get()) else {
            return false;
        };
        self.go_table = Some(table.clone());
        self.go = None;
        true
    }

    /// The GO pane in `area` above the footer in `foot`: the gene module's
    /// terms, a line each, with what they were tested against.
    pub(super) fn draw_go(&self, frame: &mut Frame, area: Rect, foot: Rect) {
        let help = help_line(&self.footer_pairs().unwrap_or_default());
        frame.render_widget(Paragraph::new(help), foot);
        let g = self.module;
        let block = panel(format!(" GO · gene module M{} ", g + 1), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let lines = match &self.go_table {
            None => vec![said("reading the GO annotations…".into(), DIM)],
            Some(Err(why)) => vec![said(why.clone(), HIGHLIGHT)],
            Some(Ok(table)) => go_lines(table, g, inner.width),
        };
        frame.render_widget(Paragraph::new(lines).wrap(Wrap { trim: false }), inner);
    }
}

/// `s` as a line, indented one column.
fn said(s: String, style: Style) -> Line<'static> {
    Line::from(Span::styled(format!(" {s}"), style))
}

/// Gene module `g`'s lines of `table`, `width` columns wide.
fn go_lines(table: &Table, g: usize, width: u16) -> Vec<Line<'static>> {
    let Some(m) = table.of(g) else {
        return vec![said("no such gene module".into(), DIM)];
    };
    let mut lines = vec![
        said(
            format!(
                "{} genes of weight at least {} ({} in all) of the {} annotated; terms at q <= {}; {}",
                m.n, table.min_weight, m.all, table.universe, table.q, table.source
            ),
            DIM,
        ),
        Line::raw(""),
    ];
    if m.hits.is_empty() {
        lines.push(said(m.why_none(table.q), DIM));
        return lines;
    }
    lines.push(said(
        format!(
            "{:>2}  {:<2} {:<46} {:>9} {:>6} {:>8}  genes",
            "#", "ns", "term", "k/K", "odds", "q"
        ),
        DIM,
    ));
    let room = (width as usize).saturating_sub(4);
    for (i, h) in m.hits.iter().enumerate() {
        let name: String = h.named().chars().take(46).collect();
        let text = format!(
            "{:>2}  {:<2} {:<46} {:>9} {:>6.2} {:>8.1e}  {}",
            i + 1,
            h.namespace,
            name,
            format!("{}/{}", h.k, h.big_k),
            h.odds,
            h.q,
            h.genes.join(",")
        );
        lines.push(said(text.chars().take(room).collect(), Style::default()));
    }
    lines.push(Line::raw(""));
    lines.push(said(
        "k: the module's genes in the term; K: the term's annotated genes".into(),
        DIM,
    ));
    lines
}

/// Text and its style, a piece of a line.
type Piece = (String, Style);

/// A term in the GO column, at most `w` wide: its name, shortened to fit,
/// and its q.
fn term(h: &Hit, w: usize) -> Vec<Piece> {
    let q = format!(" {}", super::super::go::q_said(h.q));
    let name = h.named();
    let room = w.saturating_sub(q.chars().count());
    let name = match name.chars().count() > room {
        true => (name.chars().take(room.saturating_sub(1)))
            .chain(['…'])
            .collect(),
        false => name.to_string(),
    };
    vec![(name, PLAIN), (q, DIM)]
}

/// Terms `hits` (strongest first) on two lines `w` wide: as many as fit,
/// then how many more there are.
fn packed(hits: &[Hit], w: usize) -> Vec<Vec<Piece>> {
    let terms: Vec<Vec<Piece>> = hits.iter().map(|h| term(h, w)).collect();
    (0..=terms.len())
        .rev()
        .find_map(|n| {
            let more = terms.len() - n;
            let tail = (more > 0).then(|| vec![(format!("+{more} more"), DIM)]);
            laid(terms[..n].iter().cloned().chain(tail), w, 2)
        })
        .unwrap_or_default()
}

/// `items` laid out in order on at most `lines` lines `w` wide, two
/// spaces apart; none when they do not fit.
fn laid(
    items: impl Iterator<Item = Vec<Piece>>,
    w: usize,
    lines: usize,
) -> Option<Vec<Vec<Piece>>> {
    let width = |p: &[Piece]| p.iter().map(|(t, _)| t.chars().count()).sum::<usize>();
    let mut out: Vec<Vec<Piece>> = vec![Vec::new()];
    for item in items {
        let (used, len) = (width(out.last().expect("a line")), width(&item));
        if len > w {
            return None;
        }
        let line = out.last_mut().expect("a line");
        if used == 0 {
            line.extend(item);
        } else if used + 2 + len <= w {
            line.push(("  ".into(), PLAIN));
            line.extend(item);
        } else if out.len() < lines {
            out.push(item);
        } else {
            return None;
        }
    }
    Some(out)
}
