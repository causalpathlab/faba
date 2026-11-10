//! The gene modules' GO terms, full pane on `a` (annotate; see
//! [`super::super::go`]): read on a thread of their own the first time,
//! then one gene module's terms at a time, ↑/↓ to the next, as the modules
//! list steps them.

use super::super::go::Table;
use super::{ModulesView, Pane};
use data_beans::interactive::ui::{help_line, panel, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::Rect;
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Wrap};
use ratatui::Frame;
use std::sync::{Arc, OnceLock};

/// The terms as they arrive: none until read, or why there are none.
pub(super) type Pending = Arc<OnceLock<Result<Table, String>>>;

impl ModulesView<'_> {
    /// `a`: open the GO pane on the gene module in focus, the terms read
    /// on a thread the first time.
    pub(super) fn open_go(&mut self) {
        let Some(ask) = self.go_ask.clone() else {
            self.status = Some("no GO for this view".into());
            return;
        };
        self.going = true;
        // Read once; a read that failed is tried again (its files fetched
        // since, say).
        if self.go.is_some() || matches!(self.go_table, Some(Ok(_))) {
            return;
        }
        self.go_table = None;
        let pending: Pending = Arc::default();
        self.go = Some(pending.clone());
        let source = ask.source(self.m);
        let (genes, roles) = (self.m.gene_main.clone(), self.m.gene_roles);
        std::thread::spawn(move || {
            let table = source.and_then(|s| {
                super::super::go::terms(&genes, roles, &s, ask.min_weight)
                    .map_err(|e| format!("{e:#}"))
            });
            let _ = pending.set(table);
        });
    }

    /// A key while the GO pane is open: `a` or Esc closes it, ↑/↓ (j/k)
    /// steps the gene module as the modules list does (and with it the
    /// focus there), `q` quits.
    pub(super) fn go_key(&mut self, key: KeyEvent) {
        let by = match key.code {
            KeyCode::Char('a') | KeyCode::Esc => return self.going = false,
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
        let help = help_line(&[("↑/↓", "gene module"), ("a/Esc", "close"), ("q", "quit")]);
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
                "{} genes of weight at least {} ({} in all) of the {} annotated; {}",
                m.n, table.min_weight, m.all, table.universe, table.source
            ),
            DIM,
        ),
        Line::raw(""),
    ];
    if m.hits.is_empty() {
        lines.push(said(m.why_none(), DIM));
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
        // A GMT's name column is its description (a URL, or `-`): its id
        // names it.
        let named = if h.namespace.is_empty() {
            &h.id
        } else {
            &h.name
        };
        let name: String = named.chars().take(46).collect();
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
