//! How the cell modules spread over the batches, full pane on `b` (see
//! [`crate::cluster::batches`]), and the filter on `B` that has the panes
//! show every cell module, the shared ones alone (not a batch effect), or
//! the batch-specific ones alone.

use super::grid::gray_cell;
use super::{ModulesView, Sort};
use crate::cluster::batches::{Label, CONDITION, EVEN, FEWEST};
use data_beans::interactive::ui::{help_line, panel, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::Rect;
use ratatui::layout::{Constraint, Layout};
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Wrap};
use ratatui::Frame;

/// Which cell modules the panes show.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    All,
    /// Those not a batch effect: shared, or a condition's alone.
    Shared,
    /// The batch-specific ones.
    Batch,
}

impl Filter {
    fn next(self) -> Self {
        match self {
            Filter::All => Filter::Shared,
            Filter::Shared => Filter::Batch,
            Filter::Batch => Filter::All,
        }
    }

    /// Whether a module so labelled is shown.
    fn keeps(self, label: &Label) -> bool {
        match self {
            Filter::All => true,
            Filter::Shared => !label.is_batch(),
            Filter::Batch => label.is_batch(),
        }
    }

    /// What the header says of it; none when every module is shown.
    pub fn shown(self) -> Option<&'static str> {
        match self {
            Filter::All => None,
            Filter::Shared => Some("not batch-specific"),
            Filter::Batch => Some("batch-specific"),
        }
    }
}

/// A batch column's width at most, its name cut from the front to fit.
const BATCH_W: usize = 14;

/// Why there is no batches pane.
const ONE_BATCH: &str = "the cells are of one batch";

impl ModulesView<'_> {
    /// `b`: open the batches pane, or say on the footer why there is none.
    pub(super) fn open_batches(&mut self) {
        match &self.m.batches {
            Some(_) => self.batching = true,
            None => self.status = Some(ONE_BATCH.into()),
        }
    }

    /// A key while the batches pane is open: `b` or Esc closes it, `B`
    /// filters, ↑/↓ scroll, `q` quits.
    pub(super) fn batches_key(&mut self, key: KeyEvent) {
        self.status = None;
        match key.code {
            KeyCode::Char('b') | KeyCode::Esc => self.batching = false,
            KeyCode::Char('B') => self.next_filter(),
            KeyCode::Down | KeyCode::Char('j') => {
                self.batch_from = (self.batch_from + 1).min(self.m.cell_roles.saturating_sub(1))
            }
            KeyCode::Up | KeyCode::Char('k') => self.batch_from = self.batch_from.saturating_sub(1),
            KeyCode::Char('q') => self.done = true,
            _ => {}
        }
    }

    /// `B`: the next filter that leaves a cell module to show (every one,
    /// those not batch-specific, the batch-specific), the panes back at
    /// their first.
    pub(super) fn next_filter(&mut self) {
        let m = self.m;
        let Some(batches) = &m.batches else {
            self.status = Some(ONE_BATCH.into());
            return;
        };
        let mut filter = self.filter;
        self.shown = loop {
            filter = filter.next();
            let places: Vec<usize> = (0..m.cell_roles)
                .filter(|&k| filter.keeps(&batches.spread[m.columns[k]].label))
                .collect();
            if !places.is_empty() {
                break places;
            }
            self.status = Some(match filter {
                Filter::Batch => "no cell module is batch-specific".into(),
                _ => "every cell module is batch-specific".into(),
            });
        };
        self.filter = filter;
        self.col_from = 0;
        self.col_fit = self.shown.len();
        if let Sort::Column(k) = self.sort {
            if !self.shown.contains(&k) {
                self.resort(Sort::Module, false);
            }
        }
    }

    /// The batches pane in `area` above the footer in `foot`: each cell
    /// module's share in each batch (each batch's cells weighing alike),
    /// its evenness over each condition's replicates and its label, from
    /// the row scrolled to; then the conditions and the rule.
    pub(super) fn draw_batches(&self, frame: &mut Frame, area: Rect, foot: Rect) {
        let help = match &self.status {
            Some(said) => Line::from(Span::styled(format!(" {said}"), HIGHLIGHT)),
            None => help_line(&[
                ("b/Esc", "close"),
                ("B", "filter"),
                ("↑/↓", "scroll"),
                ("q", "quit"),
            ]),
        };
        frame.render_widget(Paragraph::new(help), foot);
        let Some(batches) = &self.m.batches else {
            return;
        };
        let block = panel(
            " batches · each cell module's share in each batch, its cells by batch size ".into(),
            true,
        );
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let d = &batches.design;
        let order = d.order();
        let w = (d
            .batches
            .iter()
            .map(|b| b.chars().count())
            .max()
            .unwrap_or(0))
        .clamp(5, BATCH_W);
        // A name too long keeps its end, where a replicate token sits.
        let cut = |s: &str| -> String {
            let n = s.chars().count();
            let s: String = match n > w {
                true => std::iter::once('…')
                    .chain(s.chars().skip(n + 1 - w))
                    .collect(),
                false => s.into(),
            };
            format!("{s:>w$} ")
        };
        let mut table = vec![Line::from(
            [Span::styled(" module ", DIM)]
                .into_iter()
                .chain(order.iter().map(|&b| Span::styled(cut(&d.batches[b]), DIM)))
                .chain([Span::styled("  even  label", DIM)])
                .collect::<Vec<_>>(),
        )];
        for &c in self.m.columns.iter().skip(self.batch_from) {
            let s = &batches.spread[c];
            let top = s.per_batch.iter().copied().fold(0.0, f64::max).max(1e-12);
            let shown = self.filter.keeps(&s.label);
            let mut spans = vec![Span::styled(
                format!(" C{:<5} ", c + 1),
                if shown { Style::default() } else { DIM },
            )];
            spans.extend(order.iter().map(|&b| {
                let p = s.per_batch[b];
                let g = (40.0 + 180.0 * (p / top).clamp(0.0, 1.0)).round() as u8;
                let (fill, ink) = gray_cell(g);
                Span::styled(
                    format!("{:>w$} ", format!("{:.0}%", 100.0 * p)),
                    Style::default().bg(fill).fg(ink),
                )
            }));
            let even = s.evenness.map_or("    -".into(), |e| format!("{e:>5.2}"));
            let style = match (&s.label, shown) {
                (_, false) => DIM,
                (Label::Batch(_), true) => HIGHLIGHT,
                _ => Style::default(),
            };
            spans.push(Span::styled(format!("  {even}  {}", s.label.word()), style));
            table.push(Line::from(spans));
        }
        let mut notes: Vec<Line> = (d.listing().into_iter())
            .map(|l| Line::from(Span::styled(format!(" {l}"), DIM)))
            .collect();
        if !d.replicated() {
            notes.push(Line::from(Span::styled(
                " no condition has two replicates: name them with --batch-groups",
                HIGHLIGHT,
            )));
        }
        notes.push(Line::from(Span::styled(
            format!(
                " even: over each condition's replicates, 1 evenly, 0 in one; under {EVEN:.2} \
                 batch-specific; even with {:.0}% in one condition, that condition's; under \
                 {FEWEST:.0} cells, too few",
                100.0 * CONDITION
            ),
            DIM,
        )));
        let [top, below] = Layout::vertical([
            Constraint::Fill(1),
            Constraint::Length(notes.len() as u16 + 2),
        ])
        .areas(inner);
        frame.render_widget(Paragraph::new(table), top);
        frame.render_widget(Paragraph::new(notes).wrap(Wrap { trim: false }), below);
    }
}
