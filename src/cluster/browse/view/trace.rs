//! The fit's saved trace, full pane on `t`: each restart's ELBO by step
//! as the fitting view drew it, the restarts named in their colours, and
//! under it the verdict and its advice (see
//! [`crate::cluster::progress::advice`]).

use super::ModulesView;
use crate::cluster::progress::view::trace::{draw, restart_style};
use crate::tui::ui::{help_line, panel, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

impl ModulesView<'_> {
    /// `t`: open the trace, or say on the footer why there is none.
    pub(super) fn open_trace(&mut self) {
        match &self.m.fit.trace {
            Ok(_) => self.full = Some(super::Full::Trace),
            Err(why) => self.status = Some(why.clone()),
        }
    }

    /// A key while the trace is open: `t` or Esc closes it, `q` quits.
    pub(super) fn trace_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('t') | KeyCode::Esc => self.full = None,
            KeyCode::Char('q') => self.done = true,
            _ => {}
        }
    }

    /// The header's word on convergence, `ELBO ok` or `ELBO check`; none
    /// without a record of the fit.
    pub(super) fn badge(&self) -> Option<String> {
        self.m.fit.verdict().map(|v| v.badge())
    }

    /// The trace pane in `area` above the footer in `foot`: the plot over
    /// the steps run, then the restarts' names, where the k-means start
    /// ends, and (with a record of the fit) the verdict and its advice.
    pub(super) fn draw_trace(&self, frame: &mut Frame, area: Rect, foot: Rect) {
        let help = help_line(&[("t/Esc", "close"), ("q", "quit")]);
        frame.render_widget(Paragraph::new(help), foot);
        let Ok(t) = &self.m.fit.trace else {
            return;
        };
        let block = panel(" trace · each restart's ELBO by epoch ".into(), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let mut names = vec![Span::raw(" ")];
        names.extend(
            (0..t.restarts.len()).map(|i| Span::styled(format!("restart {i}  "), restart_style(i))),
        );
        let mut under = vec![
            Line::from(names),
            Line::from(Span::styled(
                crate::cluster::progress::vi_caption(t.warmup),
                DIM,
            )),
        ];
        // Without a record, the trace alone: no verdict to say.
        if let Some(verdict) = self.m.fit.verdict() {
            let style = if verdict.converged() { DIM } else { HIGHLIGHT };
            under.extend([
                Line::raw(""),
                Line::from(Span::styled(format!(" {}", verdict.said()), style)),
                Line::from(Span::styled(format!(" {}", verdict.advice()), DIM)),
            ]);
        }
        let [plot, below] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(under.len() as u16)])
                .areas(inner);
        draw(frame, plot, &t.restarts, (self.m.fit.span(), t.warmup));
        frame.render_widget(Paragraph::new(under), below);
    }
}
