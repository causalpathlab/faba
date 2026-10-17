//! The fit's saved trace, full pane on `t`: each chain's log joint by
//! sweep as the fitting view drew it, the chains named in their colours,
//! and under it the verdict and its advice (see
//! [`crate::cluster::progress::advice`]).

use super::ModulesView;
use crate::cluster::progress::view::trace::{chain_style, draw};
use data_beans::interactive::ui::{help_line, panel, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

impl ModulesView<'_> {
    /// `t`: open the trace, or say on the footer why there is none.
    pub(super) fn open_trace(&mut self) {
        match &self.m.fit.trace {
            Ok(_) => self.tracing = true,
            Err(why) => self.status = Some(why.clone()),
        }
    }

    /// A key while the trace is open: `t` or Esc closes it, `q` quits.
    pub(super) fn trace_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('t') | KeyCode::Esc => self.tracing = false,
            KeyCode::Char('q') => self.done = true,
            _ => {}
        }
    }

    /// The header's word on convergence, its R-hat and `ok` or `check`;
    /// none without a record of the R-hat.
    pub(super) fn badge(&self) -> Option<String> {
        self.m.fit.verdict().map(|v| v.badge())
    }

    /// The trace pane in `area` above the footer in `foot`: the plot over
    /// the sweeps run, then the chains' names, where the burn-in ends, the
    /// verdict and its advice.
    pub(super) fn draw_trace(&self, frame: &mut Frame, area: Rect, foot: Rect) {
        let help = help_line(&[("t/Esc", "close"), ("q", "quit")]);
        frame.render_widget(Paragraph::new(help), foot);
        let (Ok(t), Some(verdict)) = (&self.m.fit.trace, self.m.fit.verdict()) else {
            return;
        };
        let block = panel(" trace · each chain's log joint by sweep ".into(), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let mut names = vec![Span::raw(" ")];
        names.extend(
            (0..t.chains.len()).map(|i| Span::styled(format!("chain {i}  "), chain_style(i))),
        );
        let style = if verdict.converged() { DIM } else { HIGHLIGHT };
        let under = vec![
            Line::from(names),
            Line::from(Span::styled(
                format!(" │ the burn-in ends at sweep {}", t.burn_in),
                DIM,
            )),
            Line::raw(""),
            Line::from(Span::styled(format!(" {}", verdict.said()), style)),
            Line::from(Span::styled(format!(" {}", verdict.advice()), DIM)),
        ];
        let [plot, below] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(under.len() as u16)])
                .areas(inner);
        draw(frame, plot, &t.chains, (self.m.fit.span(), t.burn_in));
        frame.render_widget(Paragraph::new(under), below);
    }
}
