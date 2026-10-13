//! Several genes' browsers stacked top to bottom, for the genes marked in
//! the gene list: one takes the keys (Tab moves on), how the bars are drawn
//! applies to all of them, and `s` saves them all as one figure.

use data_beans::interactive::ui::{help_line, Screen};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::Line;
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use super::{browser, Exit, PileupView, Settings, Track, View};
use crate::figure::{status_line, Controls, Key};

/// What a stack says when a key needs one gene to itself.
const ONE_GENE: &str = "that works on one gene: g goes back to the list to open it alone";

/// How the genes are laid out (`|` cycles), on screen and in the figure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arrange {
    /// One under another.
    Column,
    /// Side by side.
    Row,
    /// About as many across as down.
    Grid,
}

impl Arrange {
    fn next(self) -> Self {
        match self {
            Arrange::Column => Arrange::Row,
            Arrange::Row => Arrange::Grid,
            Arrange::Grid => Arrange::Column,
        }
    }

    /// Genes to a row, of `n`.
    fn cols(self, n: usize) -> usize {
        match self {
            Arrange::Column => 1,
            Arrange::Row => n.max(1),
            Arrange::Grid => (n as f64).sqrt().ceil().max(1.0) as usize,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Arrange::Column => "one under another",
            Arrange::Row => "side by side",
            Arrange::Grid => "a grid",
        }
    }
}

pub struct StackView<'a> {
    views: Vec<PileupView<'a>>,
    /// The view the keys go to.
    focus: usize,
    arrange: Arrange,
    /// Saving the whole stack, and images for all.
    controls: Controls,
    status: Option<String>,
    exit: Option<Exit>,
}

impl<'a> StackView<'a> {
    pub fn new(views: Vec<PileupView<'a>>) -> Self {
        let controls = Controls::new(&save_name(&views));
        let mut stack = Self {
            views,
            focus: 0,
            arrange: Arrange::Column,
            controls,
            status: None,
            exit: None,
        };
        stack.mark_focus();
        stack
    }

    /// The badge says which gene has the keys.
    fn mark_focus(&mut self) {
        for (i, v) in self.views.iter_mut().enumerate() {
            v.badge = if i == self.focus {
                "▸ keys"
            } else {
                "pileup"
            };
        }
    }

    /// Every gene drawn as the one with the keys is, so they compare.
    fn share_settings(&mut self) {
        let mut s = self.views[self.focus].settings();
        s.images = Some(self.controls.images_on());
        for v in &mut self.views {
            v.apply(&s);
        }
    }

    /// The genes' figures as one, laid out as on screen.
    pub(super) fn figure(&self) -> String {
        let figures: Vec<String> = self.views.iter().map(PileupView::figure).collect();
        crate::figure::tile_svgs(&figures, self.arrange.cols(figures.len()))
    }

    pub(super) fn footer_line(&self) -> Line<'static> {
        if let Some(msg) = &self.status {
            return status_line(msg);
        }
        if let Some(line) = self.controls.footer() {
            return line;
        }
        let v = &self.views[self.focus];
        let mut keys = vec![
            ("Tab/click", "next gene"),
            ("|", self.arrange.next().name()),
            ("←/→", "bar"),
            ("n/N", "next/prev site"),
            ("=/-", "zoom"),
            ("+/_", "zoom all"),
            ("0", "whole"),
            (")", "all whole"),
            ("y", "scale"),
            (
                "Y",
                if v.y_lock.is_some() {
                    "unlock y"
                } else {
                    "lock y (all genes)"
                },
            ),
        ];
        if v.tracks.iter().any(|t| t.is_reads()) {
            keys.push(("c", v.show.next().short()));
        }
        self.controls.help_keys(&mut keys);
        keys.extend([("g/Esc", "gene list"), ("q", "quit")]);
        help_line(&keys)
    }
}

impl crate::tui::View for StackView<'_> {
    /// The apply key submits the save name being typed.
    fn takes_apply(&self) -> bool {
        self.controls.prompting()
    }
}

impl Screen for StackView<'_> {
    fn takes_mouse(&self) -> bool {
        true
    }

    /// A press gives the gene under it the keys; then, as in one gene's
    /// browser, a click on its bars moves its cursor, a drag pans it, and
    /// the wheel pans the gene under the pointer.
    fn mouse(&mut self, m: ratatui::crossterm::event::MouseEvent) -> bool {
        use ratatui::crossterm::event::MouseEventKind;
        let at = ratatui::layout::Position::new(m.column, m.row);
        let under = self.views.iter().position(|v| v.drawn.contains(at));
        let i = match m.kind {
            // Let go anywhere: the gene that was grabbed.
            MouseEventKind::Up(_) => self.focus,
            MouseEventKind::Down(_) => {
                let Some(i) = under else { return false };
                self.focus = i;
                self.status = None;
                self.mark_focus();
                i
            }
            _ => match under {
                Some(i) => i,
                None => return false,
            },
        };
        let changed = self.views[i].pointer(m);
        changed || matches!(m.kind, MouseEventKind::Down(_))
    }

    fn done(&self) -> bool {
        self.exit.is_some()
    }

    fn interrupt(&mut self) {
        self.exit = Some(Exit::Quit);
    }

    /// The annotation arriving redraws every gene.
    fn tick(&mut self) -> bool {
        // Every gene takes its annotation, so none may be skipped.
        let mut any = false;
        for v in &mut self.views {
            any |= v.tick();
        }
        any
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.status = None;
        match self.controls.key(key) {
            Key::Used => return self.share_settings(),
            Key::Save(prefix) => return self.controls.save(&self.figure(), &prefix),
            Key::Pass => {}
        }
        let n = self.views.len();
        match key.code {
            KeyCode::Tab => self.focus = (self.focus + 1) % n,
            KeyCode::Char('|') => {
                self.arrange = self.arrange.next();
                for v in &mut self.views {
                    v.plots
                        .iter_mut()
                        .for_each(crate::figure::term::PlotImage::invalidate);
                }
            }
            KeyCode::BackTab => self.focus = (self.focus + n - 1) % n,
            KeyCode::Char('q') => self.exit = Some(Exit::Quit),
            KeyCode::Char('g') | KeyCode::Esc => self.exit = Some(Exit::Genes),
            // A search or the other modality opens one gene, not a stack.
            KeyCode::Char('/' | 'm' | 'M') => self.status = Some(ONE_GENE.into()),
            // Shifted, the zoom keys act on every gene: `+` (Shift+=) in,
            // `_` (Shift+-) out, `)` (Shift+0) each whole.
            KeyCode::Char(c @ ('+' | '_' | ')')) => {
                let one = match c {
                    '+' => '=',
                    '_' => '-',
                    _ => '0',
                };
                let key = KeyEvent::new(KeyCode::Char(one), key.modifiers);
                for v in &mut self.views {
                    v.handle_key(key);
                    // Out past a gene would reload one gene's span alone.
                    if v.exit.take().is_some() {
                        self.status = Some(ONE_GENE.into());
                    }
                }
            }
            // One y axis for every gene: the tallest of them, as a locked
            // maximum only ever raises an axis (no bar is cut off).
            KeyCode::Char('Y') if self.views[self.focus].y_lock.is_none() => {
                let tallest = self
                    .views
                    .iter()
                    .filter_map(|v| v.bins().shared)
                    .fold(0.0, f64::max);
                self.views[self.focus].y_lock = (tallest > 0.0).then_some(tallest);
                self.share_settings();
            }
            _ => {
                let v = &mut self.views[self.focus];
                v.handle_key(key);
                // Zooming out past the gene reloads one gene's span.
                if v.exit.take().is_some() {
                    self.status = Some(ONE_GENE.into());
                }
                self.share_settings();
            }
        }
        self.mark_focus();
    }

    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        let help = crate::tui::wrap_help(self.footer_line(), area.width);
        let [body, footer] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(help.len() as u16)])
                .areas(area);
        let n = self.views.len();
        let cols = self.arrange.cols(n);
        let rows = n.div_ceil(cols);
        let row_areas = Layout::vertical(vec![Constraint::Ratio(1, rows as u32); rows]).split(body);
        for (r, chunk) in self.views.chunks_mut(cols).enumerate() {
            let cells = Layout::horizontal(vec![Constraint::Ratio(1, cols as u32); cols])
                .split(row_areas[r]);
            for (v, &cell) in chunk.iter_mut().zip(cells.iter()) {
                v.render_in(frame, cell, false);
            }
        }
        frame.render_widget(Paragraph::new(help), footer);
        self.controls.render_save(frame, area);
    }
}

/// The name a save suggests: `pileup_` and the first genes' names.
fn save_name(views: &[PileupView]) -> String {
    let mut names: Vec<&str> = views.iter().map(|v| v.title.as_str()).take(3).collect();
    if views.len() > 3 {
        names.push("etc");
    }
    format!("pileup_{}", names.join("_"))
}

/// Browse the genes of `views` stacked, each with its tracks, until the
/// user goes back to the gene list or quits. `settings` holds how the gene
/// before was drawn, and is left as the stack was.
pub fn show_stack(views: Vec<(View, Vec<Track>)>, settings: &mut Settings) -> anyhow::Result<Exit> {
    let views: Vec<PileupView> = views
        .into_iter()
        .map(|(view, tracks)| {
            let mut b = browser(view, tracks, None);
            b.apply(settings);
            b
        })
        .collect();
    anyhow::ensure!(!views.is_empty(), "no gene to show");
    let mut stack = StackView::new(views);
    // Text plots on screen, as one gene's browser; the saved figure is a picture.
    stack.controls = Controls::new(&save_name(&stack.views));
    if let Some(on) = settings.images {
        stack.controls.set_images(on);
    }
    stack.share_settings();
    crate::tui::run_view(&mut stack)?;
    *settings = stack.views[stack.focus].settings();
    Ok(stack.exit.unwrap_or(Exit::Quit))
}

#[cfg(test)]
#[path = "../tests/stack.rs"]
mod tests;
