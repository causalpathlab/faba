//! The view `faba modules` shows while it fits: a bar per chain, its
//! burn-in marked, with the time taken and left; under them each chain's
//! log joint as it goes; then whether a chain lags, and once summarised
//! the verdict on that and the chains' split R-hat, and what to do about
//! it (see [`super::advice`]). `q` cancels; once the fit writes, too late
//! for that, it leaves without opening the browser.

use super::advice::Verdict;
use super::Progress;
use crate::tui::{clock, on_thread, run_view, spinner_frame, Pending, View};
use data_beans::interactive::ui::{header, help_line, panel, Screen, ACCENTED, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use std::borrow::Cow;
use std::time::{Duration, Instant};
use trace::chain_style;

pub(in crate::cluster) mod trace;

/// How long the R-hat stays up once the fit is done, unless a key comes.
const HOLD: Duration = Duration::from_millis(1500);

/// The screen shown while the fit runs; done once its result is in (and
/// has been up a moment, unless it is to leave).
struct Fitting<'a, T> {
    title: &'a str,
    progress: &'a Progress,
    /// The fit's result, once it comes.
    pending: Pending<'a, anyhow::Result<T>>,
    /// A key came after the fit: open the browser now.
    skip: bool,
    /// A cancel came too late to stop the fit: leave once it is done,
    /// without opening the browser.
    leave: bool,
    status: Option<&'static str>,
    start: Instant,
}

impl<T> Fitting<'_, T> {
    /// One chain's line: its bar, `│` where the burn-in ends (burn-in
    /// filled in one colour, sampling in another), its sweeps and phase,
    /// then its time taken and left since it started (a chain queued for a
    /// thread waits), in `width` cells.
    fn chain_line(&self, i: usize, (total, warmup): (usize, usize), width: u16) -> Line<'static> {
        let done = self.progress.done(i).min(total);
        let running = self.progress.running(i);
        let phase = match running {
            None => "waiting ",
            Some(_) if done < warmup => "burn-in ",
            Some(_) => "sampling",
        };
        let unknown = Cow::Borrowed("--:--:--");
        let taken = running.map_or(unknown.clone(), |t| clock(t).into());
        let left = match running {
            Some(elapsed) if done > 0 => {
                clock(elapsed.mul_f64((total - done) as f64 / done as f64)).into()
            }
            _ => unknown,
        };
        let tail = format!(" {done}/{total} {phase}  {taken} taken  {left} left");
        // Numbered from 0 by its seed's offset, as the run record's best chain.
        let head = format!(" chain {i:<3}");
        let w = usize::from(width)
            .saturating_sub(head.len() + tail.len() + 1)
            .max(4);
        let mut spans = vec![Span::styled(head, chain_style(i))];
        spans.extend(chain_bar(done, total, warmup, w));
        spans.push(Span::styled(tail, DIM));
        Line::from(spans)
    }

    /// The lines under the plot: once every chain is done, whether one
    /// lags; once they are summarised, in its place the verdict on the
    /// R-hat and the traces, its advice, and what the fit is doing.
    fn after(
        &self,
        (total, warmup): (usize, usize),
        traces: &[Vec<(usize, f64)>],
    ) -> Vec<Line<'static>> {
        let Some(rhat) = self.progress.rhat() else {
            let done = (0..self.progress.chains()).all(|i| self.progress.done(i) >= total);
            return done
                .then(|| trace::verdict(traces, warmup))
                .flatten()
                .into_iter()
                .collect();
        };
        let verdict = Verdict::of(traces, warmup, total, rhat);
        let style = if verdict.converged() { DIM } else { HIGHLIGHT };
        let mut lines = vec![
            Line::raw(""),
            Line::from(Span::styled(format!(" {}", verdict.said()), style)),
            Line::from(Span::styled(format!(" {}", verdict.advice()), DIM)),
        ];
        lines.push(Line::from(Span::styled(
            if self.leave {
                " not opening the modules"
            } else if self.pending.out().is_some() {
                " opening the modules"
            } else {
                " writing the tables"
            },
            DIM,
        )));
        lines
    }

    /// The bars in `area`, the traces under them in what rows are left
    /// over (none when too few), then the lines after; the bars first, the
    /// lines after next, when rows run short.
    fn sampling(&self, frame: &mut Frame, area: Rect, schedule: (usize, usize)) {
        let mut bars: Vec<Line> = (0..self.progress.chains())
            .map(|i| self.chain_line(i, schedule, area.width))
            .collect();
        bars.push(Line::from(Span::styled(
            format!(" │ the burn-in ends at sweep {}", schedule.1),
            DIM,
        )));
        let traces = self.progress.traces();
        let after = self.after(schedule, &traces);
        let rows = |n: usize| u16::try_from(n).unwrap_or(u16::MAX);
        let top = rows(bars.len()).min(area.height);
        let bottom = rows(after.len()).min(area.height - top);
        let [bars_at, plot_at, after_at] = Layout::vertical([
            Constraint::Length(top),
            Constraint::Fill(1),
            Constraint::Length(bottom),
        ])
        .areas(area);
        frame.render_widget(Paragraph::new(bars), bars_at);
        trace::draw(frame, plot_at, &traces, schedule);
        frame.render_widget(Paragraph::new(after), after_at);
    }

    /// Stop the fit, or, too late for that, leave once it is written
    /// without opening the browser; before the chains, it stops once the
    /// step in hand (an input read, the filters) ends.
    fn cancel(&mut self) {
        let sampling = self.progress.schedule().is_some();
        self.status = Some(if !self.progress.cancel() {
            self.leave = true;
            "too late to cancel: leaving without the modules once the tables are written"
        } else if sampling {
            "cancelling: the chains stop and nothing is written"
        } else {
            "cancelling once the step in hand ends: nothing is written"
        });
    }
}

/// A `width`-cell bar for `done` of `total` sweeps, `│` at the end of the
/// `warmup` burn-in sweeps, the burn-in filled in one colour and the
/// sampling in another, the rest dim.
fn chain_bar(done: usize, total: usize, warmup: usize, width: usize) -> Vec<Span<'static>> {
    let cells = |n: usize| crate::tui::filled(n as u64, total as u64, width).min(width);
    let (mark, full) = (cells(warmup), cells(done));
    // The mark's cell, when the burn-in ends inside the bar.
    let at_mark = usize::from(mark > 0 && mark < width);
    // Where the sampling begins: after the mark.
    let after = (mark + at_mark).min(width);
    let burnt = full.min(mark);
    let sampled = full.max(after) - after;
    [
        ("█", burnt, ACCENTED),
        ("░", mark - burnt, DIM),
        ("│", at_mark, DIM),
        ("█", sampled, HIGHLIGHT),
        ("░", width - after - sampled, DIM),
    ]
    .into_iter()
    .filter(|&(_, n, _)| n > 0)
    .map(|(cell, n, style)| Span::styled(cell.repeat(n), style))
    .collect()
}

impl<T> View for Fitting<'_, T> {}

impl<T> Screen for Fitting<'_, T> {
    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        let [top, body, foot] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(1),
            Constraint::Length(1),
        ])
        .areas(area);
        frame.render_widget(
            header("modules", self.title, &clock(self.start.elapsed())),
            top,
        );
        let schedule = self.progress.schedule();
        let block = panel(
            if schedule.is_some() {
                " sampling "
            } else {
                " reading "
            }
            .into(),
            true,
        );
        let inner = block.inner(body);
        frame.render_widget(block, body);
        match schedule {
            None => {
                let spin = spinner_frame(self.start.elapsed());
                let line = Line::from(vec![
                    Span::styled(format!(" {spin} "), HIGHLIGHT),
                    Span::raw("reading and filtering the reads"),
                ]);
                frame.render_widget(Paragraph::new(line), inner);
            }
            Some(schedule) => self.sampling(frame, inner, schedule),
        }
        let footer = match self.status {
            Some(s) => Line::from(Span::styled(format!(" {s}"), HIGHLIGHT)),
            None => help_line(&[("q", "cancel")]),
        };
        frame.render_widget(Paragraph::new(footer), foot);
    }

    /// `q` or Esc cancels; once the fit is done any other key opens the
    /// browser.
    fn handle_key(&mut self, key: KeyEvent) {
        if matches!(key.code, KeyCode::Char('q') | KeyCode::Esc) {
            self.cancel();
        } else if self.pending.out().is_some() {
            self.skip = true;
        }
    }

    fn interrupt(&mut self) {
        self.cancel();
    }

    fn done(&self) -> bool {
        match self.pending.out() {
            None => self.pending.gone(),
            Some(Err(_)) => true,
            Some(Ok(_)) => {
                self.skip || self.leave || self.pending.since().is_some_and(|t| t >= HOLD)
            }
        }
    }

    /// Pick up the result, else redraw for the bars and the clock.
    fn tick(&mut self) -> bool {
        self.pending.poll();
        true
    }
}

/// Run `work` (a fit of `chains` chains counting its sweeps on the
/// progress it is given) on a thread, with the view under `title` until it
/// is done; its result, none when a cancel came too late to stop it (its
/// tables are written, but the browser is not to open).
pub fn watch<T: Send>(
    title: &str,
    chains: usize,
    work: impl FnOnce(&Progress) -> anyhow::Result<T> + Send,
) -> anyhow::Result<Option<T>> {
    watch_with(title, chains, work, |screen| run_view(screen))
}

/// [`watch`], the view shown by `show`; a fit that panics is an error.
fn watch_with<T: Send>(
    title: &str,
    chains: usize,
    work: impl FnOnce(&Progress) -> anyhow::Result<T> + Send,
    show: impl FnOnce(&mut Fitting<T>) -> anyhow::Result<()>,
) -> anyhow::Result<Option<T>> {
    let progress = Progress::traced(chains);
    let ((shown, leave), out) = on_thread(
        "the fit",
        || work(&progress),
        |rx| {
            let mut screen = Fitting {
                title,
                progress: &progress,
                pending: Pending::new(rx),
                skip: false,
                leave: false,
                status: None,
                start: Instant::now(),
            };
            // A view that cannot run leaves the fit to end on its own.
            let shown = show(&mut screen);
            if shown.is_err() {
                progress.cancel();
            }
            ((shown, screen.leave), screen.pending.take())
        },
    );
    shown?;
    let out = out?;
    Ok((!leave).then_some(out))
}

#[cfg(test)]
mod tests;
