//! A spinner on screen while slow work runs, so a view never leaves the
//! terminal blank while it reads.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::RecvTimeoutError;
use std::time::{Duration, Instant};

use data_beans::interactive::ui::{header, Screen, DIM, HIGHLIGHT};
use ratatui::crossterm::event::KeyEvent;
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use super::{bar_spans, on_thread, popup_frame, run_view, spinner_frame, Pending, View};

/// How long the work may take before the spinner shows, so quick work
/// never flickers the screen.
const GRACE: Duration = Duration::from_millis(150);

/// The screen shown while the work runs; done once its result arrives.
struct Busy<'a, T> {
    badge: &'a str,
    what: &'a str,
    start: Instant,
    /// The work's result, once it comes.
    pending: Pending<'a, T>,
    /// A key (or Ctrl-C) came: say it cannot stop the work.
    pressed: bool,
    /// Files done, of `total`; no bar when `total` is 0.
    done: &'a AtomicUsize,
    total: usize,
}

impl<T> View for Busy<'_, T> {}

impl<T> Screen for Busy<'_, T> {
    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        frame.render_widget(header(self.badge, "loading", ""), area);
        let elapsed = self.start.elapsed();
        let spin = spinner_frame(elapsed);
        let w = (self.what.chars().count() as u16 + 12)
            .max(if self.pressed { 56 } else { 40 })
            .min(area.width);
        let rows = if self.pressed { 5 } else { 4 };
        let inner = popup_frame(frame, area, w, rows, " working ".into());
        let secs = format!("{}s", elapsed.as_secs());
        let second = if self.total == 0 {
            Line::from(Span::styled(format!("   {secs}"), DIM))
        } else {
            let done = self.done.load(Ordering::Relaxed).min(self.total);
            let count = format!(" {done}/{} files  {secs}", self.total);
            let bar_w = (inner.width as usize).saturating_sub(count.len() + 4);
            let [on, off] = bar_spans(done as u64, self.total as u64, bar_w, HIGHLIGHT);
            Line::from(vec![Span::raw("   "), on, off, Span::styled(count, DIM)])
        };
        let mut lines = vec![
            Line::from(vec![
                Span::styled(format!(" {spin} "), HIGHLIGHT),
                Span::raw(self.what.to_string()),
            ]),
            second,
        ];
        if self.pressed {
            lines.push(Line::from(Span::styled(
                "   it cannot stop halfway; the view opens when done",
                DIM,
            )));
        }
        frame.render_widget(Paragraph::new(lines), inner);
    }

    /// Keys wait: the work cannot be stopped halfway, which the pop-up
    /// then says.
    fn handle_key(&mut self, _key: KeyEvent) {
        self.pressed = true;
    }

    fn interrupt(&mut self) {
        self.pressed = true;
    }

    fn done(&self) -> bool {
        self.pending.out().is_some() || self.pending.gone()
    }

    /// Pick up the result, else redraw for the spinner and the clock.
    fn tick(&mut self) -> bool {
        self.pending.poll();
        true
    }
}

/// Run `work`, and when `show` (and there is a terminal) and it takes more
/// than a moment, show a spinner saying `what` it is doing under `badge`
/// until it is done.
pub fn busy<T: Send>(
    show: bool,
    badge: &str,
    what: &str,
    work: impl FnOnce() -> anyhow::Result<T> + Send,
) -> anyhow::Result<T> {
    busy_counting(show, badge, what, 0, |_| work())
}

/// [`busy`], with a bar of `total` files that `work` counts off on the
/// counter it is given.
pub fn busy_counting<T: Send>(
    show: bool,
    badge: &str,
    what: &str,
    total: usize,
    work: impl FnOnce(&AtomicUsize) -> anyhow::Result<T> + Send,
) -> anyhow::Result<T> {
    let done = AtomicUsize::new(0);
    if !show || !data_beans::interactive::tui_available() {
        return work(&done);
    }
    let counter = &done;
    let (shown, out) = on_thread(
        &format!("{what}: the work"),
        || work(counter),
        |rx| {
            match rx.recv_timeout(GRACE) {
                Ok(out) => return (Ok(()), Some(out)),
                Err(RecvTimeoutError::Disconnected) => return (Ok(()), None),
                Err(RecvTimeoutError::Timeout) => {}
            }
            let mut screen = Busy {
                badge,
                what,
                start: Instant::now(),
                pending: Pending::new(rx),
                pressed: false,
                done: counter,
                total,
            };
            let shown = run_view(&mut screen);
            (shown, screen.pending.take())
        },
    );
    shown?;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc::channel;

    #[test]
    fn draws_at_any_terminal_size_and_ends_with_the_result() {
        let (tx, rx) = channel();
        let done = AtomicUsize::new(1);
        let mut b = Busy {
            badge: "pileup",
            what: "reading GENE1",
            start: Instant::now(),
            pending: Pending::new(&rx),
            pressed: false,
            done: &done,
            total: 3,
        };
        b.tick();
        assert!(!b.done());
        for (w, h) in [(10, 3), (80, 24), (200, 60)] {
            let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
            t.draw(|f| b.render(f)).unwrap();
        }
        tx.send(7).unwrap();
        b.tick();
        assert!(b.done());
        assert_eq!(b.pending.out(), Some(&7));
    }

    #[test]
    fn work_that_dies_ends_the_spinner() {
        let (tx, rx) = channel::<i32>();
        let done = AtomicUsize::new(0);
        let mut b = Busy {
            badge: "pwm",
            what: "collecting",
            start: Instant::now(),
            pending: Pending::new(&rx),
            pressed: false,
            done: &done,
            total: 0,
        };
        drop(tx);
        b.tick();
        assert!(b.done() && b.pending.out().is_none(), "not left spinning");
    }

    #[test]
    fn without_showing_it_just_runs() {
        assert_eq!(busy(false, "x", "y", || Ok(3)).unwrap(), 3);
        assert!(busy(false, "x", "y", || anyhow::bail!("no")
            as anyhow::Result<()>)
        .is_err());
    }
}
