//! A spinner on screen while slow work runs, so a view never leaves the
//! terminal blank while it reads.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{channel, Receiver, RecvTimeoutError};
use std::time::{Duration, Instant};

use data_beans::interactive::ui::{header, Screen, DIM, HIGHLIGHT};
use ratatui::crossterm::event::KeyEvent;
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use super::{bar_spans, popup_frame, run_view, View, SPINNER};

/// How long the work may take before the spinner shows, so quick work
/// never flickers the screen.
const GRACE: Duration = Duration::from_millis(150);

/// The screen shown while the work runs; done once its result arrives.
struct Busy<'a, T> {
    badge: &'a str,
    what: &'a str,
    start: Instant,
    rx: &'a Receiver<T>,
    out: Option<T>,
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
        let frames: Vec<char> = SPINNER.chars().collect();
        let spin = frames[(elapsed.as_millis() / 200) as usize % frames.len()];
        let w = (self.what.chars().count() as u16 + 12)
            .max(40)
            .min(area.width);
        let inner = popup_frame(frame, area, w, 4, " working ".into());
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
        let lines = vec![
            Line::from(vec![
                Span::styled(format!(" {spin} "), HIGHLIGHT),
                Span::raw(self.what.to_string()),
            ]),
            second,
        ];
        frame.render_widget(Paragraph::new(lines), inner);
    }

    /// Keys wait: the work cannot be stopped halfway.
    fn handle_key(&mut self, _key: KeyEvent) {}

    fn interrupt(&mut self) {}

    fn done(&self) -> bool {
        self.out.is_some()
    }

    /// Pick up the result, else redraw for the spinner and the clock.
    fn tick(&mut self) -> bool {
        if let Ok(out) = self.rx.try_recv() {
            self.out = Some(out);
        }
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
    let (tx, rx) = channel();
    let counter = &done;
    std::thread::scope(|s| {
        s.spawn(move || {
            // The receiver outlives the scope, so the send cannot fail.
            let _ = tx.send(work(counter));
        });
        match rx.recv_timeout(GRACE) {
            Ok(out) => return out,
            Err(RecvTimeoutError::Disconnected) => anyhow::bail!("{what}: the work panicked"),
            Err(RecvTimeoutError::Timeout) => {}
        }
        let mut screen = Busy {
            badge,
            what,
            start: Instant::now(),
            rx: &rx,
            out: None,
            done: &done,
            total,
        };
        run_view(&mut screen)?;
        match screen.out {
            Some(out) => out,
            None => rx
                .recv()
                .map_err(|_| anyhow::anyhow!("{what}: the work panicked"))?,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn draws_at_any_terminal_size_and_ends_with_the_result() {
        let (tx, rx) = channel();
        let done = AtomicUsize::new(1);
        let mut b = Busy {
            badge: "pileup",
            what: "reading GENE1",
            start: Instant::now(),
            rx: &rx,
            out: None,
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
        assert_eq!(b.out, Some(7));
    }

    #[test]
    fn without_showing_it_just_runs() {
        assert_eq!(busy(false, "x", "y", || Ok(3)).unwrap(), 3);
        assert!(busy(false, "x", "y", || anyhow::bail!("no")
            as anyhow::Result<()>)
        .is_err());
    }
}
