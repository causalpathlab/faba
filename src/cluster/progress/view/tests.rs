//! The fitting view drawn from a made-up progress: a bar per restart with
//! its k-means start marked, the restarts' ELBO traced under them, the
//! verdict once in, `q` cancelling, Ctrl+C stopping the restarts early.

use super::*;
use crate::cluster::progress::advice::Verdict;
use crate::tui::screen_text as screen;
use ratatui::crossterm::event::KeyModifiers;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::time::Duration;

/// The verdict of a fit that settled in 30 epochs of 60, its 3 restarts'
/// ELBO per read within 0.001.
fn settled() -> Verdict {
    Verdict::vi(true, false, (30, 60), 3, (0.001, "read"))
}

/// Three restarts of 100 steps, 40 of them the k-means start: one started
/// 10 s ago at step 20, one 7 s ago at step 70, one not yet.
fn midway() -> Progress {
    let p = Progress::new(3);
    let ago = |s| Instant::now() - Duration::from_secs(s);
    let _ = p.schedule.set((100, 40));
    let _ = p.started[0].set(ago(10));
    let _ = p.started[1].set(ago(7));
    p.done[0].0.store(20, Ordering::Relaxed);
    p.done[1].0.store(70, Ordering::Relaxed);
    p
}

fn fitting<'a>(p: &'a Progress, rx: &'a Receiver<anyhow::Result<u8>>) -> Fitting<'a, u8> {
    Fitting {
        title: "fitting the site blockmodel",
        progress: p,
        pending: Pending::new(rx),
        skip: false,
        status: None,
        start: Instant::now(),
    }
}

fn key(c: char) -> KeyEvent {
    KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE)
}

fn channel_u8() -> (Sender<anyhow::Result<u8>>, Receiver<anyhow::Result<u8>>) {
    channel()
}

/// The k-means start in one colour up to its mark, the epochs in another
/// after it, the rest dim.
#[test]
fn a_restarts_bar_marks_where_its_start_ends() {
    let text = |spans: &[Span]| spans.iter().map(|s| s.content.as_ref()).collect::<String>();
    let early = chain_bar(20, 100, 40, 10);
    assert_eq!(text(&early), "██░░│░░░░░");
    assert_eq!(early[0].style, ACCENTED);
    let late = chain_bar(70, 100, 40, 10);
    assert_eq!(text(&late), "████│██░░░");
    assert_eq!(
        late.iter().map(|s| s.style).collect::<Vec<_>>(),
        [ACCENTED, DIM, HIGHLIGHT, DIM]
    );
    assert_eq!(text(&chain_bar(0, 100, 0, 5)), "░░░░░", "no start, no mark");
    assert_eq!(text(&chain_bar(100, 100, 40, 5)), "██│██");
}

/// Each restart's bar, steps and phase, its time taken and left counted
/// from when it started, one not started as waiting, and the key to
/// cancel; at any terminal size.
#[test]
fn draws_a_bar_per_restart_with_its_start() {
    let (_tx, rx) = channel_u8();
    let p = midway();
    let mut v = fitting(&p, &rx);
    let s = screen(&mut v, 110, 10);
    for want in [
        "restart 0",
        "restart 1",
        "20/100 k-means",
        "70/100 epochs",
        "00:00:10 taken",
        "00:00:40 left",
        "00:00:07 taken",
        "00:00:03 left",
        "restart 2",
        "0/100 waiting",
        "the k-means start ends at step 40",
        "q cancel",
        "Ctrl+C stop early",
    ] {
        assert!(s.contains(want), "{want:?} in\n{s}");
    }
    // Inside the panel's borders, each restart's bar has its mark.
    for line in s.lines().filter(|l| l.contains(" restart ")) {
        let inside = line.trim_end().trim_matches('│');
        assert!(inside.contains('│'), "{line}");
    }
    for (w, h) in [(10, 3), (40, 5), (200, 60)] {
        screen(&mut v, w, h);
    }
}

/// Before the restarts start, it says it is reading; `q` then waits for the
/// input in hand, and says so.
#[test]
fn says_it_reads_before_the_restarts_start() {
    let (_tx, rx) = channel_u8();
    let p = Progress::new(2);
    let mut v = fitting(&p, &rx);
    assert!(screen(&mut v, 80, 6).contains("reading and filtering the reads"));
    v.handle_key(key('q'));
    assert!(p.cancelled());
    let s = screen(&mut v, 80, 6);
    assert!(s.contains("cancelling once the step in hand ends"), "{s}");
}

/// `q` cancels the fit, and the footer says so; once the fit writes, `q`
/// is refused, and says why.
#[test]
fn q_cancels_until_the_fit_writes() {
    let (_tx, rx) = channel_u8();
    let p = midway();
    let mut v = fitting(&p, &rx);
    v.handle_key(key('q'));
    assert!(p.cancelled());
    assert!(screen(&mut v, 100, 8).contains("cancelling"));
    assert!(!v.done(), "until the fit returns");

    let p = midway();
    p.commit().unwrap();
    let mut v = fitting(&p, &rx);
    v.handle_key(key('q'));
    assert!(!p.cancelled());
    assert!(screen(&mut v, 100, 8).contains("too late to cancel"));
}

/// The verdict shows once the fit ends; the view ends after the result
/// on a key, or at once on an error.
#[test]
fn shows_the_verdict_then_ends() {
    let (tx, rx) = channel_u8();
    let p = midway();
    p.optimised(settled()).unwrap();
    let mut v = fitting(&p, &rx);
    let s = screen(&mut v, 100, 12);
    assert!(
        s.contains("ELBO settled in 30 epochs") && s.contains("nothing to rerun"),
        "{s}"
    );
    assert!(s.contains("writing the tables"), "{s}");
    tx.send(Ok(3)).unwrap();
    v.tick();
    assert!(!v.done(), "the verdict stays up");
    v.handle_key(key('x'));
    assert!(v.done());

    let (tx, rx) = channel_u8();
    let mut v = fitting(&p, &rx);
    tx.send(Err(crate::cluster::progress::Cancelled.into()))
        .unwrap();
    v.tick();
    assert!(v.done(), "an error ends it at once");

    let (tx, rx) = channel_u8();
    let mut v = fitting(&p, &rx);
    drop(tx);
    v.tick();
    assert!(v.done(), "a fit that died ends it");
}

/// Ticks `view` until it is done, as the terminal loop would, a key
/// coming once the result is in (the verdict stays up until one).
fn tick_until_done(view: &mut Fitting<u8>) -> anyhow::Result<()> {
    while !view.done() {
        view.tick();
        if view.pending.out().is_some() {
            view.handle_key(key('x'));
        }
        std::thread::sleep(Duration::from_millis(1));
    }
    Ok(())
}

/// A fit that panics ends the view and comes back as an error, not as the
/// panic again.
#[test]
fn a_fit_that_panics_is_an_error() {
    let p = Progress::new(1);
    let out = std::panic::catch_unwind(|| {
        watch_with(
            "t",
            &p,
            true,
            |_| -> anyhow::Result<u8> { panic!("a fit") },
            tick_until_done,
        )
    });
    let err = out.expect("no panic past the view").unwrap_err();
    assert!(err.to_string().contains("panicked: a fit"), "{err}");
    let out = watch_with("t", &Progress::new(1), true, |_| Ok(7), tick_until_done);
    assert_eq!(out.unwrap(), 7);
}

/// A view that is not to stay (one fit of several) leaves once its fit is
/// in, with no key.
#[test]
fn a_view_not_to_stay_leaves_once_the_fit_is_in() {
    let out = watch_with(
        "t",
        &Progress::new(1),
        false,
        |_| Ok(7),
        |v| {
            while !v.done() {
                v.tick();
                std::thread::sleep(Duration::from_millis(1));
            }
            Ok(())
        },
    );
    assert_eq!(out.unwrap(), 7);
}

/// Once the fit is done the verdict stays up, however long, until a key.
#[test]
fn the_verdict_stays_up_until_a_key() {
    let (tx, rx) = channel_u8();
    let p = midway();
    p.optimised(settled()).unwrap();
    let mut v = fitting(&p, &rx);
    tx.send(Ok(3)).unwrap();
    v.tick();
    screen(&mut v, 100, 12);
    std::thread::sleep(Duration::from_millis(50));
    assert!(!v.done(), "the verdict stays up once the fit is done");
    v.handle_key(key('x'));
    assert!(v.done(), "a key leaves");
}

/// Ctrl+C while the restarts run stops them early, and the footer says
/// what they reached is written, which `q` after does not undo; before the
/// restarts, with nothing reached, it cancels.
#[test]
fn ctrl_c_stops_the_restarts_early() {
    let (_tx, rx) = channel_u8();
    let p = midway();
    let mut v = fitting(&p, &rx);
    v.interrupt();
    assert_eq!(p.stopped_early(), Some(vec![20, 70, 0]));
    assert!(!p.cancelled());
    let s = screen(&mut v, 100, 8);
    assert!(s.contains("stopping early"), "{s}");
    assert!(!v.done(), "until the fit returns");
    // `q` after it does not undo it: what they reached is written.
    v.handle_key(key('q'));
    assert!(!p.cancelled() && p.stopped_early().is_some());
    assert!(screen(&mut v, 100, 8).contains("stopping early already"));

    let p = Progress::new(2);
    let mut v = fitting(&p, &rx);
    v.interrupt();
    assert!(p.cancelled() && p.stopped_early().is_none());
    assert!(screen(&mut v, 80, 6).contains("cancelling once the step in hand ends"));
}

/// `q` or Ctrl+C once the fit writes cannot stop it: the view says so and
/// ends as soon as the result is in, without the hold; the result is the
/// fit's.
#[test]
fn a_stop_too_late_leaves_once_written() {
    for (stop, said) in [
        (
            (|v: &mut Fitting<u8>| v.handle_key(key('q'))) as fn(&mut Fitting<u8>),
            "too late to cancel",
        ),
        (|v: &mut Fitting<u8>| v.interrupt(), "too late to stop"),
    ] {
        let (tx, rx) = channel_u8();
        let p = midway();
        p.optimised(settled()).unwrap();
        let mut v = fitting(&p, &rx);
        stop(&mut v);
        assert!(!p.cancelled() && p.stopped_early().is_none());
        let s = screen(&mut v, 100, 12);
        assert!(s.contains(said), "{s}");
        tx.send(Ok(3)).unwrap();
        v.tick();
        assert!(v.done(), "no hold before leaving");
    }
    let p = Progress::new(1);
    let left = watch_with(
        "t",
        &p,
        true,
        |p| {
            p.commit()?;
            Ok(7)
        },
        |v| {
            while v.pending.out().is_none() {
                v.tick();
            }
            v.interrupt();
            tick_until_done(v)
        },
    );
    assert_eq!(left.unwrap(), 7, "the tables written are the fit's");
}

/// A made-up trace: steps `from..=to`, the ELBO `level` give or take a
/// fixed wobble under 1.
fn wobble(level: f64, from: usize, to: usize) -> Vec<(usize, f64)> {
    (from..=to)
        .map(|n| (n, level + ((n * 7919) % 13) as f64 / 13.0 - 0.5))
        .collect()
}

/// Three restarts of 100 steps, 40 of them the k-means start, all done:
/// each climbs from far below in its first epochs, then sits at
/// `levels[i]`.
fn finished(levels: [f64; 3]) -> Progress {
    let p = Progress::traced(3);
    let _ = p.schedule.set((100, 40));
    for (i, level) in levels.into_iter().enumerate() {
        let _ = p.started[i].set(Instant::now() - Duration::from_secs(5));
        p.done[i].0.store(100, Ordering::Relaxed);
        let mut t = p.traces[i].lock().unwrap();
        t.extend((41..=50).map(|n| (n, level - 1000.0 * 0.5f64.powi(n as i32 - 40))));
        t.extend(wobble(level, 51, 100));
    }
    p
}

/// The view drawn `w` x `h`, cell by cell.
fn drawn(v: &mut Fitting<u8>, w: u16, h: u16) -> ratatui::buffer::Buffer {
    let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
    t.draw(|f| v.render(f)).unwrap();
    t.backend().buffer().clone()
}

/// Where the buffer has a braille dot, with its style.
fn braille(buf: &ratatui::buffer::Buffer) -> Vec<(u16, u16, ratatui::style::Style)> {
    let a = buf.area;
    (a.top()..a.bottom())
        .flat_map(|y| (a.left()..a.right()).map(move |x| (x, y)))
        .filter(|&(x, y)| {
            buf[(x, y)]
                .symbol()
                .chars()
                .any(|c| ('\u{2801}'..='\u{28ff}').contains(&c))
        })
        .map(|(x, y)| (x, y, buf[(x, y)].style()))
        .collect()
}

/// Under the bars, a line per restart in its label's colour, the k-means
/// start's end a line down the plot.
#[test]
fn traces_each_restart_under_the_bars() {
    let (_tx, rx) = channel_u8();
    let p = finished([0.0, 10.0, -10.0]);
    let mut v = fitting(&p, &rx);
    let buf = drawn(&mut v, 80, 24);
    let dots = braille(&buf);
    for i in 0..3 {
        let colour = chain_style(i).fg;
        assert!(colour.is_some(), "restart {i} has a colour");
        assert!(
            dots.iter().any(|d| d.2.fg == colour),
            "restart {i}'s trace drawn"
        );
        // The bar's label is in the same colour.
        let label = (0..24).find_map(|y| {
            let row: String = (0..80).map(|x| buf[(x, y)].symbol()).collect();
            row.find(&format!("restart {i}")).map(|x| (x as u16, y))
        });
        let (x, y) = label.expect("a bar per restart");
        assert_eq!(buf[(x, y)].style().fg, colour, "restart {i}'s label");
    }
    let mark: Vec<u16> = dots
        .iter()
        .filter(|d| d.2.fg == Some(trace::MARK))
        .map(|d| d.0)
        .collect();
    assert!(mark.len() >= 3, "the start's end drawn down the plot");
    assert!(
        mark.iter().all(|&x| x == mark[0]),
        "in one column: {mark:?}"
    );
    let (left, right) = (
        dots.iter().map(|d| d.0).min(),
        dots.iter().map(|d| d.0).max(),
    );
    // The trace starts at the first epoch, right of the mark.
    assert!(
        left >= Some(mark[0]) && Some(mark[0]) < right,
        "inside the plot"
    );
    let s = screen(&mut v, 80, 24);
    for want in ["ELBO", "k-means iteration, then epoch"] {
        assert!(s.contains(want), "{want:?} in\n{s}");
    }
    assert!(
        !s.contains("lag"),
        "an ELBO's level says nothing of lagging"
    );
}

/// Too few rows for a plot: the bars stay, the plot goes; the verdict
/// keeps its line.
#[test]
fn a_small_terminal_keeps_the_bars_and_drops_the_plot() {
    let (_tx, rx) = channel_u8();
    let p = finished([0.0, 0.3, -0.2]);
    p.optimised(settled()).unwrap();
    let mut v = fitting(&p, &rx);
    for (w, h) in [(80, 12), (30, 24)] {
        let buf = drawn(&mut v, w, h);
        assert!(braille(&buf).is_empty(), "{w}x{h}: no plot");
        let s = screen(&mut v, w, h);
        assert!(s.contains("restart 2"), "{w}x{h}\n{s}");
    }
    let s = screen(&mut v, 80, 14);
    assert!(s.contains("ELBO settled in 30 epochs"), "{s}");
    for (w, h) in [(10, 3), (20, 8), (200, 60)] {
        screen(&mut v, w, h);
    }
}

/// Each of many restarts has its own colour, on its label and its line.
#[test]
fn many_chains_each_have_their_own_colour() {
    let colours: Vec<_> = (0..16).map(|i| chain_style(i).fg).collect();
    for (i, a) in colours.iter().enumerate() {
        assert!(a.is_some(), "chain {i} has a colour");
        for (j, b) in colours.iter().enumerate().skip(i + 1) {
            assert_ne!(a, b, "chains {i} and {j}");
        }
    }
}
