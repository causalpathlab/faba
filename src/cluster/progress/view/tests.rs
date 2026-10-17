//! The fitting view drawn from a made-up progress: a bar per chain with
//! its burn-in marked, the chains' log joint traced under them, the
//! R-hat once in, `q` cancelling.

use super::*;
use crate::cluster::progress::advice::lagging;
use crate::tui::screen_text as screen;
use ratatui::crossterm::event::KeyModifiers;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{channel, Receiver, Sender};

/// Three chains of 100 sweeps, 40 of them burn-in: one started 10 s ago
/// at sweep 20, one 7 s ago at sweep 70, one still queued.
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
        title: "fitting the blocks model",
        progress: p,
        pending: Pending::new(rx),
        skip: false,
        leave: false,
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

/// The burn-in in one colour up to its mark, the sampling in another after
/// it, the rest dim.
#[test]
fn a_chains_bar_marks_where_the_burn_in_ends() {
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
    assert_eq!(
        text(&chain_bar(0, 100, 0, 5)),
        "░░░░░",
        "no burn-in, no mark"
    );
    assert_eq!(text(&chain_bar(100, 100, 40, 5)), "██│██");
}

/// Each chain's bar, sweeps and phase, its time taken and left counted
/// from when it started, a queued chain as waiting, and the key to cancel;
/// at any terminal size.
#[test]
fn draws_a_bar_per_chain_with_its_burn_in() {
    let (_tx, rx) = channel_u8();
    let p = midway();
    let mut v = fitting(&p, &rx);
    let s = screen(&mut v, 110, 10);
    for want in [
        "chain 0",
        "chain 1",
        "20/100 burn-in",
        "70/100 sampling",
        "00:00:10 taken",
        "00:00:40 left",
        "00:00:07 taken",
        "00:00:03 left",
        "chain 2",
        "0/100 waiting",
        "the burn-in ends at sweep 40",
        "q cancel",
    ] {
        assert!(s.contains(want), "{want:?} in\n{s}");
    }
    // Inside the panel's borders, each chain's bar has its mark.
    for line in s.lines().filter(|l| l.contains(" chain ")) {
        let inside = line.trim_end().trim_matches('│');
        assert!(inside.contains('│'), "{line}");
    }
    for (w, h) in [(10, 3), (40, 5), (200, 60)] {
        screen(&mut v, w, h);
    }
}

/// Before the chains start, it says it is reading; `q` then waits for the
/// input in hand, and says so.
#[test]
fn says_it_reads_before_the_chains_start() {
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

/// The R-hat shows once the chains are summarised; the view ends a moment
/// after the result, or at once on a key or on an error.
#[test]
fn shows_the_r_hat_then_ends() {
    let (tx, rx) = channel_u8();
    let p = midway();
    p.set_rhat(1.01);
    let mut v = fitting(&p, &rx);
    let s = screen(&mut v, 100, 12);
    assert!(
        s.contains("chains agree (R-hat 1.01)") && s.contains("nothing to rerun"),
        "{s}"
    );
    assert!(s.contains("writing the tables"), "{s}");
    tx.send(Ok(3)).unwrap();
    v.tick();
    assert!(!v.done(), "the R-hat stays up a moment");
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

/// Ticks `view` until it is done, as the terminal loop would.
fn tick_until_done(view: &mut Fitting<u8>) -> anyhow::Result<()> {
    while !view.done() {
        view.tick();
        std::thread::sleep(Duration::from_millis(1));
    }
    Ok(())
}

/// A fit that panics ends the view and comes back as an error, not as the
/// panic again.
#[test]
fn a_fit_that_panics_is_an_error() {
    let out = std::panic::catch_unwind(|| {
        watch_with(
            "t",
            1,
            |_| -> anyhow::Result<u8> { panic!("a fit") },
            tick_until_done,
        )
    });
    let err = out.expect("no panic past the view").unwrap_err();
    assert!(err.to_string().contains("panicked: a fit"), "{err}");
    let out = watch_with("t", 1, |_| Ok(7), tick_until_done);
    assert_eq!(out.unwrap(), Some(7));
}

/// The R-hat stays up for its moment from when the result comes, however
/// long it was drawn before.
#[test]
fn the_hold_starts_when_the_result_comes() {
    let (tx, rx) = channel_u8();
    let p = midway();
    p.set_rhat(1.01);
    let mut v = fitting(&p, &rx);
    screen(&mut v, 100, 12);
    std::thread::sleep(HOLD + Duration::from_millis(50));
    tx.send(Ok(3)).unwrap();
    v.tick();
    assert!(!v.done(), "the R-hat stays up once the fit is done");
}

/// `q` or Ctrl+C once the fit writes cannot stop it, but the browser is
/// not opened after: the view ends with the result and says so.
#[test]
fn a_cancel_too_late_leaves_without_the_browser() {
    for stop in [
        |v: &mut Fitting<u8>| v.handle_key(key('q')),
        |v: &mut Fitting<u8>| v.interrupt(),
    ] {
        let (tx, rx) = channel_u8();
        let p = midway();
        p.set_rhat(1.01);
        p.commit().unwrap();
        let mut v = fitting(&p, &rx);
        stop(&mut v);
        assert!(!p.cancelled());
        let s = screen(&mut v, 100, 12);
        assert!(s.contains("not opening the modules"), "{s}");
        tx.send(Ok(3)).unwrap();
        v.tick();
        assert!(v.done(), "no hold before leaving");
        // Once the result is in, either one leaves too.
        let (tx, rx) = channel_u8();
        let mut v = fitting(&p, &rx);
        tx.send(Ok(3)).unwrap();
        v.tick();
        stop(&mut v);
        assert!(v.done() && v.leave);
    }
    let left = watch_with(
        "t",
        1,
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
    assert_eq!(left.unwrap(), None, "the browser is skipped");
}

/// A made-up trace: sweeps `from..=to`, the log joint `level` give or
/// take a fixed wobble under 1.
fn wobble(level: f64, from: usize, to: usize) -> Vec<(usize, f64)> {
    (from..=to)
        .map(|n| (n, level + ((n * 7919) % 13) as f64 / 13.0 - 0.5))
        .collect()
}

/// Three chains of 100 sweeps, 40 of them burn-in, all done: each climbs
/// from far below early in its burn-in, then sits at `levels[i]`.
fn finished(levels: [f64; 3]) -> Progress {
    let p = Progress::traced(3);
    let _ = p.schedule.set((100, 40));
    for (i, level) in levels.into_iter().enumerate() {
        let _ = p.started[i].set(Instant::now() - Duration::from_secs(5));
        p.done[i].0.store(100, Ordering::Relaxed);
        let mut t = p.traces[i].lock().unwrap();
        t.extend((1..=40).map(|n| (n, level - 1000.0 * 0.5f64.powi(n as i32))));
        t.extend(wobble(level, 41, 100));
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

/// Post-burn-in means far below the best chain's, against the chains'
/// own spread, lag: one chain of three, the lower of two, or most of them
/// stuck below one; burn-in points count for nothing; one chain, or a
/// chain with too few points, says nothing.
#[test]
fn a_chain_lags_when_far_below_the_others() {
    let low_burn_in = |level| {
        let mut t = vec![(1, -1e6), (2, -1e5)];
        t.extend(wobble(level, 3, 60));
        t
    };
    let agree = [low_burn_in(0.0), low_burn_in(0.3), low_burn_in(-0.3)];
    assert_eq!(lagging(&agree, 2), Some(Vec::new()));
    let one_low = [low_burn_in(0.0), low_burn_in(-5.0), low_burn_in(0.2)];
    assert_eq!(lagging(&one_low, 2), Some(vec![1]));
    // A wobble's sd is under 0.3: 0.7 below is within three of them.
    let near = [low_burn_in(0.0), low_burn_in(-0.5), low_burn_in(0.2)];
    assert_eq!(lagging(&near, 2), Some(Vec::new()));
    // Of two, the lower lags once three sds below, not six.
    let two = [low_burn_in(0.0), low_burn_in(-1.5)];
    assert_eq!(lagging(&two, 2), Some(vec![1]));
    // Most stuck below one: they lag, not the one above.
    let most = [low_burn_in(0.0), low_burn_in(-5.0), low_burn_in(-5.0)];
    assert_eq!(lagging(&most, 2), Some(vec![1, 2]));
    assert_eq!(lagging(&agree[..1], 2), None, "one chain");
    assert_eq!(lagging(&agree, 59), None, "one point after the burn-in");
}

/// Under the bars, a line per chain in its label's colour, the burn-in's
/// end a line down the plot; once done, whether a chain lags.
#[test]
fn traces_each_chain_under_the_bars() {
    let (_tx, rx) = channel_u8();
    let p = finished([0.0, 10.0, -10.0]);
    let mut v = fitting(&p, &rx);
    let buf = drawn(&mut v, 80, 24);
    let dots = braille(&buf);
    for i in 0..3 {
        let colour = chain_style(i).fg;
        assert!(colour.is_some(), "chain {i} has a colour");
        assert!(
            dots.iter().any(|d| d.2.fg == colour),
            "chain {i}'s trace drawn"
        );
        // The bar's label is in the same colour.
        let label = (0..24).find_map(|y| {
            let row: String = (0..80).map(|x| buf[(x, y)].symbol()).collect();
            row.find(&format!("chain {i}")).map(|x| (x as u16, y))
        });
        let (x, y) = label.expect("a bar per chain");
        assert_eq!(buf[(x, y)].style().fg, colour, "chain {i}'s label");
    }
    let mark: Vec<u16> = dots
        .iter()
        .filter(|d| d.2.fg == Some(trace::MARK))
        .map(|d| d.0)
        .collect();
    assert!(mark.len() >= 3, "the burn-in's end drawn down the plot");
    assert!(
        mark.iter().all(|&x| x == mark[0]),
        "in one column: {mark:?}"
    );
    let (left, right) = (
        dots.iter().map(|d| d.0).min(),
        dots.iter().map(|d| d.0).max(),
    );
    assert!(
        left < Some(mark[0]) && Some(mark[0]) < right,
        "inside the plot"
    );
    let s = screen(&mut v, 80, 24);
    for want in [
        "log joint",
        "sweep",
        "log joint after burn-in: chains 0, 2 lag",
    ] {
        assert!(s.contains(want), "{want:?} in\n{s}");
    }

    let p = finished([0.0, 0.3, -0.2]);
    let mut v = fitting(&p, &rx);
    let s = screen(&mut v, 80, 24);
    assert!(s.contains("log joint after burn-in: no chain lags"), "{s}");
    // While the chains run, no verdict yet.
    p.done[2].0.store(99, Ordering::Relaxed);
    assert!(!screen(&mut v, 80, 24).contains("after burn-in"));
}

/// Too few rows for a plot: the bars stay, the plot goes; the R-hat and
/// the verdict keep their lines.
#[test]
fn a_small_terminal_keeps_the_bars_and_drops_the_plot() {
    let (_tx, rx) = channel_u8();
    let p = finished([0.0, 0.3, -0.2]);
    p.set_rhat(1.01);
    let mut v = fitting(&p, &rx);
    for (w, h) in [(80, 12), (30, 24)] {
        let buf = drawn(&mut v, w, h);
        assert!(braille(&buf).is_empty(), "{w}x{h}: no plot");
        let s = screen(&mut v, w, h);
        assert!(s.contains("chain 2"), "{w}x{h}\n{s}");
    }
    let s = screen(&mut v, 80, 14);
    assert!(s.contains("chains agree (R-hat 1.01)"), "{s}");
    for (w, h) in [(10, 3), (20, 8), (200, 60)] {
        screen(&mut v, w, h);
    }
}

/// Once the R-hat is in, one line gives the verdict, the lag in it: the
/// lag line before it goes, so nothing is said twice.
#[test]
fn once_the_r_hat_is_in_one_line_gives_the_verdict() {
    let (_tx, rx) = channel_u8();
    for (levels, rhat, said) in [
        ([0.0, 0.3, -0.2], 1.3, "R-hat 1.30: chains disagree"),
        ([0.0, 0.3, -0.2], 1.01, "chains agree (R-hat 1.01)"),
        ([0.0, 10.0, -10.0], 1.01, "chains 0, 2 lag (R-hat 1.01)"),
    ] {
        let p = finished(levels);
        p.set_rhat(rhat);
        let mut v = fitting(&p, &rx);
        let s = screen(&mut v, 80, 24);
        assert!(s.contains(said), "{s}");
        assert!(!s.contains("after burn-in"), "{s}");
    }
}

/// Each of many chains has its own colour, on its label and its line.
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
