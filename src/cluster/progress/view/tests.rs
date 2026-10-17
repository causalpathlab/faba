//! The fitting view drawn from a made-up progress: a bar per chain with
//! its burn-in marked, the R-hat once in, `q` cancelling.

use super::*;
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
        s.contains("split R-hat 1.01") && s.contains("the chains agree"),
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
