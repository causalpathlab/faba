use super::*;
use std::process::Command;

fn sh(script: &str) -> Command {
    let mut c = Command::new("sh");
    c.args(["-c", script]);
    c
}

#[test]
fn lines_and_exit_code_are_captured() {
    let stopper = Stopper::default();
    let mut said = Vec::new();
    let r = run_one(
        sh("echo one >&2; printf 'two\\r' >&2; echo three >&2; exit 3"),
        &stopper,
        |s| said.push(s),
    );
    assert_eq!(r, Err(Failed::Exit("three".into())));
    let lines: Vec<String> = said
        .iter()
        .filter_map(|s| match s {
            Said::Line(l) => Some(l.clone()),
            _ => None,
        })
        .collect();
    assert_eq!(lines, ["one", "two", "three"]);
}

#[test]
fn bar_frames_are_progress() {
    let p = progress_of("##########---------- 12/40 (3s) genes").unwrap();
    assert_eq!((p.pos, p.len, p.what.as_str()), (12, 40, "genes"));
    assert_eq!(
        log_line("[00:00:05] ###--- 3/10 (1s) blocks"),
        "###--- 3/10 (1s) blocks"
    );
    assert!(progress_of("Step 1/5: gene counting").is_none());
}

#[test]
fn stop_interrupts_then_kills() {
    let stopper = std::sync::Arc::new(Stopper::default());
    let s2 = stopper.clone();
    let t = std::thread::spawn(move || run_one(sh("trap '' INT; exec sleep 30"), &s2, |_| {}));
    std::thread::sleep(std::time::Duration::from_millis(300));
    stopper.stop(); // ignored by the trap
    std::thread::sleep(std::time::Duration::from_millis(200));
    stopper.stop(); // kills
    assert_eq!(t.join().unwrap(), Err(Failed::Stopped));
}

#[test]
fn stop_twice_is_harmless_after_the_end() {
    let stopper = Stopper::default();
    assert_eq!(run_one(sh("true"), &stopper, |_| {}), Ok(()));
    stopper.stop();
    stopper.stop();
}

/// A child that ignores interrupts, and the file it makes once it does.
fn deaf(dir: &std::path::Path) -> (Command, std::path::PathBuf) {
    let ready = dir.join("ready");
    let c = sh(&format!(
        "trap '' INT; : > '{}'; exec sleep 30",
        ready.display()
    ));
    (c, ready)
}

fn wait_for(f: impl Fn() -> bool) {
    let t = std::time::Instant::now();
    while !f() {
        assert!(t.elapsed().as_secs() < 20, "timed out");
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}

#[test]
fn one_early_stop_interrupts_and_two_kill() {
    let tmp = tempfile::tempdir().unwrap();
    let (mut c, ready) = deaf(tmp.path());
    let child = c.stdin(Stdio::null()).spawn().unwrap();
    wait_for(|| ready.exists());
    let stopper = Stopper::default();
    stopper.stop(); // before the child is registered
    stopper.register(child);
    // The interrupt is ignored: still running a while after.
    std::thread::sleep(std::time::Duration::from_millis(300));
    let running = |s: &Stopper| {
        let mut c = s.child.lock().unwrap();
        c.as_mut().unwrap().try_wait().unwrap().is_none()
    };
    assert!(running(&stopper), "one stop must only interrupt");
    stopper.stop();
    wait_for(|| !running(&stopper));
}

#[test]
fn two_early_stops_kill_on_registration() {
    let tmp = tempfile::tempdir().unwrap();
    let (mut c, ready) = deaf(tmp.path());
    let child = c.stdin(Stdio::null()).spawn().unwrap();
    wait_for(|| ready.exists());
    let stopper = Stopper::default();
    stopper.stop();
    stopper.stop();
    stopper.register(child);
    wait_for(|| {
        let mut c = stopper.child.lock().unwrap();
        c.as_mut().unwrap().try_wait().unwrap().is_some()
    });
}

#[test]
fn a_good_exit_with_a_stop_pending_is_ok() {
    let stopper = Stopper::default();
    let r = run_one(
        sh("trap '' INT; echo ready >&2; sleep 0.3; exit 0"),
        &stopper,
        |s| {
            if s == Said::Line("ready".into()) {
                stopper.stop();
            }
        },
    );
    assert!(stopper.is_stopped());
    assert_eq!(r, Ok(()));
}
