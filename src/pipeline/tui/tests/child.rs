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
