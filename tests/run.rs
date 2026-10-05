//! `faba run` on the command line, without a terminal.

use std::process::Command;

fn faba() -> Command {
    Command::new(env!("CARGO_BIN_EXE_faba"))
}

#[test]
fn all_is_gone() {
    let out = faba().args(["all", "--help"]).output().unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(err.contains("unrecognized subcommand 'all'"), "{err}");
}

#[test]
fn batch_process_needs_its_inputs() {
    let out = faba().args(["run", "--batch-process"]).output().unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(err.contains("BAM") && err.contains("--gff"), "{err}");
}

#[test]
fn without_a_terminal_the_view_asks_for_batch_process() {
    let out = faba()
        .args(["run"])
        .stdin(std::process::Stdio::null())
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(err.contains("--batch-process"), "{err}");
}
