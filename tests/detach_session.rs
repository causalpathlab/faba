//! The full-screen views as detachable sessions, driven through
//! pseudo-terminals standing in for the user's: `faba pwm`'s logo (a
//! figure view, drawn as a picture where the terminal has one) starts
//! hosted, a detach leaves it running and asks for a name, an attach from
//! another terminal repaints it (its picture made again for that
//! terminal), a lost client leaves it running, and `faba kill` stops it
//! through its quit path, leaving no process or socket behind. `faba qc`'s
//! threshold view, a data-beans screen, is left by a kill as by q, a
//! changed threshold dropped without asking and nothing written. A view
//! named on the command line starts a new session while a bare start
//! offers the live ones. `--no-detach`, and a view started
//! inside a session, run plain; an error the command ends with at once
//! stays on the screen.

#![cfg(unix)]

use std::path::Path;
use std::sync::Arc;
use std::time::Duration;

use arrow::array::{ArrayRef, Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};
use legume_tui::testing::{alive, Private, Replies, Setup, Term, PROMPT};
use legume_tui::Size;
use parquet::arrow::ArrowWriter;

const EXE: &str = env!("CARGO_BIN_EXE_faba");

/// Wide enough for every footer key, the session keys at the title's end.
const SIZE: Size = Size {
    rows: 40,
    cols: 200,
    width: 2000,
    height: 800,
};

/// A private directory for faba's sessions; a session a failed test left
/// running is killed with it.
fn private() -> Private {
    Private::new(EXE).patience(Duration::from_secs(30))
}

fn term(dir: &Private, args: &[&str], env: &[(&str, &str)]) -> Term {
    dir.term_sized(args, env, SIZE)
}

/// A terminal that speaks kitty's image protocol, its cells 10 x 20
/// pixels, and faba told to ask it.
fn kitty_term(dir: &Private, args: &[&str]) -> Term {
    let setup = Setup::new()
        .size(SIZE)
        .replies(Replies::plain().kitty().cell(10, 20));
    dir.term_with(args, &[("FABA_TUI_IMAGES", "1")], setup)
}

/// A site table (`chr`, `m6a_pos`, `strand`) and a one-chromosome genome
/// for `faba pwm`: the paths of the two.
fn pwm_fixture(dir: &Path) -> (String, String) {
    let genome = dir.join("genome.fa");
    let seq = "ACGTTGCAAC".repeat(16);
    std::fs::write(&genome, format!(">chr1\n{seq}\n")).unwrap();
    let sites = dir.join("sites.parquet");
    let schema = Arc::new(Schema::new(vec![
        Field::new("chr", DataType::Utf8, false),
        Field::new("m6a_pos", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
    ]));
    let positions = [40i64, 60, 80, 100];
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"; positions.len()])) as ArrayRef,
            Arc::new(Int64Array::from(positions.to_vec())) as ArrayRef,
            Arc::new(StringArray::from(vec!["+"; positions.len()])) as ArrayRef,
        ],
    )
    .unwrap();
    let mut w = ArrowWriter::try_new(std::fs::File::create(&sites).unwrap(), schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
    let path = |p: &Path| p.to_string_lossy().into_owned();
    (path(&sites), path(&genome))
}

/// `{dir}/{name}.zarr` from named triplets.
fn write_matrix(dir: &Path, name: &str, rows: &[&str], cols: &[&str], entries: &[(u64, u64)]) {
    let path = dir.join(format!("{name}.zarr"));
    let triplets: Vec<(u64, u64, f32)> = entries.iter().map(|&(r, c)| (r, c, 3.0)).collect();
    let mut data = create_sparse_from_triplets_owned(
        triplets,
        (rows.len(), cols.len(), entries.len()),
        Some(&path.to_string_lossy()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let rows: Vec<Box<str>> = rows.iter().map(|s| (*s).into()).collect();
    let cols: Vec<Box<str>> = cols.iter().map(|s| (*s).into()).collect();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&cols);
}

/// A faba output directory in `dir` for `faba qc`: a count matrix, an m6A
/// site matrix and its site table, for one batch of four cells.
fn qc_fixture(dir: &Path) {
    let cells = ["C1", "C2", "C3", "C4"];
    let counts: Vec<(u64, u64)> = (0..2).flat_map(|r| (0..4).map(move |c| (r, c))).collect();
    write_matrix(
        dir,
        "b1_count",
        &["G1/count/spliced", "G2/count/spliced"],
        &cells,
        &counts,
    );
    let sites: Vec<(u64, u64)> = (0..4).flat_map(|r| (0..4).map(move |c| (r, c))).collect();
    write_matrix(
        dir,
        "b1_m6a_site",
        &[
            "G1/m6a/chr1:100/methylated",
            "G1/m6a/chr1:100/unmethylated",
            "G2/m6a/chr1:300/methylated",
            "G2/m6a/chr1:300/unmethylated",
        ],
        &cells,
        &sites,
    );
    let schema = Arc::new(Schema::new(vec![
        Field::new("chr", DataType::Utf8, false),
        Field::new("gene", DataType::Utf8, false),
        Field::new("primary_pos", DataType::Int64, false),
        Field::new("conversion_pos", DataType::Int64, true),
        Field::new("pv", DataType::Float32, false),
        Field::new("log_odds", DataType::Float32, true),
        Field::new("coverage", DataType::UInt64, false),
        Field::new("converted", DataType::UInt64, false),
        Field::new("control_coverage", DataType::UInt64, false),
        Field::new("control_converted", DataType::UInt64, false),
    ]));
    use arrow::array::{Float32Array, UInt64Array};
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1", "chr1"])) as ArrayRef,
            Arc::new(StringArray::from(vec!["G1", "G2"])) as ArrayRef,
            Arc::new(Int64Array::from(vec![99, 299])) as ArrayRef,
            Arc::new(Int64Array::from(vec![Some(100), Some(300)])) as ArrayRef,
            Arc::new(Float32Array::from(vec![1e-4, 1e-3])) as ArrayRef,
            Arc::new(Float32Array::from(vec![Some(3.0), Some(2.0)])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![20, 10])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![8, 3])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![20, 10])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![0, 0])) as ArrayRef,
        ],
    )
    .unwrap();
    let file = std::fs::File::create(dir.join("m6a_sites.parquet")).unwrap();
    let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
    std::fs::write(dir.join("pipeline_summary.json"), "{}").unwrap();
}

/// The files in `dir`, none when it is not there.
fn files_in(dir: &Path) -> Vec<String> {
    std::fs::read_dir(dir).map_or_else(
        |_| Vec::new(),
        |d| {
            d.flatten()
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .collect()
        },
    )
}

#[test]
fn pwm_detaches_attaches_and_is_killed_through_its_quit_path() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let (sites, genome) = pwm_fixture(work.path());
    let args = ["pwm", "-s", &sites, "-f", &genome];

    // Started in a terminal, it is hosted: this terminal is its first
    // client, and the title row ends with the session keys.
    let mut first = term(&dir, &args, &[]);
    first.wait_screen("ctrl+alt+d detach");
    first.wait_screen("ctrl+alt+s switch");
    first.wait_screen("21 positions");
    assert!(!first.cooked(), "the client's terminal is raw");
    let listed = dir.wait_listed("sequence logo of sites");
    assert!(listed.contains("faba-pwm-1\tattached"), "{listed:?}");
    let pid = dir.program_pid();
    assert_ne!(
        pid,
        first.child.id(),
        "the program runs apart from its client"
    );

    // Detached (ctrl+alt+d): the client ends, its terminal restored, the
    // program runs on; Enter leaves it unnamed.
    first.type_(b"\x1b\x04");
    first.answer(1, b"\r");
    assert_eq!(first.exit_code(), 0);
    assert!(first.cooked(), "terminal restored after detach");
    first.wait_for("Detached faba-pwm-1.");
    dir.wait_listed("faba-pwm-1\tdetached");
    assert!(alive(pid));

    // Attached from another terminal: back in the alternate screen and
    // painted whole there, keys reaching the program.
    let mut second = term(&dir, &["attach", "pwm-1"], &[]);
    second.wait_screen("21 positions · bits");
    second.wait_screen("ctrl+alt+d detach");
    assert!(second.alternate(), "back in the alternate screen");
    second.type_(b"b");
    second.wait_screen("21 positions · frequency");

    // The client lost without a word (a dropped connection): it runs on.
    second.child.kill().unwrap();
    second.child.wait().unwrap();
    drop(second);
    dir.wait_listed("faba-pwm-1\tdetached");
    assert!(alive(pid));

    // Killed while detached: it quits as q does, leaving nothing.
    let said = dir.run(&["kill", "pwm-1"]);
    assert!(said.contains("faba-pwm-1 stopped"), "{said:?}");
    dir.gone(pid);
    assert!(dir.sockets().is_empty());
    assert!(dir.run(&["sessions"]).contains("no live faba sessions"));
}

#[test]
fn a_figure_is_drawn_again_for_the_terminal_that_attaches() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let (sites, genome) = pwm_fixture(work.path());

    // A terminal with kitty's protocol: the logo is a picture.
    let mut first = kitty_term(&dir, &["pwm", "-s", &sites, "-f", &genome]);
    first.wait_screen("· image");
    first.wait_for("\x1b_G");
    dir.wait_listed("faba-pwm-1\tattached");
    let pid = dir.program_pid();
    first.type_(b"\x1b\x04");
    first.answer(1, b"\r");
    assert_eq!(first.exit_code(), 0);
    dir.wait_listed("faba-pwm-1\tdetached");

    // Another kitty terminal attaches: the attach deleted kitty's images,
    // so the picture is made and sent again.
    let mut again = kitty_term(&dir, &["attach"]);
    again.wait_screen("· image");
    again.wait_for("\x1b_G");

    // One without pictures attaches in its place: text plots there.
    let plain = term(&dir, &["attach"], &[("FABA_TUI_IMAGES", "1")]);
    plain.wait_screen("21 positions · bits");
    assert_eq!(again.exit_code(), 0, "taken over");
    let start = std::time::Instant::now();
    while plain.screen().contains("· image") && start.elapsed() < Duration::from_secs(10) {
        std::thread::sleep(Duration::from_millis(50));
    }
    assert!(!plain.screen().contains("· image"), "{}", plain.screen());

    let said = dir.run(&["kill"]);
    assert!(said.contains("faba-pwm-1 stopped"), "{said:?}");
    dir.gone(pid);
    assert!(dir.sockets().is_empty());
}

/// The row of the threshold in focus, as drawn on `screen`.
fn focused_row(screen: &str) -> String {
    screen
        .lines()
        .find(|l| l.contains("▸ "))
        .unwrap_or_default()
        .to_string()
}

#[test]
fn a_kill_leaves_the_qc_view_without_writing() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let input = work.path().join("in");
    std::fs::create_dir_all(&input).unwrap();
    qc_fixture(&input);
    let out = work.path().join("out");
    let (input_arg, out_arg) = (input.to_string_lossy(), out.to_string_lossy());

    let mut t = term(
        &dir,
        &["qc", &input_arg, "-o", &out_arg, "--no-cell-qc"],
        &[],
    );
    t.wait_screen("apply");
    t.wait_screen("ctrl+alt+d detach");
    let listed = dir.wait_listed("picking site thresholds for");
    assert!(listed.contains("faba-qc-1\tattached"), "{listed:?}");
    let pid = dir.program_pid();

    // A threshold changed and not applied: a bin along on the one in focus.
    let before = focused_row(&t.screen());
    t.type_(b"\x1b[C");
    let start = std::time::Instant::now();
    while focused_row(&t.screen()) == before && start.elapsed() < Duration::from_secs(20) {
        std::thread::sleep(Duration::from_millis(50));
    }
    let after = focused_row(&t.screen());
    assert_ne!(after, before, "the threshold moved");

    let said = dir.run(&["kill"]);
    assert!(said.contains("faba-qc-1 stopped"), "{said:?}");
    assert_eq!(t.exit_code(), 0, "the quit path, not an error");
    assert!(
        !t.text().contains("leave without writing"),
        "dropped without asking"
    );
    assert!(t.cooked(), "terminal restored after the program ended");
    dir.gone(pid);
    assert!(dir.sockets().is_empty());
    assert!(files_in(&out).is_empty(), "written: {:?}", files_in(&out));
}

#[test]
fn a_named_input_starts_a_new_session_and_a_bare_start_offers_the_live_ones() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let (sites, genome) = pwm_fixture(work.path());
    let args = ["pwm", "-s", &sites, "-f", &genome];

    let mut first = term(&dir, &args, &[]);
    first.wait_screen("21 positions");
    first.type_(b"\x1b\x04");
    first.answer(1, b"\r");
    assert_eq!(first.exit_code(), 0);
    dir.wait_listed("faba-pwm-1\tdetached");

    // Sites named: a new session, the detached one named in a hint.
    let mut named = term(&dir, &args, &[]);
    named.wait_for("[detached: faba-pwm-1");
    named.wait_screen("21 positions");
    let said = dir.run(&["kill", "pwm-2"]);
    assert!(said.contains("faba-pwm-2 stopped"), "{said:?}");
    assert_eq!(named.exit_code(), 0);

    // Nothing named: the live sessions first; Enter returns to the one.
    let mut bare = term(&dir, &["run"], &[]);
    bare.wait_screen("start new");
    bare.wait_screen("faba-pwm-1");
    bare.type_(b"\r");
    bare.wait_screen("21 positions");
    assert!(bare.alternate());
    let said = dir.run(&["kill"]);
    assert!(said.contains("faba-pwm-1 stopped"), "{said:?}");
    assert_eq!(bare.exit_code(), 0);
    assert!(dir.sockets().is_empty());
}

#[test]
fn no_detach_and_a_view_inside_a_session_run_plain() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let (sites, genome) = pwm_fixture(work.path());

    // Opted out by the flag, or started under a session's marker (as a
    // hosted program starting another view): in this process, no client
    // between, no socket, no session keys; q ends it, terminal restored.
    for (extra, env) in [
        (vec!["--no-detach"], vec![]),
        (vec![], vec![(legume_tui::SESSION_VAR, "faba-pwm-9")]),
    ] {
        let mut args = vec!["pwm", "-s", &sites, "-f", &genome];
        args.extend(&extra);
        let mut plain = term(&dir, &args, &env);
        plain.wait_screen("21 positions");
        assert!(!plain.cooked());
        assert!(dir.sockets().is_empty(), "{args:?} {env:?}");
        assert!(dir.run(&["sessions"]).contains("no live faba sessions"));
        assert!(!plain.text().contains("ctrl+alt+d"), "no session keys");
        assert!(!plain.text().contains(PROMPT));
        plain.type_(b"q");
        assert_eq!(plain.exit_code(), 0);
        assert!(plain.cooked(), "terminal restored");
    }
}

#[test]
fn an_error_the_command_ends_with_at_once_stays_on_the_screen() {
    let dir = private();
    let work = tempfile::tempdir().unwrap();
    let (_, genome) = pwm_fixture(work.path());
    let missing = work.path().join("none.parquet");
    let missing = missing.to_string_lossy();

    // Hosted, and plain: the error is said where the user looks, once,
    // and the session is gone with the program.
    for extra in [None, Some("--no-detach")] {
        let mut args = vec!["pwm", "-s", &missing, "-f", &genome];
        args.extend(extra);
        let mut t = term(&dir, &args, &[]);
        assert_eq!(t.exit_code(), 1, "{args:?}");
        t.wait_screen("Error:");
        assert!(
            t.screen().contains("No such file or directory"),
            "{}",
            t.screen()
        );
        assert!(!t.alternate() && t.cooked(), "terminal restored");
        assert_eq!(t.text().matches("Error:").count(), 1, "said once");
        assert!(dir.sockets().is_empty());
    }
}
