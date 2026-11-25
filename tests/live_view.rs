//! A view put in front of the user over MCP: `view_open` starts `faba
//! qc-view` as a detached session, a terminal attaches and the user's
//! keys arrange it, `view_state` reads back what they arranged,
//! `view_export` draws the figure the view's own `s` saves from it, and
//! the state kept beside the figure draws it again without the view.

#![cfg(all(unix, feature = "mcp"))]

use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, ChildStdout, Stdio};
use std::time::{Duration, Instant};

use legume_tui::testing::Private;
use legume_tui::Size;
use serde_json::{json, Value};

mod common;
use common::qc_fixture;

const EXE: &str = env!("CARGO_BIN_EXE_faba");

const SIZE: Size = Size {
    rows: 40,
    cols: 200,
    width: 2000,
    height: 800,
};

/// `faba mcp` on stdio, its sessions in `dir`'s folder, run in `cwd`.
struct Server {
    child: Child,
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    next: u64,
}

impl Server {
    fn start(dir: &Private, cwd: &Path) -> Self {
        let mut child = (dir.command())
            .arg("mcp")
            .current_dir(cwd)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let stdin = child.stdin.take().unwrap();
        let stdout = BufReader::new(child.stdout.take().unwrap());
        Self {
            child,
            stdin,
            stdout,
            next: 0,
        }
    }

    /// A tool's result: its first content as JSON, or its error; the
    /// whole result too.
    fn call(&mut self, name: &str, args: Value) -> (Result<Value, String>, Value) {
        self.next += 1;
        let id = self.next;
        let params = json!({ "name": name, "arguments": args });
        let msg = json!({ "jsonrpc": "2.0", "id": id, "method": "tools/call", "params": params });
        writeln!(self.stdin, "{msg}").unwrap();
        let r = loop {
            let mut line = String::new();
            assert!(
                self.stdout.read_line(&mut line).unwrap() > 0,
                "server closed"
            );
            let v: Value = serde_json::from_str(&line).unwrap();
            if v["id"] == json!(id) {
                break v["result"].clone();
            }
        };
        let first = r["content"][0]["text"].as_str().unwrap_or_default();
        let got = match r["isError"] == json!(true) {
            true => Err(first.to_string()),
            false => Ok(serde_json::from_str(first).unwrap_or_else(|_| json!(first))),
        };
        (got, r)
    }

    fn ok(&mut self, name: &str, args: Value) -> Value {
        let (got, _) = self.call(name, args.clone());
        got.unwrap_or_else(|e| panic!("{name} {args}: {e}"))
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// The view's state once `done` says it is as wanted, asked until then.
fn state_until(s: &mut Server, session: &str, done: impl Fn(&Value) -> bool) -> Value {
    let start = Instant::now();
    loop {
        let v = s.ok("view_state", json!({ "session": session }));
        if done(&v["state"]) {
            return v["state"].clone();
        }
        assert!(
            start.elapsed() < Duration::from_secs(30),
            "never arranged: {}",
            v["state"]
        );
        std::thread::sleep(Duration::from_millis(100));
    }
}

/// A PDF's streams, each inflated, in order of their bytes: the page's
/// drawing and its fonts, whatever order the writer put the fonts in.
fn pdf_streams(path: &Path) -> Vec<Vec<u8>> {
    use std::io::Read;
    let pdf = std::fs::read(path).unwrap();
    let mut streams = Vec::new();
    let mut rest = &pdf[..];
    while let Some(at) = rest.windows(7).position(|w| w == b"stream\n") {
        rest = &rest[at + 7..];
        let end = rest.windows(9).position(|w| w == b"endstream").unwrap();
        let raw = rest[..end].strip_suffix(b"\n").unwrap_or(&rest[..end]);
        let mut inflated = Vec::new();
        let stream = match flate2::read::ZlibDecoder::new(raw).read_to_end(&mut inflated) {
            Ok(_) => inflated,
            Err(_) => raw.to_vec(),
        };
        streams.push(stream);
        rest = &rest[end + 9..];
    }
    streams.sort();
    streams
}

#[test]
fn a_view_opened_over_mcp_is_arranged_by_the_user_and_drawn_as_arranged() {
    let dir = Private::new(EXE).patience(Duration::from_secs(30));
    let work = tempfile::tempdir().unwrap();
    let input = work.path().join("in");
    std::fs::create_dir_all(&input).unwrap();
    qc_fixture(&input);
    let input = input.to_string_lossy().into_owned();
    let mut s = Server::start(&dir, work.path());

    // A view that ends as it opens says why, and leaves no session.
    let bad = json!({ "view": "qc_view", "args": { "input": input, "front": "mut" } });
    let e = s.call("view_open", bad).0.unwrap_err();
    assert!(
        e.contains("ended as it opened") && e.contains("--front"),
        "{e}"
    );
    let start = Instant::now();
    while !dir.sockets().is_empty() && start.elapsed() < Duration::from_secs(10) {
        std::thread::sleep(Duration::from_millis(50));
    }
    assert!(dir.sockets().is_empty(), "{:?}", dir.sockets());

    let opened = s.ok(
        "view_open",
        json!({ "view": "qc_view", "args": { "input": input, "style": "view" }, "session": "t1" }),
    );
    let id = opened["session"].as_str().unwrap().to_string();
    assert_eq!(id, "faba-qc-view-t1");
    assert!(opened["attach"]
        .as_str()
        .unwrap()
        .ends_with(&format!("attach {id}")));
    assert_eq!(opened["state"]["view"], "qc_view");
    assert_eq!(opened["state"]["focus"], "site_max_pv");
    let listed = s.ok("view_list", json!({}));
    let one = &listed["sessions"][0];
    assert_eq!(
        (&one["session"], &one["answers"]),
        (&json!(id), &json!(true))
    );
    assert_eq!(one["attached"], false);

    // The user attaches and arranges it: a knob up (round to the last),
    // its threshold up a step, the control arm in front.
    let mut t = dir.term_sized(&["attach", &id], &[], SIZE);
    t.wait_screen("qc-view");
    // A key typed while the attach still asks the terminal what it can
    // do may be taken for its answer: `f` again until the view took it.
    let start = Instant::now();
    while s.ok("view_state", json!({ "session": id }))["state"]["front"] != "mut" {
        assert!(start.elapsed() < Duration::from_secs(30), "f never taken");
        t.type_(b"f");
        std::thread::sleep(Duration::from_secs(1));
    }
    t.type_(b"\x1b[A");
    t.type_(b"+");
    let arranged = state_until(&mut s, &id, |st| {
        st["focus"] == "site_min_cells" && st["thresholds"]["site_min_cells"] == 1.0
    });
    assert_eq!(arranged["front"], "mut");
    // `y` steps the knob's count scale, and the agent sets where the key
    // goes and the marginals' scale: `s` and the export draw them alike.
    let before = arranged["y_scale"]["cells"].clone();
    t.type_(b"y");
    let set = json!({ "key_at": "inside", "marginal_scale": "linear" });
    let after = s.ok("view_state", json!({ "session": id, "set": set }));
    assert_eq!(after["state"]["key_at"], "inside");
    let arranged = state_until(&mut s, &id, |st| st["y_scale"]["cells"] != before);
    assert_eq!(
        s.ok("view_state", json!({ "session": "t1" }))["attached"],
        true
    );

    // `s` saves the panel on screen, under its own name, in the view's
    // folder: the export of the same state is the same PDF.
    t.type_(b"s");
    t.wait_screen("save as");
    t.type_(b"\r");
    let saved = work.path().join("qc_view.pdf");
    // Read once it is written whole.
    let whole_pdf =
        |p: &Path| std::fs::read(p).is_ok_and(|b| b.trim_ascii_end().ends_with(b"%%EOF"));
    let start = Instant::now();
    while !whole_pdf(&saved) {
        assert!(
            start.elapsed() < Duration::from_secs(20),
            "`s` saved nothing"
        );
        std::thread::sleep(Duration::from_millis(50));
    }
    let prefix = work.path().join("exp").to_string_lossy().into_owned();
    let (got, whole) = s.call(
        "view_export",
        json!({ "session": id, "out": prefix, "formats": ["svg", "pdf", "png"] }),
    );
    let export = got.unwrap();
    assert_eq!(
        export["written"],
        json!([
            format!("{prefix}.svg"),
            format!("{prefix}.pdf"),
            format!("{prefix}.png")
        ])
    );
    assert_eq!(
        whole["content"][1]["type"], "image",
        "the figure as an image"
    );
    assert_eq!(
        pdf_streams(&saved),
        pdf_streams(Path::new(&format!("{prefix}.pdf"))),
        "the export is what `s` saves"
    );
    let kept: Value =
        serde_json::from_str(&std::fs::read_to_string(format!("{prefix}.view.json")).unwrap())
            .unwrap();
    for k in [
        "focus",
        "thresholds",
        "front",
        "x_scale",
        "y_scale",
        "key_at",
        "saves",
    ] {
        assert_eq!(kept[k], arranged[k], "{k}");
    }
    let again = s
        .call("view_export", json!({ "session": id, "out": prefix }))
        .0;
    assert!(again.unwrap_err().contains("not replacing"));

    // Drawn again from the kept state, without the view: the same SVG,
    // by the replay the export gave and by qc_view's `figure`.
    let replay: Vec<String> = serde_json::from_value(export["replay"].clone()).unwrap();
    let redrawn = work.path().join("again").to_string_lossy().into_owned();
    let mut argv: Vec<String> = replay[1..].to_vec();
    for w in argv.iter_mut() {
        if w.starts_with("--figure=") {
            *w = format!("--figure={redrawn}");
        }
    }
    let out = dir.command().args(&argv).output().unwrap();
    assert!(
        out.status.success(),
        "{argv:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let svg = |p: &str| std::fs::read_to_string(format!("{p}.svg")).unwrap();
    assert_eq!(svg(&redrawn), svg(&prefix), "replayed");
    let by_tool = work.path().join("tool").to_string_lossy().into_owned();
    let drawn = s.ok(
        "qc_view",
        json!({ "input": input, "style": "view", "view_state": format!("{prefix}.view.json"), "figure": by_tool, "format": ["svg"] }),
    );
    assert_eq!(drawn["figure"], json!([format!("{by_tool}.svg")]));
    assert_eq!(svg(&by_tool), svg(&prefix), "qc_view's figure");

    // The agent arranges it too: the user sees it change.
    let set = s.ok(
        "view_state",
        json!({ "session": id, "set": { "front": "wt" } }),
    );
    assert_eq!(set["state"]["front"], "wt");

    // q leaves the view and ends the session: nothing answers.
    t.type_(b"q");
    assert_eq!(t.exit_code(), 0);
    let start = Instant::now();
    loop {
        let (got, _) = s.call("view_state", json!({ "session": id }));
        match got {
            Err(e) => {
                assert!(
                    e.contains("no live session") || e.contains("no view answering"),
                    "{e}"
                );
                break;
            }
            Ok(_) => assert!(start.elapsed() < Duration::from_secs(20), "still answering"),
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    assert!(dir.sockets().is_empty());
}
