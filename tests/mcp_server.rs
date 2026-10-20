//! `faba mcp` driven over stdio, as an agent's client drives it: the tool
//! list, a read, and the refusals. A write onto an existing directory is
//! refused before any command runs, and a write the agent does not name is
//! recorded under `results/` even when the command then fails for want of
//! a BAM.

#![cfg(feature = "mcp")]

use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

const EXE: &str = env!("CARGO_BIN_EXE_faba");

struct Client {
    child: Child,
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    next: u64,
}

impl Client {
    fn start(dir: &Path) -> Self {
        let mut child = Command::new(EXE)
            .arg("mcp")
            .current_dir(dir)
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

    fn request(&mut self, method: &str, params: Value) -> Value {
        self.next += 1;
        let id = self.next;
        let msg = json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params });
        writeln!(self.stdin, "{msg}").unwrap();
        loop {
            let mut line = String::new();
            assert!(
                self.stdout.read_line(&mut line).unwrap() > 0,
                "server closed"
            );
            let v: Value = serde_json::from_str(&line).unwrap();
            if v["id"] == json!(id) {
                return v;
            }
        }
    }

    /// A tool's result: its structured content, its text, or its error.
    fn call(&mut self, name: &str, args: Value) -> Result<Value, String> {
        let r = self.request("tools/call", json!({ "name": name, "arguments": args }));
        let r = &r["result"];
        if r["isError"] == json!(true) {
            Err(r["content"][0]["text"].as_str().unwrap().to_string())
        } else if r["structuredContent"].is_object() {
            Ok(r["structuredContent"].clone())
        } else {
            Ok(r["content"][0]["text"].clone())
        }
    }
}

impl Drop for Client {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

#[test]
fn tools_read_and_refuse_to_replace() {
    let dir = tempfile::tempdir().unwrap();
    let mut c = Client::start(dir.path());

    let init = c.request("initialize", json!({ "protocolVersion": "2025-06-18" }));
    assert_eq!(init["result"]["serverInfo"]["name"], "faba");
    let tools = c.request("tools/list", json!({}));
    let names: Vec<&str> = tools["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    for tool in [
        "docs", "count", "run", "qc", "cluster", "pileup", "notebook",
    ] {
        assert!(names.contains(&tool), "{names:?}");
    }
    assert!(!names.contains(&"attach"));

    let listing = c.call("docs", json!({})).unwrap();
    assert!(listing.as_str().unwrap().contains("profiling"));
    let section = c.call("docs", json!({ "query": ["qc"] })).unwrap();
    assert!(section.as_str().unwrap().contains("qc"));
    let missing = c.call("docs", json!({ "query": ["no-such-section"] }));
    assert!(missing.unwrap_err().contains("no section"));

    std::fs::create_dir(dir.path().join("taken")).unwrap();
    let again = c.call(
        "count",
        json!({ "bam": ["no-such.bam"], "gff": "no-such.gff", "out": "taken" }),
    );
    assert!(again.unwrap_err().contains("not replacing"));
    // `out` is the only way to name an output.
    let refused = c.call(
        "count",
        json!({
            "bam": ["no-such.bam"], "gff": "no-such.gff", "out": "fresh",
            "options": ["--output=z"],
        }),
    );
    assert!(refused.unwrap_err().contains("not taken"));
    // Nor hidden in a cluster of short options.
    let refused = c.call(
        "count",
        json!({
            "bam": ["no-such.bam"], "gff": "no-such.gff", "out": "fresh",
            "options": ["-gX"],
        }),
    );
    assert!(refused.unwrap_err().contains("apart from its value"));
    assert!(!dir.path().join("fresh").exists());
}

#[test]
fn unnamed_outputs_are_recorded_and_inputs_keep_their_source() {
    let dir = tempfile::tempdir().unwrap();
    let mut c = Client::start(dir.path());
    c.request("initialize", json!({}));

    // The BAM is not there, so the command fails. The records are written
    // first: the output it would have made, and where the input came from.
    let failed = c.call(
        "count",
        json!({
            "bam": ["no-such.bam"],
            "gff": "genes.gff",
            "source": "the test's fixture",
        }),
    );
    assert!(failed.unwrap_err().contains("no-such.bam"));

    let mut recorded = String::new();
    let results = dir.path().join("results");
    for day in std::fs::read_dir(&results).unwrap() {
        let script = day.unwrap().path().join("runall.sh");
        if script.is_file() {
            recorded = std::fs::read_to_string(script).unwrap();
        }
    }
    assert!(
        recorded.contains("no-such.count") && recorded.contains("did not succeed"),
        "{recorded}"
    );
    let readme = std::fs::read_to_string(dir.path().join("data/README.md")).unwrap();
    assert!(readme.contains("the test's fixture"));
}
