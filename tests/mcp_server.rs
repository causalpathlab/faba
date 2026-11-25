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
        "docs",
        "count",
        "run",
        "qc",
        "site_plot",
        "cluster",
        "pileup",
    ] {
        assert!(names.contains(&tool), "{names:?}");
    }
    // The project's notebook is the `legume` server's alone (legume-mcp 0.6).
    for tool in ["attach", "notebook"] {
        assert!(!names.contains(&tool), "{names:?}");
    }

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

#[test]
fn a_producer_runs_as_a_background_job() {
    let dir = tempfile::tempdir().unwrap();
    let mut c = Client::start(dir.path());
    c.request("initialize", json!({}));
    let listed = c.request("tools/list", json!({}));
    let tools = listed["result"]["tools"].as_array().unwrap();
    let count = tools.iter().find(|t| t["name"] == "count").unwrap();
    assert_eq!(
        count["inputSchema"]["properties"]["background"]["type"],
        "boolean"
    );
    for builtin in ["job_status", "job_wait", "job_cancel"] {
        assert!(tools.iter().any(|t| t["name"] == builtin), "{builtin}");
    }

    // It returns at once with a job; waiting tells how it went.
    let started = c
        .call(
            "count",
            json!({ "bam": ["no-such.bam"], "gff": "genes.gff", "background": true }),
        )
        .unwrap();
    let job = started["job"].as_str().unwrap().to_string();
    assert!(started["written"]
        .as_str()
        .unwrap()
        .contains("no-such.count"));
    let ended = c
        .call("job_wait", json!({ "id": job, "seconds": 60 }))
        .unwrap();
    assert_eq!(ended["state"], "failed", "{ended}");
    assert!(
        ended["error"].as_str().unwrap().contains("no-such.bam"),
        "{ended}"
    );
}

#[test]
fn the_readme_gives_an_agent_the_servers_own_instructions() {
    // For clients that do not pass the MCP instructions on to their model:
    // the README's block is faba's own text, before the notes legume-mcp
    // adds (records, jobs).
    let dir = tempfile::tempdir().unwrap();
    let mut c = Client::start(dir.path());
    let init = c.request("initialize", json!({}));
    let told = init["result"]["instructions"].as_str().unwrap();
    let own = told
        .split("\n\nRecords (after Noble")
        .next()
        .unwrap()
        .split("\n\nA tool may run a long step")
        .next()
        .unwrap();
    let readme = include_str!("../README.md");
    assert!(
        readme.contains(own),
        "README's agent block differs from faba's instructions"
    );
}

/// A count fit of two kinds of cells in `dir`, fitted by the built binary:
/// 40 cells, half of 20 genes high in each kind.
fn count_fit(dir: &Path) -> String {
    use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};
    let (cells, genes) = (40usize, 20usize);
    let names = |n: usize, f: &dyn Fn(usize) -> String| -> Vec<Box<str>> {
        (0..n).map(|i| f(i).into_boxed_str()).collect()
    };
    let mut triplets = Vec::new();
    for c in 0..cells {
        for g in 0..genes {
            let high = (c < cells / 2) == (g < genes / 2);
            let umis = if high { 8 } else { 1 } + (c * 7 + g * 3) % 3;
            triplets.push((g as u64, c as u64, umis as f32));
        }
    }
    let path = dir.join("b1_count.zarr");
    let n = triplets.len();
    let mut data = create_sparse_from_triplets_owned(
        triplets,
        (genes, cells, n),
        Some(&path.to_string_lossy()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    data.register_row_names_vec(&names(genes, &|g| format!("G{g}/count/spliced")));
    data.register_column_names_vec(&names(cells, &|c| format!("CELL{c:02}-1")));
    let folder = format!("{}/", dir.display());
    let line = ["cluster", "--modality", "count", &folder, "--batch-process"];
    let more = ["--chains", "2", "--epochs", "20", "--min-cell-sites", "1"];
    let out = Command::new(EXE)
        .args(line.iter().chain(&more))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    dir.join("count_blockmodel.run.json")
        .to_string_lossy()
        .into_owned()
}

/// A PNG's width and height in pixels, from its header.
fn png_size(png: &[u8]) -> (u32, u32) {
    let at = |i: usize| u32::from_be_bytes(png[i..i + 4].try_into().unwrap());
    (at(16), at(20))
}

/// `cluster_view` over the server: the overview by default, as JSON, its
/// GO left unread (a GAF not there is never opened); the levels as JSON
/// rows; GO read only with `go`; and several figures in one call, each an
/// image at the size asked for.
#[test]
fn cluster_view_returns_json_levels_and_several_figures() {
    let dir = tempfile::tempdir().unwrap();
    let record = count_fit(dir.path());
    let mut c = Client::start(dir.path());
    c.request("initialize", json!({}));
    let missing = json!({ "data": record, "gaf": "no-such.gaf", "obo": "no-such.obo" });
    let overview = c.call("cluster_view", missing.clone()).unwrap();
    let o = &overview["overview"];
    assert!(
        o["verdict"] == "ELBO ok" || o["verdict"] == "ELBO check",
        "{overview:#}"
    );
    assert!(o["gene_modules"].as_u64().unwrap() >= 1, "{o:#}");
    assert!(o["cell_sizes"].as_array().is_some_and(|c| !c.is_empty()));
    let genes = o["gene_sizes"].as_array().unwrap();
    assert!(genes
        .iter()
        .any(|g| !g["top"].as_array().unwrap().is_empty()));
    assert!(overview.get("go").is_none(), "{overview:#}");

    let mut asked = missing.clone();
    asked["go"] = json!(true);
    let err = c.call("cluster_view", asked).unwrap_err();
    assert!(err.contains("no-such"), "{err}");

    let levels = c
        .call(
            "cluster_view",
            json!({ "data": record, "grid": true, "cells": true, "module": "M1",
                    "trace": true, "settings": true }),
        )
        .unwrap();
    for level in ["modules", "cells", "module", "trace", "settings"] {
        assert!(levels[level].is_object(), "no {level}: {levels:#}");
    }
    let row = &levels["modules"]["rows"][0];
    assert!(row["values"][0]["value"].is_number(), "{row:#}");

    let prefix = dir.path().join("fig").to_string_lossy().into_owned();
    let r = c.request(
        "tools/call",
        json!({ "name": "cluster_view", "arguments": {
            "data": record, "figure": prefix, "plots": ["modules", "trace"],
            "width": 4, "height": 3, "font_size": 7, "panel_only": true,
        } }),
    );
    let result = &r["result"];
    assert_ne!(result["isError"], json!(true), "{r:#}");
    let images = (result["content"].as_array().unwrap().iter())
        .filter(|c| c["type"] == "image")
        .count();
    assert_eq!(images, 2, "{result:#}");
    let figures = result["structuredContent"]["figures"].as_array().unwrap();
    assert_eq!(figures.len(), 2);
    for (f, plot) in figures.iter().zip(["modules", "trace"]) {
        assert_eq!(f["plot"], plot);
        let png = std::fs::read(f["png"].as_str().unwrap()).unwrap();
        // 4 by 3 inches, two pixels a point.
        assert_eq!(png_size(&png), (576, 432), "{plot}");
        assert!(Path::new(f["pdf"].as_str().unwrap()).is_file());
    }
}

/// A within record in `dir` over the count fit `count`: group C1's fit is
/// that count fit itself, group C2 skipped.
fn within_record(dir: &Path, count: &str) -> String {
    let record = json!({ "job": "cluster", "results": { "within": {
        "modality": "m6a", "count_record": count,
        "groups": [
            { "group": "C1", "cells": 20, "status": "fitted", "dir": ".",
              "record": "count_blockmodel.run.json" },
            { "group": "C2", "cells": 3, "status": "skipped", "why": "too few" },
        ],
    } } });
    let path = dir.join("m6a_within.run.json");
    std::fs::write(&path, record.to_string()).unwrap();
    path.to_string_lossy().into_owned()
}

/// `cluster_view` of a within record with no group and no level is its
/// listing as JSON: the count fit's modules, each group's fit, the modules
/// within each fit. A level named is the count fit's; a group opens that
/// fit alone, its overview when no level is named, as a single fit does.
#[test]
fn cluster_view_lists_a_within_record_unless_a_group_or_level_is_named() {
    let dir = tempfile::tempdir().unwrap();
    let count = count_fit(dir.path());
    let record = within_record(dir.path(), &count);
    let mut c = Client::start(dir.path());
    c.request("initialize", json!({}));

    let listing = c.call("cluster_view", json!({ "data": record })).unwrap();
    assert!(listing.get("overview").is_none(), "{listing:#}");
    assert_eq!(listing["modality"], "m6a", "{listing:#}");
    assert!(listing["modules"]["rows"].is_array(), "{listing:#}");
    let groups = listing["groups"].as_array().unwrap();
    let names: Vec<&str> = groups
        .iter()
        .map(|g| g["group"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["C1", "C2"], "{listing:#}");
    assert_eq!(groups[1]["status"], "skipped");
    let within = listing["within"].as_array().unwrap();
    assert_eq!(within.len(), 1, "{listing:#}");
    assert_eq!(within[0]["group"], "C1");
    assert!(within[0]["modules"]["rows"].is_array(), "{listing:#}");

    let level = c
        .call("cluster_view", json!({ "data": record, "cells": true }))
        .unwrap();
    assert!(level["cells"].is_object(), "{level:#}");
    assert!(level.get("groups").is_none(), "{level:#}");

    for one in [
        json!({ "data": record, "group": "C1" }),
        json!({ "data": count }),
    ] {
        let overview = c.call("cluster_view", one.clone()).unwrap();
        assert!(overview["overview"].is_object(), "{one}: {overview:#}");
        assert!(overview.get("groups").is_none(), "{one}: {overview:#}");
    }
}

/// `plots` without `figure` writes under today's results folder, named
/// after the fit, the next free name once any plot's file is there; no
/// file there is ever replaced, and every figure comes back as an image.
#[test]
fn cluster_view_plots_without_a_figure_go_to_the_results_folder() {
    let dir = tempfile::tempdir().unwrap();
    let record = count_fit(dir.path());
    let mut c = Client::start(dir.path());
    c.request("initialize", json!({}));
    let asked = json!({ "data": record, "plots": ["modules", "trace"] });
    let mut call = || {
        let r = c.request(
            "tools/call",
            json!({ "name": "cluster_view", "arguments": asked }),
        );
        let result = r["result"].clone();
        assert_ne!(result["isError"], json!(true), "{r:#}");
        let images = (result["content"].as_array().unwrap().iter())
            .filter(|c| c["type"] == "image")
            .count();
        assert_eq!(images, 2, "{result:#}");
        let figures = result["structuredContent"]["figures"].as_array().unwrap();
        let files: Vec<std::path::PathBuf> = (figures.iter())
            .flat_map(|f| [&f["pdf"], &f["png"]])
            .map(|p| p.as_str().unwrap().into())
            .collect();
        assert_eq!(files.len(), 4, "{result:#}");
        files
    };
    let name = |p: &std::path::PathBuf| p.file_name().unwrap().to_string_lossy().into_owned();
    let first = call();
    let day = first[0].parent().unwrap().canonicalize().unwrap();
    let results = dir.path().join("results").canonicalize().unwrap();
    assert_eq!(day.parent().unwrap(), results, "{}", day.display());
    let names: Vec<String> = first.iter().map(name).collect();
    assert_eq!(
        names,
        [
            "count_blockmodel.cluster_view.modules.pdf",
            "count_blockmodel.cluster_view.modules.png",
            "count_blockmodel.cluster_view.trace.pdf",
            "count_blockmodel.cluster_view.trace.png",
        ]
    );
    let kept: Vec<Vec<u8>> = first.iter().map(|p| std::fs::read(p).unwrap()).collect();

    // Again: the next name, the first figures as they were.
    let second = call();
    assert!(second
        .iter()
        .all(|p| name(p).starts_with("count_blockmodel.cluster_view-2.")));
    for (p, bytes) in first.iter().zip(&kept) {
        assert_eq!(&std::fs::read(p).unwrap(), bytes, "{}", p.display());
    }

    // One plot's file under the next name is enough to take it.
    let planted = day.join("count_blockmodel.cluster_view-3.trace.png");
    std::fs::write(&planted, b"kept").unwrap();
    let third = call();
    assert!(third
        .iter()
        .all(|p| name(p).starts_with("count_blockmodel.cluster_view-4.")));
    assert_eq!(std::fs::read(&planted).unwrap(), b"kept");
}
