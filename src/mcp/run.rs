//! Running a tool's `faba` command and shaping what it wrote into the
//! call's result.

use crate::mcp::args::Arguments;
use anyhow::{bail, ensure, Context, Result};
use legume_mcp::{Call, Content};
use serde_json::{json, Map, Value};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::SystemTime;

/// `faba -v` with `argv`. `-v` so the child's log lines come back as
/// progress. The binary is this server's own path: where the file there is
/// replaced while the server runs (a `cargo install`), macOS runs the new
/// file, and Linux, its path read from `/proc/self/exe` now naming a
/// deleted file, runs the server's own image. [`Binary::check`] refuses a
/// new file that is not the server's build.
pub(super) fn command(argv: &[String]) -> Result<Command> {
    let mut args = Vec::with_capacity(argv.len() + 1);
    args.push("-v".to_string());
    args.extend(argv.iter().cloned());
    legume_mcp::own_command(args)
}

/// The executable [`command`] runs.
fn exe() -> Result<PathBuf> {
    Ok(PathBuf::from(command(&[])?.get_program()))
}

/// What tells one file from another at a path: its modified time, size
/// and inode.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Stamp {
    modified: Option<SystemTime>,
    len: u64,
    inode: u64,
}

fn stamp(path: &Path) -> Option<Stamp> {
    let m = std::fs::metadata(path).ok()?;
    #[cfg(unix)]
    let inode = std::os::unix::fs::MetadataExt::ino(&m);
    #[cfg(not(unix))]
    let inode = 0;
    Some(Stamp {
        modified: m.modified().ok(),
        len: m.len(),
        inode,
    })
}

/// The faba a call runs, checked against this server's version and, when it
/// is the server's own file, its build: a server started from another faba
/// has another tool schema, so its calls would drop or refuse flags the
/// binary has.
pub(super) struct Binary {
    version: String,
    /// The server's own executable and its stamp when the server started.
    started: Option<(PathBuf, Stamp)>,
    /// Each executable's version, by its stamp when asked.
    seen: HashMap<PathBuf, (Option<Stamp>, String)>,
}

impl Binary {
    /// A check for a server of `version`, its start unknown.
    fn new(version: &str) -> Self {
        Self {
            version: version.into(),
            started: None,
            seen: HashMap::new(),
        }
    }

    /// A check for a server of `version` started from the executable `exe`.
    fn started_from(version: &str, exe: &Path) -> Self {
        Self {
            started: stamp(exe).map(|s| (exe.to_path_buf(), s)),
            ..Self::new(version)
        }
    }

    /// This server's, started from the executable [`command`] runs.
    pub(super) fn this_server() -> Self {
        let version = env!("CARGO_PKG_VERSION");
        match exe() {
            Ok(exe) => Self::started_from(version, &exe),
            Err(_) => Self::new(version),
        }
    }

    /// Refused when the faba a call would run is not this server's version,
    /// or is the server's own file rebuilt since the server started.
    pub(super) fn check(&mut self) -> Result<()> {
        self.check_exe(&exe()?)
    }

    /// [`Binary::check`] for the executable `exe`. Its `--version` is asked
    /// once per file: again only when the file at `exe` changes.
    fn check_exe(&mut self, exe: &Path) -> Result<()> {
        let now = stamp(exe);
        let theirs = match self.seen.get(exe) {
            Some((s, v)) if *s == now => v.clone(),
            _ => {
                let v = version_of(exe)?;
                self.seen
                    .insert(exe.to_path_buf(), (now.clone(), v.clone()));
                v
            }
        };
        let since_start = match &self.started {
            Some((path, s)) => path == exe && now.as_ref() != Some(s),
            None => false,
        };
        if theirs != self.version {
            let replaced = if since_start {
                " (that file was replaced since this server started)"
            } else {
                ""
            };
            bail!(
                "this MCP server is faba {ours}, but the faba it runs, {exe}, is faba \
                 {theirs}{replaced}: the server's tool schema is out of date for it, so \
                 arguments could be dropped or refused. Nothing was run. Ask the user to \
                 restart or reconnect the faba MCP server, then call again.",
                ours = self.version,
                exe = exe.display(),
            );
        }
        // Same version, another build (a rebuild that adds flags without a
        // bump): its schema may differ all the same.
        ensure!(
            !since_start,
            "the faba this MCP server runs, {exe}, is faba {theirs} too, but it was \
             rebuilt since this server started, so the server's tool schema may be out \
             of date for it. Nothing was run. Ask the user to restart or reconnect the \
             faba MCP server, then call again.",
            exe = exe.display(),
        );
        Ok(())
    }
}

/// The version `exe --version` prints (`faba 0.20.0` is `0.20.0`).
fn version_of(exe: &Path) -> Result<String> {
    let out = Command::new(exe)
        .arg("--version")
        .output()
        .with_context(|| format!("running {} --version", exe.display()))?;
    ensure!(
        out.status.success(),
        "{} --version failed: {}",
        exe.display(),
        String::from_utf8_lossy(&out.stderr).trim()
    );
    let text = String::from_utf8_lossy(&out.stdout);
    text.split_whitespace()
        .last()
        .map(str::to_string)
        .with_context(|| format!("{} --version printed nothing", exe.display()))
}

pub(super) fn absolute(p: &str) -> String {
    std::path::absolute(p)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| p.to_string())
}

/// A file or directory a tool is about to write: refused if it exists.
pub(super) fn fresh(path: &str) -> Result<()> {
    ensure!(
        !Path::new(path).exists(),
        "{path} exists; not replacing it: choose another `out`"
    );
    Ok(())
}

pub(super) fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

/// Where an input's name comes from: the file name, without a suffix, so a
/// chain of unnamed outputs can be read (`sample.count`).
fn stem_of(path: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "faba".into());
    name.split('.').next().unwrap_or("faba").to_string()
}

pub(super) fn named(input: &str, tool: &str) -> String {
    format!("{}.{}", stem_of(input), tool)
}

/// `source`, when the agent gave one, against each input.
pub(super) fn note(cx: &Call, a: &Arguments, paths: &[String]) {
    let Some(source) = a.str("source") else {
        return;
    };
    for p in paths {
        cx.provenance(Path::new(p), source);
    }
}

/// What a write returns: the path, the files now in a directory, and
/// anything the command printed.
pub(super) fn outcome(path: &str, printed: &str) -> Value {
    let mut v = Map::new();
    v.insert("written".into(), json!(absolute(path)));
    let files = entries(path);
    if !files.is_empty() {
        v.insert("files".into(), json!(files));
    }
    let printed = printed.trim();
    if !printed.is_empty() {
        v.insert("stdout".into(), json!(printed));
    }
    if let Some(results) = recorded_results(path) {
        v.insert("results".into(), results);
    }
    Value::Object(v)
}

/// Largest run-record `results` returned inline, in bytes of JSON.
const RESULTS_INLINE: usize = 16 * 1024;

/// The `results` of every run record in the directory `dir`, by job: what a
/// run reports beyond its files (a fit's verdict, its cluster sizes), which
/// the child logs to stderr and the call would otherwise only show as
/// progress.
fn recorded_results(dir: &str) -> Option<Value> {
    let mut m = Map::new();
    for e in std::fs::read_dir(dir).ok()?.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        let Some(job) = name.strip_suffix(".run.json") else {
            continue;
        };
        let Some(results) = std::fs::read_to_string(e.path())
            .ok()
            .and_then(|t| serde_json::from_str::<Value>(&t).ok())
            .and_then(|mut v| v.get_mut("results").map(Value::take))
            .filter(|r| !r.is_null())
        else {
            continue;
        };
        let size = results.to_string().len();
        let results = if size <= RESULTS_INLINE {
            results
        } else {
            json!(format!(
                "{size} bytes: read `results` in {}",
                absolute(&e.path().to_string_lossy())
            ))
        };
        m.insert(job.into(), results);
    }
    (!m.is_empty()).then_some(Value::Object(m))
}

/// The prefix `--figure` takes beside the table `table`, after checking the
/// PDF and PNG it writes are new; all are added to `outputs`.
pub(super) fn figure_prefix(table: &str, outputs: &mut Vec<PathBuf>) -> Result<(String, String)> {
    let prefix = table.strip_suffix(".tsv").unwrap_or(table);
    // Exactly the files a figure saved under `prefix` writes.
    let paths = crate::figure::saved_paths(prefix);
    for path in &paths {
        fresh(path)?;
    }
    outputs.extend(paths.iter().map(PathBuf::from));
    let [_, png] = paths;
    Ok((prefix.to_string(), png))
}

/// `v`, with the figure `png` (and its PDF beside it) named in it and
/// attached as an image.
pub(super) fn with_figure(mut v: Value, png: Option<&str>) -> Result<Vec<Content>> {
    let Some(png) = png else {
        return Ok(vec![Content::Json(v)]);
    };
    let pdf = format!("{}.pdf", png.strip_suffix(".png").unwrap_or(png));
    if let Some(m) = v.as_object_mut() {
        m.insert("figure".into(), json!([absolute(&pdf), absolute(png)]));
    }
    let image = std::fs::read(png).with_context(|| format!("reading {png}"))?;
    Ok(vec![Content::Json(v), Content::Png(image)])
}

/// Widest image returned for an SVG, in pixels.
const IMAGE_WIDTH: f64 = 1600.0;

/// The SVG at `path` as a PNG at most [`IMAGE_WIDTH`] pixels wide.
pub(super) fn svg_png(path: &str) -> Result<Vec<u8>> {
    let svg = std::fs::read_to_string(path)?;
    let (w, _) = crate::figure::svg_size(&svg);
    let scale = (IMAGE_WIDTH / w.max(1.0)).min(2.0) as f32;
    crate::figure::png(&svg, scale)
}

/// Most `qc_plots/` images attached to one `qc` call.
const QC_IMAGES: usize = 8;

/// A TSV with a header line as one object per row.
pub(super) fn tsv_rows(text: &str) -> Value {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return json!([]);
    };
    let keys: Vec<&str> = header.split('\t').collect();
    let rows: Vec<Value> = lines
        .filter(|l| !l.is_empty())
        .map(|l| {
            let m: Map<String, Value> = keys
                .iter()
                .zip(l.split('\t'))
                .map(|(k, x)| {
                    let x = x.parse::<i64>().map_or_else(|_| json!(x), |n| json!(n));
                    (k.to_string(), x)
                })
                .collect();
            Value::Object(m)
        })
        .collect();
    json!(rows)
}

/// [`outcome`] for a text file, with its first lines, so the table can be
/// read without a second call.
pub(super) fn text_outcome(path: &str, printed: &str) -> Value {
    let mut v = outcome(path, printed);
    let Ok(text) = std::fs::read_to_string(path) else {
        return v;
    };
    let lines: Vec<&str> = text.lines().collect();
    let head = lines
        .iter()
        .take(40)
        .copied()
        .collect::<Vec<_>>()
        .join("\n");
    let Some(m) = v.as_object_mut() else {
        return v;
    };
    m.insert("lines".into(), json!(lines.len()));
    m.insert("head".into(), json!(head));
    v
}

/// The files a write left at `path`: the file itself, or a directory's
/// immediate children. A `.zarr` tree is one entry, not every chunk.
fn entries(path: &str) -> Vec<String> {
    let p = Path::new(path);
    if p.is_file() {
        return vec![absolute(path)];
    }
    let mut v: Vec<String> = std::fs::read_dir(p)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| absolute(&e.path().to_string_lossy()))
        .collect();
    v.sort();
    v
}

/// Run `argv` so that it writes `out` (made of `what`) with `-o`; with
/// `background`, as a job whose id comes back at once.
pub(super) fn run_out(
    a: &Arguments,
    cx: &Call,
    mut argv: Vec<String>,
    out: &str,
    what: &str,
) -> Result<Value> {
    fresh(out)?;
    argv.extend(["-o".into(), out.into()]);
    if a.flag("background")? {
        let job = cx.spawn_job(command(&argv)?, &[PathBuf::from(out)], what)?;
        return Ok(json!({
            "job": job,
            "written": absolute(out),
            "note": "running as a background job: job_wait or job_status until it is done, then use `written`",
        }));
    }
    let printed = cx.run_making(command(&argv)?, &[PathBuf::from(out)], what)?;
    Ok(outcome(out, &printed))
}

/// [`run_out`] for a text file, returned with its first lines.
pub(super) fn run_text(cx: &Call, mut argv: Vec<String>, out: &str, what: &str) -> Result<Value> {
    fresh(out)?;
    argv.extend(["-o".into(), out.into()]);
    let printed = cx.run_making(command(&argv)?, &[PathBuf::from(out)], what)?;
    Ok(text_outcome(out, &printed))
}

/// The PNGs directly in `dir`, in name order; none when it cannot be read.
pub(super) fn pngs(dir: &Path) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = (std::fs::read_dir(dir).into_iter().flatten().flatten())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "png"))
        .collect();
    v.sort();
    v
}

/// The figures under `out`'s `qc_plots/`, as images: the view's, then
/// the panels under `qc_plots/panels/`, each in name order, at most
/// [`QC_IMAGES`] in all; half that each, the rest to whichever has more.
pub(super) fn qc_images(out: &str) -> Result<Vec<Content>> {
    let plots = Path::new(out).join("qc_plots");
    let (figures, panels) = (pngs(&plots), pngs(&plots.join("panels")));
    let half = QC_IMAGES / 2;
    let n_panels = panels.len().min(QC_IMAGES - figures.len().min(half));
    let n_figures = figures.len().min(QC_IMAGES - n_panels);
    (figures.iter().take(n_figures))
        .chain(panels.iter().take(n_panels))
        .map(|p| Ok(Content::Png(std::fs::read(p)?)))
        .collect()
}

#[cfg(test)]
mod tests {

    use super::{qc_images, recorded_results, tsv_rows, Binary, QC_IMAGES, RESULTS_INLINE};
    use legume_mcp::Content;
    use serde_json::json;
    use std::path::{Path, PathBuf};

    /// An executable at `dir/name` that prints `faba {version}` for
    /// `--version` and counts its runs in `dir/runs`.
    #[cfg(unix)]
    fn fake_faba(dir: &Path, name: &str, version: &str) -> PathBuf {
        use std::os::unix::fs::PermissionsExt;
        let path = dir.join(name);
        let runs = dir.join("runs");
        let script = format!(
            "#!/bin/sh\necho run >> '{}'\necho 'faba {version}'\n",
            runs.display()
        );
        std::fs::write(&path, script).unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
        path
    }

    #[cfg(unix)]
    fn runs(dir: &Path) -> usize {
        std::fs::read_to_string(dir.join("runs")).map_or(0, |t| t.lines().count())
    }

    /// A faba of another version, older or newer, is refused, naming both
    /// versions and what to do; the same version runs.
    #[cfg(unix)]
    #[test]
    fn a_faba_of_another_version_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        for theirs in ["0.1.0", "99.0.0"] {
            let exe = fake_faba(dir.path(), &format!("faba-{theirs}"), theirs);
            let err = Binary::new("0.20.0")
                .check_exe(&exe)
                .unwrap_err()
                .to_string();
            assert!(err.contains("0.20.0"), "{err}");
            assert!(err.contains(theirs), "{err}");
            assert!(err.contains("restart or reconnect"), "{err}");
            assert!(!err.contains("replaced"), "{err}");
        }
        let exe = fake_faba(dir.path(), "faba-same", "0.20.0");
        Binary::new("0.20.0").check_exe(&exe).unwrap();
    }

    /// `--version` runs once per executable until the file changes.
    #[cfg(unix)]
    #[test]
    fn the_version_is_asked_once_per_file() {
        let dir = tempfile::tempdir().unwrap();
        let exe = fake_faba(dir.path(), "faba", "0.20.0");
        let mut b = Binary::new("0.20.0");
        for _ in 0..3 {
            b.check_exe(&exe).unwrap();
        }
        assert_eq!(runs(dir.path()), 1);
        // Replaced (a new file renamed over it): asked again, and refused.
        let next = fake_faba(dir.path(), "next", "0.21.0");
        std::fs::rename(&next, &exe).unwrap();
        assert!(b.check_exe(&exe).is_err());
        assert!(b.check_exe(&exe).is_err());
        assert_eq!(runs(dir.path()), 2);
    }

    /// The server's own executable replaced since it started is said so.
    #[cfg(unix)]
    #[test]
    fn a_replaced_own_executable_is_named() {
        let dir = tempfile::tempdir().unwrap();
        let exe = fake_faba(dir.path(), "faba", "0.20.0");
        let mut b = Binary::started_from("0.20.0", &exe);
        b.check_exe(&exe).unwrap();
        let next = fake_faba(dir.path(), "next", "0.21.0");
        std::fs::rename(&next, &exe).unwrap();
        let err = b.check_exe(&exe).unwrap_err().to_string();
        assert!(err.contains("replaced since this server started"), "{err}");
        assert!(err.contains("0.21.0") && err.contains("0.20.0"), "{err}");
    }

    /// The server's own executable rebuilt at the same version since it
    /// started is refused as rebuilt; a server that did not start from it
    /// runs it.
    #[cfg(unix)]
    #[test]
    fn a_rebuild_of_the_same_version_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let exe = fake_faba(dir.path(), "faba", "0.20.0");
        let mut b = Binary::started_from("0.20.0", &exe);
        b.check_exe(&exe).unwrap();
        // Another build of 0.20.0: another file, another size.
        let next = dir.path().join("next");
        std::fs::write(&next, "#!/bin/sh\n# another build\necho 'faba 0.20.0'\n").unwrap();
        std::fs::set_permissions(&next, std::fs::metadata(&exe).unwrap().permissions()).unwrap();
        std::fs::rename(&next, &exe).unwrap();
        let err = b.check_exe(&exe).unwrap_err().to_string();
        assert!(err.contains("rebuilt since this server started"), "{err}");
        assert!(err.contains("restart or reconnect"), "{err}");
        assert!(!err.contains("replaced"), "{err}");
        Binary::new("0.20.0").check_exe(&exe).unwrap();
    }
    #[test]
    fn a_dry_run_table_reads_as_rows() {
        let rows = tsv_rows(
            "level	name	removes	only_this
cells	b1	3	
sites	m6A	10	2
",
        );
        assert_eq!(
            rows,
            json!([
                { "level": "cells", "name": "b1", "removes": 3, "only_this": "" },
                { "level": "sites", "name": "m6A", "removes": 10, "only_this": 2 },
            ])
        );
    }

    /// A run's results come back with its files, by job; a record too
    /// large to quote is pointed to.
    #[test]
    fn run_record_results_come_back() {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        std::fs::write(
            d.join("cluster.run.json"),
            json!({ "job": "cluster", "results": { "blockmodel": { "rhat_log_joint": 1.01 } } })
                .to_string(),
        )
        .unwrap();
        std::fs::write(d.join("qc.run.json"), json!({ "job": "qc" }).to_string()).unwrap();
        let big = "x".repeat(RESULTS_INLINE + 1);
        std::fs::write(
            d.join("big.run.json"),
            json!({ "results": big }).to_string(),
        )
        .unwrap();
        let r = recorded_results(&d.to_string_lossy()).unwrap();
        assert_eq!(r["cluster"]["blockmodel"]["rhat_log_joint"], json!(1.01));
        assert!(r.get("qc").is_none());
        assert!(r["big"].as_str().unwrap().contains("big.run.json"));
    }

    /// `qc`'s images are its view's figures and its panels under
    /// `qc_plots/panels/`, each in name order, half the room each and the
    /// rest to whichever has more.
    #[test]
    fn qc_returns_its_panels_with_its_figures() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().to_string_lossy().into_owned();
        let plots = dir.path().join("qc_plots");
        std::fs::create_dir_all(plots.join("panels")).unwrap();
        let png = |path: PathBuf, n: usize| std::fs::write(path, [n as u8]).unwrap();
        let figures = |n: usize| {
            for k in 0..n {
                png(plots.join(format!("f{k:02}.png")), k);
                std::fs::write(plots.join(format!("f{k:02}.pdf")), "").unwrap();
            }
        };
        let bytes = |c: &Content| match c {
            Content::Png(b) => b[0],
            _ => panic!("not an image"),
        };
        // Two panels: both come back, after the view's figures.
        figures(QC_IMAGES + 3);
        for k in 0..2 {
            png(plots.join("panels").join(format!("p{k:02}.png")), 100 + k);
        }
        let got: Vec<u8> = qc_images(&out).unwrap().iter().map(bytes).collect();
        let want: Vec<u8> = (0..QC_IMAGES as u8 - 2).chain([100, 101]).collect();
        assert_eq!(got, want);
        // Many panels: half the room each.
        for k in 2..QC_IMAGES + 2 {
            png(plots.join("panels").join(format!("p{k:02}.png")), 100 + k);
        }
        let got: Vec<u8> = qc_images(&out).unwrap().iter().map(bytes).collect();
        let half = QC_IMAGES as u8 / 2;
        let want: Vec<u8> = (0..half).chain(100..100 + half).collect();
        assert_eq!(got, want);
    }
}
