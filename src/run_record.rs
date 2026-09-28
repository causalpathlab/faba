//! `{job}.run.json`: what one run read, what it wrote, and with which options.
//!
//! Every producer (and every step of `faba all`) leaves one of these in its
//! output directory. Inputs are recorded as absolute paths, so the tools that
//! read an output directory later (`pileup`, `metagene`, `pwm`) can find the
//! annotation and genome it was made from without being told again; see
//! [`explicit_or_recorded`].
//!
//! Outputs are not listed by hand: the directory is snapshotted when the run
//! starts, and every top-level entry that is new or rewritten when it ends is
//! the run's output. A new output file cannot be left out of the record.

use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};

use anyhow::Context;
use rustc_hash::FxHashMap;
use serde::Serialize;
use serde_json::{json, Map, Value};

/// Serialize a field by its `Debug` form, for enums that carry no `Serialize`.
pub fn ser_debug<T: std::fmt::Debug, S: serde::Serializer>(v: &T, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_str(&format!("{v:?}"))
}

/// A run in progress: started by [`RunRecord::start`], written by
/// [`RunRecord::finish`].
pub struct RunRecord {
    job: String,
    /// `{job}.run.json` unless [`RunRecord::file_name`] says otherwise.
    file_name: String,
    dir: PathBuf,
    started: SystemTime,
    clock: Instant,
    /// Top-level entries of `dir` and their modification times at the start.
    before: FxHashMap<String, Option<SystemTime>>,
    inputs: Map<String, Value>,
    options: Value,
}

impl RunRecord {
    /// Snapshot `output_dir` before `job` writes to it.
    pub fn start(job: &str, output_dir: &str) -> Self {
        let dir = PathBuf::from(output_dir);
        Self {
            job: job.to_string(),
            file_name: format!("{job}.run.json"),
            before: snapshot(&dir),
            dir,
            started: SystemTime::now(),
            clock: Instant::now(),
            inputs: Map::new(),
            options: Value::Null,
        }
    }

    /// Write the record under `name` in place of `{job}.run.json`.
    pub fn file_name(mut self, name: &str) -> Self {
        self.file_name = name.to_string();
        self
    }

    /// Record one input file under `key` (skipped when `None`).
    pub fn input(mut self, key: &str, path: Option<&str>) -> Self {
        if let Some(p) = path {
            self.inputs.insert(key.into(), entry(p));
        }
        self
    }

    /// Record a list of input files under `key` (skipped when empty).
    pub fn inputs<S: AsRef<str>>(mut self, key: &str, paths: &[S]) -> Self {
        if !paths.is_empty() {
            let v: Vec<Value> = paths.iter().map(|p| entry(p.as_ref())).collect();
            self.inputs.insert(key.into(), Value::Array(v));
        }
        self
    }

    /// Record every option, defaults included, as they are now.
    pub fn options<T: Serialize>(mut self, options: &T) -> Self {
        self.options = to_value(options);
        self
    }

    /// Replace one option with the value the run ended up using.
    pub fn set_option<T: Serialize>(&mut self, key: &str, value: &T) {
        if let Value::Object(o) = &mut self.options {
            o.insert(key.into(), to_value(value));
        }
    }

    /// Write the record with the outcome of the run. A failure to write is
    /// logged, not raised: it must never cost the run's own outputs.
    pub fn finish<R>(self, outcome: &anyhow::Result<R>) -> Option<PathBuf> {
        let status = match outcome {
            Ok(_) => "ok".to_string(),
            Err(e) => format!("failed: {e:#}"),
        };
        self.write(status)
            .map_err(|e| log::warn!("could not write the run record: {e:#}"))
            .ok()
    }

    fn write(self, status: String) -> anyhow::Result<PathBuf> {
        let outputs: Vec<String> = changed(&self.dir, &self.before)
            .into_iter()
            .filter(|n| *n != self.file_name)
            .collect();
        let unix = |t: SystemTime| {
            t.duration_since(SystemTime::UNIX_EPOCH)
                .map_or(0, |d| d.as_secs())
        };
        let record = json!({
            "faba_version": env!("CARGO_PKG_VERSION"),
            "job": self.job,
            "status": status,
            "command_line": std::env::args().collect::<Vec<_>>(),
            "working_directory": std::env::current_dir().ok(),
            "started_unix_seconds": unix(self.started),
            "elapsed_seconds": self.clock.elapsed().as_secs_f64(),
            "inputs": self.inputs,
            "output_directory": absolute(&self.dir.to_string_lossy()),
            "outputs": outputs,
            "options": self.options,
        });
        std::fs::create_dir_all(&self.dir)?;
        let path = self.dir.join(&self.file_name);
        std::fs::write(&path, serde_json::to_string_pretty(&record)?)
            .with_context(|| format!("writing {}", path.display()))?;
        log::info!(
            "Wrote run record (inputs, outputs, options) to {}",
            path.display()
        );
        Ok(path)
    }
}

/// Run `body` between [`RunRecord::start`] and [`RunRecord::finish`], and
/// return its result.
pub fn recorded<R>(
    record: RunRecord,
    body: impl FnOnce() -> anyhow::Result<R>,
) -> anyhow::Result<R> {
    let outcome = body();
    record.finish(&outcome);
    outcome
}

fn to_value<T: Serialize>(v: &T) -> Value {
    serde_json::to_value(v).unwrap_or_else(|e| Value::String(format!("unserializable: {e}")))
}

/// One input as recorded: `path`, absolute with symlinks kept as given (the
/// name the run was handed), and `resolved`, the file actually read, when a
/// symlink makes the two differ.
fn entry(p: &str) -> Value {
    let path = absolute(p);
    match std::fs::canonicalize(p) {
        Ok(r) if r.as_os_str() != path.as_str() => {
            json!({"path": path, "resolved": r.to_string_lossy()})
        }
        _ => json!({"path": path}),
    }
}

/// `p` made absolute against the working directory, symlinks kept.
fn absolute(p: &str) -> String {
    std::path::absolute(p).map_or_else(|_| p.to_string(), |a| a.to_string_lossy().into_owned())
}

fn snapshot(dir: &Path) -> FxHashMap<String, Option<SystemTime>> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return FxHashMap::default();
    };
    entries
        .flatten()
        .map(|e| {
            let mtime = e.metadata().and_then(|m| m.modified()).ok();
            (e.file_name().to_string_lossy().into_owned(), mtime)
        })
        .collect()
}

/// Entries of `dir` that are new or rewritten since `before`, sorted.
fn changed(dir: &Path, before: &FxHashMap<String, Option<SystemTime>>) -> Vec<String> {
    let mut out: Vec<String> = snapshot(dir)
        .into_iter()
        .filter(|(name, mtime)| before.get(name) != Some(mtime))
        .map(|(name, _)| name)
        .collect();
    out.sort();
    out
}

/// An input `key` (`gff`, `genome`, ...) recorded by a run whose outputs sit
/// next to `near` (an output file, or the output directory itself). The most
/// recent record naming an existing file wins. Returns the input's path and
/// the record it came from.
fn find_recorded(near: &str, key: &str) -> Option<(String, PathBuf)> {
    let p = Path::new(near);
    let dir = if p.is_dir() && !is_matrix_dir(p) {
        p.to_path_buf()
    } else {
        p.parent()
            .map(|d| {
                if d.as_os_str().is_empty() {
                    Path::new(".")
                } else {
                    d
                }
            })?
            .to_path_buf()
    };
    let mut records: Vec<(SystemTime, PathBuf)> = std::fs::read_dir(&dir)
        .ok()?
        .flatten()
        .filter(|e| is_run_record(&e.file_name().to_string_lossy()))
        .map(|e| {
            let t = e
                .metadata()
                .and_then(|m| m.modified())
                .unwrap_or(SystemTime::UNIX_EPOCH);
            (t, e.path())
        })
        .collect();
    records.sort_by_key(|a| std::cmp::Reverse(a.0));
    records.into_iter().find_map(|(_, path)| {
        let text = std::fs::read_to_string(&path).ok()?;
        let v: Value = serde_json::from_str(&text).ok()?;
        let found = choose(v.get("inputs")?.get(key)?, &path)?;
        Some((found, path))
    })
}

/// The file to use for a recorded input: the file the run read, by the name
/// it was given when that name still leads there.
///
/// - `path` still resolves to `resolved` (or was never a link): `path`.
/// - `path` now resolves elsewhere (a link repointed since the run) and
///   `resolved` still exists: `resolved`, with a warning, since that is what
///   the outputs were made from.
/// - `path` is gone: `resolved`, if it exists.
/// - `resolved` is gone too, or was never recorded: `path`, if it exists.
///
/// A plain string (a record written before `resolved` existed) is a `path`.
/// For a list (`bam`), the first entry.
fn choose(v: &Value, record: &Path) -> Option<String> {
    let v = v
        .as_array()
        .map_or(v, |a| a.first().unwrap_or(&Value::Null));
    let (path, resolved) = match v {
        Value::String(p) => (Some(p.as_str()), None),
        Value::Object(o) => (
            o.get("path").and_then(Value::as_str),
            o.get("resolved").and_then(Value::as_str),
        ),
        _ => return None,
    };
    let exists = |p: &&str| Path::new(p).exists();
    let (path, resolved) = (path.filter(exists), resolved.filter(exists));
    match (path, resolved) {
        (Some(p), Some(r)) => {
            let now = std::fs::canonicalize(p).ok();
            if now.as_deref() == Some(Path::new(r)) {
                Some(p.to_string())
            } else {
                log::warn!(
                    "{p} now points elsewhere than when {} was written; \
                     using the file that run read, {r}",
                    record.display()
                );
                Some(r.to_string())
            }
        }
        (None, Some(r)) => Some(r.to_string()),
        (Some(p), None) => Some(p.to_string()),
        (None, None) => None,
    }
}

/// Whether a file name is a run record: `{job}.run.json`, or `faba all`'s
/// `pipeline_summary.json`.
pub fn is_run_record(name: &str) -> bool {
    name.ends_with(".run.json") || name == "pipeline_summary.json"
}

/// An input `key` (`gff`, `genome`, ...) recorded next to `near` (an output
/// file, or an output directory), logged as the `what` in use.
pub fn find_input(near: &str, key: &str, what: &str) -> Option<Box<str>> {
    let (path, record) = find_recorded(near, key)?;
    log::info!("{what}: {path} (recorded in {})", record.display());
    Some(path.into_boxed_str())
}

/// `explicit` when given, else the input recorded next to `near`.
pub fn explicit_or_recorded(
    explicit: Option<&str>,
    near: &str,
    key: &str,
    what: &str,
) -> Option<Box<str>> {
    explicit
        .map(Box::from)
        .or_else(|| find_input(near, key, what))
}

/// A zarr store is a directory, but it is one output, not an output directory.
fn is_matrix_dir(p: &Path) -> bool {
    p.extension().is_some_and(|e| e == "zarr")
}

#[cfg(test)]
mod tests;
