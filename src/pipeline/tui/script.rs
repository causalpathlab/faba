//! `faba_run.cmd.sh`: the exact command `faba run` starts, as a script that
//! runs it again and refuses to run over a finished run.

use std::io::Write;
use std::path::{Path, PathBuf};

pub const SCRIPT: &str = "faba_run.cmd.sh";
/// The whole-run record: its presence means the run already finished.
pub const GUARD: &str = "pipeline_summary.json";

/// `word` as one shell word.
pub fn quote(word: &str) -> String {
    let plain = !word.is_empty()
        && word
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "_@%+=:,./-".contains(c));
    if plain {
        word.to_string()
    } else {
        format!("'{}'", word.replace('\'', r"'\''"))
    }
}

/// Whether a word is a flag: `-[A-Za-z]` or `--[A-Za-z]...`.
fn is_flag(w: &str) -> bool {
    let mut c = w.chars();
    match (c.next(), c.next(), c.next()) {
        (Some('-'), Some(x), None) => x.is_ascii_alphabetic(),
        (Some('-'), Some('-'), Some(x)) => x.is_ascii_alphabetic(),
        _ => false,
    }
}

/// The program and `run --batch-process` on the first line, each BAM on its
/// own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String]) -> Vec<String> {
    let mut lines = vec![String::from("\"${FABA:-faba}\"")];
    // A flag has started a line: the words that follow, up to the next
    // flag, are its values.
    let mut in_flag = false;
    for (k, w) in argv.iter().enumerate() {
        let flag = is_flag(w);
        // `run` and `--batch-process` stay on the program line.
        let add_to_last = k < 2 || (!flag && in_flag);
        match lines.last_mut() {
            Some(last) if add_to_last => {
                last.push(' ');
                last.push_str(&quote(w));
            }
            _ => {
                lines.push(quote(w));
                in_flag |= flag;
            }
        }
    }
    lines
}

pub fn text(argv: &[String]) -> String {
    let mut s = String::new();
    s.push_str("#!/usr/bin/env bash\n");
    s.push_str(&format!(
        "# Made by `faba run` (faba {}). Run it again with: bash {SCRIPT}\n",
        env!("CARGO_PKG_VERSION")
    ));
    s.push_str("set -euo pipefail\n");
    s.push_str("cd \"$(dirname \"$0\")\"\n");
    s.push_str(&format!("if [ -f {GUARD} ]; then\n"));
    s.push_str(&format!(
        "  echo \"{GUARD} exists; move the outputs away to run again\" >&2\n"
    ));
    s.push_str("  exit 1\nfi\n");
    s.push_str(&command_lines(argv).join(" \\\n  "));
    s.push('\n');
    s
}

/// Write the script into `dir`, never over an existing one.
pub fn write(dir: &Path, argv: &[String]) -> anyhow::Result<PathBuf> {
    let path = dir.join(SCRIPT);
    let mut opts = std::fs::OpenOptions::new();
    opts.write(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        opts.mode(0o755);
    }
    let mut f = opts
        .open(&path)
        .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    f.write_all(text(argv).as_bytes())?;
    Ok(path)
}

#[cfg(test)]
#[path = "tests/script.rs"]
mod tests;
