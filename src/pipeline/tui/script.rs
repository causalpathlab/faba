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

/// Test if a word is a flag token: `-[A-Za-z]` or `--[A-Za-z]`.
fn is_flag(w: &str) -> bool {
    (w.len() == 2 && w.starts_with('-') && w.chars().nth(1).unwrap().is_ascii_alphabetic())
        || (w.len() > 2 && w.starts_with("--") && w.chars().nth(2).unwrap().is_ascii_alphabetic())
}

/// The program and `run --batch-process` on the first line, each BAM on its
/// own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String]) -> Vec<String> {
    let mut lines = vec![String::from("\"${FABA:-faba}\"")];
    let mut last_flag_idx = None;

    for (k, w) in argv.iter().enumerate() {
        let is_flag_token = is_flag(w);
        let in_head = k < 2; // `run` and `--batch-process` stay on the program line

        // Add to last line if in head or if this is a value for a flag at index >= 2
        let add_to_last =
            in_head || (!is_flag_token && last_flag_idx.map_or(false, |idx| idx >= 2));

        match lines.last_mut() {
            Some(last) if add_to_last => {
                last.push(' ');
                last.push_str(&quote(w));
            }
            _ => {
                lines.push(quote(w));
                if is_flag_token {
                    last_flag_idx = Some(k);
                }
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
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        let mut f = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o755)
            .open(&path)
            .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
        f.write_all(text(argv).as_bytes())?;
    }
    #[cfg(not(unix))]
    {
        let mut f = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
        f.write_all(text(argv).as_bytes())?;
    }
    Ok(path)
}

#[cfg(test)]
#[path = "tests/script.rs"]
mod tests;
