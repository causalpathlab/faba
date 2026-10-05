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

/// The program and `run --batch-process` on the first line, each BAM on its
/// own line, then one flag and its values per line.
pub fn command_lines(argv: &[String]) -> Vec<String> {
    let mut lines = vec![String::from("\"${FABA:-faba}\"")];
    for (k, w) in argv.iter().enumerate() {
        let starts = w.starts_with('-');
        let is_batch_process = k > 0 && argv[k - 1] == "--batch-process";
        let after_flag = k > 0 && argv[k - 1].starts_with('-') && !starts && !is_batch_process;
        let head = k < 2; // `run` and `--batch-process` stay on the program line
        match lines.last_mut() {
            Some(last) if head || after_flag => {
                last.push(' ');
                last.push_str(&quote(w));
            }
            _ => lines.push(quote(w)),
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
    let mut f = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    f.write_all(text(argv).as_bytes())?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755))?;
    }
    Ok(path)
}

#[cfg(test)]
#[path = "tests/script.rs"]
mod tests;
