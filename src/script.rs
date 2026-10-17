//! The command a run made, saved as a shell script next to its outputs
//! (`faba_run.cmd.sh`, `{modality}_cluster.cmd.sh`), so it can be read,
//! edited and run again.

use std::io::Write;
use std::path::{Path, PathBuf};

/// How a script is made: its file name, the command that made it, the
/// output whose presence means the run finished (it then refuses to run),
/// whether it runs from its own folder, how many words go on the program's
/// line, and whether it is written over an older one.
pub struct Script<'a> {
    pub name: &'a str,
    pub made_by: &'a str,
    pub guard: Option<&'a str>,
    pub cd: bool,
    pub head: usize,
    pub replace: bool,
}

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

/// The program and the first `head` words on the first line, each input on
/// its own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String], head: usize) -> Vec<String> {
    let mut lines = vec![String::from("\"${FABA:-faba}\"")];
    // A flag has started a line: the words that follow, up to the next
    // flag, are its values.
    let mut in_flag = false;
    for (k, w) in argv.iter().enumerate() {
        let flag = is_flag(w);
        let add_to_last = k < head || (!flag && in_flag);
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

/// The script running `argv`.
pub fn text(script: &Script, argv: &[String]) -> String {
    let Script {
        name,
        made_by,
        guard,
        cd,
        head,
        ..
    } = script;
    let mut s = String::new();
    s.push_str("#!/usr/bin/env bash\n");
    s.push_str(&format!(
        "# Made by `{made_by}` (faba {}). Run it again with: bash {name}\n",
        env!("CARGO_PKG_VERSION")
    ));
    s.push_str("set -euo pipefail\n");
    if *cd {
        s.push_str("cd \"$(dirname \"$0\")\"\n");
    }
    if let Some(guard) = guard {
        s.push_str(&format!("if [ -f {guard} ]; then\n"));
        s.push_str(&format!(
            "  echo \"{guard} exists; move the outputs away to run again\" >&2\n"
        ));
        s.push_str("  exit 1\nfi\n");
    }
    s.push_str(&command_lines(argv, *head).join(" \\\n  "));
    s.push('\n');
    s
}

/// Write the script into `dir`, executable; over an older one only when
/// it says so.
pub fn write(dir: &Path, script: &Script, argv: &[String]) -> anyhow::Result<PathBuf> {
    let path = dir.join(script.name);
    let mut opts = std::fs::OpenOptions::new();
    opts.write(true);
    if script.replace {
        opts.create(true).truncate(true);
    } else {
        opts.create_new(true);
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        opts.mode(0o755);
    }
    let mut f = opts
        .open(&path)
        .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    f.write_all(text(script, argv).as_bytes())?;
    Ok(path)
}
