//! The command a run made, saved as a shell script next to its outputs
//! (`faba_run.cmd.sh`, `{modality}_cluster.cmd.sh`, `{modality}_blockmodel.cmd.sh`), so it can be read,
//! edited and run again.

use clap::{Arg, ArgAction, ArgMatches};
use std::io::Write;
use std::path::{Path, PathBuf};

/// How several values of one flag reach the command line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Many {
    One,
    /// Joined into one argument: `--flag a,b`.
    Joined(char),
    /// Listed after one flag: `--flag a b`.
    Listed,
    /// The flag repeated: `--flag a --flag b`.
    Repeated,
}

/// How `a` takes several values, by its clap definition.
pub(crate) fn many_of(a: &Arg) -> Many {
    if let Some(c) = a.get_value_delimiter() {
        Many::Joined(c)
    } else if a.get_num_args().is_some_and(|r| r.max_values() > 1) {
        Many::Listed
    } else if matches!(a.get_action(), ArgAction::Append) {
        Many::Repeated
    } else {
        Many::One
    }
}

/// `flag` with values `vs`, laid out the way the flag takes them. A value
/// that starts like a flag goes in its flag's word: `--flag=-1`.
pub(crate) fn emit(flag: &str, many: Many, vs: Vec<String>) -> Vec<String> {
    let one = |v: String| {
        if v.starts_with('-') {
            vec![format!("{flag}={v}")]
        } else {
            vec![flag.to_string(), v]
        }
    };
    let Some(first) = vs.first() else {
        return Vec::new();
    };
    match many {
        Many::One => one(first.clone()),
        Many::Joined(c) => one(vs.join(&c.to_string())),
        Many::Listed => std::iter::once(flag.to_string()).chain(vs).collect(),
        Many::Repeated => vs.into_iter().flat_map(one).collect(),
    }
}

/// The words `id` was given in `m`, as typed, when it was given on the
/// command line or by the environment; none when clap took its default.
pub(crate) fn given_words(m: &ArgMatches, id: &str) -> Option<Vec<String>> {
    use clap::parser::ValueSource;
    let given = matches!(
        m.value_source(id),
        Some(ValueSource::CommandLine | ValueSource::EnvVariable)
    );
    given.then(|| {
        m.get_raw(id)
            .into_iter()
            .flatten()
            .map(|v| v.to_string_lossy().into_owned())
            .collect()
    })
}

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
