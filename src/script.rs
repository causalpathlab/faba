//! The command a run made, saved as a shell script next to its outputs
//! (`faba_run.cmd.sh`, `{modality}_cluster.cmd.sh`, `{modality}_blockmodel.cmd.sh`), so it can be read,
//! edited and run again: legume-tui's [`legume_tui::form::script`], calling
//! `"${FABA:-faba}"`.

use legume_tui::form::script::{self, Script};
use std::path::{Path, PathBuf};

/// The script `name`, made by `made_by`, calling faba: refusing to run once
/// `guard` exists, run from its own folder when `cd`, with `head` words of
/// the command on the program's line.
pub const fn faba<'a>(
    name: &'a str,
    made_by: &'a str,
    guard: Option<&'a str>,
    cd: bool,
    head: usize,
) -> Script<'a> {
    Script {
        program: "faba",
        version: env!("CARGO_PKG_VERSION"),
        made_by,
        name,
        guard,
        cd,
        env: &[],
        head,
    }
}

/// The program and the first `head` words on the first line, each input on
/// its own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String], head: usize) -> Vec<String> {
    script::command_lines(&script::program_word("faba"), argv, head)
}

/// Write the script into `dir`, executable; over an older one only when
/// `replace`.
pub fn write(dir: &Path, s: &Script, argv: &[String], replace: bool) -> anyhow::Result<PathBuf> {
    let path = dir.join(s.name);
    if replace {
        if let Err(e) = std::fs::remove_file(&path) {
            if e.kind() != std::io::ErrorKind::NotFound {
                anyhow::bail!("{}: {e}", path.display());
            }
        }
    }
    script::write(&path, s, argv)?;
    Ok(path)
}
