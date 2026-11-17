//! The command a run made, saved as a shell script next to its outputs
//! (`faba_run.cmd.sh`, `{modality}_cluster.cmd.sh`, `{modality}_blockmodel.cmd.sh`), so it can be read,
//! edited and run again: legume-tui's [`legume_tui::form::script`], calling
//! `"${FABA:-faba}"`.

use legume_tui::form::script::{self, Script};
use std::path::{Path, PathBuf};

/// The script `name`, made by `made_by`, calling this faba; the builder
/// adds a guard, the folder it runs from and the program's line.
pub const fn faba<'a>(name: &'a str, made_by: &'a str) -> Script<'a> {
    Script::new("faba", env!("CARGO_PKG_VERSION"), made_by, name)
}

/// The program and the first `head` words on the first line, each input on
/// its own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String], head: usize) -> Vec<String> {
    script::command_lines(&script::program_word("faba"), argv, head)
}

/// Write the script into `dir`, executable; over an older one only when
/// `replace`, and then as a new file renamed over it: a `bash` still
/// reading the old one (the run it started is the one rewriting it) keeps
/// reading the old text, a link to it is left as it was, and a write that
/// fails leaves the old script whole.
pub fn write(dir: &Path, s: &Script, argv: &[String], replace: bool) -> anyhow::Result<PathBuf> {
    let path = dir.join(s.name);
    if !replace {
        script::write(&path, s, argv)?;
        return Ok(path);
    }
    let fresh = dir.join(format!(".{}.new", s.name));
    // One left by a write that died before its rename.
    let _ = std::fs::remove_file(&fresh);
    script::write(&fresh, s, argv)?;
    std::fs::rename(&fresh, &path).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A script rewritten each run is a new file: what still holds the old
    /// one (a running `bash`, a link) keeps the old text.
    #[test]
    fn a_replaced_script_is_a_new_file() {
        let tmp = tempfile::tempdir().unwrap();
        let s = faba("x.cmd.sh", "faba cluster").head(1);
        let first = write(tmp.path(), &s, &["cluster".into(), "a".into()], true).unwrap();
        let kept = tmp.path().join("kept.sh");
        std::fs::hard_link(&first, &kept).unwrap();
        let again = write(tmp.path(), &s, &["cluster".into(), "b".into()], true).unwrap();
        assert!(std::fs::read_to_string(&again)
            .unwrap()
            .contains("cluster \\\n  b"));
        let old = std::fs::read_to_string(&kept).unwrap();
        assert!(old.contains("cluster \\\n  a"), "{old}");
        assert!(
            write(tmp.path(), &s, &["cluster".into()], false).is_err(),
            "not over one unasked"
        );
    }
}
