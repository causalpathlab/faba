//! `faba_run.cmd.sh`: the exact command `faba run` starts, as a script that
//! runs it again and refuses to run over a finished run (see
//! [`crate::script`]).

use crate::script::Script;
use std::path::{Path, PathBuf};

pub const SCRIPT: &str = "faba_run.cmd.sh";
/// The whole-run record: its presence means the run already finished.
pub const GUARD: &str = "pipeline_summary.json";

const RUN: Script = Script {
    name: SCRIPT,
    made_by: "faba run",
    guard: Some(GUARD),
    cd: true,
    head: 2,
    replace: false,
};

/// The program and `run --batch-process` on the first line, each BAM on its
/// own line, then one flag and all its values per line.
pub fn command_lines(argv: &[String]) -> Vec<String> {
    crate::script::command_lines(argv, RUN.head)
}

/// Write the script into `dir`, never over an existing one.
pub fn write(dir: &Path, argv: &[String]) -> anyhow::Result<PathBuf> {
    crate::script::write(dir, &RUN, argv)
}

#[cfg(test)]
#[path = "tests/script.rs"]
mod tests;
