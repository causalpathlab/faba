//! `faba qc`: the one place faba thresholds anything.
//!
//! The producers (`count`, `dartseq`, `atoi`, `apa`, `run`) are inclusive:
//! they write every called cell, every gene with a count, and every putative
//! editing site with its statistics, and they apply no p-value, effect-size or
//! reproducibility cutoff. `qc` reads such a directory, shows what each site
//! threshold keeps so the cut is chosen on the data, and writes a new,
//! filtered fileset. It never touches a BAM.

pub mod args;
pub mod cells_plot;
pub mod density;
pub mod layout;
pub mod matrix;
pub mod pairs;
pub mod path_tui;
pub mod progress;
pub mod repool;
pub mod run;
pub mod scale;
pub mod site_plot;
pub mod site_tui;
pub mod sites;
pub mod view;

#[cfg(test)]
#[path = "tests/fixtures.rs"]
mod fixtures;

pub use args::QcArgs;
pub use run::run_qc;
