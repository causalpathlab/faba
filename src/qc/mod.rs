//! `faba qc`: the one place faba thresholds anything.
//!
//! The producers (`count`, `dartseq`, `atoi`, `apa`, `all`) are inclusive:
//! they write every called cell, every gene with a count, and every putative
//! editing site with its statistics, and they apply no p-value, effect-size or
//! reproducibility cutoff. `qc` reads such a directory, shows what each site
//! threshold keeps so the cut is chosen on the data, and writes a new,
//! filtered fileset. It never touches a BAM.

pub mod args;
pub mod browser;
pub mod layout;
pub mod matrix;
pub mod path_tui;
pub mod progress;
pub mod repool;
pub mod run;
pub mod site_tui;
pub mod sites;
pub mod widgets;

pub use args::QcArgs;
pub use run::run_qc;
