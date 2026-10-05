//! `faba run` — the end-to-end pipeline that chains the per-modality
//! subcommands over one set of BAM files.
//!
//! Steps shared with the standalone subcommands (BAM index checks, mito QC,
//! gene-expression QC) live in [`crate::quant`], which the standalone
//! entries also use.

/// The `faba run` command-line surface.
pub mod args;
/// The `faba run` run. Binary entry: [`run::run_pipeline`].
pub mod run;
/// The per-modality steps, in run order.
mod steps;
/// The `faba run` setup view.
mod tui;

#[cfg(test)]
#[path = "tests/args.rs"]
mod args_tests;
