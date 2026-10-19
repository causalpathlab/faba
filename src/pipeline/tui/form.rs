//! The form of `faba run`'s flags: every flag but those the Inputs and
//! Steps screens own (see [`crate::tui::form`] for the rows).

pub use crate::tui::form::*;
use clap::{Arg, Command};

/// Flags with no help heading that are never passed: the setup view's
/// own, its session's among them. The Inputs and Steps screens own every
/// flag under [`OWNED_HEADINGS`].
pub const OWN: &[&str] = &[
    "batch-process",
    "help",
    "version",
    "verbose",
    "no-detach",
    "session",
];

/// The help headings whose flags the Inputs and Steps screens own.
pub const OWNED_HEADINGS: &[&str] = &["Inputs", "Steps"];

/// Whether the Flags screen shows `a`: not one of [`OWN`], nor under
/// [`OWNED_HEADINGS`].
fn shown(a: &Arg) -> bool {
    let own = a.get_long().is_some_and(|l| OWN.contains(&l));
    let heading = a.get_help_heading().unwrap_or("");
    !own && !OWNED_HEADINGS.contains(&heading)
}

/// The Flags screen's form for `run`, which must be built.
pub fn flags(run: &Command) -> Form {
    Form::new(run, shown)
}

#[cfg(test)]
#[path = "tests/form.rs"]
mod tests;
