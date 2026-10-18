//! A `faba cluster` run as the command that would run it again: each option
//! given or set away from its default, by the flag clap names it, with the
//! value the run had.

use super::ClusterArgs;
use clap::{ArgMatches, Args, FromArgMatches, ValueHint};
use rustc_hash::FxHashSet;
use serde_json::{Map, Value};

/// The options given to `faba cluster` (on the command line or by the
/// environment), by clap id.
#[derive(Debug, Clone, Default)]
pub struct Given(FxHashSet<String>);

/// `faba cluster`'s options alone, as clap defines them.
pub(super) fn definition() -> clap::Command {
    ClusterArgs::augment_args(clap::Command::new("cluster"))
}

impl ClusterArgs {
    /// Notes which of `faba cluster`'s options `matches` were given, for
    /// the command saved beside the run.
    pub fn record_given(&mut self, matches: &ArgMatches) {
        self.given = Given(
            definition()
                .get_arguments()
                .map(|a| a.get_id().as_str())
                .filter(|id| crate::script::given_words(matches, id).is_some())
                .map(String::from)
                .collect(),
        );
    }

    /// Whether option `id` (its clap id) was given.
    pub(super) fn given(&self, id: &str) -> bool {
        self.given.0.contains(id)
    }
}

/// `word` as the saved command takes it: a word of an option that takes a
/// path (or, as `--fit`, a pattern over paths) made absolute, so the script
/// runs from anywhere, unless it starts with a wildcard.
fn absolute(arg: &clap::Arg, word: &str) -> String {
    let path = matches!(
        arg.get_value_hint(),
        ValueHint::AnyPath | ValueHint::FilePath | ValueHint::DirPath
    );
    if path && !word.starts_with(['*', '?']) {
        crate::run_record::absolute(word)
    } else {
        word.to_string()
    }
}

/// `args`' options by field name, which is each one's clap id.
fn fields(args: &ClusterArgs) -> Map<String, Value> {
    match serde_json::to_value(args) {
        Ok(Value::Object(m)) => m,
        _ => Map::new(),
    }
}

/// The words of one option's value: none for unset, each item of a list.
fn words(v: &Value) -> Vec<String> {
    match v {
        Value::Null => Vec::new(),
        Value::String(s) => vec![s.clone()],
        Value::Array(items) => items.iter().flat_map(words).collect(),
        other => vec![other.to_string()],
    }
}

/// The words option `id` (its clap id) has in `args`, as the saved
/// command says them: none when unset.
pub(super) fn value_words(args: &ClusterArgs, id: &str) -> Vec<String> {
    fields(args).get(id).map(words).unwrap_or_default()
}

/// `args` as a `faba cluster` command: its inputs as absolute paths, then
/// each option given or away from its default, by the flag clap names it,
/// with the value this run had, paths among them absolute too.
pub(super) fn command(args: &ClusterArgs) -> Vec<String> {
    command_except(args, &[])
}

/// [`command`] without the options whose long flags are in `skip`.
pub(super) fn command_except(args: &ClusterArgs, skip: &[&str]) -> Vec<String> {
    let cmd = definition();
    let ran = fields(args);
    let default = cmd
        .clone()
        .try_get_matches_from(["cluster"])
        .ok()
        .and_then(|m| ClusterArgs::from_arg_matches(&m).ok())
        .map(|d| fields(&d))
        .unwrap_or_default();
    let mut argv = vec!["cluster".to_string()];
    argv.extend(
        args.data_files
            .iter()
            .map(|f| crate::run_record::absolute(f)),
    );
    for arg in cmd.get_arguments() {
        let id = arg.get_id().as_str();
        let (Some(long), Some(value)) = (arg.get_long(), ran.get(id)) else {
            continue;
        };
        if skip.contains(&long) {
            continue;
        }
        let changed = default.get(id) != Some(value);
        let flag = format!("--{long}");
        if !arg.get_action().takes_values() {
            // A switch is its flag alone, written when away from its default.
            if changed {
                argv.push(flag);
            }
        } else if changed || args.given(id) {
            let ws = words(value).iter().map(|w| absolute(arg, w)).collect();
            argv.extend(crate::script::emit(&flag, crate::script::many_of(arg), ws));
        }
    }
    argv
}
