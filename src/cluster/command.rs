//! A `faba cluster` run as the command that would run it again, its options
//! by the flags clap names them.

use super::ClusterArgs;
use rustc_hash::FxHashMap;

/// `faba cluster`'s options alone, to parse them or name their flags.
#[derive(clap::Parser)]
pub(super) struct Alone {
    #[command(flatten)]
    pub(super) args: ClusterArgs,
}

/// `args` as a `faba cluster` command: its inputs as absolute paths, then
/// each option that differs from its default, by the flag clap names it,
/// paths among them absolute too.
pub(super) fn command(args: &ClusterArgs) -> Vec<String> {
    let defaults = <Alone as clap::Parser>::parse_from(["faba", "x"]).args;
    let flags: FxHashMap<String, String> = <Alone as clap::CommandFactory>::command()
        .get_arguments()
        .filter_map(|a| Some((a.get_id().to_string(), a.get_long()?.to_string())))
        .collect();
    let mut argv = vec!["cluster".to_string()];
    argv.extend(
        args.data_files
            .iter()
            .map(|f| crate::run_record::absolute(f)),
    );
    let (Ok(serde_json::Value::Object(given)), Ok(defaults)) =
        (serde_json::to_value(args), serde_json::to_value(&defaults))
    else {
        return argv;
    };
    for (field, mut value) in given {
        // Paths, as the inputs are, so the script runs from anywhere: the
        // output, and a pattern that does not start with a wildcard.
        if let serde_json::Value::String(s) = &value {
            let path = field == "output" || (field == "fit" && !s.starts_with(['*', '?', '/']));
            if path {
                value = crate::run_record::absolute(s).into();
            }
        }
        let Some(flag) = flags
            .get(&field)
            .filter(|_| defaults.get(&field) != Some(&value))
        else {
            continue;
        };
        push_option(&mut argv, &format!("--{flag}"), &value);
    }
    argv
}

/// `--flag` with `value` as the command line takes it: a switch alone, a
/// list as its values, a value that starts like a flag in the same word.
fn push_option(argv: &mut Vec<String>, flag: &str, value: &serde_json::Value) {
    use serde_json::Value;
    let words: Vec<String> = match value {
        Value::Bool(true) => return argv.push(flag.to_string()),
        Value::Number(n) => vec![n.to_string()],
        Value::String(s) => vec![s.clone()],
        Value::Array(items) => items
            .iter()
            .filter_map(|v| match v {
                Value::String(s) => Some(s.clone()),
                Value::Number(n) => Some(n.to_string()),
                _ => None,
            })
            .collect(),
        _ => return,
    };
    match &words[..] {
        [one] if one.starts_with('-') => argv.push(format!("{flag}={one}")),
        _ => {
            argv.push(flag.to_string());
            argv.extend(words);
        }
    }
}
