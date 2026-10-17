//! The command `faba cluster` saves beside a run.

use crate::cluster::{command, ClusterArgs, ModelKind, ModulesArgs};

/// `words` parsed as `faba cluster` parses them, what was given recorded.
fn parse_cluster(words: &[&str]) -> ClusterArgs {
    use clap::{Args, FromArgMatches};
    let cmd = ClusterArgs::augment_args(clap::Command::new("faba"));
    let m = cmd.try_get_matches_from(words).unwrap();
    let mut args = ClusterArgs::from_arg_matches(&m).unwrap();
    args.record_given(&m);
    args
}

/// `argv`, a saved command, parsed back as `faba cluster` parses it.
fn parse_saved(argv: &[String]) -> ClusterArgs {
    let words: Vec<&str> = std::iter::once("faba")
        .chain(argv[1..].iter().map(String::as_str))
        .collect();
    parse_cluster(&words)
}

/// An option set in code, with no setter noting it, is in the saved
/// command all the same: the command is written from the values run.
#[test]
fn an_option_set_in_code_is_in_the_saved_command() {
    let mut args = parse_cluster(&["faba", "/x.zarr", "--cell-roles", "4"]);
    args.cell_roles = 9;
    args.model = ModelKind::Crp;
    let argv = command(&args);
    assert_eq!(after(&argv, "--cell-roles"), "9");
    assert_eq!(after(&argv, "--model"), "crp");
    assert_eq!(
        serde_json::to_value(parse_saved(&argv)).unwrap(),
        serde_json::to_value(&args).unwrap()
    );
}

/// The word after `flag` in `argv`.
fn after<'a>(argv: &'a [String], flag: &str) -> &'a str {
    let k = argv.iter().position(|w| w == flag);
    &argv[k.unwrap_or_else(|| panic!("no {flag} in {argv:?}")) + 1]
}

#[test]
fn the_saved_command_runs_the_same_fit_again() {
    let args = parse_cluster(&[
        "faba",
        "rel/a_m6a_site.zarr",
        "--model",
        "blocks",
        "--cell-roles",
        "6",
        "--fit",
        "*_wt_*",
        "--no-batch-correction",
        "--holdout",
        "0.1",
        "--shift-sd=-0.5",
        "-o",
        "out",
    ]);
    let argv = command(&args);
    assert_eq!(argv[0], "cluster");
    assert!(argv[1].ends_with("rel/a_m6a_site.zarr") && argv[1].starts_with('/'));
    // The output too, so the script writes where the run did from anywhere;
    // a pattern from `*` matches the absolute inputs as it is.
    assert!(
        after(&argv, "--output").starts_with('/') && after(&argv, "--output").ends_with("/out")
    );
    assert_eq!(after(&argv, "--fit"), "*_wt_*");
    // A value that starts like a flag stays in its flag's word.
    assert!(argv.contains(&"--shift-sd=-0.5".to_string()), "{argv:?}");
    // Parsed back, the same options; only the paths are absolute now.
    let mut again = parse_saved(&argv);
    again.data_files = args.data_files.clone();
    again.output = args.output.clone();
    assert_eq!(
        serde_json::to_value(&again).unwrap(),
        serde_json::to_value(&args).unwrap()
    );
    // A relative pattern is made absolute, as the inputs it matches are.
    let rel = command(&parse_cluster(&[
        "faba",
        "rel/a.zarr",
        "--fit",
        "rel/*_wt_*",
    ]));
    let fit = after(&rel, "--fit");
    assert!(fit.starts_with('/') && fit.ends_with("rel/*_wt_*"), "{fit}");
    // An absolute one is left as it is.
    let abs = command(&parse_cluster(&["faba", "/x.zarr", "--fit", "/d/*_wt_*"]));
    assert_eq!(after(&abs, "--fit"), "/d/*_wt_*");
    // Defaults are left out.
    let plain = command(&parse_cluster(&["faba", "/x.zarr"]));
    assert_eq!(plain, vec!["cluster".to_string(), "/x.zarr".to_string()]);
    // As a script: the program and `cluster` first, a flag a line.
    let script = crate::script::Script {
        name: "m6a_cluster.cmd.sh",
        made_by: "faba cluster",
        guard: None,
        cd: false,
        head: 1,
        replace: true,
    };
    let text = crate::script::text(&script, &argv);
    assert!(text.contains("\"${FABA:-faba}\" cluster \\\n"), "{text}");
    assert!(text.contains("  --cell-roles 6"), "{text}");
    assert!(text.contains("  --fit '*_wt_*'"), "{text}");
}

/// A default given on the command line is kept, as given: the script says
/// what was asked, not only what differs from today's defaults.
#[test]
fn the_saved_command_keeps_a_default_given() {
    let argv = command(&parse_cluster(&[
        "faba",
        "/x.zarr",
        "--chains",
        "4",
        "--model=link",
    ]));
    assert_eq!(after(&argv, "--chains"), "4");
    assert_eq!(after(&argv, "--model"), "link");
}

/// A path with a space is made absolute and stays one word in the script.
#[test]
fn the_saved_command_quotes_a_path_with_a_space() {
    let argv = command(&parse_cluster(&["faba", "/x y/a.zarr", "-o", "my out"]));
    assert_eq!(argv[1], "/x y/a.zarr");
    let out = after(&argv, "--output");
    assert!(out.starts_with('/') && out.ends_with("/my out"), "{out}");
    let script = crate::script::Script {
        name: "m6a_cluster.cmd.sh",
        made_by: "faba cluster",
        guard: None,
        cd: false,
        head: 1,
        replace: true,
    };
    let text = crate::script::text(&script, &argv);
    assert!(text.contains(&format!("--output '{out}'")), "{text}");
    assert!(text.contains("'/x y/a.zarr'"), "{text}");
}

/// `faba modules` clusters with `--model blocks` whatever was given, so its
/// saved command says so, beside the options given to `faba modules`.
#[test]
fn the_saved_command_of_modules_says_blocks() {
    use clap::{Args, FromArgMatches};
    let cmd = ModulesArgs::augment_args(clap::Command::new("faba"));
    let m = cmd
        .try_get_matches_from(["faba", "/x.zarr", "--chains", "2", "--refit"])
        .unwrap();
    let mut args = ModulesArgs::from_arg_matches(&m).unwrap();
    args.record_given(&m);
    let argv = command(&args.blocks());
    assert_eq!(after(&argv, "--model"), "blocks");
    assert_eq!(after(&argv, "--chains"), "2");
    // Only `faba cluster`'s own options.
    assert!(!argv.iter().any(|w| w == "--refit"), "{argv:?}");
    // Set in code over one given: the last word is the one run.
    let m = ModulesArgs::augment_args(clap::Command::new("faba"))
        .try_get_matches_from(["faba", "/x.zarr", "--model", "crp"])
        .unwrap();
    let mut args = ModulesArgs::from_arg_matches(&m).unwrap();
    args.record_given(&m);
    let argv = command(&args.blocks());
    assert_eq!(
        argv.iter().filter(|w| *w == "--model").count(),
        1,
        "{argv:?}"
    );
    assert_eq!(after(&argv, "--model"), "blocks");
}
