//! The report's Model section: the run's settings from its record, its
//! plate diagram drawn as they are set; and what `faba modules` makes of
//! that record: the tables opened as recorded, or fitted again.

use super::tables;
use crate::cluster::browse::{plan, report, title, Modules, ModulesArgs, Plan};
use crate::cluster::{args, record_name, ClusterArgs};
use crate::run_record::RunRecord;
use serde_json::Value;
use std::path::Path;

/// `words` after `faba cluster x.zarr`, parsed as it parses them.
fn cluster(words: &[&str]) -> ClusterArgs {
    let all: Vec<&str> = ["faba", "/x.zarr"].iter().chain(words).copied().collect();
    crate::cluster::tests::parse_args(&all)
}

/// The fixture's tables in `dir`, beside the record `faba cluster` leaves
/// for `words`.
fn recorded(dir: &Path, words: &[&str]) {
    tables(dir);
    let args = cluster(words);
    let mut record = RunRecord::start("cluster", &dir.to_string_lossy())
        .file_name(&record_name("m6a", args.model))
        .options(&args);
    record.set_result("model", &args::recorded(&args));
    record.finish(&Ok(()));
}

/// The page's model section for the tables in `dir`.
fn section(dir: &Path) -> Value {
    let m = Modules::load(dir, "m6a").unwrap();
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    data["model"].clone()
}

/// The section of a run of `words`.
fn run(words: &[&str]) -> Value {
    let tmp = tempfile::tempdir().unwrap();
    recorded(tmp.path(), words);
    section(tmp.path())
}

/// Whether the diagram draws `node`, and whether greyed (fixed).
fn drawn(svg: &str, node: &str) -> Option<bool> {
    let at = format!(r#"data-node="{node}""#);
    let i = svg.find(&at)?;
    Some(svg[i + at.len()..].starts_with(" data-fixed"))
}

/// The blocks preset: its shift drawn, `psi` drawn, `beta` fixed, no `p`,
/// and each cell's own level `nu_c` in the cells plate with its prior sd,
/// an edge to the reads.
#[test]
fn the_blocks_preset_is_drawn_as_set() {
    let s = run(&["--model", "blocks", "--cell-roles", "6", "--cell-sd", "0.7"]);
    assert_eq!(s["preset"], "blocks");
    assert_eq!(
        [&s["cellRoles"], &s["geneRoles"], &s["strength"]],
        [&Value::from(6), &Value::from(8), &Value::from("inf")]
    );
    assert_eq!(
        (s["alpha"].as_f64(), s["geneAlpha"].as_f64()),
        (Some(0.5), Some(0.5))
    );
    assert_eq!(s["shiftSd"], 1.0);
    assert!(s["gamma"].is_null(), "symmetric weights have no gamma");
    assert_eq!(s["cellSd"], 0.7, "the levels' prior sd");
    let svg = s["svg"].as_str().unwrap();
    assert_eq!(drawn(svg, "nu"), Some(false), "each cell's level");
    assert_eq!(drawn(svg, "sigma-nu"), Some(false));
    assert_eq!(drawn(svg, "sigma-nu-nu"), Some(false));
    assert_eq!(
        drawn(svg, "nu-y"),
        Some(false),
        "the level reaches the reads"
    );
    assert!(svg.contains("= 0.7"), "the run's cell_sd");
    assert_eq!(drawn(svg, "delta"), Some(false));
    assert_eq!(drawn(svg, "sigma-delta"), Some(false));
    assert_eq!(drawn(svg, "psi"), Some(false), "B = 8");
    assert_eq!(drawn(svg, "beta"), Some(true), "fixed at 1 / A");
    assert_eq!(drawn(svg, "gamma"), None);
    assert_eq!(drawn(svg, "p"), None, "S = inf: no p");
    assert_eq!(drawn(svg, "strength"), None);
    assert_eq!(
        drawn(svg, "mu-y"),
        Some(false),
        "reads straight from the baseline"
    );
    assert_eq!(drawn(svg, "delta-y"), Some(false));
    assert!(
        svg.contains("1…6") && svg.contains("1…8"),
        "the run's roles"
    );
    assert!(!svg.contains("http"), "nothing loaded");
}

/// The link preset: no shift, one gene role greyed, `beta` drawn, `p`
/// integrated out, and no cell level (none at finite strength).
#[test]
fn the_link_preset_is_drawn_as_set() {
    let s = run(&["--model", "link"]);
    assert_eq!(s["preset"], "link");
    assert_eq!(
        [
            &s["cellRoles"],
            &s["geneRoles"],
            &s["strength"],
            &s["gamma"]
        ],
        [
            &Value::from(8),
            &Value::from(1),
            &Value::from(5.0),
            &Value::from(1.0)
        ]
    );
    assert!(s["shiftSd"].is_null(), "shift off");
    assert!(s["cellSd"].is_null(), "no level at finite strength");
    let svg = s["svg"].as_str().unwrap();
    for part in ["nu", "sigma-nu", "sigma-nu-nu", "nu-y"] {
        assert_eq!(
            drawn(svg, part),
            None,
            "{part}: no level at finite strength"
        );
    }
    assert_eq!(drawn(svg, "delta"), None, "shift off");
    assert_eq!(drawn(svg, "sigma-delta"), None);
    assert_eq!(drawn(svg, "delta-p"), None);
    assert_eq!(drawn(svg, "psi"), Some(true), "one gene role");
    assert_eq!(drawn(svg, "beta"), Some(false), "HDP: drawn");
    assert_eq!(drawn(svg, "gamma-beta"), Some(false));
    assert_eq!(drawn(svg, "p"), Some(false));
    assert_eq!(drawn(svg, "strength"), Some(false));
    assert_eq!(drawn(svg, "mu-p"), Some(false));
    assert_eq!(drawn(svg, "p-y"), Some(false));
    assert_eq!(drawn(svg, "mu-y"), None);
    assert!(
        svg.contains("S = 5") && svg.contains("γ = 1"),
        "the run's values"
    );
}

/// One letter, one quantity: `A` the cell roles, `B` the gene roles, `S`
/// the profile strength, in the diagram and the settings line beside it.
#[test]
fn the_strength_and_the_cell_roles_have_their_own_letters() {
    for words in [&["--model", "link"][..], &["--model", "blocks"]] {
        let svg = run(words)["svg"].as_str().unwrap().to_string();
        assert!(!svg.contains("A = 5") && !svg.contains("A = ∞"), "{svg}");
        assert!(svg.contains("S = "), "{words:?}: the strength named S");
    }
    let page = report::FRAGMENT;
    assert!(page.contains(r#"["strength S", "#), "the settings line");
    assert!(
        page.contains("MODEL.cellSd"),
        "the levels' sd in the settings line"
    );
    assert!(!page.contains(r#"["strength A", "#), "the settings line");
}

/// A switch given over its preset's: the run is no preset, and drawn as set.
#[test]
fn a_switch_over_the_preset_leaves_no_preset() {
    let s = run(&["--model", "link", "--shift", "on", "--gene-roles", "3"]);
    assert!(s["preset"].is_null());
    let svg = s["svg"].as_str().unwrap();
    assert_eq!(drawn(svg, "delta-p"), Some(false));
    assert_eq!(drawn(svg, "psi"), Some(false));
    assert_eq!(drawn(svg, "nu"), None, "finite strength: no level");
    // Infinite strength over the link preset brings the cells' levels.
    let s = run(&[
        "--model",
        "link",
        "--shift",
        "on",
        "--profile-strength",
        "inf",
    ]);
    assert!(s["preset"].is_null());
    assert_eq!(s["cellSd"], 1.0);
    assert_eq!(drawn(s["svg"].as_str().unwrap(), "nu-y"), Some(false));
    // A switch given as its preset sets it: still the preset.
    assert_eq!(
        run(&["--model", "blocks", "--shift", "on"])["preset"],
        "blocks"
    );
}

/// Without the run's settings the section is left out, saying why.
#[test]
fn no_record_or_an_older_one_omits_the_section() {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    let s = section(tmp.path());
    assert!(s["svg"].is_null());
    let why = s["why"].as_str().unwrap();
    assert!(
        why.contains("no run record (m6a_blockmodel.run.json)"),
        "{s}"
    );
    let record = tmp.path().join("m6a_blockmodel.run.json");
    std::fs::write(&record, r#"{"job":"cluster","status":"ok","options":{}}"#).unwrap();
    let s = section(tmp.path());
    assert!(
        s["why"]
            .as_str()
            .unwrap()
            .contains("no blockmodel settings"),
        "{s}"
    );
    std::fs::write(&record, "{").unwrap();
    assert!(section(tmp.path())["why"]
        .as_str()
        .unwrap()
        .contains("cannot be read"));
}

/// The record `faba cluster --model crp` leaves in `dir`: its options, no
/// model.
fn crp_record(dir: &Path) {
    let args = cluster(&["--model", "crp"]);
    RunRecord::start("cluster", &dir.to_string_lossy())
        .file_name(&record_name("m6a", args.model))
        .options(&args)
        .finish(&Ok(()));
}

/// A `--model crp` run beside the block tables writes its own record, so
/// the section still draws the blockmodel's.
#[test]
fn a_later_crp_run_leaves_the_model_section() {
    let tmp = tempfile::tempdir().unwrap();
    recorded(tmp.path(), &["--model", "blocks"]);
    crp_record(tmp.path());
    let s = section(tmp.path());
    assert_eq!(s["preset"], "blocks", "{s}");
    assert!(s["svg"].is_string());
}

/// Tables from before the blockmodel had its own record: its settings in
/// the shared `{m}_cluster.run.json` are still read.
#[test]
fn an_older_shared_record_is_still_read() {
    let tmp = tempfile::tempdir().unwrap();
    recorded(tmp.path(), &["--model", "link"]);
    std::fs::rename(
        tmp.path().join("m6a_blockmodel.run.json"),
        tmp.path().join("m6a_cluster.run.json"),
    )
    .unwrap();
    assert_eq!(section(tmp.path())["preset"], "link");
}

/// `words` after `faba modules <dir>`, parsed as it parses them, what was
/// given recorded.
pub(super) fn modules(dir: &Path, words: &[&str]) -> ModulesArgs {
    use clap::{Args, FromArgMatches};
    let dir = dir.to_string_lossy();
    let all = ["faba", &dir].into_iter().chain(words.iter().copied());
    let m = ModulesArgs::augment_args(clap::Command::new("faba"))
        .try_get_matches_from(all)
        .unwrap();
    let mut args = ModulesArgs::from_arg_matches(&m).unwrap();
    args.record_given(&m);
    args
}

/// After a link run over a blocks run, `faba modules` with no model
/// options opens the tables as recorded and names the link preset in its
/// header, as it does given `--model link`; a model or switches that set
/// other settings than the record's fit again, saying why; `--model crp`
/// is refused with the tables there too.
#[test]
fn modules_opens_the_recorded_tables_unless_asked_for_other_settings() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    recorded(dir, &["--model", "blocks"]);
    recorded(dir, &["--model", "link"]);
    let at = |words: &[&str]| plan(&modules(dir, words), dir).unwrap();
    let open = Plan::Open("m6a".to_string());
    assert_eq!(at(&[]), open);
    assert_eq!(at(&["--chains", "2"]), open, "not a model option");
    assert_eq!(at(&["--model", "link"]), open);
    let head = title(dir, &Modules::load(dir, "m6a").unwrap());
    assert!(head.ends_with("m6a · link preset"), "{head}");
    let Plan::Fit(why) = at(&["--model", "blocks"]) else {
        panic!("another preset opened as it is");
    };
    assert!(
        why.contains("link preset") && why.contains("blocks preset"),
        "{why}"
    );
    let Plan::Fit(why) = at(&["--model", "link", "--cell-sd", "0.5", "--shift", "on"]) else {
        panic!("other switches opened as they are");
    };
    assert!(why.contains("custom settings"), "{why}");
    assert!(matches!(at(&["--refit"]), Plan::Fit(_)));
    let refused = plan(&modules(dir, &["--model", "crp"]), dir).unwrap_err();
    assert!(refused.to_string().contains("--model crp"), "{refused}");
    // Custom settings recorded are named so in the header.
    recorded(
        dir,
        &["--model", "link", "--shift", "on", "--gene-roles", "2"],
    );
    let head = title(dir, &Modules::load(dir, "m6a").unwrap());
    assert!(head.ends_with("custom settings"), "{head}");
}

/// Without a blockmodel record the options given cannot be compared: the
/// tables open as they are.
#[test]
fn modules_opens_unrecorded_tables_as_they_are() {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    crp_record(tmp.path());
    let args = modules(tmp.path(), &["--model", "blocks"]);
    assert_eq!(
        plan(&args, tmp.path()).unwrap(),
        Plan::Open("m6a".to_string())
    );
    let head = title(tmp.path(), &Modules::load(tmp.path(), "m6a").unwrap());
    assert!(head.ends_with("settings not recorded"), "{head}");
}

/// The page's script is valid JavaScript (when node is here to say).
#[test]
fn the_page_script_parses() {
    let code = report::FRAGMENT.rsplit_once("<script>").unwrap().1;
    let code = code.rsplit_once("</script>").unwrap().0;
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("page.js");
    std::fs::write(&path, code).unwrap();
    match std::process::Command::new("node")
        .arg("--check")
        .arg(&path)
        .output()
    {
        Ok(out) => assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        ),
        Err(_) => eprintln!("no node: the page's script is not checked"),
    }
}

/// `d` shows and hides the diagram, listed in the footer only with one.
#[test]
fn the_page_toggles_the_diagram_on_d() {
    let page = report::FRAGMENT;
    assert!(page.contains(r#"k === "d""#));
    assert!(page.contains(r#"["d", "diagram"]"#));
    assert!(page.contains(r#"id="model""#) && page.contains(r#"id="model-why""#));
}
