//! The report's Model section: the run's settings from its record, its
//! plate diagram drawn as they are set; and what `faba cluster-view` makes of
//! that record: the tables opened as they are, never fitted again.

use super::tables;
use crate::cluster::browse::{
    fit_record, open_tables, report, run_modules, title, Modules, ModulesArgs,
};
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
        .file_name(&record_name("m6a"))
        .options(&args);
    record.set_result("model", &args::recorded(&args));
    record.finish(&Ok(()));
}

/// The page's model section for the tables in `dir`.
fn section(dir: &Path) -> Value {
    let m = Modules::load(dir, "m6a", None).unwrap();
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

/// The blocks preset, as the variational fit records it: its shift drawn,
/// `psi` drawn, `beta` fixed (symmetric), no `p`, and each cell's own
/// level `nu_c` in the cells plate with its prior sd, an edge to the reads.
#[test]
fn the_blocks_preset_is_drawn_as_set() {
    let s = run(&["--cell-roles", "6", "--cell-sd", "0.7"]);
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
    assert!(s["gamma"].is_null(), "the global weights fixed");
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
    assert_eq!(drawn(svg, "beta"), Some(true), "fixed (symmetric)");
    assert!(!svg.contains("β, β′ drawn"), "neither side's weights drawn");
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

/// One letter, one quantity: `A` the cell roles, `B` the gene roles, `S`
/// the profile strength, in the diagram and the settings line beside it.
#[test]
fn the_strength_and_the_cell_roles_have_their_own_letters() {
    let s = run(&[]);
    let svg = s["svg"].as_str().unwrap().to_string();
    assert!(!svg.contains("A = ∞"), "{svg}");
    assert!(svg.contains("S = ∞"), "the strength named S");
    let page = report::FRAGMENT;
    assert!(page.contains(r#"["strength S", "#), "the settings line");
    assert!(
        page.contains("MODEL.cellSd"),
        "the levels' sd in the settings line"
    );
    assert!(!page.contains(r#"["strength A", "#), "the settings line");
}

/// Gene roles other than the preset's: the run is no preset, and drawn as
/// set.
#[test]
fn a_switch_over_the_preset_leaves_no_preset() {
    let s = run(&["--gene-roles", "3"]);
    assert!(s["preset"].is_null());
    assert!(s["svg"].as_str().unwrap().contains("1…3"));
    assert_eq!(run(&["--gene-roles", "8"])["preset"], "blocks");
}

/// Without the run's settings the section is left out, saying why.
#[test]
fn no_record_omits_the_section() {
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

/// A count fit's record (its tables' modality `count`) holds the count
/// model's settings, not the blockmodel's: the section is left out, saying
/// so plainly, and the header names no preset.
#[test]
fn a_count_fit_draws_no_plate_and_claims_no_preset() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    super::count_tables(dir);
    let mut record = RunRecord::start("cluster", &dir.to_string_lossy())
        .file_name("count_blockmodel.run.json")
        .options(&cluster(&["--modality", "count"]));
    record.set_result("model", &crate::cluster::count::Settings::default());
    // As `faba cluster --modality count` records it: the tables' modality
    // says which kind of fit this is.
    record.set_result("tables", &serde_json::json!({"modality": "count"}));
    record.finish(&Ok(()));
    let m = Modules::load(dir, "count", None).unwrap();
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let s = &data["model"];
    assert!(s["svg"].is_null() && s["preset"].is_null(), "{s}");
    let why = s["why"].as_str().unwrap();
    assert!(why.starts_with("a count fit: x_cg ~ Poisson("), "{why}");
    assert!(why.ends_with("no plate drawn yet"), "{why}");
    let head = title(dir, &m);
    assert!(head.ends_with("count · count blockmodel"), "{head}");
    assert!(
        !head.contains("preset") && !head.contains("not recorded"),
        "{head}"
    );
}

/// `words` after `faba cluster-view <dir>`, parsed as `faba` parses them:
/// the arguments, or clap's refusal as `faba` says it.
fn parsed(dir: &Path, words: &[&str]) -> Result<ModulesArgs, String> {
    let dir = dir.to_string_lossy();
    let line = ["faba", "cluster-view", &dir];
    let line = line.into_iter().chain(words.iter().copied());
    let m = crate::matches(&mut crate::faba_command(), line).map_err(|e| e.to_string())?;
    match crate::parse(&m).unwrap().commands {
        crate::Commands::ClusterView(args) => Ok(args),
        _ => panic!("not cluster-view"),
    }
}

/// [`parsed`], given words it takes.
fn modules(dir: &Path, words: &[&str]) -> ModulesArgs {
    parsed(dir, words).unwrap()
}

/// What `faba cluster-view` says of `words`, which it does not take: an
/// option of `faba cluster` (the model's, the reads' filters) refused
/// saying to run `faba cluster` with it, or one of the MCMC fits' saying
/// it is gone, with or without tables to open.
fn refused(dir: &Path, words: &[&str]) {
    let err = parsed(dir, words).unwrap_err();
    let flag = words[0].split('=').next().unwrap();
    let fitting = err.contains("never fits")
        && err.contains("`faba cluster`")
        && err.contains("faba cluster --modality count");
    let removed = crate::cluster::args::is_removed(flag) && err.contains("removed in faba 0.19");
    assert!(
        (fitting || removed) && err.contains(flag),
        "{words:?}: {err}"
    );
}

/// A count fit is opened as it is: `faba cluster-view --modality count`
/// opens its tables, and has no model option, saying to run the command
/// that fits it; with no count tables it says to fit them first.
#[test]
fn a_count_fit_is_opened_as_it_is_never_fitted_again() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let at = |words: &[&str]| {
        let args = modules(dir, words);
        open_tables(dir, args.modality.as_deref())
    };
    let missing = at(&["--modality", "count"]).unwrap_err().to_string();
    assert_eq!(
        missing,
        format!(
            "no count tables in {}: run `faba cluster --modality count` first",
            dir.display()
        )
    );
    super::count_tables(dir);
    assert_eq!(at(&["--modality", "count"]).unwrap(), "count");
    assert_eq!(at(&["--modality", "Count"]).unwrap(), "count");
    for words in [
        &["--chains", "2", "--modality", "count"][..],
        &["--cell-roles", "12", "--modality", "count"],
        &["--count-genes", "500", "--modality", "count"],
        &["--model", "blocks", "--modality", "count"],
    ] {
        refused(dir, words);
    }
    // `--refit` is no option of either command: refused plainly.
    let err = parsed(dir, &["--refit"]).unwrap_err();
    assert!(
        err.contains("--refit") && !err.contains("never fits"),
        "{err}"
    );
}

/// `faba cluster-view` never fits. Given a directory with no block tables
/// it says to run `faba cluster` first; given any option that sets the
/// blockmodel, the recorded settings among them, it refuses, naming `faba
/// cluster`, with or without the tables. Given none, it opens the tables
/// as they are, the recorded preset in its header.
#[test]
fn cluster_view_never_fits() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let view = |words: &[&str]| {
        let all: Vec<&str> = words.iter().copied().chain(["--batch-process"]).collect();
        run_modules(&modules(dir, &all))
    };
    let refusals = [
        &["--model", "link"][..],
        &["--model", "crp"],
        &["--cell-roles", "4"],
        &["--shift", "on"],
        &["--profile-strength", "inf"],
        // A sampler's option sets a fit too: refused, not ignored.
        &["--chains", "2"],
    ];
    let missing = view(&[]).unwrap_err().to_string();
    assert_eq!(
        missing,
        format!(
            "no block tables in {}: run `faba cluster` first",
            dir.display()
        )
    );
    let missing = view(&["--modality", "ATOI"]).unwrap_err().to_string();
    assert_eq!(
        missing,
        format!(
            "no atoi block tables in {}: run `faba cluster` first",
            dir.display()
        )
    );
    refusals.iter().for_each(|w| refused(dir, w));
    recorded(dir, &[]);
    view(&[]).unwrap();
    refusals.iter().for_each(|w| refused(dir, w));
    let head = title(dir, &Modules::load(dir, "m6a", None).unwrap());
    assert!(head.ends_with("m6a · blocks preset"), "{head}");
    // Custom settings recorded are named so in the header.
    recorded(dir, &["--gene-roles", "2"]);
    let head = title(dir, &Modules::load(dir, "m6a", None).unwrap());
    assert!(head.ends_with("custom settings"), "{head}");
}

/// A level asked for by name needs no --batch-process (it is printed as
/// text either way) and is refused beside --html, the page of every level.
#[test]
fn a_level_is_printed_without_batch_process_but_not_with_html() {
    let tmp = tempfile::tempdir().unwrap();
    for level in [
        &["--module", "M1"][..],
        &["--gene", "G1"],
        &["--cells"],
        &["--batches"],
        &["--trace"],
        &["--settings"],
    ] {
        assert!(parsed(tmp.path(), level).is_ok(), "{level:?}");
        let with: Vec<&str> = level.iter().copied().chain(["--batch-mode"]).collect();
        assert!(parsed(tmp.path(), &with).is_ok(), "{with:?}");
        let page: Vec<&str> = level.iter().copied().chain(["--html", "p.html"]).collect();
        let err = parsed(tmp.path(), &page).unwrap_err();
        assert!(err.contains("--html"), "{page:?}: {err}");
    }
    assert!(parsed(tmp.path(), &[]).is_ok());
}

/// A fit's run record says, by its contents, which fit to open: its
/// modality from its `tables` (else from its name), refused for an MCMC
/// fit of faba before 0.19; another command's record is an input as
/// before. It opens
/// alone, never with other inputs, and never fits again.
#[test]
fn a_run_record_names_the_fit_to_open() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let path = |name: &str| -> Box<str> { dir.join(name).to_string_lossy().into() };
    let fit = |p: &str| fit_record(&[p.into()]);
    // No `tables`: the modality in its name.
    recorded(dir, &[]);
    let record = path(&record_name("m6a"));
    let got = fit(&record).unwrap().unwrap();
    assert_eq!((got.dir.as_path(), got.modality.as_str()), (dir, "m6a"));
    // The record's own `tables` say the modality, whatever its name.
    let mut json: Value =
        serde_json::from_str(&std::fs::read_to_string(&*record).unwrap()).unwrap();
    json["results"]["tables"] = crate::cluster::run::sites::tables_entry("m6a");
    std::fs::write(dir.join("renamed.run.json"), json.to_string()).unwrap();
    assert_eq!(
        fit(&path("renamed.run.json")).unwrap().unwrap().modality,
        "m6a"
    );
    // An MCMC fit's record (a crp fit's, a Gibbs-sampled blockmodel's) is
    // refused, by what it holds; another command's is no fit, alone or
    // with other inputs, as before.
    for (name, mcmc) in [
        (
            "m6a_cluster.run.json",
            r#"{"job":"cluster","results":{"posterior":{},"cluster_sizes":[]}}"#,
        ),
        (
            "m6a_old.run.json",
            r#"{"job":"cluster","results":{"blockmodel":{}}}"#,
        ),
    ] {
        std::fs::write(dir.join(name), mcmc).unwrap();
        let err = fit(&path(name)).unwrap_err().to_string();
        assert!(err.contains("MCMC fit of faba before 0.19"), "{err}");
    }
    std::fs::write(dir.join("pipeline.run.json"), r#"{"job":"run"}"#).unwrap();
    assert!(fit(&path("pipeline.run.json")).unwrap().is_none());
    let pipeline_and_more = [path("pipeline.run.json"), "x.zarr".into()];
    assert!(fit_record(&pipeline_and_more).unwrap().is_none());
    let dir_only: Box<str> = dir.to_string_lossy().into();
    assert!(fit(&dir_only).unwrap().is_none());
    assert!(fit_record(&[record.clone(), dir_only]).is_err());
    run_modules(&modules(Path::new(&*record), &["--batch-process"])).unwrap();
    let elsewhere = &["-o", "elsewhere", "--batch-process"];
    let err = run_modules(&modules(Path::new(&*record), elsewhere)).unwrap_err();
    let err = err.to_string();
    assert!(
        err.contains("as it is") && err.contains("`faba cluster`"),
        "{err}"
    );
    for words in [&["--cell-roles", "4"][..], &["--fit", "*_wt_*"]] {
        refused(Path::new(&*record), words);
    }
    let other = run_modules(&modules(
        Path::new(&*record),
        &["--modality", "atoi", "--batch-process"],
    ))
    .unwrap_err();
    assert!(other.to_string().contains("not --modality atoi"), "{other}");
    run_modules(&modules(
        Path::new(&*record),
        &["--modality", "M6A", "--batch-process"],
    ))
    .unwrap();
}

/// Without a blockmodel record the tables open as they are, the header
/// saying their settings are not recorded.
#[test]
fn modules_opens_unrecorded_tables_as_they_are() {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    assert_eq!(open_tables(tmp.path(), None).unwrap(), "m6a");
    let head = title(tmp.path(), &Modules::load(tmp.path(), "m6a", None).unwrap());
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
