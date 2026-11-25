//! `cluster` and `cluster_view`: fitting cells by their reads, and
//! viewing a fit.

use crate::mcp::args::{push_flag, push_opt, Arguments};
use crate::mcp::run::{absolute, command, exists, figure_prefix, fresh, named, run_out};
use anyhow::{bail, ensure, Context, Result};
use legume_mcp::{Call, Content};
use serde_json::{json, Value};
use std::path::Path;

/// `cluster`'s `batch_map`, as `--batch-map` takes each pair: from an
/// object of `SITES: COUNT`, or a list of `SITES=COUNT`.
fn batch_map(a: &Arguments) -> Result<Vec<String>> {
    match a.0.get("batch_map") {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::Object(m)) => (m.iter())
            .map(|(sites, count)| match count {
                Value::String(c) => Ok(format!("{sites}={c}")),
                Value::Null => Ok(format!("{sites}=")),
                other => bail!("`batch_map` {sites}: the count fit's batch, a string, not {other}"),
            })
            .collect(),
        Some(_) => a.strings("batch_map"),
    }
}

pub(in crate::mcp) fn cluster(a: &Arguments, cx: &Call) -> Result<Value> {
    let data = a.paths("data")?;
    // The second round reads the matrices beside the count fit's when none
    // are named; `--within` itself is among the typed options.
    let within = a.path_opt("within")?;
    let (from, tool, what) = match (&within, data.first()) {
        (Some(record), first) => (
            first.unwrap_or(record),
            "within",
            format!("clusters within the cell modules of {record}"),
        ),
        (None, Some(first)) => (first, "cluster", format!("clusters of {}", data.join(", "))),
        (None, None) => bail!("`data` is required (or `within`, a count fit's run record)"),
    };
    let out = a.out_or(cx, (&named(from, tool), ""), exists)?;
    let mut argv = vec!["cluster".into()];
    argv.extend(data.iter().cloned());
    for pair in batch_map(a)? {
        // As one word: a batch named `-x` is not read as a flag.
        argv.push(format!("--batch-map={pair}"));
    }
    argv.extend(a.options("cluster")?);
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn cluster_view(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.path("data")?;
    let html = a.path_opt("html")?;
    if html.is_some() {
        no_levels_with_html(a)?;
    }
    // What cannot be opened as it is is refused before the command runs;
    // the pages are those `html` writes, a within record's one a fit.
    let opens = crate::cluster::opens(
        &[data.as_str().into()],
        a.str("modality"),
        a.str("group"),
        html.as_deref().map(Path::new),
    )?;
    let pages = opens.pages;
    let mut argv = vec!["cluster-view".into(), data.clone()];
    // The levels, GO's files and cutoffs and the fit opened are among the
    // typed options. A fitting option is refused by the command itself, naming `faba cluster`.
    argv.extend(a.options("cluster_view")?);
    argv.push("--batch-process".into());
    // Several figures may be asked for without a prefix: today's results
    // folder has one.
    let figure = match a.path_opt("figure")? {
        None if a.0.get("plots").is_some_and(|p| !p.is_null()) => {
            let taken = |p: &str| (figures(a, p).into_iter().flatten()).any(|(_, f)| drawn(&f));
            Some(a.out_or(cx, (&named(&data, "cluster_view"), ""), taken)?)
        }
        figure => figure,
    };
    if let Some(figure) = figure {
        return draw_figures(a, cx, &argv, &data, &figure);
    }
    if let Some(html) = html {
        for page in &pages {
            fresh(&page.to_string_lossy())?;
        }
        push_flag(&mut argv, "--html", &html);
        let what = format!("modules of {data} as a page");
        cx.run_making(command(&argv)?, &pages, &what)?;
        let written: Vec<String> = pages
            .iter()
            .map(|p| absolute(&p.to_string_lossy()))
            .collect();
        return Ok(vec![Content::Json(match &written[..] {
            [page] => json!({ "written": page }),
            pages => json!({ "written": pages[0], "pages": pages }),
        })]);
    }
    argv.extend(levels(a, opens.whole));
    let printed = cx.run_reading(command(&argv)?)?;
    let v: Value = serde_json::from_str(&printed)
        .with_context(|| format!("`faba cluster-view --json` printed no JSON: {printed}"))?;
    Ok(vec![Content::Json(v)])
}

/// The levels `cluster_view` takes, each among the typed options.
const LEVELS: [&str; 9] = [
    "overview", "grid", "module", "gene", "cells", "batches", "go", "trace", "settings",
];

/// Whether `a` asks for `level`: a switch set, or a module or gene named.
fn asked(a: &Arguments, level: &str) -> bool {
    !matches!(
        a.0.get(level),
        None | Some(Value::Null | Value::Bool(false))
    )
}

/// Refuses a level with `html`, whose page shows every level, before the
/// command runs, so no level flag reaches the page's command line.
fn no_levels_with_html(a: &Arguments) -> Result<()> {
    let given: Vec<String> = (LEVELS.into_iter())
        .filter(|level| asked(a, level))
        .map(|level| format!("`{level}`"))
        .collect();
    ensure!(
        given.is_empty(),
        "{} cannot be given with `html`: the page shows every level; leave the level \
         arguments out, or drop `html` to get them as JSON",
        given.join(", ")
    );
    Ok(())
}

/// `--json`, so the levels asked for are printed as JSON; when none is,
/// with the overview of one fit, or alone for a within record's listing
/// when it is opened `whole`. GO is read only for `go`.
fn levels(a: &Arguments, whole: bool) -> Vec<String> {
    let mut argv = vec!["--json".to_string()];
    if !whole && !LEVELS.into_iter().any(|level| asked(a, level)) {
        argv.push("--overview".into());
    }
    argv
}

/// The figures asked for (`plots`, else `plot`, else the modules grid),
/// each with the prefix it is written under: `figure` for one, and
/// `{figure}.{plot}` each for several. `plot` and `plots` together, or a
/// plot named twice, are refused.
fn figures(a: &Arguments, figure: &str) -> Result<Vec<(String, String)>> {
    let plots = a.strings("plots")?;
    ensure!(
        plots.is_empty() || a.str("plot").is_none(),
        "give `plot` or `plots`, not both"
    );
    let plots = match (plots.is_empty(), a.str("plot")) {
        (false, _) => plots,
        (true, plot) => vec![plot.unwrap_or("modules").to_string()],
    };
    for (i, p) in plots.iter().enumerate() {
        ensure!(!plots[..i].contains(p), "`plots` names {p} twice");
    }
    Ok(match &plots[..] {
        [one] => vec![(one.clone(), figure.to_string())],
        several => (several.iter())
            .map(|p| (p.clone(), format!("{figure}.{p}")))
            .collect(),
    })
}

/// Each figure `a` asks for of `data` under `figure` (see [`figures`]),
/// drawn by `argv` with its `--figure` and `--plot`: the paths written,
/// and each as an image. Every file is checked new before the first is
/// drawn.
fn draw_figures(
    a: &Arguments,
    cx: &Call,
    argv: &[String],
    data: &str,
    figure: &str,
) -> Result<Vec<Content>> {
    let mut plans = Vec::new();
    for (plot, prefix) in figures(a, figure)? {
        // The prefix checked is the one written under.
        let mut made = Vec::new();
        let (prefix, png) = figure_prefix(&prefix, &mut made)?;
        plans.push((plot, prefix, png, made));
    }
    let (mut written, mut images) = (Vec::new(), Vec::new());
    for (plot, prefix, png, made) in plans {
        let mut argv = argv.to_vec();
        push_flag(&mut argv, "--figure", &prefix);
        push_flag(&mut argv, "--plot", &plot);
        push_opt(&mut argv, "--look", a.str("look"));
        let what = format!("the {plot} figure of {data}");
        cx.run_making(command(&argv)?, &made, &what)?;
        let pdf = absolute(&made[0].to_string_lossy());
        written.push(json!({ "plot": plot, "pdf": pdf, "png": absolute(&png) }));
        let image = std::fs::read(&png).with_context(|| format!("reading {png}"))?;
        images.push(Content::Png(image));
    }
    let listed = Content::Json(json!({ "figures": written }));
    Ok(std::iter::once(listed).chain(images).collect())
}

/// Whether a file a figure saved under `prefix` writes is there.
fn drawn(prefix: &str) -> bool {
    crate::figure::saved_paths(prefix).iter().any(|p| exists(p))
}

#[cfg(test)]
mod tests {
    use crate::mcp::args::{known_keys, typed_flags, with_typed, Arguments};

    use crate::mcp::schema::tools;

    use super::{batch_map, figures, levels, no_levels_with_html};
    use crate::mcp::schema::properties;
    use serde_json::json;
    use std::path::{Path, PathBuf};

    /// `cluster_view` has `faba cluster-view` print JSON, the overview
    /// when no level is named (a within record opened whole: its listing,
    /// no level at all); `json` is the tool's, not offered and refused in
    /// `options`. Every level and the figures' look are arguments, the
    /// levels as their flags; what the tool asks for parses as the command
    /// does.
    #[test]
    fn cluster_view_asks_for_json_and_the_overview_unless_a_level_is_named() {
        let argv = |args: serde_json::Value| {
            let a = Arguments(args);
            let mut argv = typed_flags("cluster_view", &a).unwrap();
            argv.extend(levels(&a, false));
            argv
        };
        let whole = |args| levels(&Arguments(args), true);
        assert_eq!(argv(json!({})), ["--json", "--overview"]);
        assert_eq!(argv(json!({ "grid": false })), ["--json", "--overview"]);
        assert_eq!(whole(json!({})), ["--json"]);
        assert_eq!(whole(json!({ "cells": true })), ["--json"]);
        assert_eq!(
            argv(json!({ "module": "M2" })),
            ["--module", "M2", "--json"]
        );
        assert_eq!(
            argv(json!({ "go": true, "grid": true, "overview": true })),
            ["--overview", "--grid", "--go", "--json"]
        );
        for args in [
            json!({}),
            json!({ "gene": "GENE1", "trace": true, "settings": true }),
        ] {
            let line = ["faba", "cluster-view", "out/", "--batch-process"];
            let all = line.iter().copied().map(String::from).chain(argv(args));
            assert!(crate::matches(&mut crate::faba_command(), all).is_ok());
        }
        let props = properties("cluster_view").unwrap();
        assert!(props.get("json").is_none(), "the tool sets --json itself");
        for key in [
            "overview",
            "grid",
            "plots",
            "font_size",
            "xlab",
            "ylab",
            "panel_only",
            "palette",
            "width",
            "height",
        ] {
            assert!(props.contains_key(key), "cluster_view: no `{key}`");
        }
        let a = Arguments(json!({ "options": ["--json"] }));
        let err = a.options("cluster_view").unwrap_err().to_string();
        assert!(err.contains("adds it itself"), "{err}");
    }

    /// With `html` every level asked for is refused by name, saying the
    /// page shows them all; a level switched off or null is not asked for.
    #[test]
    fn cluster_view_refuses_a_level_with_html() {
        let check = |args| no_levels_with_html(&Arguments(args)).map_err(|e| e.to_string());
        assert!(check(json!({ "html": "p.html" })).is_ok());
        assert!(check(json!({ "go": false, "module": null })).is_ok());
        let err = check(json!({ "go": true, "gene": "GENE1" })).unwrap_err();
        assert!(
            err.starts_with("`gene`, `go` cannot be given with `html`"),
            "{err}"
        );
        assert!(err.contains("the page shows every level"), "{err}");
    }

    /// One figure is written under `figure`, several each under
    /// `{figure}.{plot}`; `plot` and `plots` together, or a plot twice,
    /// are refused.
    #[test]
    fn several_plots_are_written_each_under_its_own_prefix() {
        let plans = |args| figures(&Arguments(args), "out/fig").map_err(|e| e.to_string());
        let plan = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
            (pairs.iter())
                .map(|(p, f)| (p.to_string(), f.to_string()))
                .collect()
        };
        assert_eq!(plans(json!({})).unwrap(), plan(&[("modules", "out/fig")]));
        assert_eq!(
            plans(json!({ "plot": "trace" })).unwrap(),
            plan(&[("trace", "out/fig")])
        );
        assert_eq!(
            plans(json!({ "plots": ["modules", "trace", "batches"] })).unwrap(),
            plan(&[
                ("modules", "out/fig.modules"),
                ("trace", "out/fig.trace"),
                ("batches", "out/fig.batches"),
            ])
        );
        let both = plans(json!({ "plot": "trace", "plots": ["modules"] }));
        assert!(both.unwrap_err().contains("not both"));
        let twice = plans(json!({ "plots": ["trace", "trace"] }));
        assert!(twice.unwrap_err().contains("twice"));
    }
    /// Every blockmodel option, under any of its flags, passes through
    /// `options` to `faba cluster-view`, which refuses it as clap does,
    /// saying to run `faba cluster` with it (a viewer never fits), or, for
    /// the MCMC fits' flags, that they are gone.
    #[test]
    fn cluster_view_refuses_every_option_that_would_fit() {
        let removed = [
            &["--model", "link"][..],
            &["--prior-strength", "2"],
            &["--shift", "on"],
            &["--global", "on"],
            &["--profile-strength", "inf"],
            &["--community-gamma", "3"],
            &["--sweeps", "50"],
            &["--burn-in", "10"],
        ];
        for option in [
            &["--cell-roles", "4"][..],
            &["--communities", "4"],
            &["--gene-roles", "3"],
            &["--membership-alpha=2"],
            &["--gene-alpha", "1"],
            &["--shift-sd", "1"],
            &["--cell-sd", "1"],
            &["--chains", "2"],
            &["--epochs", "5"],
        ]
        .into_iter()
        .chain(removed)
        {
            let a = Arguments(json!({ "options": option }));
            let passed = a.options("cluster_view").unwrap();
            let line = ["faba", "cluster-view", "out/", "--batch-process"];
            let all = line.iter().copied().map(String::from).chain(passed);
            let err = crate::matches(&mut crate::faba_command(), all)
                .unwrap_err()
                .to_string();
            let said = match removed.contains(&option) {
                true => err.contains("removed in faba 0.19"),
                false => err.contains("never fits") && err.contains("`faba cluster`"),
            };
            assert!(said, "{option:?}: {err}");
        }
    }

    /// `cluster_view` opens a fit's record or its tables, never anything
    /// it would have to fit first: another command's record or site
    /// matrices without tables beside them, or a directory without tables,
    /// are refused.
    #[test]
    fn cluster_view_opens_only_what_needs_no_fit() {
        let tmp = tempfile::tempdir().unwrap();
        let d = tmp.path();
        let at = |name: &str| d.join(name).to_string_lossy().into_owned();
        let opens = |data: String, group| crate::cluster::opens(&[data.into()], None, group, None);
        std::fs::write(at("pipeline.run.json"), r#"{"job":"run"}"#).unwrap();
        std::fs::write(at("x.zarr"), "").unwrap();
        for refused in [at("pipeline.run.json"), at("x.zarr"), at("")] {
            let err = opens(refused.clone(), None).unwrap_err().to_string();
            assert!(err.contains("run `faba cluster` first"), "{refused}: {err}");
        }
        std::fs::write(at("m6a_block_sites.parquet"), "").unwrap();
        assert!(opens(at(""), None).unwrap().pages.is_empty());
        let fit = r#"{"job":"cluster","results":{"model":{},"tables":{"modality":"m6a"}}}"#;
        std::fs::write(at("m6a_blockmodel.run.json"), fit).unwrap();
        let one = opens(at("m6a_blockmodel.run.json"), None).unwrap();
        assert!(one.pages.is_empty() && !one.whole, "{one:?}");
        // One fit has no groups to name.
        let err = opens(at("m6a_blockmodel.run.json"), Some("C1")).unwrap_err();
        assert!(err.to_string().contains("within record"), "{err}");
    }

    /// A within record opened whole names its fitted groups, whose pages
    /// `html` writes beside the count fit's; with `group` it is one fit,
    /// refused when that group has none.
    #[test]
    fn cluster_view_opens_a_within_record_whole_or_by_group() {
        let tmp = tempfile::tempdir().unwrap();
        let d = tmp.path();
        let sub = d.join("within/m6a/C1");
        std::fs::create_dir_all(&sub).unwrap();
        let fit = r#"{"job":"cluster","results":{"model":{},"tables":{"modality":"m6a"}}}"#;
        std::fs::write(sub.join("m6a_blockmodel.run.json"), fit).unwrap();
        std::fs::write(sub.join("m6a_block_sites.parquet"), b"").unwrap();
        let record = json!({ "job": "cluster", "results": { "within": {
            "modality": "m6a", "count_record": d.join("count_blockmodel.run.json"),
            "groups": [
                { "group": "C1", "cells": 80, "status": "fitted", "dir": "within/m6a/C1",
                  "record": "m6a_blockmodel.run.json" },
                { "group": "C2", "cells": 9, "status": "skipped", "why": "too few" },
            ],
        } } });
        let path = d.join("m6a_within.run.json").to_string_lossy().into_owned();
        std::fs::write(&path, record.to_string()).unwrap();
        let html = Path::new("out/fit.html");
        let opens = |group| crate::cluster::opens(&[path.as_str().into()], None, group, Some(html));
        let whole = opens(None).unwrap();
        assert!(whole.whole);
        assert_eq!(
            whole.pages,
            [PathBuf::from("out/fit.html"), "out/fit.C1.html".into()]
        );
        let one = opens(Some("c1")).unwrap();
        assert!(!one.whole);
        assert_eq!(one.pages, [html]);
        let err = opens(Some("C2")).unwrap_err().to_string();
        assert!(err.contains("C2 has no fit to open"), "{err}");
        let err = opens(Some("C9")).unwrap_err().to_string();
        assert!(err.contains("no group C9"), "{err}");
    }

    /// The second round needs no `data`: its matrices are those beside the
    /// count fit's, and its output is named after the count fit's record.
    #[test]
    fn cluster_within_a_count_fit_needs_no_data() {
        let args = json!({ "within": "out/count_blockmodel.run.json", "within_min_cells": 20 });
        let a = Arguments(with_typed("cluster", args).unwrap());
        assert_eq!(
            a.options("cluster").unwrap(),
            [
                "--within",
                "out/count_blockmodel.run.json",
                "--within-min-cells",
                "20"
            ]
        );
        assert!(tools().iter().any(|t| {
            t.name == "cluster" && t.input_schema["properties"]["within"].is_object()
        }));
    }

    /// `batch_map` is an object of `SITES: COUNT` or a list of
    /// `SITES=COUNT`, each pair one `--batch-map`; anything else is refused.
    #[test]
    fn cluster_batch_map_is_an_object_or_a_list() {
        let map = |v| batch_map(&Arguments(json!({ "batch_map": v })));
        assert_eq!(
            map(json!({ "/x/solo": "solo", "b": "" })).unwrap(),
            ["/x/solo=solo", "b="]
        );
        assert_eq!(map(json!(["-b=c/b"])).unwrap(), ["-b=c/b"]);
        assert!(map(json!({ "b": 3 })).is_err());
        assert!(map(json!(7)).is_err());
        assert!(batch_map(&Arguments(json!({}))).unwrap().is_empty());
        let cluster = tools().into_iter().find(|t| t.name == "cluster").unwrap();
        assert!(cluster.input_schema["properties"]["batch_map"].is_object());
        assert!(cluster.description.contains("`batch_map`"));
    }

    /// The MCMC fits' arguments are refused by name, `model` blocks taken;
    /// the fit's own are not touched.
    #[test]
    fn cluster_refuses_the_mcmc_arguments() {
        let refused = |args| known_keys("cluster", &args).map_err(|e| e.to_string());
        let err = refused(json!({ "data": ["d"], "sweeps": 100, "model": "link", "thin": 2 }))
            .unwrap_err();
        assert!(err.starts_with("`model`, `sweeps`, `thin` "), "{err}");
        assert!(err.contains("removed in faba 0.19"), "{err}");
        assert!(refused(json!({ "data": ["d"], "model": "blocks", "epochs": 5 })).is_ok());
        assert!(refused(json!({ "profile_strength": 5 })).is_err());
        assert!(
            refused(json!({ "sweeps": null, "holdout": null })).is_ok(),
            "null: not given"
        );
    }

    /// `go_q` is in `cluster_view`'s schema as a number in (0, 1], and
    /// reaches `faba cluster-view` as the `--go-q` it parses.
    #[test]
    fn cluster_view_takes_the_go_q_cutoff() {
        let specs = tools();
        let spec = specs.iter().find(|s| s.name == "cluster_view").unwrap();
        let q = &spec.input_schema["properties"]["go_q"];
        assert_eq!(q["type"], "number", "{q}");
        assert_eq!(
            (q["exclusiveMinimum"].as_f64(), q["maximum"].as_f64()),
            (Some(0.0), Some(1.0))
        );
        let flags = typed_flags("cluster_view", &Arguments(json!({ "go_q": 0.1 }))).unwrap();
        let argv = ["faba", "cluster-view", "fit.run.json"]
            .into_iter()
            .map(String::from)
            .chain(flags);
        let cli = <crate::Cli as clap::Parser>::try_parse_from(argv).unwrap();
        let crate::Commands::ClusterView(args) = cli.commands else {
            panic!("not cluster-view");
        };
        assert_eq!(args.shown.go_ask().q, 0.1);
    }
}
