//! `cluster` and `cluster_view`: fitting cells by their reads, and
//! viewing a fit.

use crate::mcp::args::{push_flag, push_opt, Arguments};
use crate::mcp::run::{
    absolute, command, exists, figure_prefix, fresh, named, run_out, with_figure,
};
use anyhow::{bail, Result};
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
    push_opt(&mut argv, "--modality", a.str("modality"));
    push_opt(&mut argv, "--fit", a.str("fit"));
    for pair in batch_map(a)? {
        // As one word: a batch named `-x` is not read as a flag.
        argv.push(format!("--batch-map={pair}"));
    }
    argv.extend(a.numbers(&[
        "cell_roles",
        "gene_roles",
        "chains",
        "svi_epochs",
        "batch_size",
        "epochs",
        "tol",
    ])?);
    argv.extend(a.options("cluster")?);
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn cluster_view(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.path("data")?;
    let html = a.path_opt("html")?;
    // What cannot be opened as it is is refused before the command runs;
    // the pages are those `html` writes, a within record's one a fit.
    let pages = crate::cluster::opens(
        &[data.as_str().into()],
        a.str("modality"),
        a.str("group"),
        html.as_deref().map(Path::new),
    )?;
    let mut argv = vec!["cluster-view".into(), data.clone()];
    push_opt(&mut argv, "--modality", a.str("modality"));
    push_opt(&mut argv, "--group", a.str("group"));
    push_opt(&mut argv, "--batch-groups", a.str("batch_groups"));
    for file in ["gaf", "obo", "gmt"] {
        if let Some(path) = a.path_opt(file)? {
            push_flag(&mut argv, &format!("--{file}"), &path);
        }
    }
    argv.extend(a.numbers(&["go_min_weight"])?);
    if let Some(gtf) = a.path_opt("gtf")? {
        push_flag(&mut argv, "--gtf", &gtf);
    }
    // A fitting option is refused by the command itself, naming `faba cluster`.
    argv.extend(a.options("cluster_view")?);
    argv.push("--batch-process".into());
    push_opt(&mut argv, "--module", a.str("module"));
    push_opt(&mut argv, "--gene", a.str("gene"));
    for level in ["cells", "batches", "go", "trace", "settings"] {
        if a.flag(level)? {
            argv.push(format!("--{level}"));
        }
    }
    if let Some(figure) = a.path_opt("figure")? {
        let mut outputs = Vec::new();
        // The prefix checked is the one written under.
        let (prefix, png) = figure_prefix(&figure, &mut outputs)?;
        push_flag(&mut argv, "--figure", &prefix);
        push_opt(&mut argv, "--plot", a.str("plot"));
        push_opt(&mut argv, "--look", a.str("look"));
        let what = format!("a figure of {data}");
        cx.run_making(command(&argv)?, &outputs, &what)?;
        return with_figure(json!({}), Some(&png));
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
    let text = cx.run_reading(command(&argv)?)?;
    Ok(vec![Content::Text(text)])
}

#[cfg(test)]
mod tests {
    use crate::mcp::args::{known_keys, with_typed, Arguments};

    use crate::mcp::schema::tools;

    use super::batch_map;
    use serde_json::json;
    use std::path::{Path, PathBuf};
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
        assert!(opens(at(""), None).unwrap().is_empty());
        let fit = r#"{"job":"cluster","results":{"model":{},"tables":{"modality":"m6a"}}}"#;
        std::fs::write(at("m6a_blockmodel.run.json"), fit).unwrap();
        assert!(opens(at("m6a_blockmodel.run.json"), None)
            .unwrap()
            .is_empty());
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
        assert_eq!(
            opens(None).unwrap(),
            [PathBuf::from("out/fit.html"), "out/fit.C1.html".into()]
        );
        assert_eq!(opens(Some("c1")).unwrap(), [html]);
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
}
