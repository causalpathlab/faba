//! The HTML report.

use super::loaded;
use crate::cluster::browse::{gene_models, report};

#[test]
fn the_report_holds_the_modules_and_loads_nothing() {
    let m = loaded();
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("report.html");
    let model = crate::site_analysis::miami::genemodel::GeneModel {
        chr: "chr1".into(),
        lo: 50,
        hi: 401,
        forward: true,
        exons: vec![(50, 120), (280, 401)],
        symbol: "GENE1".into(),
        key: "GENE1".into(),
    };
    let models = gene_models(&m, &[model]);
    assert_eq!(models.len(), 1, "matched by key on its chromosome");
    report::write(&m, &models, &path).unwrap();
    let html = std::fs::read_to_string(&path).unwrap();
    assert!(!html.contains("http://") && !html.contains("https://"));
    let at = html
        .find(r#"<script id="data" type="application/json">"#)
        .unwrap();
    let json = &html[at..];
    let json = &json[json.find('>').unwrap() + 1..json.find("</script>").unwrap()];
    let data: serde_json::Value = serde_json::from_str(json).unwrap();
    let modules = data["modules"].as_array().unwrap();
    assert_eq!(modules[0]["name"], "M2", "by shift range");
    // GENE1 on its model: positions within it, exons from its start.
    let gene1 = &modules[1]["genes"][0];
    assert_eq!(
        (gene1["name"].as_str(), gene1["lo"].as_i64()),
        (Some("GENE1"), Some(50))
    );
    assert_eq!(gene1["exons"][0], serde_json::json!([0, 70]));
    let site = &gene1["sites"][0];
    let on_chr = site["pos"]
        .as_i64()
        .zip(gene1["lo"].as_i64())
        .map(|(p, lo)| p + lo);
    assert_eq!((site["pos"].as_i64(), on_chr), (Some(50), Some(100)));
    assert_eq!(site["c"][0], serde_json::json!([8.0, 10.0]));
    assert_eq!(site["s"], serde_json::json!([20, 20, 9]));
    assert_eq!(data["on"], "methylated");
}

/// The page's data for `m`, GENE1 on a model.
fn data_json(m: &crate::cluster::browse::Modules) -> String {
    let model = crate::site_analysis::miami::genemodel::GeneModel {
        chr: "chr1".into(),
        lo: 50,
        hi: 401,
        forward: true,
        exons: vec![(50, 120), (280, 401)],
        symbol: "GENE1".into(),
        key: "GENE1".into(),
    };
    let models = gene_models(m, &[model]);
    serde_json::to_string(&report::data(m, &models)).unwrap()
}

/// The page's data is in the panes' column order: the same tables with
/// their cell roles numbered the other way give the same data, only the
/// cell modules' names following their roles.
#[test]
fn the_reports_data_is_in_column_order() {
    let head = r#"{"title":"gene modules × cell modules · m6a","modality":"m6a","on":"methylated","off":"unmethylated","cellNames":"#;
    let tail = r#","cellLabels":["1.6","1.4"],"allLabel":"3.0","diagonal":["M2","M1"],"modules":[{"name":"M2","size":1,"range":1.5,"genes":[{"name":"GENE3","chr":"chr1","lo":-40,"len":100,"exons":[],"sites":[{"pos":50,"c":[[30.0,40.0],[10.0,40.0]],"s":[40,40,12]}]}]},{"name":"M1","size":2,"range":0.2,"genes":[{"name":"GENE1","chr":"chr1","lo":50,"len":350,"exons":[[0,70],[230,351]],"sites":[{"pos":50,"c":[[8.0,10.0],[1.0,10.0]],"s":[20,20,9]},{"pos":250,"c":[[5.0,10.0],[0.0,0.0]],"s":[10,10,4]}]},{"name":"GENE2","chr":"chr1","lo":0,"len":100,"exons":[],"sites":[{"pos":50,"c":[[1.0,2.0],[0.0,0.0]],"w":[[1.0,2.0],[3.0,4.0]],"s":[4,2,3]}]}]}],"model":{"why":"no run record (m6a_blockmodel.run.json) beside the tables"},"trace":{"why":"no m6a_trace.parquet beside the tables: older fits kept no trace (rerun with --refit)"}}"#;
    for (m, names) in [loaded(), super::swapped()]
        .iter()
        .zip([r#"["C1","C2"]"#, r#"["C2","C1"]"#])
    {
        assert_eq!(data_json(m), format!("{head}{names}{tail}"));
    }
}

/// The page written holds its title, the browser's fragment, and the data
/// as [`report::data`] makes it.
#[test]
fn the_page_holds_its_title_fragment_and_data() {
    for m in [loaded(), super::swapped()] {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("report.html");
        let models = Default::default();
        report::write(&m, &models, &path).unwrap();
        let html = std::fs::read_to_string(&path).unwrap();
        assert!(html.contains(&format!("<title>faba modules · {}</title>", m.modality)));
        let script = report::FRAGMENT
            .rsplit_once("<script>")
            .map(|(_, code)| code)
            .unwrap();
        assert!(html.contains(script), "the fragment's code as it is");
        let at = html
            .find(r#"<script id="data" type="application/json">"#)
            .unwrap();
        let json = &html[at..];
        let json = &json[json.find('>').unwrap() + 1..json.find("</script>").unwrap()];
        let back: serde_json::Value = serde_json::from_str(json).unwrap();
        let data = serde_json::to_value(report::data(&m, &models)).unwrap();
        assert_eq!(back, data);
    }
}

/// The page's trace: each chain's points as saved, in its terminal colour,
/// the burn-in and sweeps to mark, and the verdict's words and advice.
#[test]
fn the_report_draws_the_trace_with_its_verdict() {
    let m = super::trace::load(1.31, true);
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let t = &data["trace"];
    let chains = t["chains"].as_array().unwrap();
    assert_eq!(chains.len(), 2);
    let saved = super::trace::chains([0.0, 0.5]);
    for (i, c) in chains.iter().enumerate() {
        assert_eq!(c["name"], format!("chain {i}"));
        assert!(c["color"].as_str().unwrap().starts_with('#'), "{c}");
        let points: Vec<(usize, f64)> = serde_json::from_value(c["points"].clone()).unwrap();
        assert_eq!(points.len(), saved[i].len());
        for (p, q) in points.iter().zip(&saved[i]) {
            assert!(p.0 == q.0 && (p.1 - q.1).abs() < 1e-3, "{p:?} {q:?}");
        }
    }
    assert_ne!(chains[0]["color"], chains[1]["color"]);
    assert_eq!(
        (t["sweeps"].as_u64(), t["burnIn"].as_u64()),
        (Some(200), Some(20))
    );
    assert_eq!(t["rhat"].as_f64(), Some(1.31));
    assert_eq!(t["converged"], false);
    assert_eq!(t["said"], "R-hat 1.31: chains disagree");
    assert_eq!(
        t["advice"],
        "run longer (--sweeps 400) and check the traces"
    );
}

/// Older tables kept no trace: the section says why instead.
#[test]
fn without_a_trace_the_report_says_why() {
    let m = super::trace::load(1.02, false);
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let why = data["trace"]["why"].as_str().unwrap();
    assert!(why.contains("no m6a_trace.parquet"), "{why}");
    assert!(data["trace"].get("chains").is_none());
}

/// The Model section's own parts (heading, settings, diagram, legend) sit
/// in one block, before the trace, that the page shows only with a
/// diagram: a trace without one shows no empty box above it.
#[test]
fn the_models_parts_hide_together_apart_from_the_trace() {
    let page = report::FRAGMENT;
    let body = page
        .find(r#"<div id="model-body" hidden>"#)
        .expect("one block");
    let trace = page.find(r#"<div class="trace" id="trace""#).unwrap();
    for part in [
        "model-toggle",
        "model-settings",
        "model-diagram",
        "model-legend",
    ] {
        let at = page.find(&format!(r#"id="{part}""#)).unwrap();
        assert!(body < at && at < trace, "{part} in the block");
    }
    assert!(page.contains(r#"getElementById("model-body").hidden = false"#));
}
