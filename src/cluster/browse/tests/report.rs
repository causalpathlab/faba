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
