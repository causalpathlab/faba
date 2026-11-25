//! The levels as JSON (`--json`): each level an object under its name, its
//! rows objects with named fields; the overview, which reads no GO unless
//! asked; GO's terms as the text table lists them.

use super::super::json::levels;
use super::super::text::{printed, Shown};
use super::batches::batched;
use super::trace::load;
use super::{counted, loaded};
use serde_json::{json, Value};

fn shown(f: impl FnOnce(&mut Shown)) -> Shown {
    let mut s = Shown::default();
    f(&mut s);
    s
}

/// `v[key]` as an array, saying what was there when it is not one.
fn rows<'a>(v: &'a Value, key: &str) -> &'a Vec<Value> {
    v[key]
        .as_array()
        .unwrap_or_else(|| panic!("no `{key}` rows in {v:#}"))
}

/// Each of `fields` is in `row`.
fn has(row: &Value, fields: &[&str]) {
    for f in fields {
        assert!(row.get(f).is_some(), "no `{f}` in {row:#}");
    }
}

/// No level asked: the modules grid, a row a gene module with a value and
/// reads in each cell module, as the text's first level has them.
#[test]
fn the_modules_grid_comes_back_as_rows() {
    let m = loaded();
    let v = levels(&m, &Shown::default()).unwrap();
    assert_eq!(v.as_object().unwrap().len(), 1, "{v:#}");
    let grid = &v["modules"];
    assert_eq!(
        (grid["unit"].as_str(), grid["total"].as_str()),
        (Some("%"), Some("reads"))
    );
    let text = printed(&m, &Shown::default()).unwrap();
    let rows = rows(grid, "rows");
    assert_eq!(rows.len(), text.len() - 1);
    for (row, line) in rows.iter().zip(&text[1..]) {
        has(row, &["module", "main_genes", "range", "values"]);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(row["module"], cols[0]);
        assert_eq!(row["main_genes"].to_string(), cols[1]);
        let values = row["values"].as_array().unwrap();
        assert_eq!(values.len(), m.cell_roles);
        for (k, cell) in values.iter().enumerate() {
            has(cell, &["cell_module", "value", "reads"]);
            assert_eq!(cell["cell_module"], format!("C{}", m.columns[k] + 1));
            let reads = cell["reads"].as_f64().unwrap();
            assert_eq!(format!("{reads:.1}"), cols[4 + 2 * k]);
        }
    }
    // Asked for beside another level, as `grid`.
    let both = levels(&m, &shown(|s| (s.grid, s.cells) = (true, true))).unwrap();
    assert_eq!(both["modules"], *grid);
    assert!(both["cells"].is_object());
}

/// Every level asked for comes back under its name, its rows carrying
/// their fields.
#[test]
fn each_level_comes_back_as_rows_with_its_fields() {
    let m = loaded();
    let v = levels(
        &m,
        &shown(|s| {
            s.cells = true;
            s.module = Some("M1".into());
            s.gene = Some("GENE1".into());
            s.trace = true;
            s.settings = true;
        }),
    )
    .unwrap();
    for (row, size) in rows(&v["cells"], "rows").iter().zip(&m.columns) {
        has(row, &["module", "cells"]);
        assert_eq!(row["cells"].as_f64(), Some(m.cell_sizes[*size]));
    }
    let module = &v["module"];
    assert_eq!(module["module"], "M1");
    has(module, &["main_genes", "sites", "range", "unit", "total"]);
    for row in rows(module, "rows") {
        has(row, &["gene", "site", "values"]);
        assert!(row["site"].as_str().unwrap().starts_with("chr1:"), "{row}");
    }
    let gene = &v["gene"];
    assert_eq!(
        (gene["gene"].as_str(), gene["main_module"].as_str()),
        (Some("GENE1"), Some("M1"))
    );
    assert_eq!(gene["chr"], "chr1");
    let sites = rows(gene, "rows");
    assert_eq!(sites.len(), 2);
    for row in sites {
        has(
            row,
            &["site", "values", "converted", "unconverted", "cells"],
        );
    }
    // No record beside the fixture's tables: no trace, no verdict, saying why.
    let trace = &v["trace"];
    assert!(
        trace["restarts"].as_array().unwrap().is_empty(),
        "{trace:#}"
    );
    assert!(trace["no_trace"].as_str().is_some() && trace["verdict"].is_null());
    assert!(trace["unjudged"].as_str().is_some(), "{trace:#}");
    assert!(v["settings"]["model"].as_str().is_some(), "{v:#}");

    // A fit with a record: each restart's points and the verdict.
    let traced = load(false, true);
    let t = &levels(&traced, &shown(|s| s.trace = true)).unwrap()["trace"];
    let restarts = rows(t, "restarts");
    assert_eq!(restarts.len(), 2);
    for point in restarts[0].as_array().unwrap() {
        has(point, &["step", "elbo"]);
    }
    assert_eq!(
        (t["verdict"].as_str(), t["converged"].as_bool()),
        (Some("ELBO check"), Some(false))
    );
    has(t, &["said", "advice", "steps", "warmup"]);

    // A count fit's values are UMIs per cell, beside each module's cells.
    let c = counted();
    let g = &levels(&c, &shown(|s| s.gene = Some("G3".into()))).unwrap()["gene"];
    assert_eq!(
        (g["unit"].as_str(), g["total"].as_str()),
        (Some("UMI/cell"), Some("cells"))
    );
    let row = &rows(g, "rows")[0];
    assert!(
        row["site"].is_null(),
        "a count fit's gene has no place: {row}"
    );
    let at = |cell: usize| {
        let k = c.columns.iter().position(|&x| x == cell).unwrap();
        row["values"][k]["value"].as_f64().unwrap()
    };
    assert_eq!((at(0), at(1)), (0.0, 3.0), "{row}");
}

/// The batches: a row a cell module, its label, evenness and share in each
/// batch by name; the conditions; one batch, said.
#[test]
fn the_batches_come_back_as_rows() {
    let m = batched();
    let b = &levels(&m, &shown(|s| s.batches = true)).unwrap()["batches"];
    assert_eq!(
        b["batches"],
        json!(["rep1_a", "rep2_a", "rep1_b", "rep2_b"])
    );
    let rows = rows(b, "rows");
    assert_eq!(rows.len(), m.cell_roles);
    for row in rows {
        has(row, &["module", "label", "evenness", "shares"]);
        assert_eq!(row["shares"].as_object().unwrap().len(), 4, "{row}");
    }
    assert!(
        rows.iter().any(|r| r["label"] == "batch: rep1_a"),
        "{rows:?}"
    );
    assert!(
        b["conditions"].as_array().is_some_and(|c| !c.is_empty()),
        "{b:#}"
    );
    let one = &levels(&loaded(), &shown(|s| s.batches = true)).unwrap()["batches"];
    assert!(
        one["batches"].is_null() && one["note"].as_str().is_some(),
        "{one:#}"
    );
}

/// The overview: the verdict, how many gene and cell modules, each cell
/// module's size and each gene module's size and top genes by weight. It
/// reads no GO: a GAF that is not there is never opened, which `go` would.
#[test]
fn the_overview_says_the_fit_at_a_glance_and_reads_no_go() {
    let m = load(true, true);
    let missing = shown(|s| {
        s.overview = true;
        s.gaf = Some("/no/such/annotations.gaf".into());
        s.obo = Some("/no/such/go.obo".into());
    });
    let v = levels(&m, &missing).unwrap();
    assert_eq!(v.as_object().unwrap().len(), 1, "the overview alone: {v:#}");
    let o = &v["overview"];
    assert_eq!(o["verdict"], "ELBO ok");
    assert_eq!(
        (o["gene_modules"].as_u64(), o["cell_modules"].as_u64()),
        (Some(3), Some(2))
    );
    let cells = rows(o, "cell_sizes");
    assert_eq!(cells.len(), 2);
    for row in cells {
        has(row, &["module", "cells"]);
    }
    let genes = rows(o, "gene_sizes");
    assert_eq!(genes.len(), 3);
    // Gene module 1: GENE1 (0.9) before GENE2 (0.6), its main genes.
    assert_eq!(genes[0]["module"], "M1");
    assert_eq!(genes[0]["genes"], 2);
    let top: Vec<&str> = (genes[0]["top"].as_array().unwrap().iter())
        .map(|g| g["gene"].as_str().unwrap())
        .collect();
    assert_eq!(top, ["GENE1", "GENE2"]);
    let w = genes[0]["top"][0]["weight"].as_f64().unwrap();
    assert!((w - 0.9).abs() < 1e-6, "{w}");
    // As text, alike: no GO read.
    let text = printed(&m, &missing).unwrap();
    assert!(
        text[0].starts_with("overview: m6a fit, ELBO ok"),
        "{text:?}"
    );
    // With `go`, the GAF is read, and its absence refused.
    let go = shown(|s| {
        s.overview = true;
        s.go = true;
        s.gaf = missing.gaf.clone();
        s.obo = missing.obo.clone();
    });
    let e = levels(&m, &go).unwrap_err().to_string();
    assert!(e.contains("/no/such/"), "{e}");
}

/// GO's terms: per gene module the terms the text table lists, each with
/// its id, name, namespace, counts, odds, p, q and genes, the q cutoff and
/// weight used, and why a module lists none.
#[test]
fn go_comes_back_as_the_text_tables_terms() {
    let mut m = loaded();
    m.gene_main = (0..40)
        .map(|i| {
            (
                format!("ENSG{i:011}_G{i}").into(),
                usize::from(i >= 20),
                0.9,
            )
        })
        .collect();
    let tmp = tempfile::tempdir().unwrap();
    let set = |name: &str, genes: std::ops::Range<usize>| {
        let g: Vec<String> = genes.map(|i| format!("G{i}")).collect();
        format!("{name}\t-\t{}", g.join("\t"))
    };
    let gmt = tmp.path().join("sets.gmt");
    let sets = [set("UP", 0..12), set("REST", 12..40)];
    std::fs::write(&gmt, sets.join("\n")).unwrap();
    let ask = shown(|s| {
        s.go = true;
        s.gmt = Some(gmt.to_string_lossy().into());
    });
    let go = &levels(&m, &ask).unwrap()["go"];
    assert_eq!(go["go_q"], 0.05);
    assert_eq!(go["min_weight"], 0.5);
    has(go, &["source", "genes", "universe"]);
    // The text's rows, a term each: module, genes, rank, term, namespace,
    // name, k, K, odds, p, q, overlap.
    let text = printed(&m, &ask).unwrap();
    let mut listed = Vec::new();
    let modules = rows(go, "modules");
    assert_eq!(modules.len(), m.gene_roles);
    for module in modules {
        has(module, &["module", "genes", "of", "terms", "why_none"]);
        let terms = rows(module, "terms");
        assert_eq!(
            terms.is_empty(),
            module["why_none"].is_string(),
            "{module:#}"
        );
        for (r, t) in terms.iter().enumerate() {
            has(
                t,
                &[
                    "id",
                    "name",
                    "namespace",
                    "k",
                    "n",
                    "K",
                    "N",
                    "odds",
                    "p",
                    "q",
                    "genes",
                ],
            );
            let genes: Vec<&str> = (t["genes"].as_array().unwrap().iter())
                .map(|g| g.as_str().unwrap())
                .collect();
            listed.push(
                [
                    module["module"].as_str().unwrap().to_string(),
                    format!("{}/{}", module["genes"], module["of"]),
                    (r + 1).to_string(),
                    t["id"].as_str().unwrap().into(),
                    t["namespace"].as_str().unwrap().into(),
                    t["name"].as_str().unwrap().into(),
                    t["k"].to_string(),
                    t["K"].to_string(),
                    format!("{:.2}", t["odds"].as_f64().unwrap()),
                    format!("{:.2e}", t["p"].as_f64().unwrap()),
                    format!("{:.2e}", t["q"].as_f64().unwrap()),
                    genes.join(","),
                ]
                .join("\t"),
            );
        }
    }
    let table: Vec<&String> = text[2..].iter().filter(|l| !l.contains("\t-\t(")).collect();
    assert!(!listed.is_empty());
    assert_eq!(listed.iter().collect::<Vec<_>>(), table);
}
