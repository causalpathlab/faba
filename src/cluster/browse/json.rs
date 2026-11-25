//! The browser's levels as JSON, for `--json` and the MCP tool: the same
//! levels [`super::text`] prints, each an object under its name, its rows
//! objects with named fields and its numbers numbers, so a program reads
//! them without taking text apart. GO is read only for `go`.

use super::text::{gene_module, place, sampled, top_genes, Shown};
use super::{Count, Modules};
use serde_json::{json, Map, Value};

/// What `shown` asks for, a key a level (`overview`, `modules`, `cells`,
/// `module`, `gene`, `batches`, `go`, `trace`, `settings`); the modules
/// (level 1) when it asks for none.
pub fn levels(m: &Modules, shown: &Shown) -> anyhow::Result<Value> {
    let mut out = Map::new();
    let mut put = |key: &str, v: Value| {
        out.insert(key.into(), v);
    };
    if !shown.any() || shown.grid {
        put("modules", modules(m));
    }
    if shown.overview {
        put("overview", overview(m));
    }
    if shown.cells {
        put("cells", cells(m));
    }
    if let Some(name) = &shown.module {
        put("module", module(m, gene_module(m, name)?));
    }
    if let Some(name) = &shown.gene {
        put("gene", gene(m, name)?);
    }
    if shown.batches {
        put("batches", batches(m));
    }
    if shown.go {
        let ask = shown.go_ask();
        let source = ask.source(m).map_err(anyhow::Error::msg)?;
        put(
            "go",
            super::go::table(m, &source, ask.min_weight, ask.q)?.json(),
        );
    }
    if shown.trace {
        put("trace", trace(m));
    }
    if shown.settings {
        put("settings", settings(m));
    }
    Ok(Value::Object(out))
}

/// `x` to four decimals, as the rows carry it.
pub(super) fn round(x: f64) -> f64 {
    (x * 1e4).round() / 1e4
}

/// What a row's values are and what is beside them: `%` and `reads` for a
/// site fit, `UMI/cell` and `cells` for a count fit.
fn units(m: &Modules) -> (&'static str, &'static str) {
    match m.per_cell() {
        true => ("UMI/cell", "cells"),
        false => ("%", "reads"),
    }
}

/// A value of `r`: its methylated % (a count fit's UMIs per cell); none
/// without reads.
fn value(m: &Modules, r: &Count) -> Option<f64> {
    let scale = if m.per_cell() { 1.0 } else { 100.0 };
    r.fraction().map(|f| round(scale * f))
}

/// `reads`, in the panes' order, a cell module each: its name, its value
/// and what is beside it (`totals`, by column; else each count's own).
fn values(m: &Modules, reads: &[Count], totals: Option<&[f64]>) -> Value {
    let (_, total) = units(m);
    let row = (m.columns.iter().zip(reads).enumerate()).map(|(k, (c, r))| {
        let beside = totals.map_or(r.total, |t| t[k]);
        let mut v = Map::new();
        v.insert("cell_module".into(), json!(format!("C{}", c + 1)));
        v.insert("value".into(), json!(value(m, r)));
        v.insert(total.into(), json!(round(beside)));
        Value::Object(v)
    });
    Value::Array(row.collect())
}

/// The modules (level 1): a row a gene module with sites, in the order of
/// their shifts' range, its main genes, that range, and its value in each
/// cell module (a count fit's beside each cell module's size).
pub fn modules(m: &Modules) -> Value {
    let (unit, total) = units(m);
    let sizes: Vec<f64> = m.columns.iter().map(|&c| m.cell_sizes[c]).collect();
    let totals = m.per_cell().then_some(sizes.as_slice());
    let rows: Vec<Value> = (m.order.iter())
        .map(|&g| {
            json!({
                "module": format!("M{}", g + 1),
                "main_genes": m.gene_sizes[g],
                "range": round(m.shift_range(g)),
                "values": values(m, m.block(g), totals),
            })
        })
        .collect();
    json!({
        "unit": unit,
        "total": total,
        "cell_modules": cell_names(m),
        "rows": rows,
    })
}

/// The cell modules' names, in the panes' order.
fn cell_names(m: &Modules) -> Vec<String> {
    m.columns.iter().map(|c| format!("C{}", c + 1)).collect()
}

/// Each cell module's posterior size, in the panes' order.
fn cell_sizes(m: &Modules) -> Vec<Value> {
    (m.columns.iter())
        .map(|&c| json!({ "module": format!("C{}", c + 1), "cells": m.cell_sizes[c] }))
        .collect()
}

/// The fit at a glance: its verdict, its gene and cell modules, each cell
/// module's size, and each gene module's genes with its top ones by
/// weight. GO is not read.
fn overview(m: &Modules) -> Value {
    let v = m.fit.verdict();
    let genes: Vec<Value> = (0..m.gene_roles)
        .map(|g| {
            let top: Vec<Value> = (top_genes(m, g).into_iter())
                .map(|(gene, w)| json!({ "gene": gene, "weight": round(f64::from(w)) }))
                .collect();
            json!({ "module": format!("M{}", g + 1), "genes": m.gene_sizes[g], "top": top })
        })
        .collect();
    json!({
        "modality": m.modality,
        "verdict": v.map(|v| v.badge()),
        "said": v.map(|v| v.said()),
        "unjudged": v.is_none().then(|| m.fit.unjudged()),
        "gene_modules": m.gene_roles,
        "cell_modules": m.cell_roles,
        "approximate_sizes": m.approximate.sizes,
        "cell_sizes": cell_sizes(m),
        "gene_sizes": genes,
    })
}

/// The cell modules: each its posterior size.
fn cells(m: &Modules) -> Value {
    json!({ "approximate": m.approximate.sizes, "rows": cell_sizes(m) })
}

/// Gene module `g` (level 2): its genes, those with the most reads first,
/// and each gene's sites with their values in each cell module.
fn module(m: &Modules, g: usize) -> Value {
    let (unit, total) = units(m);
    let rows: Vec<Value> = (m.genes(g).iter().flat_map(Clone::clone))
        .map(|i| {
            let member = &m.members[g][i];
            let site = &m.sites[member.site];
            json!({
                "gene": site.gene,
                "site": at(site),
                "values": values(m, &member.reads, None),
            })
        })
        .collect();
    json!({
        "module": format!("M{}", g + 1),
        "main_genes": m.gene_sizes[g],
        "sites": m.members[g].len(),
        "range": round(m.shift_range(g)),
        "unit": unit,
        "total": total,
        "rows": rows,
    })
}

/// A site's place, `chr:pos`; none for a count fit's gene.
fn at(site: &crate::cluster::read::Site) -> Option<String> {
    (!site.chr.is_empty()).then(|| place(site))
}

/// Gene `name` (levels 3 and 4): its main gene module, and its sites on
/// the chromosome with most of them, each with its reads in each cell
/// module over every gene module and, when the tables have them, its
/// counts over every cell read.
fn gene(m: &Modules, name: &str) -> anyhow::Result<Value> {
    let (key, main) = m
        .find_gene(name)
        .ok_or_else(|| anyhow::anyhow!("no gene {name} in the modules"))?;
    let pileup = m.pileup(key).expect("a gene with sites");
    let (unit, total) = units(m);
    let rows: Vec<Value> = (pileup.sites.clone())
        .map(|s| {
            let mut row = json!({ "site": at(&m.sites[s]), "values": values(m, m.whole(s), None) });
            if let Some(t) = &m.stats {
                row["converted"] = json!(t.on[s]);
                row["unconverted"] = json!(t.off[s]);
                row["cells"] = json!(t.cells[s]);
            }
            row
        })
        .collect();
    Ok(json!({
        "gene": key,
        "main_module": format!("M{}", main + 1),
        "chr": (!pileup.chr.is_empty()).then_some(pileup.chr),
        "note": super::gene::note(m, pileup.chr, &pileup.dropped),
        "unit": unit,
        "total": total,
        "rows": rows,
    }))
}

/// How each cell module spreads over the batches: a row a module, in the
/// panes' order, its label, evenness and share in each batch by name; the
/// conditions; none, said, for cells of one batch.
fn batches(m: &Modules) -> Value {
    let Some(b) = &m.batches else {
        return json!({ "batches": null, "note": "the cells are of one batch" });
    };
    let d = &b.design;
    let order = d.order();
    let rows: Vec<Value> = (m.columns.iter())
        .map(|&c| {
            let s = &b.spread[c];
            let shares: Map<String, Value> = (order.iter())
                .map(|&i| (d.batches[i].clone(), json!(round(s.per_batch[i]))))
                .collect();
            json!({
                "module": format!("C{}", c + 1),
                "label": s.label.word(),
                "batch_specific": s.label.is_batch(),
                "evenness": s.evenness.map(round),
                "shares": shares,
            })
        })
        .collect();
    let batches: Vec<&String> = order.iter().map(|&i| &d.batches[i]).collect();
    json!({
        "even": crate::cluster::batches::EVEN,
        "batches": batches,
        "rows": rows,
        "conditions": d.listing(),
    })
}

/// The fit's trace and verdict: about ten points of each restart's ELBO,
/// the verdict and what to do about it, or why there is none.
fn trace(m: &Modules) -> Value {
    let fit = &m.fit;
    let (restarts, warmup, why) = match &fit.trace {
        Ok(t) => {
            let restarts: Vec<Value> = (t.restarts.iter())
                .map(|r| {
                    let points = sampled(r).map(|(step, v)| json!({ "step": step, "elbo": v }));
                    Value::Array(points.collect())
                })
                .collect();
            (restarts, Some(t.warmup), None)
        }
        Err(why) => (Vec::new(), None, Some(why.clone())),
    };
    let v = fit.verdict();
    json!({
        "steps": fit.span(),
        "warmup": warmup,
        "restarts": restarts,
        "no_trace": why,
        "verdict": v.map(|v| v.badge()),
        "converged": v.map(|v| v.converged()),
        "said": v.map(|v| v.said()),
        "advice": v.map(|v| v.advice()).filter(|a| !a.is_empty()),
        "unjudged": v.is_none().then(|| fit.unjudged()),
    })
}

/// The settings the fit ran with, as its run record holds them, or why
/// there are none.
fn settings(m: &Modules) -> Value {
    let (model, settings) = super::text::recorded(m);
    json!({ "model": model, "settings": settings })
}
