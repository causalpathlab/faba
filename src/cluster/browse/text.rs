//! The browser's levels as text, for `--batch-process` and the MCP tool:
//! the modules (level 1, see [`modules`]), one gene module's genes
//! and sites (level 2), one gene's sites in every cell module and over
//! every gene module (levels 3 and 4), the cell modules, the fit's trace and
//! verdict, and its settings. Each a list of lines, tab-separated where a
//! line is a row, so a reader (or a program) can take the columns apart.

use super::Modules;

/// What the batch printout shows instead of the modules (level 1): each
/// level the browser opens, named on the command line, each given with
/// --batch-process alone (the browser shows every level itself).
#[derive(clap::Args, Debug, Default, Clone, PartialEq)]
pub struct Shown {
    /// A gene module's genes and sites (`M3` or `3`).
    #[arg(
        long = "module",
        value_name = "M",
        requires = "batch_process",
        help = "With --batch-process: print gene module M's genes and sites (level 2), as M3 or 3"
    )]
    pub module: Option<String>,

    /// A gene's sites in every cell module (by key or symbol).
    #[arg(
        long = "gene",
        value_name = "GENE",
        requires = "batch_process",
        help = "With --batch-process: print a gene's sites in each cell module, over every gene module (levels 3 and 4), by key or symbol"
    )]
    pub gene: Option<String>,

    #[arg(
        long = "cells",
        requires = "batch_process",
        help = "With --batch-process: print the cell modules and their sizes"
    )]
    pub cells: bool,

    #[arg(
        long = "trace",
        requires = "batch_process",
        help = "With --batch-process: print each chain's trace, the verdict and what to do about it"
    )]
    pub trace: bool,

    #[arg(
        long = "settings",
        requires = "batch_process",
        help = "With --batch-process: print the settings the fit ran with"
    )]
    pub settings: bool,
}

impl Shown {
    /// Whether any level beside the modules is asked for.
    pub fn any(&self) -> bool {
        self.module.is_some() || self.gene.is_some() || self.cells || self.trace || self.settings
    }
}

/// What `shown` asks for, each part under a heading, in the order the
/// browser goes: the cell modules, a gene module, a gene, the trace, the
/// settings; the modules (level 1) when it asks for none.
pub fn printed(m: &Modules, shown: &Shown) -> anyhow::Result<Vec<String>> {
    if !shown.any() {
        return Ok(modules(m));
    }
    let mut out = Vec::new();
    let mut part = |lines: Vec<String>| {
        if !out.is_empty() {
            out.push(String::new());
        }
        out.extend(lines);
    };
    if shown.cells {
        part(cells(m));
    }
    if let Some(name) = &shown.module {
        part(module(m, gene_module(m, name)?));
    }
    if let Some(name) = &shown.gene {
        part(gene(m, name)?);
    }
    if shown.trace {
        part(trace(m));
    }
    if shown.settings {
        part(settings(m));
    }
    Ok(out)
}

/// The columns' heads: `C{k} %` and `C{k} reads` for each cell module, in
/// the panes' order (`C{k} UMI/cell` and `C{k} cells` for a count fit).
fn heads(m: &Modules) -> Vec<String> {
    let (unit, total) = if m.per_cell() {
        ("UMI/cell", "cells")
    } else {
        ("%", "reads")
    };
    m.columns
        .iter()
        .flat_map(|c| [format!("C{} {unit}", c + 1), format!("C{} {total}", c + 1)])
        .collect()
}

/// The modules (level 1): per gene module, in the order of their shifts'
/// range, its main genes, that range, and each cell module's value and
/// reads (a count fit: its UMIs per cell and gene, a module's genes
/// averaged, and its cells).
pub fn modules(m: &Modules) -> Vec<String> {
    let mut head = vec!["module".to_string(), "main genes".into(), "range".into()];
    head.extend(heads(m));
    let mut out = vec![head.join("\t")];
    for &g in &m.order {
        let mut row = vec![
            format!("M{}", g + 1),
            m.gene_sizes[g].to_string(),
            format!("{:.2}", m.shift_range(g)),
        ];
        row.extend(values(m.block(g)));
        out.push(row.join("\t"));
    }
    out
}

/// A row's values, `reads` in the panes' order: each cell module's value
/// (`-` without reads) and its reads (cells), to a tenth.
fn values(reads: &[super::Count]) -> Vec<String> {
    reads
        .iter()
        .flat_map(|r| {
            [
                r.number().unwrap_or_else(|| "-".into()),
                format!("{:.1}", r.total),
            ]
        })
        .collect()
}

/// Gene module `name` (`M3`, `m3` or `3`) as an index, when there is one.
fn gene_module(m: &Modules, name: &str) -> anyhow::Result<usize> {
    let n: usize = name
        .trim_start_matches(['M', 'm'])
        .parse()
        .map_err(|_| anyhow::anyhow!("--module takes a gene module, as M3 or 3, not {name}"))?;
    anyhow::ensure!(
        (1..=m.gene_roles).contains(&n),
        "no gene module M{n}: the fit has M1 to M{}",
        m.gene_roles
    );
    Ok(n - 1)
}

/// The cell modules: each its posterior size (its fitted cells'
/// memberships summed) and, for each gene module, its reads there (level
/// 1's column, turned on its side).
fn cells(m: &Modules) -> Vec<String> {
    let mut out = vec![format!(
        "cell modules ({}), in the panes' order: each its posterior size",
        m.cell_roles
    )];
    out.push("module\tcells".into());
    for &c in &m.columns {
        out.push(format!("C{}\t{}", c + 1, m.cells_label(m.cell_sizes[c])));
    }
    out
}

/// Gene module `g` (level 2): its genes, those with the most reads first,
/// and each gene's sites (for a count fit, the gene itself) with their
/// values in each cell module.
fn module(m: &Modules, g: usize) -> Vec<String> {
    let mut out = vec![format!(
        "gene module M{}: {} main genes, {} sites with reads, shifts ranging {:.2} over the cell modules",
        g + 1,
        m.gene_sizes[g],
        m.members[g].len(),
        m.shift_range(g)
    )];
    let mut head = vec!["gene".to_string(), "site".to_string()];
    head.extend(heads(m));
    out.push(head.join("\t"));
    for span in m.genes(g) {
        for i in span.clone() {
            let member = &m.members[g][i];
            let site = &m.sites[member.site];
            let mut row = vec![site.gene.to_string(), place(site)];
            row.extend(values(&member.reads));
            out.push(row.join("\t"));
        }
    }
    out
}

/// A site's place, `chr:pos`; `-` for a count fit's gene, which has none.
fn place(site: &crate::cluster::read::Site) -> String {
    if site.chr.is_empty() {
        "-".into()
    } else {
        format!("{}:{}", site.chr, site.pos)
    }
}

/// Gene `name` (by key or symbol, ASCII case aside; levels 3 and 4): its
/// main gene module, and its sites on the chromosome with most of them,
/// each with its reads in each cell module over every gene module (its
/// full pileup) and, when the tables have them, its counts over every cell
/// read.
fn gene(m: &Modules, name: &str) -> anyhow::Result<Vec<String>> {
    let (key, main) = m
        .find_gene(name)
        .ok_or_else(|| anyhow::anyhow!("no gene {name} in the modules"))?;
    let pileup = m.pileup(key).expect("a gene with sites");
    // A count fit's gene is one feature, with no place of its own.
    let mut head = format!("gene {key}: main gene module M{}", main + 1);
    if !pileup.chr.is_empty() {
        head.push_str(&format!("; {} sites on {}", pileup.sites.len(), pileup.chr));
    }
    if let Some(note) = super::gene::note(m, pileup.chr, &pileup.dropped) {
        head.push_str(&format!("; {note}"));
    }
    let mut out = vec![head];
    let stats = m.stats.as_ref();
    let mut head = vec!["site".to_string()];
    head.extend(heads(m));
    if stats.is_some() {
        head.extend(["converted".into(), "unconverted".into(), "cells".into()]);
    }
    out.push(head.join("\t"));
    for s in pileup.sites {
        let mut row = vec![place(&m.sites[s])];
        row.extend(values(m.whole(s)));
        if let Some(t) = stats {
            row.extend([t.on[s], t.off[s], t.cells[s]].map(|x| x.to_string()));
        }
        out.push(row.join("\t"));
    }
    Ok(out)
}

/// The fit's trace and verdict: each chain's log joint (a count fit's
/// log-likelihood) at a few sweeps, the verdict and what to do about it.
fn trace(m: &Modules) -> Vec<String> {
    let fit = &m.fit;
    let mut out = vec![format!("trace: {} sweeps per chain", fit.span())];
    match &fit.trace {
        Ok(t) => {
            out.push(format!("burn-in {} sweeps", t.burn_in));
            out.push("chain\tsweep\tlog joint".into());
            for (i, chain) in t.chains.iter().enumerate() {
                // About ten points a chain, the last always.
                let step = (chain.len() / 10).max(1);
                let picked = chain
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j % step == 0 || j + 1 == chain.len());
                for (_, (sweep, v)) in picked {
                    out.push(format!("{i}\t{sweep}\t{v:.1}"));
                }
            }
        }
        Err(why) => out.push(format!("no trace: {why}")),
    }
    match fit.verdict() {
        Some(v) => {
            out.push(format!("verdict: {}", v.said()));
            let advice = v.advice();
            if !advice.is_empty() {
                out.push(format!("advice: {advice}"));
            }
        }
        None => out.push("verdict: none (no R-hat recorded)".into()),
    }
    out
}

/// The settings the fit ran with, as its run record holds them (named by
/// their preset, or as a count fit's), or that there is no record.
fn settings(m: &Modules) -> Vec<String> {
    let (model, json) = match &m.run {
        Ok(r) => (format!("the {}", r.label()), serde_json::to_value(r).ok()),
        Err(_) if m.per_cell() => ("the count blockmodel".into(), m.count_settings.clone()),
        Err(why) => (why.clone(), None),
    };
    let mut out = vec![format!("settings: {model}")];
    if let Some(json) = json {
        let json = serde_json::to_string_pretty(&json).unwrap_or_default();
        out.extend(json.lines().map(String::from));
    }
    out
}
