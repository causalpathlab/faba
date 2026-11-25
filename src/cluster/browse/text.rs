//! The browser's levels as text, for `--batch-process`: the fit at a
//! glance (see [`overview`]), the modules (level 1, see [`modules`]), one
//! gene module's genes and sites (level 2), one gene's sites in every cell
//! module and over every gene module (levels 3 and 4), the cell modules,
//! the fit's trace and verdict, and its settings. Each a list of lines,
//! tab-separated where a line is a row, so a reader can take the columns
//! apart; `--json` prints the same levels as JSON (see [`super::json`]).

use super::Modules;

/// What the printout shows instead of the modules (level 1): each level
/// the browser opens, named on the command line. Asking for one prints it
/// as text, as --batch-process does; with --html, which writes the page of
/// every level, it is refused.
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct Shown {
    #[arg(
        long = "overview",
        conflicts_with = "html",
        help = "Print the fit at a glance: its verdict, its gene and cell modules' sizes and each gene module's top genes by weight (GO only with --go)"
    )]
    pub overview: bool,

    #[arg(
        long = "grid",
        conflicts_with = "html",
        help = "Print the modules (level 1), as no level asked for does, beside the levels named"
    )]
    pub grid: bool,

    #[arg(
        long = "json",
        conflicts_with = "html",
        help = "Print the levels as one JSON object, a key a level, its rows objects with named fields, instead of tab-separated text"
    )]
    pub json: bool,

    /// A gene module's genes and sites (`M3` or `3`).
    #[arg(
        long = "module",
        value_name = "M",
        conflicts_with = "html",
        help = "Print gene module M's genes and sites (level 2), as M3 or 3"
    )]
    pub module: Option<String>,

    /// A gene's sites in every cell module (by key or symbol).
    #[arg(
        long = "gene",
        value_name = "GENE",
        conflicts_with = "html",
        help = "Print a gene's sites in each cell module, over every gene module (levels 3 and 4), by key or symbol"
    )]
    pub gene: Option<String>,

    #[arg(
        long = "cells",
        conflicts_with = "html",
        help = "Print the cell modules and their sizes"
    )]
    pub cells: bool,

    #[arg(
        long = "trace",
        conflicts_with = "html",
        help = "Print each restart's ELBO trace, the verdict and what to do about it"
    )]
    pub trace: bool,

    #[arg(
        long = "go",
        conflicts_with = "html",
        help = "Print the GO terms (or --gmt gene sets) over-represented in each gene module"
    )]
    pub go: bool,

    #[arg(
        long = "gaf",
        value_name = "FILE",
        value_hint = clap::ValueHint::FilePath,
        help = "GO annotations for --go and the browser's a (default: lupin's cached GO Consortium file for the genes' species)"
    )]
    pub gaf: Option<String>,

    #[arg(
        long = "obo",
        value_name = "FILE",
        value_hint = clap::ValueHint::FilePath,
        help = "The GO graph for --gaf, terms propagated up it (default: lupin's cached go-basic.obo)"
    )]
    pub obo: Option<String>,

    #[arg(
        long = "gmt",
        value_name = "FILE",
        value_hint = clap::ValueHint::FilePath,
        conflicts_with = "gaf",
        help = "Gene sets to test instead of GO, a GMT file"
    )]
    pub gmt: Option<String>,

    #[arg(
        long = "go-min-weight",
        value_name = "WEIGHT",
        default_value_t = 0.5,
        help = "A gene counts in its main gene module for --go when its weight there is at least this"
    )]
    pub go_min_weight: f64,

    #[arg(
        long = "batches",
        conflicts_with = "html",
        help = "Print how each cell module spreads over the batches: shared, one condition's, or batch-specific"
    )]
    pub batches: bool,

    #[arg(
        long = "settings",
        conflicts_with = "html",
        help = "Print the settings the fit ran with"
    )]
    pub settings: bool,
}

impl Default for Shown {
    /// As the command line has it: no level, GO's weight cut at 0.5.
    fn default() -> Self {
        use clap::Parser;
        #[derive(Parser)]
        struct Line {
            #[command(flatten)]
            shown: Shown,
            /// What the levels conflict with in `faba cluster-view`.
            #[arg(id = "html", long = "html", hide = true)]
            _html: Option<String>,
        }
        Line::parse_from(["shown"]).shown
    }
}

impl Shown {
    /// Whether any level beside the modules is asked for.
    pub fn any(&self) -> bool {
        self.overview
            || self.grid
            || self.module.is_some()
            || self.gene.is_some()
            || self.cells
            || self.batches
            || self.go
            || self.trace
            || self.settings
    }
}

/// What `shown` asks for, each part under a heading, in the order the
/// browser goes: the overview, the modules, the cell modules, a gene
/// module, a gene, the batches, the GO terms, the trace, the settings; the
/// modules (level 1) when it asks for none.
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
    if shown.overview {
        part(overview(m));
    }
    if shown.grid {
        part(modules(m));
    }
    if shown.cells {
        part(cells(m));
    }
    if let Some(name) = &shown.module {
        part(module(m, gene_module(m, name)?));
    }
    if let Some(name) = &shown.gene {
        part(gene(m, name)?);
    }
    if shown.batches {
        part(batches(m));
    }
    if shown.go {
        part(go_lines(m, shown)?);
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
        if m.per_cell() {
            // A count module's UMIs per gene and cell, beside each cell
            // module's size (the block's own total sums its genes' cells).
            let sizes = m.columns.iter().map(|&c| m.cells_label(m.cell_sizes[c]));
            row.extend(
                m.block(g)
                    .iter()
                    .zip(sizes)
                    .flat_map(|(r, n)| [r.number().unwrap_or_else(|| "-".into()), n]),
            );
        } else {
            row.extend(values(m.block(g)));
        }
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

/// The fit at a glance: its verdict, its gene and cell modules, each cell
/// module's size, and each gene module's genes with its top ones by
/// weight. GO is not read.
fn overview(m: &Modules) -> Vec<String> {
    let verdict = match m.fit.verdict() {
        Some(v) => format!("{}: {}", v.badge(), v.said()),
        None => format!("no verdict ({})", m.fit.unjudged()),
    };
    let mut out = vec![
        format!(
            "overview: {} fit, {verdict}; {} gene modules, {} cell modules",
            m.modality, m.gene_roles, m.cell_roles
        ),
        "module\tcells".into(),
    ];
    for &c in &m.columns {
        out.push(format!("C{}\t{}", c + 1, m.cells_label(m.cell_sizes[c])));
    }
    out.push("module\tgenes\ttop genes by weight".into());
    for g in 0..m.gene_roles {
        let top: Vec<String> = (top_genes(m, g).iter())
            .map(|(gene, w)| format!("{gene} ({w:.2})"))
            .collect();
        out.push(format!(
            "M{}\t{}\t{}",
            g + 1,
            m.gene_sizes[g],
            top.join(", ")
        ));
    }
    out
}

/// Genes the overview names a gene module by, at most.
const TOP_GENES: usize = 5;

/// Gene module `g`'s main genes with the most weight there, the most first
/// (by key on a tie), at most [`TOP_GENES`]: each its symbol and weight.
pub(super) fn top_genes(m: &Modules, g: usize) -> Vec<(Box<str>, f32)> {
    let mut genes: Vec<&(Box<str>, usize, f32)> = m
        .gene_main
        .iter()
        .filter(|(_, role, _)| *role == g)
        .collect();
    genes.sort_by(|a, b| b.2.total_cmp(&a.2).then_with(|| a.0.cmp(&b.0)));
    (genes.into_iter().take(TOP_GENES))
        .map(|(key, _, w)| (crate::site_analysis::pileup::query_symbol(key), *w))
        .collect()
}

/// Gene module `name` (`M3`, `m3` or `3`) as an index, when there is one.
pub(super) fn gene_module(m: &Modules, name: &str) -> anyhow::Result<usize> {
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
pub(super) fn place(site: &crate::cluster::read::Site) -> String {
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

/// The GO terms (or gene sets) of each gene module (see [`super::go`]).
fn go_lines(m: &Modules, shown: &Shown) -> anyhow::Result<Vec<String>> {
    let ask = shown.go_ask();
    let source = ask.source(m).map_err(anyhow::Error::msg)?;
    Ok(super::go::table(m, &source, ask.min_weight)?.lines())
}

impl Shown {
    /// What `--go` and the browser's `a` ask of GO.
    pub fn go_ask(&self) -> super::go::Ask {
        super::go::Ask {
            gaf: self.gaf.clone(),
            gmt: self.gmt.clone(),
            obo: self.obo.clone(),
            min_weight: self.go_min_weight,
        }
    }
}

/// How each cell module spreads over the batches (see
/// [`crate::cluster::batches`]): a line a module, in the panes' order, its
/// label, evenness and share in each batch; then each condition's batches.
fn batches(m: &Modules) -> Vec<String> {
    let Some(b) = &m.batches else {
        return vec!["batches: the cells are of one batch".into()];
    };
    let d = &b.design;
    let mut out = vec![format!(
        "batches: each cell module's share in each batch (its cells by batch size); \
         under {:.2} evenness over a condition's replicates it is batch-specific",
        crate::cluster::batches::EVEN
    )];
    let order = d.order();
    out.push(
        ["module", "label", "evenness"]
            .into_iter()
            .map(String::from)
            .chain(order.iter().map(|&i| d.batches[i].clone()))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    for &c in &m.columns {
        let s = &b.spread[c];
        let even = s.evenness.map_or("-".into(), |e| format!("{e:.2}"));
        let shares = order.iter().map(|&i| format!("{:.3}", s.per_batch[i]));
        let row: Vec<String> = [format!("C{}", c + 1), s.label.word(), even]
            .into_iter()
            .chain(shares)
            .collect();
        out.push(row.join("\t"));
    }
    out.extend(d.listing().into_iter().map(|l| format!("condition {l}")));
    out
}

/// The fit's trace and verdict: each restart's ELBO at a few steps, the
/// verdict and what to do about it.
fn trace(m: &Modules) -> Vec<String> {
    let fit = &m.fit;
    let mut out = vec![format!("trace: {} steps per restart", fit.span())];
    match &fit.trace {
        Ok(t) => {
            out.push(format!(
                "k-means start {} steps, then epochs of variational inference",
                t.warmup
            ));
            out.push("restart\tstep\tELBO".into());
            push_points(&mut out, &t.restarts);
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
        None => out.push(format!("verdict: none ({})", fit.unjudged())),
    }
    out
}

/// About ten points of each of `restarts`, the last always, a line each.
fn push_points(out: &mut Vec<String>, restarts: &[Vec<(usize, f64)>]) {
    for (i, trace) in restarts.iter().enumerate() {
        for (step, v) in sampled(trace) {
            out.push(format!("{i}\t{step}\t{v:.1}"));
        }
    }
}

/// About ten points of `trace`, the last always.
pub(super) fn sampled(trace: &[(usize, f64)]) -> impl Iterator<Item = (usize, f64)> + '_ {
    let every = (trace.len() / 10).max(1);
    (trace.iter().enumerate())
        .filter(move |&(j, _)| j % every == 0 || j + 1 == trace.len())
        .map(|(_, &p)| p)
}

/// The settings the fit ran with, as its run record holds them (named by
/// their preset, or as a count fit's), or that there is no record.
fn settings(m: &Modules) -> Vec<String> {
    let (model, json) = recorded(m);
    let mut out = vec![format!("settings: {model}")];
    if let Some(json) = json {
        let json = serde_json::to_string_pretty(&json).unwrap_or_default();
        out.extend(json.lines().map(String::from));
    }
    out
}

/// The model the fit ran with, named (by its preset, or as a count fit's),
/// and its settings as its run record holds them; why not, and none,
/// without a record.
pub(super) fn recorded(m: &Modules) -> (String, Option<serde_json::Value>) {
    match &m.run {
        Ok(r) => (format!("the {}", r.label()), serde_json::to_value(r).ok()),
        // A count fit's settings are not the blockmodel's: kept as recorded.
        Err(_) if m.count_settings.is_some() => {
            ("the count blockmodel".into(), m.count_settings.clone())
        }
        Err(why) => (why.to_string(), None),
    }
}
