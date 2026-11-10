//! GO terms (or any gene sets) over-represented in a fit's gene modules.
//!
//! Each gene module's genes (those whose main gene role it is, at a weight
//! of at least `--go-min-weight`) are tested against the fit's own genes
//! that some term annotates (the universe), term by term: the
//! hypergeometric upper tail of the overlap, with Haldane's odds ratio, and
//! BH over every (module, term) pair tested: every term of 5 to 500 of the
//! universe's genes against every module of at least 5, whatever their
//! overlap (lest the family depend on the outcome). GO annotations (a GAF) are propagated up `is_a` and `part_of` (the
//! GO graph, an OBO), their NAS and ND rows left out (an author's statement
//! nobody traced, no data); a GMT's sets are taken as they are. A module's
//! terms are listed by p, a term left out when its overlap is a term's
//! already listed (a parent with the same genes).
//!
//! The fit's genes are matched to the annotations' by symbol or Ensembl id
//! first; a gene not matched so then by one of its aliases, only onto a fit
//! gene no other matched and an alias no other gene has, so one gene never
//! takes another's terms.

use super::Modules;
use data_beans::aux::gene_sets::{read_gaf, read_gmt, GafOpts, GeneSets};
use data_beans::aux::ontology::Ontology;
use faba::hypothesis_tests::ln_fisher_exact_greater;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// The sizes, in the universe's genes, of the terms tested.
const SIZES: std::ops::RangeInclusive<usize> = 5..=500;
/// A module tested has at least this many genes in the universe.
const LEAST_MODULE: usize = 5;
/// Terms listed per module at most, and the q they need.
const TOP: usize = 10;
pub(super) const Q: f64 = 0.05;
/// GAF evidence left out (see the module).
const DROPPED_EVIDENCE: [&str; 2] = ["NAS", "ND"];

/// Where the gene sets come from: GO annotations with the GO graph, or a
/// GMT.
#[derive(Debug, Clone, PartialEq)]
pub enum Source {
    Go { gaf: PathBuf, obo: PathBuf },
    Gmt(PathBuf),
}

impl Source {
    /// The source named: a GAF (with its OBO), else a GMT; with neither,
    /// the GO Consortium's annotations for the genes' species from lupin's
    /// cache when they are there. Why not, else.
    pub fn of(
        gaf: Option<&str>,
        gmt: Option<&str>,
        obo: Option<&str>,
        genes: &[(Box<str>, usize, f32)],
    ) -> Result<Self, String> {
        let cache = lupin_cache();
        let in_cache = |name: &str| cache.as_ref().map(|d| d.join(name)).filter(|p| p.is_file());
        let obo = obo.map(PathBuf::from).or_else(|| in_cache("go-basic.obo"));
        match (gaf, gmt) {
            (Some(_), Some(_)) => Err("give --gaf or --gmt, not both".into()),
            (None, Some(gmt)) => Ok(Source::Gmt(gmt.into())),
            (Some(gaf), None) => Ok(Source::Go {
                gaf: gaf.into(),
                obo: obo.ok_or("--gaf needs --obo, the GO graph (go-basic.obo)")?,
            }),
            (None, None) => {
                let file = species_gaf(genes).ok_or(
                    "the genes are not human or mouse (Ensembl ids): give --gaf and --obo, or --gmt",
                )?;
                match (in_cache(file), obo) {
                    (Some(gaf), Some(obo)) => Ok(Source::Go { gaf, obo }),
                    _ => Err(format!(
                        "no {file} and go-basic.obo in lupin's cache: fetch them (`lupin data \
                         --fetch`) or give --gaf and --obo"
                    )),
                }
            }
        }
    }

    /// In a few words, for the output's head.
    pub fn said(&self) -> String {
        match self {
            Source::Go { gaf, .. } => format!("GO annotations {}", name(gaf)),
            Source::Gmt(gmt) => format!("gene sets {}", name(gmt)),
        }
    }
}

fn name(p: &Path) -> String {
    p.file_name()
        .map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into())
}

/// lupin's cache of data files, where lupin keeps them (`dirs::cache_dir`):
/// `~/Library/Caches/lupin` on macOS; else `$XDG_CACHE_HOME/lupin` when that
/// is an absolute path, `~/.cache/lupin` otherwise.
fn lupin_cache() -> Option<PathBuf> {
    let home = PathBuf::from(std::env::var_os("HOME")?);
    if cfg!(target_os = "macos") {
        return Some(home.join("Library/Caches/lupin"));
    }
    let xdg = std::env::var_os("XDG_CACHE_HOME").map(PathBuf::from);
    Some(match xdg.filter(|x| x.is_absolute()) {
        Some(x) => x.join("lupin"),
        None => home.join(".cache/lupin"),
    })
}

/// The GO Consortium's annotation file for the genes' species, by their
/// Ensembl ids (`ENSG` human, `ENSMUSG` mouse); none for another.
fn species_gaf(genes: &[(Box<str>, usize, f32)]) -> Option<&'static str> {
    let count = |p: &str| genes.iter().filter(|g| g.0.starts_with(p)).count();
    let (human, mouse) = (count("ENSG0"), count("ENSMUSG"));
    match human.max(mouse) * 2 > genes.len() {
        true if human >= mouse => Some("goa_human.gaf.gz"),
        true => Some("mgi.gaf.gz"),
        false => None,
    }
}

/// What a call asks GO of: the files named (none: found, see
/// [`Source::of`]) and the weight a gene needs in its main module.
#[derive(Debug, Clone, PartialEq)]
pub struct Ask {
    pub gaf: Option<String>,
    pub gmt: Option<String>,
    pub obo: Option<String>,
    pub min_weight: f64,
}

impl Ask {
    /// The source for `m`'s genes.
    pub fn source(&self, m: &Modules) -> Result<Source, String> {
        Source::of(
            self.gaf.as_deref(),
            self.gmt.as_deref(),
            self.obo.as_deref(),
            &m.gene_main,
        )
    }
}

/// One term over-represented in a module.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    pub id: String,
    pub name: String,
    /// BP, MF or CC (a GO term's namespace); empty for a GMT's set.
    pub namespace: &'static str,
    /// The module's genes in the term, of the module's (`n`), of the
    /// term's (`big_k`) and the universe's (`big_n`).
    pub k: usize,
    pub n: usize,
    pub big_k: usize,
    pub big_n: usize,
    pub odds: f64,
    /// `ln p` and the BH q.
    pub ln_p: f64,
    pub q: f64,
    /// The overlap's genes' symbols, the module's strongest first.
    pub genes: Vec<String>,
}

/// One gene module's terms: its genes before and after the weight cut (in
/// the universe), and the terms listed.
#[derive(Debug, Clone, PartialEq)]
pub struct ModuleTerms {
    pub all: usize,
    pub n: usize,
    pub hits: Vec<Hit>,
}

/// Every gene module's terms, and what they were tested against.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    pub source: String,
    /// The fit's genes, and those some term annotates (the universe).
    pub genes: usize,
    pub universe: usize,
    pub min_weight: f64,
    /// Every gene module, by its number from 0.
    pub modules: Vec<ModuleTerms>,
}

/// A gene's symbol: its key less its Ensembl id (`ENSG..._DPM1`: `DPM1`;
/// `ENSG..._Y_RNA`: `Y_RNA`), as the browser reads it.
fn symbol(key: &str) -> Box<str> {
    crate::site_analysis::pileup::query_symbol(key)
}

/// A gene's Ensembl id, its version left out; none for a key without one.
fn ensembl(key: &str) -> Option<&str> {
    let id = key.split('_').next().filter(|id| id.starts_with("ENS"))?;
    Some(id.split('.').next().unwrap_or(id))
}

/// The gene sets of `source`, read.
fn read(source: &Source) -> anyhow::Result<(GeneSets, FxHashMap<Box<str>, &'static str>)> {
    match source {
        Source::Gmt(gmt) => Ok((read_gmt(&gmt.to_string_lossy())?, FxHashMap::default())),
        Source::Go { gaf, obo } => {
            let onto = Ontology::load_obo(&obo.to_string_lossy())?;
            let kept = without_dropped_evidence(gaf)?;
            let raw = read_gaf(&kept.path().to_string_lossy(), &GafOpts { no_iea: false })?;
            Ok((raw.into_gene_sets(Some(&onto)), namespaces(obo)?))
        }
    }
}

/// `gaf` without its [`DROPPED_EVIDENCE`] rows, as a file of its own that
/// lasts as long as the handle.
fn without_dropped_evidence(gaf: &Path) -> anyhow::Result<tempfile::NamedTempFile> {
    let reader = legume_numeric::matrix::common_io::open_buf_reader(&gaf.to_string_lossy())?;
    let mut out = tempfile::Builder::new().suffix(".gaf").tempfile()?;
    {
        let mut w = std::io::BufWriter::new(out.as_file_mut());
        for line in reader.lines() {
            let line = line?;
            let evidence = line.split('\t').nth(6).unwrap_or("");
            if !line.starts_with('!') && DROPPED_EVIDENCE.contains(&evidence) {
                continue;
            }
            writeln!(w, "{line}")?;
        }
        w.flush()?;
    }
    Ok(out)
}

/// Each GO term's namespace, BP, MF or CC, from the OBO.
fn namespaces(obo: &Path) -> anyhow::Result<FxHashMap<Box<str>, &'static str>> {
    let reader = legume_numeric::matrix::common_io::open_buf_reader(&obo.to_string_lossy())?;
    let mut out = FxHashMap::default();
    let mut id: Option<Box<str>> = None;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('[') {
            id = None;
        } else if let Some(x) = line.strip_prefix("id: ") {
            id = Some(x.trim().into());
        } else if let Some(ns) = line.strip_prefix("namespace: ") {
            let short = match ns.trim() {
                "biological_process" => "BP",
                "molecular_function" => "MF",
                "cellular_component" => "CC",
                _ => continue,
            };
            if let Some(id) = id.take() {
                out.insert(id, short);
            }
        }
    }
    Ok(out)
}

/// Each annotated gene (by its key in `sets`) matched to a row of `genes`:
/// by symbol or Ensembl id, then, for those left, by an alias unique to it
/// onto a row no other took (see the module).
fn matched(sets: &GeneSets, genes: &[(Box<str>, usize, f32)]) -> FxHashMap<Box<str>, usize> {
    let mut by_name: FxHashMap<String, usize> = FxHashMap::default();
    for (row, (key, ..)) in genes.iter().enumerate() {
        by_name
            .entry(symbol(key).to_ascii_lowercase())
            .or_insert(row);
        if let Some(id) = ensembl(key) {
            by_name.entry(id.to_ascii_lowercase()).or_insert(row);
        }
    }
    let mut out = FxHashMap::default();
    let mut taken: FxHashSet<usize> = FxHashSet::default();
    // Every gene of every set is a key of `gene_aliases` (GAF and GMT alike).
    let mut keys: Vec<&Box<str>> = sets.gene_aliases.keys().collect();
    keys.sort();
    for gene in keys {
        if let Some(&row) = by_name.get(&gene.to_ascii_lowercase()) {
            if taken.insert(row) {
                out.insert(gene.clone(), row);
            }
        }
    }
    // An alias names one gene alone to count.
    let mut owners: FxHashMap<String, usize> = FxHashMap::default();
    for aliases in sets.gene_aliases.values() {
        for a in aliases {
            *owners.entry(a.to_ascii_lowercase()).or_default() += 1;
        }
    }
    let mut left: Vec<&Box<str>> = sets
        .gene_aliases
        .keys()
        .filter(|g| !out.contains_key(*g))
        .collect();
    left.sort();
    for gene in left {
        let hit = sets.gene_aliases[gene]
            .iter()
            .map(|a| a.to_ascii_lowercase())
            .find_map(|a| {
                let row = *by_name.get(&a)?;
                (owners.get(&a) == Some(&1) && !taken.contains(&row)).then_some(row)
            });
        if let Some(row) = hit {
            taken.insert(row);
            out.insert(gene.clone(), row);
        }
    }
    out
}

/// The terms of `m`'s gene modules from `source`, each module's genes those
/// whose main role it is at a weight of at least `min_weight`.
pub fn table(m: &Modules, source: &Source, min_weight: f64) -> anyhow::Result<Table> {
    terms(&m.gene_main, m.gene_roles, source, min_weight)
}

/// [`table`] of the genes `genes` (key, main role, its weight) over
/// `roles` gene roles.
pub fn terms(
    genes: &[(Box<str>, usize, f32)],
    roles: usize,
    source: &Source,
    min_weight: f64,
) -> anyhow::Result<Table> {
    let (sets, ns) = read(source)?;
    let rows = matched(&sets, genes);
    // Each term's rows; the universe: every row some term has.
    let mut term_rows: Vec<(&Box<str>, Vec<usize>)> = (sets.term_genes.iter())
        .map(|(t, g)| {
            let mut r: Vec<usize> = g.iter().filter_map(|g| rows.get(g).copied()).collect();
            r.sort_unstable();
            r.dedup();
            (t, r)
        })
        .filter(|(_, r)| !r.is_empty())
        .collect();
    term_rows.sort_by(|a, b| a.0.cmp(b.0));
    let mut universe = vec![false; genes.len()];
    term_rows
        .iter()
        .flat_map(|(_, r)| r)
        .for_each(|&r| universe[r] = true);
    let big_n = universe.iter().filter(|&&u| u).count();
    anyhow::ensure!(
        big_n >= 10 && big_n * 10 >= genes.len(),
        "{}: only {big_n} of the fit's {} genes annotated: are they the same species, \
         and their names symbols or Ensembl ids?",
        source.said(),
        genes.len()
    );
    let tested: Vec<&(&Box<str>, Vec<usize>)> = term_rows
        .iter()
        .filter(|(_, r)| SIZES.contains(&r.len()))
        .collect();
    // Every module's genes in the universe, before and after the cut.
    let mut modules: Vec<(usize, Vec<usize>)> = vec![(0, Vec::new()); roles];
    for (row, (_, role, w)) in genes.iter().enumerate() {
        if universe[row] && *role < roles {
            modules[*role].0 += 1;
            if f64::from(*w) >= min_weight {
                modules[*role].1.push(row);
            }
        }
    }
    // Every (module, term) tested: its overlap and ln p.
    let mut all: Vec<(usize, usize, usize, f64)> = Vec::new();
    for (g, (_, members)) in modules.iter().enumerate() {
        let n = members.len();
        if n < LEAST_MODULE {
            continue;
        }
        let mut is = vec![false; genes.len()];
        members.iter().for_each(|&r| is[r] = true);
        for (t, (_, r)) in tested.iter().enumerate() {
            let k = r.iter().filter(|&&x| is[x]).count();
            let [a, b, c, d] = cells(k, n, r.len(), big_n);
            let ln_p = match k {
                0 => 0.0,
                _ => ln_fisher_exact_greater(a, b, c, d),
            };
            all.push((g, t, k, ln_p));
        }
    }
    let q = bh(&all.iter().map(|x| x.3).collect::<Vec<_>>());
    let mut out: Vec<ModuleTerms> = (modules.iter())
        .map(|(all, members)| ModuleTerms {
            all: *all,
            n: members.len(),
            hits: Vec::new(),
        })
        .collect();
    let mut order: Vec<usize> = (0..all.len())
        .filter(|&i| q[i] <= Q && all[i].2 >= 2)
        .collect();
    order.sort_by(|&a, &b| all[a].3.total_cmp(&all[b].3).then(all[a].1.cmp(&all[b].1)));
    // Each module's terms by p, a term whose overlap one listed has left out.
    let mut seen: Vec<Vec<Vec<usize>>> = vec![Vec::new(); roles];
    for i in order {
        let (g, t, k, ln_p) = all[i];
        let (id, r) = tested[t];
        let module = &mut out[g];
        if module.hits.len() >= TOP {
            continue;
        }
        // The overlap: both lists ascend.
        let members = &modules[g].1;
        let overlap: Vec<usize> = (r.iter().copied())
            .filter(|x| members.binary_search(x).is_ok())
            .collect();
        if seen[g].contains(&overlap) {
            continue;
        }
        let (n, big_k) = (module.n, r.len());
        let [a, b, c, d] = cells(k, n, big_k, big_n).map(|x| x as f64);
        let mut strongest = overlap.clone();
        strongest.sort_by(|&x, &y| genes[y].2.total_cmp(&genes[x].2).then(x.cmp(&y)));
        module.hits.push(Hit {
            id: id.to_string(),
            name: sets
                .names
                .get(*id)
                .map_or_else(|| id.to_string(), |n| n.to_string()),
            namespace: ns.get(*id).copied().unwrap_or(""),
            k,
            n,
            big_k,
            big_n,
            odds: (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)),
            ln_p,
            q: q[i],
            genes: strongest
                .iter()
                .take(5)
                .map(|&r| symbol(&genes[r].0).to_string())
                .collect(),
        });
        seen[g].push(overlap);
    }
    Ok(Table {
        source: source.said(),
        genes: genes.len(),
        universe: big_n,
        min_weight,
        modules: out,
    })
}

/// The 2 x 2 table of a module of `n` genes against a term of `big_k`, `k`
/// shared, in a universe of `big_n`: in both, the module's alone, the
/// term's alone, neither.
fn cells(k: usize, n: usize, big_k: usize, big_n: usize) -> [u64; 4] {
    [k, n - k, big_k - k, big_n - big_k - n + k].map(|x| x as u64)
}

/// Benjamini-Hochberg q-values of `ln_p` (natural-log p-values), in their
/// order.
fn bh(ln_p: &[f64]) -> Vec<f64> {
    let m = ln_p.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| ln_p[a].total_cmp(&ln_p[b]));
    let mut q = vec![1.0; m];
    let mut least = 1.0f64;
    for (rank, &i) in order.iter().enumerate().rev() {
        let v = (ln_p[i].exp() * m as f64 / (rank + 1) as f64).min(1.0);
        least = least.min(v);
        q[i] = least;
    }
    q
}

impl Table {
    /// The table as text: a head of what was tested, then a line per
    /// (module, term), tab-separated; a module with none says why.
    pub fn lines(&self) -> Vec<String> {
        let mut out = vec![
            format!(
                "go: {}; {} of the fit's {} genes annotated (the universe); a module's genes \
                 are those of main role weight at least {}; terms of {}..={} genes; q: BH over \
                 every (module, term) tested",
                self.source,
                self.universe,
                self.genes,
                self.min_weight,
                SIZES.start(),
                SIZES.end()
            ),
            "module\tgenes\trank\tterm\tnamespace\tname\tk\tK\todds\tp\tq\toverlap".into(),
        ];
        for (role, m) in self.modules.iter().enumerate() {
            let head = format!("M{}\t{}/{}", role + 1, m.n, m.all);
            if m.hits.is_empty() {
                out.push(format!("{head}\t-\t({})", m.why_none()));
            }
            for (r, h) in m.hits.iter().enumerate() {
                out.push(format!(
                    "{head}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2e}\t{:.2e}\t{}",
                    r + 1,
                    h.id,
                    h.namespace,
                    h.name,
                    h.k,
                    h.big_k,
                    h.odds,
                    h.ln_p.exp(),
                    h.q,
                    h.genes.join(",")
                ));
            }
        }
        out
    }

    /// Gene role `role`'s terms.
    pub fn of(&self, role: usize) -> Option<&ModuleTerms> {
        self.modules.get(role)
    }
}

impl ModuleTerms {
    /// Why it lists no term.
    pub fn why_none(&self) -> String {
        match self.n < LEAST_MODULE {
            true => format!("fewer than {LEAST_MODULE} genes to test"),
            false => format!("no term at q <= {Q} (BH over every module and term)"),
        }
    }
}

#[cfg(test)]
mod tests;
