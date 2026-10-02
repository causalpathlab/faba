//! Row names, gene and region queries, and the selector built from them.

use super::*;

/// Parse a faba row name `gene_key/modality/detail`. `detail` is either
/// `chr:pos` (site output, e.g. `ENSG00000000003_GENE3/m6A/chr13:32350000`)
/// or a bare component ordinal (mixture output, e.g.
/// `ENSG00000000003_GENE3/m6A/0`). Returns `(gene_part, chr, x)` where
/// `chr` is empty for mixture rows and `x` is the genomic position or the
/// component ordinal used as the pileup x-coordinate.
/// Returns `(gene, modality, chr, pos)`. `modality` is the middle
/// `/`-delimited token (e.g. `m6A`), used by the figure's
/// `--top-modality` filter; `chr` is empty for mixture (component) rows.
///
/// Site rows written with a channel suffix (`.../chr:pos/methylated`) pile
/// up their converted channel only: the suffix is dropped from `methylated`
/// and `edited` rows, and `unmethylated` / `unedited` rows yield `None`.
pub(super) fn parse_row_name_full(name: &str) -> Option<(&str, &str, &str, i64)> {
    parse_row_channel(name)
        .filter(|row| row.4 != Some(false))
        .map(|(gene, modality, chr, pos, _)| (gene, modality, chr, pos))
}

/// [`parse_row_name_full`] for either channel: whether the row is the
/// modality's first (converted) channel, `None` for rows without one.
/// Rows split with [`parse_feature_row`], so a gene symbol may hold `/`.
pub(super) fn parse_row_channel(name: &str) -> Option<(&str, &str, &str, i64, Option<bool>)> {
    let row = parse_feature_row(name)?;
    let (detail, converted) = match row.subunit {
        Some(detail) => {
            let (on, off) =
                channels(row.modality).or_else(|| channels(&row.modality.to_ascii_lowercase()))?;
            let converted = match row.channel {
                c if c == on => true,
                c if c == off => false,
                _ => return None,
            };
            (detail, Some(converted))
        }
        None => (row.channel, None),
    };
    if let Some((chr, pos_str)) = detail.split_once(':') {
        let pos = pos_str.parse::<i64>().ok()?;
        Some((row.gene, row.modality, chr, pos, converted))
    } else {
        // Mixture rows carry a component ordinal, not a chromosome.
        let component = detail.parse::<i64>().ok()?;
        Some((row.gene, row.modality, "", component, converted))
    }
}

/// Relaxed gene matching, consistent with the data_beans::aux
/// `FeatureNameKind::Gene` canonicalization used for cross-file row
/// alignment. A row matches when its `gene_part` shares any `_`-split
/// component with the query, or agrees on the canonical gene symbol
/// (last `_`-delimited component, Cell Ranger feature-type suffix
/// stripped) — so both a symbol (`GENE1`) and an Ensembl ID
/// (`ENSG00000000001`) query resolve the same row. Case-insensitive.
/// `query_sym` is the pre-canonicalized query symbol.
pub(crate) fn gene_matches(query: &str, query_sym: &str, gene_part: &str) -> bool {
    // A position or locus (a BAF row's `chr:pos`) is not a gene name: only
    // the whole of it matches, so an alt contig's `_` pieces never do.
    if genomic_data::coordinates::is_region(gene_part) {
        return gene_part.eq_ignore_ascii_case(query);
    }
    // Allocation-free component check first — it directly covers symbol and
    // Ensembl-ID queries (and subsumes a full-composite match). Fall back to
    // the suffix-stripping canonicalizer only when the components miss.
    gene_part.split('_').any(|c| c.eq_ignore_ascii_case(query))
        || FeatureNameKind::Gene { delim: '_' }
            .canonicalize(gene_part)
            .eq_ignore_ascii_case(query_sym)
}

/// Canonical query symbol used by [`gene_matches`].
pub(crate) fn query_symbol(query: &str) -> Box<str> {
    FeatureNameKind::Gene { delim: '_' }.canonicalize(query)
}

/// A `chr:lb-ub` genomic window. Bounds are inclusive.
pub(crate) struct Region {
    pub(super) chr: Box<str>,
    pub(super) lb: i64,
    pub(super) ub: i64,
}

/// Parse `chr:lb-ub` (e.g. `chr17:1000-2000`). Reversed bounds are
/// swapped so `lb <= ub` always holds.
pub(super) fn parse_region(spec: &str) -> anyhow::Result<Region> {
    let bad = || anyhow::anyhow!("region '{}' must be formatted chr:lb-ub", spec);
    let (chr, range) = spec.split_once(':').ok_or_else(bad)?;
    let (lb, ub) = range.split_once('-').ok_or_else(bad)?;
    let lb: i64 = lb.trim().parse().map_err(|_| bad())?;
    let ub: i64 = ub.trim().parse().map_err(|_| bad())?;
    let (lb, ub) = if lb <= ub { (lb, ub) } else { (ub, lb) };
    let chr = chr.trim();
    if chr.is_empty() {
        return Err(bad());
    }
    Ok(Region {
        chr: chr.into(),
        lb,
        ub,
    })
}

/// Something typed into a `/` search or the gene list.
pub(crate) enum Query {
    /// A window, or a single position (`true`).
    Locus(Region, bool),
    Gene(Box<str>),
}

/// Parse `chr:start-end`, `chr:pos` (commas allowed) or a gene name.
pub(crate) fn parse_query(text: &str) -> Option<Query> {
    let t: String = text
        .chars()
        .filter(|c| !c.is_whitespace() && *c != ',')
        .collect();
    if t.is_empty() {
        return None;
    }
    match t.split_once(':') {
        Some((_, range)) if range.contains('-') => {
            parse_region(&t).ok().map(|r| Query::Locus(r, false))
        }
        Some((chr, pos)) => pos
            .parse::<i64>()
            .ok()
            .filter(|_| !chr.is_empty())
            .map(|p| {
                let r = Region {
                    chr: chr.into(),
                    lb: p,
                    ub: p,
                };
                Query::Locus(r, true)
            }),
        None => Some(Query::Gene(t.into())),
    }
}

/// Combined gene + region row selector. A row is selected when it
/// matches any requested gene OR falls inside any requested region
/// (union), so callers can pass either or both.
pub(crate) struct Selector {
    pub(super) genes: Vec<Box<str>>,
    pub(super) gene_syms: Vec<Box<str>>,
    pub(super) regions: Vec<Region>,
    /// Match `genes` as whole row keys rather than by the relaxed rules.
    pub(super) exact: bool,
}

impl Selector {
    pub(crate) fn build(genes: &[Box<str>], regions: &[Box<str>]) -> anyhow::Result<Self> {
        let genes: Vec<Box<str>> = genes
            .iter()
            .map(|g| g.trim())
            .filter(|g| !g.is_empty())
            .map(Into::into)
            .collect();
        let gene_syms: Vec<Box<str>> = genes.iter().map(|g| query_symbol(g)).collect();
        let regions: Vec<Region> = regions
            .iter()
            .map(|r| r.trim())
            .filter(|r| !r.is_empty())
            .map(parse_region)
            .collect::<anyhow::Result<_>>()?;
        if genes.is_empty() && regions.is_empty() {
            anyhow::bail!("provide at least one of --genes or --regions");
        }
        Ok(Self {
            genes,
            gene_syms,
            regions,
            exact: false,
        })
    }

    /// Exactly the gene whose row key is `gene_part`, as a picked catalog
    /// entry names it.
    pub(crate) fn exact(gene_part: &str) -> Self {
        Self {
            genes: vec![gene_part.into()],
            gene_syms: Vec::new(),
            regions: Vec::new(),
            exact: true,
        }
    }

    pub(crate) fn matches_gene(&self, gene_part: &str) -> bool {
        if self.exact {
            return self.genes.iter().any(|g| g.as_ref() == gene_part);
        }
        self.genes
            .iter()
            .zip(&self.gene_syms)
            .any(|(g, sym)| gene_matches(g, sym, gene_part))
    }

    fn matches_region(&self, chr: &str, pos: i64) -> bool {
        // Mixture rows carry no chromosome (empty `chr`); they can never sit
        // inside a region, so guard explicitly rather than relying on
        // `parse_region` having rejected empty region chromosomes.
        !chr.is_empty()
            && self
                .regions
                .iter()
                .any(|r| chr_eq(&r.chr, chr) && pos >= r.lb && pos <= r.ub)
    }

    /// Row is kept when it matches any gene or any region.
    pub(super) fn selects(&self, gene_part: &str, chr: &str, pos: i64) -> bool {
        self.matches_gene(gene_part) || self.matches_region(chr, pos)
    }

    /// Short description of the active selection for log/error messages.
    pub(super) fn describe(&self) -> String {
        let mut parts = Vec::new();
        if !self.genes.is_empty() {
            let g: Vec<&str> = self.genes.iter().map(|g| g.as_ref()).collect();
            parts.push(format!("genes [{}]", g.join(",")));
        }
        if !self.regions.is_empty() {
            let r: Vec<String> = self
                .regions
                .iter()
                .map(|r| format!("{}:{}-{}", r.chr, r.lb, r.ub))
                .collect();
            parts.push(format!("regions [{}]", r.join(",")));
        }
        parts.join(" + ")
    }
}

/// Human-readable label for the set of matched genes. A single gene is
/// shown verbatim; multiple matches collapse to `N genes: a,b,...`.
pub(super) fn summarize_genes(distinct: &FxHashMap<Box<str>, usize>) -> Box<str> {
    if distinct.len() == 1 {
        return distinct.keys().next().unwrap().clone();
    }
    let mut names: Vec<&str> = distinct.keys().map(|k| k.as_ref()).collect();
    names.sort_unstable();
    let shown = names.iter().take(5).copied().collect::<Vec<_>>().join(",");
    format!(
        "{} genes: {}{}",
        names.len(),
        shown,
        if names.len() > 5 { ",..." } else { "" }
    )
    .into()
}

/// Label for the chromosome axis. `component` when matched rows carry no
/// chromosome (mixture output), the single chromosome when all matches
/// agree, else `*` for a multi-chromosome aggregate.
pub(super) fn summarize_chr(matched_chrs: &[Box<str>]) -> Box<str> {
    let mut chrs: Vec<&str> = matched_chrs
        .iter()
        .map(|c| c.as_ref())
        .filter(|c| !c.is_empty())
        .collect();
    chrs.sort_unstable();
    chrs.dedup();
    match chrs.len() {
        0 => "component".into(),
        1 => chrs[0].into(),
        _ => "*".into(),
    }
}
