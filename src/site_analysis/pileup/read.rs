//! Reading the selected rows: sparse matrices, site parquets, read depth.

use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Default)]
struct PosAgg {
    pub(super) sum: f64,
    pub(super) nnz: usize,
}

pub(super) struct MatrixGeneData {
    pub(super) gene: Box<str>,
    /// The matched gene keys, sorted.
    pub(super) genes: Vec<Box<str>>,
    pub(super) chr: Box<str>,
    pub(super) modality: Box<str>,
    pub(super) total: Option<Vec<(i64, f64)>>,
    pub(super) positions: Vec<(i64, f64)>,
}

/// Gene boundaries and per-site signals from the sites parquet.
pub(super) struct SiteAnnotation {
    pub(super) gene_start: i64,
    pub(super) gene_stop: i64,
    pub(super) positions: Vec<(i64, f64)>,
    pub(super) num_sites: usize,
}

/// How a modality's two channels are named on screen: A-to-I reads are
/// converted or not, m6A reads methylated or not.
pub(crate) fn channel_names(modality: &str) -> (&'static str, &'static str) {
    if modality.eq_ignore_ascii_case(ATOI) {
        ("converted", "unconverted")
    } else {
        (METHYLATED, UNMETHYLATED)
    }
}

pub(super) fn read_site_annotation(
    site_file: &str,
    selector: &Selector,
    site_signal: &SiteSignal,
) -> anyhow::Result<SiteAnnotation> {
    let file = std::fs::File::open(site_file)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut gene_start: Option<i64> = None;
    let mut gene_stop: Option<i64> = None;
    let mut positions: Vec<(i64, f64)> = Vec::new();
    let mut distinct_genes: FxHashMap<Box<str>, usize> = FxHashMap::default();

    // Resolve signal column names once
    let signal_cols: &[&str] = match site_signal {
        SiteSignal::WtCoverage => &["wt_a", "wt_t", "wt_g", "wt_c"],
        SiteSignal::MutCoverage => &["mut_a", "mut_t", "mut_g", "mut_c"],
        _ => &[],
    };

    for batch in reader {
        let batch = batch?;

        let gene_col = batch
            .column_by_name("gene")
            .ok_or_else(|| anyhow::anyhow!("missing 'gene' column in {}", site_file))?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| anyhow::anyhow!("'gene' column is not a string array"))?;

        let pos_col = batch
            .column_by_name("primary_pos")
            .ok_or_else(|| anyhow::anyhow!("missing 'primary_pos' column in {}", site_file))?
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or_else(|| anyhow::anyhow!("'primary_pos' column is not Int64"))?;

        // Optional: enables `--regions` filtering on the parquet track.
        let chr_col = batch
            .column_by_name("chr")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>());

        let start_col = batch
            .column_by_name("gene_start")
            .and_then(|c| c.as_any().downcast_ref::<Int64Array>());

        let stop_col = batch
            .column_by_name("gene_stop")
            .and_then(|c| c.as_any().downcast_ref::<Int64Array>());

        // Resolve signal columns once per batch
        let pv_col = batch
            .column_by_name("pv")
            .and_then(|c| c.as_any().downcast_ref::<Float32Array>());

        let base_cols: Vec<Option<&UInt64Array>> = signal_cols
            .iter()
            .map(|name| {
                batch
                    .column_by_name(name)
                    .and_then(|c| c.as_any().downcast_ref::<UInt64Array>())
            })
            .collect();

        for i in 0..batch.num_rows() {
            let gene_val = gene_col.value(i);
            let pos = pos_col.value(i);
            let chr_val = chr_col.map(|c| c.value(i)).unwrap_or("");
            if !selector.selects(gene_val, chr_val, pos) {
                continue;
            }

            *distinct_genes.entry(gene_val.into()).or_insert(0) += 1;

            // Update gene boundaries
            if let (Some(sc), Some(tc)) = (start_col, stop_col) {
                let gs = sc.value(i);
                let gt = tc.value(i);
                gene_start = Some(gene_start.map_or(gs, |v: i64| v.min(gs)));
                gene_stop = Some(gene_stop.map_or(gt, |v: i64| v.max(gt)));
            }

            let value = match site_signal {
                SiteSignal::Count => 1.0,
                SiteSignal::WtCoverage | SiteSignal::MutCoverage => {
                    let mut total = 0u64;
                    for col in base_cols.iter().flatten() {
                        total += col.value(i);
                    }
                    total as f64
                }
                SiteSignal::NegLog10Pv => {
                    if let Some(pv) = pv_col {
                        let p = pv.value(i);
                        if p > 0.0 {
                            -(p as f64).log10()
                        } else {
                            300.0
                        }
                    } else {
                        0.0
                    }
                }
            };

            positions.push((pos, value));
        }
    }

    if positions.is_empty() {
        return Err(anyhow::anyhow!(
            "no sites matching {} in {}",
            selector.describe(),
            site_file
        ));
    }

    // Aggregate every matched gene's sites; gene boundaries below widen to
    // the span across all of them rather than erroring on ambiguity.
    if distinct_genes.len() > 1 {
        info!(
            "{} matched {} genes in parquet; aggregating all sites",
            selector.describe(),
            distinct_genes.len()
        );
    }

    positions.sort_unstable_by_key(|(pos, _)| *pos);

    let gs = gene_start.unwrap_or_else(|| positions.first().unwrap().0);
    let gt = gene_stop.unwrap_or_else(|| positions.last().unwrap().0);

    info!(
        "loaded {} sites from parquet, gene boundaries: {}-{}",
        positions.len(),
        gs,
        gt
    );

    let num_sites = positions.len();
    Ok(SiteAnnotation {
        gene_start: gs,
        gene_stop: gt,
        positions,
        num_sites,
    })
}

/// Matrix positions grouped by cell type (the figure's per-panel top
/// track). With `membership = None` every cell falls into a single
/// synthetic `""` group, which is exactly the all-cells aggregate the
/// ASCII path uses — so [`read_matrix_positions`] is a thin wrapper.
pub(super) struct GroupedMatrix {
    pub(super) gene: Box<str>,
    pub(super) chr: Box<str>,
    /// The first matched row's modality (e.g. `m6a`).
    pub(super) modality: Box<str>,
    /// celltype label -> sorted `(pos, value)`
    pub(super) by_group: FxHashMap<Box<str>, Vec<(i64, f64)>>,
    /// Both channels summed per position (all cells), when asked for and
    /// the matrices have an unconverted channel.
    pub(super) total: Option<Vec<(i64, f64)>>,
    /// Converted rows the selection matched, over all files.
    pub(super) matched: usize,
    /// The matched gene keys, sorted.
    pub(super) genes: Vec<Box<str>>,
}

/// One file's share of [`read_matrix_positions_grouped`], merged with the
/// other files' in file order.
#[derive(Default)]
struct FileRead {
    by_group: FxHashMap<Box<str>, FxHashMap<i64, PosAgg>>,
    distinct_genes: FxHashMap<Box<str>, usize>,
    matched_chrs: Vec<Box<str>>,
    matched: usize,
    first_modality: Option<Box<str>>,
    totals: FxHashMap<i64, f64>,
    unconverted_rows: usize,
}

impl FileRead {
    /// Add `other`, a later file's share: per group and position the sums
    /// add, and the first file's modality stays.
    fn merge(&mut self, other: FileRead) {
        if self.matched == 0 && self.unconverted_rows == 0 && self.by_group.is_empty() {
            // Nothing yet: take it whole rather than copy it over.
            *self = other;
            return;
        }
        for (grp, aggs) in other.by_group {
            let into = self.by_group.entry(grp).or_default();
            for (pos, agg) in aggs {
                let a = into.entry(pos).or_default();
                a.sum += agg.sum;
                a.nnz += agg.nnz;
            }
        }
        for (gene, n) in other.distinct_genes {
            *self.distinct_genes.entry(gene).or_insert(0) += n;
        }
        self.matched_chrs.extend(other.matched_chrs);
        self.matched += other.matched;
        if self.first_modality.is_none() {
            self.first_modality = other.first_modality;
        }
        for (pos, v) in other.totals {
            *self.totals.entry(pos).or_insert(0.0) += v;
        }
        self.unconverted_rows += other.unconverted_rows;
    }
}

/// Read the rows `selector` picks from `data_file`.
fn read_one_file(
    data_file: &str,
    selector: &Selector,
    membership: Option<&CellMembership>,
    modality_filter: Option<&[Box<str>]>,
    with_total: bool,
) -> anyhow::Result<FileRead> {
    let mut part = FileRead::default();
    let (backend, resolved_path) = resolve_backend_file(data_file, None)?;
    let data = open_sparse_matrix(&resolved_path, &backend)?;

    let row_names = data.row_names()?;

    let mut matched_rows: Vec<(usize, i64, bool)> = Vec::new();
    for (idx, name) in row_names.iter().enumerate() {
        if let Some((gene_part, modality, chr, pos, converted)) = parse_row_channel(name) {
            // Unconverted rows count only toward the total.
            let is_converted = converted != Some(false);
            if !is_converted && !with_total {
                continue;
            }
            if let Some(mf) = modality_filter {
                if !mf.iter().any(|m| m.eq_ignore_ascii_case(modality)) {
                    continue;
                }
            }
            if selector.selects(gene_part, chr, pos) {
                if is_converted {
                    part.first_modality.get_or_insert_with(|| modality.into());
                    *part.distinct_genes.entry(gene_part.into()).or_insert(0) += 1;
                    part.matched_chrs.push(chr.into());
                } else {
                    part.unconverted_rows += 1;
                }
                matched_rows.push((idx, pos, is_converted));
            }
        }
    }

    if matched_rows.is_empty() {
        return Ok(part);
    }
    part.matched += matched_rows.iter().filter(|r| r.2).count();

    // Column index -> cell type (None to drop). Only needed when
    // stratifying; the all-cells path skips reading column names.
    let col_groups: Option<Vec<Option<Box<str>>>> = match membership {
        None => None,
        Some(m) => {
            let col_names = data.column_names()?;
            Some(
                col_names
                    .iter()
                    .map(|bc| m.matches_barcode(&CellBarcode::Barcode(Arc::from(bc.as_ref()))))
                    .collect(),
            )
        }
    };

    let local_to_pos: Vec<(i64, bool)> = matched_rows.iter().map(|r| (r.1, r.2)).collect();
    let row_indices: Vec<usize> = matched_rows.iter().map(|r| r.0).collect();
    let (_nrow, _ncol, triplets) = data.read_triplets_by_rows(row_indices)?;

    // All-cells: pre-seed every matched position so zero-signal sites
    // still appear (axis markers), matching the legacy behavior.
    if membership.is_none() {
        let g = part.by_group.entry("".into()).or_default();
        for &(pos, is_converted) in &local_to_pos {
            if is_converted {
                g.entry(pos).or_default();
            }
        }
    }

    for (row, col, val) in &triplets {
        let local_idx = *row as usize;
        if local_idx >= local_to_pos.len() || *val == 0.0 {
            continue;
        }
        let group: Box<str> = match &col_groups {
            None => "".into(),
            Some(cg) => match cg.get(*col as usize).and_then(|o| o.clone()) {
                Some(ct) => ct,
                None => continue,
            },
        };
        let (pos, is_converted) = local_to_pos[local_idx];
        if with_total {
            *part.totals.entry(pos).or_insert(0.0) += *val as f64;
        }
        if !is_converted {
            continue;
        }
        let agg = part
            .by_group
            .entry(group)
            .or_default()
            .entry(pos)
            .or_default();
        agg.sum += *val as f64;
        agg.nnz += 1;
    }
    Ok(part)
}

pub(super) fn read_matrix_positions_grouped(
    data_files: &[Box<str>],
    selector: &Selector,
    signal: &PileupSignal,
    membership: Option<&CellMembership>,
    top_modality: &[Box<str>],
    with_total: bool,
    done: Option<&AtomicUsize>,
) -> anyhow::Result<GroupedMatrix> {
    // Only the rows of these modalities, when some are named.
    let modality_filter = (!top_modality.is_empty()).then_some(top_modality);

    // The files are independent: read them together, then merge in file
    // order. Multiple input files (e.g. replicates) merge per genomic
    // position; gene/chr labels reflect the union.
    let parts: Vec<anyhow::Result<FileRead>> = {
        use rayon::prelude::*;
        data_files
            .par_iter()
            .map(|f| {
                let part = read_one_file(f, selector, membership, modality_filter, with_total);
                if let Some(done) = done {
                    done.fetch_add(1, Ordering::Relaxed);
                }
                part
            })
            .collect()
    };
    let mut all = FileRead::default();
    for part in parts {
        all.merge(part?);
    }
    let FileRead {
        by_group,
        distinct_genes,
        matched_chrs,
        matched: total_matched,
        first_modality,
        totals,
        unconverted_rows,
    } = all;

    // The selection can span several genes (and chromosomes); rather than
    // erroring on ambiguity we pile them together and label the aggregate.
    let gene: Box<str> = summarize_genes(&distinct_genes);
    let mut genes: Vec<Box<str>> = distinct_genes.keys().cloned().collect();
    genes.sort_unstable();
    let chr: Box<str> = summarize_chr(&matched_chrs);

    info!(
        "matched {} rows across {} gene(s) in {} file(s) for {} [{}]",
        total_matched,
        distinct_genes.len(),
        data_files.len(),
        selector.describe(),
        gene
    );

    let by_group = by_group
        .into_iter()
        .map(|(grp, pos_agg)| {
            let mut positions: Vec<(i64, f64)> = pos_agg
                .iter()
                .map(|(&pos, agg)| {
                    let value = match signal {
                        PileupSignal::Sum | PileupSignal::Log10Sum => agg.sum,
                        PileupSignal::Nnz => agg.nnz as f64,
                    };
                    (pos, value)
                })
                .collect();
            positions.sort_unstable_by_key(|(pos, _)| *pos);
            (grp, positions)
        })
        .collect();

    // No unconverted rows (a matrix without that channel): no total.
    let total = (unconverted_rows > 0).then(|| {
        let mut total: Vec<(i64, f64)> = totals.into_iter().collect();
        total.sort_unstable_by_key(|p| p.0);
        total
    });
    Ok(GroupedMatrix {
        gene,
        chr,
        modality: first_modality.unwrap_or_default(),
        by_group,
        total,
        matched: total_matched,
        genes,
    })
}

/// The selection's converted positions over `data_files` (and with
/// `with_total` both channels summed), or `None` when no row matched.
pub(super) fn read_matrix_positions(
    data_files: &[Box<str>],
    selector: &Selector,
    signal: &PileupSignal,
    with_total: bool,
    done: Option<&AtomicUsize>,
) -> anyhow::Result<Option<MatrixGeneData>> {
    let grouped =
        read_matrix_positions_grouped(data_files, selector, signal, None, &[], with_total, done)?;
    if grouped.matched == 0 {
        return Ok(None);
    }
    // membership = None yields exactly one synthetic "" group.
    let positions = grouped
        .by_group
        .into_iter()
        .next()
        .map(|(_, p)| p)
        .unwrap_or_default();
    Ok(Some(MatrixGeneData {
        gene: grouped.gene,
        genes: grouped.genes,
        chr: grouped.chr,
        modality: grouped.modality,
        positions,
        total: grouped.total,
    }))
}

/// Parse a `_depth` row name `chr:start-end`, by the shared strict locus rule.
pub(super) fn parse_depth_row(name: &str) -> Option<(&str, i64, i64)> {
    genomic_data::coordinates::split_interval(name)
}

/// Read depth over `chr:lo-hi` from `_depth` matrices: each overlapping
/// bin `(start, end, reads)`, reads summed over cells and files, sorted.
pub(super) fn read_depth(
    files: &[Box<str>],
    chr: &str,
    (lo, hi): (i64, i64),
) -> anyhow::Result<Vec<(i64, i64, f64)>> {
    let mut bins: FxHashMap<(i64, i64), f64> = FxHashMap::default();
    for file in files {
        let data = crate::qc::matrix::open_matrix(file)?;
        let names = data.row_names()?;
        let rows: Vec<(usize, (i64, i64))> = names
            .iter()
            .enumerate()
            .filter_map(|(i, n)| parse_depth_row(n).map(|(c, s, e)| (i, c, s, e)))
            .filter(|&(_, c, s, e)| chr_eq(c, chr) && e > lo && s <= hi)
            .map(|(i, _, s, e)| (i, (s, e)))
            .collect();
        if rows.is_empty() {
            continue;
        }
        let (_, _, triplets) = data.read_triplets_by_rows(rows.iter().map(|r| r.0).collect())?;
        for bin in rows.iter().map(|r| r.1) {
            bins.entry(bin).or_insert(0.0);
        }
        for (row, _, value) in triplets {
            if let Some(&(_, bin)) = rows.get(row as usize) {
                *bins.entry(bin).or_insert(0.0) += value as f64;
            }
        }
    }
    let mut out: Vec<(i64, i64, f64)> = bins.into_iter().map(|((s, e), v)| (s, e, v)).collect();
    out.sort_unstable_by_key(|b| b.0);
    Ok(out)
}
