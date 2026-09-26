//! The editing site tables (`m6a_sites.parquet`, `atoi_sites.parquet`): read
//! them, decide each site under the [`SiteFilterArgs`] thresholds, and write
//! the kept and the dropped rows back out, the dropped ones with a `reason`.
//!
//! The producers write every putative site; this is the only place a site is
//! ever rejected, and it rejects on stored columns, so a different cut is a
//! rerun of `faba qc`, never of a BAM scan.

use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, Float32Array, Int64Array, StringArray, UInt64Array,
};
use arrow::compute::{concat_batches, filter_record_batch};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::{RecordBatch, RecordBatchReader};
use data_beans::aux::feature_rows::{parse_feature_row, EDITED, METHYLATED};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::{FxHashMap, FxHashSet};

use super::args::SiteFilterArgs;
use crate::editing::io::write_record_batch;

/// One editing modality's site table, with everything a cut reads computed
/// once, per site, at load time.
pub struct SiteTable {
    pub modality: Box<str>,
    pub batch: RecordBatch,
    /// `{chr}:{pos}`, the subunit the `_site` matrix rows carry: the
    /// conversion position for m6A, the editing position for A-to-I, i.e.
    /// [`crate::editing::ConversionSite::conversion_pos`].
    pub key: Vec<Box<str>>,
    /// Dense gene id per site, for counting genes without rehashing strings.
    pub gene_id: Vec<u32>,
    pub n_genes: usize,
    pub pv: Vec<f32>,
    pub coverage: Vec<u64>,
    pub converted: Vec<u64>,
    pub control_coverage: Vec<u64>,
    /// `converted / coverage` on the signal arm.
    pub edit_ratio: Vec<f32>,
    /// Signal edit rate over control edit rate; `+inf` at a clean control.
    pub fold: Vec<f32>,
    /// The RAW cross-product log odds ratio `ln((a_w·u_m)/(u_w·a_m))` from the
    /// stored counts, the value the retired producer guard tested. NOT the
    /// parquet's `log_odds`: that one is Haldane-corrected for reporting, and
    /// thresholding it would re-impose a minimum signal rate at thin controls
    /// (with `a_m = 0` a corrected guard passes only when the signal odds exceed
    /// `0.5/(n_control + 0.5)`). Raw, a clean control reads `+inf` and passes,
    /// and a variant converting equally in both arms reads exactly `0.0`.
    pub raw_log_odds: Vec<f64>,
}

impl SiteTable {
    pub fn len(&self) -> usize {
        self.key.len()
    }

    pub fn is_empty(&self) -> bool {
        self.key.is_empty()
    }

    /// Whether this modality carries a control arm (m6A does, A-to-I does not).
    pub fn has_control(&self) -> bool {
        &*self.modality == data_beans::aux::feature_rows::M6A
    }

    /// The converted-read channel of this modality's `_site` rows.
    pub fn converted_channel(&self) -> &'static str {
        if self.has_control() {
            METHYLATED
        } else {
            EDITED
        }
    }

    /// Kept cells per site, aligned to the table, from the pooled per-key
    /// counts of [`accumulate_site_cells`]; a key with no row counts 0.
    pub fn cells_per_site(&self, acc: &FxHashMap<Box<str>, usize>) -> Vec<usize> {
        self.key
            .iter()
            .map(|k| acc.get(k).copied().unwrap_or(0))
            .collect()
    }
}

fn col<'a, T: 'static>(batch: &'a RecordBatch, name: &str) -> anyhow::Result<&'a T> {
    batch
        .column_by_name(name)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "site table lacks the `{name}` column; it was written by an older faba \
                 (rerun the producer, which now writes every putative site with its statistics)"
            )
        })?
        .as_any()
        .downcast_ref::<T>()
        .ok_or_else(|| anyhow::anyhow!("`{name}` column has an unexpected type"))
}

pub fn read_site_table(path: &str, modality: &str) -> anyhow::Result<SiteTable> {
    let file = std::fs::File::open(path)?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;
    let schema = reader.schema();
    let batches: Vec<RecordBatch> = reader.collect::<Result<_, _>>()?;
    let batch = concat_batches(&schema, &batches)?;
    let n = batch.num_rows();

    let chr: &StringArray = col(&batch, "chr")?;
    let gene: &StringArray = col(&batch, "gene")?;
    let primary: &Int64Array = col(&batch, "primary_pos")?;
    let conversion: &Int64Array = col(&batch, "conversion_pos")?;
    let pv: &Float32Array = col(&batch, "pv")?;
    let coverage: &UInt64Array = col(&batch, "coverage")?;
    let converted: &UInt64Array = col(&batch, "converted")?;
    let control_coverage: &UInt64Array = col(&batch, "control_coverage")?;
    let control_converted: &UInt64Array = col(&batch, "control_converted")?;

    let mut key = Vec::with_capacity(n);
    let mut gene_id = Vec::with_capacity(n);
    let mut gene_index: FxHashMap<&str, u32> = FxHashMap::default();
    let mut edit_ratio = Vec::with_capacity(n);
    let mut fold = Vec::with_capacity(n);
    let mut raw_log_odds = Vec::with_capacity(n);
    for i in 0..n {
        let pos = if conversion.is_null(i) {
            primary.value(i)
        } else {
            conversion.value(i)
        };
        key.push(format!("{}:{}", chr.value(i), pos).into());
        let next = gene_index.len() as u32;
        gene_id.push(*gene_index.entry(gene.value(i)).or_insert(next));

        let (a_w, n_w) = (converted.value(i), coverage.value(i));
        let (a_m, n_m) = (control_converted.value(i), control_coverage.value(i));
        let rate_w = if n_w == 0 {
            0.0
        } else {
            a_w as f32 / n_w as f32
        };
        let rate_m = if n_m == 0 {
            0.0
        } else {
            a_m as f32 / n_m as f32
        };
        edit_ratio.push(rate_w);
        fold.push(if rate_m <= 0.0 {
            f32::INFINITY
        } else {
            rate_w / rate_m
        });
        raw_log_odds.push(faba::hypothesis_tests::log_odds_ratio(
            a_w,
            n_w.saturating_sub(a_w),
            a_m,
            n_m.saturating_sub(a_m),
        ));
    }
    let n_genes = gene_index.len();
    Ok(SiteTable {
        modality: modality.into(),
        key,
        gene_id,
        n_genes,
        pv: pv.values().to_vec(),
        coverage: coverage.values().to_vec(),
        converted: converted.values().to_vec(),
        control_coverage: control_coverage.values().to_vec(),
        edit_ratio,
        fold,
        raw_log_odds,
        batch,
    })
}

/// Why a site was dropped, written to `{modality}_sites_dropped.parquet`.
/// The first failing check in a fixed order is recorded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DropReason {
    Coverage,
    Converted,
    EditRatio,
    Fold,
    LogOdds,
    Pvalue,
    Cells,
}

impl DropReason {
    pub fn label(&self) -> &'static str {
        match self {
            DropReason::Coverage => "coverage",
            DropReason::Converted => "converted",
            DropReason::EditRatio => "edit_ratio",
            DropReason::Fold => "fold",
            DropReason::LogOdds => "log_odds",
            DropReason::Pvalue => "pvalue",
            DropReason::Cells => "cells",
        }
    }
}

/// One threshold of [`SiteFilterArgs`]. The variants are in the order
/// [`SiteFilterArgs::reason`] checks them, so `c as usize` indexes
/// [`Criterion::ALL`] and the lowest failing one is the recorded reason.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criterion {
    MinCoverage,
    MinConverted,
    MinEditRatio,
    MaxEditRatio,
    MinFold,
    MinLogOdds,
    MaxPv,
    MinCells,
}

impl Criterion {
    pub const ALL: [Criterion; 8] = [
        Criterion::MinCoverage,
        Criterion::MinConverted,
        Criterion::MinEditRatio,
        Criterion::MaxEditRatio,
        Criterion::MinFold,
        Criterion::MinLogOdds,
        Criterion::MaxPv,
        Criterion::MinCells,
    ];

    /// The `faba qc` flag that sets it.
    pub fn flag(self) -> &'static str {
        match self {
            Criterion::MinCoverage => "--site-min-coverage",
            Criterion::MinConverted => "--site-min-converted",
            Criterion::MinEditRatio => "--site-min-edit-ratio",
            Criterion::MaxEditRatio => "--site-max-edit-ratio",
            Criterion::MinFold => "--site-min-fold",
            Criterion::MinLogOdds => "--site-min-log-odds",
            Criterion::MaxPv => "--site-max-pv",
            Criterion::MinCells => "--site-min-cells",
        }
    }

    /// The written `reason`; both edit-ratio bounds are `edit_ratio`.
    pub fn drop_reason(self) -> DropReason {
        match self {
            Criterion::MinCoverage => DropReason::Coverage,
            Criterion::MinConverted => DropReason::Converted,
            Criterion::MinEditRatio | Criterion::MaxEditRatio => DropReason::EditRatio,
            Criterion::MinFold => DropReason::Fold,
            Criterion::MinLogOdds => DropReason::LogOdds,
            Criterion::MaxPv => DropReason::Pvalue,
            Criterion::MinCells => DropReason::Cells,
        }
    }

    /// Whether the threshold is an upper bound on the raw value.
    pub fn is_max(self) -> bool {
        matches!(self, Criterion::MaxPv | Criterion::MaxEditRatio)
    }

    /// Whether the check runs at all: fold and log odds need a control arm,
    /// cells a `_site` matrix.
    pub fn applies(self, t: &SiteTable, has_cells: bool) -> bool {
        match self {
            Criterion::MinFold | Criterion::MinLogOdds => t.has_control(),
            Criterion::MinCells => has_cells,
            _ => true,
        }
    }

    /// Whether site `i` fails this check under `f`.
    pub fn fails(
        self,
        f: &SiteFilterArgs,
        t: &SiteTable,
        i: usize,
        n_cells: Option<usize>,
    ) -> bool {
        match self {
            Criterion::MinCoverage => t.coverage[i] + t.control_coverage[i] < f.site_min_coverage,
            Criterion::MinConverted => t.converted[i] < f.site_min_converted,
            Criterion::MinEditRatio => t.edit_ratio[i] < f.site_min_edit_ratio,
            Criterion::MaxEditRatio => t.edit_ratio[i] > f.site_max_edit_ratio,
            Criterion::MinFold => t.fold[i] < f.site_min_fold,
            Criterion::MinLogOdds => t.raw_log_odds[i] < f.site_min_log_odds as f64,
            Criterion::MaxPv => t.pv[i] > f.site_max_pv,
            Criterion::MinCells => n_cells.is_some_and(|n| n < f.site_min_cells),
        }
    }

    /// The raw value this check tests at site `i`.
    pub fn raw(self, t: &SiteTable, i: usize, n_cells: Option<usize>) -> f64 {
        match self {
            Criterion::MinCoverage => (t.coverage[i] + t.control_coverage[i]) as f64,
            Criterion::MinConverted => t.converted[i] as f64,
            Criterion::MinEditRatio | Criterion::MaxEditRatio => t.edit_ratio[i] as f64,
            Criterion::MinFold => t.fold[i] as f64,
            Criterion::MinLogOdds => t.raw_log_odds[i],
            Criterion::MaxPv => t.pv[i] as f64,
            Criterion::MinCells => n_cells.unwrap_or(0) as f64,
        }
    }

    pub fn get(self, f: &SiteFilterArgs) -> f64 {
        match self {
            Criterion::MinCoverage => f.site_min_coverage as f64,
            Criterion::MinConverted => f.site_min_converted as f64,
            Criterion::MinEditRatio => f.site_min_edit_ratio as f64,
            Criterion::MaxEditRatio => f.site_max_edit_ratio as f64,
            Criterion::MinFold => f.site_min_fold as f64,
            Criterion::MinLogOdds => f.site_min_log_odds as f64,
            Criterion::MaxPv => f.site_max_pv as f64,
            Criterion::MinCells => f.site_min_cells as f64,
        }
    }

    /// Set the threshold to `raw`, rounded so that a site whose value is
    /// exactly `raw` stays kept: a lower bound rounds down, an upper bound up.
    pub fn set(self, f: &mut SiteFilterArgs, raw: f64) {
        let as_f32 = |x: f64| {
            let y = x as f32;
            if self.is_max() && (y as f64) < x {
                y.next_up()
            } else if !self.is_max() && (y as f64) > x {
                y.next_down()
            } else {
                y
            }
        };
        let as_count = |x: f64| x.max(0.0).ceil();
        match self {
            Criterion::MinCoverage => f.site_min_coverage = as_count(raw) as u64,
            Criterion::MinConverted => f.site_min_converted = as_count(raw) as u64,
            Criterion::MinEditRatio => f.site_min_edit_ratio = as_f32(raw),
            Criterion::MaxEditRatio => f.site_max_edit_ratio = as_f32(raw),
            Criterion::MinFold => f.site_min_fold = as_f32(raw),
            Criterion::MinLogOdds => f.site_min_log_odds = as_f32(raw),
            Criterion::MaxPv => f.site_max_pv = as_f32(raw),
            Criterion::MinCells => f.site_min_cells = as_count(raw) as usize,
        }
    }

    pub fn permissive(self) -> f64 {
        self.get(&SiteFilterArgs::permissive())
    }

    /// Whether `f` sets this threshold where it keeps everything.
    pub fn is_off(self, f: &SiteFilterArgs) -> bool {
        let (v, off) = (self.get(f), self.permissive());
        if self.is_max() {
            v >= off
        } else {
            v <= off
        }
    }
}

impl SiteFilterArgs {
    /// Thresholds that keep everything; `qc-report` starts from it and moves
    /// one knob at a time.
    pub fn permissive() -> Self {
        SiteFilterArgs {
            site_max_pv: 1.0,
            site_min_log_odds: f32::NEG_INFINITY,
            site_min_coverage: 0,
            site_min_converted: 0,
            site_min_edit_ratio: 0.0,
            site_max_edit_ratio: f32::INFINITY,
            site_min_fold: 0.0,
            site_min_cells: 0,
        }
    }

    /// Parse `--site-*` flags as `faba qc` would, defaults filling the rest.
    /// Panics (through clap) on a flag it does not know.
    pub fn parse_flags<'a>(flags: impl IntoIterator<Item = &'a str>) -> Self {
        #[derive(clap::Parser)]
        struct Flags {
            #[command(flatten)]
            site: SiteFilterArgs,
        }
        <Flags as clap::Parser>::parse_from(std::iter::once("faba").chain(flags)).site
    }

    /// The command-line defaults, as `faba qc` uses them with no `--site-*`
    /// flag given.
    pub fn default_values() -> Self {
        Self::parse_flags([])
    }

    /// Verdict for site `i`: the first failing [`Criterion`] in check order.
    /// `n_cells` is the site's kept-cell count, or `None` when no `_site`
    /// matrix exists to read it from (the cells check is then skipped).
    pub fn reason(&self, t: &SiteTable, i: usize, n_cells: Option<usize>) -> Option<DropReason> {
        Criterion::ALL
            .into_iter()
            .find(|c| c.applies(t, n_cells.is_some()) && c.fails(self, t, i, n_cells))
            .map(Criterion::drop_reason)
    }

    /// Verdicts for every site; `n_cells` is aligned to the table.
    pub fn reasons(&self, t: &SiteTable, n_cells: Option<&[usize]>) -> Vec<Option<DropReason>> {
        (0..t.len())
            .map(|i| self.reason(t, i, n_cells.map(|c| c[i])))
            .collect()
    }
}

/// Per-site kept-cell counts from one `_site` matrix: the nnz (over kept
/// cells) of each site's converted-channel row, keyed by the `{chr}:{pos}`
/// subunit. Summed into `acc` so batches pool.
pub fn accumulate_site_cells(
    row_names: &[Box<str>],
    row_nnz: &[usize],
    converted_channel: &str,
    acc: &mut FxHashMap<Box<str>, usize>,
) {
    for (name, &nnz) in row_names.iter().zip(row_nnz) {
        let Some(row) = parse_feature_row(name) else {
            continue;
        };
        if row.channel != converted_channel {
            continue;
        }
        if let Some(sub) = row.subunit {
            match acc.get_mut(sub) {
                Some(v) => *v += nnz,
                None => {
                    acc.insert(sub.into(), nnz);
                }
            }
        }
    }
}

/// Ascending row ids of a `_site` matrix whose subunit is a kept site (both
/// channels of a kept site survive together); none when no site table exists.
pub fn site_matrix_rows(row_names: &[Box<str>], kept: Option<&FxHashSet<Box<str>>>) -> Vec<usize> {
    let Some(kept) = kept else {
        return Vec::new();
    };
    row_names
        .iter()
        .enumerate()
        .filter(|(_, name)| {
            parse_feature_row(name)
                .and_then(|r| r.subunit)
                .is_some_and(|s| kept.contains(s))
        })
        .map(|(i, _)| i)
        .collect()
}

/// Write the kept rows to `kept_path` and the dropped rows, plus a `reason`
/// column, to `dropped_path`. Returns `(n_kept, n_dropped)`.
pub fn write_site_tables(
    t: &SiteTable,
    reasons: &[Option<DropReason>],
    kept_path: &str,
    dropped_path: &str,
) -> anyhow::Result<(usize, usize)> {
    let keep: Vec<bool> = reasons.iter().map(|r| r.is_none()).collect();
    let drop: Vec<bool> = keep.iter().map(|k| !k).collect();

    let kept = filter_record_batch(&t.batch, &BooleanArray::from(keep))?;
    write_record_batch(&kept, kept_path)?;

    let dropped = filter_record_batch(&t.batch, &BooleanArray::from(drop))?;
    let labels: Vec<&str> = reasons.iter().flatten().map(|r| r.label()).collect();
    let reason_array = Arc::new(StringArray::from(labels)) as ArrayRef;
    let mut fields: Vec<Field> = dropped
        .schema()
        .fields()
        .iter()
        .map(|f| f.as_ref().clone())
        .collect();
    fields.push(Field::new("reason", DataType::Utf8, false));
    let mut columns = dropped.columns().to_vec();
    columns.push(reason_array);
    let dropped = RecordBatch::try_new(Arc::new(Schema::new(fields)), columns)?;
    write_record_batch(&dropped, dropped_path)?;

    Ok((kept.num_rows(), dropped.num_rows()))
}
