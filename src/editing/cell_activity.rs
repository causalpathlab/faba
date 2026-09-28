//! Null-cell calling: which WT cells actually carry APOBEC1-YTH editing.
//!
//! Expression-based cell QC passes cells with no functional fusion protein. This
//! stage (QC, not a test) cuts as aggressively as possible while the discarded
//! pool stays indistinguishable from the catalytically-dead control.

use crate::common::*;
use genomic_data::gff::GffRecordMap;
use genomic_data::sam::CellBarcode;
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;

pub mod scan;

/// A cell's genome-wide editing tally (converted out of covered reads), for either arm.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CellActivity {
    /// Reads showing the conversion (C→T forward, G→A reverse) at motif candidates.
    pub edited: u64,
    /// Reads covering a motif candidate with either the reference or converted base.
    pub covered: u64,
    /// Same, at same-channel bases away from motif candidates (promiscuous,
    /// m6A-independent deamination). Reported only, never used in the cut.
    pub bg_edited: u64,
    pub bg_covered: u64,
}

impl CellActivity {
    pub fn add(&mut self, edited: u64, covered: u64) {
        self.edited += edited;
        self.covered += covered;
    }

    pub fn add_background(&mut self, edited: u64, covered: u64) {
        self.bg_edited += edited;
        self.bg_covered += covered;
    }

    /// Promiscuous (non-motif) conversion rate, the APOBEC1-YTH activity proxy.
    pub fn background_rate(&self) -> f64 {
        if self.bg_covered == 0 {
            0.0
        } else {
            self.bg_edited as f64 / self.bg_covered as f64
        }
    }

    /// Observed conversion rate; `0.0` when the cell has no coverage (never NaN).
    pub fn rate(&self) -> f64 {
        if self.covered == 0 {
            0.0
        } else {
            self.edited as f64 / self.covered as f64
        }
    }
}

impl std::ops::AddAssign<&CellActivity> for CellActivity {
    fn add_assign(&mut self, other: &Self) {
        self.edited += other.edited;
        self.covered += other.covered;
        self.bg_edited += other.bg_edited;
        self.bg_covered += other.bg_covered;
    }
}

/// Per-cell tallies for one library, built by [`scan`].
pub type ActivityTally = FxHashMap<CellBarcode, CellActivity>;

/// Competent cells per signal BAM path (barcodes are library-local); `Arc` so the
/// per-gene rayon loop shares rather than clones the set.
pub type CompetentCells = FxHashMap<Box<str>, Arc<FxHashSet<CellBarcode>>>;

/// Control (catalytically-dead arm) cells as counts only; the null needs no identity.
/// A-to-I has no equivalent: endogenous ADAR is graded, with no null population.
pub type ControlCells = [CellActivity];

/// Coverage-weighted pooled rate, so a shallow cell cannot outvote a deep one.
pub fn pooled_rate<'a>(cells: impl Iterator<Item = &'a CellActivity>) -> f64 {
    let (mut e, mut n) = (0u64, 0u64);
    for c in cells {
        e += c.edited;
        n += c.covered;
    }
    if n == 0 {
        0.0
    } else {
        e as f64 / n as f64
    }
}

/// Closed-form method-of-moments fit of `y_i ~ BetaBinomial(n_i, mean, rho)`
/// (Kleinman 1973, JASA), using `E[S] = m(1−m)·[(k−1) + ρ·(N − Σnᵢ²/N)]` for
/// `S = Σ nᵢ(p̂ᵢ − m)²`. `ρ` matters: inflation `1+(n−1)ρ` is large for deep cells.
///
/// Returns `(mean, rho)`, ρ clamped to `[0, 1)`; `rho = 0` with < 3 covered cells.
pub fn fit_betabinom_mom(counts: &[(u64, u64)]) -> (f64, f64) {
    let live: Vec<(f64, f64)> = counts
        .iter()
        .filter(|(_, n)| *n > 0)
        .map(|(y, n)| (*y as f64, *n as f64))
        .collect();
    let k = live.len();
    let total_n: f64 = live.iter().map(|(_, n)| n).sum();
    if total_n <= 0.0 {
        return (0.0, 0.0);
    }
    let mean = live.iter().map(|(y, _)| y).sum::<f64>() / total_n;
    if k < 3 || mean <= 0.0 || mean >= 1.0 {
        return (mean, 0.0);
    }
    let s: f64 = live.iter().map(|(y, n)| n * (y / n - mean).powi(2)).sum();
    let sum_sq: f64 = live.iter().map(|(_, n)| n * n).sum();
    let denom = total_n - sum_sq / total_n;
    if denom <= 0.0 {
        return (mean, 0.0);
    }
    let rho = (s / (mean * (1.0 - mean)) - (k as f64 - 1.0)) / denom;
    (mean, rho.clamp(0.0, 0.99))
}

/// Equal-count strata labels `0..n_strata` over `values`; ties break by index.
pub fn quantile_strata(values: &[f64], n_strata: usize, min_per: usize) -> Vec<usize> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }
    let k = n_strata.min(n / min_per.max(1)).max(1);
    if k <= 1 {
        return vec![0; n];
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        values[a]
            .partial_cmp(&values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    let mut label = vec![0usize; n];
    for (rank, &i) in order.iter().enumerate() {
        label[i] = rank * k / n;
    }
    label
}

/// Allowed discarded-pool rate as a multiple of the control; `1.0` is the
/// parameter-free, conservative point (discarded pool equals control).
pub const DEFAULT_REJECT_TOLERANCE: f64 = 1.0;

/// Minimum coverage for a cell to be scored at all.
pub const DEFAULT_SCAN_MIN_COVERAGE: u64 = 50;

/// Knobs for [`call_competent_cells`].
#[derive(Clone, Copy, Debug)]
pub struct NullCallOpts {
    /// Minimum coverage for a cell to be scored; below it the cell is rejected unscored.
    pub min_coverage: u64,
    /// Number of coverage strata (equal-count).
    pub n_strata: usize,
    /// Minimum cells per stratum; fewer ⇒ fewer strata.
    pub min_per_stratum: usize,
    /// Max discarded-pool rate as a multiple of the control rate.
    pub reject_tolerance: f64,
    /// Never reject more than this fraction (guards against a pathological control).
    pub max_reject_frac: f64,
    /// If > 0, keep WT cells scoring above the `1 - control_tail` quantile of
    /// depth-matched control cells instead of using [`Self::reject_tolerance`].
    pub control_tail: f64,
}

impl Default for NullCallOpts {
    fn default() -> Self {
        Self {
            min_coverage: DEFAULT_SCAN_MIN_COVERAGE,
            n_strata: 12,
            min_per_stratum: 50,
            reject_tolerance: DEFAULT_REJECT_TOLERANCE,
            max_reject_frac: 0.95,
            control_tail: 0.0,
        }
    }
}

/// The outcome, including everything worth logging as QC.
#[derive(Clone, Debug)]
pub struct NullCellCall {
    /// Cells that edit; handed to discovery as `CellMembership`.
    pub selected: FxHashSet<CellBarcode>,
    pub n_scored: usize,
    /// Pooled conversion rate of the kept cells.
    pub selected_rate: f64,
    /// Pooled conversion rate of the discarded cells.
    pub rejected_rate: f64,
    /// Pooled conversion rate of the control arm.
    pub control_rate: f64,
    /// Percentile of depth-matched control cells the weakest kept cell exceeds.
    pub control_percentile: f64,
}

impl NullCellCall {
    pub fn n_selected(&self) -> usize {
        self.selected.len()
    }

    pub fn kept_frac(&self) -> f64 {
        if self.n_scored == 0 {
            0.0
        } else {
            self.n_selected() as f64 / self.n_scored as f64
        }
    }

    /// The QC number: ≈1 means the discarded cells look like the dead enzyme.
    pub fn rejected_over_control(&self) -> f64 {
        self.rejected_rate / self.control_rate
    }
}

/// Beta-binomial standardized deviate of a cell's rate from its stratum's control
/// mean (not a p-value, which would underflow and tie the top of the ranking).
fn stratum_score(a: &CellActivity, mean: f64, rho: f64) -> f64 {
    let n = a.covered as f64;
    if n <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let m = mean.clamp(1e-9, 1.0 - 1e-9);
    let var = m * (1.0 - m) / n * (1.0 + (n - 1.0).max(0.0) * rho);
    if var <= 0.0 {
        return f64::NEG_INFINITY;
    }
    (a.rate() - m) / var.sqrt()
}

/// Call competent WT cells: score each against its coverage stratum's control
/// null, rank, and take the most aggressive cut whose rejected pool is within
/// tolerance of the control.
pub fn call_competent_cells(
    wt: &ActivityTally,
    control: &ControlCells,
    opts: &NullCallOpts,
) -> NullCellCall {
    // Same cells, same order everywhere, for the `n_wt + i` offset below.
    let ctrl: Vec<&CellActivity> = control
        .iter()
        .filter(|c| c.covered >= opts.min_coverage)
        .collect();
    let control_rate = pooled_rate(ctrl.iter().copied());

    let mut cells: Vec<(&CellBarcode, &CellActivity)> = wt
        .iter()
        .filter(|(_, a)| a.covered >= opts.min_coverage)
        .collect();
    cells.sort_by(|a, b| a.0.cmp(b.0));
    if cells.is_empty() || control_rate <= 0.0 {
        // No null to calibrate against: refuse to cut.
        return NullCellCall {
            selected: wt.keys().cloned().collect(),
            n_scored: 0,
            control_percentile: f64::NAN,
            selected_rate: pooled_rate(wt.values()),
            rejected_rate: 0.0,
            control_rate,
        };
    }
    let n_wt = cells.len();

    let ctrl: Vec<&CellActivity> = control
        .iter()
        .filter(|c| c.covered >= opts.min_coverage)
        .collect();

    // One coverage scale for both arms; the null is fit on control only.
    let mut depths: Vec<f64> = cells.iter().map(|(_, a)| a.covered as f64).collect();
    depths.extend(ctrl.iter().map(|c| c.covered as f64));
    let strata = quantile_strata(&depths, opts.n_strata, opts.min_per_stratum);
    let n_strata = strata.iter().copied().max().unwrap_or(0) + 1;
    let mut by_stratum: Vec<Vec<(u64, u64)>> = vec![Vec::new(); n_strata];
    for (i, c) in ctrl.iter().enumerate() {
        by_stratum[strata[n_wt + i]].push((c.edited, c.covered));
    }
    let null_by_stratum: Vec<(f64, f64)> = by_stratum
        .iter()
        .map(|counts| {
            // Too few cells for a dispersion: use the pooled control rate.
            if counts.len() >= 3 {
                fit_betabinom_mom(counts)
            } else {
                (control_rate, 0.0)
            }
        })
        .collect();

    // Control scores, so the cut can be reported (or placed) on the control's scale.
    let control_scores: Vec<f64> = {
        let mut v: Vec<f64> = ctrl
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let k = strata.get(n_wt + i).copied().unwrap_or(0);
                let (m, rho) = null_by_stratum
                    .get(k)
                    .copied()
                    .unwrap_or((control_rate, 0.0));
                stratum_score(c, m, rho)
            })
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v
    };

    // Rank by the depth-adjusted deviate, best first; ties by barcode.
    let mut ranked: Vec<(f64, &CellBarcode, &CellActivity)> = cells
        .iter()
        .enumerate()
        .map(|(i, (cb, a))| {
            let (m, rho) = null_by_stratum[strata[i]];
            (stratum_score(a, m, rho), *cb, *a)
        })
        .collect();
    ranked.sort_by(|x, y| {
        y.0.partial_cmp(&x.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(x.1.cmp(y.1))
    });

    // Suffix sums: the rejected pool for a cut at `k` is ranks `k..n`.
    let n = ranked.len();
    let mut suf_e = vec![0u64; n + 1];
    let mut suf_n = vec![0u64; n + 1];
    for i in (0..n).rev() {
        suf_e[i] = suf_e[i + 1] + ranked[i].2.edited;
        suf_n[i] = suf_n[i + 1] + ranked[i].2.covered;
    }
    let ratio_at = |k: usize| -> f64 {
        if suf_n[k] == 0 {
            0.0
        } else {
            (suf_e[k] as f64 / suf_n[k] as f64) / control_rate
        }
    };

    // `min_keep` guards against rejecting everything.
    let min_keep = ((1.0 - opts.max_reject_frac) * n as f64).ceil() as usize;
    let mut cut = n;
    if opts.control_tail > 0.0 && !control_scores.is_empty() {
        let idx = ((1.0 - opts.control_tail) * control_scores.len() as f64).floor() as usize;
        let threshold = control_scores[idx.min(control_scores.len() - 1)];
        cut = ranked
            .partition_point(|(s, _, _)| *s > threshold)
            .max(min_keep.max(1));
    } else {
        for k in min_keep.max(1)..=n {
            if ratio_at(k) <= opts.reject_tolerance {
                cut = k;
                break;
            }
        }
    }

    // Where the chosen cut lands among depth-matched control cells (descriptive).
    let control_percentile = match ranked.get(cut.saturating_sub(1)) {
        Some((boundary, _, _)) if !control_scores.is_empty() => {
            let below = control_scores.partition_point(|v| v < boundary);
            100.0 * below as f64 / control_scores.len() as f64
        }
        _ => f64::NAN,
    };

    let selected: FxHashSet<CellBarcode> = ranked[..cut]
        .iter()
        .map(|(_, cb, _)| (*cb).clone())
        .collect();
    let sel_rate = pooled_rate(ranked[..cut].iter().map(|(_, _, a)| *a));
    let rej_rate = pooled_rate(ranked[cut..].iter().map(|(_, _, a)| *a));

    NullCellCall {
        selected,
        n_scored: n,
        control_percentile,
        selected_rate: sel_rate,
        rejected_rate: rej_rate,
        control_rate,
    }
}

#[cfg(test)]
mod tests;

/// Knobs shared by `faba dartseq` and `faba all`. No off switch: the scan no-ops
/// when there is nothing to calibrate against (see [`call_and_report`]).
#[derive(Args, Debug, Clone, serde::Serialize)]
pub struct CellScanArgs {
    /// Keep WT cells above this upper tail of depth-matched control cells (0: off)
    #[arg(long = "cell-scan-control-tail", default_value_t = 0.0)]
    pub cell_scan_control_tail: f64,

    /// Diagnostic only: m6A reader genes to summarise per cell, by symbol
    #[arg(long = "reader-genes", value_delimiter = ',')]
    pub reader_genes: Vec<String>,

    /// Diagnostic only: m6A writer/eraser genes to summarise per cell, by symbol
    #[arg(long = "writer-genes", value_delimiter = ',')]
    pub writer_genes: Vec<String>,

    /// Override the data-driven cut (multiple of the control rate)
    #[arg(
        long = "cell-scan-tolerance",
        help = "Override the data-driven cut (multiple of the control rate the discarded pool may reach)",
        long_help = "Override the data-driven cut.\n\
                     \n\
                     Unset (the default), the cut is decided from the data:\n\
                     cells are dropped up to the point where the discarded pool's conversion rate EQUALS the control's,\n\
                     so nothing demonstrably real is thrown away. Nothing is tuned.\n\
                     \n\
                     Set it to cut deeper and concentrate the kept pool.\n\
                     e.g. `1.5` permits the discarded pool to edit 50% above the control.\n\
                     Every run logs `dropped/control` and the equivalent control percentile,\n\
                     so the cost of raising it is visible. Raising it trades purity for recall;\n\
                     the default is the defensible choice, not the maximal one.\n\
                     \n\
                     `--cell-scan-control-tail` takes precedence over this if both are set."
    )]
    pub cell_scan_tolerance: Option<f64>,

    /// Minimum conversion-site coverage for a cell to be scored
    #[arg(long = "cell-scan-min-coverage", default_value_t = DEFAULT_SCAN_MIN_COVERAGE)]
    pub cell_scan_min_coverage: u64,

    /// Drop null cells from the OUTPUT matrices too, not just from site calling
    #[arg(
        long = "quantify-competent-only",
        help = "Restrict the m6A matrices to editing-competent cells as well"
    )]
    pub quantify_competent_only: bool,
}

impl CellScanArgs {
    fn opts(&self) -> NullCallOpts {
        NullCallOpts {
            min_coverage: self.cell_scan_min_coverage,
            reject_tolerance: self.cell_scan_tolerance.unwrap_or(DEFAULT_REJECT_TOLERANCE),
            control_tail: self.cell_scan_control_tail,
            ..NullCallOpts::default()
        }
    }
}

/// Tally, call and audit per-cell editing; return the competent set per signal
/// BAM, or `None` (keep every cell) when there is nothing to calibrate against.
pub fn call_and_report(
    gff_map: &GffRecordMap,
    params: &crate::editing::pipeline::ConversionParams,
    args: &CellScanArgs,
    gene_qc: Option<&crate::quant::GeneCountQc>,
    output_dir: &str,
    label: &str,
) -> anyhow::Result<Option<Arc<CompetentCells>>> {
    let qc_cells = gene_qc.map(|q| &q.cells_by_batch);
    let matrix_paths: Vec<Box<str>> = gene_qc
        .map(|q| q.matrix_by_batch.values().cloned().collect())
        .unwrap_or_default();
    let matrix_paths = matrix_paths.as_slice();
    // No control arm (A-to-I, or m6A without `--control-bam`): no null to calibrate.
    if params.mut_bam_files.is_empty() {
        info!("cell scan: no control arm; keeping every cell");
        return Ok(None);
    }

    // Fail early on a bad `--genome` rather than panic in a rayon worker.
    crate::data::util_htslib::load_fasta_index(&params.genome_file)?;

    // All control libraries pool into one null (pairing is not expressible in
    // general); per-library rates are logged so a hot control is visible.
    let mut control = Vec::new();
    let mut per_library = Vec::new();
    for bam in &params.mut_bam_files {
        let mut tally = scan::scan_cell_activity(gff_map, params, bam, "control")?;
        retain_qc_cells(&mut tally, qc_cells, bam);
        let rate = pooled_rate(
            tally
                .values()
                .filter(|c| c.covered >= args.cell_scan_min_coverage),
        );
        per_library.push((bam.clone(), rate, tally.len()));
        control.extend(tally.into_values());
    }
    if per_library.len() > 1 {
        for (bam, rate, n) in &per_library {
            info!(
                "control library {bam}: {n} cells, C->U {:.4}%",
                100.0 * rate
            );
        }
        let rates: Vec<f64> = per_library.iter().map(|(_, r, _)| *r).collect();
        let (lo, hi) = rates
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
        if lo > 0.0 && hi / lo > 1.25 {
            log::warn!(
                "control libraries differ {:.2}x in background editing ({:.4}% to {:.4}%); \
                 they are POOLED into one null, so the hotter library raises the bar for \
                 every signal library",
                hi / lo,
                100.0 * lo,
                100.0 * hi
            );
        }
    }
    if control.is_empty() {
        info!("cell scan: control tallied no cells; keeping every cell");
        return Ok(None);
    }
    let opts = args.opts();

    // One call per signal library; batch names from the same list as the matrices.
    let quant = params.quant_bam_files();
    let quant_batches = uniq_batch_names(&quant)?;
    let batch_of: FxHashMap<&str, &str> = quant
        .iter()
        .zip(quant_batches.iter())
        .map(|(f, b)| (f.as_ref(), b.as_ref()))
        .collect();
    let mut by_batch = CompetentCells::default();
    for bam in params.wt_bam_files.iter() {
        let batch = batch_of.get(bam.as_ref()).copied().unwrap_or(bam.as_ref());
        let mut signal = scan::scan_cell_activity(gff_map, params, bam, "signal")?;
        let observed = signal.len();
        retain_qc_cells(&mut signal, qc_cells, bam);
        if signal.is_empty() {
            // Left out of the map, this library keeps every cell; say so.
            log::warn!(
                "cell scan ({batch}): no cells tallied (missing {} tag?); this library \
                 stays UNDILUTED while the others are filtered",
                params.cell_barcode_tag
            );
            continue;
        }
        let scored = signal
            .values()
            .filter(|a| a.covered >= opts.min_coverage)
            .count();
        // The cell-QC cascade per stage, since each stage can fail independently.
        info!(
            "cell QC cascade ({batch}): {observed} barcodes observed -> {qc} passed \
             droplet/complexity/mito QC ({qc_pct:.1}%) -> {scored} had coverage >= {floor} \
             ({scored_pct:.1}% of QC) -> see the competence call below",
            observed = observed,
            qc = signal.len(),
            qc_pct = 100.0 * signal.len() as f64 / observed.max(1) as f64,
            scored = scored,
            floor = opts.min_coverage,
            scored_pct = 100.0 * scored as f64 / signal.len().max(1) as f64,
        );
        if log::log_enabled!(log::Level::Info) {
            log_rate_histogram(&signal, &control, &opts);
        }
        let call = call_competent_cells(&signal, &control, &opts);
        info!(
            "cell QC ({batch}): competent {}/{} scored cells ({:.1}%); editing {:.4}% \
             kept vs {:.4}% dropped vs {:.4}% control; dropped/control {:.2}; \
             cut sits at control p{:.0}",
            call.n_selected(),
            call.n_scored,
            100.0 * call.kept_frac(),
            100.0 * call.selected_rate,
            100.0 * call.rejected_rate,
            100.0 * call.control_rate,
            call.rejected_over_control(),
            call.control_percentile
        );
        let rule = if args.cell_scan_control_tail > 0.0 {
            format!(
                "control p{:.0}",
                100.0 * (1.0 - args.cell_scan_control_tail)
            )
        } else if let Some(t) = args.cell_scan_tolerance {
            format!("tolerance {t} (overriding the data-driven cut)")
        } else {
            "data-driven (discarded pool == control)".to_string()
        };
        info!("cell QC ({batch}): cut placed by {rule}");
        if call.rejected_over_control() > opts.reject_tolerance {
            info!(
                "cell QC ({batch}): dropped cells still edit {:.2}x the control -- the \
                 cut is limited by --cell-scan-tolerance, not by the data",
                call.rejected_over_control()
            );
        }

        let path = format!("{output_dir}/{batch}_{label}_cell_qc.tsv.gz");
        write_audit(&path, &signal, &call, opts.min_coverage)?;
        info!("cell QC audit -> {path}");

        // Diagnostic only; runs after the cut and never feeds it.
        for (lbl, genes) in [
            ("m6A readers", &args.reader_genes),
            ("m6A writers/erasers", &args.writer_genes),
        ] {
            if let Err(e) = log_family_expression(lbl, genes, matrix_paths, &signal, &call.selected)
            {
                log::warn!("{lbl}: could not summarise expression ({e})");
            }
        }

        if call.n_selected() == 0 {
            info!("cell QC ({batch}): selected no cells; keeping every cell instead");
            continue;
        }
        by_batch.insert(bam.clone(), Arc::new(call.selected));
    }

    Ok((!by_batch.is_empty()).then(|| Arc::new(by_batch)))
}

/// Cells for the m6A matrices: `Some` only with `--quantify-competent-only` and a
/// competence call; otherwise the caller uses the full QC set.
pub fn cells_for_quantification(
    competent: Option<&Arc<CompetentCells>>,
    qc_cells: Option<&FxHashMap<Box<str>, FxHashSet<CellBarcode>>>,
    competent_only: bool,
) -> Option<FxHashMap<Box<str>, FxHashSet<CellBarcode>>> {
    match (competent, qc_cells) {
        (Some(comp), Some(qc)) if competent_only => Some(intersect_for_quantification(qc, comp)),
        _ => None,
    }
}

/// Intersect each signal BAM's QC cells with its competent cells; BAMs without a
/// call (including controls) keep every QC cell.
pub fn intersect_for_quantification(
    qc_cells: &FxHashMap<Box<str>, FxHashSet<CellBarcode>>,
    competent: &CompetentCells,
) -> FxHashMap<Box<str>, FxHashSet<CellBarcode>> {
    qc_cells
        .iter()
        .map(|(bam, cells)| match competent.get(bam) {
            Some(keep) => {
                let kept: FxHashSet<CellBarcode> = cells
                    .iter()
                    .filter(|c| keep.contains(*c))
                    .cloned()
                    .collect();
                info!(
                    "quantification ({bam}): restricted to {} of {} QC cells",
                    kept.len(),
                    cells.len()
                );
                (bam.clone(), kept)
            }
            // A control BAM (or a signal BAM with no call) keeps every QC cell.
            None => (bam.clone(), cells.clone()),
        })
        .collect()
}

/// Keep only cells that passed upstream cell QC, so ambient droplets are never
/// scored. A BAM with no QC entry is left unfiltered.
fn retain_qc_cells(
    tally: &mut ActivityTally,
    qc_cells: Option<&FxHashMap<Box<str>, FxHashSet<CellBarcode>>>,
    bam_file: &str,
) {
    let Some(keep) = qc_cells.and_then(|m| m.get(bam_file)) else {
        return;
    };
    tally.retain(|cb, _| keep.contains(cb));
}

/// ASCII histogram of per-cell rate, signal vs control, on a logit axis (rates
/// are small, so a linear axis would pile them into the leftmost bins).
fn log_rate_histogram(signal: &ActivityTally, control: &ControlCells, opts: &NullCallOpts) {
    const BINS: usize = 24;
    const WIDTH: usize = 28;
    let logit = |p: f64| -> f64 {
        let p = p.clamp(1e-6, 1.0 - 1e-6);
        (p / (1.0 - p)).ln()
    };
    let grab = |it: &mut dyn Iterator<Item = &CellActivity>| -> Vec<f64> {
        it.filter(|c| c.covered >= opts.min_coverage && c.edited > 0)
            .map(|c| logit(c.rate()))
            .collect()
    };
    let sig = grab(&mut signal.values());
    let ctl = grab(&mut control.iter());
    info!("per-cell C->U rate, logit axis (signal | control), cells with >=1 conversion:");
    info!("{:>10}  {:>28} | {:<28}", "rate", "SIGNAL", "CONTROL");
    two_series_histogram(&sig, &ctl, BINS, WIDTH, |centre| {
        format!("{:>9.4}%", 100.0 / (1.0 + (-centre).exp()))
    });
}

/// Does this feature row belong to `symbol`? Accepts raw `{id}_{symbol}` and
/// canonicalised `{symbol}` names; field-anchored, so `GENE1` never matches `GENE1L`.
fn feature_is_gene(feature: &str, symbol: &str) -> bool {
    data_beans::aux::feature_rows::parse_feature_row(feature)
        .is_some_and(|r| r.gene == symbol || r.gene.rsplit('_').next() == Some(symbol))
}

/// Diagnostic only: per-10k-UMI expression of a gene family, competent vs null
/// cells. Nothing here feeds the cut.
fn log_family_expression(
    label: &str,
    symbols: &[String],
    matrix_paths: &[Box<str>],
    scanned: &ActivityTally,
    selected: &FxHashSet<CellBarcode>,
) -> anyhow::Result<()> {
    use data_beans::aux::data_loading::{read_data_on_shared_rows, ReadSharedRowsArgs};
    use data_beans::qc::collect_column_stat_across_vec;
    use data_beans::sparse_io_vector::ColumnAlignment;
    if symbols.is_empty() {
        return Ok(());
    }
    // `--valid-cells` skips the QC stage that builds these matrices.
    if matrix_paths.is_empty() {
        log::warn!(
            "{label}: {symbols:?} requested, but there is no gene-count matrix to read \
             (--valid-cells skips the gene-count QC stage that builds it); skipping"
        );
        return Ok(());
    }
    let loaded = read_data_on_shared_rows(ReadSharedRowsArgs {
        data_files: matrix_paths.to_vec(),
        column_alignment: ColumnAlignment::Disjoint,
        // The null side is every unselected barcode; gating would change the contrast.
        keep_empty_barcodes: true,
        ..Default::default()
    })?;
    let backend = &loaded.data;
    let features = backend.row_names()?;
    let rows: Vec<usize> = features
        .iter()
        .enumerate()
        .filter(|(_, f)| symbols.iter().any(|g| feature_is_gene(f, g)))
        .map(|(i, _)| i)
        .collect();
    if rows.is_empty() {
        log::warn!("{label}: none of {symbols:?} matched a feature row; skipping");
        return Ok(());
    }
    let (_, fam, _, _) = collect_column_stat_across_vec(backend, Some(&rows), None)?.to_f32_vecs();
    let (_, tot, _, _) = collect_column_stat_across_vec(backend, None, None)?.to_f32_vecs();
    let barcodes = backend.column_names()?;
    let mut kept: Vec<f32> = Vec::new();
    let mut dropped: Vec<f32> = Vec::new();
    for (i, cb) in barcodes.iter().enumerate() {
        // Strip the Disjoint loader's `@{batch}` suffix; the scan keys on bare barcodes.
        let bare = cb.split_once('@').map_or(cb.as_ref(), |(b, _)| b);
        let key = CellBarcode::Barcode(bare.into());
        // Scanned library only; control cells would contrast conditions instead.
        if !scanned.contains_key(&key) {
            continue;
        }
        let t = *tot.get(i).unwrap_or(&0.0);
        if t <= 0.0 {
            continue;
        }
        let v = *fam.get(i).unwrap_or(&0.0) / t * 1e4;
        if selected.contains(&key) {
            kept.push(v)
        } else {
            dropped.push(v)
        }
    }
    // An empty side means matrix columns and scan barcodes name different cells.
    if kept.is_empty() || dropped.is_empty() {
        log::warn!(
            "{label}: cannot compare — {} competent / {} null cells matched the matrix. \
             Matrix column e.g. {:?}; scanned barcode e.g. {:?}",
            kept.len(),
            dropped.len(),
            barcodes.first().map(|b| b.as_ref()),
            scanned.keys().next()
        );
        return Ok(());
    }
    let (mk, md) = (
        legume_numeric::matrix::utils::median(&kept),
        legume_numeric::matrix::utils::median(&dropped),
    );
    info!(
        "{label} ({} genes, {} rows): per-10k median {:.1} in competent cells vs {:.1} in \
         null cells{}",
        symbols.len(),
        rows.len(),
        mk,
        md,
        if (mk - md).abs() <= 0.15 * mk.max(1e-9) {
            "  -- indistinguishable, so competence is not explained by this family"
        } else {
            ""
        }
    );
    side_by_side(&kept, &dropped, "competent", "null");
    Ok(())
}

/// Two distributions on one shared log axis, per-10k units; zeros are counted
/// separately rather than clamped onto the axis.
fn side_by_side(a: &[f32], b: &[f32], la: &str, lb: &str) {
    let split = |v: &[f32]| -> (usize, Vec<f64>) {
        let pos: Vec<f64> = v
            .iter()
            .filter(|x| **x > 0.0)
            .map(|x| (*x as f64).ln())
            .collect();
        (v.len() - pos.len(), pos)
    };
    let ((za, x), (zb, y)) = (split(a), split(b));
    info!(
        "  zero expression: {}/{} {la} vs {}/{} {lb}",
        za,
        a.len(),
        zb,
        b.len()
    );
    info!("{:>10}  {:>22} | {:<22}", "per-10k", la, lb);
    two_series_histogram(&x, &y, 16, 22, |c| format!("{:>10.1}", c.exp()));
}

/// Two pre-transformed series binned on one shared axis, side by side; `axis`
/// renders a bin centre in display units.
fn two_series_histogram(
    a: &[f64],
    b: &[f64],
    bins: usize,
    width: usize,
    axis: impl Fn(f64) -> String,
) {
    if a.is_empty() || b.is_empty() {
        return;
    }
    let (lo, hi) = a
        .iter()
        .chain(b.iter())
        .fold((f64::MAX, f64::MIN), |(x, y), &v| (x.min(v), y.max(v)));
    // Also catches NaN.
    if hi <= lo || !hi.is_finite() || !lo.is_finite() {
        return;
    }
    let bin = |v: f64| -> usize { (((v - lo) / (hi - lo) * bins as f64) as usize).min(bins - 1) };
    let (mut ha, mut hb) = (vec![0usize; bins], vec![0usize; bins]);
    for v in a {
        ha[bin(*v)] += 1;
    }
    for v in b {
        hb[bin(*v)] += 1;
    }
    let peak = ha
        .iter()
        .chain(hb.iter())
        .copied()
        .max()
        .unwrap_or(1)
        .max(1);
    let bar = |n: usize| {
        let k = (n * width).div_ceil(peak).min(width);
        format!("{}{}", "#".repeat(k), " ".repeat(width - k))
    };
    for i in (0..bins).rev() {
        let centre = lo + (i as f64 + 0.5) * (hi - lo) / bins as f64;
        info!(
            "{}  {} | {}  {:>5} {:>5}",
            axis(centre),
            bar(ha[i]),
            bar(hb[i]),
            ha[i],
            hb[i]
        );
    }
}

/// Per-cell audit: every scored cell, its tally, and whether it was kept.
fn write_audit(
    path: &str,
    signal: &ActivityTally,
    call: &NullCellCall,
    min_coverage: u64,
) -> anyhow::Result<()> {
    let mut rows: Vec<_> = signal.iter().collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    // `scored` separates "rejected as null" from "too little coverage to assess".
    let mut lines: Vec<Box<str>> = Vec::with_capacity(rows.len() + 1);
    lines.push(
        "barcode\tedited\tcovered\trate\tbg_edited\tbg_covered\tbg_rate\tscored\tkept".into(),
    );
    for (cb, act) in rows {
        lines.push(
            format!(
                "{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{}\t{}",
                cb,
                act.edited,
                act.covered,
                act.rate(),
                act.bg_edited,
                act.bg_covered,
                act.background_rate(),
                u8::from(act.covered >= min_coverage),
                u8::from(call.selected.contains(cb))
            )
            .into(),
        );
    }
    // `write_lines` propagates write errors (a `GzEncoder` dropped would swallow them).
    write_lines(&lines, path)
}
