//! Shared quantification layer for the BAM subcommands: output paths, gene keys,
//! UMI tags, cell/gene/mito QC gates, channel rows and matrix sinks.

use crate::common::*;
use crate::data::conversion::*;
use crate::data::util_htslib::*;
use crate::gene_count::splice::{count_read_per_gene_splice, format_gene_key, CountReadOpts};

use dashmap::DashMap as HashMap;
use data_beans::zarr_io::finalize_zarr_output;
use genomic_data::gff::GffRecordMap;
use genomic_data::sam::CellBarcode;
use rustc_hash::FxHashMap;

/// Backend paths: `write_path` is written, `target_path` is final (`.zarr.zip` if zipping).
pub struct BackendOutputPath {
    pub write_path: Box<str>,
    pub target_path: Box<str>,
}

impl BackendOutputPath {
    pub fn new(output_dir: &str, name: &str, backend: &SparseIoBackend, zip: bool) -> Self {
        match backend {
            SparseIoBackend::HDF5 => {
                let p: Box<str> = format!("{}/{}.h5", output_dir, name).into_boxed_str();
                Self {
                    write_path: p.clone(),
                    target_path: p,
                }
            }
            SparseIoBackend::Zarr if zip => Self {
                write_path: format!("{}/{}.zarr", output_dir, name).into_boxed_str(),
                target_path: format!("{}/{}.zarr.zip", output_dir, name).into_boxed_str(),
            },
            SparseIoBackend::Zarr => {
                let p: Box<str> = format!("{}/{}.zarr", output_dir, name).into_boxed_str();
                Self {
                    write_path: p.clone(),
                    target_path: p,
                }
            }
        }
    }

    /// Zip the staging `.zarr` into its `.zarr.zip` target (no-op otherwise).
    pub fn finalize(&self) -> anyhow::Result<()> {
        finalize_zarr_output(&self.write_path, &self.target_path)
    }
}

/// Check BAM indices for all files
pub fn check_all_bam_indices(bam_files: &[Box<str>]) -> anyhow::Result<()> {
    for bam_file in bam_files {
        info!("checking .bai file for {}...", bam_file);
        check_bam_index(bam_file, None)?;
    }
    Ok(())
}

/// Gene key `{gene_id}_{symbol}`, the first part of `{gene_key}/{modality}/{detail}`.
pub fn create_gene_key_function(
    gff_map: &GffRecordMap,
) -> impl Fn(&BedWithGene) -> Box<str> + Send + Sync + '_ {
    |x: &BedWithGene| -> Box<str> {
        gff_map
            .get(&x.gene)
            .map(|gff| format!("{}_{}", gff.gene_id, gff.gene_name))
            .unwrap_or_else(|| format!("{}", x.gene))
            .into_boxed_str()
    }
}

/// Push a channel-last `feature_row` of `count`, skipping zeros; `subunit` is a site or component.
pub fn push_channel_row(
    triplets: &mut Vec<(CellBarcode, Box<str>, f32)>,
    cb: &CellBarcode,
    gene: &str,
    modality: &str,
    channel: &str,
    subunit: Option<&str>,
    count: usize,
) {
    if count > 0 {
        triplets.push((
            cb.clone(),
            data_beans::aux::feature_rows::feature_row(gene, modality, channel, subunit),
            count as f32,
        ));
    }
}

/// Pool conversion stats per (cell, gene) into `pos_channel`/`neg_channel` rows, skipping zeros.
pub fn summarize_stats_two_channel<F>(
    stats: &[(CellBarcode, BedWithGene, ConversionData)],
    gene_key_func: F,
    modality: &str,
    pos_channel: &str,
    neg_channel: &str,
) -> TripletsRowsCols
where
    F: Fn(&BedWithGene) -> Box<str> + Send + Sync,
{
    let combined: HashMap<(CellBarcode, Box<str>), ConversionData> = HashMap::default();
    stats.par_iter().for_each(|(cb, bed, dat)| {
        let key = (cb.clone(), gene_key_func(bed));
        combined.entry(key).or_default().add_assign(dat);
    });

    let mut triplets: Vec<(CellBarcode, Box<str>, f32)> = Vec::with_capacity(combined.len() * 2);
    for ((cb, gene), dat) in combined {
        push_channel_row(
            &mut triplets,
            &cb,
            &gene,
            modality,
            pos_channel,
            None,
            dat.converted,
        );
        push_channel_row(
            &mut triplets,
            &cb,
            &gene,
            modality,
            neg_channel,
            None,
            dat.unconverted,
        );
    }
    format_data_triplets(triplets)
}

/// Per-site form of [`summarize_stats_two_channel`]; the subunit is `{chr}:{start}` (0-based).
pub fn summarize_stats_per_site<F>(
    stats: &[(CellBarcode, BedWithGene, ConversionData)],
    gene_key_func: F,
    modality: &str,
    pos_channel: &str,
    neg_channel: &str,
) -> TripletsRowsCols
where
    F: Fn(&BedWithGene) -> Box<str> + Send + Sync,
{
    let combined: HashMap<(CellBarcode, Box<str>, Box<str>), ConversionData> = HashMap::default();
    stats.par_iter().for_each(|(cb, bed, dat)| {
        let gene = gene_key_func(bed);
        let site: Box<str> = format!("{}:{}", bed.chr, bed.start).into();
        let key = (cb.clone(), gene, site);
        combined.entry(key).or_default().add_assign(dat);
    });

    let entries: Vec<_> = combined.into_iter().collect();

    let mut triplets: Vec<(CellBarcode, Box<str>, f32)> = Vec::with_capacity(entries.len() * 2);
    for ((cb, gene, site), dat) in &entries {
        push_channel_row(
            &mut triplets,
            cb,
            gene,
            modality,
            pos_channel,
            Some(site.as_ref()),
            dat.converted,
        );
        push_channel_row(
            &mut triplets,
            cb,
            gene,
            modality,
            neg_channel,
            Some(site.as_ref()),
            dat.unconverted,
        );
    }
    format_data_triplets(triplets)
}

/// Gene-count QC: genes pooled across batches, cells kept per batch (barcodes are per library).
pub struct GeneCountQc {
    pub gene_ids: rustc_hash::FxHashSet<GeneId>,
    /// Valid cells keyed by BAM path (basenames collide across libraries).
    pub cells_by_batch: rustc_hash::FxHashMap<Box<str>, rustc_hash::FxHashSet<CellBarcode>>,
    /// Persisted gene-count matrix per BAM path; empty on the `--valid-cells` path.
    pub matrix_by_batch: rustc_hash::FxHashMap<Box<str>, Box<str>>,
}

/// Optional matrix write for [`run_gene_count_qc`]; QC artifacts are always written.
pub struct GeneMatrixSink<'a> {
    pub backend: &'a SparseIoBackend,
    pub zip: bool,
}

/// Gene key of a feature name: `GENE1/count/spliced` gives `GENE1`.
#[inline]
pub fn extract_gene_key(feat: &str) -> &str {
    feat.rfind("/count/")
        .map(|pos| &feat[..pos])
        .unwrap_or(feat)
}

/// One batch's per-gene `(nnz, total)` (spliced + unspliced) and spliced-only passing cells.
pub struct BatchQc {
    pub gene_stats: FxHashMap<Box<str>, (usize, f64)>,
    pub passing_cells: rustc_hash::FxHashSet<CellBarcode>,
}

/// Compute [`BatchQc`]: gene stats include unspliced; cell calling is spliced-only (Cell Ranger).
pub fn batch_qc(
    spliced: &[(CellBarcode, Box<str>, f32)],
    unspliced: &[(CellBarcode, Box<str>, f32)],
    want_total: bool,
    cell_min_genes: usize,
    cell_call: &crate::cell_qc::CellCallParams,
) -> BatchQc {
    let spliced_pairs: rustc_hash::FxHashSet<(&CellBarcode, &str)> = spliced
        .par_iter()
        .map(|(cb, feat, _)| (cb, extract_gene_key(feat)))
        .collect();

    let mut cell_nnz: FxHashMap<&CellBarcode, usize> = FxHashMap::default();
    for &(cb, _) in &spliced_pairs {
        *cell_nnz.entry(cb).or_default() += 1;
    }

    let mut gene_cell_pairs = spliced_pairs;
    gene_cell_pairs.extend(
        unspliced
            .iter()
            .map(|(cb, feat, _)| (cb, extract_gene_key(feat))),
    );

    let mut gene_nnz: FxHashMap<&str, usize> = FxHashMap::default();
    for &(_, gk) in &gene_cell_pairs {
        *gene_nnz.entry(gk).or_default() += 1;
    }

    let gene_total: FxHashMap<&str, f64> = if want_total {
        let mut m: FxHashMap<&str, f64> = FxHashMap::default();
        for (_, feat, v) in spliced.iter().chain(unspliced.iter()) {
            *m.entry(extract_gene_key(feat)).or_default() += *v as f64;
        }
        m
    } else {
        FxHashMap::default()
    };

    let gene_stats: FxHashMap<Box<str>, (usize, f64)> = gene_nnz
        .into_iter()
        .map(|(gk, nnz)| {
            let total = if want_total {
                gene_total.get(gk).copied().unwrap_or(0.0)
            } else {
                0.0
            };
            (Box::from(gk), (nnz, total))
        })
        .collect();

    // `cell_min_genes` always applies; the `Nnz` policy keeps that raw superset.
    let nnz_cells: rustc_hash::FxHashSet<CellBarcode> = cell_nnz
        .into_iter()
        .filter(|(_, n)| *n >= cell_min_genes)
        .map(|(cb, _)| cb.clone())
        .collect();

    let passing_cells: rustc_hash::FxHashSet<CellBarcode> =
        if cell_call.filter == crate::cell_qc::CellFilter::Nnz {
            nnz_cells
        } else {
            let counts = crate::cell_qc::CellCounts::from_triplets(spliced, &[]);
            crate::cell_qc::call_cells(&counts, cell_call)
                .into_iter()
                .filter(|cb| nnz_cells.contains(cb))
                .collect()
        };

    BatchQc {
        gene_stats,
        passing_cells,
    }
}

/// Gene thresholds on `(nnz, total)` stats; `gene_min_counts == 0` disables the count floor.
pub fn passing_genes_from_stats(
    gene_stats: &FxHashMap<Box<str>, (usize, f64)>,
    gene_min_cells: usize,
    gene_min_counts: usize,
) -> rustc_hash::FxHashSet<Box<str>> {
    let min_counts = gene_min_counts as f64;
    gene_stats
        .iter()
        .filter(|(_, (nnz, _))| *nnz >= gene_min_cells)
        .filter(|(_, (_, total))| gene_min_counts == 0 || *total >= min_counts)
        .map(|(gk, _)| gk.clone())
        .collect()
}

/// Add one batch's gene stats into `pooled`; exact, since each cell is in one library.
pub fn accumulate_gene_stats(
    pooled: &mut FxHashMap<Box<str>, (usize, f64)>,
    batch_stats: FxHashMap<Box<str>, (usize, f64)>,
) {
    for (gene_key, (nnz, total)) in batch_stats {
        let e = pooled.entry(gene_key).or_default();
        e.0 += nnz;
        e.1 += total;
    }
}

/// Splice-aware gene counting and QC over every batch; the one gene-counting loop.
pub fn run_gene_count_qc(gff_file: &str, req: &GeneQcRequest) -> anyhow::Result<GeneCountQc> {
    info!("=== Gene expression QC ===");

    let bam_files = req.bam_files;
    let persist = req.persist.as_ref();
    let batch_names = uniq_batch_names(bam_files)?;

    let all_records = read_gff_record_vec(gff_file)?;
    let gene_map = build_gene_map(&all_records, Some(&FeatureType::Gene))?;
    let exon_map = build_exon_intervals(&all_records);
    let exon_intervals: FxHashMap<GeneId, Vec<(i64, i64)>> = exon_map.into_iter().collect();

    let gff_map = GffRecordMap::from_map(gene_map);
    info!("Loaded {} genes for expression QC", gff_map.len());

    let records = gff_map.records();

    let gene_key_to_id: FxHashMap<Box<str>, GeneId> = records
        .iter()
        .map(|rec| (format_gene_key(rec), rec.gene_id.clone()))
        .collect();

    let gate = GeneGate::new(&records, req.gene_type, req.mito.clone());

    let mut expressed_gene_ids: rustc_hash::FxHashSet<GeneId> = rustc_hash::FxHashSet::default();
    let mut cells_by_batch: rustc_hash::FxHashMap<Box<str>, rustc_hash::FxHashSet<CellBarcode>> =
        rustc_hash::FxHashMap::default();
    let mut matrix_by_batch: rustc_hash::FxHashMap<Box<str>, Box<str>> =
        rustc_hash::FxHashMap::default();
    // Pooled, then thresholded once, so the gene set is partition-invariant.
    let mut pooled_gene_stats: FxHashMap<Box<str>, (usize, f64)> = FxHashMap::default();

    for (bam_file, batch_name) in bam_files.iter().zip(batch_names.iter()) {
        let njobs = records.len() as u64;
        info!(
            "Counting genes (splice-aware) in {} ({} genes)",
            bam_file, njobs
        );

        let results: Vec<_> = records
            .par_iter()
            .progress_with(new_progress_bar(njobs))
            .map_init(crate::data::bam_io::BamReaderCache::new, |cache, rec| {
                count_read_per_gene_splice(cache, bam_file, rec, &exon_intervals, req.count)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        let mut spliced_triplets = Vec::new();
        let mut unspliced_triplets = Vec::new();
        for r in results {
            spliced_triplets.extend(r.spliced);
            unspliced_triplets.extend(r.unspliced);
        }

        info!(
            "{} spliced, {} unspliced triplets",
            spliced_triplets.len(),
            unspliced_triplets.len()
        );

        let bq = qc_one_batch(
            &spliced_triplets,
            &unspliced_triplets,
            &gate,
            req.gene_min_counts,
            req.cell_min_genes,
            &req.cell_call,
            Some(QcArtifacts {
                dir: req.output_dir,
                batch_name,
            }),
        )?;

        // Count QC sees every gene; the gate then narrows to the biotype and drops mito.
        let passing_genes: rustc_hash::FxHashSet<Box<str>> =
            passing_genes_from_stats(&bq.gene_stats, req.gene_min_cells, req.gene_min_counts)
                .into_iter()
                .filter(|gk| gate.quantify(gk))
                .collect();
        accumulate_gene_stats(&mut pooled_gene_stats, bq.gene_stats);

        info!(
            "{}: {} genes, {} cells passed QC",
            batch_name,
            passing_genes.len(),
            bq.passing_cells.len()
        );

        if let Some(sink) = persist {
            let keep = |t: Vec<(CellBarcode, Box<str>, f32)>| -> Vec<_> {
                t.into_par_iter()
                    .filter(|(cb, feat, _)| {
                        passing_genes.contains(extract_gene_key(feat))
                            && bq.passing_cells.contains(cb)
                    })
                    .collect::<Vec<_>>()
            };
            let spliced = keep(spliced_triplets);
            let unspliced = keep(unspliced_triplets);

            let UnionNames {
                col_names,
                cell_to_index,
                row_names,
                feature_to_index,
            } = collect_union_names(&spliced, &unspliced);

            let out = BackendOutputPath::new(
                req.output_dir,
                &format!("{}_count", batch_name),
                sink.backend,
                sink.zip,
            );
            let merged: Vec<_> = spliced.into_iter().chain(unspliced).collect();
            format_data_triplets_shared(
                merged,
                &feature_to_index,
                &cell_to_index,
                row_names,
                col_names,
            )
            .to_backend(&out.write_path)?;
            out.finalize()?;
            info!(
                "{}: wrote spliced + unspliced to {}",
                batch_name, out.target_path
            );
            matrix_by_batch.insert(bam_file.clone(), out.target_path);
        }

        cells_by_batch.insert(bam_file.clone(), bq.passing_cells);
    }

    let passing_genes =
        passing_genes_from_stats(&pooled_gene_stats, req.gene_min_cells, req.gene_min_counts);
    for gene_key in passing_genes.iter().filter(|gk| gate.quantify(gk)) {
        if let Some(gene_id) = gene_key_to_id.get(gene_key) {
            expressed_gene_ids.insert(gene_id.clone());
        }
    }
    write_qc_genes(req.output_dir, &expressed_gene_ids)?;

    let total_cells: usize = cells_by_batch.values().map(|s| s.len()).sum();
    info!(
        "Gene QC summary: {} genes, {} cells passed across {} BAM files",
        expressed_gene_ids.len(),
        total_cells,
        bam_files.len()
    );

    Ok(GeneCountQc {
        gene_ids: expressed_gene_ids,
        cells_by_batch,
        matrix_by_batch,
    })
}

/// Resolve `--no-umi-dedup` / `--umi-tag` to the dedup tag (`None` = no dedup).
pub fn resolve_umi_tag(no_umi_dedup: bool, umi_tag: &str) -> Option<&[u8]> {
    if no_umi_dedup {
        None
    } else {
        Some(umi_tag.as_bytes())
    }
}

/// Inputs for [`resolve_gene_qc`], built by each modality runner from its own args.
pub struct GeneQcRequest<'a> {
    pub bam_files: &'a [Box<str>],
    /// BAM tags and read admission, matching what the modality's pileup admits.
    pub count: CountReadOpts<'a>,
    /// GFF for the recompute path; `None` skips recompute (reuse can still run).
    pub gff_file: Option<&'a str>,
    /// QC artifact directory, always written (the matrix is optional, see `persist`).
    pub output_dir: &'a str,
    /// Biotype to quantify (`""` = all); never affects cell calling.
    pub gene_type: &'a str,
    pub gene_min_cells: usize,
    pub gene_min_counts: usize,
    pub cell_min_genes: usize,
    pub cell_call: crate::cell_qc::CellCallParams,
    /// Mito cell filter and MT-gene exclusion policy.
    pub mito: MitoQcParams,
    pub valid_cells_file: Option<&'a str>,
    pub valid_genes_file: Option<&'a str>,
    pub skip_gene_qc: bool,
    /// Persist the QC gene-count matrix; `None` skips the write.
    pub persist: Option<GeneMatrixSink<'a>>,
}

/// Reuse `--valid-cells` or recompute QC; an empty `gene_ids` means keep all genes.
pub fn resolve_gene_qc(req: &GeneQcRequest) -> anyhow::Result<Option<GeneCountQc>> {
    if let Some(dir) = req.valid_cells_file {
        let cells_by_batch = load_valid_cells_dir(dir, req.bam_files)?;
        let gene_ids = match req.valid_genes_file {
            Some(gf) => load_valid_genes(gf)?,
            None => rustc_hash::FxHashSet::default(),
        };
        Ok(Some(GeneCountQc {
            gene_ids,
            cells_by_batch,
            matrix_by_batch: rustc_hash::FxHashMap::default(),
        }))
    } else if !req.skip_gene_qc {
        match req.gff_file {
            Some(gff_file) => Ok(Some(run_gene_count_qc(gff_file, req)?)),
            None => Ok(None),
        }
    } else {
        Ok(None)
    }
}

/// [`resolve_gene_qc`], then retain `gff_map` to the passing genes when a gene filter exists.
pub fn resolve_modality_gene_qc(
    gff_map: &mut GffRecordMap,
    req: &GeneQcRequest,
) -> anyhow::Result<Option<GeneCountQc>> {
    let qc = resolve_gene_qc(req)?;
    if let Some(ref qc) = qc {
        if !qc.gene_ids.is_empty() {
            gff_map.retain_by_ids(&qc.gene_ids);
            info!("After gene QC: {} genes retained", gff_map.len());
        }
    }
    Ok(qc)
}

/// Write retained barcodes to `{dir}/{batch}_cells.tsv.gz`, one per line.
pub fn write_qc_cells(
    dir: &str,
    batch_name: &str,
    cells: &rustc_hash::FxHashSet<CellBarcode>,
) -> anyhow::Result<()> {
    let mut lines: Vec<Box<str>> = cells
        .iter()
        .map(|c| c.to_string().into_boxed_str())
        .collect();
    lines.sort();
    let path = format!("{}/{}_cells.tsv.gz", dir, batch_name);
    write_lines(&lines, &path)?;
    info!("wrote {} retained cells to {}", lines.len(), path);
    Ok(())
}

/// Write retained gene ids to `{dir}/genes_kept.tsv.gz` (pooled, not per batch).
pub fn write_qc_genes(dir: &str, gene_ids: &rustc_hash::FxHashSet<GeneId>) -> anyhow::Result<()> {
    let mut lines: Vec<Box<str>> = gene_ids
        .iter()
        .map(|g| g.to_string().into_boxed_str())
        .collect();
    lines.sort();
    let path = format!("{}/genes_kept.tsv.gz", dir);
    write_lines(&lines, &path)?;
    info!("wrote {} retained genes to {}", lines.len(), path);
    Ok(())
}

/// Default mitochondrial chromosomes (clap default and [`MitoQcParams::default`]).
pub const MITO_CHR_DEFAULT: &str = "chrM,chrMT,MT,M";

/// Mitochondrial QC flags shared by every subcommand that does gene-count QC.
#[derive(clap::Args, Debug, Clone, serde::Serialize)]
pub struct MitoQcArgs {
    #[arg(
        long = "mito-chr",
        default_value = MITO_CHR_DEFAULT,
        help = "Mitochondrial chromosome name(s) (comma-separated)",
        long_help = "Genes on these chromosomes are treated as mitochondrial:\n\
                     excluded from the quantified gene set (unless --keep-mito) and summarized in the per-cell MT-fraction QC.\n\
                     Matched case-insensitively against the GFF seqname."
    )]
    pub mito_chr: Box<str>,

    #[arg(
        long = "keep-mito",
        default_value_t = false,
        help = "Keep mitochondrial genes in the quantified gene set",
        long_help = "By default mitochondrial genes are dropped from what is quantified.\n\
                     Their per-cell MT fraction is still reported as QC.\n\
                     Use this flag to retain them."
    )]
    pub keep_mito: bool,

    #[arg(
        long = "max-mito-frac",
        default_value_t = 0.0,
        help = "Max MT fraction per cell: >0 = fixed cutoff; 0 = elbow cutoff",
        long_help = "Cells whose mitochondrial UMI fraction exceeds the cutoff are removed during QC.\n\
                     A value > 0 is a fixed cutoff;\n\
                     the default 0 uses a data-driven elbow cutoff on the MT% distribution (drops the high-MT burst tail).\n\
                     See --no-mito-cell-qc to disable."
    )]
    pub max_mito_frac: f64,

    #[arg(
        long = "no-mito-cell-qc",
        default_value_t = false,
        help = "Disable MT cell QC (report MT% only, drop no cells)",
        long_help = "Report per-cell MT% but drop no cells.\n\
                     Mitochondrial genes are still excluded from the quantified gene set unless --keep-mito."
    )]
    pub no_mito_cell_qc: bool,
}

impl Default for MitoQcArgs {
    fn default() -> Self {
        let d = MitoQcParams::default();
        Self {
            mito_chr: d.mito_chr,
            keep_mito: d.keep_mito,
            max_mito_frac: d.max_mito_frac,
            no_mito_cell_qc: d.no_mito_cell_qc,
        }
    }
}

impl MitoQcArgs {
    /// Resolve to [`MitoQcParams`].
    pub fn params(&self) -> MitoQcParams {
        MitoQcParams {
            mito_chr: self.mito_chr.clone(),
            keep_mito: self.keep_mito,
            max_mito_frac: self.max_mito_frac,
            no_mito_cell_qc: self.no_mito_cell_qc,
        }
    }
}

/// Resolved, clap-free form of [`MitoQcArgs`].
#[derive(Debug, Clone)]
pub struct MitoQcParams {
    pub mito_chr: Box<str>,
    pub keep_mito: bool,
    /// Fixed per-cell MT-fraction cutoff; 0 uses the data-driven elbow.
    pub max_mito_frac: f64,
    pub no_mito_cell_qc: bool,
}

impl Default for MitoQcParams {
    fn default() -> Self {
        Self {
            mito_chr: MITO_CHR_DEFAULT.into(),
            keep_mito: false,
            max_mito_frac: 0.0,
            no_mito_cell_qc: false,
        }
    }
}

/// Genes a run quantifies (not mito-excluded, in the biotype); cell calling never uses it.
pub struct GeneGate {
    mito: MitoQcParams,
    mito_keys: rustc_hash::FxHashSet<Box<str>>,
    selected: Option<rustc_hash::FxHashSet<Box<str>>>,
}

impl GeneGate {
    /// Resolve against the annotation; an empty `gene_type` keeps all biotypes.
    pub fn new(records: &[GffRecord], gene_type: &str, mito: MitoQcParams) -> Self {
        let mito_keys = mito_gene_keys(records, &mito.mito_chr);
        info!(
            "{} mitochondrial gene(s) on {} ({})",
            mito_keys.len(),
            mito.mito_chr,
            if mito.keep_mito {
                "kept in matrix"
            } else {
                "excluded from matrix"
            }
        );

        let selected = if gene_type.is_empty() {
            info!("Gene biotype filter: OFF (all biotypes quantified)");
            None
        } else {
            let target: genomic_data::gff::GeneType = Box::<str>::from(gene_type).into();
            let keys: rustc_hash::FxHashSet<Box<str>> = records
                .iter()
                .filter(|r| r.gene_type == target)
                .map(format_gene_key)
                .collect();
            info!(
                "Gene biotype filter: quantifying {} '{}' genes (QC keeps all biotypes)",
                keys.len(),
                gene_type
            );
            Some(keys)
        };

        Self::from_keys(mito, mito_keys, selected)
    }

    /// Build from an already-resolved MT key set.
    pub fn from_keys(
        mito: MitoQcParams,
        mito_keys: rustc_hash::FxHashSet<Box<str>>,
        selected: Option<rustc_hash::FxHashSet<Box<str>>>,
    ) -> Self {
        Self {
            mito,
            mito_keys,
            selected,
        }
    }

    /// The mitochondrial QC policy.
    pub fn mito(&self) -> &MitoQcParams {
        &self.mito
    }

    /// `gene_key`s on mitochondrial chromosomes (the MT-fraction numerator).
    pub fn mito_keys(&self) -> &rustc_hash::FxHashSet<Box<str>> {
        &self.mito_keys
    }

    /// Whether a gene survives to quantification.
    pub fn quantify(&self, gene_key: &str) -> bool {
        (self.mito.keep_mito || !self.mito_keys.contains(gene_key))
            && self.selected.as_ref().is_none_or(|s| s.contains(gene_key))
    }
}

/// `gene_key`s on the comma-separated `mito_chr_spec` chromosomes (case-insensitive).
pub fn mito_gene_keys(
    records: &[GffRecord],
    mito_chr_spec: &str,
) -> rustc_hash::FxHashSet<Box<str>> {
    let mito: rustc_hash::FxHashSet<String> = mito_chr_spec
        .split(',')
        .map(|s| s.trim().to_ascii_lowercase())
        .filter(|s| !s.is_empty())
        .collect();
    records
        .iter()
        .filter(|r| mito.contains(&r.seqname.as_ref().to_ascii_lowercase()))
        .map(format_gene_key)
        .collect()
}

/// Per-cell `(mito, total)` UMI over passing cells. Serial on purpose: a parallel `f32`
/// reduce is order-dependent, so a cell on the cutoff could flip between runs.
pub fn mito_cell_stats(
    triplet_sets: &[&[(CellBarcode, Box<str>, f32)]],
    passing_cells: &rustc_hash::FxHashSet<CellBarcode>,
    mito_keys: &rustc_hash::FxHashSet<Box<str>>,
) -> FxHashMap<CellBarcode, (f32, f32)> {
    let has_mito = !mito_keys.is_empty();
    let mut stats: FxHashMap<CellBarcode, (f32, f32)> = FxHashMap::default();
    for set in triplet_sets {
        for (cb, feat, val) in set.iter() {
            if !passing_cells.contains(cb) {
                continue;
            }
            // `get_mut` first to avoid cloning the `Arc` barcode per triplet.
            let e = match stats.get_mut(cb) {
                Some(e) => e,
                None => stats.entry(cb.clone()).or_insert((0.0, 0.0)),
            };
            e.1 += *val; // total
            if has_mito && mito_keys.contains(extract_gene_key(feat)) {
                e.0 += *val; // mito
            }
        }
    }
    stats
}

/// Write `{dir}/{batch}_mt_qc.tsv.gz` (`barcode total_umi mt_umi mt_frac`, sorted).
pub fn write_mt_qc(
    dir: &str,
    batch_name: &str,
    stats: &FxHashMap<CellBarcode, (f32, f32)>,
) -> anyhow::Result<()> {
    let mut rows: Vec<(&CellBarcode, f32, f32)> =
        stats.iter().map(|(cb, (m, t))| (cb, *m, *t)).collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    let mut lines: Vec<Box<str>> = Vec::with_capacity(rows.len() + 1);
    lines.push("barcode\ttotal_umi\tmt_umi\tmt_frac".into());
    for (cb, mt, tot) in rows {
        let frac = if tot > 0.0 { mt / tot } else { 0.0 };
        lines.push(format!("{}\t{}\t{}\t{:.6}", cb, tot, mt, frac).into());
    }
    let path = format!("{}/{}_mt_qc.tsv.gz", dir, batch_name);
    write_lines(&lines, &path)?;
    info!("wrote MT QC ({} cells) to {}", lines.len() - 1, path);
    Ok(())
}

/// Elbow on ascending MT fractions (rank farthest from the end-to-end chord). `None` if too
/// few cells, flat, or the elbow is in the lower half (never cut a majority).
pub fn mito_elbow_cutoff(sorted_fracs: &[f64]) -> Option<f64> {
    let n = sorted_fracs.len();
    if n < 50 {
        return None;
    }
    let (ymin, ymax) = (sorted_fracs[0], sorted_fracs[n - 1]);
    let span = ymax - ymin;
    if span <= 1e-9 {
        return None; // flat distribution (e.g. no mito genes)
    }
    let xn = (n - 1) as f64;
    let (mut best_i, mut best_d) = (0usize, f64::NEG_INFINITY);
    for (i, &f) in sorted_fracs.iter().enumerate() {
        let x = i as f64 / xn;
        let y = (f - ymin) / span;
        let d = (x - y).abs();
        if d > best_d {
            best_d = d;
            best_i = i;
        }
    }
    if best_i < n / 2 {
        return None;
    }
    Some(sorted_fracs[best_i])
}

/// Drop high-MT cells: none if `disable`, else `max_frac` if > 0, else [`mito_elbow_cutoff`].
pub fn apply_mito_filter(
    passing_cells: rustc_hash::FxHashSet<CellBarcode>,
    stats: &FxHashMap<CellBarcode, (f32, f32)>,
    max_frac: f64,
    disable: bool,
) -> rustc_hash::FxHashSet<CellBarcode> {
    let frac_of = |cb: &CellBarcode| -> f64 {
        match stats.get(cb) {
            Some((m, t)) if *t > 0.0 => *m as f64 / *t as f64,
            _ => 0.0,
        }
    };
    let mut fracs: Vec<f64> = passing_cells.iter().map(&frac_of).collect();
    if fracs.is_empty() {
        return passing_cells;
    }
    fracs.sort_by(|a, b| a.total_cmp(b));
    let median = fracs[fracs.len() / 2];
    let mean = fracs.iter().sum::<f64>() / fracs.len() as f64;
    info!(
        "MT fraction over {} cells: median {:.3}, mean {:.3}, max {:.3}",
        fracs.len(),
        median,
        mean,
        fracs[fracs.len() - 1]
    );
    if disable {
        return passing_cells; // report-only
    }
    let (cutoff, kind) = if max_frac > 0.0 {
        (Some(max_frac), "fixed")
    } else {
        (mito_elbow_cutoff(&fracs), "elbow")
    };
    let Some(cutoff) = cutoff else {
        info!("MT QC: no clear high-MT burst population; no cells dropped");
        return passing_cells;
    };
    let kept: rustc_hash::FxHashSet<CellBarcode> = passing_cells
        .into_iter()
        .filter(|cb| frac_of(cb) <= cutoff)
        .collect();
    info!(
        "MT QC: dropped {} cells with MT fraction > {:.3} ({})",
        fracs.len() - kept.len(),
        cutoff,
        kind
    );
    kept
}

/// Destination for a batch's `{batch}_mt_qc.tsv.gz` and `{batch}_cells.tsv.gz`.
pub struct QcArtifacts<'a> {
    pub dir: &'a str,
    pub batch_name: &'a str,
}

/// QC one batch (cell call, MT metric, mito filter, artifacts); the one place passing
/// cells are decided. Cell calling is mito-blind, and the MT fraction uses the full
/// pre-gate counts (gating genes first would move the cutoff).
pub fn qc_one_batch(
    spliced: &[(CellBarcode, Box<str>, f32)],
    unspliced: &[(CellBarcode, Box<str>, f32)],
    gate: &GeneGate,
    gene_min_counts: usize,
    cell_min_genes: usize,
    cell_call: &crate::cell_qc::CellCallParams,
    out: Option<QcArtifacts>,
) -> anyhow::Result<BatchQc> {
    let bq = batch_qc(
        spliced,
        unspliced,
        gene_min_counts > 0,
        cell_min_genes,
        cell_call,
    );

    let mt_stats = mito_cell_stats(&[spliced, unspliced], &bq.passing_cells, gate.mito_keys());
    let mito = gate.mito();
    if let Some(ref out) = out {
        write_mt_qc(out.dir, out.batch_name, &mt_stats)?;
    }
    let passing_cells = apply_mito_filter(
        bq.passing_cells,
        &mt_stats,
        mito.max_mito_frac,
        mito.no_mito_cell_qc,
    );
    if let Some(ref out) = out {
        write_qc_cells(out.dir, out.batch_name, &passing_cells)?;
    }

    Ok(BatchQc {
        gene_stats: bq.gene_stats,
        passing_cells,
    })
}

/// Load per-batch `{batch}_cells.tsv.gz` from `dir`; a missing file leaves the batch unfiltered.
pub fn load_valid_cells_dir(
    dir: &str,
    bam_files: &[Box<str>],
) -> anyhow::Result<rustc_hash::FxHashMap<Box<str>, rustc_hash::FxHashSet<CellBarcode>>> {
    let batch_names = uniq_batch_names(bam_files)?;
    let mut out: rustc_hash::FxHashMap<Box<str>, rustc_hash::FxHashSet<CellBarcode>> =
        rustc_hash::FxHashMap::default();
    for (bam_file, batch) in bam_files.iter().zip(batch_names.iter()) {
        let path = format!("{}/{}_cells.tsv.gz", dir, batch);
        if !std::path::Path::new(&path).exists() {
            log::warn!(
                "--valid-cells: no file for batch '{}' ({}); not filtered",
                batch,
                path
            );
            continue;
        }
        let cells: rustc_hash::FxHashSet<CellBarcode> = read_lines(&path)?
            .into_iter()
            .filter(|s| s.as_ref() != ".")
            .map(|s| CellBarcode::Barcode(std::sync::Arc::from(s.as_ref())))
            .collect();
        info!(
            "--valid-cells: loaded {} cells for batch '{}'",
            cells.len(),
            batch
        );
        out.insert(bam_file.clone(), cells);
    }
    Ok(out)
}

/// Load `genes_kept.tsv.gz` (one gene id per line, shared across batches).
pub fn load_valid_genes(path: &str) -> anyhow::Result<rustc_hash::FxHashSet<GeneId>> {
    let genes: rustc_hash::FxHashSet<GeneId> = read_lines(path)?
        .into_iter()
        .filter(|s| s.as_ref() != ".")
        .map(|s| GeneId::Ensembl(s.as_ref().into()))
        .collect();
    info!("--valid-genes: loaded {} genes from {}", genes.len(), path);
    Ok(genes)
}

#[cfg(test)]
mod tests;
