//! `faba run` command-line surface.
//!
//! One struct for the whole pipeline: each step reads the subset it needs and
//! builds the standalone subcommand's own args from it, so a chained run and a
//! hand-run subcommand cannot drift apart.

use crate::common::*;

#[derive(Args, Debug, serde::Serialize)]
pub struct PipelineArgs {
    // Inputs: required only under --batch-process (see `check_batch`)
    #[arg(
        value_delimiter = ',',
        help = "Input BAM files (comma-separated)",
        long_help = "Comma-separated BAM files used across every modality: gene counting, ATOI,\n\
                     APA and m6A quantification."
    )]
    pub bam_files: Vec<Box<str>>,

    #[arg(short = 'g', long = "gff", help = "Gene annotation (GFF) file")]
    pub gff_file: Option<Box<str>>,

    #[arg(
        short = 'f',
        long = "genome",
        help = "Reference genome FASTA file (.fa/.fasta, must be indexed)"
    )]
    pub genome_file: Option<Box<str>>,

    #[arg(
        short = 'o',
        long = "output",
        help = "Output directory (flat structure)"
    )]
    pub output: Option<Box<str>>,

    #[arg(
        long = "control-bam",
        alias = "mut",
        alias = "control",
        value_delimiter = ',',
        help = "Control BAM files (catalytically-dead YTHmut) for the m6A contrast",
        long_help = "Comma-separated control (catalytically-dead YTHmut) BAM files.\n\
                     m6A is called by a WT-vs-MUT contrast:\n\
                     the signal arm is the positional BAMs MINUS these controls.\n\
                     That split is used only for m6A site discovery;\n\
                     otherwise these controls are quantified like positional samples.\n\
                     Their SNP, gene, ATOI, APA and m6A per-cell matrices are produced too,\n\
                     with cells frozen per control BAM in step 1. Optional,\n\
                     but the m6A (DART) step is skipped without it:\n\
                     m6A cannot be separated from genomic C/T variation without a control."
    )]
    pub control_bam_files: Vec<Box<str>>,

    #[arg(
        long = "batch-process",
        default_value_t = false,
        help = "Run straight through with the flags given, without the setup view",
        long_help = "Run straight through with the flags given, without the setup view.\n\
                     The BAMs, -g, -f and -o are then required. Without it, `faba run`\n\
                     opens a view to pick the inputs, steps and flags, which needs stdin\n\
                     and stdout on a terminal."
    )]
    pub batch_process: bool,

    ///////////////////////
    // Shared parameters //
    ///////////////////////
    #[arg(
        help_heading = "Common",
        long,
        default_value = "CB",
        help = "Cell barcode tag"
    )]
    pub cell_barcode_tag: Box<str>,

    #[arg(
        help_heading = "Common",
        long,
        default_value = "GX",
        help = "Gene barcode tag"
    )]
    pub gene_barcode_tag: Box<str>,

    #[arg(
        help_heading = "Common",
        long,
        value_enum,
        default_value = "zarr",
        help = "Sparse matrix backend (zarr or hdf5)"
    )]
    // `SparseIoBackend` is data-beans' and does not implement `Serialize`; it is a plain
    // enum, so its `Debug` form ("Zarr" / "Hdf5") is exactly what belongs in the summary.
    #[serde(serialize_with = "crate::run_record::ser_debug")]
    pub backend: SparseIoBackend,

    #[arg(
        help_heading = "Common",
        long = "no-zip",
        default_value_t = true,
        action = clap::ArgAction::SetFalse,
        help = "Keep a `.zarr` directory instead of producing a `.zarr.zip` archive",
        long_help = "Keep a `.zarr` directory instead of producing a `.zarr.zip` archive.\n\
                     Zarr backend only; no effect on hdf5."
    )]
    pub zip: bool,

    #[arg(
        help_heading = "Common",
        long,
        alias = "threads",
        default_value_t = 16,
        help = "Maximum number of threads"
    )]
    pub max_threads: usize,

    ///////////////////////////////
    // Gene expression filtering //
    ///////////////////////////////
    #[arg(
        help_heading = "count",
        long,
        default_value_t = 1,
        help = "Minimum cells per gene; 1 = drop only empty rows",
        long_help = "Minimum cells per gene, for gene filtering.\n\
                     The default of 1 drops only genes with no counts at all;\n\
                     an opinionated floor belongs to `faba qc --row-nnz-cutoff`.\n\
                     Matches the standalone subcommands."
    )]
    pub gene_min_cells: usize,

    #[arg(
        help_heading = "count",
        long,
        default_value_t = 0,
        help = "Minimum UMI per gene; 0 = off",
        long_help = "Minimum UMI per gene, for gene filtering. 0 turns it off,\n\
                     matching Cell Ranger, which keeps every gene and feature."
    )]
    pub gene_min_counts: usize,

    #[arg(
        help_heading = "count",
        long,
        default_value_t = 1,
        help = "Minimum detected genes (nnz) per cell; 1 = drop only empty columns",
        long_help = "Minimum detected genes, as nnz, per cell.\n\
                     Cells below it are dropped from gene counts,\n\
                     and from every downstream modality.\n\
                     \n\
                     The default of 1 drops only cells with no counts at all,\n\
                     so cell calling is pure Cell Ranger EmptyDrops/OrdMag;\n\
                     an opinionated floor belongs to `faba qc`."
    )]
    pub cell_min_genes: usize,

    #[command(flatten, next_help_heading = "count")]
    pub cell_qc: crate::cell_qc::CellQcArgs,

    //////////////////////////////////////////
    // Gene biotype (quantification subset) //
    //////////////////////////////////////////
    #[arg(
        help_heading = "count",
        long,
        default_value = "",
        help = "Gene biotype to quantify; empty keeps all",
        long_help = "Gene biotype to quantify. Empty is the default, and keeps all biotypes.\n\
                     Pass a value to restrict: protein_coding, lncRNA, pseudogene.\n\
                     \n\
                     QC and cell-calling always use ALL biotypes.\n\
                     Only the quantified gene set is restricted.\n\
                     That set is gene counts plus ATOI, APA and m6A."
    )]
    pub gene_type: Box<str>,

    //////////////////////
    // Mitochondrial QC //
    //////////////////////
    #[command(flatten, next_help_heading = "count")]
    pub mito_qc: crate::quant::MitoQcArgs,

    ////////////////////////////////////////////////////
    // Shared read-quality filters (ATOI / m6A / SNP) //
    ////////////////////////////////////////////////////
    #[arg(
        help_heading = "Common",
        long,
        default_value_t = 20,
        help = "Minimum base quality for editing/SNP base calls (ATOI/m6A/SNP)"
    )]
    pub min_base_quality: u8,

    #[arg(
        help_heading = "Common",
        long,
        default_value_t = 20,
        help = "Minimum mapping quality (MAPQ) for every read the pipeline admits",
        long_help = "Minimum mapping quality (MAPQ) for every read the pipeline admits.\n\
                     Applies to the gene counts that freeze the cell set,\n\
                     and to the editing/SNP pileups (ATOI/m6A/SNP).\n\
                     One knob on purpose: cells are then called on the same\n\
                     alignments the site tests later count."
    )]
    pub min_mapping_quality: u8,

    #[arg(
        help_heading = "APA",
        long = "no-apa-pdui",
        default_value_t = false,
        help = "Skip the APA PDUI (proximal/distal count) matrix output ({batch}_apa)"
    )]
    pub no_apa_pdui: bool,

    /////////////////////
    // ATOI parameters //
    /////////////////////
    // Shared with `faba atoi` by const, so the same BAM cannot produce a
    // different A-to-I site list depending on which command ran it.
    #[arg(
        help_heading = "ATOI",
        long,
        default_value_t = crate::editing::pipeline::DEFAULT_ATOI_MIN_COVERAGE,
        help = "Minimum coverage (ref + alt) for an ATOI site to be written (matches `faba atoi`)"
    )]
    pub atoi_min_coverage: usize,

    #[arg(
        help_heading = "ATOI",
        long,
        default_value_t = crate::editing::pipeline::DEFAULT_ATOI_MIN_CONVERSION,
        help = "Minimum A-to-G (alt) reads for an ATOI site to be written (matches `faba atoi`)"
    )]
    pub atoi_min_conversion: usize,

    ///////////////////////////////////////////////////////
    // Editing statistical null (ATOI only; m6A is a contrast) //
    ///////////////////////////////////////////////////////
    #[arg(
        help_heading = "ATOI",
        long = "edit-error-rate",
        alias = "error-rate",
        default_value_t = 0.01,
        help = "Sequencing-error rate ε: the beta-binomial null mean",
        long_help = "Sequencing-error rate ε. It is the beta-binomial null mean.\n\
                     The edited fraction is tested against it. The test is reference-anchored,\n\
                     with no control sample."
    )]
    pub edit_error_rate: f64,

    #[arg(
        help_heading = "ATOI",
        long = "edit-overdispersion",
        alias = "overdispersion",
        default_value_t = 0.1,
        help = "Beta-binomial overdispersion ρ for the editing null (0 ⇒ binomial)"
    )]
    pub edit_overdispersion: f64,

    ////////////////////
    // APA parameters //
    ////////////////////
    #[arg(
        help_heading = "APA",
        long,
        default_value_t = 10,
        help = "Minimum coverage for APA detection"
    )]
    pub apa_min_coverage: usize,

    #[arg(
        help_heading = "APA",
        long = "apa-max-sites",
        default_value_t = 20,
        help = "Cap candidate poly-A sites per UTR; 0 = unlimited",
        long_help = "Cap on candidate poly-A sites per UTR, for APA BIC selection.\n\
                     Sites are ranked top-N by coverage; 0 is unlimited.\n\
                     This bounds EM cost on long 3'UTRs."
    )]
    pub apa_max_sites: usize,

    #[arg(
        help_heading = "APA",
        long = "apa-em-pdui",
        default_value_t = false,
        help = "Use the full SCAPE EM for PDUI",
        long_help = "Use the full SCAPE EM for PDUI.\n\
                     The default is a fast top-2 nearest-site assignment. The EM is slower.\n\
                     --mixture also forces it."
    )]
    pub apa_em_pdui: bool,

    #[arg(
        help_heading = "APA",
        long,
        default_value_t = 10,
        help = "Minimum poly(A) tail length"
    )]
    pub polya_min_tail_length: usize,

    /////////////////////
    // DART parameters //
    /////////////////////
    // Shared with `faba dartseq` by const, not by convention — see
    // `DEFAULT_M6A_MIN_COVERAGE` for why these two must not be allowed to drift,
    // and the long_help on `faba dartseq --min-coverage` for the measured cost of
    // the current values.
    #[arg(
        help_heading = "m6A",
        long,
        default_value_t = crate::editing::pipeline::DEFAULT_M6A_MIN_COVERAGE,
        help = "Minimum total reads (signal + control) for an m6A site to be written (matches `faba dartseq`)"
    )]
    pub m6a_min_coverage: usize,

    #[arg(
        help_heading = "m6A",
        long,
        default_value_t = crate::editing::pipeline::DEFAULT_M6A_MIN_CONVERSION,
        help = "Minimum converted (C->T) signal reads for an m6A site to be written (matches `faba dartseq`)"
    )]
    pub m6a_min_conversion: usize,

    #[command(flatten, next_help_heading = "m6A")]
    pub cell_scan: crate::editing::cell_activity::CellScanArgs,

    ////////////////////////////////////////////////////////
    // Mixture model weighting (shared by m6A and A-to-I) //
    ////////////////////////////////////////////////////////
    #[arg(
        help_heading = "m6A",
        long = "mixture-weight",
        value_enum,
        default_value_t = crate::editing::pipeline::MixtureWeightMode::Posterior,
        help = "How to weight each (cell, site) observation in the mixture EM"
    )]
    pub mixture_weight: crate::editing::pipeline::MixtureWeightMode,

    // α = β = 1 (uniform Beta) is what `faba dartseq` and `faba atoi` use. The
    // pipeline used 1e-4, a near-improper prior that leaves a 1-of-1 site at
    // almost full weight — the very thing posterior weighting exists to damp.
    // The help string already said "(default: 1.0)" while the code said 1e-4,
    // so the restated default is dropped here: clap prints the real one.
    #[arg(
        help_heading = "m6A",
        long = "mixture-prior-alpha",
        default_value_t = 1.0,
        help = "Beta prior α for posterior-rate weighting"
    )]
    pub mixture_prior_alpha: f32,

    #[arg(
        help_heading = "m6A",
        long = "mixture-prior-beta",
        default_value_t = 1.0,
        help = "Beta prior β for posterior-rate weighting"
    )]
    pub mixture_prior_beta: f32,

    #[arg(
        help_heading = "m6A",
        long = "drop-single-component",
        default_value_t = false,
        help = "Drop genes with a single mixture component across m6A/ATOI/APA"
    )]
    pub drop_single_component: bool,

    #[arg(
        help_heading = "m6A",
        long = "mixture",
        default_value_t = false,
        help = "Also produce the per-gene component-mixture matrices",
        long_help = "Also produce the per-gene component-mixture matrices.\n\
                     They come from an EM fit, which is slow.\n\
                     \n\
                     This is off by default.\n\
                     Only the gene-level and per-site matrices are produced then.\n\
                     \n\
                     For m6A / A-to-I this SKIPS the 1-D Gaussian mixture EM entirely when off.\n\
                     For APA it also runs the SCAPE poly-A EM and writes `_apa_mixture`;\n\
                     when off, PDUI takes proximal and distal sites from a fast split of read\n\
                     positions, with no EM."
    )]
    pub mixture: bool,

    ////////////////////
    // SNP parameters //
    ////////////////////
    #[arg(
        long = "known-snps",
        help = "Known SNP sites VCF/BCF/Parquet to force-call",
        long_help = "Path to known SNP sites. Accepts:\n\
                     - VCF/BCF (.vcf, .vcf.gz, .bcf): standard variant calls\n\
                     - Parquet (.parquet): output from a previous `faba snp` run\n\
                     When provided,\n\
                     force-calls genotypes at these positions in addition to de novo discovery."
    )]
    pub known_snps: Option<Box<str>>,

    #[arg(
        help_heading = "SNP",
        long,
        default_value_t = 5,
        help = "Minimum depth for SNP calling"
    )]
    pub snp_min_depth: usize,

    #[arg(
        help_heading = "SNP",
        long,
        default_value_t = 20.0,
        help = "Minimum genotype quality (Phred) to emit a call"
    )]
    pub snp_min_gq: f32,

    #[arg(
        help_heading = "SNP",
        long,
        default_value_t = 10,
        help = "Minimum coverage for de novo SNP discovery"
    )]
    pub snp_min_coverage: usize,

    #[arg(
        help_heading = "SNP",
        long,
        default_value_t = 3,
        help = "Minimum alt allele reads for SNP discovery"
    )]
    pub snp_min_alt_count: usize,

    #[arg(
        help_heading = "SNP",
        long,
        default_value_t = 0.1,
        help = "Minimum alt allele frequency for SNP discovery"
    )]
    pub snp_min_alt_freq: f64,

    //////////////////////////////////////////////
    // UMI deduplication (applies to all steps) //
    //////////////////////////////////////////////
    #[arg(
        help_heading = "Common",
        long = "umi-tag",
        default_value = "UB",
        help = "UMI barcode BAM tag for deduplication (all steps)"
    )]
    pub umi_tag: Box<str>,

    #[arg(
        help_heading = "Common",
        long = "no-umi-dedup",
        default_value_t = false,
        help = "Disable UMI deduplication (for bulk data without UMIs)"
    )]
    pub no_umi_dedup: bool,

    //////////////////
    // Step control //
    //////////////////
    #[arg(long, default_value_t = false, help = "Skip SNP genotyping step")]
    pub skip_snp: bool,

    #[arg(
        long = "skip-count",
        alias = "skip-genes",
        default_value_t = false,
        help = "Skip the gene counting / cell calling step"
    )]
    pub skip_count: bool,

    #[arg(long, default_value_t = false, help = "Skip ATOI detection step")]
    pub skip_atoi: bool,

    #[arg(long, default_value_t = false, help = "Skip APA quantification step")]
    pub skip_apa: bool,

    #[arg(long, default_value_t = false, help = "Skip m6A detection step")]
    pub skip_m6a: bool,

    #[arg(
        long = "depth-resolution-kb",
        help = "Also bin per-cell read depth at this resolution, in KILOBASES",
        long_help = "Bin per-cell read depth genome-wide at this resolution, in kilobases.\n\
                     \n\
                     Omit the flag and the depth step is skipped entirely.\n\
                     It is opt-in because it costs a full extra pass over every BAM.\n\
                     \n\
                     Depth reuses the cells frozen by gene counting,\n\
                     so its columns line up with the other modalities.\n\
                     Rows are named {chr}:{start}-{end}, the one modality not keyed by gene,\n\
                     because a bin is not a gene.\n\
                     \n\
                     The grid is anchored at each chromosome's start and does not depend on how the work was chunked.\n\
                     Only the last bin of a chromosome is short,\n\
                     as in the standard CNV binning tools."
    )]
    pub depth_resolution_kb: Option<f32>,
}

impl PipelineArgs {
    /// Under `--batch-process`, every input the run needs, or what is missing.
    pub fn check_batch(&self) -> anyhow::Result<()> {
        let mut missing = Vec::new();
        if self.bam_files.is_empty() {
            missing.push("BAM files");
        }
        if self.gff_file.is_none() {
            missing.push("-g/--gff");
        }
        if self.genome_file.is_none() {
            missing.push("-f/--genome");
        }
        if self.output.is_none() {
            missing.push("-o/--output");
        }
        anyhow::ensure!(
            missing.is_empty(),
            "`faba run --batch-process` needs {}",
            missing.join(", ")
        );
        Ok(())
    }

    /// The annotation; call after [`Self::check_batch`].
    pub fn gff(&self) -> &str {
        self.gff_file.as_deref().unwrap_or_default()
    }

    /// The genome; call after [`Self::check_batch`].
    pub fn genome(&self) -> &str {
        self.genome_file.as_deref().unwrap_or_default()
    }

    /// The output directory; call after [`Self::check_batch`].
    pub fn out(&self) -> &str {
        self.output.as_deref().unwrap_or_default()
    }
}
