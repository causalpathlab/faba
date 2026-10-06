use crate::data::cell_membership::CellMembership;
use crate::site_analysis::miami::bin::BinEdges;
use crate::site_analysis::miami::depth::read_depth_binned;
use crate::site_analysis::miami::genemodel::{
    load_gene_models, load_gene_models_where, models_extent, GeneModel,
};
use crate::site_analysis::miami::palette::Palette;
use crate::site_analysis::miami::render::{render_miami, FigOpts, PanelData};
use arrow::array::{Array, Float32Array, Int64Array, StringArray, UInt64Array};
use clap::Args;
use data_beans::aux::feature_names::FeatureNameKind;
use data_beans::aux::feature_rows::{channels, parse_feature_row, ATOI, METHYLATED, UNMETHYLATED};
use data_beans::hdf5_io::resolve_backend_file;
use data_beans::sparse_io::open_sparse_matrix;
use genomic_data::bed::Bed;
use genomic_data::coordinates::chr_eq;
use genomic_data::sam::CellBarcode;
use log::info;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::sync::Arc;

mod ascii;
mod figure;
mod read;
mod select;

pub(crate) use ascii::*;
use figure::*;
pub(crate) use read::*;
pub(crate) use select::*;

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum PileupSignal {
    /// Sum of values across all cells at each position
    Sum,
    /// Number of non-zero cells at each position
    Nnz,
    /// log10(1 + sum) at each position
    Log10Sum,
}

impl PileupSignal {
    fn name(&self) -> &'static str {
        match self {
            PileupSignal::Sum => "sum",
            PileupSignal::Nnz => "nnz",
            PileupSignal::Log10Sum => "log10-sum",
        }
    }
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum SiteSignal {
    /// Number of sites per bin
    Count,
    /// Sum of wild-type base coverage per bin
    WtCoverage,
    /// Sum of mutant base coverage per bin
    MutCoverage,
    /// Sum of -log10(p-value) per bin
    NegLog10Pv,
}

impl SiteSignal {
    fn name(&self) -> &'static str {
        match self {
            SiteSignal::Count => "count",
            SiteSignal::WtCoverage => "wt-coverage",
            SiteSignal::MutCoverage => "mut-coverage",
            SiteSignal::NegLog10Pv => "-log10(pv)",
        }
    }
}

/// Figure output format when only one is wanted. Omit to get the default
/// SVG + PDF pair.
#[derive(clap::ValueEnum, Clone, Debug)]
pub enum FigFormat {
    Svg,
    Pdf,
}

#[derive(Args, Debug, Clone)]
pub struct PileupArgs {
    #[arg(
        num_args = 1..,
        value_name = "MATRIX|DIR",
        help = "Site matrices (zarr or h5) from faba output, or the output directory (or a run record in it); browse for one when left out",
        long_help = "Site matrices (zarr or h5) from faba output: `{batch}_m6a_site` or\n\
                     `{batch}_atoi_site`. Multiple files, e.g. replicates via a shell glob,\n\
                     are aggregated per genomic position into a single track.\n\
                     \n\
                     Or a faba output directory (from `faba run`, a producer, or `faba qc`),\n\
                     or a run record in it (`*.run.json`, `pipeline_summary.json`): its site\n\
                     matrices of --modality are taken, one track per batch, with its\n\
                     `{modality}_sites.parquet` for -s and its recorded GFF for the genes.\n\
                     \n\
                     Left out, a browser asks: mark site matrices to pile up together,\n\
                     in one track, or choose an output directory."
    )]
    pub data_files: Vec<Box<str>>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi",
        help = "Which sites to take from an output directory (default: m6a, else atoi)"
    )]
    modality: Option<Box<str>>,

    #[arg(
        short = 'q',
        long = "genes",
        visible_alias = "gene",
        value_delimiter = ',',
        help = "Genes to pile up: comma-separated symbols (`GENE1,GENE2`) or Ensembl IDs",
        long_help = "Genes to pile up: comma-separated symbols (`GENE1,GENE2`) or Ensembl IDs,\n\
                     case-insensitive. Uses the shared relaxed gene-name scheme;\n\
                     all matched genes are aggregated into one pileup."
    )]
    genes: Vec<Box<str>>,

    #[arg(
        long = "regions",
        visible_alias = "region",
        value_delimiter = ',',
        help = "Genomic regions to pile up: comma-separated `chr:lb-ub`",
        long_help = "Genomic regions to pile up:\n\
                     comma-separated `chr:lb-ub` (`chr17:1000-2000,chr1:50-99`).\n\
                     Selects rows by position, with or without `--genes`.\n\
                     Without either, the browser opens a gene list to pick from;\n\
                     with --batch-process, at least one is required."
    )]
    regions: Vec<Box<str>>,

    #[arg(
        short = 's',
        long = "sites-parquet",
        help = "Site-level parquet file (from dartseq or atoi output) for the second track"
    )]
    site_file: Option<Box<str>>,

    #[arg(
        long,
        value_enum,
        default_value = "sum",
        help = "Signal aggregation mode for sparse matrix track"
    )]
    signal: PileupSignal,

    #[arg(
        long = "site-signal",
        value_enum,
        default_value = "wt-coverage",
        help = "Signal for site-level parquet track"
    )]
    site_signal: SiteSignal,

    #[arg(
        short = 'n',
        long = "bins",
        default_value_t = 80,
        help = "Number of bins along the gene body"
    )]
    num_bins: usize,

    #[arg(
        long = "height",
        default_value_t = 20,
        help = "Height of ASCII plot in terminal rows (per track)"
    )]
    plot_height: usize,

    #[arg(short, long, help = "Output TSV file path (optional)")]
    output: Option<Box<str>>,

    #[arg(long, help = "Suppress ASCII plot")]
    quiet: bool,

    #[arg(
        long = "depth",
        value_name = "FILE",
        num_args = 1..,
        help = "`_depth` matrices for a read-depth row in the browser: each bin's reads summed over cells and files"
    )]
    depth_files: Vec<Box<str>>,

    #[arg(
        long = "track",
        value_name = "LABEL=PATTERN",
        help = "Group input files into a labelled track (repeatable): each file joins the first track whose PATTERN (`*`, `?`) matches its path",
        long_help = "Group the input files into labelled tracks, e.g.\n\
                     `--track wt='*_wt_*' --track mut='*_mut_*'`. Each file joins the first\n\
                     track whose PATTERN matches its path (`*` any run, `?` one character);\n\
                     files matching none are an error. Without `--track`, every file is\n\
                     pooled into one `matrix` track. Tracks are drawn, written and browsed\n\
                     separately; files within a track are pooled per position."
    )]
    tracks: Vec<Box<str>>,

    #[arg(
        long = "batch-process",
        visible_alias = "batch-mode",
        default_value_t = false,
        help = "Print the ASCII pileup (and write -o) instead of browsing it full screen",
        long_help = "Print the ASCII pileup (and write -o) instead of browsing it full screen.\n\
                     By default `faba pileup` opens a browser: pan, zoom and jump between\n\
                     sites, with a gene list to pick from when no --genes/--regions are given.\n\
                     Without a terminal it prints, as this flag does."
    )]
    batch_process: bool,

    #[arg(
        long = "separate",
        default_value_t = false,
        help = "Give each named matrix its own track, named after its batch, rather than one together",
        long_help = "Give each named matrix its own track, named after its batch, rather than\n\
                     summing them into one. The browser's file picker sets this with Tab.\n\
                     --track groups them otherwise."
    )]
    separate: bool,

    /// The browser is the default now; kept so older command lines still run.
    #[arg(
        short = 'I',
        long = "interactive",
        hide = true,
        default_value_t = false
    )]
    interactive: bool,

    ///////////////////////
    // Miami figure mode //
    ///////////////////////
    // Providing any of `--gtf`, `--bam`, `--format`, `--svg`, or `--png`
    // switches `pileup` from the ASCII histogram to a faceted Miami plot
    // (epi sites up / gene model / read depth down, one panel per cell
    // type). Otherwise the existing ASCII/TSV behavior is unchanged.
    #[arg(
        long = "gtf",
        help = "Gene annotation GTF/GFF for the middle gene-model track (exons, introns)",
        long_help = "Gene annotation GTF/GFF for the middle gene-model track (exons, introns, strand).\n\
                     Enables figure mode with --batch-process; in the browser, draws the genes in view.\n\
                     Without it, the browser and the figure use the GFF recorded in the\n\
                     `*.run.json` next to the first matrix, when there is one."
    )]
    gtf: Option<Box<str>>,

    #[arg(
        long = "bam",
        num_args = 1..,
        help = "BAM file(s) for the bottom read-depth track",
        long_help = "BAM file(s) for the bottom read-depth track. Repeatable;\n\
                     replicates are pooled. Enables figure mode."
    )]
    bam_files: Vec<Box<str>>,

    #[arg(
        long = "cell-membership",
        visible_alias = "membership",
        help = "Cell barcode -> cell type membership (TSV/CSV/Parquet)",
        long_help = "Cell barcode -> cell type membership (TSV/CSV/Parquet).\n\
                     Panels are faceted by cell type;\n\
                     without it the figure is a single \"all cells\" panel."
    )]
    cell_membership_file: Option<Box<str>>,

    #[arg(
        long = "membership-barcode-col",
        default_value_t = 0,
        help = "Column index of the cell barcode in the membership file"
    )]
    membership_barcode_col: usize,

    #[arg(
        long = "membership-celltype-col",
        default_value_t = 1,
        help = "Column index of the cell type in the membership file"
    )]
    membership_celltype_col: usize,

    #[arg(
        long = "exact-barcode-match",
        default_value_t = false,
        help = "Require exact barcode matching",
        long_help = "Require exact barcode matching.\n\
                     By default, membership barcodes match as prefixes of BAM or matrix barcodes,\n\
                     which handles \"-1\" suffixes."
    )]
    exact_barcode_match: bool,

    #[arg(
        long = "cell-barcode-tag",
        default_value = "CB",
        help = "BAM tag holding the cell barcode (read-depth track)"
    )]
    cell_barcode_tag: Box<str>,

    #[arg(
        long = "top-modality",
        value_delimiter = ',',
        help = "Restrict the top track to these modalities (e.g. `m6A,A-to-I`)",
        long_help = "Restrict the top track to these modalities (e.g. `m6A,A-to-I`).\n\
                     Empty = all matrix rows.\n\
                     Matched case-insensitively against the `gene/MODALITY/detail` row name."
    )]
    top_modality: Vec<Box<str>>,

    #[arg(
        long = "out",
        help = "Output prefix for figure files (`<prefix>.miami.{svg,pdf,png}`)",
        long_help = "Output prefix for figure files (`<prefix>.miami.{svg,pdf,png}`).\n\
                     Defaults to the gene label."
    )]
    out: Option<Box<str>>,

    #[arg(
        long = "format",
        value_enum,
        help = "Emit only this format. Omit for the default SVG + PDF"
    )]
    format: Option<FigFormat>,

    #[arg(
        long = "svg",
        default_value_t = false,
        help = "Also write the SVG (always written unless `--format pdf`)"
    )]
    svg: bool,

    #[arg(
        long = "png",
        default_value_t = false,
        help = "Also write a flattened PNG"
    )]
    png: bool,

    #[arg(long = "no-pdf", default_value_t = false, help = "Skip PDF output")]
    no_pdf: bool,

    #[arg(
        long = "fig-width",
        default_value_t = 8.0,
        help = "Figure width in inches"
    )]
    fig_width: f32,

    #[arg(
        long = "dpi",
        default_value_t = 300,
        help = "Figure resolution (dots per inch)"
    )]
    dpi: u32,

    #[arg(
        long = "palette",
        value_enum,
        default_value = "auto",
        help = "Qualitative color palette for cell-type panels"
    )]
    palette: Palette,

    #[arg(
        long = "raster-threshold",
        default_value_t = 300,
        help = "Rasterize the per-site dot layer once a panel exceeds this many sites",
        long_help = "Rasterize the per-site dot layer once a panel exceeds this many sites.\n\
                     That keeps SVG/PDF size bounded; axes and areas stay vector."
    )]
    raster_threshold: usize,
}

/// Whether `text` matches `pattern`, `*` standing for any run of
/// characters and `?` for one.
fn wildcard(pattern: &str, text: &str) -> bool {
    let (p, t): (Vec<char>, Vec<char>) = (pattern.chars().collect(), text.chars().collect());
    let (mut pi, mut ti, mut star, mut mark) = (0, 0, None, 0);
    while ti < t.len() {
        if pi < p.len() && (p[pi] == '?' || p[pi] == t[ti]) {
            pi += 1;
            ti += 1;
        } else if pi < p.len() && p[pi] == '*' {
            star = Some(pi);
            mark = ti;
            pi += 1;
        } else if let Some(sp) = star {
            pi = sp + 1;
            mark += 1;
            ti = mark;
        } else {
            return false;
        }
    }
    p[pi..].iter().all(|&c| c == '*')
}

/// A matrix's own track name: its batch, else its file name.
fn track_name(file: &str) -> Box<str> {
    let name = crate::qc::layout::file_name(file);
    crate::qc::layout::matrix_batch(&name).unwrap_or(name)
}

/// One labelled group of input matrices.
struct TrackFiles {
    label: Box<str>,
    files: Vec<Box<str>>,
}

/// The input files grouped by `--track`, in the order the tracks were given;
/// one `matrix` track of every file without it.
fn track_files(files: &[Box<str>], specs: &[Box<str>]) -> anyhow::Result<Vec<TrackFiles>> {
    if specs.is_empty() {
        return Ok(vec![TrackFiles {
            label: "matrix".into(),
            files: files.to_vec(),
        }]);
    }
    let mut groups: Vec<(TrackFiles, &str)> = Vec::new();
    for spec in specs {
        let (label, pattern) = spec
            .split_once('=')
            .filter(|(l, p)| !l.is_empty() && !p.is_empty())
            .ok_or_else(|| anyhow::anyhow!("--track {spec}: expected LABEL=PATTERN"))?;
        let track = TrackFiles {
            label: label.into(),
            files: Vec::new(),
        };
        groups.push((track, pattern));
    }
    let mut unmatched = Vec::new();
    for f in files {
        match groups.iter_mut().find(|(_, p)| wildcard(p, f)) {
            Some((g, _)) => g.files.push(f.clone()),
            None => unmatched.push(f.as_ref()),
        }
    }
    anyhow::ensure!(
        unmatched.is_empty(),
        "no --track pattern matches {}",
        unmatched.join(", ")
    );
    for (g, p) in &groups {
        anyhow::ensure!(
            !g.files.is_empty(),
            "--track {}={p} matches no input file",
            g.label
        );
    }
    Ok(groups.into_iter().map(|(g, _)| g).collect())
}

/// One labelled matrix track: the converted channel's sorted
/// `(position, value)` pairs, and, when loaded, both channels summed.
struct MatrixTrack {
    label: Box<str>,
    converted: Vec<(i64, f64)>,
    total: Option<Vec<(i64, f64)>>,
}

/// One selection loaded for drawing: each labelled matrix track (empty when
/// none of its files hold the selection), the site annotation, the extent.
struct Loaded {
    gene: Box<str>,
    chr: Box<str>,
    modality: Box<str>,
    tracks: Vec<MatrixTrack>,
    sites: Option<SiteAnnotation>,
    extent: (i64, i64),
    /// The genes whose models to draw; `None` for a searched locus, which
    /// draws every gene in view.
    keys: Option<Vec<Box<str>>>,
    /// `--depth` bins over the extent.
    depth: Vec<(i64, i64, f64)>,
    /// Where the browser's cursor starts; the first site when `None`.
    cursor: Option<i64>,
}

impl Loaded {
    /// This view, its cursor starting `at`.
    fn at(self, at: i64) -> Self {
        Self {
            cursor: Some(at),
            ..self
        }
    }
}

/// Load `selector` for drawing. With `totals`, each matrix track also sums
/// both channels, for the converted-over-total bars (not for `nnz`, where
/// adding cells across channels would count cells twice). A `locus` fixes
/// the extent to a searched window.
fn load(
    args: &PileupArgs,
    groups: &[TrackFiles],
    selector: &Selector,
    totals: bool,
    locus: Option<(i64, i64)>,
    done: Option<&std::sync::atomic::AtomicUsize>,
) -> anyhow::Result<Loaded> {
    let totals = totals && !matches!(args.signal, PileupSignal::Nnz);
    let mut label: Option<(Box<str>, Box<str>, Box<str>)> = None;
    let mut keys: Vec<Box<str>> = Vec::new();
    let mut tracks = Vec::with_capacity(groups.len());
    // The tracks' files are independent: read them together, then take
    // them in order.
    let read: Vec<_> = {
        use rayon::prelude::*;
        groups
            .par_iter()
            .map(|g| read_matrix_positions(&g.files, selector, &args.signal, totals, done))
            .collect()
    };
    for (g, read) in groups.iter().zip(read) {
        let (converted, total) = match read? {
            Some(m) => {
                label.get_or_insert((m.gene, m.chr, m.modality));
                keys.extend(m.genes);
                (m.positions, m.total)
            }
            None => (Vec::new(), None),
        };
        tracks.push(MatrixTrack {
            label: g.label.clone(),
            converted,
            total,
        });
    }
    let Some((gene, chr, modality)) = label else {
        anyhow::bail!(
            "no rows matching {} in {} file(s)",
            selector.describe(),
            args.data_files.len()
        );
    };
    let sites = args
        .site_file
        .as_ref()
        .map(|sf| read_site_annotation(sf, selector, &args.site_signal))
        .transpose()?;
    let extent = match &sites {
        Some(sa) => (sa.gene_start, sa.gene_stop),
        None => {
            let all = tracks.iter().flat_map(|t| t.converted.iter().map(|x| x.0));
            let lo = all.clone().min().unwrap_or(0);
            (lo, all.max().unwrap_or(lo))
        }
    };
    let extent = locus.unwrap_or(extent);
    let depth = read_depth(&args.depth_files, &chr, extent)?;
    Ok(Loaded {
        gene,
        chr,
        modality,
        tracks,
        sites,
        extent,
        keys: locus.is_none().then(|| {
            keys.sort_unstable();
            keys.dedup();
            keys
        }),
        depth,
        cursor: None,
    })
}

/// [`load`], for the browser (`browser`, with the totals) under a spinner
/// saying `what` is read.
fn loading(
    browser: bool,
    args: &PileupArgs,
    groups: &[TrackFiles],
    selector: &Selector,
    locus: Option<(i64, i64)>,
    what: &str,
) -> anyhow::Result<Loaded> {
    let files: usize = groups.iter().map(|g| g.files.len()).sum();
    let text = format!("reading {what}");
    crate::tui::busy_counting(browser, "pileup", &text, files, |done| {
        load(args, groups, selector, browser, locus, Some(done))
    })
}

/// The gene models to draw: those of the matched gene `keys`, by key or
/// failing that by symbol (a GTF may version its gene ids differently), or
/// every model for a searched locus (`None`).
pub(crate) fn genes_to_draw<'a>(
    genes: &'a [GeneModel],
    keys: Option<&[Box<str>]>,
) -> Vec<&'a GeneModel> {
    let Some(keys) = keys else {
        return genes.iter().collect();
    };
    let symbols: Vec<Box<str>> = keys.iter().map(|k| query_symbol(k)).collect();
    genes
        .iter()
        .filter(|g| {
            keys.contains(&g.key) || symbols.iter().any(|s| g.symbol.eq_ignore_ascii_case(s))
        })
        .collect()
}

/// Browse `loaded` full screen, with `status` on the footer first and the
/// annotation's `genes` drawn under the tracks.
fn browse(
    args: &PileupArgs,
    loaded: &Loaded,
    genes: Option<&SharedModels>,
    status: Option<String>,
    switchable: bool,
) -> anyhow::Result<tui::Exit> {
    let is_log = matches!(args.signal, PileupSignal::Log10Sum);
    let (on, off) = channel_names(&loaded.modality);
    let stacked_name = format!("{on} / {off} ({})", args.signal.name());
    let mut tracks: Vec<tui::Track> = loaded
        .tracks
        .iter()
        .map(|t| {
            let name = match t.total {
                Some(_) => &stacked_name,
                None => args.signal.name(),
            };
            let track = tui::Track::single(&t.label, name, &t.converted, is_log);
            match &t.total {
                Some(total) => track.with_total(total),
                None => track,
            }
        })
        .collect();
    if !args.depth_files.is_empty() {
        tracks.push(tui::Track::depth("depth", &loaded.depth));
    }
    let view = tui::View {
        title: &loaded.gene,
        chr: &loaded.chr,
        extent: loaded.extent,
        on,
        off,
        keys: loaded.keys.as_deref(),
        genes: genes.cloned(),
        cursor: loaded.cursor,
        switchable,
    };
    tui::show_pileup(view, tracks, status)
}

/// Bases shown around a single searched position.
const POSITION_FLANK: i64 = 5_000;

/// Where a search led.
enum Found {
    View(Loaded),
    /// Several genes match: show the list filtered by the query.
    List,
}

/// The inputs' genes, read from their row names the first time asked.
struct Catalog {
    files: Vec<Box<str>>,
    genes: std::cell::OnceCell<Vec<picker::GeneEntry>>,
}

impl Catalog {
    fn new(files: Vec<Box<str>>) -> Self {
        Self {
            files,
            genes: std::cell::OnceCell::new(),
        }
    }

    fn get(&self) -> anyhow::Result<&[picker::GeneEntry]> {
        if let Some(genes) = self.genes.get() {
            return Ok(genes);
        }
        let what = format!("reading genes from {} file(s)", self.files.len());
        let files = &self.files;
        let list = crate::tui::busy(true, "pileup", &what, || picker::gene_catalog(files))?;
        anyhow::ensure!(!list.is_empty(), "no site rows in the input files");
        Ok(self.genes.get_or_init(|| list))
    }
}

/// Load `chr:lo-hi` for the browser, every gene in it drawn. `Err` carries
/// a message for the footer.
fn load_locus(
    args: &PileupArgs,
    groups: &[TrackFiles],
    chr: &str,
    lo: i64,
    hi: i64,
) -> Result<Loaded, String> {
    let locus = format!("{chr}:{lo}-{hi}");
    let selector = Selector::build(&[], &[locus.as_str().into()]).map_err(|e| e.to_string())?;
    loading(true, args, groups, &selector, Some((lo, hi)), &locus)
        .map_err(|_| format!("no sites in {locus}"))
}

/// Open what `query` names: a locus anywhere, or a gene of the catalog.
/// `Err` carries a message for the footer.
fn search(
    args: &PileupArgs,
    groups: &[TrackFiles],
    catalog: &Catalog,
    query: &str,
) -> Result<Found, String> {
    match parse_query(query) {
        Some(Query::Locus(r, single)) => {
            let (lo, hi) = if single {
                ((r.lb - POSITION_FLANK).max(0), r.ub + POSITION_FLANK)
            } else {
                (r.lb, r.ub)
            };
            load_locus(args, groups, &r.chr, lo, hi).map(Found::View)
        }
        Some(Query::Gene(g)) => {
            let genes = catalog.get().map_err(|e| e.to_string())?;
            let sym = query_symbol(&g);
            let hits: Vec<&picker::GeneEntry> = genes
                .iter()
                .filter(|e| gene_matches(&g, &sym, &e.gene))
                .collect();
            match hits.as_slice() {
                [] => Err(format!("no gene matches {g}")),
                [one] => loading(
                    true,
                    args,
                    groups,
                    &Selector::exact(&one.gene),
                    None,
                    &one.gene,
                )
                .map(Found::View)
                .map_err(|e| e.to_string()),
                _ => Ok(Found::List),
            }
        }
        None => Err("nothing to search for".into()),
    }
}

/// The annotation's gene models, filled by a thread so the browser opens
/// at once; views pick them up when they arrive.
pub(crate) type SharedModels = std::sync::Arc<std::sync::OnceLock<Result<Vec<GeneModel>, String>>>;

/// Start reading `gtf` on a thread; `None` without an annotation.
fn start_gene_models(gtf: Option<&str>) -> Option<SharedModels> {
    let gtf = gtf?.to_string();
    let models = SharedModels::default();
    let slot = models.clone();
    std::thread::spawn(move || {
        // Logs are paused under a view, so a failure is shown on its footer.
        let read = load_gene_models_where(&gtf, |_| true).map_err(|e| format!("{gtf}: {e}"));
        let _ = slot.set(read);
    });
    Some(models)
}

/// One modality's inputs, resolved, with their gene list read when first
/// asked for.
struct Source {
    args: PileupArgs,
    groups: Vec<TrackFiles>,
    catalog: Catalog,
    /// The sites' modality, once known.
    modality: Option<Box<str>>,
}

impl Source {
    fn new(base: &PileupArgs) -> anyhow::Result<Self> {
        let (args, groups) = base.resolved()?;
        let catalog = Catalog::new(args.data_files.clone());
        Ok(Self {
            args,
            groups,
            catalog,
            modality: base.modality.clone(),
        })
    }
}

/// `M`: the other modality's sources (resolved the first time), and `here`
/// reopened in them: the same gene, else the same span. The index of the
/// sources and the view, or why not, for the footer.
fn switch_modality(
    base: &PileupArgs,
    sources: &mut Vec<Source>,
    shown: usize,
    here: &Loaded,
) -> Result<(usize, Loaded), String> {
    sources[shown]
        .modality
        .get_or_insert_with(|| here.modality.clone());
    let to = crate::qc::layout::other_site_modality(&here.modality)
        .ok_or_else(|| format!("no other sites to switch {} to", here.modality))?;
    let i = match sources
        .iter()
        .position(|s| s.modality.as_deref() == Some(to))
    {
        Some(i) => i,
        None => {
            let other = PileupArgs {
                modality: Some(to.into()),
                ..base.clone()
            };
            sources.push(Source::new(&other).map_err(|e| e.to_string())?);
            sources.len() - 1
        }
    };
    let Source { args, groups, .. } = &sources[i];
    let reopened = match here.keys.as_deref() {
        Some([key]) => {
            loading(true, args, groups, &Selector::exact(key), None, key).map_err(|e| e.to_string())
        }
        _ => {
            let (lo, hi) = here.extent;
            load_locus(args, groups, &here.chr, lo, hi)
        }
    };
    let l = reopened.map_err(|_| format!("no {to} sites here"))?;
    Ok((i, l))
}

/// Browse `source` (resolved from `base`) interactively, starting from
/// `first` (a selection given on the command line) or from the gene list;
/// `M` switches between the m6A and A-to-I sites of the output folders
/// among `base`'s inputs.
fn interactive(base: &PileupArgs, source: Source, first: Option<Loaded>) -> anyhow::Result<()> {
    let genes = start_gene_models(source.args.annotation().as_deref());
    let switchable = base
        .data_files
        .iter()
        .any(|g| crate::run_record::output_dir_of(g).is_some());
    let mut sources = vec![source];
    let mut shown = 0;
    let mut filter = String::new();
    let mut current = first;
    let mut status: Option<String> = None;
    loop {
        let Source {
            args,
            groups,
            catalog,
            ..
        } = &sources[shown];
        let loaded = match current.take() {
            Some(l) => l,
            None => {
                let genes = catalog.get()?;
                let mut list = picker::GenePicker::new(genes);
                list.set_filter(&filter);
                list.set_status(status.take());
                let choice = list.pick()?;
                filter = list.filter().to_string();
                let query = match choice {
                    picker::Choice::Quit => return Ok(()),
                    picker::Choice::Gene(i) => {
                        let gene = &genes[i].gene;
                        let selector = Selector::exact(gene);
                        current = Some(loading(true, args, groups, &selector, None, gene)?);
                        continue;
                    }
                    picker::Choice::Locus(q) => q,
                };
                match search(args, groups, catalog, &query) {
                    Ok(Found::View(l)) => l,
                    Ok(Found::List) => continue,
                    Err(msg) => {
                        status = Some(msg);
                        continue;
                    }
                }
            }
        };
        match browse(args, &loaded, genes.as_ref(), status.take(), switchable)? {
            tui::Exit::Quit => return Ok(()),
            tui::Exit::Genes => {}
            tui::Exit::Locus { lo, hi, at } => {
                match load_locus(args, groups, &loaded.chr, lo, hi) {
                    // Wider, but where it was: the same title and cursor.
                    Ok(l) => {
                        current = Some(Loaded {
                            gene: loaded.gene.clone(),
                            ..l.at(at)
                        })
                    }
                    Err(msg) => {
                        status = Some(msg);
                        current = Some(loaded);
                    }
                }
            }
            tui::Exit::Modality { at } => {
                let here = loaded.at(at);
                match switch_modality(base, &mut sources, shown, &here) {
                    Ok((i, l)) => {
                        shown = i;
                        current = Some(l.at(at));
                    }
                    Err(msg) => {
                        status = Some(msg);
                        current = Some(here);
                    }
                }
            }
            tui::Exit::Search(q) => match search(args, groups, catalog, &q) {
                Ok(Found::View(l)) => current = Some(l),
                Ok(Found::List) => filter = q,
                Err(msg) => {
                    status = Some(msg);
                    current = Some(loaded);
                }
            },
        }
    }
}

impl PileupArgs {
    /// These arguments with each output directory (or run record) among
    /// the inputs replaced by its site matrices of `--modality` and its
    /// site table as `-s` when none was given; and the tracks: `--track`'s,
    /// else with a directory or `--separate` one per batch and one per file
    /// named, else one of them all.
    fn resolved(&self) -> anyhow::Result<(PileupArgs, Vec<TrackFiles>)> {
        let mut out = self.clone();
        let mut files = Vec::new();
        let mut own: Vec<TrackFiles> = Vec::new();
        let mut from_dir = false;
        for f in &self.data_files {
            let Some(dir) = crate::run_record::output_dir_of(f) else {
                files.push(f.clone());
                own.push(TrackFiles {
                    label: track_name(f),
                    files: vec![f.clone()],
                });
                continue;
            };
            from_dir = true;
            let sites =
                crate::site_analysis::output_dir::output_sites(&dir, self.modality.as_deref())?;
            anyhow::ensure!(
                !sites.matrices.is_empty(),
                "{} has no {}_site matrices to pile up",
                dir.display(),
                sites.modality
            );
            info!(
                "{}: {} {}_site matrices",
                dir.display(),
                sites.matrices.len(),
                sites.modality
            );
            for (batch, path) in sites.matrices {
                files.push(path.clone());
                own.push(TrackFiles {
                    label: batch,
                    files: vec![path],
                });
            }
            if out.site_file.is_none() {
                out.site_file = sites.site_table;
            }
        }
        out.data_files = files;
        let groups = if (from_dir || self.separate) && self.tracks.is_empty() {
            own
        } else {
            track_files(&out.data_files, &self.tracks)?
        };
        Ok((out, groups))
    }

    /// `--gtf`, or else the GFF the first matrix was made from, found in its
    /// run record.
    fn annotation(&self) -> Option<Box<str>> {
        let first = self.data_files.first()?;
        crate::run_record::explicit_or_recorded(self.gtf.as_deref(), first, "gff", "gene models")
    }
}

pub fn run_pileup(args: &PileupArgs) -> anyhow::Result<()> {
    let mut base = args.clone();
    if base.data_files.is_empty() {
        let ask = crate::site_analysis::input_picker::ask_inputs(
            "pileup",
            args.batch_process,
            "site matrices",
            crate::qc::layout::is_site_matrix_name,
            true,
        );
        let Some(chosen) = ask? else { return Ok(()) };
        base.data_files = chosen.paths;
        base.separate = chosen.separate;
    }
    let base = &base;
    let source = Source::new(base)?;
    let (args, groups) = (&source.args, &source.groups);
    // A figure flag makes the Miami figure; otherwise the browser, unless
    // asked to print (or there is no terminal to browse in), when `--gtf`
    // makes the figure too. In the browser `--gtf` draws the genes in view.
    let figure_flags = !args.bam_files.is_empty() || args.format.is_some() || args.svg || args.png;
    let browse_it = !args.batch_process && !figure_flags;
    let terminal = data_beans::interactive::tui_available();
    if browse_it && !terminal {
        log::info!("no terminal to browse in: printing the pileup, as --batch-process does");
    }
    let interactive = browse_it && terminal;
    let figure_mode = !interactive && (args.gtf.is_some() || figure_flags);
    if interactive && args.genes.is_empty() && args.regions.is_empty() {
        return self::interactive(base, source, None);
    }
    let selector = Selector::build(&args.genes, &args.regions)?;
    if figure_mode {
        if !args.tracks.is_empty() {
            log::warn!("--track applies to the ASCII pileup, not the figure; ignoring it");
        }
        return run_miami_figure(args, &selector);
    }

    let what = selector.describe();
    let loaded = loading(interactive, args, groups, &selector, None, &what)?;
    if interactive {
        // Straight to the browser; -o still writes the binned table.
        write_pileup(args, &loaded, false)?;
        return self::interactive(base, source, Some(loaded));
    }
    write_pileup(args, &loaded, true)
}

/// With `print`, print `loaded` as ASCII histograms (unless `--quiet`); and
/// write the binned table to `-o`.
fn write_pileup(args: &PileupArgs, loaded: &Loaded, print: bool) -> anyhow::Result<()> {
    let (min_pos, max_pos) = loaded.extent;
    let max_sites = loaded
        .tracks
        .iter()
        .map(|t| t.converted.len())
        .chain(loaded.sites.iter().map(|sa| sa.num_sites))
        .max()
        .unwrap_or(0);
    let effective_bins = args.num_bins.min(max_sites.max(1));
    let is_log = matches!(args.signal, PileupSignal::Log10Sum);

    let mut pileups: Vec<BinnedPileup> = loaded
        .tracks
        .iter()
        .map(|t| BinnedPileup {
            gene: loaded.gene.clone(),
            chr: loaded.chr.clone(),
            bins: bin_positions_with_extent(&t.converted, effective_bins, min_pos, max_pos, is_log),
            sites: distinct_positions(&t.converted),
            min_pos,
            max_pos,
            num_sites: t.converted.len(),
            track_label: t.label.clone(),
            signal_name: args.signal.name(),
        })
        .collect();
    if let Some(sa) = &loaded.sites {
        pileups.push(BinnedPileup {
            gene: loaded.gene.clone(),
            chr: loaded.chr.clone(),
            bins: bin_positions_with_extent(&sa.positions, effective_bins, min_pos, max_pos, false),
            sites: distinct_positions(&sa.positions),
            min_pos,
            max_pos,
            num_sites: sa.num_sites,
            track_label: "sites".into(),
            signal_name: args.site_signal.name(),
        });
    }

    if print && !args.quiet {
        for p in &pileups {
            print_vertical_histogram(p, args.plot_height);
        }
    }

    if let Some(ref output) = args.output {
        let tracks: Vec<&BinnedPileup> = pileups.iter().collect();
        write_pileup_tsv(&tracks, output)?;
        info!("wrote pileup TSV to {}", output);
    }

    Ok(())
}

mod picker;

mod tui;

#[cfg(test)]
mod tests;
