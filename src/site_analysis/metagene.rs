//! Metagene profiles over a simplified 5'UTR / CDS / 3'UTR transcript.
//!
//! Follows MetaPlotR, so a difference from a published profile is one of data,
//! not procedure (Olarerin-George & Jaffrey, Bioinformatics 33:1563, 2017,
//! <https://doi.org/10.1093/bioinformatics/btx002>).

use super::site_io::*;
use clap::{Args, ValueEnum};
use genomic_data::gff::*;
use genomic_data::sam::Strand;
use genomic_data::transcript::{
    build_transcript_models, elect_longest_isoform, merge_intervals, TranscriptModel,
};
use log::info;
use rustc_hash::FxHashMap;
use std::io::Write;

/// Which isoforms carry the sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum IsoformPolicy {
    /// One transcript per gene, the longest spliced one. MetaPlotR's procedure.
    Longest,
    /// Every coding transcript. A site in several isoforms is counted in each.
    All,
}

#[derive(Args, Debug)]
pub struct MetageneArgs {
    #[arg(
        short = 's',
        long = "sites",
        required = true,
        help = "Site-level parquet file (from dartseq, atoi or apa output)"
    )]
    site_file: Box<str>,

    #[arg(
        short = 'g',
        long = "gff",
        help = "GFF annotation file (default: the one recorded next to --sites)",
        long_help = "GFF annotation file.\n\
                     Without it, the GFF recorded in the `*.run.json` next to --sites is used."
    )]
    gff_file: Option<Box<str>>,

    #[arg(
        short = 'n',
        long = "bins",
        default_value_t = 200,
        value_parser = clap::value_parser!(u32).range(1..=1_000_000),
        help = "Total bins across the metagene (default: 200)",
        long_help = "Total number of bins across the 5'UTR, CDS and 3'UTR.\n\
                     \n\
                     MetaPlotR plots its metagene with 200 breaks, hence the default.\n\
                     Bins are split between the three regions in proportion,\n\
                     by each region's MEDIAN spliced length over the assigned sites.\n\
                     The median rather than the maximum, which one gene would set:\n\
                     titin's merged CDS is 114,586 nt against a median of 1,347.\n\
                     A region that has sites always keeps at least one bin.\n\
                     The split depends on the annotation and on the sites,\n\
                     so compare the shape of two profiles rather than their bar widths."
    )]
    num_bins: u32,

    #[arg(
        long = "isoforms",
        value_enum,
        default_value = "longest",
        help = "Which isoforms sites are placed on",
        long_help = "Which isoforms sites are placed on.\n\
                     \n\
                     `longest` keeps one transcript per gene, the longest spliced one.\n\
                     That is MetaPlotR's stated procedure. `all` keeps every coding transcript.\n\
                     A site inside several isoforms is then counted once per isoform,\n\
                     which is what MetaPlotR's own distance table does."
    )]
    isoforms: IsoformPolicy,

    #[arg(
        long = "include-non-coding",
        help = "Also profile non-coding genes, as a separate ncRNA track",
        long_help = "Also profile non-coding genes, as a separate ncRNA track.\n\
                     \n\
                     This has no counterpart in MetaPlotR,\n\
                     which profiles coding transcripts only.\n\
                     A non-coding gene has no start or stop codon to split on.\n\
                     Its whole body becomes one undivided track on its own [0,1] axis,\n\
                     and its density is normalized within that track alone."
    )]
    include_non_coding: bool,

    #[arg(short, long, required = true, help = "Output TSV file path")]
    output: Box<str>,

    #[arg(
        long = "dist-measures",
        help = "Also write MetaPlotR's per-site distance table to this path",
        long_help = "Also write the per-site distance table,\n\
                     one row per site and transcript it was placed on.\n\
                     \n\
                     Column names match MetaPlotR's own `dist_measures` output, so its\n\
                     `visualize_metagenes.R` runs on this file unmodified.\n\
                     That turns \"our profile looks like theirs\" into \"their script, our data\"."
    )]
    dist_measures: Option<Box<str>>,

    #[arg(long = "print", help = "Print ASCII histogram to stderr")]
    print_histogram: bool,

    #[arg(
        short = 'I',
        long = "interactive",
        default_value_t = false,
        help = "After writing, browse the profile full screen (needs a terminal)"
    )]
    interactive: bool,

    #[arg(
        long = "max-width",
        default_value_t = 60,
        value_parser = clap::value_parser!(u32).range(1..),
        help = "Maximum width of ASCII histogram"
    )]
    max_width: u32,
}

/// TSV feature labels in report order; an output contract (no apostrophes).
const FEATURE_LABELS: [&str; 4] = ["5UTR", "CDS", "3UTR", "ncRNA"];

/// On-screen region names (the TSV's `FEATURE_LABELS` avoid apostrophes).
pub(crate) const REGION_NAMES: [&str; 4] = ["5'UTR", "CDS", "3'UTR", "ncRNA"];

/// Region indices into [`FEATURE_LABELS`], and the base of each region's
/// MetaPlotR coordinate: 5'UTR spans [0,1), CDS [1,2), 3'UTR [2,3).
const UTR5: usize = 0;
const CDS: usize = 1;
const UTR3: usize = 2;
const NCRNA: usize = 3;

/// One interval of one region of one transcript.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct IndexedInterval {
    start: i64,
    stop: i64,
    strand: Strand,
    /// Spliced length of the same region lying genomically before this interval.
    cum_before: i64,
    /// Spliced length of the region this interval belongs to.
    total_len: i64,
    /// Which of [`FEATURE_LABELS`] this interval is.
    region: usize,
    /// Index into the model table; `u32::MAX` for the ncRNA track.
    model: u32,
}

impl IndexedInterval {
    /// 0-based offset of `pos` along the spliced region, read 5'->3'.
    fn relative_pos(&self, pos: i64) -> i64 {
        let rel_genomic = self.cum_before + (pos - self.start);
        let rel = match self.strand {
            Strand::Forward => rel_genomic,
            Strand::Backward => self.total_len - 1 - rel_genomic,
        };
        rel.clamp(0, (self.total_len - 1).max(0))
    }

    /// This interval's placement of `pos`, ready to bin.
    fn place(&self, site: u32, pos: i64) -> SiteAssignment {
        SiteAssignment {
            site,
            model: (self.region != NCRNA).then_some(self.model),
            region: self.region,
            rel: self.relative_pos(pos),
            total_len: self.total_len,
        }
    }
}

/// One chromosome's intervals, sorted by start.
struct ChromIntervals {
    intervals: Vec<IndexedInterval>,
    /// Largest `stop` among `intervals[..=i]`; lets a backward scan stop early.
    max_stop: Vec<i64>,
}

/// Per-chromosome sorted interval index over every region of every transcript.
struct RegionIndex {
    by_chr: FxHashMap<Box<str>, ChromIntervals>,
}

impl RegionIndex {
    fn build(models: &[TranscriptModel], non_coding: &[NonCodingBody]) -> Self {
        let mut by_chr: FxHashMap<Box<str>, Vec<IndexedInterval>> = FxHashMap::default();

        for (mi, m) in models.iter().enumerate() {
            let chrom = by_chr.entry(m.seqname.clone()).or_default();
            for (region, intervals, total_len) in [
                (UTR5, &m.utr5, m.utr5_size),
                (CDS, &m.cds, m.cds_size),
                (UTR3, &m.utr3, m.utr3_size),
            ] {
                let mut cum_before = 0;
                for &(start, stop) in intervals.iter() {
                    chrom.push(IndexedInterval {
                        start,
                        stop,
                        strand: m.strand,
                        cum_before,
                        total_len,
                        region,
                        model: mi as u32,
                    });
                    cum_before += stop - start + 1;
                }
            }
        }

        for body in non_coding.iter() {
            let chrom = by_chr.entry(body.seqname.clone()).or_default();
            let total_len: i64 = body.intervals.iter().map(|&(s, e)| e - s + 1).sum();
            let mut cum_before = 0;
            for &(start, stop) in body.intervals.iter() {
                chrom.push(IndexedInterval {
                    start,
                    stop,
                    strand: body.strand,
                    cum_before,
                    total_len,
                    region: NCRNA,
                    model: u32::MAX,
                });
                cum_before += stop - start + 1;
            }
        }

        let by_chr = by_chr
            .into_iter()
            .map(|(chr, mut intervals)| {
                // Total order, without run-dependent `model`, for reproducibility.
                intervals.sort_by_key(|iv| {
                    (
                        iv.start,
                        iv.stop,
                        iv.region,
                        iv.total_len,
                        iv.cum_before,
                        iv.strand,
                    )
                });
                let mut running = i64::MIN;
                let max_stop = intervals
                    .iter()
                    .map(|iv| {
                        running = running.max(iv.stop);
                        running
                    })
                    .collect();
                (
                    chr,
                    ChromIntervals {
                        intervals,
                        max_stop,
                    },
                )
            })
            .collect();

        RegionIndex { by_chr }
    }

    /// Every interval on `strand` containing `position` (1-based GFF coords),
    /// as MetaPlotR's `intersectBed -wo -s`.
    fn find_all(&self, chr: &str, position: i64, strand: Strand, out: &mut Vec<IndexedInterval>) {
        out.clear();
        let Some(chrom) = self.by_chr.get(chr) else {
            return;
        };
        // Every interval left of `idx` starts at or before `position`.
        let idx = chrom.intervals.partition_point(|iv| iv.start <= position);
        for i in (0..idx).rev() {
            if chrom.max_stop[i] < position {
                break;
            }
            let iv = chrom.intervals[i];
            if position <= iv.stop && iv.strand == strand {
                out.push(iv);
            }
        }
    }
}

/// One site placed on one transcript, independent of the bin grid.
struct SiteAssignment {
    /// Index into the site list.
    site: u32,
    /// `None` on the ncRNA track, which has no transcript.
    model: Option<u32>,
    region: usize,
    /// 0-based offset along the spliced region and its length (exact binning).
    rel: i64,
    total_len: i64,
}

impl SiteAssignment {
    /// Bin within a track of `nbins` bins.
    fn bin(&self, nbins: usize) -> usize {
        let total = self.total_len.max(1) as usize;
        ((self.rel as usize) * nbins / total).min(nbins.saturating_sub(1))
    }

    /// MetaPlotR's coordinate: 5'UTR [0,1), CDS [1,2), 3'UTR [2,3); ncRNA [0,1).
    fn rel_location(&self) -> f64 {
        let base = if self.region == NCRNA {
            0.0
        } else {
            self.region as f64
        };
        base + self.rel as f64 / self.total_len.max(1) as f64
    }
}

/// Collect every (site, transcript) placement.
fn assign_sites(sites: &[GenomicSite], index: &RegionIndex) -> (Vec<SiteAssignment>, usize) {
    let mut out = Vec::new();
    let mut hits = Vec::new();
    let mut unassigned = 0usize;

    for (si, site) in sites.iter().enumerate() {
        // Sites use 0-based positions; GFF uses 1-based.
        let gff_pos = site.position + 1;
        index.find_all(site.chr.as_ref(), gff_pos, site.strand, &mut hits);
        if hits.is_empty() {
            unassigned += 1;
            continue;
        }
        // ncRNA is a fallback, never a second placement of a coding site.
        let on_coding = hits.iter().any(|iv| iv.region != NCRNA);
        for iv in hits.iter() {
            if on_coding && iv.region == NCRNA {
                continue;
            }
            out.push(iv.place(si as u32, gff_pos));
        }
    }
    (out, unassigned)
}

/// MetaPlotR's display widths: each UTR's median size relative to the CDS's.
/// Medians are site-weighted, as in `visualize_metagenes.R`.
struct ScaleFactors {
    /// Each region's median size, doubled so an even-n median stays exact.
    twice_median: [i64; 3],
    utr5_sf: f64,
    utr3_sf: f64,
}

impl ScaleFactors {
    /// Place a within-region fraction on MetaPlotR's rescaled axis (CDS width 1).
    fn rescale(&self, region: usize, within: f64) -> f64 {
        match region {
            UTR5 => 1.0 - self.utr5_sf * (1.0 - within),
            CDS => 1.0 + within,
            UTR3 => 2.0 + self.utr3_sf * within,
            _ => within, // ncRNA: its own [0,1] axis
        }
    }

    /// Stand-in when no site is coding; the coding tracks are then zero-wide.
    fn none() -> Self {
        Self {
            twice_median: [0; 3],
            utr5_sf: 1.0,
            utr3_sf: 1.0,
        }
    }

    /// Median region sizes, for reporting only.
    fn median(&self) -> [f64; 3] {
        [
            self.twice_median[0] as f64 / 2.0,
            self.twice_median[1] as f64 / 2.0,
            self.twice_median[2] as f64 / 2.0,
        ]
    }
}

/// Median, doubled so it stays an integer for exact bin allocation.
fn twice_median(values: &mut [i64]) -> i64 {
    if values.is_empty() {
        return 0;
    }
    let n = values.len();
    let (lo, mid, _) = values.select_nth_unstable(n / 2);
    if n % 2 == 1 {
        2 * *mid
    } else {
        // `lo` holds the n/2 smallest; its max is the lower middle value.
        *mid + *lo.iter().max().expect("n >= 2 when n is even and non-zero")
    }
}

fn scale_factors(
    assignments: &[SiteAssignment],
    models: &[TranscriptModel],
) -> Option<ScaleFactors> {
    let mut per_region: [Vec<i64>; 3] = Default::default();
    for a in assignments.iter() {
        let Some(mi) = a.model else {
            continue; // ncRNA: no transcript, no region sizes
        };
        let m = &models[mi as usize];
        per_region[UTR5].push(m.utr5_size);
        per_region[CDS].push(m.cds_size);
        per_region[UTR3].push(m.utr3_size);
    }
    let mut m = [0i64; 3];
    for r in 0..3 {
        m[r] = twice_median(&mut per_region[r]);
    }
    if m[CDS] == 0 {
        // No coding axis; not an error, as the ncRNA track needs none.
        return None;
    }
    Some(ScaleFactors {
        twice_median: m,
        utr5_sf: m[UTR5] as f64 / m[CDS] as f64,
        utr3_sf: m[UTR3] as f64 / m[CDS] as f64,
    })
}

/// Split `n` bins between the three regions in proportion to their medians,
/// by exact largest remainder. A region with sites keeps at least one bin,
/// or its sites would silently vanish from the counts and the density.
fn allocate_bins(n: usize, m: &[i64; 3]) -> [usize; 3] {
    let total: i64 = m.iter().sum();
    if total <= 0 || n == 0 {
        return [0, n, 0];
    }
    // i128 so a large --bins times a genomic length cannot overflow.
    let total = total as i128;
    let mut out = [0usize; 3];
    let mut rem = [(0i128, 0usize); 3];
    let mut used = 0usize;
    for r in 0..3 {
        let exact = m[r] as i128 * n as i128;
        out[r] = (exact / total) as usize;
        used += out[r];
        rem[r] = (exact % total, r);
    }
    // Largest remainder first; ties to the wider region, then to the earlier.
    rem.sort_by_key(|&(f, r)| (std::cmp::Reverse(f), std::cmp::Reverse(m[r]), r));
    for &(_, r) in rem.iter().take(n.saturating_sub(used)) {
        out[r] += 1;
    }

    // Floor every represented region at one bin, paying from the widest.
    for r in 0..3 {
        if m[r] > 0 && out[r] == 0 {
            let donor = (0..3)
                .filter(|&d| out[d] > 1)
                .max_by_key(|&d| out[d])
                .unwrap_or(r);
            if donor != r {
                out[donor] -= 1;
                out[r] = 1;
            }
        }
    }
    out
}

/// The bin count of every track, decided in one constructor.
struct BinGrid([usize; 4]);

impl BinGrid {
    /// Placements per bin, one row per track, each adding its weight.
    fn tally<'a>(
        &self,
        weighted: impl Iterator<Item = (&'a SiteAssignment, usize)>,
    ) -> [Vec<usize>; 4] {
        let mut counts: [Vec<usize>; 4] =
            std::array::from_fn(|region| vec![0usize; self.0[region]]);
        for (a, w) in weighted {
            let track = &mut counts[a.region];
            let width = track.len();
            if width > 0 {
                track[a.bin(width)] += w;
            }
        }
        counts
    }

    fn new(n: usize, scale: Option<&ScaleFactors>, include_non_coding: bool) -> Self {
        let coding = match scale {
            Some(s) => allocate_bins(n, &s.twice_median),
            None => [0usize; 3],
        };
        // ncRNA: whole budget on its own axis, or 0 (no TSV rows) if not asked.
        BinGrid([
            coding[UTR5],
            coding[CDS],
            coding[UTR3],
            if include_non_coding { n } else { 0 },
        ])
    }
}

pub struct GeneFeatureHistogram {
    /// One row of bins per track; each row's length is its bin count.
    counts: [Vec<usize>; 4],
    scale: ScaleFactors,
}

impl GeneFeatureHistogram {
    /// Tally every placement, once the grid has fixed the bin widths.
    fn accumulate(grid: &BinGrid, scale: ScaleFactors, assignments: &[SiteAssignment]) -> Self {
        let counts = grid.tally(assignments.iter().map(|a| (a, 1)));
        GeneFeatureHistogram { counts, scale }
    }

    /// Rescaled coordinate spanned by one bin, as MetaPlotR's plot draws it.
    fn bin_edges(&self, region: usize, i: usize) -> (f64, f64) {
        let b = self.counts[region].len().max(1) as f64;
        (
            self.scale.rescale(region, i as f64 / b),
            self.scale.rescale(region, (i + 1) as f64 / b),
        )
    }

    pub fn print(&self, max_width: usize) {
        let nmax = self
            .counts
            .iter()
            .flat_map(|c| c.iter())
            .cloned()
            .max()
            .unwrap_or(0);
        if nmax == 0 {
            eprintln!("(no sites mapped to gene features)");
            return;
        }
        let scale = nmax.div_ceil(max_width);
        for (region, data) in self.counts.iter().enumerate() {
            for &n in data.iter() {
                let n1 = n.div_ceil(scale);
                let n0 = max_width.saturating_sub(n1);
                eprintln!(
                    "{:<6}{}{} {}",
                    FEATURE_LABELS[region],
                    "*".repeat(n1),
                    " ".repeat(n0),
                    n
                );
            }
        }
    }

    pub fn to_tsv(&self, file_path: &str) -> anyhow::Result<()> {
        let mut writer = legume_numeric::matrix::common_io::open_buf_writer(file_path)?;
        // Contract: scripts read the first three columns positionally.
        writeln!(
            writer,
            "#feature\tgenomic_bin\tcount\tbin_start\tbin_end\tfrac\tdensity"
        )?;

        // Coding regions share one density; ncRNA normalizes within itself.
        let coding_total: usize = self.counts[..3].iter().flat_map(|c| c.iter()).sum();
        let nc_total: usize = self.counts[NCRNA].iter().sum();

        for (region, data) in self.counts.iter().enumerate() {
            let total = if region == NCRNA {
                nc_total
            } else {
                coding_total
            };
            for (i, &n) in data.iter().enumerate() {
                let (lo, hi) = self.bin_edges(region, i);
                let width = hi - lo;
                let (frac, density) = if total == 0 || width <= 0.0 {
                    (0.0, 0.0)
                } else {
                    let f = n as f64 / total as f64;
                    (f, f / width)
                };
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                    FEATURE_LABELS[region], i, n, lo, hi, frac, density
                )?;
            }
        }
        writer.flush()?;
        Ok(())
    }
}

/// MetaPlotR's `*.dist.measures.txt` schema, so its `visualize_metagenes.R`
/// runs on this file unmodified.
///
/// The first fourteen columns are `rel_and_abs_dist_calc.pl`'s, in its order;
/// `strand` and `rescaled_location` are appended. The six `_st`/`_end` columns
/// are absolute distances `mrna_pos - endpoint` in 1-based spliced coordinates
/// running 5'->3' (so `utr3_st` is the signed distance from the stop codon);
/// a missing region prints `NA`. `coord` is 1-based, as MetaPlotR reads it
/// from a BED `end`; the site parquet stores 0-based positions.
const DIST_MEASURES_HEADER: &str = "chr\tcoord\tgene_name\trefseqID\trel_location\t\
     utr5_st\tutr5_end\tcds_st\tcds_end\tutr3_st\tutr3_end\t\
     utr5_size\tcds_size\tutr3_size\tstrand\trescaled_location";

fn write_dist_measures(
    path: &str,
    sites: &[GenomicSite],
    assignments: &[SiteAssignment],
    models: &[TranscriptModel],
    scale: &ScaleFactors,
) -> anyhow::Result<()> {
    let mut w = legume_numeric::matrix::common_io::open_buf_writer(path)?;
    writeln!(w, "{}", DIST_MEASURES_HEADER)?;

    for a in assignments.iter() {
        let Some(mi) = a.model else {
            continue; // ncRNA has no transcript, and no MetaPlotR counterpart
        };
        let m = &models[mi as usize];
        let site = &sites[a.site as usize];
        let rel_location = a.rel_location();

        // 1-based position along the mature transcript.
        let preceding = match a.region {
            UTR5 => 0,
            CDS => m.utr5_size,
            _ => m.utr5_size + m.cds_size,
        };
        let mrna_pos = preceding + a.rel + 1;

        // Region boundaries in the same frame, inclusive on both ends.
        let bounds = [
            (m.utr5_size > 0).then_some((1, m.utr5_size)),
            (m.cds_size > 0).then_some((m.utr5_size + 1, m.utr5_size + m.cds_size)),
            (m.utr3_size > 0).then_some((
                m.utr5_size + m.cds_size + 1,
                m.utr5_size + m.cds_size + m.utr3_size,
            )),
        ];
        let mut abs = String::new();
        for b in bounds.iter() {
            match b {
                Some((st, end)) => {
                    abs.push_str(&format!("{}\t{}\t", mrna_pos - st, mrna_pos - end))
                }
                None => abs.push_str("NA\tNA\t"),
            }
        }

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.6}\t{}{}\t{}\t{}\t{}\t{:.6}",
            site.chr,
            site.position + 1,
            m.gene_name,
            m.transcript_id,
            rel_location,
            abs,
            m.utr5_size,
            m.cds_size,
            m.utr3_size,
            m.strand,
            scale.rescale(a.region, rel_location - a.region as f64)
        )?;
    }
    w.flush()?;
    Ok(())
}

/// One non-coding gene's merged exons.
struct NonCodingBody {
    seqname: Box<str>,
    strand: Strand,
    intervals: Vec<(i64, i64)>,
}

/// Merged exons (not the gene span, so introns are excluded) per non-coding gene.
fn non_coding_bodies(records: &[GffRecord]) -> Vec<NonCodingBody> {
    // Keyed on sequence name too, so pseudoautosomal copies stay apart.
    let mut by_gene: FxHashMap<(GeneId, Box<str>), NonCodingBody> = FxHashMap::default();
    for rec in records.iter() {
        if rec.gene_type == GeneType::CodingGene
            || rec.feature_type != FeatureType::Exon
            || rec.stop < rec.start
        {
            continue;
        }
        by_gene
            .entry((rec.gene_id.clone(), rec.seqname.clone()))
            .or_insert_with(|| NonCodingBody {
                seqname: rec.seqname.clone(),
                strand: rec.strand,
                intervals: Vec::new(),
            })
            .intervals
            .push((rec.start, rec.stop));
    }
    by_gene
        .into_values()
        .map(|mut b| {
            merge_intervals(&mut b.intervals);
            b
        })
        .collect()
}

pub fn run_metagene(args: &MetageneArgs) -> anyhow::Result<()> {
    let sites = read_sites(&args.site_file)?;
    let gff_file = crate::run_record::explicit_or_recorded(
        args.gff_file.as_deref(),
        &args.site_file,
        "gff",
        "annotation",
    )
    .ok_or_else(|| {
        anyhow::anyhow!(
            "-g/--gff is required: no run record next to {} names a GFF",
            args.site_file
        )
    })?;
    let records = read_gff_record_vec(&gff_file)?;

    let models = build_transcript_models(&records);
    if models.is_empty() {
        // Models are built from `exon` records; name that cause if missing.
        let has_exon = records.iter().any(|r| r.feature_type == FeatureType::Exon);
        let has_cds = records.iter().any(|r| r.feature_type == FeatureType::CDS);
        if has_cds && !has_exon {
            anyhow::bail!(
                "{} has CDS records but no `exon` records, and the transcript model is \
                 built from exons. GENCODE, Ensembl and RefSeq all emit exon lines; a \
                 CDS/UTR-only or hand-subsetted annotation does not.",
                gff_file
            );
        }
        anyhow::bail!(
            "no coding transcript could be built from {}. Coding transcripts need \
             `exon` and `CDS` records carrying both `gene_type`/`gene_biotype` \
             protein_coding and a `transcript_id` attribute.",
            gff_file
        );
    }
    let n_models = models.len();
    let models = match args.isoforms {
        IsoformPolicy::Longest => elect_longest_isoform(models),
        IsoformPolicy::All => models,
    };
    let non_coding = if args.include_non_coding {
        non_coding_bodies(&records)
    } else {
        Vec::new()
    };
    drop(records);

    info!(
        "{} coding transcripts, {} kept under --isoforms {:?}",
        n_models,
        models.len(),
        args.isoforms
    );

    let index = RegionIndex::build(&models, &non_coding);

    // Placement first: bin widths depend on the site-weighted medians.
    let (assignments, unassigned) = assign_sites(&sites, &index);
    let scale = scale_factors(&assignments, &models);
    let grid = BinGrid::new(
        args.num_bins as usize,
        scale.as_ref(),
        args.include_non_coding,
    );
    let nbins = grid.0;

    if scale.is_none() {
        if nbins[NCRNA] == 0 {
            anyhow::bail!(
                "no site was placed on a coding transcript, so there is nothing to profile. \
                 MetaPlotR's bin widths are all relative to the median CDS, so a coding axis \
                 needs at least one coding assignment. Pass --include-non-coding to profile \
                 the ncRNA track instead, or check that the GFF and the sites use the same \
                 chromosome names."
            );
        }
        info!("no coding assignment: writing the ncRNA track only");
    }
    let scale = scale.unwrap_or_else(ScaleFactors::none);

    let n_rows: usize = assignments.iter().filter(|a| a.model.is_some()).count();
    info!(
        "sites {} | assigned rows {} | unassigned {} ({:.2}%)",
        sites.len(),
        n_rows,
        unassigned,
        100.0 * unassigned as f64 / sites.len().max(1) as f64
    );
    let median = scale.median();
    info!(
        "site-weighted medians 5'UTR/CDS/3'UTR = {:.1}/{:.1}/{:.1} nt; SF5 = {:.4}, SF3 = {:.4}",
        median[0], median[1], median[2], scale.utr5_sf, scale.utr3_sf
    );
    info!(
        "bins 5'UTR/CDS/3'UTR/ncRNA = {}/{}/{}/{}",
        nbins[0], nbins[1], nbins[2], nbins[3]
    );

    if let Some(path) = args.dist_measures.as_ref() {
        write_dist_measures(path, &sites, &assignments, &models, &scale)?;
        info!("wrote per-site distance table to {}", path);
    }

    let histogram = GeneFeatureHistogram::accumulate(&grid, scale, &assignments);
    histogram.to_tsv(&args.output)?;
    info!("wrote metagene histogram to {}", args.output);

    if args.print_histogram {
        histogram.print(args.max_width as usize);
    }
    if args.interactive {
        let title = crate::qc::layout::file_name(&args.site_file);
        crate::figure::term::when_terminal(|| tui::show_metagene(&title, &histogram))?;
    }

    Ok(())
}

/// Longest coding transcript per gene, indexed for placing sites.
pub struct MetaModels {
    models: Vec<TranscriptModel>,
    index: RegionIndex,
}

impl MetaModels {
    /// Whether no coding transcript could be built.
    pub fn is_empty(&self) -> bool {
        self.models.is_empty()
    }

    pub(crate) fn from_records(records: &[GffRecord]) -> Self {
        let models = elect_longest_isoform(build_transcript_models(records));
        let index = RegionIndex::build(&models, &[]);
        Self { models, index }
    }

    /// Place `sites` on a coding metagene of `n_bins` whose widths all of them
    /// fix. `None` when no site lands on a coding transcript.
    pub fn layout(&self, sites: &[GenomicSite], n_bins: usize) -> Option<MetaLayout> {
        let (assignments, unassigned) = assign_sites(sites, &self.index);
        let scale = scale_factors(&assignments, &self.models)?;
        let grid = BinGrid::new(n_bins, Some(&scale), false);
        Some(MetaLayout {
            grid,
            assignments,
            unassigned,
        })
    }
}

/// Every site placed once on a fixed axis; [`Self::counts`] re-tallies any
/// subset without recomputing region widths.
pub struct MetaLayout {
    grid: BinGrid,
    assignments: Vec<SiteAssignment>,
    /// Sites on no coding transcript.
    pub unassigned: usize,
}

impl MetaLayout {
    /// Bins per coding region: 5'UTR, CDS, 3'UTR.
    pub fn region_bins(&self) -> [usize; 3] {
        [self.grid.0[UTR5], self.grid.0[CDS], self.grid.0[UTR3]]
    }

    /// Per bin (5'UTR, CDS, 3'UTR), the sum of `weight(site index)` over
    /// placements; a site on two isoforms counts on each.
    pub fn counts(&self, weight: impl Fn(usize) -> usize) -> Vec<usize> {
        let weighted = self
            .assignments
            .iter()
            .map(|a| (a, weight(a.site as usize)));
        let [utr5, cds, utr3, _] = self.grid.tally(weighted);
        [utr5, cds, utr3].concat()
    }
}

mod tui;

#[cfg(test)]
mod tests;
