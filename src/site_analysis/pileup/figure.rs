//! The faceted Miami figure.

use super::*;

/// Faceted Miami-plot path: stratified matrix epi sites (top), GTF gene
/// model (middle), and BAM read depth (bottom), one panel per track and
/// cell type.
pub(super) fn run_miami_figure(
    args: &PileupArgs,
    groups: &[TrackFiles],
    selector: &Selector,
) -> anyhow::Result<()> {
    // Optional cell-type membership for faceting. allow_prefix = !exact.
    let membership = match &args.cell_membership_file {
        Some(p) => Some(CellMembership::from_file(
            p,
            args.membership_barcode_col,
            args.membership_celltype_col,
            !args.exact_barcode_match,
        )?),
        None => None,
    };

    // Top track: stratified matrix epi sites.
    let groups = figure_tracks(args, groups)?;
    let sites = figure_sites(args, &groups, selector, membership.as_ref())?;

    // Optional parquet refines gene bounds in figure mode.
    let site_annotation = args
        .site_file
        .as_ref()
        .map(|sf| read_site_annotation(sf, selector, &args.site_signal))
        .transpose()?;

    // Middle track: gene model(s) from GTF.
    let models = match args.annotation() {
        Some(gtf) => load_gene_models(&gtf, selector)?,
        None => Vec::new(),
    };

    // Shared extent = union of gene-model footprint, all matrix sites, and
    // parquet bounds — so the whole model and every site are visible.
    let mut lo = i64::MAX;
    let mut hi = i64::MIN;
    if let Some((_, mlo, mhi)) = models_extent(&models) {
        lo = lo.min(mlo);
        hi = hi.max(mhi);
    }
    for panel in &sites.panels {
        for &(p, _) in &panel.sites {
            lo = lo.min(p);
            hi = hi.max(p);
        }
    }
    if let Some(sa) = &site_annotation {
        lo = lo.min(sa.gene_start);
        hi = hi.max(sa.gene_stop);
    }
    if lo > hi {
        anyhow::bail!(
            "no genomic coordinates to plot for {} (need positional matrix rows, a GTF, or a sites parquet)",
            selector.describe()
        );
    }

    let edges = BinEdges::new(lo, hi, args.num_bins);

    // Region chromosome: prefer the matrix chr (same origin as the BAM),
    // else the GTF gene's chr.
    let region_chr: Box<str> = if sites.chr.as_ref() != "*"
        && sites.chr.as_ref() != "component"
        && !sites.chr.is_empty()
    {
        sites.chr.clone()
    } else if let Some(m) = models.first() {
        m.chr.clone()
    } else {
        sites.chr.clone()
    };

    // Bottom track: read depth from BAM, stratified by cell type, each
    // track's from the BAMs its pattern matches.
    let region = Bed {
        chr: region_chr.clone(),
        start: lo,
        stop: hi + 1,
    };
    let depth = track_bams(&groups, &args.bam_files)
        .iter()
        .map(|bams| {
            if bams.is_empty() {
                return Ok(FxHashMap::default());
            }
            read_depth_binned(
                bams,
                &region,
                &edges,
                &args.cell_barcode_tag,
                membership.as_ref(),
            )
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let is_log = matches!(args.signal, PileupSignal::Log10Sum);
    let panels = panel_data(sites.panels, &depth, edges.num_bins, is_log);

    // Output: SVG only (PNG/PDF previously came from plot-utils). `--svg`
    // asks for what is written, so it is not one of those ignored.
    if args.format.is_some() || args.png || args.no_pdf {
        log::warn!("Miami figure writes SVG only; `--format` / `--png` / `--no-pdf` are ignored");
    }

    let out_prefix: Box<str> = args.out.clone().unwrap_or_else(|| slug(&sites.gene).into());

    let top_label: Box<str> = if args.top_modality.is_empty() {
        "epi sites".into()
    } else {
        args.top_modality
            .iter()
            .map(|m| m.as_ref())
            .collect::<Vec<_>>()
            .join("/")
            .into()
    };

    let title: Box<str> = format!(
        "{}  {}:{}-{}",
        sites.gene,
        region_chr,
        fmt_thousands(lo),
        fmt_thousands(hi)
    )
    .into();

    let (width_in, height_in) = args.size.inches()?;
    let opts = FigOpts {
        out_prefix,
        width_in: width_in.unwrap_or(8.0) as f32,
        height_in: height_in.map(|h| h as f32),
        dpi: args.dpi,
        palette: args.palette.clone(),
        title,
        top_label,
    };

    let n = render_miami(&panels, &models, &edges, &opts)?;
    info!(
        "rendered Miami plot: {} panel(s), {} file(s)",
        panels.len(),
        n
    );
    Ok(())
}

/// The figure's tracks: `groups` as `--track` gave them, else every file
/// pooled into one (a track per batch, as a directory or `--separate`
/// gives the browser, would draw a panel per batch and cell type).
pub(super) fn figure_tracks(
    args: &PileupArgs,
    groups: &[TrackFiles],
) -> anyhow::Result<Vec<TrackFiles>> {
    if !args.tracks.is_empty() {
        return Ok(groups.to_vec());
    }
    if groups.len() > 1 {
        info!(
            "the figure pools its {} tracks into one; --track gives it a panel per track",
            groups.len()
        );
    }
    track_files(&args.data_files, &[])
}

/// Each track's BAMs: a pattern track's those its pattern matches (each
/// BAM in the first track that matches it), a pooled track's all of them.
pub(super) fn track_bams(groups: &[TrackFiles], bams: &[Box<str>]) -> Vec<Vec<Box<str>>> {
    let mut out = vec![Vec::new(); groups.len()];
    for b in bams {
        let track = groups
            .iter()
            .position(|g| g.pattern.as_deref().is_none_or(|p| wildcard(p, b)));
        match track {
            Some(i) => out[i].push(b.clone()),
            None => log::warn!("--bam {b} matches no --track pattern; its reads are not drawn"),
        }
    }
    if !bams.is_empty() {
        for (g, b) in groups.iter().zip(&out) {
            if b.is_empty() {
                info!("track {} matches no --bam: no read depth drawn", g.label);
            }
        }
    }
    out
}

/// One figure panel's top-track sites.
pub(super) struct FigurePanel {
    /// What the panel is titled; empty for all cells.
    pub(super) label: Box<str>,
    /// The index of the track the panel is of.
    pub(super) track: usize,
    /// The cell type whose read depth the panel draws; empty for all cells.
    pub(super) celltype: Box<str>,
    /// The cell type's index among the figure's, its colour.
    pub(super) hue: usize,
    /// Sorted `(position, value)` pairs.
    pub(super) sites: Vec<(i64, f64)>,
}

/// The figure's top-track sites, panel by panel, and the selection's gene
/// and chromosome over every track.
pub(super) struct FigureSites {
    pub(super) gene: Box<str>,
    pub(super) chr: Box<str>,
    pub(super) panels: Vec<FigurePanel>,
}

/// The selection's sites, a panel per track of `groups` and cell type of
/// `membership` (sorted), or per track alone without it. One track's
/// panels are named by cell type only, all cells when there is none.
pub(super) fn figure_sites(
    args: &PileupArgs,
    groups: &[TrackFiles],
    selector: &Selector,
    membership: Option<&CellMembership>,
) -> anyhow::Result<FigureSites> {
    let celltypes: Vec<Box<str>> = match membership {
        Some(m) => {
            let mut v = m.cell_types();
            v.sort();
            v
        }
        None => vec!["".into()],
    };
    // The tracks' files are independent: read them together, then build
    // the panels in order.
    let read: Vec<_> = {
        use rayon::prelude::*;
        groups
            .par_iter()
            .map(|g| {
                read_matrix_positions_grouped(
                    &g.files,
                    selector,
                    &args.signal,
                    membership,
                    &args.top_modality,
                    false,
                    None,
                )
            })
            .collect()
    };
    // The figure is named by every track's genes and chromosomes.
    let mut genes: FxHashMap<Box<str>, usize> = FxHashMap::default();
    let mut chrs: Vec<Box<str>> = Vec::new();
    let mut panels = Vec::with_capacity(groups.len() * celltypes.len());
    for (track, (g, grouped)) in groups.iter().zip(read).enumerate() {
        let mut grouped = grouped?;
        for gene in grouped.genes {
            *genes.entry(gene).or_insert(0) += 1;
        }
        chrs.append(&mut grouped.chrs);
        for (hue, ct) in celltypes.iter().enumerate() {
            let label: Box<str> = match (groups.len(), ct.is_empty()) {
                (1, _) => ct.clone(),
                (_, true) => g.label.clone(),
                _ => format!("{}: {ct}", g.label).into(),
            };
            panels.push(FigurePanel {
                label,
                track,
                celltype: ct.clone(),
                hue,
                sites: grouped.by_group.remove(ct).unwrap_or_default(),
            });
        }
    }
    anyhow::ensure!(
        !genes.is_empty(),
        "no rows matching {} in {} file(s)",
        selector.describe(),
        args.data_files.len()
    );
    let (gene, chr) = (summarize_genes(&genes), summarize_chr(&chrs));
    Ok(FigureSites { gene, chr, panels })
}

/// The panels to draw: each with its track's read depth of its cell type
/// from `depth` (by track; none drawn when absent), its sites on a log
/// scale with `is_log`.
pub(super) fn panel_data(
    panels: Vec<FigurePanel>,
    depth: &[FxHashMap<Box<str>, Vec<f64>>],
    num_bins: usize,
    is_log: bool,
) -> Vec<PanelData> {
    panels
        .into_iter()
        .map(|panel| {
            let epi_sites = if is_log {
                panel
                    .sites
                    .into_iter()
                    .map(|(p, v)| (p, (1.0 + v).log10()))
                    .collect()
            } else {
                panel.sites
            };
            let depth_bins = depth
                .get(panel.track)
                .and_then(|d| d.get(&panel.celltype))
                .cloned()
                .unwrap_or_else(|| vec![0.0; num_bins]);
            PanelData {
                celltype: panel.label,
                hue: panel.hue,
                epi_sites,
                depth_bins,
            }
        })
        .collect()
}

/// Filesystem-safe slug for a default figure output prefix.
pub(super) fn slug(s: &str) -> String {
    let out: String = s
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    if out.is_empty() {
        "miami".to_string()
    } else {
        out
    }
}
