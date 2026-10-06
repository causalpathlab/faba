//! The faceted Miami figure.

use super::*;

/// Faceted Miami-plot path: stratified matrix epi sites (top), GTF gene
/// model (middle), and BAM read depth (bottom), one panel per cell type.
pub(super) fn run_miami_figure(args: &PileupArgs, selector: &Selector) -> anyhow::Result<()> {
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
    let grouped = read_matrix_positions_grouped(
        &args.data_files,
        selector,
        &args.signal,
        membership.as_ref(),
        &args.top_modality,
        false,
        None,
    )?;
    anyhow::ensure!(
        grouped.matched > 0,
        "no rows matching {} in {} file(s)",
        selector.describe(),
        args.data_files.len()
    );

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
    for positions in grouped.by_group.values() {
        for &(p, _) in positions {
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
    let region_chr: Box<str> = if grouped.chr.as_ref() != "*"
        && grouped.chr.as_ref() != "component"
        && !grouped.chr.is_empty()
    {
        grouped.chr.clone()
    } else if let Some(m) = models.first() {
        m.chr.clone()
    } else {
        grouped.chr.clone()
    };

    // Bottom track: read depth from BAM, stratified by cell type.
    let depth_by_group = if args.bam_files.is_empty() {
        FxHashMap::default()
    } else {
        let region = Bed {
            chr: region_chr.clone(),
            start: lo,
            stop: hi + 1,
        };
        read_depth_binned(
            &args.bam_files,
            &region,
            &edges,
            &args.cell_barcode_tag,
            membership.as_ref(),
        )?
    };

    // Panel order: membership cell types (sorted) or one all-cells panel.
    let celltypes: Vec<Box<str>> = match &membership {
        Some(m) => {
            let mut v = m.cell_types();
            v.sort();
            v
        }
        None => vec!["".into()],
    };

    let is_log = matches!(args.signal, PileupSignal::Log10Sum);
    let mut panels: Vec<PanelData> = Vec::with_capacity(celltypes.len());
    for ct in &celltypes {
        let raw = grouped.by_group.get(ct).cloned().unwrap_or_default();
        let epi_sites = if is_log {
            raw.into_iter()
                .map(|(p, v)| (p, (1.0 + v).log10()))
                .collect()
        } else {
            raw
        };
        let depth_bins = depth_by_group
            .get(ct)
            .cloned()
            .unwrap_or_else(|| vec![0.0; edges.num_bins]);
        panels.push(PanelData {
            celltype: ct.clone(),
            epi_sites,
            depth_bins,
        });
    }

    // Output: SVG only (PNG/PDF previously came from plot-utils).
    if args.format.is_some() || args.png || args.svg || args.no_pdf {
        log::warn!(
            "Miami figure writes SVG only; `--format` / `--png` / `--svg` / `--no-pdf` are ignored"
        );
    }

    let out_prefix: Box<str> = args
        .out
        .clone()
        .unwrap_or_else(|| slug(&grouped.gene).into());

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
        grouped.gene,
        region_chr,
        fmt_thousands(lo),
        fmt_thousands(hi)
    )
    .into();

    let opts = FigOpts {
        out_prefix,
        width_in: args.fig_width,
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
