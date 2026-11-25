//! The faceted Miami figure.

use super::*;

/// The figure the browser saves, of the gene's sites (or the whole gene,
/// `--window gene`): its rows on the page
/// and in the look asked for, a pair of alike tracks (`--track WT=...
/// --track MUT=...`) one above the other's mirror on one scale, under the
/// gene model; written to `{out}.miami.svg` (and PDF, PNG as asked). With
/// several genes (`-q A,B`, `--gene-file`), a figure each,
/// `{out}_{SYMBOL}`, and with `--row` or `--grid` one of them all; the
/// files a stack of those genes saves in the browser.
pub(super) fn run_rows_figure(args: &PileupArgs, groups: &[TrackFiles]) -> anyhow::Result<()> {
    let sheet = crate::figure::style::Sheet::new(
        crate::figure::size::Figure::Pileup,
        &args.size,
        &args.look,
    )?;
    let genes = args.gene_list()?;
    let many = genes.len() > 1;
    anyhow::ensure!(
        !many || args.regions.is_empty(),
        "several genes draw a figure each: give --regions without them"
    );
    let selectors: Vec<Selector> = match many {
        true => (genes.iter())
            .map(|g| Selector::build(std::slice::from_ref(g), &[]))
            .collect::<anyhow::Result<_>>()?,
        false => vec![Selector::build(&genes, &args.regions)?],
    };
    let combine = args.combine(selectors.len())?;
    let loaded = (selectors.iter())
        .map(|s| load(args, groups, s, true, None, None))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let models = match args.annotation() {
        Some(gtf) => {
            let keep = |key: &str| selectors.iter().any(|s| s.matches_gene(key));
            let models = SharedModels::default();
            let _ = models.set(Ok(load_gene_models_where(&gtf, keep)?));
            Some(models)
        }
        None => None,
    };
    let names: Vec<String> = loaded.iter().map(|l| stacked_name(args, l)).collect();
    let settings = args.figure_settings();
    let columns = args.num_bins.unwrap_or(tui::SCREEN_BARS);
    let mut views: Vec<_> = (loaded.iter().zip(&names))
        .map(|(l, name)| {
            let view = view_of(l, models.as_ref(), false);
            let at = (args.figure_window(), columns);
            tui::figure_view(view, tracks_of(args, l, name), &settings, at, args.key)
        })
        .collect();
    let files = tui::export(&mut views, &sheet, combine)?;
    let prefix: String = match &args.out {
        Some(p) => p.to_string(),
        None if many => "pileup".into(),
        None => slug(&query_symbol(&loaded[0].gene)),
    };
    let mut written = Vec::new();
    for (k, (suffix, svg)) in files.iter().enumerate() {
        // One gene's figure keeps its name of old.
        let name = match (k, many || combine.is_some()) {
            (0, false) => format!("{prefix}.miami"),
            _ => format!("{prefix}{suffix}"),
        };
        written.extend(write_figure(args, svg, &name)?);
    }
    info!("wrote {}", written.join(", "));
    Ok(())
}

/// `svg` under `prefix`: the SVG (unless `--format pdf`), the PDF (unless
/// `--no-pdf` or `--format svg`) and with `--png` a PNG; the paths written.
fn write_figure(args: &PileupArgs, svg: &str, prefix: &str) -> anyhow::Result<Vec<String>> {
    let mut written = Vec::new();
    if !matches!(args.format, Some(FigFormat::Pdf)) {
        let path = format!("{prefix}.svg");
        std::fs::write(&path, svg)?;
        written.push(path);
    }
    if !args.no_pdf && !matches!(args.format, Some(FigFormat::Svg)) {
        written.extend(crate::figure::save(svg, prefix)?);
    }
    if args.png {
        let png = format!("{prefix}.png");
        std::fs::write(&png, crate::figure::png(svg, args.dpi as f32 / 72.0)?)?;
        written.push(png);
    }
    Ok(written)
}

impl PileupArgs {
    /// The genes asked for: `-q`'s, then `--gene-file`'s, each once.
    pub(super) fn gene_list(&self) -> anyhow::Result<Vec<Box<str>>> {
        let mut genes = self.genes.clone();
        if let Some(file) = &self.gene_file {
            let text = std::fs::read_to_string(file.as_ref())
                .map_err(|e| anyhow::anyhow!("--gene-file {file}: {e}"))?;
            genes.extend(genes_in(&text));
        }
        let mut seen = std::collections::HashSet::new();
        genes.retain(|g| seen.insert(g.to_ascii_uppercase()));
        Ok(genes)
    }

    /// How the figure draws, as the browser's state: the y scale, both
    /// channels, the pair mirrored, the difference row as asked, the key
    /// where asked.
    pub(super) fn figure_settings(&self) -> tui::Settings {
        tui::Settings {
            y_scale: Some(self.yscale.into()),
            show: Some(crate::site_analysis::show::Show::Both),
            mirror: Some(true),
            difference: Some(self.difference),
            key_at: Some(self.key_at),
            pad: self.pad.map(i64::from),
            ..tui::Settings::default()
        }
    }

    /// What the figure spans: `--window`, else the sites.
    pub(super) fn figure_window(&self) -> tui::Window {
        self.window.unwrap_or(tui::Window::Sites)
    }

    /// The figure of all `n` genes `--row` or `--grid` asks for.
    pub(super) fn combine(&self, n: usize) -> anyhow::Result<Option<tui::Combine>> {
        if self.row {
            return Ok(Some(tui::Combine {
                cols: n,
                name: "row",
            }));
        }
        let Some(grid) = &self.grid else {
            return Ok(None);
        };
        let parsed = grid.split_once(['x', 'X']).and_then(|(c, r)| {
            Some((
                c.trim().parse::<usize>().ok()?,
                r.trim().parse::<usize>().ok()?,
            ))
        });
        let Some((cols, rows)) = parsed.filter(|&(c, r)| c > 0 && r > 0) else {
            anyhow::bail!("--grid {grid}: write it as COLSxROWS, e.g. 4x2");
        };
        anyhow::ensure!(
            cols * rows >= n,
            "--grid {grid} holds {} panels: {n} genes need more",
            cols * rows
        );
        Ok(Some(tui::Combine { cols, name: "grid" }))
    }
}

/// The genes a gene file names: words split at lines, commas and spaces,
/// what follows a `#` on a line left out.
pub(super) fn genes_in(text: &str) -> Vec<Box<str>> {
    text.lines()
        .map(|l| l.split('#').next().unwrap_or_default())
        .flat_map(|l| l.split([',', ' ', '\t']))
        .map(str::trim)
        .filter(|w| !w.is_empty())
        .map(Box::from)
        .collect()
}

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

    let edges = BinEdges::new(lo, hi, args.num_bins.unwrap_or(80));

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
        palette: args.celltype_palette.clone(),
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
