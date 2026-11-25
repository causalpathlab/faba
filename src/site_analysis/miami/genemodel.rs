//! Gene-model middle track for the Miami figure.
//!
//! Loads the selected gene(s) from a GTF/GFF — gene footprint + strand
//! (via `build_gene_map`) and merged exons (via `build_exon_intervals`) —
//! and draws a band of exon rectangles on an intron line with
//! strand-direction chevrons. The x-axis is purely positional (never
//! reversed on the minus strand); only the chevrons flip direction, so
//! the model stays aligned with the matrix/depth tracks above and below.

use super::bin::BinEdges;
use crate::figure::{Anchor, Canvas, INK};
use crate::site_analysis::pileup::Selector;
use genomic_data::gff::{
    build_exon_intervals, build_gene_map, parse_gff, FeatureType, GeneSymbol, GffRecord,
};
use genomic_data::sam::Strand;

/// One gene resolved from the GTF.
#[derive(Clone)]
pub struct GeneModel {
    pub chr: Box<str>,
    /// 0-based half-open gene footprint `[lo, hi)`.
    pub lo: i64,
    pub hi: i64,
    pub forward: bool,
    /// 0-based half-open merged exons, sorted.
    pub exons: Vec<(i64, i64)>,
    pub symbol: Box<str>,
    /// `{gene_id}_{symbol}`, as matrix rows name the gene.
    pub key: Box<str>,
}

/// Load the gene model(s) matching `selector` from the GTF. Matching
/// reuses the same `{gene_id}_{symbol}` key + relaxed canonicalizer the
/// matrix rows use, so a `-q GENE1` query resolves the GTF gene too.
pub fn load_gene_models(gtf: &str, selector: &Selector) -> anyhow::Result<Vec<GeneModel>> {
    load_gene_models_where(gtf, |key| selector.matches_gene(key))
}

/// Every gene model of the GTF whose `{gene_id}_{symbol}` key passes `keep`.
pub fn load_gene_models_where(
    gtf: &str,
    keep: impl Fn(&str) -> bool,
) -> anyhow::Result<Vec<GeneModel>> {
    gene_models_from_records(&read_gene_and_exon_records(gtf)?, keep)
}

/// Every gene model in `records` (at least their `gene` and `exon` lines)
/// whose `{gene_id}_{symbol}` key passes `keep`.
pub fn gene_models_from_records(
    records: &[GffRecord],
    keep: impl Fn(&str) -> bool,
) -> anyhow::Result<Vec<GeneModel>> {
    let gene_map = build_gene_map(records, Some(&FeatureType::Gene))?;
    let exon_map = build_exon_intervals(records);

    let mut out = Vec::new();
    for entry in gene_map.iter() {
        let gene_id = entry.key();
        let rec = entry.value();

        let (symbol, gene_key): (Box<str>, Box<str>) = match &rec.gene_name {
            GeneSymbol::Symbol(s) if !s.is_empty() => {
                (s.clone(), format!("{}_{}", gene_id, s).into())
            }
            _ => {
                let id: Box<str> = gene_id.to_string().into();
                (id.clone(), id)
            }
        };

        if !keep(&gene_key) {
            continue;
        }

        let exons = exon_map
            .get(gene_id)
            .map(|e| e.value().clone())
            .unwrap_or_default();

        out.push(GeneModel {
            chr: rec.seqname.clone(),
            lo: rec.start - 1, // GTF 1-based inclusive -> 0-based half-open
            hi: rec.stop,
            forward: matches!(rec.strand, Strand::Forward),
            exons,
            symbol,
            key: gene_key,
        });
    }
    // Stable order (smallest start first) for deterministic rendering.
    out.sort_by_key(|g| (g.lo, g.hi));
    Ok(out)
}

/// The `gene` and `exon` records of a GTF/GFF.
fn read_gene_and_exon_records(gtf: &str) -> anyhow::Result<Vec<GffRecord>> {
    read_records_of(gtf, &["gene", "Gene", "exon"])
}

/// The records of a GTF/GFF whose feature (column 3) is one of `features`,
/// read line by line: only those lines are split and parsed.
pub fn read_records_of(gtf: &str, features: &[&str]) -> anyhow::Result<Vec<GffRecord>> {
    use std::io::BufRead;
    let reader = legume_numeric::matrix::common_io::open_buf_reader(gtf)
        .map_err(|e| anyhow::anyhow!("opening {gtf}: {e}"))?;
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        if !line
            .split('\t')
            .nth(2)
            .is_some_and(|f| features.contains(&f))
        {
            continue;
        }
        let words = line.split('\t').map(|w| Box::from(w.trim())).collect();
        records.extend(parse_gff(words));
    }
    Ok(records)
}

/// Genomic extent spanned by a set of gene models (0-based half-open),
/// or `None` if empty / no usable coordinates.
pub fn models_extent(models: &[GeneModel]) -> Option<(Box<str>, i64, i64)> {
    let mut it = models.iter();
    let first = it.next()?;
    let mut lo = first.lo;
    let mut hi = first.hi;
    for g in it {
        lo = lo.min(g.lo);
        hi = hi.max(g.hi);
    }
    Some((first.chr.clone(), lo, hi))
}

/// The share of a gene model's band its exons take: thin blocks on the
/// intron line, not a bar across the band. (The models carry merged exons
/// only, no UTRs, so every exon is drawn alike.)
pub const EXON_SHARE: f64 = 1.0 / 3.0;

/// A gene model's exons, and its intron line, chevrons and break marks:
/// lighter than the bars' ink, so the model is a quiet guide beside the
/// data.
pub const EXON: &str = "#9a9a9a";
pub const INTRON: &str = "#b0b0b0";

/// Draw one gene model into `[x_left, x_left + plot_w]`, centered on
/// `y_mid` in a band `band_h` tall: exons as boxes [`EXON_SHARE`] of it on
/// an intron line with strand
/// chevrons, the symbol just left of the gene.
pub fn draw_gene_model(
    c: &mut Canvas,
    g: &GeneModel,
    edges: &BinEdges,
    x_left: f64,
    plot_w: f64,
    y_mid: f64,
    band_h: f64,
) {
    let x0 = draw_gene_body(c, g, edges, (x_left, plot_w), y_mid, band_h);
    let lx = (x0 - 2.0).max(x_left);
    c.text(
        lx,
        y_mid + band_h * 0.35,
        &g.symbol,
        band_h * 0.9,
        Anchor::End,
        INK,
    );
}

/// [`draw_gene_model`] without the symbol, where a title names the gene;
/// where the gene starts across.
pub fn draw_gene_body(
    c: &mut Canvas,
    g: &GeneModel,
    edges: &BinEdges,
    (x_left, plot_w): (f64, f64),
    y_mid: f64,
    band_h: f64,
) -> f64 {
    let x_of = |pos: i64| edges.x_px(pos, x_left as f32, plot_w as f32) as f64;
    let (x0, x1) = (x_of(g.lo), x_of(g.hi));
    let exon_h = band_h * EXON_SHARE;

    // Intron line (gene footprint), then strand chevrons along it, each
    // thinning with the band so a small model is not drawn coarse.
    let stroke = (band_h * 0.06).clamp(0.25, 0.8);
    c.line(x0, y_mid, x1, y_mid, INTRON, stroke);
    let step = band_h * 2.5;
    let head = exon_h * 0.35;
    let mut x = x0 + step * 0.5;
    while x < x1 - 1.0 {
        let (tip, base) = if g.forward {
            (x + head, x - head)
        } else {
            (x - head, x + head)
        };
        let chevron = [(base, y_mid - head), (tip, y_mid), (base, y_mid + head)];
        c.polyline(&chevron, INTRON, stroke * 0.8);
        x += step;
    }

    let top = y_mid - exon_h / 2.0;
    for &(es, ee) in &g.exons {
        if ee <= edges.min_pos || es > edges.max_pos {
            continue;
        }
        let (ex0, ex1) = (x_of(es), x_of(ee));
        c.rect(ex0, top, (ex1 - ex0).max(0.8), exon_h, EXON);
    }
    // A gene running past the window is cut there: a break mark, two
    // short slashes across its line, says so at that end.
    let (left, right) = cut_ends(g, edges);
    let cut = [(left, x0, 1.0), (right, x1, -1.0)];
    for (_, at, inward) in cut.into_iter().filter(|c| c.0) {
        let (rise, lean) = (exon_h * 0.9, exon_h * 0.3);
        for k in [0.0, 1.0] {
            let x = at + inward * (lean + k * 1.6 * lean);
            c.line(
                x - lean,
                y_mid + rise,
                x + lean,
                y_mid - rise,
                INTRON,
                stroke,
            );
        }
    }
    x0
}

/// Whether `g` runs past the window `edges` span, left and right.
pub fn cut_ends(g: &GeneModel, edges: &BinEdges) -> (bool, bool) {
    (g.lo < edges.min_pos, g.hi - 1 > edges.max_pos)
}

/// [`draw_gene_model`] as SVG shapes to place in a larger drawing.
pub fn gene_model_svg(
    g: &GeneModel,
    edges: &BinEdges,
    x_left: f32,
    plot_w: f32,
    y_mid: f32,
    band_h: f32,
) -> String {
    let mut c = Canvas::layer();
    let f = |v: f32| v as f64;
    draw_gene_model(&mut c, g, edges, f(x_left), f(plot_w), f(y_mid), f(band_h));
    c.into_body()
}

#[cfg(test)]
mod tests;
