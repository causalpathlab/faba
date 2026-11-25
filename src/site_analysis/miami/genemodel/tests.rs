use super::*;
use crate::site_analysis::pileup::Selector;

fn write_gtf(forward: bool) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    // Unique per call: two tests both write the `forward=true` fixture, so a
    // pid+strand name collides under parallel execution (one test removes the
    // file the other is mid-read). The counter keeps every path distinct.
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let strand = if forward { "+" } else { "-" };
    // gene 100..550 (1-based) with two exons.
    let gtf = format!(
        "chr1\tT\tgene\t100\t550\t.\t{s}\t.\tgene_id \"ENSG1\"; gene_name \"GENEA\"\n\
             chr1\tT\texon\t100\t200\t.\t{s}\t.\tgene_id \"ENSG1\"; gene_name \"GENEA\"\n\
             chr1\tT\texon\t400\t550\t.\t{s}\t.\tgene_id \"ENSG1\"; gene_name \"GENEA\"\n",
        s = strand
    );
    let p = std::env::temp_dir().join(format!(
        "miami_genemodel_{}_{}.gtf",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::write(&p, gtf).unwrap();
    p
}

#[test]
fn loads_and_matches_by_symbol() {
    let gtf = write_gtf(true);
    let sel = Selector::build(&["GENEA".into()], &[]).unwrap();
    let models = load_gene_models(gtf.to_str().unwrap(), &sel).unwrap();
    std::fs::remove_file(&gtf).ok();
    assert_eq!(models.len(), 1);
    let g = &models[0];
    assert_eq!(g.lo, 99); // 100 (1-based) -> 99 (0-based)
    assert_eq!(g.hi, 550);
    assert!(g.forward);
    assert_eq!(g.exons.len(), 2);
    assert_eq!(g.exons[0], (99, 200));
}

#[test]
fn svg_endpoints_and_arrow_direction() {
    let gtf_fwd = write_gtf(true);
    let sel = Selector::build(&["GENEA".into()], &[]).unwrap();
    let fwd = load_gene_models(gtf_fwd.to_str().unwrap(), &sel).unwrap();
    std::fs::remove_file(&gtf_fwd).ok();
    let edges = BinEdges::new(fwd[0].lo, fwd[0].hi, 80);
    let svg = gene_model_svg(&fwd[0], &edges, 10.0, 80.0, 100.0, 8.0);
    assert!(svg.contains("<rect"));
    assert!(svg.contains("<polyline"));
    assert!(svg.contains("GENEA"));

    // forward chevron tip is to the right of its base (tipx > basex).
    let gtf_rev = write_gtf(false);
    let rev = load_gene_models(gtf_rev.to_str().unwrap(), &sel).unwrap();
    std::fs::remove_file(&gtf_rev).ok();
    assert!(!rev[0].forward);
}

/// Every `attr="v"` value of the elements named `tag` in `svg`.
fn attrs(svg: &str, tag: &str, attr: &str) -> Vec<f64> {
    let key = format!(" {attr}=\"");
    svg.lines()
        .filter(|l| l.trim_start().starts_with(&format!("<{tag} ")))
        .filter_map(|l| {
            let at = l.find(&key)? + key.len();
            l[at..].split('"').next()?.parse().ok()
        })
        .collect()
}

/// The exons take about a third of the band, in a gray lighter than the
/// bars' ink; the intron line and the strand chevrons thin with the band
/// and are lighter still, so the model is a quiet guide above the data.
#[test]
fn exons_take_a_modest_share_of_the_band_and_the_lines_scale_with_it() {
    let g = GeneModel {
        chr: "chr1".into(),
        lo: 0,
        hi: 1000,
        forward: true,
        exons: vec![(0, 100), (400, 500), (900, 1000)],
        symbol: "GENE1".into(),
        key: "ID1_GENE1".into(),
    };
    let edges = BinEdges::new(0, 1000, 80);
    for band in [6.0, 9.0, 12.0] {
        let svg = gene_model_svg(&g, &edges, 0.0, 200.0, 20.0, band as f32);
        let heights = attrs(&svg, "rect", "height");
        assert_eq!(heights.len(), 3, "{svg}");
        for h in heights {
            let share = h / band;
            assert!((0.3..=0.37).contains(&share), "exon {h} of band {band}");
        }
        // Lighter than the bars' ink: each channel above it.
        let grey = |c: &str| u8::from_str_radix(&c[1..3], 16).unwrap();
        for fill in svg.lines().filter(|l| l.starts_with("<rect")) {
            let c = fill.split("fill=\"").nth(1).unwrap();
            assert!(grey(c) > grey(crate::figure::BAR) + 0x30, "{fill}");
        }
        let stroke = |l: &str| l.split("stroke=\"").nth(1).map(grey);
        for l in svg
            .lines()
            .filter(|l| l.starts_with("<line") || l.starts_with("<polyline"))
        {
            assert!(stroke(l).unwrap() >= grey(EXON), "{l}");
        }
        let line = attrs(&svg, "line", "stroke-width");
        assert_eq!(line.len(), 1, "{svg}");
        assert!(
            line[0] <= band * 0.1 + 1e-9,
            "intron {} of band {band}",
            line[0]
        );
        for w in attrs(&svg, "polyline", "stroke-width") {
            assert!(w <= line[0] + 1e-9, "chevron {w} past the intron line");
        }
    }
    // A small band's chevrons are small and few: each no taller than the
    // exons, spaced at least two of their own widths apart.
    let svg = gene_model_svg(&g, &edges, 0.0, 60.0, 20.0, 6.0);
    let chevrons = svg.matches("<polyline").count();
    assert!(
        (1..=5).contains(&chevrons),
        "{chevrons} chevrons over 60 pt"
    );
}

/// A window inside the gene clips its model there: nothing drawn past
/// the plot, and a break mark (two short slashes across the line) at each
/// end the gene runs past; a gene the window holds whole has none.
#[test]
fn a_gene_cut_by_the_window_is_clipped_and_its_cut_ends_marked() {
    let g = GeneModel {
        chr: "chr1".into(),
        lo: 0,
        hi: 1000,
        forward: true,
        exons: vec![(0, 100), (400, 500), (900, 1000)],
        symbol: "GENE1".into(),
        key: "ID1_GENE1".into(),
    };
    let slashes = |svg: &str| {
        svg.lines()
            .filter(|l| l.starts_with("<line"))
            .filter(|l| {
                let v = |k: &str| {
                    let at = l.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
                    l[at..].split('"').next().unwrap().parse::<f64>().unwrap()
                };
                (v("y1") - v("y2")).abs() > 1e-9
            })
            .count()
    };
    let (left, w) = (10.0, 200.0);
    let whole = gene_model_svg(&g, &BinEdges::new(0, 999, 80), left, w, 20.0, 9.0);
    assert_eq!(slashes(&whole), 0, "{whole}");
    for (lo, hi, marks) in [(300, 999, 2), (0, 600, 2), (300, 600, 4)] {
        let svg = gene_model_svg(&g, &BinEdges::new(lo, hi, 80), left, w, 20.0, 9.0);
        assert_eq!(slashes(&svg), marks, "{lo}-{hi}: {svg}");
        for x in attrs(&svg, "rect", "x")
            .into_iter()
            .chain(attrs(&svg, "line", "x1"))
        {
            assert!(
                x >= left as f64 - 3.0 && x <= (left + w) as f64 + 3.0,
                "{x}"
            );
        }
        // Only the exons in the window are drawn.
        assert_eq!(
            attrs(&svg, "rect", "x").len(),
            1 + usize::from(lo == 0) + usize::from(hi == 999)
        );
    }
}
