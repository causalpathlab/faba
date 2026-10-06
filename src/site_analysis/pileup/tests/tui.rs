use super::*;
use crate::site_analysis::miami::genemodel::GeneModel;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

/// Sites every 97 bp over a 100 kb gene body, plus a dense cluster.
fn positions() -> Vec<(i64, f64)> {
    let mut v: Vec<(i64, f64)> = (0..1000).map(|i| (1_000_000 + i * 97, 1.0)).collect();
    v.extend((0..50).map(|i| (1_050_000 + i, 5.0)));
    v.sort_by_key(|p| p.0);
    v
}

fn sites() -> Vec<(i64, f64)> {
    vec![(1_010_000, 2.0), (1_050_010, 7.5), (1_090_000, 1.0)]
}

const EXTENT: (i64, i64) = (1_000_000, 1_100_000);

fn view<'a>(m: &'a [(i64, f64)], s: &'a [(i64, f64)]) -> PileupView<'a> {
    let tracks = vec![
        Track::single("matrix", "sum", m, false),
        Track::single("sites", "count", s, false),
    ];
    PileupView::new("GENE1", "chr1", tracks, EXTENT)
}

fn press(v: &mut PileupView, code: KeyCode) {
    v.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn screen(v: &mut PileupView, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    let buf = term.backend().buffer().clone();
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn bins_hold_exactly_the_window() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    for _ in 0..3 {
        press(&mut v, KeyCode::Char('+'));
        let (lo, hi) = v.window;
        let want: f64 = m
            .iter()
            .filter(|p| p.0 >= lo && p.0 <= hi)
            .map(|p| p.1)
            .sum();
        let got: f64 = v.tracks[0].bin(&v.edges()).0.iter().sum();
        assert!((got - want).abs() < 1e-9, "{lo}-{hi}: {got} vs {want}");
    }
    let log = Track::single("m", "log10-sum", &m, true);
    let edges = BinEdges::new(1_050_000, 1_050_049, 1);
    assert!((log.bin(&edges).0[0] - (1.0f64 + 250.0).log10()).abs() < 1e-9);
}

#[test]
fn zoom_keeps_the_cursor_and_stops_at_a_base_per_bar() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    press(&mut v, KeyCode::Char('n'));
    press(&mut v, KeyCode::Char('n'));
    let c = v.cursor;
    for _ in 0..20 {
        press(&mut v, KeyCode::Char('+'));
        assert!(v.window.0 <= c && c <= v.window.1, "cursor left the window");
    }
    assert!(v.window.1 - v.window.0 >= v.columns as i64 - 1);
    assert_eq!(v.bin_width(), 1);
    for _ in 0..20 {
        press(&mut v, KeyCode::Char('-'));
    }
    assert_eq!(v.window, EXTENT);
    // Past the whole extent, `-` asks for twice its span around it.
    assert_eq!(v.exit, Some(Exit::Locus(950_000, 1_150_000)));
    v.exit = None;
    press(&mut v, KeyCode::Char('+'));
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.window, EXTENT);
}

#[test]
fn site_jumps_and_moves_stay_in_the_extent() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    assert_eq!(v.cursor, 1_000_000, "starts at the first site");
    press(&mut v, KeyCode::Char('n'));
    assert_eq!(v.cursor, 1_000_097);
    press(&mut v, KeyCode::Char('p'));
    assert_eq!(v.cursor, 1_000_000);
    press(&mut v, KeyCode::Char('p'));
    assert_eq!(v.cursor, 1_000_000, "no site before the first");
    for _ in 0..500 {
        press(&mut v, KeyCode::Right);
    }
    assert_eq!(v.cursor, EXTENT.1);
    // Zoomed in, moving past the window pans it.
    for _ in 0..6 {
        press(&mut v, KeyCode::Char('+'));
    }
    for _ in 0..500 {
        press(&mut v, KeyCode::Left);
        assert!(v.window.0 <= v.cursor && v.cursor <= v.window.1);
    }
}

#[test]
fn renders_both_tracks_with_coordinates() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    let text = screen(&mut v, 110, 26);
    assert!(text.contains("matrix · sum"), "{text}");
    assert!(text.contains("sites · count"), "{text}");
    assert!(text.contains("1,000,000"), "{text}");
    assert!(text.contains("chr1:1,000,000"), "{text}");
    assert!(text.contains("site(s) in bar"), "{text}");
    // One track, and a tiny terminal, must not panic.
    let one = PileupView::new(
        "GENE1",
        "chr1",
        vec![Track::single("matrix", "sum", &m, false)],
        EXTENT,
    );
    let mut one = one;
    screen(&mut one, 12, 6);
    press(&mut one, KeyCode::Char('g'));
    assert_eq!(one.exit, Some(Exit::Genes), "g goes back to the gene list");
}

#[test]
fn saves_the_window() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    press(&mut v, KeyCode::Char('+'));
    let svg = v.figure();
    assert!(svg.contains("matrix · sum") && svg.contains("sites · count"));
    assert!(svg.contains("bp per bar"));
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("pile1.png");
    press(&mut v, KeyCode::Char('s'));
    for _ in 0..40 {
        press(&mut v, KeyCode::Backspace);
    }
    for ch in prefix.to_str().unwrap().chars() {
        press(&mut v, KeyCode::Char(ch));
    }
    press(&mut v, KeyCode::Enter);
    assert!(dir.path().join("pile1.pdf").exists() && dir.path().join("pile1.png").exists());
}

#[test]
fn draws_through_the_image_path() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    v.controls = Controls::new("x").with_picker(ratatui_image::picker::Picker::halfblocks());
    screen(&mut v, 160, 24);
    assert_eq!(v.plots.len(), 2, "one image per track");
    press(&mut v, KeyCode::Char('+'));
    let text = screen(&mut v, 160, 24);
    assert!(text.contains("image/text"));
}

fn type_search(v: &mut PileupView, q: &str) {
    press(v, KeyCode::Char('/'));
    for ch in q.chars() {
        press(v, KeyCode::Char(ch));
    }
    press(v, KeyCode::Enter);
}

#[test]
fn search_moves_within_the_view_or_leaves_it() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    type_search(&mut v, "chr1:1,050,000-1,050,040");
    assert!(v.exit.is_none());
    assert!(v.window.0 <= 1_050_000 && v.window.1 >= 1_050_040);
    assert!(v.window.1 - v.window.0 < 1_000, "zoomed to the window");
    type_search(&mut v, "chr1:1090000");
    assert_eq!(v.cursor, 1_090_000);
    assert!(v.exit.is_none());

    type_search(&mut v, "chr2:5-10");
    assert_eq!(
        v.exit,
        Some(Exit::Search("chr2:5-10".into())),
        "another chromosome"
    );
    let mut v = view(&m, &s);
    type_search(&mut v, "GENE2");
    assert_eq!(v.exit, Some(Exit::Search("GENE2".into())));
    let mut v = view(&m, &s);
    press(&mut v, KeyCode::Char('/'));
    press(&mut v, KeyCode::Char('x'));
    press(&mut v, KeyCode::Esc);
    assert!(
        v.exit.is_none() && !v.search.active(),
        "Esc drops the search"
    );
    v.status = Some("no gene matches X".into());
    assert!(screen(&mut v, 100, 24).contains("no gene matches X"));
}

#[test]
fn stacked_tracks_draw_front_over_total() {
    let m = positions();
    let total: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 3.0)).collect();
    let matrix =
        Track::single("wt", "methylated / unmethylated (sum)", &m, false).with_total(&total);
    let mut v = PileupView::new("GENE1", "chr1", vec![matrix], EXTENT);
    // Converted reads first, in the accent.
    let text = screen(&mut v, 110, 26);
    assert!(text.contains("wt · methylated reads"), "{text}");
    let (front, behind) = v.tracks[0].bin(&v.edges());
    let behind = behind.expect("stacked");
    assert!(front.iter().zip(&behind).all(|(f, b)| f <= b));
    assert_eq!(v.bins().shown[0].values, front);
    assert!(v.figure().contains(crate::figure::ACCENT));
    // `c`: the unconverted reads, total less converted.
    press(&mut v, KeyCode::Char('c'));
    let shown = &v.bins().shown[0];
    let less: Vec<f64> = behind.iter().zip(&front).map(|(b, f)| b - f).collect();
    assert_eq!(shown.values, less);
    assert!(!shown.accent && shown.front.is_none());
    assert!(screen(&mut v, 110, 26).contains("wt · unmethylated reads"));
    // Both: converted in front of the total; the readout gives front/total.
    press(&mut v, KeyCode::Char('c'));
    let shown = &v.bins().shown[0];
    assert_eq!(
        (shown.values.clone(), shown.front.clone()),
        (behind.clone(), Some(front.clone()))
    );
    assert!(v.readout(&v.bins()).to_string().contains('/'));
    // Sites: one per distinct position, in the bar that holds it.
    press(&mut v, KeyCode::Char('c'));
    let sites: f64 = v.bins().shown[0].values.iter().sum();
    assert_eq!(sites as usize, distinct_positions(&m).len());
    assert!(screen(&mut v, 110, 26).contains("wt · sites"));
    press(&mut v, KeyCode::Char('c'));
    assert_eq!(v.show, Show::Converted, "and round again");
}

#[test]
fn contrast_measures() {
    let d = Measure::Difference;
    assert_eq!(d.of((5.0, 10.0), (2.0, 10.0)), Some(30.0));
    assert_eq!(
        d.of((1.0, 10.0), (2.0, 4.0)).map(|v| v.round()),
        Some(-40.0)
    );
    assert_eq!(d.of((1.0, 0.0), (2.0, 4.0)), None, "no reads, no value");
    let f = Measure::Log2Fold;
    let up = f.of((8.0, 10.0), (2.0, 10.0)).unwrap();
    assert!(up > 1.0 && f.of((2.0, 10.0), (8.0, 10.0)).unwrap() == -up);
}

#[test]
fn contrast_row_compares_the_first_two_tracks() {
    let (m, s) = (positions(), sites());
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    let quarter: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 4.0)).collect();
    let a = Track::single("wt", "sum", &m, false).with_total(&half);
    let b = Track::single("mut", "sum", &m, false).with_total(&quarter);
    let mut v = PileupView::new(
        "GENE1",
        "chr1",
        vec![a, b, Track::single("sites", "count", &s, false)],
        EXTENT,
    );
    assert_eq!(
        v.rows(),
        vec![Row::Contrast(DIFFERENCE), Row::Mirror(0, 1), Row::Track(2)]
    );
    let text = screen(&mut v, 110, 30);
    assert!(
        text.contains("wt vs mut · methylated fraction difference"),
        "{text}"
    );
    assert!(text.contains("d difference/fold"), "{text}");
    // 1/2 vs 1/4 methylated wherever there are reads: +25 pp.
    let values = v.contrast_values(DIFFERENCE, &v.bins());
    assert!(values.iter().flatten().all(|&x| (x - 25.0).abs() < 1e-9));
    press(&mut v, KeyCode::Char('d'));
    assert_eq!(v.contrast.map(|c| c.measure), Some(Measure::Log2Fold));
    assert!(screen(&mut v, 110, 30).contains("log2 fold"));
    assert!(v.figure().contains("wt vs mut"));

    let single = PileupView::new(
        "GENE1",
        "chr1",
        vec![Track::single("m", "sum", &m, false)],
        EXTENT,
    );
    assert!(single.contrast.is_none(), "nothing to compare");
}

/// The first two tracks compared by difference, as a view starts.
const DIFFERENCE: Contrast = Contrast {
    a: 0,
    b: 1,
    measure: Measure::Difference,
};

fn gene(symbol: &str, lo: i64, hi: i64, forward: bool) -> GeneModel {
    GeneModel {
        chr: "chr1".into(),
        lo,
        hi,
        forward,
        exons: vec![(lo, lo + 500), (hi - 800, hi)],
        symbol: symbol.into(),
        key: format!("ID1_{symbol}").into(),
    }
}

#[test]
fn genes_row_stacks_overlapping_genes() {
    let (m, s) = (positions(), sites());
    let genes = [
        gene("GENE1", 1_000_000, 1_060_000, true),
        gene("GENE2", 1_040_000, 1_090_000, false),
        gene("GENE3", 1_095_000, 1_120_000, true),
        gene("GENE4", 5_000_000, 5_010_000, true),
    ];
    let mut v = view(&m, &s);
    v.genes = genes.to_vec();
    screen(&mut v, 110, 34);
    let lanes = v.gene_lanes();
    let lane_of = |sym: &str| {
        lanes
            .iter()
            .find(|(_, g)| &*g.symbol == sym)
            .map(|(l, _)| *l)
    };
    assert_eq!(lane_of("GENE1"), Some(0));
    assert_eq!(lane_of("GENE2"), Some(1), "overlaps GENE1");
    assert_eq!(lane_of("GENE4"), None, "outside the window");
    assert_eq!(v.rows().last(), Some(&Row::Genes));
    let text = screen(&mut v, 110, 34);
    assert!(text.contains("GENE1") && text.contains("GENE2"), "{text}");
    assert!(
        text.contains('━') && text.contains('›') && text.contains('‹'),
        "{text}"
    );
    let svg = v.figure();
    assert!(
        svg.contains("GENE3") && svg.contains("<rect"),
        "gene models in the figure"
    );
}

#[test]
fn contrast_and_titles_name_the_channels() {
    let m = positions();
    let total: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    let a = Track::single("wt", "sum", &m, false).with_total(&total);
    let b = Track::single("mut", "sum", &m, false).with_total(&total);
    let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
    v.on = crate::site_analysis::pileup::channel_names("atoi").0;
    assert!(screen(&mut v, 110, 30).contains("converted fraction difference"));
}

#[test]
fn depth_track_shows_the_bin_under_each_column() {
    let ranges = [(1_000_000, 1_050_000, 40.0), (1_050_000, 1_100_000, 10.0)];
    let edges = BinEdges::new(1_000_000, 1_100_000, 10);
    assert_eq!(
        ranges_per_column(&ranges, &edges),
        vec![40.0, 40.0, 40.0, 40.0, 40.0, 10.0, 10.0, 10.0, 10.0, 10.0]
    );
    let gap = [(1_000_000, 1_010_000, 5.0)];
    assert_eq!(ranges_per_column(&gap, &edges)[5], 0.0, "no bin, no depth");

    let m = positions();
    let tracks = vec![
        Track::single("m", "sum", &m, false),
        Track::depth("depth", &ranges),
    ];
    let mut v = PileupView::new("GENE1", "chr1", tracks, EXTENT);
    let text = screen(&mut v, 110, 30);
    assert!(text.contains("depth · reads per depth bin"), "{text}");
    assert!(v.figure().contains("depth"));
}

#[test]
fn genes_arriving_late_show_without_a_key() {
    let m = positions();
    let models = crate::site_analysis::pileup::SharedModels::default();
    let mut v = PileupView::new(
        "ENSG1_GENE1",
        "chr1",
        vec![Track::single("m", "sum", &m, false)],
        EXTENT,
    );
    v.pending_genes = Some(models.clone());
    screen(&mut v, 110, 30);
    assert!(!v.rows().contains(&Row::Genes), "nothing yet");
    assert!(!v.tick(), "no redraw while waiting");
    let _ = models.set(Ok(vec![
        gene("GENE1", 1_000_000, 1_060_000, true),
        gene("GENE2", 1_070_000, 1_090_000, true),
    ]));
    assert!(v.tick(), "arrival redraws without a key");
    let text = screen(&mut v, 110, 30);
    assert!(v.rows().contains(&Row::Genes));
    assert!(
        text.contains("GENE1") && !text.contains("GENE2"),
        "only the opened gene: {text}"
    );
}

#[test]
fn the_genes_model_widens_the_view_to_its_tss_and_tes() {
    let m = positions();
    let models = crate::site_analysis::pileup::SharedModels::default();
    let mut v = PileupView::new(
        "ENSG1_GENE1",
        "chr1",
        vec![Track::single("m", "sum", &m, false)],
        EXTENT,
    );
    v.pending_genes = Some(models.clone());
    let _ = models.set(Ok(vec![
        gene("GENE1", 950_000, 1_150_000, true),
        gene("GENE2", 1_050_000, 1_300_000, true),
    ]));
    v.tick();
    assert_eq!(v.extent, (950_000, 1_149_999), "GENE1 only, TSS to TES");
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.window, v.extent, "zooming all the way out shows it");
}

#[test]
fn read_tracks_share_one_scale() {
    let m = positions();
    let big: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 10.0)).collect();
    let depth = [(1_000_000, 1_100_000, 1e9)];
    let tracks = vec![
        Track::single("a", "sum", &m, false),
        Track::single("b", "sum", &big, false),
        Track::depth("depth", &depth),
    ];
    let mut v = PileupView::new("GENE1", "chr1", tracks, EXTENT);
    screen(&mut v, 110, 30);
    let tallest_b = v.tracks[1]
        .bin(&v.edges())
        .0
        .into_iter()
        .fold(0.0, f64::max);
    assert_eq!(
        v.bins().shared,
        Some(tallest_b),
        "the larger track sets it; depth does not"
    );
    // Both read panels label the same top in the figure.
    let svg = v.figure();
    let top = data_beans::interactive::ui::compact(tallest_b);
    assert!(svg.matches(&format!(">{top}<")).count() >= 2, "{top}");
}

#[test]
fn like_tracks_share_a_mirrored_row_until_split() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    let a = Track::single("wt", "sum", &m, false).with_total(&half);
    let b = Track::single("mut", "sum", &half, false).with_total(&half);
    let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
    let text = screen(&mut v, 160, 30);
    assert!(
        text.contains("wt above, mut below · methylated reads"),
        "{text}"
    );
    assert!(text.contains("m split"), "{text}");
    assert!(v.figure().contains("wt above, mut below"));
    // Bars grow both ways from the zero line.
    assert!(text.contains('▀') || text.lines().any(|l| l.contains('█')));

    press(&mut v, KeyCode::Char('m'));
    assert_eq!(
        v.rows(),
        vec![Row::Contrast(DIFFERENCE), Row::Track(0), Row::Track(1)]
    );
    assert!(screen(&mut v, 160, 30).contains("m mirror"));
}

#[test]
fn unlike_tracks_are_not_mirrored() {
    let (m, s) = (positions(), sites());
    let v = view(&m, &s);
    assert_eq!(v.rows(), vec![Row::Track(0), Row::Track(1)]);
    let depth = [(1_000_000, 1_100_000, 5.0)];
    let v = PileupView::new(
        "GENE1",
        "chr1",
        vec![
            Track::single("wt", "sum", &m, false),
            Track::depth("depth", &depth),
        ],
        EXTENT,
    );
    assert_eq!(v.rows(), vec![Row::Track(0), Row::Track(1)]);
}
