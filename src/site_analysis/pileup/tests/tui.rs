use super::*;
use crate::figure::style::{Sheet, StyleArgs};
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

/// A figure's width over its height.
fn aspect((w, h): (f64, f64)) -> f64 {
    w / h
}

#[test]
fn a_saved_window_has_the_views_shape_even_on_a_short_terminal() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    let (cw, ch) = crate::figure::CELL;
    for (w, h) in [(120, 40), (120, 12), (60, 50)] {
        screen(&mut v, w, h);
        let on_screen = aspect((v.drawn.width as f64 * cw, v.drawn.height as f64 * ch));
        let saved = aspect(crate::figure::svg_size(&v.figure()));
        assert!(
            (saved - on_screen).abs() < 1e-6,
            "{w}x{h}: {saved} vs {on_screen}"
        );
    }
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
    assert_eq!(
        v.exit,
        Some(Exit::Locus {
            lo: 950_000,
            hi: 1_150_000,
            at: c
        }),
        "the cursor kept"
    );
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
    press(&mut v, KeyCode::Char('N'));
    assert_eq!(v.cursor, 1_000_000);
    press(&mut v, KeyCode::Char('N'));
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
    assert!(dir.path().join("pile1.pdf").exists() && !dir.path().join("pile1.png").exists());
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
    let mut a = Track::single("wt", "sum", &m, false).with_total(&half);
    // As `faba pileup` flags two alike tracks.
    a.mirror_with_next = true;
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
        text.contains('█') && text.contains('›') && text.contains('‹'),
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
    assert_eq!(
        v.window, v.extent,
        "not zoomed yet: it opens on the whole gene"
    );
    press(&mut v, KeyCode::Char('+'));
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
    let top = crate::tui::ui::compact(tallest_b);
    assert!(svg.matches(&format!(">{top}<")).count() >= 2, "{top}");
}

#[test]
fn like_tracks_share_a_mirrored_row_until_split() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    let mut a = Track::single("wt", "sum", &m, false).with_total(&half);
    // As `faba pileup` flags two alike tracks.
    a.mirror_with_next = true;
    let b = Track::single("mut", "sum", &half, false).with_total(&half);
    let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
    let text = screen(&mut v, 160, 30);
    assert!(
        text.contains("wt above, mut below · methylated reads"),
        "{text}"
    );
    assert!(text.contains("v split"), "{text}");
    assert!(v.figure().contains("wt above, mut below"));
    // Bars grow both ways from the zero line.
    assert!(text.contains('▀') || text.lines().any(|l| l.contains('█')));

    press(&mut v, KeyCode::Char('v'));
    assert_eq!(
        v.rows(),
        vec![Row::Contrast(DIFFERENCE), Row::Track(0), Row::Track(1)]
    );
    assert!(screen(&mut v, 160, 30).contains("v mirror"));
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

#[test]
fn shift_m_switches_modality_only_where_it_can() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    press(&mut v, KeyCode::Char('n'));
    press(&mut v, KeyCode::Char('M'));
    assert_eq!(v.exit, None, "named files: nothing to switch to");
    v.switchable = true;
    assert!(screen(&mut v, 200, 24).contains("m6A/A-to-I"));
    press(&mut v, KeyCode::Char('M'));
    assert_eq!(v.exit, Some(Exit::Modality { at: v.cursor }));
}

#[test]
fn every_footer_key_shows_at_80_columns() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    let a = Track::single("wt", "sum", &m, false).with_total(&half);
    let b = Track::single("mut", "sum", &half, false).with_total(&half);
    let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
    v.switchable = true;
    let screen = screen(&mut v, 80, 30);
    let missing = crate::tui::missing_keys(&v.footer_line(v.rows_hidden()), &screen);
    assert!(missing.is_empty(), "{missing:?}\n{screen}");
}

#[test]
fn enter_stays_esc_goes_back_to_the_genes_and_q_quits() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    press(&mut v, KeyCode::Enter);
    assert!(v.exit.is_none(), "Enter never quits");
    press(&mut v, KeyCode::Esc);
    assert!(matches!(v.exit, Some(Exit::Genes)));
    let mut v = view(&m, &s);
    press(&mut v, KeyCode::Char('q'));
    assert!(matches!(v.exit, Some(Exit::Quit)));
}

#[test]
fn from_cluster_view_esc_names_it_and_no_search_is_offered() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    let said = help_line(&v.footer_pairs(false)).to_string();
    assert!(
        said.contains("g/Esc gene list") && said.contains("/ search"),
        "{said}"
    );
    v.from_cluster = true;
    let said = help_line(&v.footer_pairs(false)).to_string();
    assert!(said.contains("g/Esc cluster-view"), "{said}");
    assert!(!said.contains("/ search"), "{said}");
    // Nor is it bound: the search is cluster-view's.
    press(&mut v, KeyCode::Char('/'));
    assert!(!v.search.active(), "/ opened a search from cluster-view");
}

#[test]
fn question_mark_lists_the_keys_but_types_in_the_search() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    let keys = crate::tui::help::shown(&mut v, 120, 40).expect("? opens the keys");
    for key in [
        "←/→ h/l",
        "PgUp/PgDn",
        "n/N",
        "+/=",
        "0",
        "/",
        "Y",
        "g/Esc",
        "q",
    ] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
    press(&mut v, KeyCode::Char('/'));
    assert!(crate::tui::help::shown(&mut v, 120, 40).is_none());
    assert_eq!(v.search.text(), Some("?"), "? is typed");
}

#[test]
fn how_bars_draw_carries_to_the_next_gene() {
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    screen(&mut v, 100, 24);
    press(&mut v, KeyCode::Char('y'));
    press(&mut v, KeyCode::Char('Y'));
    let kept = v.settings();
    assert!(kept.y_lock.is_some());
    let mut next = view(&m, &s);
    next.apply(&kept);
    assert_eq!(next.y_scale, v.y_scale, "the y scale carries over");
    assert_eq!(next.y_lock, v.y_lock, "and the locked y maximum");
    assert_eq!(next.bins().shared, v.y_lock);
    press(&mut next, KeyCode::Char('Y'));
    assert!(next.y_lock.is_none(), "Y again unlocks");
}

#[test]
fn the_browser_draws_text_and_saves_turned_position_labels() {
    let (m, s) = (positions(), sites());
    let view = View {
        title: "GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: 0,
    };
    let tracks = vec![Track::single("matrix", "sum", &m, false)];
    let mut b = browser(view, tracks, None);
    assert!(b.controls.images().is_none(), "no image mode in pileup");
    let text = screen(&mut b, 100, 24);
    assert!(!text.contains("image/text"), "{text}");
    let _ = s;
    assert!(b.figure().contains("rotate(-90)"), "position labels turned");
}

#[test]
fn only_the_bars_grab_and_any_button_clicks() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind};
    let (m, s) = (positions(), sites());
    let mut v = view(&m, &s);
    press(&mut v, KeyCode::Char('+'));
    screen(&mut v, 100, 30);
    let at = |kind, column, row| MouseEvent {
        kind,
        column,
        row,
        modifiers: KeyModifiers::NONE,
    };
    let (x, before) = (v.chart_x + 40, v.window);
    // A drag from the header row does not pan.
    v.mouse(at(MouseEventKind::Down(MouseButton::Left), x, 0));
    v.mouse(at(MouseEventKind::Up(MouseButton::Left), x - 10, 0));
    assert_eq!(v.window, before);
    // A right click on the bars moves the cursor.
    let y = v.chart_rows.0 + 2;
    let cursor = v.cursor;
    v.mouse(at(MouseEventKind::Down(MouseButton::Right), x, y));
    v.mouse(at(MouseEventKind::Up(MouseButton::Right), x, y));
    assert_ne!(v.cursor, cursor);
    // The wheel at the start of the span changes nothing, and says so.
    press(&mut v, KeyCode::Char('0'));
    assert!(!v.mouse(at(MouseEventKind::ScrollUp, x, y)));
}

#[test]
fn paired_tracks_mirror_in_pairs_and_the_margin_keeps_its_own_scale() {
    // All cells (the sum) and one group, each converted then unconverted.
    let all_on = vec![(1_010_000, 30.0), (1_050_000, 40.0)];
    let all_off = vec![(1_010_000, 10.0), (1_050_000, 20.0)];
    let one_on = vec![(1_010_000, 3.0), (1_050_000, 4.0)];
    let one_off = vec![(1_010_000, 1.0), (1_050_000, 2.0)];
    let view = View {
        title: "GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: 0,
    };
    let tracks = vec![
        Track::single("all: methylated", "reads", &all_on, false)
            .on_own_scale()
            .converted_above(),
        Track::single("all: unmethylated", "reads", &all_off, false).on_own_scale(),
        Track::single("C1: methylated", "reads", &one_on, false).converted_above(),
        Track::single("C1: unmethylated", "reads", &one_off, false),
    ];
    let mut b = browser(view, tracks, None);
    assert_eq!(b.rows(), vec![Row::Mirror(0, 1), Row::Mirror(2, 3)]);
    let bins = b.bins();
    // The shared scale is the group's, never the margin's.
    assert_eq!(bins.shared, Some(4.0));
    assert_eq!(b.row_max(Row::Mirror(0, 1), &bins), None, "its own");
    assert_eq!(b.row_max(Row::Mirror(2, 3), &bins), Some(4.0));
    let accent: Vec<bool> = bins.shown.iter().map(|s| s.accent).collect();
    assert_eq!(accent, vec![true, false, true, false], "converted above");
    // Split (`v`), each track a row, the margin still on its own scale.
    press(&mut b, KeyCode::Char('v'));
    assert_eq!(b.rows(), (0..4).map(Row::Track).collect::<Vec<_>>());
    let bins = b.bins();
    assert_eq!(b.row_max(Row::Track(1), &bins), None);
    assert_eq!(b.row_max(Row::Track(3), &bins), Some(4.0));
    // Locked (`Y`), the lock is the group's and the margin stays its own.
    press(&mut b, KeyCode::Char('Y'));
    let bins = b.bins();
    assert_eq!(b.row_max(Row::Track(0), &bins), None);
    screen(&mut b, 120, 40);
}

#[test]
fn a_track_flagged_to_mirror_draws_with_the_next_its_accent_on_top() {
    let m = positions();
    let tracks = vec![
        Track::single("alone", "count", &m, false),
        Track::single("up", "sum", &m, false).converted_above(),
        Track::single("down", "sum", &m, false),
    ];
    let mut v = PileupView::new("GENE1", "chr1", tracks, EXTENT);
    assert_eq!(v.rows(), vec![Row::Track(0), Row::Mirror(1, 2)]);
    let accent: Vec<bool> = v.bins().shown.iter().map(|s| s.accent).collect();
    assert_eq!(
        accent,
        vec![false, true, false],
        "the upper half in the accent"
    );
    let text = screen(&mut v, 160, 30);
    assert!(text.contains("up above, down below"), "{text}");
    assert!(text.contains("v split"), "{text}");
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(v.rows(), (0..3).map(Row::Track).collect::<Vec<_>>());

    // Tracks are taken as given: the first two, though alike, stay apart.
    let tracks = vec![
        Track::single("a", "sum", &m, false),
        Track::single("b", "sum", &m, false).converted_above(),
        Track::single("c", "sum", &m, false),
    ];
    let v = PileupView::new("GENE1", "chr1", tracks, EXTENT);
    assert_eq!(v.rows(), vec![Row::Track(0), Row::Mirror(1, 2)]);
}

/// The last track has none after it to share a row with: its flag is
/// dropped, so a view of it alone has nothing for `v` to split.
#[test]
fn the_last_track_shares_no_row() {
    let m = positions();
    let tracks = vec![
        Track::single("a", "sum", &m, false),
        Track::single("b", "sum", &m, false).converted_above(),
    ];
    let v = PileupView::new("GENE1", "chr1", tracks, EXTENT);
    assert_eq!(v.rows(), vec![Row::Track(0), Row::Track(1)]);
    assert!(!v.mirrorable());
}

/// Labels for a margin and `n` groups, each converted then unconverted.
fn group_labels(n: usize) -> Vec<String> {
    std::iter::once("all cells".to_string())
        .chain((1..=n).map(|k| format!("CT{k}")))
        .flat_map(|g| [format!("{g}: methylated"), format!("{g}: unmethylated")])
        .collect()
}

/// A mirrored row per pair of `labels`, the first on its own scale.
fn group_tracks<'a>(labels: &'a [String], m: &'a [(i64, f64)]) -> Vec<Track<'a>> {
    labels
        .chunks(2)
        .enumerate()
        .flat_map(|(k, pair)| {
            let up = Track::single(&pair[0], "reads", m, false).converted_above();
            let down = Track::single(&pair[1], "reads", m, false);
            match k {
                0 => [up.on_own_scale(), down.on_own_scale()],
                _ => [up, down],
            }
        })
        .collect()
}

/// Each mirrored row's title line and the lines of its panel below it.
fn panels(s: &str) -> Vec<(String, Vec<String>)> {
    let lines: Vec<&str> = s.lines().collect();
    let tops: Vec<usize> = (0..lines.len())
        .filter(|&y| lines[y].contains(" above, "))
        .collect();
    tops.iter()
        .enumerate()
        .map(|(k, &y)| {
            let end = tops.get(k + 1).copied().unwrap_or(lines.len());
            let body = lines[y + 1..end].iter().map(|l| l.to_string()).collect();
            (lines[y].to_string(), body)
        })
        .collect()
}

#[test]
fn many_mirrored_rows_each_keep_room_to_draw_and_the_rest_scroll() {
    let m = positions();
    let labels = group_labels(12);
    let mut v = PileupView::new("GENE1", "chr1", group_tracks(&labels, &m), EXTENT);
    let s = screen(&mut v, 200, 50);
    let shown = panels(&s);
    assert!((2..13).contains(&shown.len()), "{s}");
    for (title, body) in &shown {
        let bars = body.iter().any(|l| l.contains('█') || l.contains('▄'));
        let axis = body.iter().any(|l| l.contains('└'));
        assert!(bars && axis, "{title} draws nothing: {s}");
    }
    assert!(shown[0].0.contains("all cells"), "the margin on top: {s}");
    let hidden = 13 - shown.len();
    let last = &shown[shown.len() - 1].0;
    assert!(last.contains(&format!("{hidden}›")), "{last}");
    assert!(!s.contains('‹'), "none above yet: {s}");
    assert!(s.contains("[/]"), "the footer names the keys: {s}");
    // `]` scrolls the groups, the margin kept on top.
    press(&mut v, KeyCode::Char(']'));
    let s = screen(&mut v, 200, 50);
    let after = panels(&s);
    assert!(after[0].0.contains("all cells"), "{s}");
    assert!(after[1].0.contains("‹1 CT2"), "{s}");
    let last = &after[after.len() - 1].0;
    assert!(last.contains(&format!("{}›", hidden - 1)), "{last}");
    // No further than the last group just filling the rows.
    for _ in 0..20 {
        press(&mut v, KeyCode::Char(']'));
    }
    let s = screen(&mut v, 200, 50);
    let end = panels(&s);
    assert_eq!(end.len(), shown.len(), "{s}");
    assert!(
        end[end.len() - 1].0.contains("CT12") && !s.contains('›'),
        "{s}"
    );
    press(&mut v, KeyCode::Char('['));
    assert!(screen(&mut v, 200, 50).contains("1›"));
}

#[test]
fn a_short_track_list_shows_every_row_without_scrolling() {
    let m = positions();
    let labels = group_labels(1);
    let mut v = PileupView::new("GENE1", "chr1", group_tracks(&labels, &m), EXTENT);
    let s = screen(&mut v, 200, 50);
    let shown = panels(&s);
    assert_eq!(shown.len(), 2, "{s}");
    assert!(
        !s.contains('‹') && !s.contains('›') && !s.contains("[/]"),
        "{s}"
    );
    // The rows share the height as before: as many lines each.
    assert!(shown[0].1.len().abs_diff(shown[1].1.len()) <= 1, "{s}");
    press(&mut v, KeyCode::Char(']'));
    assert_eq!(screen(&mut v, 200, 50), s, "nothing to scroll");
}

#[test]
fn the_view_opens_scrolled_to_the_row_asked_for() {
    let m = positions();
    let labels = group_labels(12);
    let view = View {
        title: "GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: 8,
    };
    let mut b = browser(view, group_tracks(&labels, &m), None);
    let s = screen(&mut b, 200, 50);
    let shown = panels(&s);
    assert!(shown[1].0.contains("‹3 CT4"), "{s}");
    // Split, a row a track, it still starts at the same group.
    press(&mut b, KeyCode::Char('v'));
    let s = screen(&mut b, 200, 50);
    let titles: Vec<&str> = s.lines().filter(|l| l.starts_with('╭')).collect();
    assert!(titles[0].contains("all cells: methylated ·"), "{s}");
    assert!(titles[1].contains("‹7 CT4: methylated ·"), "{s}");
}

/// A browser over a margin and 12 groups, opened at track `from`.
fn groups_from<'a>(labels: &'a [String], m: &'a [(i64, f64)], from: usize) -> PileupView<'a> {
    let view = View {
        title: "GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: from,
    };
    browser(view, group_tracks(labels, m), None)
}

#[test]
fn a_view_opened_split_starts_at_the_group_asked_for() {
    let (m, labels) = (positions(), group_labels(12));
    let mut b = groups_from(&labels, &m, 8);
    // Split by the gene before, as it opens.
    b.apply(&Settings {
        mirror: Some(false),
        ..Settings::default()
    });
    let s = screen(&mut b, 200, 50);
    let titles: Vec<&str> = s.lines().filter(|l| l.starts_with('╭')).collect();
    assert!(titles[1].contains("CT4: methylated ·"), "{s}");
    // Mirrored again, the same group's row.
    b.apply(&Settings {
        mirror: Some(true),
        ..Settings::default()
    });
    let s = screen(&mut b, 200, 50);
    assert!(panels(&s)[1].0.contains("‹3 CT4"), "{s}");
}

#[test]
fn the_last_rows_are_reached_at_every_height() {
    let (m, labels) = (positions(), group_labels(12));
    // Widths at which naming `[/]` takes the footer another line.
    for w in [60, 85, 160] {
        for h in 20..60 {
            let mut v = PileupView::new("GENE1", "chr1", group_tracks(&labels, &m), EXTENT);
            screen(&mut v, w, h);
            for _ in 0..30 {
                press(&mut v, KeyCode::Char(']'));
            }
            screen(&mut v, w, h);
            let rows = v.rows();
            assert_eq!(v.shown_rows(&rows).last(), rows.last(), "{w}x{h}");
        }
    }
}

/// `v` as the check draws it, zoomed in on a site mid-gene, so that
/// every move, zoom and jump moves.
fn mid_gene(mut v: PileupView) -> PileupView {
    screen(&mut v, 160, 50);
    press(&mut v, KeyCode::Char('+'));
    press(&mut v, KeyCode::Char('n'));
    press(&mut v, KeyCode::Char('n'));
    v
}

#[test]
fn the_keys_listed_are_the_keys_bound() {
    use crate::tui::help::assert_lists_its_keys;
    let (m, s) = (positions(), sites());
    assert_lists_its_keys("two tracks", || mid_gene(view(&m, &s)));
    assert_lists_its_keys("from cluster-view", || {
        let mut v = view(&m, &s);
        v.from_cluster = true;
        mid_gene(v)
    });
    // Reads to show, two alike tracks to compare and mirror, and another
    // modality to switch to.
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v * 2.0)).collect();
    assert_lists_its_keys("reads, compared, switchable", || {
        let mut a = Track::single("wt", "sum", &m, false).with_total(&half);
        a.mirror_with_next = true;
        let b = Track::single("mut", "sum", &half, false).with_total(&half);
        let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
        v.switchable = true;
        mid_gene(v)
    });
    // Both channels drawn, so the key and where it goes.
    assert_lists_its_keys("both channels, keyed", || {
        let mut a = Track::single("wt", "sum", &m, false).with_total(&half);
        a.mirror_with_next = true;
        let b = Track::single("mut", "sum", &half, false).with_total(&half);
        let mut v = PileupView::new("GENE1", "chr1", vec![a, b], EXTENT);
        v.show = Show::Both;
        mid_gene(v)
    });
    // Rows hidden, scrolled down one so `[` and `]` both scroll.
    let labels = group_labels(12);
    assert_lists_its_keys("rows hidden", || {
        let tracks = group_tracks(&labels, &m);
        let mut v = mid_gene(PileupView::new("GENE1", "chr1", tracks, EXTENT));
        press(&mut v, KeyCode::Char(']'));
        v
    });
}

/// The README's picture of the browser: made-up reads at sites over a
/// two-gene stretch, three batches (two samples mirrored, a control), the
/// gene models under them. Written where `FABA_SCREENSHOTS_DIR` names.
#[test]
#[ignore = "writes the README's screenshot; run by hand"]
fn readme_screenshots_pileup() {
    const LO: i64 = 2_400_000;
    const HI: i64 = 2_436_000;
    // A fixed stream of numbers in [0, 1), so the picture never changes.
    let mut state = 0x2545_f491_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        (state >> 33) as f64 / f64::from(1u32 << 31)
    };
    let sites: Vec<i64> = (0..140)
        .map(|_| LO + 600 + (next().powf(0.6) * (HI - LO - 1200) as f64) as i64)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    // Coverage and the share converted, per site, by batch: the samples
    // alike but for a few sites, the control seldom converted.
    let batch = |gain: f64, next: &mut dyn FnMut() -> f64| {
        let mut converted = Vec::new();
        let mut total = Vec::new();
        for &p in &sites {
            let depth = (8.0 + 90.0 * next() * (p - LO) as f64 / (HI - LO) as f64).round();
            let share = (gain * next() * next()).min(1.0);
            converted.push((p, (depth * share).round()));
            total.push((p, depth));
        }
        (converted, total)
    };
    let (s1, s1_total) = batch(0.9, &mut next);
    let (s2, s2_total) = batch(0.55, &mut next);
    let (ctrl1, ctrl1_total) = batch(0.08, &mut next);
    let name = "methylated / unmethylated (sum)";
    let mut tracks = vec![
        Track::single("s1", name, &s1, false).with_total(&s1_total),
        Track::single("s2", name, &s2, false).with_total(&s2_total),
        Track::single("ctrl1", name, &ctrl1, false).with_total(&ctrl1_total),
    ];
    tracks[0].mirror_with_next = true;
    let mut v = PileupView::new("GENE1", "chr1", tracks, (LO, HI));
    let model = |symbol: &str, lo: i64, hi: i64, forward: bool, exons: &[(i64, i64)]| GeneModel {
        chr: "chr1".into(),
        lo,
        hi,
        forward,
        exons: exons.to_vec(),
        symbol: symbol.into(),
        key: format!("ID1_{symbol}").into(),
    };
    v.genes = vec![
        model(
            "GENE1",
            LO + 400,
            HI - 900,
            true,
            &[
                (LO + 400, LO + 1_300),
                (LO + 7_800, LO + 8_100),
                (LO + 14_000, LO + 14_250),
                (LO + 21_500, LO + 21_900),
                (LO + 27_000, HI - 900),
            ],
        ),
        model(
            "GENE2",
            HI - 6_000,
            HI + 9_000,
            false,
            &[(HI - 6_000, HI - 3_800), (HI + 2_000, HI + 2_600)],
        ),
    ];
    v.keys = None;
    let mut term = Terminal::new(TestBackend::new(160, 44)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    crate::tui::screenshot::save("pileup", term.backend().buffer());
}

/// A WT and a MUT track of GENE1 (its row name a gene id and symbol), each
/// over a total, drawn as `faba pileup --track` draws them: one mirrored
/// row over the gene model, both channels, no contrast row.
fn paper_view<'a>(
    wt: &'a [(i64, f64)],
    mutant: &'a [(i64, f64)],
    total: &'a [(i64, f64)],
) -> PileupView<'a> {
    let mut a = Track::single("WT", "sum", wt, false).with_total(total);
    a.mirror_with_next = true;
    let b = Track::single("MUT", "sum", mutant, false).with_total(total);
    let mut v = PileupView::new("ID1_GENE1", "chr1", vec![a, b], EXTENT);
    v.genes = vec![gene("GENE1", 1_000_000, 1_100_000, true)];
    v.contrast = None;
    v.show = Show::Both;
    v
}

/// The figure `--fig-width w --fig-height h` (inches) asks for, in `look`.
fn paper_figure(v: &mut PileupView, (w, h): (f64, f64), look: StyleArgs) -> anyhow::Result<String> {
    let size = crate::figure::SizeArgs::new(Some(w), Some(h));
    let sheet = Sheet::new(crate::figure::size::Figure::Pileup, &size, &look)?;
    v.columns = 80;
    v.figure_on(&sheet)
}

const PAPER: StyleArgs = StyleArgs {
    font_size: Some(7.0),
    panel_only: true,
    palette: crate::figure::style::Palette::Colour,
};

/// The values with each doubled, as a total over them.
fn doubled(m: &[(i64, f64)]) -> Vec<(i64, f64)> {
    m.iter().map(|&(p, v)| (p, v * 2.0)).collect()
}

/// A paper's panel, 1 by 2.5 inches at 7 points: exactly that size, every
/// label at 7 points, titled by the symbol (not the row's gene id), the
/// gene model drawn once, WT above and MUT below one zero line on one y
/// scale, its axis titled, and the colours keyed only when asked (`--key`).
#[test]
fn a_paper_panel_is_one_mirror_over_one_gene_model() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v / 2.0)).collect();
    let total = doubled(&m);
    let mut v = paper_view(&m, &half, &total);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    assert_eq!(crate::figure::svg_size(&svg), (72.0, 180.0));
    let labels = crate::figure::labels(&svg);
    assert!(labels.iter().all(|l| l.0 == "7"), "{labels:?}");
    let texts: Vec<&str> = labels.iter().map(|l| l.1.as_str()).collect();
    // The symbol as the title, once: no gene id, no label beside the model.
    let named = texts.iter().filter(|&&t| t == "GENE1").count();
    assert_eq!(named, 1, "{texts:?}");
    assert!(!svg.contains("ID1_"), "{svg}");
    let bold = |l: &&str| l.contains(r#"font-weight="bold""#) && l.contains(">GENE1<");
    assert!(svg.lines().any(|l| bold(&l)), "{svg}");
    // One gene model: one intron line, thin as its band.
    let intron = crate::site_analysis::miami::genemodel::INTRON;
    let line = |l: &&str| l.starts_with("<line") && l.contains(intron);
    assert_eq!(svg.lines().filter(line).count(), 1, "{svg}");
    // One plot: one zero, the one scale's top above and below it.
    assert_eq!(texts.iter().filter(|&&t| t == "0").count(), 1, "{texts:?}");
    let top = crate::tui::ui::compact(v.bins().shared.unwrap());
    let tops = texts.iter().filter(|&&t| t == top).count();
    assert_eq!(tops, 2, "{top}: {texts:?}");
    let y_of = |s: &str| labels.iter().find(|l| l.1 == s).and_then(|l| l.2).unwrap();
    assert!(y_of("WT") < y_of("0") && y_of("0") < y_of("MUT"));
    // Its y axis is titled; a panel alone has no key unless asked for.
    assert!(svg.contains(">reads<"), "{svg}");
    assert!(!texts.contains(&"methylated"), "{texts:?}");
    v.with_key = true;
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let labels = crate::figure::labels(&svg);
    let texts: Vec<&str> = labels.iter().map(|l| l.1.as_str()).collect();
    assert!(texts.contains(&"methylated") && texts.contains(&"unmethylated"));
}

/// The position ticks at 1 inch wide and 7 points do not overlap, nor
/// leave the page; some are drawn.
#[test]
fn ticks_are_thinned_to_fit_a_narrow_panel() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ticks: Vec<_> = crate::figure::text_boxes(&svg)
        .into_iter()
        .filter(|(_, s)| mb_labels(&svg).iter().any(|t| t.0 == *s))
        .collect();
    assert!(ticks.len() >= 2, "{ticks:?}");
    let em = 7.0 / 8.0;
    for (k, ((l, t, r, b), s)) in ticks.iter().enumerate() {
        assert!(*l >= 0.0 && *r <= 72.0 / em, "{s} leaves the page");
        for ((l2, t2, r2, b2), s2) in &ticks[k + 1..] {
            let apart = r <= l2 || r2 <= l || b <= t2 || b2 <= t;
            assert!(apart, "{s} and {s2} overlap");
        }
    }
}

/// A short page is taken when the rows fit it; one they cannot is refused
/// naming a size that is drawn.
#[test]
fn a_small_page_is_taken_when_it_fits() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let svg = paper_figure(&mut v, (2.2, 2.4), StyleArgs::default()).unwrap();
    assert_eq!(crate::figure::svg_size(&svg), (2.2 * 72.0, 2.4 * 72.0));
    let e = paper_figure(&mut v, (1.0, 0.8), PAPER)
        .unwrap_err()
        .to_string();
    let inches = |s: &str| s.parse::<f64>().unwrap();
    let (w, rest) = e.split_once(" by ").unwrap();
    let w = inches(w.rsplit(' ').next().unwrap());
    let h = inches(rest.split(' ').next().unwrap());
    assert!(h > 0.8 && h < 2.0, "{e}");
    assert!(paper_figure(&mut v, (w, h), PAPER).is_ok(), "{e}");
}

/// The titles of the panels on `s`, top to bottom.
fn titles(s: &str) -> Vec<&str> {
    s.lines().filter(|l| l.starts_with('╭')).collect()
}

/// The gene model is on top on screen, as the saved figure draws it; `D`
/// hides the WT vs MUT difference row on screen and in the figure saved
/// (and the figure the CLI draws from the same state), and shows it again;
/// the next gene opens as this one was left.
#[test]
fn genes_on_top_and_d_hides_the_difference_row_in_the_view_and_its_figure() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    v.contrast = Some(DIFFERENCE);
    let s = screen(&mut v, 110, 40);
    let t = titles(&s);
    assert!(t[0].contains("genes"), "{s}");
    assert!(t[1].contains("WT vs MUT"), "{s}");
    assert!(v.figure().contains("WT vs MUT"));
    press(&mut v, KeyCode::Char('D'));
    let s = screen(&mut v, 110, 40);
    assert!(!s.contains("WT vs MUT"), "{s}");
    assert!(titles(&s)[0].contains("genes"), "{s}");
    assert!(!v.figure().contains("WT vs MUT"));
    let size = crate::figure::SizeArgs::new(Some(3.0), Some(4.0));
    let sheet = Sheet::new(
        crate::figure::size::Figure::Pileup,
        &size,
        &StyleArgs::default(),
    );
    let svg = v.figure_on(&sheet.unwrap()).unwrap();
    assert_eq!(svg.matches(">0<").count(), 1, "one plot, the mirror: {svg}");
    // The next gene opens without it too.
    let mut next = paper_view(&m, &m, &total);
    next.contrast = Some(DIFFERENCE);
    next.apply(&v.settings());
    assert!(!screen(&mut next, 110, 40).contains("WT vs MUT"));
    press(&mut v, KeyCode::Char('D'));
    assert!(screen(&mut v, 110, 40).contains("WT vs MUT"));
}

/// `K` steps where a saved figure puts its key (inside, right, alone),
/// which the header says; alone, `s` saves it as a file of its own.
#[test]
fn k_steps_where_the_key_goes() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    assert!(screen(&mut v, 110, 40).contains("key inside"));
    for want in ["key right", "key alone", "key inside"] {
        press(&mut v, KeyCode::Char('K'));
        assert!(screen(&mut v, 110, 40).contains(want), "{want}");
        // Alone, `s` saves the key on a page of its own beside the figure.
        let saved: Vec<String> = v.saved().into_iter().map(|f| f.0).collect();
        let alone = want == "key alone";
        assert_eq!(saved.len(), 1 + usize::from(alone), "{want}: {saved:?}");
        assert_eq!(saved.contains(&"_key".to_string()), alone);
    }
}

/// A figure drawn without a view takes the bars the browser shows the
/// whole gene in on a 110-column terminal: the same bases per bar.
#[test]
fn a_figure_has_the_browsers_bases_per_bar() {
    let m = positions();
    let total = doubled(&m);
    let mut on_screen = paper_view(&m, &m, &total);
    screen(&mut on_screen, 110, 40);
    assert_eq!(on_screen.columns, SCREEN_BARS);
    let view = View {
        title: "ID1_GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: 0,
    };
    let mut a = Track::single("WT", "sum", &m, false).with_total(&total);
    a.mirror_with_next = true;
    let b = Track::single("MUT", "sum", &m, false).with_total(&total);
    let at = (Window::Gene, SCREEN_BARS);
    let v = figure_view(view, vec![a, b], &Settings::default(), at, false);
    assert_eq!(v.bin_width(), on_screen.bin_width());
    let svg = v.figure_on(&Sheet::OWN).unwrap();
    let said = format!(">{} bp per bar<", on_screen.bin_width());
    assert!(svg.contains(&said) || svg.contains(&format!("{} bp per bar<", on_screen.bin_width())));
}

/// A tight page's position ticks are two or three: flat where they fit
/// across, else turned; never more than three.
#[test]
fn a_tight_page_has_two_or_three_ticks_flat_where_they_fit() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let wide = paper_figure(&mut v, (5.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&wide);
    assert!((2..=3).contains(&ticks.len()), "{ticks:?}");
    assert!(ticks.iter().all(|t| !t.1), "flat on a wide page: {ticks:?}");
    let narrow = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&narrow);
    assert!((2..=3).contains(&ticks.len()), "{ticks:?}");
}

/// The y scale the view shows is the figure's: the same top, the bars
/// drawn anew on each.
#[test]
fn the_y_scale_of_the_view_draws_the_figure() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let mut drawn = Vec::new();
    for (scale, name) in [
        (Scale::Linear, "linear"),
        (Scale::Sqrt, "sqrt"),
        (Scale::Log, "log"),
    ] {
        v.y_scale = scale;
        assert!(screen(&mut v, 110, 40).contains(&format!("y {name}")));
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        if scale == Scale::Linear {
            let top = crate::tui::ui::compact(v.bins().shared.unwrap());
            assert!(svg.contains(&format!(">{top}<")), "{name}: {top}");
        }
        assert!(!drawn.contains(&svg), "{name} draws as another scale did");
        drawn.push(svg);
    }
}

/// The y labels of a mirror in `svg`, as written: those right-aligned
/// against its axis (where its zero is), top to bottom.
fn y_labels(svg: &str) -> Vec<String> {
    let labels = crate::figure::labels(svg);
    let x_of = |l: &str| {
        l.split(" x=\"")
            .nth(1)
            .and_then(|r| r.split('"').next())
            .map(String::from)
    };
    let texts = || svg.lines().filter(|l| l.starts_with("<text"));
    let axis = texts().find(|l| l.ends_with(">0</text>")).and_then(x_of);
    let right = |l: &&str| l.contains(r#"text-anchor="end""#) && x_of(l) == axis;
    let ends: Vec<&str> = texts().filter(right).collect();
    let mut out: Vec<(f64, String)> = labels
        .into_iter()
        .filter(|l| ends.iter().any(|e| e.contains(&format!(">{}<", l.1))))
        .filter_map(|l| Some((l.2?, l.1)))
        .collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.into_iter().map(|l| l.1).collect()
}

/// A log10 axis labels powers of ten up each arm of the mirror (10, 100,
/// 1k, ...), a square-root one round counts (1, 2 or 5 times a power of
/// ten), each the same up as down, around one zero.
#[test]
fn log10_and_sqrt_axes_label_round_values_on_both_arms() {
    let m: Vec<(i64, f64)> = positions().iter().map(|&(p, v)| (p, v * 2_000.0)).collect();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let top = v.bins().shared.unwrap();
    assert!(top > 1_000.0, "{top}");
    for scale in [Scale::Log, Scale::Sqrt] {
        v.y_scale = scale;
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        let ys = y_labels(&svg);
        let zero = ys.iter().position(|l| l == "0").unwrap();
        let (up, down) = (&ys[..zero], &ys[zero + 1..]);
        assert!(up.len() >= 2, "{scale:?}: {ys:?}");
        let mirrored: Vec<&String> = down.iter().rev().collect();
        assert_eq!(up.iter().collect::<Vec<_>>(), mirrored, "{scale:?}: {ys:?}");
        let round = |l: &str| {
            let n = l.trim_end_matches(['k', 'M']);
            let lead = n.trim_end_matches('0');
            match scale {
                Scale::Log => l == "10" || l == "100" || lead == "1",
                _ => ["1", "2", "5"].contains(&lead),
            }
        };
        assert!(up.iter().all(|l| round(l)), "{scale:?}: {ys:?}");
    }
}

/// The key goes where asked: under the plots, beside them (the plots
/// narrower), or on a page of its own, `_key`, which a figure of the genes
/// side by side shares; a combined figure draws it once.
#[test]
fn the_key_goes_inside_right_or_alone() {
    let m = positions();
    let total = doubled(&m);
    let look = StyleArgs::default();
    let plot_right = |svg: &str| {
        let boxes = crate::figure::text_boxes(svg);
        let key = boxes.iter().find(|(_, s)| s == "methylated").map(|b| b.0);
        key
    };
    let mut v = paper_view(&m, &m, &total);
    let inside = paper_figure(&mut v, (3.0, 3.0), look).unwrap();
    let (l, t, ..) = plot_right(&inside).unwrap();
    v.key_at = KeyAt::Right;
    let right = paper_figure(&mut v, (3.0, 3.0), look).unwrap();
    let (l2, t2, ..) = plot_right(&right).unwrap();
    assert!(
        l2 > l + 50.0 && t2 < t,
        "beside the plots, not under: {l},{t} vs {l2},{t2}"
    );
    v.key_at = KeyAt::Alone;
    let alone = paper_figure(&mut v, (3.0, 3.0), look).unwrap();
    assert!(plot_right(&alone).is_none());
    let sheet = Sheet::new(
        crate::figure::size::Figure::Pileup,
        &crate::figure::SizeArgs::new(Some(1.0), Some(2.5)),
        &PAPER,
    )
    .unwrap();
    let mut two = vec![paper_view(&m, &m, &total), paper_view(&m, &m, &total)];
    two[1].title = "ID2_GENE2".into();
    two[1].genes = vec![gene("GENE2", 1_000_000, 1_100_000, true)];
    for (key_at, with_key, names) in [
        (KeyAt::Alone, false, vec!["_GENE1", "_GENE2", "_row"]),
        (KeyAt::Alone, true, vec!["_GENE1", "_GENE2", "_row", "_key"]),
        (KeyAt::Inside, true, vec!["_GENE1", "_GENE2", "_row"]),
    ] {
        for v in &mut two {
            v.key_at = key_at;
            v.with_key = with_key;
        }
        let row = Combine {
            cols: 2,
            name: "row",
        };
        let files = export(&mut two, &sheet, Some(row)).unwrap();
        let got: Vec<&str> = files.iter().map(|f| f.0.as_str()).collect();
        assert_eq!(got, names, "{key_at:?}");
        let combined = &files[2].1;
        let keys = combined.matches(">unmethylated<").count();
        let want = usize::from(with_key && key_at != KeyAt::Alone);
        assert_eq!(keys, want, "{key_at:?}: one key for the row");
        // Each panel at the page asked for, side by side.
        let (w, _) = crate::figure::svg_size(combined);
        assert_eq!(w, 144.0);
    }
}

/// A bar narrower than a hairline when the page is small is drawn a
/// hairline wide, centred on its place, so a peak an inch-wide panel
/// holds in one bar still shows; a bar of a few reads, a hairline tall.
#[test]
fn a_bar_is_never_thinner_nor_shorter_than_a_hairline() {
    let m = positions();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    v.columns = 400;
    let size = crate::figure::SizeArgs::new(Some(1.0), Some(2.5));
    let sheet = Sheet::new(crate::figure::size::Figure::Pileup, &size, &PAPER).unwrap();
    let svg = v.figure_on(&sheet).unwrap();
    let widths: Vec<f64> = svg
        .lines()
        .filter(|l| l.starts_with("<rect") && l.contains("#cfcfcf"))
        .filter_map(|l| {
            let at = l.find(" width=\"")? + 8;
            l[at..].split('"').next()?.parse().ok()
        })
        .collect();
    assert!(!widths.is_empty());
    assert!(
        widths.iter().all(|&w| w >= crate::figure::MIN_BAR - 1e-9),
        "{widths:?}"
    );
    // A few reads beside a tall peak show as the browser shows them, a
    // hairline tall, not lost.
    let mut few = positions();
    few.push((1_099_000, 0.001));
    let total = doubled(&few);
    let mut v = paper_view(&few, &few, &total);
    v.columns = 400;
    let svg = v.figure_on(&sheet).unwrap();
    let heights: Vec<f64> = svg
        .lines()
        .filter(|l| l.starts_with("<rect") && l.contains("#cfcfcf"))
        .filter_map(|l| {
            let at = l.find(" height=\"")? + 9;
            l[at..].split('"').next()?.parse().ok()
        })
        .collect();
    let least = heights.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(least >= crate::figure::MIN_BAR - 1e-9, "{least}");
}

/// Genes exported together take one room for their position labels, so
/// their plots line up whatever each gene's coordinates: in their own
/// files and side by side.
#[test]
fn panels_in_a_row_line_their_plots_up() {
    let m = positions();
    let total = doubled(&m);
    let far: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p + 120_000_000, v)).collect();
    let far_total = doubled(&far);
    let mut two = vec![
        paper_view(&m, &m, &total),
        paper_view(&far, &far, &far_total),
    ];
    two[1].extent = (121_000_000, 121_100_000);
    two[1].window = two[1].extent;
    two[1].genes = vec![gene("GENE2", 121_000_000, 121_100_000, true)];
    two[1].title = "ID2_GENE2".into();
    let size = crate::figure::SizeArgs::new(Some(1.0), Some(2.5));
    let sheet = Sheet::new(crate::figure::size::Figure::Pileup, &size, &PAPER).unwrap();
    let row = Combine {
        cols: 2,
        name: "row",
    };
    let files = export(&mut two, &sheet, Some(row)).unwrap();
    // The zero line of each: the x axis drawn across at the mirror's middle.
    let zero = |svg: &str| -> f64 {
        let l = svg
            .lines()
            .find(|l| l.starts_with("<text") && l.ends_with(">0</text>"))
            .unwrap();
        let at = l.find(" y=\"").unwrap() + 4;
        l[at..].split('"').next().unwrap().parse().unwrap()
    };
    assert!(
        (zero(&files[0].1) - zero(&files[1].1)).abs() < 1e-6,
        "the zero lines differ"
    );
    // And so in the figure of both.
    let zeros: Vec<&str> = (files[2].1.lines())
        .filter(|l| l.starts_with("<text") && l.ends_with(">0</text>"))
        .collect();
    assert_eq!(zeros.len(), 2);
    assert_eq!(zero(zeros[0]), zero(zeros[1]));
}

/// On a small page the difference row's y title fits along its plot, as
/// the reads' does: a short name, not the browser's long one run past it.
#[test]
fn a_small_pages_y_titles_fit_their_plots() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v / 2.0)).collect();
    let total = doubled(&m);
    let mut v = paper_view(&m, &half, &total);
    v.contrast = Some(DIFFERENCE);
    let svg = paper_figure(&mut v, (1.5, 3.2), PAPER).unwrap();
    let boxes = crate::figure::text_boxes(&svg);
    let y_of = |s: &str| boxes.iter().find(|b| b.1 == s).map(|b| b.0 .1).unwrap();
    let labels: Vec<&str> = boxes.iter().map(|b| b.1.as_str()).collect();
    let top = labels
        .iter()
        .find(|l| l.starts_with('+') && **l != "+0.0")
        .unwrap();
    let bottom = labels.iter().find(|l| l.starts_with('-')).unwrap();
    let title = boxes.iter().find(|b| b.1.starts_with('Δ')).unwrap();
    let (t, b) = (title.0 .1, title.0 .3);
    assert!(
        t >= y_of(top) - 1e-6 && b <= y_of(bottom) + 8.0,
        "{title:?} past {top} {}..{bottom} {}",
        y_of(top),
        y_of(bottom)
    );
}

/// Sites only in the middle of the gene: WT and MUT over a total, from
/// 1,040,000 to 1,060,000 of the 100 kb gene.
fn middle_sites() -> Vec<(i64, f64)> {
    (0..200)
        .map(|i| (1_040_000 + i * 100, 1.0 + (i % 7) as f64))
        .collect()
}

/// `--window sites` and `0`: the window spans the sites of any track (a
/// read or its total), by default centred on them, their span and a
/// quarter of it either side but 1 kb at least, kept within the gene;
/// with a pad, from the first to the last and the pad either side; a base
/// a bar at least; `0` steps from the whole gene to it and back; a gene
/// with no sites stays whole.
#[test]
fn the_sites_window_spans_the_sites_and_its_pad_and_0_steps_to_it() {
    let m = middle_sites();
    // The total reaches a site past the reads': it counts.
    let mut total = doubled(&m);
    total.push((1_070_000, 1.0));
    let mut v = paper_view(&m, &m, &total);
    v.columns = 80;
    // 30 kb of sites: a quarter of it either side.
    assert_eq!(v.sites_window(), (1_032_500, 1_077_500));
    // 400 bases of sites: 1 kb centred on them.
    let short: Vec<(i64, f64)> = (0..5).map(|i| (1_050_000 + i * 100, 2.0)).collect();
    let mut s = paper_view(&short, &short, &short);
    s.columns = 80;
    assert_eq!(s.sites_window(), (1_049_700, 1_050_700));
    // At the gene's start: kept within it, still 1 kb.
    let first: Vec<(i64, f64)> = (0..3).map(|i| (1_000_100 + i * 100, 2.0)).collect();
    let mut s = paper_view(&first, &first, &first);
    s.columns = 80;
    assert_eq!(s.sites_window(), (1_000_000, 1_001_000));
    v.pad = Some(10);
    assert_eq!(v.sites_window(), (1_039_990, 1_070_010));
    // `0` from the whole gene to the sites, and back.
    assert_eq!(v.window, EXTENT);
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.window, (1_039_990, 1_070_010));
    assert!(v.window.0 <= v.cursor && v.cursor <= v.window.1);
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.window, EXTENT);
    // Zoomed in, `0` shows the whole gene first.
    press(&mut v, KeyCode::Char('+'));
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.window, EXTENT);
    // Past the gene's ends the pad stops there.
    v.pad = Some(80_000);
    assert_eq!(v.sites_window(), EXTENT);
    // One site: a base a bar.
    let one = [(1_050_000, 3.0)];
    let mut v = paper_view(&one, &one, &one);
    v.columns = 80;
    v.pad = Some(0);
    let (lo, hi) = v.sites_window();
    assert!(
        hi - lo >= 80 && lo <= 1_050_000 && 1_050_000 <= hi,
        "{lo}-{hi}"
    );
    // None: the whole gene.
    let mut v = paper_view(&[], &[], &[]);
    v.columns = 80;
    assert_eq!(v.sites_window(), EXTENT);
}

/// The figure `--batch-process` draws at the sites is the browser's
/// window after `0`: the same span, the same bases per bar on a
/// 110-column terminal, and the gene model cut there, its cut ends marked.
#[test]
fn a_figure_at_the_sites_is_the_browsers_window_after_0() {
    let m = middle_sites();
    let total = doubled(&m);
    let mut on_screen = paper_view(&m, &m, &total);
    screen(&mut on_screen, 110, 40);
    press(&mut on_screen, KeyCode::Char('0'));
    assert!(screen(&mut on_screen, 110, 40).contains("(sites)"));
    let view = View {
        title: "ID1_GENE1",
        chr: "chr1",
        extent: EXTENT,
        on: "methylated",
        off: "unmethylated",
        keys: None,
        genes: None,
        cursor: None,
        switchable: false,
        from_cluster: false,
        track_from: 0,
    };
    let mut a = Track::single("WT", "sum", &m, false).with_total(&total);
    a.mirror_with_next = true;
    let b = Track::single("MUT", "sum", &m, false).with_total(&total);
    let at = (Window::Sites, SCREEN_BARS);
    let one_row = Settings {
        difference: Some(false),
        ..Settings::default()
    };
    let mut v = figure_view(view, vec![a, b], &one_row, at, false);
    v.genes = vec![gene("GENE1", 1_000_000, 1_100_000, true)];
    assert_eq!(v.window, on_screen.window);
    assert_eq!(v.bin_width(), on_screen.bin_width());
    let span = (v.window.1 - v.window.0) as u64;
    assert_eq!(v.bin_width() as u64, span.div_ceil(SCREEN_BARS as u64));
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    // Two slashes at either cut end.
    assert_eq!(slashes(&svg), 4, "{svg}");
    // `--window gene` draws it whole, uncut.
    v.set_window(Window::Gene);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    assert_eq!(slashes(&svg), 0, "{svg}");
    assert_eq!(v.window, EXTENT);
}

/// The slanted lines of `svg`: a cut gene's break marks.
fn slashes(svg: &str) -> usize {
    let slanted = |l: &&str| {
        let v = |k: &str| {
            let at = l.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
            l[at..].split('"').next().unwrap().parse::<f64>().unwrap()
        };
        (v("y1") - v("y2")).abs() > 1e-9 && (v("x1") - v("x2")).abs() > 1e-9
    };
    svg.lines()
        .filter(|l| l.starts_with("<line"))
        .filter(slanted)
        .count()
}

/// A tight page's ticks are in Mb at round positions, two or three, with
/// just the decimals that tell them apart, and the chromosome under them
/// names the unit; flat on a wide page.
#[test]
fn a_tight_pages_ticks_are_round_and_in_mb() {
    let m = middle_sites();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    v.columns = 80;
    v.window = (1_039_646, 1_041_355);
    let labels = |t: Vec<(usize, String)>| t.into_iter().map(|t| t.1).collect::<Vec<_>>();
    assert_eq!(labels(v.mb_ticks(2)), ["1.040", "1.041"]);
    assert_eq!(labels(v.mb_ticks(3)), ["1.0400", "1.0405", "1.0410"]);
    v.window = (1_000_000, 1_100_000);
    assert_eq!(labels(v.mb_ticks(3)), ["1.00", "1.05", "1.10"]);
    // Each where its position falls.
    let ticks = v.mb_ticks(3);
    assert_eq!(ticks[1].0, v.edges().col_of(1_050_000));
    let wide = paper_figure(&mut v, (5.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&wide);
    assert!((2..=3).contains(&ticks.len()), "{ticks:?}");
    assert!(ticks.iter().all(|t| !t.1), "flat on a wide page: {ticks:?}");
    assert!(wide.contains(">chr1 (Mb)<"), "{wide}");
    assert!(!wide.contains(",0"), "no full coordinates: {wide}");
}

/// The position labels of a figure in Mb: digits and a point, and
/// whether each is turned.
fn mb_labels(svg: &str) -> Vec<(String, bool)> {
    svg.lines()
        .filter(|l| l.starts_with("<text"))
        .filter_map(|l| {
            let s = &l[l.find('>')? + 1..l.rfind("</text>")?];
            let mb = s.chars().all(|c| c.is_ascii_digit() || c == '.') && s.contains('.');
            mb.then(|| (s.to_string(), l.contains("rotate(-90)")))
        })
        .collect()
}

/// At an inch and 7 points the position labels leave most of the panel to
/// the plot: they take less room under it than full coordinates turned.
#[test]
fn an_inch_wide_panels_ticks_leave_the_plot_its_height() {
    let m = middle_sites();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    v.set_window(Window::Sites);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&svg);
    assert!((2..=3).contains(&ticks.len()), "{ticks:?}");
    // The x axis's title is the lowest text; the bottom of the plot, the
    // MUT name just above it.
    let boxes = crate::figure::text_boxes(&svg);
    let y = |s: &str| boxes.iter().find(|b| b.1 == s).unwrap().0;
    let under = y("chr1 (Mb)").3 - y("MUT").3;
    let height = 180.0 / (7.0 / 8.0);
    assert!(under / height < 0.25, "{under} of {height} under the plot");
}

/// In gray the converted reads are near black over their total in light
/// gray, as the qc panels draw kept over all: drawn so and keyed so.
#[test]
fn in_gray_the_converted_reads_are_near_black_over_a_light_total() {
    let m = positions();
    let half: Vec<(i64, f64)> = m.iter().map(|&(p, v)| (p, v / 2.0)).collect();
    let total = doubled(&m);
    let mut v = paper_view(&m, &half, &total);
    v.with_key = true;
    let gray = StyleArgs {
        palette: crate::figure::style::Palette::Gray,
        ..PAPER
    };
    let front = crate::figure::style::Palette::Gray.roles().front;
    assert_eq!(front, crate::figure::style::Palette::Gray.roles().bar);
    for columns in [80, 400] {
        let svg = {
            let size = crate::figure::SizeArgs::new(Some(1.0), Some(2.5));
            let sheet = Sheet::new(crate::figure::size::Figure::Pileup, &size, &gray).unwrap();
            v.columns = columns;
            v.figure_on(&sheet).unwrap()
        };
        let fills = |c: &str| svg.matches(&format!(r#"fill="{c}""#)).count();
        // Every bar's total behind, its converted reads in front, and the
        // key's two swatches.
        assert!(fills(front) > 2, "{columns}: no near-black bars");
        assert!(
            fills(crate::figure::FAINT) > 2,
            "{columns}: no light totals"
        );
        assert!(!svg.contains("#b3b3b3"), "{columns}: mid gray over light");
    }
}

/// A symbol too long for its panel's middle moves over to stay on the
/// page, and one too long for the panel shrinks to fit it.
#[test]
fn a_long_title_stays_inside_its_panel() {
    let m = middle_sites();
    let total = doubled(&m);
    let em = 7.0 / 8.0;
    for symbol in ["GENE12345", "GENEWWWWWMMMMM12345"] {
        let mut v = paper_view(&m, &m, &total);
        v.genes = vec![gene(symbol, 1_000_000, 1_100_000, true)];
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        let line = (svg.lines())
            .find(|l| l.contains(r#"font-weight="bold""#) && l.contains(&format!(">{symbol}<")))
            .unwrap_or_else(|| panic!("{symbol}: {svg}"));
        let num = |k: &str| {
            let at = line.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
            line[at..]
                .split('"')
                .next()
                .unwrap()
                .parse::<f64>()
                .unwrap()
        };
        let half = crate::figure::bold_width(symbol, num("font-size")) / 2.0;
        let (l, r) = (num("x") - half, num("x") + half);
        assert!(
            l >= 0.0 && r <= 72.0 / em,
            "{symbol}: {l}..{r} past 0..{}",
            72.0 / em
        );
    }
}

/// Where a gene at full coordinates starts: past a hundred Mb, so its
/// position labels are as long as a real chromosome's.
const FAR: i64 = 127_740_000;

/// Sites from `FAR + 1_100` over 400 bases, as a short gene's are.
fn far_sites() -> Vec<(i64, f64)> {
    (0..20)
        .map(|i| (FAR + 1_100 + i * 20, 5.0 + (i % 4) as f64 * 40.0))
        .collect()
}

/// [`paper_view`] of GENE{k}, 5 kb from [`FAR`], its window `window`.
fn far_view<'a>(
    m: &'a [(i64, f64)],
    total: &'a [(i64, f64)],
    k: usize,
    window: (i64, i64),
) -> PileupView<'a> {
    let mut v = paper_view(m, m, total);
    v.title = format!("ID{k}_GENE{k}");
    v.extent = (FAR, FAR + 5_000);
    v.window = window;
    v.cursor = window.0;
    v.genes = vec![gene(&format!("GENE{k}"), FAR, FAR + 5_000, true)];
    v.columns = 80;
    v
}

/// The panels nested in a combined figure, in order.
fn nested(svg: &str) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for l in svg.lines() {
        if l.starts_with("<svg x=") {
            out.push(String::new());
        }
        if let Some(p) = out.last_mut() {
            p.push_str(l);
            p.push('\n');
        }
    }
    out
}

/// A sheet of `w` by `h` inches in `look`.
fn sheet_of((w, h): (f64, f64), look: StyleArgs) -> Sheet {
    let size = crate::figure::SizeArgs::new(Some(w), Some(h));
    Sheet::new(crate::figure::size::Figure::Pileup, &size, &look).unwrap()
}

/// In a row or grid only the first panel of each row titles its y axis;
/// the others give that room to their plots, so at full coordinates on an
/// inch at 7 points their position labels lie flat. Each gene's own file
/// keeps its title.
#[test]
fn only_a_rows_first_panel_titles_its_y_axis_and_the_rest_lie_flat() {
    let m = far_sites();
    let total = doubled(&m);
    let window = (FAR + 900, FAR + 2_100);
    let mut views: Vec<_> = (1..=3).map(|k| far_view(&m, &total, k, window)).collect();
    let sheet = sheet_of((1.0, 2.5), PAPER);
    let grid = Combine {
        cols: 2,
        name: "grid",
    };
    let files = export(&mut views, &sheet, Some(grid)).unwrap();
    for (name, svg) in files.iter().filter(|f| f.0.starts_with("_GENE")) {
        assert!(svg.contains(">reads<"), "{name} keeps its y title");
    }
    let combined = &files.iter().find(|f| f.0 == "_grid").unwrap().1;
    let panels = nested(combined);
    let titled: Vec<bool> = panels.iter().map(|p| p.contains(">reads<")).collect();
    assert_eq!(titled, [true, false, true]);
    let ticks = mb_labels(&panels[1]);
    assert!(ticks.len() >= 2, "{ticks:?}");
    assert!(ticks.iter().all(|t| !t.1), "flat: {ticks:?}");
    assert!(ticks.iter().all(|t| t.0.starts_with("127.7")), "{ticks:?}");
}

/// The y gridlines of `svg`'s mirror: each y label's line, 3 above its
/// baseline, laid out.
fn gridlines(svg: &str) -> Vec<f64> {
    let ys = y_labels(svg);
    crate::figure::text_boxes(svg)
        .into_iter()
        .filter(|(_, s)| ys.contains(s))
        .map(|((_, _, _, b), _)| b - 0.25 * 8.0 - 3.0)
        .collect()
}

/// The bars of `svg` as boxes, left, top, right, bottom, laid out: every
/// rectangle but the page.
fn bar_boxes(svg: &str) -> Vec<(f64, f64, f64, f64)> {
    let num = |l: &str, k: &str| {
        let at = l.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
        l[at..].split('"').next().unwrap().parse::<f64>().unwrap()
    };
    (svg.lines())
        .filter(|l| l.starts_with("<rect") && !l.contains("#ffffff"))
        .map(|l| {
            let (x, y) = (num(l, "x"), num(l, "y"));
            (x, y, x + num(l, "width"), y + num(l, "height"))
        })
        .collect()
}

/// The arms' names keep clear of the y labels' lines and of the bars, on
/// every scale, even where the tallest bars stand at the left and a round
/// label sits at the axis's top.
#[test]
fn the_arms_names_keep_clear_of_the_tick_labels_and_the_bars() {
    // Tall bars under the left corner, to a top of a power of ten.
    let tall = |i: i64| (5..40).contains(&i);
    let m: Vec<(i64, f64)> = (0..60)
        .map(|i| (1_000_000 + i * 50, if tall(i) { 1_000.0 } else { 3.0 }))
        .collect();
    let mut v = paper_view(&m, &m, &m);
    v.window = (1_000_000, 1_008_000);
    for scale in [Scale::Linear, Scale::Sqrt, Scale::Log] {
        v.y_scale = scale;
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        let boxes = crate::figure::text_boxes(&svg);
        let lines = gridlines(&svg);
        assert!(lines.len() >= 3, "{scale:?}: {lines:?}");
        let bars = bar_boxes(&svg);
        for name in ["WT", "MUT"] {
            let (l, t, r, b) = boxes.iter().find(|x| x.1 == name).unwrap().0;
            for g in &lines {
                assert!(
                    *g <= t || *g >= b,
                    "{scale:?}: {name} {t}..{b} on a line at {g}"
                );
            }
            for &(bl, bt, br, bb) in &bars {
                let apart = r <= bl || br <= l || b <= bt || bb <= t;
                assert!(apart, "{scale:?}: {name} on a bar {bl},{bt},{br},{bb}");
            }
        }
    }
}

/// Where tall bars stand at both ends of an arm, its name goes between
/// them, clear of them and of the y labels' lines.
#[test]
fn an_arms_name_goes_between_bars_at_both_its_ends() {
    let ends = |i: i64| !(10..50).contains(&i);
    let m: Vec<(i64, f64)> = (0..60)
        .map(|i| (1_000_000 + i * 50, if ends(i) { 1_000.0 } else { 3.0 }))
        .collect();
    let mut v = paper_view(&m, &m, &m);
    v.window = (1_000_000, 1_003_000);
    for scale in [Scale::Linear, Scale::Sqrt, Scale::Log] {
        v.y_scale = scale;
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        let boxes = crate::figure::text_boxes(&svg);
        let bars = bar_boxes(&svg);
        for name in ["WT", "MUT"] {
            let (l, t, r, b) = boxes.iter().find(|x| x.1 == name).unwrap().0;
            for g in gridlines(&svg) {
                assert!(
                    g <= t || g >= b,
                    "{scale:?}: {name} {t}..{b} on a line at {g}"
                );
            }
            for &(bl, bt, br, bb) in &bars {
                let apart = r <= bl || br <= l || b <= bt || bb <= t;
                assert!(apart, "{scale:?}: {name} on a bar {bl},{bt},{br},{bb}");
            }
        }
    }
}

/// The counts of a y axis's labels, as written (`1.5k` is 1500).
fn counts(labels: &[String]) -> Vec<f64> {
    labels
        .iter()
        .map(|l| match l.strip_suffix('k') {
            Some(k) => k.parse::<f64>().unwrap() * 1e3,
            None => l.parse::<f64>().unwrap(),
        })
        .collect()
}

/// A square-root axis labels counts evenly spread, two or more, a round
/// step apart from zero; a log10 one labels its top when that is a power
/// of ten.
#[test]
fn sqrt_labels_are_evenly_spread_and_a_log_top_of_a_power_of_ten_is_labelled() {
    let tops = [1_100.0, 2_400.0, 26_000.0, 32_000.0, 260_000.0];
    for (top, high) in tops.into_iter().flat_map(|t| [(t, 2.5), (t, 2.0)]) {
        let m: Vec<(i64, f64)> = (0..40)
            .map(|i| (1_000_000 + i * 2_500, if i == 7 { top } else { 10.0 }))
            .collect();
        let mut v = paper_view(&m, &m, &m);
        v.y_scale = Scale::Sqrt;
        let svg = paper_figure(&mut v, (1.0, high), PAPER).unwrap();
        let ys = y_labels(&svg);
        let zero = ys.iter().position(|l| l == "0").unwrap();
        let mut up = counts(&ys[..zero]);
        up.reverse();
        assert!(up.len() >= 2, "{top} at {high}: {ys:?}");
        let step = up[0];
        for (k, c) in up.iter().enumerate() {
            assert!(
                (c - step * (k + 1) as f64).abs() < 1e-6,
                "{top} at {high}: {ys:?}"
            );
        }
    }
    let m: Vec<(i64, f64)> = (0..40)
        .map(|i| (1_000_000 + i * 2_500, if i == 7 { 1_000.0 } else { 3.0 }))
        .collect();
    let mut v = paper_view(&m, &m, &m);
    v.y_scale = Scale::Log;
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ys = y_labels(&svg);
    assert_eq!(ys.first().map(String::as_str), Some("1k"), "{ys:?}");
    assert_eq!(ys.last().map(String::as_str), Some("1k"), "{ys:?}");
}

/// The page keeps a margin at its right: the last panel's title, shrunk
/// to fit, stays clear of the edge.
#[test]
fn a_title_keeps_a_margin_at_the_pages_right() {
    let m = middle_sites();
    let total = doubled(&m);
    let em = 7.0 / 8.0;
    let mut v = paper_view(&m, &m, &total);
    let symbol = "GENEWWWWWMMMMM12345";
    v.genes = vec![gene(symbol, 1_000_000, 1_100_000, true)];
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let line = (svg.lines())
        .find(|l| l.contains(r#"font-weight="bold""#) && l.contains(&format!(">{symbol}<")))
        .unwrap();
    let num = |k: &str| {
        let at = line.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
        line[at..]
            .split('"')
            .next()
            .unwrap()
            .parse::<f64>()
            .unwrap()
    };
    let right = num("x") + crate::figure::bold_width(symbol, num("font-size")) / 2.0;
    assert!(right <= 72.0 / em - 2.0, "{right} of {}", 72.0 / em);
}

/// Where three round ticks crowd and no round step gives exactly two,
/// the first and last of the three lie flat before any are turned.
#[test]
fn two_ticks_from_the_ends_of_three_lie_flat_before_turning() {
    let m = far_sites();
    let total = doubled(&m);
    let mut v = far_view(&m, &total, 1, (FAR + 1_500, FAR + 7_000));
    v.y_untitled = true;
    let labels = |t: Vec<(usize, String)>| t.into_iter().map(|t| t.1).collect::<Vec<_>>();
    assert_eq!(labels(v.mb_ticks(3)), ["127.742", "127.744", "127.746"]);
    assert!(v.mb_ticks(2).len() < 2, "{:?}", v.mb_ticks(2));
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&svg);
    let want = [
        ("127.742".to_string(), false),
        ("127.746".to_string(), false),
    ];
    assert_eq!(ticks, want);
}

/// An untitled panel's plot starts just past the y labels it draws: on a
/// square-root or log axis its round labels, not its top's own count,
/// which it does not write.
#[test]
fn an_untitled_plot_starts_past_the_labels_it_draws() {
    let m: Vec<(i64, f64)> = (0..40)
        .map(|i| (1_000_000 + i * 2_500, if i == 7 { 2_300.0 } else { 10.0 }))
        .collect();
    let mut v = paper_view(&m, &m, &m);
    v.y_untitled = true;
    for scale in [Scale::Sqrt, Scale::Log] {
        v.y_scale = scale;
        let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
        let ys = y_labels(&svg);
        assert!(!ys.contains(&"2.3k".to_string()), "{scale:?}: {ys:?}");
        let lefts = crate::figure::text_boxes(&svg)
            .into_iter()
            .filter(|(_, s)| ys.contains(s))
            .map(|(b, _)| b.0);
        let least = lefts.fold(f64::INFINITY, f64::min);
        assert!(
            least <= 6.0,
            "{scale:?}: the labels start at {least}: {ys:?}"
        );
    }
}

/// A position label's point is reckoned at a point's width, not a
/// digit's: two labels that fit only so lie flat.
#[test]
fn a_position_labels_point_takes_a_points_width() {
    let m = far_sites();
    let total = doubled(&m);
    let mut v = far_view(&m, &total, 1, (FAR + 950, FAR + 2_630));
    v.y_untitled = true;
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    let ticks = mb_labels(&svg);
    let want = [
        ("127.741".to_string(), false),
        ("127.742".to_string(), false),
    ];
    assert_eq!(ticks, want);
}

/// Where bars cover an arm end to end, its name goes where it overlaps
/// them least, over a white halo; a name clear of the bars has none.
#[test]
fn a_name_over_bars_has_a_halo_and_a_clear_one_none() {
    let halo = |svg: &str, name: &str| {
        let l = (svg.lines())
            .find(|l| l.starts_with("<text") && l.ends_with(&format!(">{name}</text>")))
            .unwrap()
            .to_string();
        l.contains(r#"paint-order="stroke""#)
    };
    let full: Vec<(i64, f64)> = (0..80).map(|i| (1_000_000 + i * 50, 1_000.0)).collect();
    let mut v = paper_view(&full, &full, &full);
    v.window = (1_000_000, 1_004_000);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    assert!(halo(&svg, "WT") && halo(&svg, "MUT"), "{svg}");
    let m = middle_sites();
    let total = doubled(&m);
    let mut v = paper_view(&m, &m, &total);
    let svg = paper_figure(&mut v, (1.0, 2.5), PAPER).unwrap();
    assert!(!halo(&svg, "WT") && !halo(&svg, "MUT"), "{svg}");
}
