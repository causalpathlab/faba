use super::*;
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
        Track {
            label: "matrix",
            signal: "sum",
            positions: m,
            log: false,
        },
        Track {
            label: "sites",
            signal: "count",
            positions: s,
            log: false,
        },
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
        let got: f64 = v.tracks[0].bin(&v.edges()).iter().sum();
        assert!((got - want).abs() < 1e-9, "{lo}-{hi}: {got} vs {want}");
    }
    let log = Track {
        label: "m",
        signal: "log10-sum",
        positions: &m,
        log: true,
    };
    let edges = BinEdges::new(1_050_000, 1_050_049, 1);
    assert!((log.bin(&edges)[0] - (1.0f64 + 250.0).log10()).abs() < 1e-9);
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
        vec![Track {
            label: "matrix",
            signal: "sum",
            positions: &m,
            log: false,
        }],
        EXTENT,
    );
    let mut one = one;
    screen(&mut one, 12, 6);
}
