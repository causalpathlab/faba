use super::*;
use crate::site_analysis::metagene::ScaleFactors;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

/// 5'UTR 4 bins, CDS 10, 3'UTR 6, and an optional ncRNA track of 5.
fn hist(non_coding: bool) -> GeneFeatureHistogram {
    GeneFeatureHistogram {
        counts: [
            vec![1, 2, 3, 4],
            (0..10).map(|i| 10 + i).collect(),
            vec![30, 20, 10, 5, 2, 1],
            if non_coding {
                vec![7, 0, 3, 0, 1]
            } else {
                Vec::new()
            },
        ],
        scale: ScaleFactors {
            twice_median: [200, 2000, 1200],
            utr5_sf: 0.1,
            utr3_sf: 0.6,
        },
    }
}

fn press(v: &mut MetageneView, code: KeyCode) {
    v.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn screen(v: &mut MetageneView, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    let buf = term.backend().buffer().clone();
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn merging_keeps_regions_and_totals() {
    let h = hist(false);
    for merge in 1..12 {
        let b = bars(&h, Track::Coding, merge);
        let total: usize = h.counts[..3].iter().flatten().sum();
        assert_eq!(
            b.iter().map(|b| b.count).sum::<usize>(),
            total,
            "merge {merge}"
        );
        assert!(b.iter().all(|b| b.last < h.counts[b.region].len()));
        // Never across a boundary: every region starts a fresh bar.
        for r in [UTR5, CDS, UTR3] {
            assert!(b.iter().any(|b| b.region == r && b.first == 0));
        }
    }
}

#[test]
fn fitting_merge_fits() {
    let h = hist(false);
    for width in [1, 3, 5, 8, 13, 20, 40] {
        let m = fitting_merge(&h, Track::Coding, width);
        let n = bars(&h, Track::Coding, m).len();
        assert!(n <= width || m == 10, "width {width} merge {m}: {n} bars");
        if m > 1 {
            assert!(
                bars(&h, Track::Coding, m - 1).len() > width,
                "not the smallest"
            );
        }
    }
}

#[test]
fn the_cursor_survives_merging() {
    let h = hist(false);
    let mut v = MetageneView::new("x", &h);
    screen(&mut v, 80, 16);
    for _ in 0..6 {
        press(&mut v, KeyCode::Right);
    }
    let (region, bin) = v.cursor;
    assert_eq!(region, CDS);
    press(&mut v, KeyCode::Char(']'));
    press(&mut v, KeyCode::Char(']'));
    let bars = v.bars();
    let b = bars[v.cursor_bar(&bars)];
    assert_eq!(b.region, region);
    assert!((b.first..=b.last).contains(&bin));
    press(&mut v, KeyCode::Char('a'));
    assert_eq!(v.merge, None);
}

#[test]
fn tab_switches_to_the_non_coding_track_only_when_there_is_one() {
    let h = hist(false);
    let mut v = MetageneView::new("x", &h);
    press(&mut v, KeyCode::Tab);
    assert_eq!(v.track, Track::Coding);

    let h = hist(true);
    let mut v = MetageneView::new("x", &h);
    press(&mut v, KeyCode::Tab);
    assert_eq!(v.track, Track::NonCoding);
    assert_eq!(v.cursor, (NCRNA, 0));
    press(&mut v, KeyCode::Tab);
    assert_eq!(v.track, Track::Coding);
}

#[test]
fn renders_region_labels_and_the_cursor_readout() {
    let h = hist(true);
    let mut v = MetageneView::new("sites", &h);
    let s = screen(&mut v, 100, 18);
    for name in ["5'UTR", "CDS", "3'UTR"] {
        assert!(s.contains(name), "{name}\n{s}");
    }
    assert!(s.contains("bin 0"), "{s}");
    assert!(s.contains("1 sites"), "{s}");
    assert!(s.contains("median nt 100/1000/600"), "{s}");
    // Narrow terminals merge to fit, and tiny ones do not panic.
    screen(&mut v, 22, 10);
    assert!(v.shown_merge > 1, "14 columns cannot hold 20 bins");
    screen(&mut v, 8, 4);
}
