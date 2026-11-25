use super::*;
use crate::site_analysis::metagene::ScaleFactors;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

/// 5'UTR 4 bins, CDS 10, 3'UTR 6, and an optional ncRNA track of 5.
fn hist(non_coding: bool) -> GeneFeatureHistogram {
    GeneFeatureHistogram {
        reads: None,
        counts: [
            vec![1.0, 2.0, 3.0, 4.0],
            (0..10).map(|i| 10.0 + i as f64).collect(),
            vec![30.0, 20.0, 10.0, 5.0, 2.0, 1.0],
            if non_coding {
                vec![7.0, 0.0, 3.0, 0.0, 1.0]
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
        let b = bars(&h, Track::Coding, merge, Show::Sites);
        let total: f64 = h.counts[..3].iter().flatten().sum();
        assert_eq!(
            b.iter().map(|b| b.count).sum::<f64>(),
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
        let n = bars(&h, Track::Coding, m, Show::Sites).len();
        assert!(n <= width || m == 10, "width {width} merge {m}: {n} bars");
        if m > 1 {
            assert!(
                bars(&h, Track::Coding, m - 1, Show::Sites).len() > width,
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
    press(&mut v, KeyCode::Char('0'));
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

#[test]
fn saves_the_profile() {
    let h = hist(true);
    let mut v = MetageneView::new("sites", &h);
    screen(&mut v, 60, 16);
    let svg = v.figure(crate::figure::style::Sheet::OWN);
    for name in ["5&apos;UTR", "CDS", "3&apos;UTR"] {
        assert!(svg.contains(name), "{name}");
    }
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("mg1");
    press(&mut v, KeyCode::Char('s'));
    for _ in 0..20 {
        press(&mut v, KeyCode::Backspace);
    }
    for ch in prefix.to_str().unwrap().chars() {
        press(&mut v, KeyCode::Char(ch));
    }
    press(&mut v, KeyCode::Enter);
    assert!(dir.path().join("mg1.pdf").exists() && !dir.path().join("mg1.png").exists());
}

#[test]
fn draws_through_the_image_path() {
    let h = hist(true);
    let mut v = MetageneView::new("x", &h);
    v.controls = Controls::new("x").with_picker(ratatui_image::picker::Picker::halfblocks());
    let s = screen(&mut v, 80, 16);
    assert!(s.contains("image"));
    press(&mut v, KeyCode::Right);
    screen(&mut v, 80, 16);
    press(&mut v, KeyCode::Char('i'));
    assert!(screen(&mut v, 80, 16).contains("· text"));
}

/// [`hist`] with each bin's reads: two converted and one unconverted per site.
fn hist_with_reads() -> GeneFeatureHistogram {
    let mut h = hist(false);
    let times = |k: f64| h.counts.clone().map(|r| r.iter().map(|v| v * k).collect());
    h.reads = Some([times(2.0), times(1.0)]);
    h
}

#[test]
fn reads_show_converted_first_and_c_cycles() {
    let h = hist_with_reads();
    let mut v = MetageneView::new("x", &h);
    let text = screen(&mut v, 120, 24);
    assert!(text.contains("converted reads per region"), "{text}");
    let total = |v: &MetageneView| v.bars().iter().map(|b| b.count).sum::<f64>();
    let sites: f64 = h.counts[..3].iter().flatten().sum();
    assert_eq!(total(&v), 2.0 * sites, "converted reads");
    press(&mut v, KeyCode::Char('c'));
    assert_eq!(total(&v), sites, "unconverted reads");
    press(&mut v, KeyCode::Char('c'));
    let bars = v.bars();
    assert_eq!(total(&v), 3.0 * sites, "both: the total");
    assert!(bars.iter().all(|b| b.front == Some(2.0 * b.count / 3.0)));
    press(&mut v, KeyCode::Char('c'));
    assert_eq!(total(&v), sites, "the sites");
    assert!(screen(&mut v, 120, 24).contains("c: converted"));
}

#[test]
fn without_reads_it_counts_sites_and_c_does_nothing() {
    let h = hist(false);
    let mut v = MetageneView::new("x", &h);
    assert_eq!(v.show, Show::Sites);
    press(&mut v, KeyCode::Char('c'));
    assert_eq!(v.show, Show::Sites);
    assert!(!screen(&mut v, 120, 24).contains("c: "));
}

#[test]
fn the_text_plot_names_regions_with_no_ticks_and_a_colour_key() {
    let h = hist_with_reads();
    let mut v = MetageneView::new("x", &h);
    let text = screen(&mut v, 120, 24);
    assert!(!text.contains('┴'), "no tick under a region name:\n{text}");
    for name in ["5'UTR", "CDS", "3'UTR"] {
        assert!(
            text.contains(&format!("█ {name}")),
            "{name} in the key:\n{text}"
        );
    }
}

#[test]
fn every_footer_key_shows_at_80_columns() {
    let h = hist(true);
    let mut v = MetageneView::new("sites", &h);
    let s = screen(&mut v, 80, 24);
    for key in [
        "Home/End first/last",
        "+/- finer/coarser",
        "[/] coarser/finer",
        "0 fit",
        "Tab track",
        "s save",
        "q quit",
        "? keys",
    ] {
        assert!(s.contains(key), "{key}\n{s}");
    }
    assert!(!s.contains("q/Esc"), "Esc goes back, never quits:\n{s}");
    assert!(!s.contains("c reads"), "no reads, no c:\n{s}");
    let h = hist_with_reads();
    let mut v = MetageneView::new("sites", &h);
    let s = screen(&mut v, 80, 24);
    assert!(s.contains("c reads"), "{s}");
}

#[test]
fn esc_does_not_quit_but_q_does() {
    let h = hist(false);
    let mut v = MetageneView::new("sites", &h);
    press(&mut v, KeyCode::Esc);
    assert!(!v.done(), "Esc never quits a top-level view");
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done());
}

#[test]
fn zero_fits_a_does_not_and_end_takes_the_last_bar() {
    let h = hist(false);
    let mut v = MetageneView::new("sites", &h);
    screen(&mut v, 100, 18);
    press(&mut v, KeyCode::Char('-'));
    assert!(v.merge.is_some());
    press(&mut v, KeyCode::Char('a'));
    assert!(v.merge.is_some(), "a is no longer fit");
    press(&mut v, KeyCode::Char('0'));
    assert!(v.merge.is_none());
    press(&mut v, KeyCode::End);
    let bars = v.bars();
    assert_eq!(v.cursor_bar(&bars), bars.len() - 1);
    press(&mut v, KeyCode::Home);
    assert_eq!(v.cursor_bar(&v.bars()), 0);
}

#[test]
fn question_mark_lists_every_key_but_types_in_the_save_prompt() {
    let h = hist_with_reads();
    let mut v = MetageneView::new("sites", &h);
    let text = crate::tui::help::shown(&mut v, 100, 30).expect("? opens the keys");
    for key in ["←/→ h/l", "Home", "End", "c", "y", "+/=", "s", "q"] {
        assert!(text.contains(key), "{key}:\n{text}");
    }
    press(&mut v, KeyCode::Char('s'));
    assert!(
        crate::tui::help::shown(&mut v, 100, 30).is_none(),
        "? types"
    );
}

#[test]
fn enter_stays_and_plus_is_finer() {
    let h = hist(true);
    let mut v = MetageneView::new("sites", &h);
    screen(&mut v, 100, 18);
    press(&mut v, KeyCode::Enter);
    assert!(!v.done(), "Enter never quits");
    press(&mut v, KeyCode::Char('-'));
    let coarser = v.merge;
    press(&mut v, KeyCode::Char('+'));
    assert!(
        v.merge < coarser,
        "+ undoes -, finer: {:?} < {coarser:?}",
        v.merge
    );
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done());
}

/// A paper panel's size is taken as asked, however small, as `faba
/// metagene` reads it: the profile is drawn at that size, its text and
/// margins shrunk with the font, every label on the page; none refused.
#[test]
fn a_panel_sized_page_is_not_refused() {
    let size = crate::figure::SizeArgs::new(Some(1.8), Some(1.4));
    let style = crate::figure::style::StyleArgs::default();
    let sheet =
        crate::figure::style::Sheet::new(crate::figure::size::Figure::Metagene, &size, &style)
            .unwrap();
    let svg = super::figure("TITLE1", &hist(true), sheet);
    assert_eq!(crate::figure::svg_size(&svg), (129.6, 100.8));
    let labels = crate::figure::labels(&svg);
    assert!(labels.iter().any(|(_, t, _)| t.contains("CDS")), "{svg}");
    for (_, text, y) in &labels {
        assert!(
            y.is_none_or(|y| (0.0..=100.8).contains(&y)),
            "{text} at {y:?}"
        );
    }
}

/// At 1.8 by 1.4 inches the profile is exactly that size, every label at
/// the font asked for and on the page; the panel alone keeps its bars,
/// axes and their titles only.
#[test]
fn a_panel_is_drawn_exactly_at_the_size_and_font_asked() {
    use crate::figure::style::{Sheet, StyleArgs};
    let h = hist(true);
    let page = (Some(1.8 * 72.0), Some(1.4 * 72.0));
    let at = |panel_only| {
        let look = StyleArgs {
            font_size: Some(6.0),
            panel_only,
            ..Default::default()
        };
        super::figure("TITLE1", &h, Sheet { page, look })
    };
    for svg in [at(false), at(true)] {
        assert_eq!(crate::figure::svg_size(&svg), (129.6, 100.8));
        let labels = crate::figure::labels(&svg);
        assert!(!labels.is_empty());
        for (size, text, y) in &labels {
            assert_eq!(size, "6", "{text}");
            assert!(
                y.is_none_or(|y| (0.0..=100.8).contains(&y)),
                "{text} at {y:?}"
            );
        }
    }
    let panel = at(true);
    for gone in ["TITLE1", "per region", "per bar"] {
        assert!(!panel.contains(gone), "{gone}\n{panel}");
    }
    assert!(panel.contains("CDS") && panel.contains("metagene position"));
    assert!(at(false).contains("bin(s) per bar"));
}

/// A size asked for is the figure's, whatever the view's shape.
#[test]
fn a_figure_is_laid_out_for_a_size_asked_for() {
    let h = hist(true);
    let svg = super::figure(
        "sites",
        &h,
        crate::figure::style::Sheet::of((Some(252.0), Some(144.0))),
    );
    assert_eq!(crate::figure::svg_size(&svg), (252.0, 144.0));
}

#[test]
fn a_saved_figure_keeps_the_views_shape() {
    let h = hist(true);
    let mut v = MetageneView::new("sites", &h);
    let size = crate::figure::svg_size;
    screen(&mut v, 150, 40);
    let (cw, ch) = crate::figure::CELL;
    assert_eq!(
        size(&v.figure(crate::figure::style::Sheet::OWN)),
        (150.0 * cw, 40.0 * ch)
    );
    screen(&mut v, 90, 50);
    assert_eq!(
        size(&v.figure(crate::figure::style::Sheet::OWN)),
        (90.0 * cw, 50.0 * ch),
        "a taller view, a taller page"
    );
}

/// The profile at its own size, no font asked for, is drawn as it always
/// was.
#[test]
fn the_default_profile_is_unchanged() {
    let h = hist(true);
    let svg = super::figure("sites", &h, crate::figure::style::Sheet::OWN);
    assert_eq!(crate::figure::fingerprint(&svg), 6676292108214670719);
}

/// Each with a bar mid-track at a merge of its own, so every move moves.
#[test]
fn the_keys_listed_are_the_keys_bound() {
    for (what, h) in [
        ("sites and an ncRNA track", hist(true)),
        ("reads", hist_with_reads()),
    ] {
        crate::tui::help::assert_lists_its_keys(what, || {
            let mut v = MetageneView::new("x", &h);
            screen(&mut v, 160, 50);
            press(&mut v, KeyCode::Char('-'));
            press(&mut v, KeyCode::Char('-'));
            press(&mut v, KeyCode::Right);
            v
        });
    }
}
