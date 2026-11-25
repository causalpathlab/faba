use super::*;
use crate::data::dna::Dna;
use crate::figure::style::{AxisTitles, Palette, StyleArgs};
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

fn count(a: usize, c: usize, g: usize, t: usize) -> DnaBaseCount {
    let mut x = DnaBaseCount::new();
    x.add(Some(&Dna::A), a);
    x.add(Some(&Dna::C), c);
    x.add(Some(&Dna::G), g);
    x.add(Some(&Dna::T), t);
    x
}

/// Window 2: uniform, A-rich, pure A at the site, C-rich, empty.
fn pwm() -> Vec<DnaBaseCount> {
    vec![
        count(10, 10, 10, 10),
        count(30, 4, 3, 3),
        count(40, 0, 0, 0),
        count(2, 30, 4, 4),
        count(0, 0, 0, 0),
    ]
}

fn press(v: &mut PwmView, code: KeyCode) {
    v.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn screen(v: &mut PwmView, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    let buf = term.backend().buffer().clone();
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn information_content_bounds() {
    let v = PwmView::new("x", &pwm(), 2);
    assert!(
        v.columns[0].bits().abs() < 1e-12,
        "uniform carries no information"
    );
    assert!(
        (v.columns[2].bits() - 2.0).abs() < 1e-12,
        "a pure column carries 2 bits"
    );
    assert_eq!(v.columns[4].bits(), 0.0, "an empty column carries none");
    let f = v.columns[1].freqs();
    assert!((f.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert_eq!(v.columns[1].counts, [30, 4, 3, 3]);
}

#[test]
fn allotted_rows_add_up() {
    for h in [
        [0.5, 0.25, 0.15, 0.1],
        [0.0; 4],
        [1.0 / 3.0; 4],
        [0.9, 0.05, 0.03, 0.02],
    ] {
        for rows in [3.0, 7.5, 10.0, 21.0] {
            let got = allot(&h, rows);
            let want = (h.iter().sum::<f64>() * rows).round() as usize;
            assert_eq!(got.iter().sum::<usize>(), want, "{h:?} {rows}");
        }
    }
}

#[test]
fn keys_move_the_cursor_and_switch_mode() {
    let mut v = PwmView::new("x", &pwm(), 2);
    assert_eq!(v.rel(v.cursor), 0, "starts at the site");
    for _ in 0..10 {
        press(&mut v, KeyCode::Right);
    }
    assert_eq!(v.cursor, 4);
    press(&mut v, KeyCode::Home);
    assert_eq!(v.cursor, 0);
    press(&mut v, KeyCode::Char('0'));
    assert_eq!(v.rel(v.cursor), 0);
    press(&mut v, KeyCode::Char('b'));
    assert_eq!(v.mode, Mode::Freq);
    assert!(!v.done());
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done());
}

#[test]
fn renders_a_logo_with_the_site_labelled() {
    let mut v = PwmView::new("sites", &pwm(), 2);
    let s = screen(&mut v, 100, 20);
    assert!(s.contains('A'), "{s}");
    assert!(s.contains("position +0"), "{s}");
    assert!(s.contains("IC 2.000 bits"), "{s}");
    // In bits the uniform column is empty; in frequency it is full.
    let bits_letters = s.matches(['A', 'C', 'G', 'T']).count();
    press(&mut v, KeyCode::Char('b'));
    let freq_letters = screen(&mut v, 100, 20)
        .matches(['A', 'C', 'G', 'T'])
        .count();
    assert!(freq_letters > bits_letters);
}

#[test]
fn a_wide_window_scrolls_to_the_cursor() {
    let wide: Vec<DnaBaseCount> = (0..201).map(|i| count(i, 1, 1, 1)).collect();
    let mut v = PwmView::new("x", &wide, 100);
    screen(&mut v, 40, 12);
    press(&mut v, KeyCode::End);
    screen(&mut v, 40, 12);
    assert!(v.offset > 0 && v.offset <= v.cursor);
    // A tiny terminal must not panic.
    screen(&mut v, 8, 4);
}

/// A paper panel's size is taken as asked, however many positions: the
/// letters fill what room there is.
#[test]
fn a_logo_of_any_width_is_not_refused() {
    use crate::figure::style::{Sheet, StyleArgs};
    let wide: Vec<DnaBaseCount> = (0..201).map(|i| count(i, 1, 1, 1)).collect();
    let size = crate::figure::SizeArgs::new(Some(1.8), Some(1.4));
    let sheet = Sheet::new(
        crate::figure::size::Figure::Logo,
        &size,
        &StyleArgs::default(),
    )
    .unwrap();
    let svg = super::figure("sites", &wide, 100, sheet, &AxisTitles::default());
    assert_eq!(crate::figure::svg_size(&svg), (129.6, 100.8));
}

/// At 1.8 by 1.4 inches the logo is exactly that size, every label at the
/// font asked for and on the page (the letters fill their boxes); the
/// panel alone has no title.
#[test]
fn a_logo_panel_is_drawn_exactly_at_the_size_and_font_asked() {
    use crate::figure::style::{Sheet, StyleArgs};
    let page = (Some(1.8 * 72.0), Some(1.4 * 72.0));
    let at = |panel_only| {
        let look = StyleArgs {
            font_size: Some(5.5),
            panel_only,
            ..Default::default()
        };
        super::figure(
            "TITLE1",
            &pwm(),
            2,
            Sheet { page, look },
            &AxisTitles::default(),
        )
    };
    for svg in [at(false), at(true)] {
        assert_eq!(crate::figure::svg_size(&svg), (129.6, 100.8));
        for (size, text, y) in crate::figure::labels(&svg) {
            assert_eq!(size, "5.5", "{text}");
            assert!(
                y.is_none_or(|y| (0.0..=100.8).contains(&y)),
                "{text} at {y:?}"
            );
        }
        assert!(svg.contains(">A<"), "the letters are drawn");
    }
    assert!(at(false).contains("TITLE1"));
    let panel = at(true);
    assert!(!panel.contains("TITLE1"), "{panel}");
    // The y axis 56 points in and 10 down at 8-point text, scaled by 5.5/8.
    assert!(
        panel.contains(r#"<g transform="scale(0.6875)">"#),
        "{panel}"
    );
    assert!(
        panel.contains(r#"x1="56.00" y1="10.00" x2="56.00""#),
        "{panel}"
    );
    assert!(panel.contains("position relative to the site"));
}

/// A size asked for is the logo's, whatever the view's shape.
#[test]
fn the_logo_is_laid_out_for_a_size_asked_for() {
    let svg = super::figure(
        "sites",
        &pwm(),
        2,
        crate::figure::style::Sheet::of((Some(504.0), Some(216.0))),
        &AxisTitles::default(),
    );
    assert_eq!(crate::figure::svg_size(&svg), (504.0, 216.0));
}

#[test]
fn saves_the_logo() {
    let mut v = PwmView::new("sites", &pwm(), 2);
    let svg = v.figure(crate::figure::style::Sheet::OWN);
    assert!(svg.contains(">A<") && svg.contains(crate::figure::ACCENT));
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("logo1");
    press(&mut v, KeyCode::Char('s'));
    for _ in 0..20 {
        press(&mut v, KeyCode::Backspace);
    }
    for ch in prefix.to_str().unwrap().chars() {
        press(&mut v, KeyCode::Char(ch));
    }
    press(&mut v, KeyCode::Enter);
    assert!(dir.path().join("logo1.pdf").exists() && !dir.path().join("logo1.png").exists());
    assert!(!v.done());
}

#[test]
fn big_letters_are_drawn_from_bitmaps() {
    assert!(bitmap(0, 3, 0) && bitmap(0, 3, 4), "A's crossbar");
    assert!(!bitmap(3, 6, 0) && bitmap(3, 6, 2), "T's stem");
    let mut v = PwmView::new("x", &pwm(), 2);
    let tall = screen(&mut v, 60, 40);
    assert!(tall.contains('█'), "{tall}");
}

#[test]
fn draws_through_the_image_path() {
    let mut v = PwmView::new("x", &pwm(), 2);
    v.controls = Controls::new("x").with_picker(ratatui_image::picker::Picker::halfblocks());
    screen(&mut v, 60, 20);
    press(&mut v, KeyCode::Char('i'));
    assert!(v.controls.images().is_none());
    let text = screen(&mut v, 60, 20);
    assert!(text.contains("image/text"));
}

#[test]
fn every_footer_key_shows_at_60_columns() {
    let mut v = PwmView::new("x", &pwm(), 2);
    let s = screen(&mut v, 60, 20);
    for key in [
        "←/→ position",
        "Home/End ends",
        "0 site",
        "v bits/frequency",
        "s save",
        "q quit",
        "? keys",
    ] {
        assert!(s.contains(key), "{key}\n{s}");
    }
    assert!(!s.contains("b bits"), "b is an alias, not shown:\n{s}");
}

#[test]
fn enter_and_esc_stay_q_leaves() {
    let mut v = PwmView::new("x", &pwm(), 2);
    screen(&mut v, 60, 20);
    press(&mut v, KeyCode::Enter);
    assert!(!v.done(), "Enter never quits");
    press(&mut v, KeyCode::Esc);
    assert!(!v.done(), "Esc never quits a top-level view");
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done());
}

#[test]
fn v_switches_bits_and_frequency_as_b_still_does() {
    let mut v = PwmView::new("x", &pwm(), 2);
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(v.mode, Mode::Freq);
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(v.mode, Mode::Bits);
    press(&mut v, KeyCode::Char('b'));
    assert_eq!(v.mode, Mode::Freq);
}

#[test]
fn question_mark_lists_every_key_but_types_in_the_save_prompt() {
    let mut v = PwmView::new("x", &pwm(), 2);
    let text = crate::tui::help::shown(&mut v, 80, 24).expect("? opens the keys");
    for key in ["←/→ h/l", "Home", "End", "0", "v b", "s", "q"] {
        assert!(text.contains(key), "{key}:\n{text}");
    }
    press(&mut v, KeyCode::Char('s'));
    assert!(crate::tui::help::shown(&mut v, 80, 24).is_none(), "? types");
}

#[test]
fn a_saved_logo_keeps_the_views_shape() {
    let mut v = PwmView::new("x", &pwm(), 2);
    let size = crate::figure::svg_size;
    let (cw, ch) = crate::figure::CELL;
    // The logo's own box on screen: letters at most a few cells wide, so a
    // wide terminal leaves room beside it that the figure does not take.
    for (w, h) in [(120, 30), (120, 40), (40, 30)] {
        screen(&mut v, w, h);
        assert!(v.drawn.width < w, "{w}x{h}: {:?}", v.drawn);
        let (fw, fh) = size(&v.figure(crate::figure::style::Sheet::OWN));
        let on_screen = v.drawn.width as f64 * cw / (v.drawn.height as f64 * ch);
        assert!((fw / fh - on_screen).abs() < 1e-6, "{w}x{h}: {fw}x{fh}");
    }
}

#[test]
fn a_long_motif_keeps_readable_letters() {
    // 41 positions saved from a narrow terminal: at least 8 points a letter.
    let wide: Vec<DnaBaseCount> = (0..41).map(|_| count(5, 1, 1, 1)).collect();
    let mut v = PwmView::new("x", &wide, 20);
    screen(&mut v, 60, 20);
    let (w, _) = crate::figure::svg_size(&v.figure(crate::figure::style::Sheet::OWN));
    assert!(w >= 76.0 + 8.0 * 41.0, "{w}");
}

/// The logo at its own size, no font asked for, is drawn as it always was.
#[test]
fn the_default_logo_is_unchanged() {
    let svg = super::figure("sites", &pwm(), 2, Sheet::OWN, &AxisTitles::default());
    assert_eq!(crate::figure::fingerprint(&svg), 4045831948546200367);
}

/// A panel `w` by `h` inches with every label at `font` points, alone,
/// in `palette`.
fn panel(w: f64, h: f64, font: f64, palette: Palette) -> Sheet {
    Sheet {
        page: (Some(w * 72.0), Some(h * 72.0)),
        look: StyleArgs {
            font_size: Some(font),
            panel_only: true,
            palette,
        },
    }
}

/// `--xlab` and `--ylab` title the logo's axes in place of its own.
#[test]
fn the_axis_titles_asked_replace_the_logos_own() {
    let sheet = panel(2.4, 1.8, 7.0, Palette::Colour);
    let labels = |titles: &AxisTitles| -> Vec<String> {
        let svg = super::figure("sites", &pwm(), 2, sheet, titles);
        (crate::figure::labels(&svg).into_iter())
            .map(|l| l.1)
            .collect()
    };
    let own = labels(&AxisTitles::default());
    for t in ["position relative to the site", "bits"] {
        assert!(own.iter().any(|l| l == t), "{t}: {own:?}");
    }
    let titles = AxisTitles {
        xlab: Some("position relative to the site (nt)".into()),
        ylab: Some("information (bits)".into()),
    };
    let asked = labels(&titles);
    for t in ["position relative to the site (nt)", "information (bits)"] {
        assert!(asked.iter().any(|l| l == t), "{t}: {asked:?}");
    }
    for t in ["position relative to the site", "bits"] {
        assert!(!asked.iter().any(|l| l == t), "{t}: {asked:?}");
    }
}

/// Each letter's fill, by letter, in the logo `svg`.
fn letter_fills(svg: &str) -> Vec<(char, String)> {
    (svg.lines()
        .filter(|l| l.contains(" scale(") && l.starts_with("<text")))
    .map(|l| {
        let fill = l.split(" fill=\"").nth(1).unwrap().split('"').next();
        let ch = l[l.find('>').unwrap() + 1..].chars().next().unwrap();
        (ch, fill.unwrap().to_string())
    })
    .collect()
}

/// `--palette gray` draws the letters in grays, each base its own shade,
/// and the site's mark in ink; colour, the default, in the logo colours.
#[test]
fn a_gray_logo_tells_its_letters_apart_by_shade() {
    let at = |palette| {
        let sheet = panel(2.4, 1.8, 7.0, palette);
        super::figure("sites", &pwm(), 2, sheet, &AxisTitles::default())
    };
    let colour = at(Palette::Colour);
    for (ch, fill) in letter_fills(&colour) {
        let b = BASES.iter().position(|&c| c == ch).unwrap();
        assert_eq!(fill, base_fill(b, Palette::Colour), "{ch}");
    }
    assert!(colour.contains(ACCENT));
    let gray = at(Palette::Gray);
    let fills = letter_fills(&gray);
    assert!(fills.len() >= 4, "{fills:?}");
    let mut shades: Vec<(char, String)> = Vec::new();
    for (ch, fill) in fills {
        let (r, g, b) = crate::figure::rgb(&fill);
        assert!(r == g && g == b, "{ch} in {fill}: not a gray");
        match shades.iter().find(|s| s.0 == ch) {
            Some(s) => assert_eq!(s.1, fill, "{ch} in one shade"),
            None => shades.push((ch, fill)),
        }
    }
    assert_eq!(shades.len(), 4, "{shades:?}");
    for (i, a) in shades.iter().enumerate() {
        for b in &shades[i + 1..] {
            let (x, y) = (crate::figure::rgb(&a.1).0, crate::figure::rgb(&b.1).0);
            assert!((x - y).abs() >= 40.0, "{a:?} and {b:?} too close");
        }
    }
    assert!(!gray.contains(ACCENT), "no accent in gray");
}

/// The position labels under the logo `svg`, as the boxes they take: the
/// numbers on the lowest line of numbers (the y axis' are above it).
fn position_labels(svg: &str) -> Vec<((f64, f64, f64, f64), String)> {
    let numbers: Vec<_> = (crate::figure::text_boxes(svg).into_iter())
        .filter(|(_, s)| s.parse::<i64>().is_ok())
        .collect();
    let low = numbers.iter().map(|(b, _)| b.3).fold(f64::MIN, f64::max);
    (numbers.into_iter())
        .filter(|(b, _)| (b.3 - low).abs() < 0.01)
        .collect()
}

/// At 2 inches across and 7 points, the window's position labels thin
/// (every 5, then every 10) so none runs into another; 0 always shows,
/// and with room every one does.
#[test]
fn position_labels_thin_to_keep_apart() {
    let window: Vec<DnaBaseCount> = (0..21).map(|i| count(10 + i, 5, 5, 5)).collect();
    let at = |w: f64, font: f64| {
        let sheet = panel(w, 1.8, font, Palette::Colour);
        position_labels(&super::figure(
            "sites",
            &window,
            10,
            sheet,
            &AxisTitles::default(),
        ))
    };
    for (w, font) in [(2.0, 7.0), (2.4, 7.0), (1.5, 7.0), (2.0, 8.0)] {
        let ticks = at(w, font);
        let shown: Vec<&str> = ticks.iter().map(|t| t.1.as_str()).collect();
        assert!(
            shown.contains(&"0") && shown.len() >= 3,
            "{w} in: {shown:?}"
        );
        for (i, (a, s)) in ticks.iter().enumerate() {
            for (b, z) in &ticks[i + 1..] {
                let apart = a.2 <= b.0 || b.2 <= a.0;
                assert!(apart, "{w} in: {s} at {a:?} runs into {z} at {b:?}");
            }
        }
    }
    assert_eq!(at(2.0, 7.0).len(), 5, "every 5 at 2 inches");
    assert_eq!(at(7.0, 7.0).len(), 21, "every one with room");
}

#[test]
fn the_keys_listed_are_the_keys_bound() {
    // Off the site, so `0`, Home and End each move.
    crate::tui::help::assert_lists_its_keys("a position off the site", || {
        let mut v = PwmView::new("x", &pwm(), 2);
        press(&mut v, KeyCode::Left);
        v
    });
}

/// The README's picture of the logo view: made-up counts around 2,000
/// sites, a motif at the site in a window of 10 either side, a leaning
/// base elsewhere. Written where `FABA_SCREENSHOTS_DIR` names.
#[test]
#[ignore = "writes the README's screenshot; run by hand"]
fn readme_screenshots_pwm() {
    let motif: [[usize; 4]; 5] = [
        [820, 120, 700, 360],
        [900, 40, 1060, 0],
        [2000, 0, 0, 0],
        [0, 2000, 0, 0],
        [880, 620, 60, 440],
    ];
    let pwm: Vec<DnaBaseCount> = (-10i64..=10)
        .map(
            |rel| match usize::try_from(rel + 2).ok().and_then(|i| motif.get(i)) {
                Some(&[a, c, g, t]) => count(a, c, g, t),
                None => {
                    // A leaning base, stronger nearer the site.
                    let i = (rel + 10) as usize;
                    let lean = 1100 - 80 * rel.unsigned_abs() as usize;
                    let mut c = [
                        300 + i * 37 % 120,
                        260 + i * 53 % 90,
                        280,
                        320 - i * 29 % 70,
                    ];
                    c[[0, 3, 2, 0, 3, 1][i % 6]] += lean;
                    count(c[0], c[1], c[2], c[3])
                }
            },
        )
        .collect();
    let mut v = PwmView::new("demo_m6a_sites", &pwm, 10);
    let mut term = Terminal::new(TestBackend::new(160, 40)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    crate::tui::screenshot::save("pwm", term.backend().buffer());
}
