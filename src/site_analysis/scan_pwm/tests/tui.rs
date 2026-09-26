use super::*;
use crate::data::dna::Dna;
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
    press(&mut v, KeyCode::Char('m'));
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
    press(&mut v, KeyCode::Char('m'));
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
