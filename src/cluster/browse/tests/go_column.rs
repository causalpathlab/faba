//! The gene modules' GO terms beside the grid: read on the first `a`, a
//! column of each module's top terms on its own row, `a` again to hide or
//! show it, `G` for the full pane, and no column where there is no room.

use super::super::go::Ask;
use super::*;

/// [`loaded`] with 40 genes `ENSG..._G{i}`: gene module 1 the first 20, 2
/// the next 4, 3 the rest; and a GMT whose `TERM1_SET_OF_GENES` to `TERM8_...` each hold 12 of
/// gene module 1's (each a different 12) and `OTHER` gene module 3's, so
/// gene module 1 has eight terms and 2 none (none of its genes annotated).
fn with_go(dir: &std::path::Path) -> (Modules, Ask) {
    let mut m = loaded();
    m.gene_main = (0..40)
        .map(|i| {
            let role = match i {
                0..20 => 0,
                20..24 => 1,
                _ => 2,
            };
            (format!("ENSG{i:011}_G{i}").into(), role, 0.9)
        })
        .collect();
    let set = |name: String, genes: std::ops::Range<usize>| {
        let g: Vec<String> = genes.map(|i| format!("G{i}")).collect();
        format!("{name}\t-\t{}", g.join("\t"))
    };
    let mut lines: Vec<String> = (0..8)
        .map(|t| set(format!("TERM{}_SET_OF_GENES", t + 1), t..t + 12))
        .collect();
    lines.push(set("OTHER".into(), 24..40));
    let gmt = dir.join("sets.gmt");
    std::fs::write(&gmt, lines.join("\n")).unwrap();
    let ask = Ask {
        gaf: None,
        gmt: Some(gmt.to_string_lossy().into()),
        obo: None,
        min_weight: 0.5,
    };
    (m, ask)
}

/// The screen once GO has been read, `w` by `h`.
fn read(v: &mut ModulesView, w: u16, h: u16) -> String {
    for _ in 0..500 {
        if Screen::tick(v) {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    screen(v, w, h)
}

/// The two lines of gene module `name`'s row in the grid.
fn row_of(s: &str, name: &str) -> String {
    let lines: Vec<&str> = s.lines().collect();
    let label = format!("{name:<4}");
    let y = (lines.iter())
        .position(|l| l.contains(&format!("│{label}")))
        .unwrap_or_else(|| panic!("no row {name}:\n{s}"));
    format!("{}\n{}", lines[y], lines[y + 1])
}

#[test]
fn a_reads_go_then_shows_each_modules_terms_on_its_row() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    let mut v = ModulesView::new("test", &m).with_go(ask);
    let s = screen(&mut v, 140, 30);
    assert!(
        !s.contains("GO (q") && !s.contains("TERM"),
        "no column yet:\n{s}"
    );
    assert!(s.contains("a GO column") && s.contains("G GO terms"), "{s}");
    press(&mut v, KeyCode::Char('a'));
    assert!(screen(&mut v, 140, 30).contains("reading GO..."));
    let s = read(&mut v, 140, 30);
    assert!(s.contains("GO (q <= 0.05)"), "{s}");
    // Gene module 1's terms, strongest first, as many as fit, then how
    // many more; gene module 2 none.
    let m1 = row_of(&s, "M1");
    assert!(m1.contains("TERM1_") && m1.contains("TERM2_"), "{m1}");
    let more = m1.split('+').nth(1).expect("+N more");
    let n: usize = more.split(' ').next().unwrap().parse().unwrap();
    let shown = (1..=8)
        .filter(|t| m1.contains(&format!("TERM{t}_")))
        .count();
    assert!(n > 0 && shown + n == 8, "{shown} shown, {n} more:\n{m1}");
    assert!(more.starts_with(&format!("{n} more")), "{m1}");
    let m2 = row_of(&s, "M2");
    assert!(
        m2.contains("none at q <= 0.05") && !m2.contains("TERM"),
        "{m2}"
    );
    // `a` again hides it, and again shows it, nothing read anew.
    press(&mut v, KeyCode::Char('a'));
    let s = screen(&mut v, 140, 30);
    assert!(!s.contains("GO (q") && !s.contains("TERM1"), "{s}");
    press(&mut v, KeyCode::Char('a'));
    let s = screen(&mut v, 140, 30);
    assert!(
        s.contains("GO (q <= 0.05)") && !s.contains("reading"),
        "{s}"
    );
}

#[test]
fn g_opens_the_full_pane_and_esc_or_g_closes_it() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    let mut v = ModulesView::new("test", &m).with_go(ask);
    v.module = 0;
    press(&mut v, KeyCode::Char('G'));
    let s = read(&mut v, 140, 30);
    assert!(
        s.contains("GO · gene module M1") && s.contains(" TERM1_SET_OF_GENES "),
        "{s}"
    );
    press(&mut v, KeyCode::Esc);
    let s = screen(&mut v, 140, 30);
    assert!(s.contains("gene modules ·") && !s.contains("GO (q"), "{s}");
    press(&mut v, KeyCode::Char('G'));
    assert!(screen(&mut v, 140, 30).contains("GO · gene module"));
    press(&mut v, KeyCode::Char('G'));
    assert!(!v.done() && screen(&mut v, 140, 30).contains("gene modules ·"));
    let keys = crate::tui::help::shown(&mut v, 140, 50).expect("? lists the keys");
    for key in ["a", "G", "v"] {
        assert!(
            keys.lines().any(|l| l.contains(&format!(" {key} "))),
            "{key}:\n{keys}"
        );
    }
}

#[test]
fn a_screen_too_narrow_hides_the_column_and_says_so() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    let mut v = ModulesView::new("test", &m).with_go(ask);
    press(&mut v, KeyCode::Char('a'));
    let s = read(&mut v, 96, 30);
    assert!(!s.contains("GO (q") && !s.contains("TERM1"), "{s}");
    assert!(s.contains("no room for the GO column"), "{s}");
    // Narrower than wide, the terms' names are shortened to fit.
    let s = screen(&mut v, 112, 30);
    assert!(s.contains("GO (q"), "{s}");
    assert!(row_of(&s, "M1").contains('…'), "{s}");
}

/// A read that failed is tried again on `a`, and the column stays on to
/// show what the retry reads.
#[test]
fn a_after_a_failed_read_retries_with_the_column_on() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, mut ask) = with_go(tmp.path());
    let later = tmp.path().join("later.gmt");
    let real = ask.gmt.replace(later.to_string_lossy().into()).unwrap();
    let mut v = ModulesView::new("test", &m).with_go(ask);
    press(&mut v, KeyCode::Char('a'));
    let s = read(&mut v, 140, 30);
    assert!(
        s.contains("GO (q") && !s.contains("TERM1"),
        "the error:\n{s}"
    );
    // The file is there now: `a` reads it, and its terms show.
    std::fs::copy(real, &later).unwrap();
    press(&mut v, KeyCode::Char('a'));
    let s = read(&mut v, 140, 30);
    assert!(row_of(&s, "M1").contains("TERM1_"), "{s}");
}

/// With no room for the column, every `a` says so; `G` still opens GO.
#[test]
fn each_a_without_room_says_so() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    let mut v = ModulesView::new("test", &m).with_go(ask);
    press(&mut v, KeyCode::Char('a'));
    let s = read(&mut v, 96, 30);
    assert!(s.contains("no room for the GO column"), "{s}");
    // Each a still turns the column off and on, saying which.
    press(&mut v, KeyCode::Char('a'));
    assert!(screen(&mut v, 96, 30).contains("GO column off"));
    press(&mut v, KeyCode::Char('a'));
    assert!(screen(&mut v, 96, 30).contains("no room for the GO column"));
    press(&mut v, KeyCode::Char('a'));
    // Off when the screen widens: the column the user hid stays hidden.
    assert!(!screen(&mut v, 260, 30).contains("GO (q <="));
    press(&mut v, KeyCode::Char('G'));
    assert!(screen(&mut v, 96, 30).contains("GO · gene module"));
}
