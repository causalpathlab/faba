use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn key(code: KeyCode) -> KeyEvent {
    KeyEvent::new(code, KeyModifiers::NONE)
}

fn panel() -> String {
    let values = [0.0, 3.0, 10.0, 4.0, 1.0];
    let front = [0.0, 1.0, 8.0, 4.0, 0.0];
    let mut c = Canvas::new(300.0, 160.0);
    Bars {
        values: &values,
        front: Some(&front),
        accent: &|i| i < 2,
        colour: &|i| if i == 3 { "#0072b2" } else { BAR },
        y_scale: Scale::Log,
        y_max: None,
        ticks: vec![(0, "0".into()), (4, "L4 & <x>".into())],
        pointer: Some(2),
        marks: vec![1, 3],
        dividers: vec![2],
        title: "PANEL1".into(),
        x_title: "value".into(),
        y_title: "n".into(),
    }
    .draw(&mut c, 0.0, 0.0, 300.0, 160.0);
    c.glyph(10.0, 10.0, 20.0, 40.0, 'A', ACCENT);
    c.finish()
}

#[test]
fn a_divider_sits_on_the_bar_edge_and_drops_tick_stubs() {
    // Five bars over a 236-wide plot from x = 52: bar 2's left edge is 146.4.
    let svg = panel();
    assert!(
        svg.contains(r#"x1="146.40" y1="20.00" x2="146.40""#),
        "{svg}"
    );
    // With dividers, tick labels name spans: no solid stub under them. The
    // stub of the tick at bar 4 would sit at its centre, x = 264.4.
    assert!(!svg.contains(r#"x1="264.40""#));
}

#[test]
fn a_bar_takes_its_own_colour_unless_accented() {
    let svg = panel();
    // Bar 3 is not accented, so it takes its colour; bars 0 and 1 stay accent.
    assert!(svg.contains(r##"fill="#0072b2""##), "{svg}");
    assert!(svg.contains(&format!(r#"fill="{ACCENT}""#)));
}

#[test]
fn svg_is_well_formed_and_escaped() {
    let svg = panel();
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("L4 &amp; &lt;x&gt;"));
    assert!(svg.contains(ACCENT));
    let tree = usvg::Tree::from_str(&svg, &usvg::Options::default()).unwrap();
    assert_eq!(tree.size().width(), 300.0);
}

#[test]
fn save_writes_pdf_and_png() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("fig1.pdf");
    let paths = save(&panel(), prefix.to_str().unwrap()).unwrap();
    assert_eq!(paths.len(), 2);
    let pdf = std::fs::read(&paths[0]).unwrap();
    assert!(pdf.starts_with(b"%PDF"), "{}", paths[0]);
    let png = std::fs::read(&paths[1]).unwrap();
    assert!(png.starts_with(b"\x89PNG"), "{}", paths[1]);
    assert!(paths[0].ends_with("fig1.pdf") && paths[1].ends_with("fig1.png"));
}

#[test]
fn prompt_confirms_edits_and_cancels() {
    let mut p = SavePrompt::new("view1");
    assert!(!p.active() && p.footer().is_none());
    p.open();
    assert!(p.active());
    for _ in 0..5 {
        p.handle(key(KeyCode::Backspace));
    }
    for ch in "out/a".chars() {
        assert_eq!(p.handle(key(KeyCode::Char(ch))), None);
    }
    assert_eq!(p.handle(key(KeyCode::Enter)).as_deref(), Some("out/a"));
    assert!(!p.active());
    p.report(Ok(vec!["out/a.pdf".into()]));
    assert!(p.footer().is_some());
    p.dismiss();
    assert!(p.footer().is_none());

    p.open();
    assert!(p.footer().is_some(), "the last name is offered again");
    p.handle(key(KeyCode::Esc));
    assert!(!p.active());
}

#[test]
fn empty_shapes_are_not_written() {
    let mut c = Canvas::new(100.0, 200.0);
    c.rect(1.0, 1.0, 10.0, 0.0, INK);
    c.rect(1.0, 1.0, 0.0, 10.0, INK);
    c.rect(1.0, 1.0, 10.0, f64::NAN, INK);
    let values = [Some(0.0), Some(1.5), None, Some(-0.5), Some(0.0)];
    Diverging {
        values: &values,
        ticks: Vec::new(),
        pointer: None,
        title: String::new(),
        x_title: String::new(),
        y_title: String::new(),
        label: &|v| format!("{v:+.1}"),
    }
    .draw(&mut c, 0.0, 0.0, 100.0, 200.0);
    let svg = c.finish();
    assert!(!svg.contains(r#"height="0.00""#), "{svg}");
    assert!(!svg.contains("NaN"), "{svg}");
    // Two non-zero bars plus the background.
    assert_eq!(svg.matches("<rect").count(), 3, "{svg}");
}
