use super::*;
use data_beans::aux::feature_rows::M6A;
use ratatui::crossterm::event::KeyModifiers;

fn press(p: &mut PathPicker, code: KeyCode) {
    crate::tui::feed(p, KeyEvent::new(code, KeyModifiers::NONE));
}

fn type_text(p: &mut PathPicker, text: &str) {
    for ch in text.chars() {
        press(p, KeyCode::Char(ch));
    }
}

/// `root/{plain, prof/{inner}}`, with `prof` a faba output directory.
fn tree() -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    std::fs::create_dir_all(root.join("plain")).unwrap();
    std::fs::create_dir_all(root.join("prof/inner")).unwrap();
    std::fs::write(root.join(format!("prof/{M6A}_sites.parquet")), b"").unwrap();
    tmp
}

fn select(p: &mut PathPicker, name: &str) {
    assert!(p.browser.list.select(name));
}

#[test]
fn browse_choose_an_input_and_name_the_output() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut p = PathPicker::new(root.clone(), None, None);
    let names: Vec<&str> = p.browser.list.shown().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "plain", "prof"]);
    assert!(
        p.browser.list.shown().nth(2).unwrap().tagged
            && !p.browser.list.shown().nth(1).unwrap().tagged
    );

    // Not a faba directory: refused, still browsing.
    select(&mut p, "plain");
    press(&mut p, KeyCode::Char(' '));
    assert!(p.error.is_some() && matches!(p.step, Step::Input));

    // Into `prof` and back up lands on it again.
    select(&mut p, "prof");
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.browser.cwd, root.join("prof"));
    press(&mut p, KeyCode::Left);
    assert_eq!(p.browser.cwd, root);
    assert_eq!(p.browser.highlighted().unwrap().name, "prof");

    press(&mut p, KeyCode::Char(' '));
    let suggested = format!("{}_qc", root.join("prof").display());
    assert!(matches!(&p.step, Step::Output { line, .. } if line.text() == Some(&*suggested)));

    // A non-empty output is refused and stays open for editing.
    let Step::Output { line, .. } = &mut p.step else {
        unreachable!()
    };
    line.open("");
    type_text(&mut p, &root.join("prof").to_string_lossy());
    press(&mut p, KeyCode::Enter);
    assert!(p.error.is_some() && !p.done());
    let Step::Output { line, .. } = &mut p.step else {
        unreachable!()
    };
    line.open("");
    let out = root.join("new").to_string_lossy().into_owned();
    type_text(&mut p, &out);
    press(&mut p, KeyCode::Enter);
    assert_eq!(
        p.decision,
        Some(Some((
            root.join("prof").to_string_lossy().into_owned(),
            out
        )))
    );
}

#[test]
fn only_what_is_missing_is_asked() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let input = root.join("prof");
    let out = root.join("out").to_string_lossy().into_owned();

    // Both given: nothing to ask.
    let p = PathPicker::new(root.clone(), Some(input.clone()), Some(out.clone()));
    assert!(p.done());

    // The input given: straight to the output, and nothing listed.
    let p = PathPicker::new(root.clone(), Some(input.clone()), None);
    assert!(matches!(p.step, Step::Output { .. }) && p.browser.list.is_empty());

    // An output given but not empty: asked again, with why.
    let p = PathPicker::new(
        root.clone(),
        Some(input.clone()),
        Some(input.to_string_lossy().into()),
    );
    assert!(matches!(p.step, Step::Output { .. }) && p.error.is_some());

    // The output given: browse, and choosing finishes.
    let mut p = PathPicker::new(input.clone(), None, Some(out.clone()));
    press(&mut p, KeyCode::Char('.'));
    assert!(
        p.decision.is_none() && p.browser.list.find.text == ".",
        "letters find"
    );
    press(&mut p, KeyCode::Esc);
    assert!(p.decision.is_none(), "Esc clears the find first");
    crate::tui::feed(
        &mut p,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    assert_eq!(
        p.decision,
        Some(Some((input.to_string_lossy().into_owned(), out)))
    );

    // Esc cancels.
    let mut p = PathPicker::new(root, None, None);
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.decision, Some(None));
}

#[test]
fn draws_at_any_terminal_size() {
    use ratatui::backend::TestBackend;
    use ratatui::Terminal;
    let tmp = tree();
    let mut p = PathPicker::new(tmp.path().to_path_buf(), None, None);
    for (w, h) in [(120, 40), (30, 8), (10, 3)] {
        let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
        term.draw(|f| p.render(f)).unwrap();
    }
    let mut term = Terminal::new(TestBackend::new(120, 40)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(screen.contains("prof/") && screen.contains("faba output"));
}
