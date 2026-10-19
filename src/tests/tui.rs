use super::*;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

#[test]
fn the_apply_key_is_ctrl_r_and_nothing_else() {
    let k = |c, m| KeyEvent::new(c, m);
    let ctrl_shift = KeyModifiers::CONTROL | KeyModifiers::SHIFT;
    assert!(is_apply(&k(KeyCode::Char('r'), KeyModifiers::CONTROL)));
    assert!(is_apply(&k(KeyCode::Char('R'), ctrl_shift)));
    assert!(is_apply(&apply_key()));
    for m in [
        KeyModifiers::NONE,
        KeyModifiers::SHIFT,
        KeyModifiers::CONTROL,
    ] {
        assert!(!is_apply(&k(KeyCode::Enter, m)), "{m:?}+Enter");
    }
    assert!(!is_apply(&k(KeyCode::Char('j'), KeyModifiers::CONTROL)));
    for c in ['G', 'g', 'y', 'r', 'R'] {
        assert!(!is_apply(&k(KeyCode::Char(c), KeyModifiers::NONE)), "{c}");
    }
    // Typed into a find, it types nothing.
    let mut f = Find::default();
    assert!(!f.key(&apply_key()) && f.text.is_empty());
}

#[test]
fn modified_keys_are_stray_but_the_apply_key_and_capitals_are_not() {
    let k = |c, m| KeyEvent::new(c, m);
    for m in [
        KeyModifiers::SHIFT,
        KeyModifiers::ALT,
        KeyModifiers::SUPER,
        KeyModifiers::CONTROL,
    ] {
        assert!(is_stray(&k(KeyCode::Enter, m)), "{m:?}+Enter");
    }
    assert!(is_stray(&k(KeyCode::Char('s'), KeyModifiers::CONTROL)));
    assert!(is_stray(&k(KeyCode::Char('j'), KeyModifiers::CONTROL)));
    assert!(is_stray(&k(KeyCode::Char('q'), KeyModifiers::ALT)));
    assert!(!is_stray(&k(KeyCode::Char('r'), KeyModifiers::CONTROL)));
    assert!(!is_stray(&k(KeyCode::Enter, KeyModifiers::NONE)));
    assert!(!is_stray(&k(KeyCode::Char('R'), KeyModifiers::SHIFT)));
    assert!(!is_stray(&k(KeyCode::BackTab, KeyModifiers::SHIFT)));
}

#[test]
fn centered_fits_inside_small_areas() {
    let r = centered(Rect::new(0, 0, 10, 4), 40, 12);
    assert_eq!((r.width, r.height), (10, 4));
    assert_eq!(first_visible(9, 10, 4), 6);
}

#[test]
fn output_folders_are_numbered_and_checked() {
    let tmp = tempfile::tempdir().unwrap();
    let base = tmp.path();
    assert_eq!(next_free(base, "out"), base.join("out"));
    std::fs::create_dir(base.join("out")).unwrap();
    assert_eq!(next_free(base, "out"), base.join("out2"));
    let out = base.join("out").to_string_lossy().into_owned();
    assert_eq!(output_problem(&out), None, "an empty folder is fine");
    std::fs::write(base.join("out/x"), b"").unwrap();
    assert!(output_problem(&out)
        .unwrap()
        .contains("already contains files"));
    let file = base.join("out/x").to_string_lossy().into_owned();
    assert!(output_problem(&file).unwrap().contains("is a file"));
    assert!(output_problem("").is_some());
}

#[test]
fn files_are_numbered_before_their_extension() {
    let tmp = tempfile::tempdir().unwrap();
    let base = tmp.path();
    assert_eq!(next_free_file(base, "r", ".html"), base.join("r.html"));
    std::fs::write(base.join("r.html"), b"").unwrap();
    assert_eq!(next_free_file(base, "r", ".html"), base.join("r2.html"));
    std::fs::write(base.join("r2.html"), b"").unwrap();
    assert_eq!(next_free_file(base, "r", ".html"), base.join("r3.html"));
}

/// A view that keeps what it is handed.
#[derive(Default)]
struct Keeps {
    apply: bool,
    /// The chords it takes.
    chords: Vec<KeyEvent>,
    keys: Vec<KeyEvent>,
    stray: usize,
}

impl Screen for Keeps {
    fn render(&mut self, _: &mut Frame) {}
    fn handle_key(&mut self, key: KeyEvent) {
        self.keys.push(key);
    }
    fn interrupt(&mut self) {}
    fn done(&self) -> bool {
        false
    }
}

impl View for Keeps {
    fn stray_enter(&mut self) {
        self.stray += 1;
    }
    fn takes_apply(&self) -> bool {
        self.apply
    }
    fn takes_chord(&self, key: &KeyEvent) -> bool {
        self.chords.contains(key)
    }
}

#[test]
fn feed_hands_on_only_what_a_view_takes() {
    let k = |c, m| KeyEvent::new(c, m);
    for apply in [false, true] {
        let mut v = Keeps {
            apply,
            ..Keeps::default()
        };
        feed(&mut v, apply_key());
        feed(&mut v, k(KeyCode::Char('s'), KeyModifiers::CONTROL));
        feed(&mut v, k(KeyCode::Enter, KeyModifiers::SHIFT));
        feed(&mut v, k(KeyCode::Char('r'), KeyModifiers::NONE));
        let want = if apply {
            vec![apply_key(), k(KeyCode::Char('r'), KeyModifiers::NONE)]
        } else {
            vec![k(KeyCode::Char('r'), KeyModifiers::NONE)]
        };
        assert_eq!(v.keys, want, "apply taken: {apply}");
        assert_eq!(v.stray, 1);
    }
}

#[test]
fn a_view_gets_only_the_chords_it_names() {
    let k = |c, m| KeyEvent::new(c, m);
    let ctrl_e = k(KeyCode::Char('e'), KeyModifiers::CONTROL);
    let ctrl_s = k(KeyCode::Char('s'), KeyModifiers::CONTROL);
    let ctrl_enter = k(KeyCode::Enter, KeyModifiers::CONTROL);
    // Taking none, as every view does by default: all dropped.
    let mut none = Keeps::default();
    for key in [ctrl_e, ctrl_s, ctrl_enter] {
        feed(&mut none, key);
    }
    assert!(none.keys.is_empty());
    assert_eq!(none.stray, 1, "Ctrl+Enter still says what applies");
    // Naming Ctrl+E: that one comes, the others are dropped as before.
    let mut named = Keeps {
        chords: vec![ctrl_e],
        ..Keeps::default()
    };
    for key in [ctrl_e, ctrl_s, ctrl_enter] {
        feed(&mut named, key);
    }
    assert_eq!(named.keys, vec![ctrl_e]);
    assert_eq!(named.stray, 1);
}

#[test]
fn a_plain_axis_keeps_its_line_and_labels_but_not_its_ticks() {
    let area = Rect::new(0, 0, 8, 2);
    let mut buf = ratatui::buffer::Buffer::empty(area);
    buf.set_string(0, 0, "──┴──┴──", ratatui::style::Style::new());
    buf.set_string(0, 1, "  CDS   ", ratatui::style::Style::new());
    plain_axis(&mut buf, area);
    let row = |y| (0..8).map(|x| buf[(x, y)].symbol()).collect::<String>();
    assert_eq!(row(0), "────────");
    assert_eq!(row(1), "  CDS   ");
}

#[test]
fn help_wraps_whole_keys_within_the_width() {
    use data_beans::interactive::ui::help_line;
    let keys = [
        ("type", "find"),
        ("↑/↓", "move"),
        ("Space", "mark file / choose folder"),
    ];
    let line = help_line(&keys);
    assert_eq!(
        wrap_help(line.clone(), 200).len(),
        1,
        "all on one line when it fits"
    );
    let lines = wrap_help(line, 40);
    assert_eq!(lines.len(), 2);
    assert!(
        lines.iter().all(|l| l.width() <= 40),
        "each within the width"
    );
    assert_eq!(
        lines[1].spans[1].content, "Space",
        "a key moves with what it does"
    );
}

#[test]
fn digits_count_a_number_as_it_is_written() {
    for n in [
        0,
        9,
        10,
        99,
        100,
        99_999,
        100_000,
        -1,
        -10,
        i64::MAX,
        i64::MIN,
    ] {
        assert_eq!(digits(n), n.to_string().len(), "{n}");
    }
}

/// Work on a thread comes back with its result, whether the wait took it
/// or left it; work that panics is an error saying with what.
#[test]
fn work_on_a_thread_comes_back_or_says_why_it_panicked() {
    let (waited, out) = on_thread("w", || Ok(3), |rx| (1, rx.recv().ok()));
    assert_eq!((waited, out.unwrap()), (1, 3));
    let (_, out) = on_thread("w", || Ok(4), |_| ((), None));
    assert_eq!(out.unwrap(), 4);
    let owned = String::from("an owned message");
    let (_, out) = on_thread(
        "the work",
        || -> anyhow::Result<u8> { panic!("{owned}") },
        |_| ((), None),
    );
    let err = out.unwrap_err().to_string();
    assert_eq!(err, "the work panicked: an owned message");
    let (_, out) = on_thread(
        "the work",
        || -> anyhow::Result<u8> { panic!("a fixed message") },
        |_| ((), None),
    );
    assert_eq!(
        out.unwrap_err().to_string(),
        "the work panicked: a fixed message"
    );
}

/// The spinner turns a frame every 200 ms, round and round.
#[test]
fn the_spinner_turns_with_the_time() {
    let frames: Vec<char> = SPINNER.chars().collect();
    let at = |ms| spinner_frame(std::time::Duration::from_millis(ms));
    assert_eq!(at(0), frames[0]);
    assert_eq!(at(199), frames[0]);
    assert_eq!(at(200), frames[1]);
    assert_eq!(at(200 * frames.len() as u64), frames[0]);
}

#[test]
fn a_stop_is_known_through_the_context_it_gathers() {
    let e = anyhow::Error::from(Stopped).context("picking the inputs");
    assert!(is_stopped(&e));
    assert!(!is_stopped(&anyhow::anyhow!("no such file")));
}

/// Row `y` of what `draw` drew `w` wide.
fn row(w: u16, y: u16, draw: impl FnOnce(&mut Frame)) -> String {
    let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, 3)).unwrap();
    t.draw(draw).unwrap();
    let buf = t.backend().buffer();
    (0..w).map(|x| buf[(x, y)].symbol()).collect()
}

#[test]
fn the_session_keys_end_the_title_row_whole() {
    let keys = vec![
        ("ctrl+alt+d".to_string(), "detach".to_string()),
        ("ctrl+alt+s".to_string(), "switch".to_string()),
    ];
    let title = "T".repeat(200);
    for w in [60, 120] {
        let top = row(w, 0, |f| {
            f.render_widget(Paragraph::new(title.as_str()), f.area());
            draw_session_keys(f, f.area(), &keys);
        });
        let shown = top.trim_end();
        assert!(shown.ends_with("ctrl+alt+s switch"), "{top:?}");
        assert!(shown.contains("ctrl+alt+d detach"), "{top:?}");
        assert!(top.starts_with('T'), "the title gives way: {top:?}");
    }
    // Plain, no keys: the title row as the view drew it.
    let top = row(40, 0, |f| {
        f.render_widget(Paragraph::new(title.as_str()), f.area());
        draw_session_keys(f, f.area(), &[]);
    });
    assert_eq!(top, "T".repeat(40));
}

/// A view that notes how it was asked to stop.
#[derive(Default)]
struct Stoppable {
    interrupted: bool,
    quit: bool,
}

impl Screen for Stoppable {
    fn render(&mut self, _: &mut Frame) {}
    fn handle_key(&mut self, _: KeyEvent) {}
    fn interrupt(&mut self) {
        self.interrupted = true;
    }
    fn done(&self) -> bool {
        self.interrupted || self.quit
    }
}

impl View for Stoppable {}

/// [`Stoppable`] with a quit path of its own.
#[derive(Default)]
struct Quits(Stoppable);

impl Screen for Quits {
    fn render(&mut self, _: &mut Frame) {}
    fn handle_key(&mut self, _: KeyEvent) {}
    fn interrupt(&mut self) {
        self.0.interrupt();
    }
    fn done(&self) -> bool {
        self.0.done()
    }
}

impl View for Quits {
    fn quit(&mut self) {
        self.0.quit = true;
    }
}

#[test]
fn a_kill_takes_ctrl_c_s_path_unless_the_view_has_its_own() {
    let mut v = Stoppable::default();
    View::quit(&mut v);
    assert!(v.interrupted, "by default as Ctrl-C");
    let mut q = Quits::default();
    View::quit(&mut q);
    assert!(q.0.quit && !q.0.interrupted);
}

#[test]
fn outside_a_kill_a_framed_view_takes_ctrl_c_as_ever() {
    let session = legume_tui::Session::plain(&TOOL);
    let mut q = Quits::default();
    let mut framed = Framed {
        view: &mut q,
        session: &session,
        keys: data_beans::interactive::session_keys(),
    };
    framed.interrupt();
    assert!(framed.keys.is_empty(), "plain: no session keys");
    assert!(q.0.interrupted && !q.0.quit);
}
