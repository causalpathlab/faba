use super::*;
use inputs::InputsFocus;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use ratatui::Terminal;

fn press(a: &mut App, code: KeyCode) {
    crate::tui::feed(a, KeyEvent::new(code, KeyModifiers::NONE));
}

fn app_with_bams() -> (tempfile::TempDir, App) {
    let tmp = tempfile::tempdir().unwrap();
    for f in ["sample_A.bam", "control_A.bam", "genes.gff", "genome.fa"] {
        std::fs::write(tmp.path().join(f), b"").unwrap();
    }
    let mut a = App::new(run_cmd(), tmp.path().to_path_buf());
    for n in ["control_A.bam", "sample_A.bam"] {
        assert!(a.inputs.bams.list.select(n));
        press(&mut a, KeyCode::Char(' '));
    }
    assert!(a.inputs.bams.list.select("control_A.bam"));
    press(&mut a, KeyCode::Char(' '));
    a.inputs.gff = Some(tmp.path().join("genes.gff"));
    a.inputs.genome = Some(tmp.path().join("genome.fa"));
    a.inputs.output = tmp.path().join("faba_out").to_string_lossy().into_owned();
    // Off the BAM list, where letters type into its find.
    a.inputs.focus = InputsFocus::Rows;
    a.refresh();
    (tmp, a)
}

#[test]
fn tab_walks_the_screens_but_run_needs_a_job() {
    let (_t, mut a) = app_with_bams();
    a.inputs.focus = InputsFocus::Bams;
    // The BAM list, its rows, then each screen; Shift+Tab walks back.
    press(&mut a, KeyCode::Tab);
    assert_eq!((a.page, a.inputs.focus), (Page::Inputs, InputsFocus::Rows));
    for want in [Page::Steps, Page::Flags, Page::Inputs] {
        press(&mut a, KeyCode::Tab);
        assert_eq!(a.page, want);
    }
    assert_eq!(a.inputs.focus, InputsFocus::Bams);
    press(&mut a, KeyCode::BackTab);
    assert_eq!(a.page, Page::Flags);
    press(&mut a, KeyCode::BackTab);
    press(&mut a, KeyCode::BackTab);
    assert_eq!((a.page, a.inputs.focus), (Page::Inputs, InputsFocus::Rows));
    press(&mut a, KeyCode::Char('4'));
    assert_eq!(a.page, Page::Inputs, "no run yet");
}

#[test]
fn the_preview_lists_the_command_and_its_problems() {
    let (_t, mut a) = app_with_bams();
    let argv = a.argv();
    assert_eq!(&argv[..2], ["run", "--batch-process"]);
    assert!(argv.iter().any(|w| w == "--control-bam"));
    assert!(argv.windows(2).any(|w| w == ["-o", "."]));
    assert!(a.problems().is_empty(), "{:?}", a.problems());
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    assert!(a.preview);
    press(&mut a, KeyCode::Esc);
    assert!(!a.preview);
    a.inputs.gff = None;
    assert!(a.problems().iter().any(|p| p.contains("GFF")));
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    assert!(a.job.is_none(), "problems block the start");
}

fn ctrl_r(a: &mut App) {
    crate::tui::feed(a, KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL));
}

#[test]
fn only_ctrl_r_previews_and_starts() {
    let (tmp, mut a) = app_with_bams();
    fake_program(&mut a, tmp.path(), "exit 0");
    assert!(a.problems().is_empty(), "{:?}", a.problems());
    // `G` and `y` do nothing on any screen; nor does plain Enter.
    for page in ['1', '2', '3'] {
        press(&mut a, KeyCode::Char(page));
        for code in [KeyCode::Char('G'), KeyCode::Char('y')] {
            press(&mut a, code);
            assert!(!a.preview && a.job.is_none(), "{page} {code:?}");
        }
    }
    press(&mut a, KeyCode::Char('2'));
    press(&mut a, KeyCode::Enter);
    assert!(!a.preview);
    let hint = a.note.clone().unwrap_or_default();
    assert!(
        hint.contains("c copies")
            && hint.contains("p prints it and leaves")
            && hint.contains("--batch-process"),
        "{hint}"
    );
    ctrl_r(&mut a);
    assert!(a.preview);
    for code in [KeyCode::Char('G'), KeyCode::Char('y'), KeyCode::Enter] {
        press(&mut a, code);
        assert!(a.preview && a.job.is_none(), "{code:?} does not start");
    }
    assert!(a.note.as_deref().is_some_and(|n| n.contains("c copies")));
    ctrl_r(&mut a);
    assert!(a.job.is_some() && !a.preview, "Ctrl+R starts");
    a.job
        .as_mut()
        .unwrap()
        .handle
        .take()
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn shift_and_alt_enter_are_not_enter() {
    let (_t, mut a) = app_with_bams();
    press(&mut a, KeyCode::Char('3'));
    let i = a.flag_row().unwrap();
    let before = a.form.fields[i].value.clone();
    for m in [KeyModifiers::SHIFT, KeyModifiers::ALT] {
        crate::tui::feed(&mut a, KeyEvent::new(KeyCode::Enter, m));
        assert!(
            a.editing.is_none() && !a.preview,
            "{m:?}+Enter opens nothing"
        );
        assert_eq!(
            a.form.fields[i].value, before,
            "{m:?}+Enter toggles nothing"
        );
        assert!(a.note.as_deref().is_some_and(|n| n.contains("Ctrl+R")));
    }
    // Ctrl+letters are not letters: Ctrl+Q does not quit.
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('q'), KeyModifiers::CONTROL),
    );
    assert!(!a.quit);
    // Only Ctrl+R goes: not Ctrl+J, nor Ctrl+Enter.
    for code in [KeyCode::Char('j'), KeyCode::Enter] {
        crate::tui::feed(&mut a, KeyEvent::new(code, KeyModifiers::CONTROL));
        assert!(!a.preview, "Ctrl+{code:?}");
    }
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    assert!(a.preview);
}

#[test]
fn q_asks_before_losing_a_setup() {
    let (_t, mut a) = app_with_bams();
    press(&mut a, KeyCode::Char('q'));
    assert!(!a.quit && a.note.as_deref().is_some_and(|n| n.contains("q again")));
    // Any other key forgets the question.
    press(&mut a, KeyCode::Tab);
    press(&mut a, KeyCode::Char('q'));
    assert!(!a.quit);
    press(&mut a, KeyCode::Char('q'));
    assert!(a.quit);
    // With nothing picked, q leaves at once.
    let tmp = tempfile::tempdir().unwrap();
    let mut b = App::new(run_cmd(), tmp.path().to_path_buf());
    b.inputs.focus = InputsFocus::Rows;
    press(&mut b, KeyCode::Char('q'));
    assert!(b.quit);
}

#[test]
fn the_output_folder_is_asked_for_never_assumed() {
    let (tmp, mut a) = app_with_bams();
    a.inputs.output.clear();
    assert!(a.problems().iter().any(|p| p.contains("no output folder")));
    let s = screen(&mut a, 200, 20);
    assert!(s.contains("Enter to name it"), "{s}");
    // The apply key opens the output line, holding a suggestion.
    press(&mut a, KeyCode::Char('3'));
    ctrl_r(&mut a);
    assert!(!a.preview && a.page == Page::Inputs);
    let suggested = tmp.path().join("faba_out").to_string_lossy().into_owned();
    match &a.editing {
        Some((Target::Output, line)) => assert_eq!(line.text(), Some(suggested.as_str())),
        _ => panic!("the output line is open"),
    }
    // Esc leaves it unnamed; nothing goes on.
    press(&mut a, KeyCode::Esc);
    assert!(a.inputs.output().is_empty() && !a.preview);
    // Naming it goes on to the preview the apply key asked for.
    ctrl_r(&mut a);
    press(&mut a, KeyCode::Enter);
    assert_eq!(a.inputs.output(), suggested);
    assert!(a.preview && a.problems().is_empty(), "{:?}", a.problems());
    // `p` asks too, then prints with the folder named.
    let mut b = app_with_bams().1;
    b.inputs.output.clear();
    press(&mut b, KeyCode::Char('p'));
    assert!(!b.quit && b.editing.is_some());
    for c in "/x".chars() {
        press(&mut b, KeyCode::Char(c));
    }
    press(&mut b, KeyCode::Enter);
    assert!(b.quit);
    assert!(b
        .printed
        .as_deref()
        .is_some_and(|l| l.contains("faba_out/x")));
}

#[test]
fn letters_in_the_bam_list_find_rather_than_command() {
    let (_t, mut a) = app_with_bams();
    a.inputs.focus = InputsFocus::Bams;
    for c in "qcp1b".chars() {
        press(&mut a, KeyCode::Char(c));
    }
    assert!(!a.quit && a.page == Page::Inputs && a.note.is_none());
    assert_eq!(a.inputs.bams.list.find.text, "qcp1b");
    assert!(a.inputs.bams.list.is_empty());
    let s = screen(&mut a, 200, 20);
    assert!(
        s.contains("find: qcp1b") && s.contains("nothing matches"),
        "{s}"
    );
    // Esc clears it; the commands work again off the list.
    press(&mut a, KeyCode::Esc);
    for c in "sample".chars() {
        press(&mut a, KeyCode::Char(c));
    }
    assert_eq!(a.inputs.bams.highlighted().unwrap().name, "sample_A.bam");
    // The file pop-up finds the same way.
    press(&mut a, KeyCode::Tab);
    press(&mut a, KeyCode::Enter);
    for c in "genes".chars() {
        press(&mut a, KeyCode::Char(c));
    }
    let (_, b) = a.picking.as_ref().unwrap();
    assert_eq!(b.list.len(), 1);
    press(&mut a, KeyCode::Enter);
    assert!(a.picking.is_none() && a.inputs.gff.is_some());
    press(&mut a, KeyCode::Char('q'));
    assert!(!a.quit, "q asks first with BAMs picked");
}

#[test]
fn a_folder_with_files_is_refused_as_the_output() {
    let (tmp, mut a) = app_with_bams();
    a.inputs.output.clear();
    ctrl_r(&mut a);
    // The BAMs' own folder has files in it: asked again, with why.
    for _ in 0..40 {
        press(&mut a, KeyCode::Backspace);
    }
    for c in tmp.path().to_string_lossy().chars() {
        press(&mut a, KeyCode::Char(c));
    }
    press(&mut a, KeyCode::Enter);
    assert!(a.inputs.output().is_empty() && !a.preview);
    assert!(matches!(a.editing, Some((Target::Output, _))));
    let s = screen(&mut a, 200, 20);
    assert!(s.contains("already contains files"), "{s}");
}

#[test]
fn a_download_shows_its_progress_in_a_pop_up_esc_hides_d_shows() {
    let (tmp, mut a) = app_with_bams();
    let f = fetch::Fetch::idle("a reference", tmp.path().join("ref"));
    if let Ok(mut s) = f.status.lock() {
        s.step = "genome".into();
        s.bytes = Some((412 << 20, Some(846 << 20)));
    }
    a.fetch = Some(f);
    a.fetch_shown = true;
    let s = screen(&mut a, 120, 30);
    assert!(
        s.contains("downloading a reference") && s.contains("412 / 846 MB"),
        "{s}"
    );
    press(&mut a, KeyCode::Esc);
    assert!(!a.fetch_shown && a.fetch.is_some());
    press(&mut a, KeyCode::Char('d'));
    assert!(a.fetch_shown);
}

fn mouse(a: &mut App, kind: ratatui::crossterm::event::MouseEventKind, hit: Hit) {
    let _ = screen(a, 160, 30);
    let (column, row) = a
        .hits
        .spot(&hit)
        .unwrap_or_else(|| panic!("{hit:?} not drawn"));
    let event = ratatui::crossterm::event::MouseEvent {
        kind,
        column,
        row,
        modifiers: KeyModifiers::NONE,
    };
    a.mouse(event);
}

#[test]
fn clicks_pick_screens_panes_and_rows_and_the_wheel_moves() {
    use ratatui::crossterm::event::{MouseButton, MouseEventKind as M};
    let click = M::Down(MouseButton::Left);
    let (_t, mut a) = app_with_bams();
    mouse(&mut a, click, Hit::Bams(Some(1)));
    assert_eq!(a.inputs.focus, InputsFocus::Bams);
    assert_eq!(a.inputs.bams.list.at, 1);
    mouse(&mut a, click, Hit::Rows(Some(3)));
    assert_eq!((a.inputs.focus, a.inputs.row), (InputsFocus::Rows, 3));
    // The wheel over a pane focuses it and moves its cursor.
    mouse(&mut a, M::ScrollDown, Hit::Bams(None));
    assert_eq!(
        (a.inputs.focus, a.inputs.bams.list.at),
        (InputsFocus::Bams, 2)
    );
    mouse(&mut a, click, Hit::Tab(Page::Steps));
    assert_eq!(a.page, Page::Steps);
    mouse(&mut a, click, Hit::Step(2));
    assert_eq!(a.steps.at, 2);
    mouse(&mut a, click, Hit::Tab(Page::Flags));
    mouse(&mut a, click, Hit::Flag(3));
    assert_eq!(a.flags_at, 3);
    // No Run screen without a run.
    mouse(&mut a, click, Hit::Tab(Page::Run));
    assert_eq!(a.page, Page::Flags);
    // A pop-up takes no clicks.
    a.preview = true;
    mouse(&mut a, click, Hit::Tab(Page::Inputs));
    assert_eq!(a.page, Page::Flags);
}

#[test]
fn ctrl_r_and_the_buttons_preview_and_cancel_with_the_folder_as_named() {
    use ratatui::crossterm::event::{MouseButton, MouseEventKind as M};
    let click = M::Down(MouseButton::Left);
    let (tmp, mut a) = app_with_bams();
    crate::tui::feed(
        &mut a,
        KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL),
    );
    assert!(a.preview, "Ctrl+R previews");
    // The preview shows the folder as named, not the script's `.`.
    let s = screen(&mut a, 200, 40);
    let out = tmp.path().join("faba_out");
    assert!(s.contains(&format!("-o {}", out.display())), "{s}");
    assert!(s.contains("[ ▶ start ]") && s.contains("[ cancel ]"), "{s}");
    mouse(&mut a, click, Hit::Cancel);
    assert!(!a.preview);
    mouse(&mut a, click, Hit::Preview);
    assert!(a.preview);
    // The copied command names the folder too.
    assert!(a.argv_out().contains(&out.to_string_lossy().into_owned()));
}

#[test]
fn a_run_that_ends_says_so_in_a_pop_up_once() {
    for (body, title, why) in [
        ("exit 0", "faba run finished", None),
        (
            "echo 'Error: no reads' >&2; exit 1",
            "faba run failed",
            Some("no reads"),
        ),
    ] {
        let (tmp, mut a) = app_with_bams();
        fake_program(&mut a, tmp.path(), body);
        a.start().unwrap();
        let job = a.job.as_mut().unwrap();
        job.handle.take().unwrap().join().unwrap();
        assert!(a.tick() && a.ended_shown, "{body}");
        let notice = a.take_notice().unwrap_or_default();
        assert_eq!(
            format!(" {notice} "),
            format!(" {title} "),
            "told to whoever is away"
        );
        assert!(a.take_notice().is_none(), "once");
        let s = screen(&mut a, 160, 30);
        assert!(s.contains(title) && s.contains("[ ok ]"), "{s}");
        assert!(why.is_none_or(|w| s.contains(w)), "{s}");
        press(&mut a, KeyCode::Enter);
        assert!(!a.ended_shown && a.page == Page::Run);
        a.tick();
        assert!(!a.ended_shown, "told once");
    }
}

#[test]
fn d_offers_references_only_on_request() {
    let (_t, mut a) = app_with_bams();
    let s = screen(&mut a, 200, 20);
    assert!(!s.contains("download a reference"));
    a.inputs.gff = None;
    let s = screen(&mut a, 200, 20);
    assert!(s.contains("d downloads"), "{s}");
    press(&mut a, KeyCode::Char('d'));
    let s = screen(&mut a, 200, 30);
    assert!(
        s.contains("download a reference") && s.contains("GENCODE"),
        "{s}"
    );
    press(&mut a, KeyCode::Esc);
    assert!(a.catalogue.is_none());
}

#[test]
fn c_copies_the_command_from_every_screen() {
    let (_t, mut a) = app_with_bams();
    for page in ['1', '2', '3'] {
        press(&mut a, KeyCode::Char(page));
        press(&mut a, KeyCode::Char('c'));
        assert!(!a.preview);
        assert!(a.note.as_deref().is_some_and(|n| n.contains("copied")));
    }
    let s = screen(&mut a, 200, 20);
    assert!(s.contains("copied"), "{s}");
    press(&mut a, KeyCode::BackTab);
    let s = screen(&mut a, 200, 20);
    assert!(s.contains("c copy the command") && !s.contains("/G"), "{s}");
}

#[test]
fn a_bad_flag_value_is_a_problem() {
    let (_t, mut a) = app_with_bams();
    let i = a
        .form
        .fields
        .iter()
        .position(|f| f.long == "max-threads")
        .unwrap();
    a.form.fields[i].set_text("many");
    assert!(a.problems().iter().any(|p| p.contains("max-threads")));
}

#[test]
fn starting_writes_the_script_and_runs_the_child() {
    let (tmp, mut a) = app_with_bams();
    // A fake faba: prints a line and a bar frame, then exits.
    let fake = tmp.path().join("fake-faba");
    std::fs::write(
        &fake,
        "#!/usr/bin/env bash\necho 'Step 1/5: gene counting' >&2\nprintf '##-- 1/2 (1s) genes\\r' >&2\n",
    )
    .unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    a.program = fake.clone();
    a.start().unwrap();
    let out = PathBuf::from(a.inputs.output());
    assert!(out.join(script::SCRIPT).exists());
    assert_eq!(a.page, Page::Run);
    let job = a.job.as_mut().unwrap();
    job.handle.take().unwrap().join().unwrap();
    let log = job.log.lock().unwrap();
    assert!(log.lines.iter().any(|l| l.contains("gene counting")));
    assert_eq!(log.ended, Some(Ok(())));
}

#[test]
fn every_page_draws_at_any_size() {
    let (_t, mut a) = app_with_bams();
    for page in [Page::Inputs, Page::Steps, Page::Flags] {
        a.page = page;
        for (w, h) in [(140, 44), (60, 16), (20, 6)] {
            let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
            term.draw(|f| a.render(f)).unwrap();
            a.preview = true;
            term.draw(|f| a.render(f)).unwrap();
            a.preview = false;
        }
    }
    let mut term = Terminal::new(TestBackend::new(140, 44)).unwrap();
    a.page = Page::Inputs;
    term.draw(|f| a.render(f)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(screen.contains("[fg]") && screen.contains("[bg]"));
}

#[test]
fn the_run_page_and_pop_ups_draw_at_any_size() {
    let (tmp, mut a) = app_with_bams();
    let fake = tmp.path().join("fake-faba");
    std::fs::write(
        &fake,
        "#!/usr/bin/env bash\nfor i in 1 2 3; do echo \"Step $i/5: part $i\" >&2; done\nexit 3\n",
    )
    .unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    a.program = fake;
    a.start().unwrap();
    a.job
        .as_mut()
        .unwrap()
        .handle
        .take()
        .unwrap()
        .join()
        .unwrap();
    assert!(matches!(
        a.job.as_ref().unwrap().log.lock().unwrap().ended,
        Some(Err(child::Failed::Exit(_)))
    ));
    press(&mut a, KeyCode::Up);
    for (w, h) in [(140, 44), (60, 16), (20, 6)] {
        let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
        term.draw(|f| a.render(f)).unwrap();
    }
    // The file pop-up and a line input, over the Inputs screen.
    press(&mut a, KeyCode::Char('1'));
    press(&mut a, KeyCode::Enter);
    assert!(a.picking.is_some());
    for (w, h) in [(140, 44), (20, 6)] {
        let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
        term.draw(|f| a.render(f)).unwrap();
    }
    press(&mut a, KeyCode::Esc);
    a.inputs.row = 3;
    press(&mut a, KeyCode::Enter);
    assert!(matches!(a.editing, Some((Target::Output, _))));
    let mut term = Terminal::new(TestBackend::new(20, 6)).unwrap();
    term.draw(|f| a.render(f)).unwrap();
}

#[test]
fn flags_keys_edit_the_highlighted_row() {
    let (_t, mut a) = app_with_bams();
    press(&mut a, KeyCode::Char('3'));
    assert_eq!(a.page, Page::Flags);
    let i = a
        .form
        .fields
        .iter()
        .position(|f| f.long == "max-threads")
        .unwrap();
    a.flags_at = a.visible_flags().iter().position(|&k| k == i).unwrap();
    assert_eq!(a.flag_row(), Some(i));
    press(&mut a, KeyCode::Enter);
    for c in "7".chars() {
        press(&mut a, KeyCode::Char(c));
    }
    press(&mut a, KeyCode::Enter);
    assert!(a.editing.is_none());
    assert_ne!(a.form.fields[i].value, a.form.fields[i].default);
    press(&mut a, KeyCode::Char('r'));
    assert!(a.form.fields[i].is_default());
    press(&mut a, KeyCode::Char('/'));
    for c in "zzz-none".chars() {
        press(&mut a, KeyCode::Char(c));
    }
    press(&mut a, KeyCode::Enter);
    assert!(a.visible_flags().is_empty() && a.flag_row().is_none());
    let mut term = Terminal::new(TestBackend::new(60, 16)).unwrap();
    term.draw(|f| a.render(f)).unwrap();
}

#[test]
fn q_is_refused_while_a_run_goes() {
    let (tmp, mut a) = app_with_bams();
    let fake = tmp.path().join("fake-faba");
    std::fs::write(&fake, "#!/usr/bin/env bash\nsleep 30\n").unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    a.program = fake;
    a.start().unwrap();
    press(&mut a, KeyCode::Char('q'));
    assert!(!a.quit && a.note.is_some(), "refused while running");
    // s asks, a second s stops.
    press(&mut a, KeyCode::Char('s'));
    assert!(a.job.as_ref().unwrap().asking);
    press(&mut a, KeyCode::Char('s'));
    let job = a.job.as_mut().unwrap();
    assert!(job.stopper.is_stopped());
    job.stopper.stop();
    job.handle.take().unwrap().join().unwrap();
    assert_eq!(
        job.log.lock().unwrap().ended,
        Some(Err(child::Failed::Stopped))
    );
    press(&mut a, KeyCode::Char('q'));
    assert!(a.quit, "q leaves once the run is over");
}

#[test]
fn flags_on_the_command_line_prefill_the_view() {
    use clap::FromArgMatches;
    let tmp = tempfile::tempdir().unwrap();
    std::fs::write(tmp.path().join("sample_A.bam"), b"").unwrap();
    let bam = tmp
        .path()
        .join("sample_A.bam")
        .to_string_lossy()
        .into_owned();
    let cmd = run_cmd();
    let m = cmd
        .clone()
        .try_get_matches_from([
            "run",
            &bam,
            "-g",
            "genes.gff",
            "--max-threads",
            "3",
            "--skip-apa",
        ])
        .unwrap();
    let args = crate::pipeline::args::PipelineArgs::from_arg_matches(&m).unwrap();
    let mut a = App::new(cmd, tmp.path().to_path_buf());
    a.prefill(&m, &args);
    let argv = a.argv();
    assert!(argv.contains(&bam));
    assert!(argv
        .windows(2)
        .any(|w| w[0] == "-g" && w[1].ends_with("genes.gff")));
    assert!(argv.windows(2).any(|w| w == ["--max-threads", "3"]));
    assert!(argv.contains(&"--skip-apa".to_string()));
}

#[test]
fn one_bam_under_two_spellings_is_one_pick() {
    use clap::FromArgMatches;
    let tmp = tempfile::tempdir().unwrap();
    std::fs::write(tmp.path().join("sample_A.bam"), b"").unwrap();
    std::fs::create_dir(tmp.path().join("sub")).unwrap();
    let spelled = format!("{}/sub/../sample_A.bam", tmp.path().display());
    let cmd = run_cmd();
    let m = cmd.clone().try_get_matches_from(["run", &spelled]).unwrap();
    let args = crate::pipeline::args::PipelineArgs::from_arg_matches(&m).unwrap();
    let mut a = App::new(cmd, tmp.path().to_path_buf());
    a.prefill(&m, &args);
    assert_eq!(a.inputs.picked.len(), 1);
    assert!(a.inputs.bams.list.select("sample_A.bam"));
    let shown = a.inputs.bams.cwd.join("sample_A.bam");
    assert_eq!(a.inputs.role_of(&shown), Some(Role::Fg));
    press(&mut a, KeyCode::Char(' '));
    assert_eq!(
        a.inputs.role_of(&shown),
        Some(Role::Bg),
        "Space makes it bg"
    );
    press(&mut a, KeyCode::Char(' '));
    assert!(a.inputs.picked.is_empty(), "then removes it");
    press(&mut a, KeyCode::Char(' '));
    assert_eq!(a.inputs.picked.len(), 1);
}

/// The screen as text, drawn at `w` x `h`.
fn screen(a: &mut App, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| a.render(f)).unwrap();
    term.backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect()
}

/// A stand-in for faba running `body`, as the app's program.
fn fake_program(a: &mut App, dir: &std::path::Path, body: &str) {
    let fake = dir.join("fake-faba");
    std::fs::write(&fake, format!("#!/usr/bin/env bash\n{body}\n")).unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    a.program = fake;
}

#[test]
fn the_footer_offers_q_only_when_it_works() {
    let (tmp, mut a) = app_with_bams();
    fake_program(&mut a, tmp.path(), "exec sleep 30");
    a.start().unwrap();
    for page in ['4', '1'] {
        press(&mut a, KeyCode::Char(page));
        let s = screen(&mut a, 160, 20);
        assert!(!s.contains("q quit"), "no `q quit` while running");
    }
    press(&mut a, KeyCode::Char('4'));
    assert!(screen(&mut a, 160, 20).contains("refused"));
    // Asked to stop, then the run ends some other way: the question goes.
    press(&mut a, KeyCode::Char('s'));
    assert!(a.job.as_ref().unwrap().asking);
    let job = a.job.as_mut().unwrap();
    job.stopper.stop();
    job.stopper.stop();
    job.handle.take().unwrap().join().unwrap();
    let s = screen(&mut a, 160, 20);
    assert!(!a.job.as_ref().unwrap().asking);
    assert!(!s.contains("again stops") && s.contains("quit"), "{s}");
}

#[test]
fn the_log_scrolls_no_further_than_its_top() {
    let (tmp, mut a) = app_with_bams();
    fake_program(
        &mut a,
        tmp.path(),
        "for i in $(seq 1 30); do echo \"line $i\" >&2; done",
    );
    a.start().unwrap();
    a.job
        .as_mut()
        .unwrap()
        .handle
        .take()
        .unwrap()
        .join()
        .unwrap();
    let s = screen(&mut a, 80, 12);
    assert!(s.contains("line 30"));
    for _ in 0..100 {
        press(&mut a, KeyCode::Up);
    }
    let s = screen(&mut a, 80, 12);
    assert!(s.contains("line 1 ") && !s.contains("line 30"), "{s}");
    // One Down moves at once: no presses lost past the top.
    press(&mut a, KeyCode::Down);
    let s = screen(&mut a, 80, 12);
    assert!(!s.contains("line 1 ") && s.contains("line 2 "), "{s}");
}

#[test]
fn the_preview_shows_problems_above_a_long_command() {
    let (tmp, mut a) = app_with_bams();
    for k in 0..60 {
        let path = tmp.path().join(format!("extra_{k:02}.bam"));
        std::fs::write(&path, b"").unwrap();
        a.inputs.picked.push(Picked {
            path,
            role: Role::Fg,
        });
    }
    a.inputs.gff = None;
    a.preview = true;
    // Problems are checked after each key, not at each draw.
    a.refresh();
    let s = screen(&mut a, 120, 24);
    assert!(s.contains("GFF") && s.contains("saved as"), "{s}");
    a.inputs.gff = Some(tmp.path().join("genes.gff"));
    a.refresh();
    let s = screen(&mut a, 120, 24);
    assert!(s.contains("no problems") && s.contains("saved as"), "{s}");
}

#[test]
fn the_check_is_cached_until_a_key() {
    let (_t, mut a) = app_with_bams();
    a.refresh();
    assert!(a.checked.problems.is_empty(), "{:?}", a.checked.problems);
    a.inputs.gff = None;
    assert!(a.checked.problems.is_empty(), "a draw does not re-check");
    press(&mut a, KeyCode::Char('3'));
    assert!(a.checked.problems.iter().any(|p| p.contains("GFF")));
}

#[test]
fn verbose_is_passed_on_to_the_run_and_its_script() {
    use clap::FromArgMatches;
    let (tmp, mut a) = app_with_bams();
    let bam = tmp.path().join("sample_A.bam");
    let m = run_cmd()
        .try_get_matches_from(["run", "-v", &bam.to_string_lossy()])
        .unwrap();
    let args = crate::pipeline::args::PipelineArgs::from_arg_matches(&m).unwrap();
    let mut b = App::new(run_cmd(), tmp.path().to_path_buf());
    b.prefill(&m, &args);
    b.inputs.picked.clone_from(&a.inputs.picked);
    b.inputs.gff = a.inputs.gff.take();
    b.inputs.genome = a.inputs.genome.take();
    b.inputs.output = std::mem::take(&mut a.inputs.output);
    assert_eq!(b.argv().last().map(String::as_str), Some("-v"));
    assert!(b.problems().is_empty(), "{:?}", b.problems());
    fake_program(&mut b, tmp.path(), "echo \"$@\" > args.txt");
    b.start().unwrap();
    b.job
        .as_mut()
        .unwrap()
        .handle
        .take()
        .unwrap()
        .join()
        .unwrap();
    let out = PathBuf::from(b.inputs.output());
    let script = std::fs::read_to_string(out.join(script::SCRIPT)).unwrap();
    assert!(script.lines().any(|l| l.trim() == "-v"), "{script}");
    let given = std::fs::read_to_string(out.join("args.txt")).unwrap();
    assert!(given.trim_end().ends_with(" -v"), "{given}");
}

#[test]
fn p_prints_the_command_and_leaves() {
    let (tmp, mut a) = app_with_bams();
    let out = tmp.path().join("my out");
    a.inputs.output = out.to_string_lossy().into_owned();
    press(&mut a, KeyCode::Char('3'));
    press(&mut a, KeyCode::Char('p'));
    assert!(a.quit && a.job.is_none());
    assert!(!out.exists(), "nothing is written");
    let d = tmp.path().display();
    let want = format!(
        "faba run --batch-process {d}/sample_A.bam --control-bam {d}/control_A.bam \
         -g {d}/genes.gff -f {d}/genome.fa -o '{d}/my out'"
    );
    assert_eq!(a.printed.as_deref(), Some(want.as_str()));
}

#[test]
fn p_is_ignored_while_a_run_goes() {
    let (tmp, mut a) = app_with_bams();
    fake_program(&mut a, tmp.path(), "exec sleep 30");
    a.start().unwrap();
    for page in ['4', '1'] {
        press(&mut a, KeyCode::Char(page));
        press(&mut a, KeyCode::Char('p'));
        assert!(!a.quit && a.printed.is_none(), "page {page}");
    }
    let job = a.job.as_mut().unwrap();
    job.stopper.stop();
    job.stopper.stop();
    job.handle.take().unwrap().join().unwrap();
    press(&mut a, KeyCode::Char('p'));
    assert!(
        a.quit && a.printed.is_some(),
        "p works once the run is over"
    );
}
