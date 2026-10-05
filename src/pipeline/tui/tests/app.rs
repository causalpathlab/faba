use super::*;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use ratatui::Terminal;

fn run_cmd() -> clap::Command {
    use clap::{CommandFactory, Parser};
    #[derive(Parser)]
    #[command(name = "run")]
    struct Run {
        #[command(flatten)]
        a: crate::pipeline::args::PipelineArgs,
    }
    let mut c = Run::command();
    c.build();
    c
}

fn press(a: &mut App, code: KeyCode) {
    a.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn app_with_bams() -> (tempfile::TempDir, App) {
    let tmp = tempfile::tempdir().unwrap();
    for f in ["sample_A.bam", "control_A.bam", "genes.gff", "genome.fa"] {
        std::fs::write(tmp.path().join(f), b"").unwrap();
    }
    let mut a = App::new(run_cmd(), tmp.path().to_path_buf());
    for n in ["control_A.bam", "sample_A.bam"] {
        a.inputs.bams.at = a
            .inputs
            .bams
            .entries
            .iter()
            .position(|e| e.name == n)
            .unwrap();
        press(&mut a, KeyCode::Char(' '));
    }
    a.inputs.bams.at = a
        .inputs
        .bams
        .entries
        .iter()
        .position(|e| e.name == "control_A.bam")
        .unwrap();
    press(&mut a, KeyCode::Char('b'));
    a.inputs.gff = Some(tmp.path().join("genes.gff"));
    a.inputs.genome = Some(tmp.path().join("genome.fa"));
    (tmp, a)
}

#[test]
fn tab_walks_the_screens_but_run_needs_a_job() {
    let (_t, mut a) = app_with_bams();
    for want in [Page::Steps, Page::Flags, Page::Inputs] {
        press(&mut a, KeyCode::Tab);
        assert_eq!(a.page, want);
    }
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
    a.handle_key(KeyEvent::new(KeyCode::Enter, KeyModifiers::SHIFT));
    assert!(a.preview);
    press(&mut a, KeyCode::Esc);
    assert!(!a.preview);
    a.inputs.gff = None;
    assert!(a.problems().iter().any(|p| p.contains("GFF")));
    press(&mut a, KeyCode::Char('G'));
    press(&mut a, KeyCode::Char('y'));
    assert!(a.job.is_none(), "problems block the start");
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
    press(&mut a, KeyCode::Char('l'));
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
fn base64_matches_the_standard() {
    assert_eq!(base64(b""), "");
    assert_eq!(base64(b"f"), "Zg==");
    assert_eq!(base64(b"fo"), "Zm8=");
    assert_eq!(base64(b"foo"), "Zm9v");
    assert_eq!(base64(b"foobar"), "Zm9vYmFy");
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
    let at = a
        .inputs
        .bams
        .entries
        .iter()
        .position(|e| e.name == "sample_A.bam")
        .unwrap();
    a.inputs.bams.at = at;
    let shown = a.inputs.bams.cwd.join("sample_A.bam");
    assert_eq!(a.inputs.role_of(&shown), Some(Role::Fg));
    press(&mut a, KeyCode::Char(' '));
    assert!(a.inputs.picked.is_empty(), "Space removes it");
    press(&mut a, KeyCode::Char(' '));
    press(&mut a, KeyCode::Char(' '));
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
    let s = screen(&mut a, 120, 24);
    assert!(s.contains("GFF") && s.contains("saved as"), "{s}");
    a.inputs.gff = Some(tmp.path().join("genes.gff"));
    let s = screen(&mut a, 120, 24);
    assert!(s.contains("no problems") && s.contains("saved as"), "{s}");
}

#[test]
fn verbose_is_passed_on_to_the_run_and_its_script() {
    use clap::{CommandFactory, FromArgMatches};
    let (tmp, mut a) = app_with_bams();
    let bam = tmp.path().join("sample_A.bam");
    let line = ["faba", "run", "-v", &bam.to_string_lossy()].map(String::from);
    let (run_cmd, m) = crate::pipeline::run::run_matches(crate::Cli::command(), line).unwrap();
    let args = crate::pipeline::args::PipelineArgs::from_arg_matches(&m).unwrap();
    let mut b = App::new(run_cmd, tmp.path().to_path_buf());
    b.prefill(&m, &args);
    b.inputs.picked.clone_from(&a.inputs.picked);
    b.inputs.gff = a.inputs.gff.take();
    b.inputs.genome = a.inputs.genome.take();
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
