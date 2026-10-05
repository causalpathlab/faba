//! `faba run`'s setup view: pick the inputs, steps and flags, preview the
//! exact command, save it as a script and run it with its log on screen.

mod child;
mod draw;
mod form;
mod inputs;
mod keys;
mod script;
mod steps;

use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use data_beans::interactive::ui::Screen;
use ratatui::crossterm::event::KeyEvent;
use ratatui::Frame;

use crate::figure::LineInput;
use crate::pipeline::args::PipelineArgs;
use crate::tui::ShiftEnter;
use child::{Failed, Progress, Said, Stopper};
use form::{Form, OWN};
use inputs::{Browser, FileRow, Inputs, Picked, Role};
use steps::Steps;

/// Lines of the child's log kept for the Run screen.
const LOG_CAP: usize = 2000;

/// The screens, in tab order.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Page {
    Inputs,
    Steps,
    Flags,
    Run,
}

/// What an open line input fills in.
pub enum Target {
    Output,
    DepthKb,
    /// A row of the form, by index.
    Flag(usize),
    Find,
}

/// What the running child has said so far.
#[derive(Default)]
pub struct Log {
    /// The last [`LOG_CAP`] lines.
    pub lines: VecDeque<String>,
    /// Its progress bar's latest frame.
    pub progress: Option<Progress>,
    /// The last `Step ...` line.
    pub step: String,
    /// How it ended, once it has.
    pub ended: Option<Result<(), Failed>>,
}

/// A run under way, or over.
pub struct Job {
    pub log: Arc<Mutex<Log>>,
    pub stopper: Arc<Stopper>,
    pub handle: Option<std::thread::JoinHandle<()>>,
    pub started: Instant,
    pub out: PathBuf,
    /// Lines from the end; 0 follows the tail.
    pub scroll: usize,
    /// `s` pressed once: asking.
    pub asking: bool,
    /// Log lines the Run screen showed when last drawn: how far up it scrolls.
    pub rows: std::cell::Cell<usize>,
}

impl Job {
    /// The child's thread is still going.
    pub fn running(&self) -> bool {
        self.handle.as_ref().is_some_and(|h| !h.is_finished())
    }
}

/// What keeps the run from starting, worked out after each key rather than
/// at each draw: clap's check parses the whole command.
#[derive(Default)]
pub struct Checked {
    /// clap's complaint about the command, if any.
    pub complaint: Option<String>,
    /// Everything that keeps the run from starting.
    pub problems: Vec<String>,
}

pub struct App {
    pub page: Page,
    pub inputs: Inputs,
    pub steps: Steps,
    pub form: Form,
    pub run_cmd: clap::Command,
    pub preview: bool,
    pub editing: Option<(Target, LineInput)>,
    /// A file row's pop-up browser: the row and the browser.
    pub picking: Option<(FileRow, Browser)>,
    /// The highlighted row among the visible flags.
    pub flags_at: usize,
    pub advanced: bool,
    pub find: String,
    pub note: Option<String>,
    pub job: Option<Job>,
    pub shift_enter: ShiftEnter,
    /// The program the run starts: this binary, or a stand-in in tests.
    pub program: PathBuf,
    pub quit: bool,
    /// `faba -v`: passed on to the run.
    pub verbose: bool,
    /// The run's end has been drawn.
    end_drawn: bool,
    /// [`App::check`] as of the last key.
    pub checked: Checked,
}

impl App {
    pub fn new(run_cmd: clap::Command, cwd: PathBuf) -> App {
        // The browser's `..` needs an absolute path to climb.
        let cwd = std::path::absolute(&cwd).unwrap_or(cwd);
        let mut app = App {
            page: Page::Inputs,
            inputs: Inputs::new(cwd),
            steps: Steps::default(),
            form: Form::new(&run_cmd, OWN),
            run_cmd,
            preview: false,
            editing: None,
            picking: None,
            flags_at: 0,
            advanced: false,
            find: String::new(),
            note: None,
            job: None,
            shift_enter: ShiftEnter::default(),
            program: std::env::current_exe().unwrap_or_else(|_| "faba".into()),
            quit: false,
            verbose: false,
            end_drawn: false,
            checked: Checked::default(),
        };
        app.refresh();
        app
    }

    /// Take what the command line gave: its flags, steps and inputs.
    pub fn prefill(&mut self, m: &clap::ArgMatches, args: &PipelineArgs) {
        self.form.prefill(m);
        self.steps.prefill(args);
        // `faba`'s own flag, present when the command is built inside `faba`.
        self.verbose = matches!(m.try_get_one::<bool>("verbose"), Ok(Some(true)));
        let controls: Vec<Picked> = args
            .control_bam_files
            .iter()
            .map(|b| Picked::new(Path::new(&**b), Role::Bg))
            .collect();
        for b in &args.bam_files {
            let p = Picked::new(Path::new(&**b), Role::Fg);
            let is_control = controls.iter().any(|c| c.path == p.path);
            if !is_control && self.inputs.role_of(&p.path).is_none() {
                self.inputs.picked.push(p);
            }
        }
        for c in controls {
            if self.inputs.role_of(&c.path).is_none() {
                self.inputs.picked.push(c);
            }
        }
        for (row, given) in [
            (FileRow::Gff, &args.gff_file),
            (FileRow::Genome, &args.genome_file),
            (FileRow::KnownSnps, &args.known_snps),
        ] {
            self.inputs.set(row, given.as_deref().map(Path::new));
        }
        if let Some(o) = &args.output {
            self.inputs.output = o.to_string();
        }
        self.refresh();
    }

    /// Indices of the form rows shown: advanced ones only when asked for,
    /// and only those matching the find text.
    pub fn visible_flags(&self) -> Vec<usize> {
        self.form
            .fields
            .iter()
            .enumerate()
            .filter(|(_, f)| (!f.advanced || self.advanced) && f.long.contains(&self.find))
            .map(|(i, _)| i)
            .collect()
    }

    /// The highlighted flag's row in the form.
    pub fn flag_row(&self) -> Option<usize> {
        let v = self.visible_flags();
        v.get(self.flags_at.min(v.len().saturating_sub(1))).copied()
    }

    /// The command after the program, run in the output folder.
    pub fn argv(&self) -> Vec<String> {
        let has_bg = !self.inputs.bg().is_empty();
        let mut v: Vec<String> = vec!["run".into(), "--batch-process".into()];
        v.extend(self.inputs.argv());
        v.extend(["-o".into(), ".".into()]);
        v.extend(self.steps.argv(has_bg));
        v.extend(self.form.argv());
        if self.verbose {
            v.push("-v".into());
        }
        v
    }

    /// clap's complaint and everything that keeps the run from starting,
    /// as things stand now.
    pub fn check(&self) -> Checked {
        let complaint = form::check(&self.run_cmd, &self.argv()).err();
        let mut problems = self.inputs.problems();
        problems.extend(self.steps.problems());
        if let Some(c) = &complaint {
            problems.push(match form::blamed(c, &self.form.fields) {
                Some(l) => format!("--{l}: {c}"),
                None => c.clone(),
            });
        }
        Checked {
            complaint,
            problems,
        }
    }

    /// Everything that keeps the run from starting, worked out afresh.
    pub fn problems(&self) -> Vec<String> {
        self.check().problems
    }

    /// Bring [`App::checked`] up to date.
    pub fn refresh(&mut self) {
        self.checked = self.check();
    }

    /// A run is going.
    pub fn running(&self) -> bool {
        self.job.as_ref().is_some_and(Job::running)
    }

    /// Save the script in the output folder and start the run there.
    pub fn start(&mut self) -> anyhow::Result<()> {
        if self.running() {
            self.note = Some("a run is going already".into());
            return Ok(());
        }
        if let Some(p) = self.problems().first() {
            self.note = Some(format!("cannot start: {p}"));
            return Ok(());
        }
        let out = PathBuf::from(self.inputs.output());
        // Pin it: the suggestion moves on once the folder exists.
        self.inputs.output = out.to_string_lossy().into_owned();
        std::fs::create_dir_all(&out).map_err(|e| anyhow::anyhow!("{}: {e}", out.display()))?;
        let argv = self.argv();
        script::write(&out, &argv)?;

        let log = Arc::new(Mutex::new(Log::default()));
        let stopper = Arc::new(Stopper::default());
        let mut cmd = std::process::Command::new(&self.program);
        cmd.args(&argv).current_dir(&out).env(
            "RUST_LOG",
            std::env::var("RUST_LOG").unwrap_or_else(|_| "info".into()),
        );
        let handle = {
            let (log, stopper) = (Arc::clone(&log), Arc::clone(&stopper));
            std::thread::spawn(move || {
                let result = child::run_one(cmd, &stopper, |said| {
                    let Ok(mut log) = log.lock() else { return };
                    match said {
                        Said::Line(l) => {
                            if l.starts_with("Step ") {
                                log.step.clone_from(&l);
                            }
                            log.lines.push_back(l);
                            while log.lines.len() > LOG_CAP {
                                log.lines.pop_front();
                            }
                        }
                        Said::Progress(p) => log.progress = Some(p),
                    }
                });
                if let Ok(mut log) = log.lock() {
                    log.ended = Some(result);
                }
            })
        };
        self.job = Some(Job {
            log,
            stopper,
            handle: Some(handle),
            started: Instant::now(),
            out,
            scroll: 0,
            asking: false,
            rows: std::cell::Cell::new(0),
        });
        self.end_drawn = false;
        self.page = Page::Run;
        self.preview = false;
        Ok(())
    }

    /// `q` or Ctrl-C: leave, unless a run is going.
    fn leave(&mut self) {
        if self.running() {
            self.note = Some("a run is going: stop it with s first".into());
        } else {
            self.quit = true;
            self.shift_enter.release();
        }
    }

    /// Put `text` on the terminal's clipboard (OSC 52).
    fn copy(&mut self, text: &str) {
        let _ = crossterm::execute!(
            std::io::stdout(),
            crossterm::clipboard::CopyToClipboard::to_clipboard_from(text)
        );
        self.note = Some("copied the command to the clipboard".into());
    }
}

impl Screen for App {
    fn render(&mut self, frame: &mut Frame) {
        self.shift_enter.arm();
        if let Some(job) = self.job.as_mut().filter(|j| !j.running()) {
            // Nothing is left to stop.
            job.asking = false;
            self.end_drawn = true;
        }
        draw::draw(self, frame);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.key(key);
    }

    fn interrupt(&mut self) {
        self.leave();
    }

    fn done(&self) -> bool {
        self.quit
    }

    fn tick(&mut self) -> bool {
        match &self.job {
            Some(j) => j.running() || !self.end_drawn,
            None => false,
        }
    }
}

/// Open the setup view in `faba run`'s flags, pre-filled from a command
/// line when one is given.
pub(crate) fn run_view(
    run_cmd: clap::Command,
    prefill: Option<(&clap::ArgMatches, &PipelineArgs)>,
) -> anyhow::Result<()> {
    let mut app = App::new(run_cmd, std::env::current_dir()?);
    if let Some((m, args)) = prefill {
        app.prefill(m, args);
    }
    app.shift_enter = ShiftEnter::wanted();
    let shown = data_beans::interactive::ui::run_screen(&mut app);
    app.shift_enter.release();
    if let Some(job) = app.job.as_mut() {
        if job.running() {
            job.stopper.stop();
            job.stopper.stop();
        }
        if let Some(h) = job.handle.take() {
            let _ = h.join();
        }
    }
    shown
}

/// `faba`'s built `run` subcommand, as `main` hands it to the view.
#[cfg(test)]
pub(crate) fn run_cmd() -> clap::Command {
    crate::faba_command()
        .find_subcommand("run")
        .cloned()
        .expect("faba has a run command")
}

#[cfg(test)]
#[path = "tui/tests/app.rs"]
mod tests;
