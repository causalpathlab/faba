//! `faba run`'s setup view: pick the inputs, steps and flags, preview the
//! exact command, save it as a script and run it with its log on screen.

pub mod child;
mod draw;
pub mod form;
pub mod inputs;
pub(crate) mod keys;
pub mod script;
pub mod steps;

use std::collections::VecDeque;
use std::path::PathBuf;
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
use inputs::{Browser, Inputs, Picked, Role};
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

pub struct App {
    pub page: Page,
    pub inputs: Inputs,
    pub steps: Steps,
    pub form: Form,
    pub run_cmd: clap::Command,
    pub preview: bool,
    pub editing: Option<(Target, LineInput)>,
    /// A file row's pop-up browser: the row and the browser.
    pub picking: Option<(usize, Browser)>,
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
}

impl App {
    pub fn new(run_cmd: clap::Command, cwd: PathBuf) -> App {
        // The browser's `..` needs an absolute path to climb.
        let cwd = std::path::absolute(&cwd).unwrap_or(cwd);
        App {
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
        }
    }

    /// Take what the command line gave: its flags, steps and inputs.
    pub fn prefill(&mut self, m: &clap::ArgMatches, args: &PipelineArgs) {
        self.form.prefill(m);
        self.steps.prefill(args);
        // `faba`'s own flag, present when the command is built inside `faba`.
        self.verbose = matches!(m.try_get_one::<bool>("verbose"), Ok(Some(true)));
        let abs = |s: &str| inputs::normalize(std::path::Path::new(s));
        let controls: Vec<PathBuf> = args.control_bam_files.iter().map(|s| abs(s)).collect();
        for b in &args.bam_files {
            let path = abs(b);
            if !controls.contains(&path) && self.inputs.role_of(&path).is_none() {
                self.inputs.picked.push(Picked {
                    path,
                    role: Role::Fg,
                });
            }
        }
        for path in controls {
            if self.inputs.role_of(&path).is_none() {
                self.inputs.picked.push(Picked {
                    path,
                    role: Role::Bg,
                });
            }
        }
        self.inputs.gff = args.gff_file.as_deref().map(abs);
        self.inputs.genome = args.genome_file.as_deref().map(abs);
        self.inputs.known_snps = args.known_snps.as_deref().map(abs);
        if let Some(o) = &args.output {
            self.inputs.output = o.to_string();
        }
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

    /// clap's complaint about the command, if any.
    pub fn clap_complaint(&self) -> Option<String> {
        form::check(&self.run_cmd, &self.argv()).err()
    }

    /// Everything that keeps the run from starting.
    pub fn problems(&self) -> Vec<String> {
        let mut v = self.inputs.problems();
        v.extend(self.steps.problems());
        if let Some(c) = self.clap_complaint() {
            v.push(match form::blamed(&c, &self.form.fields) {
                Some(l) => format!("--{l}: {c}"),
                None => c,
            });
        }
        v
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
        use std::io::Write;
        let mut out = std::io::stdout();
        let _ = write!(out, "\x1b]52;c;{}\x07", base64(text.as_bytes()));
        let _ = out.flush();
        self.note = Some("copied the command to the clipboard".into());
    }
}

/// Standard base64, padded.
fn base64(bytes: &[u8]) -> String {
    const ABC: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut s = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            chunk.get(1).copied().unwrap_or(0),
            chunk.get(2).copied().unwrap_or(0),
        ];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for k in 0..4 {
            if k <= chunk.len() {
                s.push(ABC[(n >> (18 - 6 * k) & 63) as usize] as char);
            } else {
                s.push('=');
            }
        }
    }
    s
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
pub fn run_view(
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

#[cfg(test)]
#[path = "tui/tests/app.rs"]
mod tests;
