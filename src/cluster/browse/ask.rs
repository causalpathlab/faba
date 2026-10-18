//! The form `faba modules` asks the blockmodel's settings in before it
//! fits in a terminal: a row per setting of [`ROWS`], from `faba cluster`'s
//! clap definition (see [`crate::tui::form`]), prefilled from the flags
//! given, else the preset's; checked after every key as `faba cluster`
//! checks them. Ctrl+R fits with them, Esc leaves without, Ctrl+C quits.

use super::super::args::{self, Global, Switch};
use super::super::command::{command_except, definition, value_words};
use super::super::unified::Weights;
use super::{ClusterArgs, ModelKind};
use crate::figure::{Edit, LineInput};
use crate::tui::form::{blamed, parse, row_line, Form, Kind};
use crate::tui::{is_apply, moved, window, View, APPLY_KEYS};
use clap::{Command, FromArgMatches};
use data_beans::interactive::ui::{
    header, help_line, input_line, panel, Screen, ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::Line;
use ratatui::widgets::{Paragraph, Wrap};
use ratatui::Frame;

/// The settings asked, by long flag, in the order shown.
pub(super) const ROWS: [&str; 16] = [
    "model",
    "cell-roles",
    "gene-roles",
    "chains",
    "sweeps",
    "burn-in",
    "thin",
    "seed",
    "profile-strength",
    "global",
    "shift",
    "cell-sd",
    "membership-alpha",
    "gene-alpha",
    "shift-sd",
    "holdout",
];

/// Longest text a row takes.
const MAX_TYPED: usize = 64;

/// The rows a preset sets, each with its value there: what `base` fits as
/// `model` with none of them given.
fn preset(base: &ClusterArgs, model: ModelKind) -> [(&'static str, String); 4] {
    let s = args::blockmodel(&ClusterArgs {
        model,
        gene_roles: None,
        shift: None,
        global: None,
        profile_strength: None,
        ..base.clone()
    });
    let shift = if s.shift_sd.is_some() {
        Switch::On
    } else {
        Switch::Off
    };
    let global = match s.weights {
        Weights::Hdp(_) => Global::Hdp,
        Weights::Symmetric => Global::Symmetric,
    };
    [
        ("gene-roles", s.gene_roles.to_string()),
        ("shift", name(shift)),
        ("global", name(global)),
        ("profile-strength", s.strength.to_string()),
    ]
}

/// A value-enum's name, as its flag takes it.
fn name(v: impl clap::ValueEnum) -> String {
    v.to_possible_value()
        .map(|p| p.get_name().to_string())
        .unwrap_or_default()
}

/// What the form was left with.
pub(super) enum Asked {
    /// Ctrl+R: fit with these.
    Fit(Box<ClusterArgs>),
    /// Esc or `q`: no fit; on to what the footer named.
    Back,
    /// Ctrl+C: no fit, and nothing after.
    Quit,
}

/// Whether `a` and `b` read as the same number.
fn same_number(a: &str, b: &str) -> bool {
    matches!((a.parse::<f64>(), b.parse::<f64>()), (Ok(x), Ok(y)) if x == y)
}

/// The form and what came of it.
pub(super) struct FitForm {
    /// The settings the form starts from: the inputs and every option not
    /// asked are taken from here as they are.
    base: ClusterArgs,
    /// `faba cluster`'s definition, built: the rows' source and their check.
    cmd: Command,
    form: Form,
    /// Per row: given on the command line, or changed here; a preset
    /// switch refills only the rows that are not.
    touched: Vec<bool>,
    at: usize,
    /// The row being typed into.
    editing: Option<LineInput>,
    /// The settings the rows stand for, or why they cannot be fitted.
    checked: Result<ClusterArgs, String>,
    title: String,
    /// Why the fit is to run.
    why: String,
    note: Option<String>,
    /// What Esc leaves to, as the footer says it.
    leave: String,
    outcome: Option<Asked>,
}

impl FitForm {
    /// The form over `base` (a blockmodel's settings), under `title`, the
    /// fit asked for because of `why`; Esc leaves to `leave`.
    pub(super) fn new(base: &ClusterArgs, title: &str, why: &str, leave: &str) -> Self {
        let mut cmd = definition();
        cmd.build();
        let mut form = Form::only(&cmd, &ROWS);
        let preset = preset(base, base.model);
        let mut touched = Vec::new();
        for f in &mut form.fields {
            if let Kind::Choice(values) = &mut f.kind {
                // Unset leaves it to the preset, which the row shows; the
                // mixture writes no modules.
                values.retain(|v| !v.is_empty() && *v != ModelKind::Crp.name());
            }
            if f.long == "model" {
                f.default = ModelKind::Blocks.name();
            }
            if let Some((_, v)) = preset.iter().find(|(l, _)| *l == f.long) {
                f.default.clone_from(v);
            }
            let given = base.given(&f.id);
            f.value = if f.long == "model" {
                base.model.name()
            } else if given {
                let words = value_words(base, &f.id).join(",");
                // `1.0` given where the default reads `1` is the default.
                if same_number(&words, &f.default) {
                    f.default.clone()
                } else {
                    words
                }
            } else {
                f.default.clone()
            };
            touched.push(given);
        }
        let mut view = FitForm {
            base: base.clone(),
            cmd,
            form,
            touched,
            at: 0,
            editing: None,
            checked: Err(String::new()),
            title: title.to_string(),
            why: why.to_string(),
            note: None,
            leave: leave.to_string(),
            outcome: None,
        };
        view.refresh();
        view
    }

    /// The value of row `long`, as it stands.
    #[cfg(test)]
    pub(super) fn value(&self, long: &str) -> &str {
        self.form.get(long).map_or("", |f| f.value.as_str())
    }

    /// The settings Ctrl+R took; none when the form was left without.
    #[cfg(test)]
    pub(super) fn chosen(&self) -> Option<ClusterArgs> {
        match &self.outcome {
            Some(Asked::Fit(args)) => Some(*args.clone()),
            _ => None,
        }
    }

    /// What the form was left with; none while it is open.
    #[cfg(test)]
    pub(super) fn outcome(&self) -> Option<&Asked> {
        self.outcome.as_ref()
    }

    /// The settings `form`'s rows stand for, `faba cluster`'s checks
    /// passed; the complaint if not.
    fn settings(&self, form: &Form) -> Result<ClusterArgs, String> {
        let mut argv = command_except(&self.base, &ROWS);
        for (f, &touched) in form.fields.iter().zip(&self.touched) {
            // The model always: `faba cluster`'s default is not ours.
            if touched || !f.is_default() || f.long == "model" {
                argv.extend(f.line());
            }
        }
        let m = parse(&self.cmd, &argv)?;
        let mut args = ClusterArgs::from_arg_matches(&m).map_err(|e| e.to_string())?;
        args.record_given(&m);
        args::check(&args).map_err(|e| format!("{e:#}"))?;
        Ok(args)
    }

    /// Check the rows again, with the text being typed in its row.
    fn refresh(&mut self) {
        let typed = self.editing.as_ref().and_then(LineInput::text);
        self.checked = match typed {
            Some(text) => {
                let mut form = self.form.clone();
                form.fields[self.at].set_text(text);
                self.settings(&form)
            }
            None => self.settings(&self.form),
        };
    }

    /// The model row was changed: the rows its preset sets follow it,
    /// unless given or changed here.
    fn follow_preset(&mut self) {
        let model = match self.form.get("model").map(|f| f.value.as_str()) {
            Some(v) if v == ModelKind::Link.name() => ModelKind::Link,
            _ => ModelKind::Blocks,
        };
        for (long, value) in preset(&self.base, model) {
            let Some(i) = self.form.fields.iter().position(|f| f.long == long) else {
                continue;
            };
            let f = &mut self.form.fields[i];
            if !self.touched[i] {
                f.value.clone_from(&value);
            }
            f.default = value;
        }
    }

    /// Row `at`'s choice stepped by `d`.
    fn cycle(&mut self, d: isize) {
        let f = &mut self.form.fields[self.at];
        if !matches!(f.kind, Kind::Choice(_)) {
            return;
        }
        f.cycle(d);
        self.touched[self.at] = true;
        if f.long == "model" {
            self.follow_preset();
        }
    }

    /// The text being typed, set in its row.
    fn commit(&mut self, text: &str) {
        self.editing = None;
        self.form.fields[self.at].set_text(text);
        self.touched[self.at] = true;
    }

    fn key_editing(&mut self, key: KeyEvent) {
        let Some(line) = self.editing.as_mut() else {
            return;
        };
        match line.handle(key) {
            Edit::Typing => {}
            Edit::Cancelled => self.editing = None,
            Edit::Submitted(text) => self.commit(&text),
        }
    }

    fn key_rows(&mut self, key: KeyEvent) {
        let last = self.form.fields.len().saturating_sub(1);
        let nav = match key.code {
            KeyCode::Up | KeyCode::Char('k') => Some(-1),
            KeyCode::Down | KeyCode::Char('j') => Some(1),
            KeyCode::PageUp => Some(-10),
            KeyCode::PageDown => Some(10),
            _ => None,
        };
        if let Some(d) = nav {
            self.at = moved(self.at, d, last);
            return;
        }
        let choice = matches!(self.form.fields[self.at].kind, Kind::Choice(_));
        match key.code {
            KeyCode::Home => self.at = 0,
            KeyCode::End => self.at = last,
            KeyCode::Left => self.cycle(-1),
            KeyCode::Right | KeyCode::Char(' ') => self.cycle(1),
            KeyCode::Enter if choice => self.cycle(1),
            KeyCode::Enter => {
                let mut line = LineInput::new(MAX_TYPED);
                line.open(&self.form.fields[self.at].value);
                self.editing = Some(line);
            }
            KeyCode::Char('r') => {
                self.form.fields[self.at].reset();
                self.touched[self.at] = false;
                if self.form.fields[self.at].long == "model" {
                    self.follow_preset();
                }
            }
            KeyCode::Esc | KeyCode::Char('q') => self.outcome = Some(Asked::Back),
            _ => {}
        }
    }

    /// The rows, the one in focus marked, a reason beside the row it
    /// blames.
    fn rows(&self) -> Vec<Line<'static>> {
        let complaint = self.checked.as_ref().err();
        let blamed = complaint.and_then(|c| blamed(c, &self.form.fields));
        let typed = self.editing.as_ref().and_then(LineInput::text);
        self.form
            .fields
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let shown = match (typed, &f.kind) {
                    (Some(t), _) if i == self.at => format!("{t}▏"),
                    (_, Kind::Choice(_)) => format!("‹ {} ›", f.shown()),
                    _ => f.shown(),
                };
                let said = complaint
                    .map(String::as_str)
                    .filter(|_| blamed == Some(f.long.as_str()));
                row_line(f, i == self.at, 20, PLAIN, shown, said)
            })
            .collect()
    }

    /// Whether the fit can start, or why not.
    fn status(&self) -> Line<'static> {
        if let Some(note) = &self.note {
            return Line::styled(note.clone(), ACCENTED);
        }
        match &self.checked {
            Ok(args) => Line::styled(
                format!(
                    "{APPLY_KEYS} fits the {} model with these settings",
                    args.model.name()
                ),
                DIM,
            ),
            Err(why) => Line::styled(format!("cannot fit: {why}"), ACCENTED),
        }
    }

    fn footer(&self) -> Line<'static> {
        if let Some(text) = self.editing.as_ref().and_then(LineInput::text) {
            let prompt = format!("--{} ", self.form.fields[self.at].long);
            return input_line(&prompt, text, &[("Enter", "set"), ("Esc", "cancel")]);
        }
        let fit = if self.checked.is_ok() {
            "fit"
        } else {
            "fit (blocked)"
        };
        help_line(&[
            ("↑/↓", "move"),
            ("←/→/Space", "choice"),
            ("Enter", "type"),
            ("r", "reset"),
            (APPLY_KEYS, fit),
            ("Esc", &self.leave),
            ("Ctrl+C", "quit"),
        ])
    }
}

impl View for FitForm {
    fn stray_enter(&mut self) {
        self.note = Some(format!(
            "only a plain Enter acts here: Enter types or steps a row; {APPLY_KEYS} fits"
        ));
    }

    fn takes_apply(&self) -> bool {
        true
    }
}

impl Screen for FitForm {
    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        let help = crate::tui::wrap_help(self.footer(), area.width);
        let [top, body, foot] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(help.len() as u16),
        ])
        .areas(area);
        frame.render_widget(header("modules", &self.title, &self.why), top);
        let changed = self.form.changed();
        let block = panel(format!(" fit settings · {changed} changed "), true);
        let inner = block.inner(body);
        frame.render_widget(block, body);
        let [list, about, status] = Layout::vertical([
            Constraint::Fill(1),
            Constraint::Length(5),
            Constraint::Length(2),
        ])
        .areas(inner);
        let rows = list.height as usize;
        let lines = self.rows();
        frame.render_widget(Paragraph::new(window(lines, self.at, rows)), list);
        if let Some(f) = self.form.fields.get(self.at) {
            let mut text = vec![Line::styled(format!("--{}", f.long), HIGHLIGHT)];
            text.extend(f.long_help.lines().map(|l| Line::raw(l.trim().to_string())));
            frame.render_widget(Paragraph::new(text).wrap(Wrap { trim: true }), about);
        }
        frame.render_widget(
            Paragraph::new(self.status()).wrap(Wrap { trim: true }),
            status,
        );
        frame.render_widget(Paragraph::new(help), foot);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.note = None;
        if is_apply(&key) {
            // The text being typed counts: set it, then fit.
            if let Some(text) = self.editing.as_ref().and_then(LineInput::text) {
                let text = text.to_string();
                self.commit(&text);
                self.refresh();
            }
            match &self.checked {
                Ok(args) => self.outcome = Some(Asked::Fit(Box::new(args.clone()))),
                Err(why) => self.note = Some(format!("{APPLY_KEYS} is blocked: {why}")),
            }
        } else if self.editing.is_some() {
            self.key_editing(key);
        } else {
            self.key_rows(key);
        }
        self.refresh();
    }

    fn interrupt(&mut self) {
        self.outcome = Some(Asked::Quit);
    }

    fn done(&self) -> bool {
        self.outcome.is_some()
    }
}

/// Ask the settings to fit `base` with in the form, under `title`, the fit
/// asked for because of `why`, Esc leaving to `leave`.
pub(super) fn ask(
    base: &ClusterArgs,
    title: &str,
    why: &str,
    leave: &str,
) -> anyhow::Result<Asked> {
    let mut view = FitForm::new(base, title, why, leave);
    crate::tui::run_view(&mut view)?;
    Ok(view.outcome.unwrap_or(Asked::Quit))
}

/// The settings to fit: `args` as they are unless `interactive`, else what
/// `ask` asks.
pub(super) fn chosen(
    args: ClusterArgs,
    interactive: bool,
    ask: impl FnOnce(&ClusterArgs) -> anyhow::Result<Asked>,
) -> anyhow::Result<Asked> {
    if interactive {
        ask(&args)
    } else {
        Ok(Asked::Fit(Box::new(args)))
    }
}
