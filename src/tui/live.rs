//! A view's state, read and written from outside it: what the user
//! arranged in a live session (the knob in focus, the scales, the look),
//! for an agent to read back and to draw figures from, and the same state
//! kept as a file (`--view-state`) to open a view, or draw its figure
//! without one, as it was.
//!
//! A view that runs in a named session ([`legume_tui::Session::name`])
//! answers on a socket of its own beside the session's, `{id}.ctl` in
//! [`legume_tui::runtime_dir`]: one JSON line asked ([`Request`]), one
//! JSON line answered, `{"ok": ...}` or `{"error": "..."}`. The socket's
//! thread hands each request to the view's loop, which answers it between
//! keys ([`Control::serve`], from the view's tick), so the state read is
//! the one on screen and a figure is drawn as `s` draws it.

use std::io::{BufRead, BufReader, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::Duration;

use anyhow::Context;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::figure::size::{Figure, SizeArgs};
use crate::figure::style::{Palette, Sheet};
use crate::figure::Format;

/// What a view answers: its state as the user left it, set from a file,
/// and the figure its `s` saves.
pub trait Live {
    /// The view's name, as MCP's `view` names it (`qc_view`).
    fn view_name(&self) -> &'static str;

    /// What the user arranged, as JSON: the view's own keys (see each
    /// view's state module).
    fn view_state(&self) -> Value;

    /// Arrange the view as `state` says (as [`Live::view_state`] wrote
    /// it); a key left out stays as it is, one it does not know is passed
    /// over.
    fn set_view_state(&mut self, state: &Value) -> anyhow::Result<()>;

    /// The size and style `s` saves at, and the figure it is (whose floor
    /// a size asked for is held to).
    fn own_sheet(&self) -> (Sheet, Figure);

    /// The figure `s` saves, as SVG, laid out for `sheet`.
    fn saved_figure(&self, sheet: Sheet) -> anyhow::Result<String>;
}

/// The sheet `view`'s figure is drawn on with `look` over its own.
pub fn sheet_of(view: &(impl Live + ?Sized), look: &Look) -> anyhow::Result<Sheet> {
    let (own, figure) = view.own_sheet();
    look.over(own, figure)
}

/// `sheet` as a [`Look`] says it: each side in inches (none, the
/// figure's own), the font (none, fit), the plot alone, the palette.
pub fn sheet_json(sheet: &Sheet) -> Value {
    let inches = |p: Option<f64>| p.map(|p| p / crate::figure::size::PT);
    json!({
        "width": inches(sheet.page.0),
        "height": inches(sheet.page.1),
        "font_size": sheet.look.font_size,
        "panel_only": sheet.look.panel_only,
        "palette": sheet.look.palette,
    })
}

/// A figure's size and style asked for over a view's own: each side in
/// inches, the font in points, the plot alone, the palette; what is not
/// given is the view's.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Look {
    pub width: Option<f64>,
    pub height: Option<f64>,
    pub font_size: Option<f64>,
    pub panel_only: Option<bool>,
    pub palette: Option<Palette>,
}

impl Look {
    /// `own` with what is asked for over it, checked as the command line
    /// checks it for `figure`.
    pub fn over(&self, own: Sheet, figure: Figure) -> anyhow::Result<Sheet> {
        let mut look = own.look;
        look.font_size = self.font_size.or(look.font_size);
        look.panel_only = self.panel_only.unwrap_or(look.panel_only);
        look.palette = self.palette.unwrap_or(look.palette);
        look.check()?;
        let size = SizeArgs {
            fig_width: self.width,
            fig_height: self.height,
        };
        let asked = figure.page(&size)?;
        let page = (asked.0.or(own.page.0), asked.1.or(own.page.1));
        Ok(Sheet { page, look })
    }
}

/// What is asked of a live view.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
pub enum Request {
    /// Its state ([`Live::view_state`]), with the view's name and the
    /// command line it was opened with.
    State,
    /// Arrange it as `state` says ([`Live::set_view_state`]); its state
    /// after.
    Set { state: Value },
    /// Write the figure `s` saves, at `look`, as `{prefix}.{format}` for
    /// each of `formats`, and its state beside it ([`state_file`]); none
    /// of them may exist.
    Export {
        prefix: String,
        formats: Vec<Format>,
        #[serde(default)]
        look: Look,
    },
}

/// The file a figure's state is kept in beside it: `{prefix}.view.json`.
pub fn state_file(prefix: &str) -> String {
    format!("{prefix}.view.json")
}

/// A `--view-state` file: a JSON object as [`Request::State`] answers.
pub fn read_state(path: &str) -> anyhow::Result<Value> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let v: Value = serde_json::from_str(&text).with_context(|| format!("{path}: not JSON"))?;
    anyhow::ensure!(v.is_object(), "{path}: a view's state is a JSON object");
    Ok(v)
}

/// `state` checked to be for `view`: the name it was written by, if it
/// names one, is `view`'s.
pub fn check_view(state: &Value, view: &str) -> anyhow::Result<()> {
    match state.get("view").and_then(Value::as_str) {
        Some(v) if v != view => anyhow::bail!("this state is {v}'s, not {view}'s"),
        _ => Ok(()),
    }
}

/// `view`'s state with its name and the command line it runs with (the
/// arguments after the program's name), so a file of it says what it is
/// and what it was opened on.
pub fn snapshot(view: &(impl Live + ?Sized)) -> Value {
    let mut state = view.view_state();
    if let Some(m) = state.as_object_mut() {
        m.insert("view".into(), json!(view.view_name()));
        let command: Vec<String> = std::env::args().skip(1).collect();
        m.insert("command".into(), json!(command));
    }
    state
}

/// The socket session `id`'s view answers on.
pub fn control_path(id: &str) -> anyhow::Result<PathBuf> {
    Ok(legume_tui::runtime_dir()?.join(format!("{id}.ctl")))
}

/// A request on its way to the view's loop, and where its answer goes.
struct Asked {
    request: Request,
    reply: mpsc::Sender<Value>,
}

/// How long a request waits for the view's loop to answer it: a figure
/// drawn at a large size takes a while.
const ANSWER: Duration = Duration::from_secs(120);

/// A live view's socket, answered by its loop; removed when dropped.
pub struct Control {
    path: PathBuf,
    asked: mpsc::Receiver<Asked>,
}

impl Control {
    /// The socket of the session this view runs in, when it runs in a
    /// named one; none otherwise, or when it cannot be made (said in the
    /// log: the view runs on without it).
    pub fn start() -> Option<Self> {
        let id = legume_tui::current()?.name()?;
        Self::bind(&id)
            .inspect_err(|e| log::warn!("{id}: no control socket ({e:#}); the view runs on"))
            .ok()
    }

    fn bind(id: &str) -> anyhow::Result<Self> {
        Self::at(control_path(id)?)
    }

    /// Answering on `path`.
    fn at(path: PathBuf) -> anyhow::Result<Self> {
        // A program of this id killed outright left its socket: ids are
        // unique among the live sessions, so it is no one's.
        let _ = std::fs::remove_file(&path);
        let listener = UnixListener::bind(&path)?;
        let (tx, asked) = mpsc::channel();
        std::thread::Builder::new()
            .name("view-control".into())
            .spawn(move || listen(&listener, &tx))?;
        Ok(Self { path, asked })
    }

    /// Answer every request waiting, from the view's loop; whether one
    /// changed the view, which is then drawn again.
    pub fn serve(&self, view: &mut (impl Live + ?Sized)) -> bool {
        let mut changed = false;
        while let Ok(Asked { request, reply }) = self.asked.try_recv() {
            changed |= matches!(request, Request::Set { .. });
            let answer = match answer(view, request) {
                Ok(v) => json!({ "ok": v }),
                Err(e) => json!({ "error": format!("{e:#}") }),
            };
            let _ = reply.send(answer);
        }
        changed
    }
}

impl Drop for Control {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// Arrange `view` as `state` says, all or nothing: when a key is refused,
/// the view is put back as it was before (and drawn again), and the
/// refusal returned.
pub fn set_all(view: &mut (impl Live + ?Sized), state: &Value) -> anyhow::Result<()> {
    let before = view.view_state();
    view.set_view_state(state).inspect_err(|_| {
        // What the view itself wrote sets again; nothing is lost if not.
        let _ = view.set_view_state(&before);
    })
}

/// What the view says to `request`.
fn answer(view: &mut (impl Live + ?Sized), request: Request) -> anyhow::Result<Value> {
    match request {
        Request::State => Ok(snapshot(view)),
        Request::Set { state } => {
            check_view(&state, view.view_name())?;
            set_all(view, &state)?;
            Ok(snapshot(view))
        }
        Request::Export {
            prefix,
            formats,
            look,
        } => {
            anyhow::ensure!(!formats.is_empty(), "no format asked for");
            let kept = state_file(&prefix);
            anyhow::ensure!(
                !Path::new(&kept).exists(),
                "{kept} exists; not replacing it"
            );
            let sheet = sheet_of(view, &look)?;
            let svg = view.saved_figure(sheet)?;
            let written = crate::figure::save_as(&svg, &prefix, &formats)?;
            let state = snapshot(view);
            std::fs::write(&kept, serde_json::to_string_pretty(&state)? + "\n")?;
            Ok(json!({
                "written": written,
                "state_file": kept,
                "sheet": sheet_json(&sheet),
                "state": state,
            }))
        }
    }
}

/// The socket's thread: each connection asked in turn, until the view
/// is gone.
fn listen(listener: &UnixListener, tx: &mpsc::Sender<Asked>) {
    for stream in listener.incoming() {
        let Ok(stream) = stream else { continue };
        if !take(stream, tx) {
            return;
        }
    }
}

/// One connection: its request read, handed on, answered; false once the
/// view's loop is gone.
fn take(stream: UnixStream, tx: &mpsc::Sender<Asked>) -> bool {
    let _ = stream.set_read_timeout(Some(Duration::from_secs(10)));
    let mut line = String::new();
    if BufReader::new(&stream).read_line(&mut line).is_err() {
        return true;
    }
    let (reply, answered) = mpsc::channel();
    let answer = match serde_json::from_str::<Request>(&line) {
        Err(e) => json!({ "error": format!("not a request: {e}") }),
        Ok(request) => {
            if tx.send(Asked { request, reply }).is_err() {
                let _ = say(&stream, &json!({ "error": "the view has ended" }));
                return false;
            }
            answered.recv_timeout(ANSWER).unwrap_or_else(
                |_| json!({ "error": "the view did not answer in time: busy, or ending" }),
            )
        }
    };
    let _ = say(&stream, &answer);
    true
}

fn say(mut stream: &UnixStream, v: &Value) -> std::io::Result<()> {
    writeln!(stream, "{v}")?;
    stream.flush()
}

/// Ask session `id`'s view `request`; what it answered, or why not.
pub fn ask(id: &str, request: &Request) -> anyhow::Result<Value> {
    ask_at(&control_path(id)?, id, request)
}

/// [`ask`] on the socket at `path`, session `id`'s.
fn ask_at(path: &Path, id: &str, request: &Request) -> anyhow::Result<Value> {
    let stream = UnixStream::connect(path).map_err(|e| match e.kind() {
        std::io::ErrorKind::NotFound | std::io::ErrorKind::ConnectionRefused => anyhow::anyhow!(
            "{id} has no view answering: still opening, ended, or a view that keeps no state"
        ),
        _ => anyhow::Error::new(e).context(format!("asking {id}")),
    })?;
    stream.set_read_timeout(Some(ANSWER + Duration::from_secs(5)))?;
    say(&stream, &serde_json::to_value(request)?)?;
    let mut line = String::new();
    BufReader::new(&stream)
        .read_line(&mut line)
        .with_context(|| format!("{id} did not answer"))?;
    let v: Value =
        serde_json::from_str(&line).with_context(|| format!("{id} answered {line:?}"))?;
    match (v.get("ok"), v.get("error").and_then(Value::as_str)) {
        (Some(ok), _) => Ok(ok.clone()),
        (None, Some(e)) => anyhow::bail!("{id}: {e}"),
        _ => anyhow::bail!("{id} answered {line:?}"),
    }
}

#[cfg(test)]
#[path = "../tests/tui/live.rs"]
mod tests;
