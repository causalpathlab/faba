//! Pieces the terminal views share: centred pop-ups, lists that keep their
//! selection in view and narrow as you type, the key rules every view is fed
//! by (the apply key above all), typed paths, and the rules for an output
//! folder.

pub mod browser;
pub mod busy;
pub mod form;

pub use busy::{busy, busy_counting};

/// Frames of the spinners the workspace draws.
pub const SPINNER: &str = "⠁⠂⠄⡀⢀⠠⠐⠈";

use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, TryRecvError};
use std::time::{Duration, Instant};

use data_beans::interactive::ui::{panel, Screen, Tui, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent};
use ratatui::layout::Rect;
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Clear, Paragraph};
use ratatui::Frame;

/// Whether to write and stop rather than open the view: when `asked`
/// (`--batch-process`), or when there is no terminal to show it in.
pub fn batch(asked: bool, command: &str) -> bool {
    if asked {
        return true;
    }
    let terminal = data_beans::interactive::tui_available();
    if !terminal {
        log::info!("no terminal for the view: {command} runs as with --batch-process");
    }
    !terminal
}

/// A `w` x `h` rectangle centred in `area`, shrunk to fit it.
pub fn centered(area: Rect, w: u16, h: u16) -> Rect {
    let (w, h) = (w.min(area.width), h.min(area.height));
    Rect::new(
        area.x + (area.width - w) / 2,
        area.y + (area.height - h) / 2,
        w,
        h,
    )
}

/// Clear a `w` x `h` pop-up centred in `area`, frame it in a panel titled
/// `title`, and return the inside.
pub fn popup_frame(frame: &mut Frame, area: Rect, w: u16, h: u16, title: String) -> Rect {
    let rect = centered(area, w, h);
    frame.render_widget(Clear, rect);
    let block = panel(title, true);
    let inner = block.inner(rect);
    frame.render_widget(block, rect);
    inner
}

/// A row of buttons to click, `[ label ]` each, two apart: the line, and
/// each button's column (from the line's start) and width.
pub fn buttons(labels: &[&str]) -> (Line<'static>, Vec<(u16, u16)>) {
    let style = HIGHLIGHT.add_modifier(ratatui::style::Modifier::REVERSED);
    let mut spans = vec![Span::raw(" ")];
    let mut spots = Vec::new();
    let mut x = 1;
    for (k, label) in labels.iter().enumerate() {
        if k > 0 {
            spans.push(Span::raw("  "));
            x += 2;
        }
        let text = format!("[ {label} ]");
        let w = text.chars().count() as u16;
        spans.push(Span::styled(text, style));
        spots.push((x, w));
        x += w;
    }
    (Line::from(spans), spots)
}

/// `lines` in a pop-up sized to them, at least 40 wide, centred in `area`;
/// the inside, where they went.
pub fn popup(frame: &mut Frame, area: Rect, title: &str, lines: Vec<Line<'static>>) -> Rect {
    let width = lines.iter().map(Line::width).max().unwrap_or(0) as u16 + 4;
    let inner = popup_frame(
        frame,
        area,
        width.max(40),
        lines.len() as u16 + 2,
        title.into(),
    );
    frame.render_widget(Paragraph::new(lines), inner);
    inner
}

/// A pop-up with `buttons` on its first line (where a long body cannot push
/// them out of sight), then `body`; each button's spot goes into `hits`.
pub fn button_popup<T: Clone>(
    frame: &mut Frame,
    area: Rect,
    title: &str,
    buttons_: &[(&str, T)],
    body: Vec<Line<'static>>,
    hits: &Hits<T>,
) -> Rect {
    let labels: Vec<&str> = buttons_.iter().map(|(l, _)| *l).collect();
    let (line, spots) = buttons(&labels);
    let mut lines = vec![line, Line::raw("")];
    lines.extend(body);
    let inner = popup(frame, area, title, lines);
    for (&(dx, w), (_, what)) in spots.iter().zip(buttons_) {
        hits.add(Rect::new(inner.x + dx, inner.y, w, 1), what.clone());
    }
    inner
}

/// A tab line drawn in `area`: the names, `active` lit, then `button` if
/// any; each name's and the button's spot goes into `hits`.
pub fn tab_bar<T: Clone>(
    area: Rect,
    names: &[&str],
    active: usize,
    button: Option<(&str, T)>,
    hits: &Hits<T>,
    tab: impl Fn(usize) -> T,
) -> Line<'static> {
    let mut spans = vec![Span::raw(" ")];
    let mut x = area.x + 1;
    for (i, name) in names.iter().enumerate() {
        let style = if i == active {
            HIGHLIGHT
        } else {
            data_beans::interactive::ui::DIM
        };
        let w = name.chars().count() as u16;
        hits.add(Rect::new(x, area.y, w, 1), tab(i));
        spans.push(Span::styled(format!("{name}  "), style));
        x += w + 2;
    }
    if let Some((label, what)) = button {
        let (line, spots) = buttons(&[label]);
        if let Some(&(dx, w)) = spots.first() {
            hits.add(Rect::new(x + dx, area.y, w, 1), what);
        }
        spans.extend(line.spans);
    }
    Line::from(spans)
}

/// A footer `line` as lines no wider than `width`, so no key is cut off:
/// a help line ([`help_line`], [`input_line`]) is a lead span, then a key
/// and what it does, span by span, and a key is never parted from what it
/// does. A line that fits is left as it is.
///
/// [`help_line`]: data_beans::interactive::ui::help_line
/// [`input_line`]: data_beans::interactive::ui::input_line
pub fn wrap_help(line: Line<'static>, width: u16) -> Vec<Line<'static>> {
    let width = usize::from(width);
    if line.width() <= width {
        return vec![line];
    }
    // A key (drawn highlighted) starts a group with what follows it.
    let mut groups: Vec<Vec<Span<'static>>> = vec![Vec::new()];
    for span in line.spans {
        if span.style == HIGHLIGHT {
            groups.push(Vec::new());
        }
        groups.last_mut().expect("one group at least").push(span);
    }
    let mut out = Vec::new();
    let mut row: Vec<Span<'static>> = Vec::new();
    let mut used = 0;
    for group in groups {
        let w: usize = group.iter().map(Span::width).sum();
        if used + w > width && used > 1 {
            out.push(Line::from(std::mem::take(&mut row)));
            row.push(Span::raw(" "));
            used = 1;
        }
        used += w;
        row.extend(group);
    }
    out.push(Line::from(row));
    out
}

/// The keys and descriptions of `help` that `screen` does not show, as a
/// test reads a drawn terminal: empty when every one is there.
#[cfg(test)]
pub fn missing_keys(help: &Line, screen: &str) -> Vec<String> {
    help.spans
        .iter()
        .map(|s| s.content.trim().to_string())
        .filter(|t| !t.is_empty() && !screen.contains(t.as_str()))
        .collect()
}

/// A `width`-cell bar for `done` of `total`: the filled part in `on`, the
/// rest dim.
pub fn bar_spans(done: u64, total: u64, width: usize, on: Style) -> [Span<'static>; 2] {
    let full = filled(done, total, width);
    [
        Span::styled("█".repeat(full), on),
        Span::styled("░".repeat(width.saturating_sub(full)), DIM),
    ]
}

/// `d` as `hh:mm:ss`.
pub fn clock(d: std::time::Duration) -> String {
    let s = d.as_secs();
    format!("{:02}:{:02}:{:02}", s / 3600, s / 60 % 60, s % 60)
}

/// How many of `width` cells a bar shows filled for `done` of `total`.
pub fn filled(done: u64, total: u64, width: usize) -> usize {
    match total {
        0 => 0,
        _ => (done.min(total) as f64 / total as f64 * width as f64).round() as usize,
    }
}

/// The [`SPINNER`] frame `elapsed` into a spin: one every 200 ms.
pub fn spinner_frame(elapsed: Duration) -> char {
    let frames = SPINNER.chars().count();
    SPINNER
        .chars()
        .nth((elapsed.as_millis() / 200) as usize % frames)
        .expect("a frame within the spinner")
}

/// The result of work on another thread, as a view picks it up between
/// frames from the channel it comes on.
pub struct Pending<'a, T> {
    rx: &'a Receiver<T>,
    /// The result, and when it came.
    out: Option<(T, Instant)>,
    /// The work ended without a result: it panicked.
    gone: bool,
}

impl<'a, T> Pending<'a, T> {
    /// Waiting on `rx`.
    pub fn new(rx: &'a Receiver<T>) -> Self {
        Self {
            rx,
            out: None,
            gone: false,
        }
    }

    /// Pick up the result, if it has come.
    pub fn poll(&mut self) {
        match self.rx.try_recv() {
            Ok(out) => self.out = Some((out, Instant::now())),
            Err(TryRecvError::Disconnected) => self.gone = true,
            Err(TryRecvError::Empty) => {}
        }
    }

    /// The result, once it came.
    pub fn out(&self) -> Option<&T> {
        self.out.as_ref().map(|(out, _)| out)
    }

    /// The time since the result came; none before.
    pub fn since(&self) -> Option<Duration> {
        self.out.as_ref().map(|(_, at)| at.elapsed())
    }

    /// Whether the work ended without a result: it panicked.
    pub fn gone(&self) -> bool {
        self.out.is_none() && self.gone
    }

    /// The result itself, once it came.
    pub fn take(self) -> Option<T> {
        self.out.map(|(out, _)| out)
    }
}

/// The screen of `view` drawn `w` x `h`, as text.
#[cfg(test)]
pub fn screen_text(view: &mut impl Screen, w: u16, h: u16) -> String {
    let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
    t.draw(|f| view.render(f)).unwrap();
    let buf = t.backend().buffer();
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

/// Run `work` on a thread while `meanwhile` waits here on the channel its
/// result comes on; what `meanwhile` returns, and the work's result: the
/// one `meanwhile` took, else the one left on the channel, else (the work
/// panicked) an error saying `what` panicked, with the panic's message.
pub fn on_thread<T: Send, R>(
    what: &str,
    work: impl FnOnce() -> anyhow::Result<T> + Send,
    meanwhile: impl FnOnce(&Receiver<anyhow::Result<T>>) -> (R, Option<anyhow::Result<T>>),
) -> (R, anyhow::Result<T>) {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::scope(|s| {
        let worker = s.spawn(move || {
            // The receiver outlives the scope, so the send cannot fail.
            let _ = tx.send(work());
        });
        let (waited, taken) = meanwhile(&rx);
        // Joined here, a panic is the work's error: the scope raises it
        // again only for threads it joins itself.
        let joined = worker.join();
        let out = match (taken.or_else(|| rx.try_recv().ok()), joined) {
            (Some(out), _) => out,
            (None, Err(panic)) => Err(anyhow::anyhow!("{what} panicked: {}", panic_text(&*panic))),
            (None, Ok(())) => Err(anyhow::anyhow!("{what} ended without a result")),
        };
        (waited, out)
    })
}

/// A panic's message, when it is text.
fn panic_text(panic: &(dyn std::any::Any + Send)) -> &str {
    panic
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| panic.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("no message")
}

/// The x axis of a text plot in `area` without its tick marks: the axis
/// line runs through where the plot put a `┴` under each label. For labels
/// that name a span (a metagene's regions), where a tick would read as a
/// position.
pub fn plain_axis(buf: &mut ratatui::buffer::Buffer, area: Rect) {
    for y in area.top()..area.bottom() {
        for x in area.left()..area.right() {
            if buf[(x, y)].symbol() == "┴" {
                buf[(x, y)].set_symbol("─");
            }
        }
    }
}

/// The marker of a list's highlighted row.
pub fn marker(on: bool) -> Span<'static> {
    Span::styled(if on { "▸ " } else { "  " }, HIGHLIGHT)
}

/// The `rows` lines of `lines` around line `at`.
pub fn window(lines: Vec<Line<'static>>, at: usize, rows: usize) -> Vec<Line<'static>> {
    let first = first_visible(at, lines.len(), rows);
    lines.into_iter().skip(first).take(rows).collect()
}

/// The first of `rows` visible items of a `len`-long list that keeps item
/// `at` near the middle.
pub fn first_visible(at: usize, len: usize, rows: usize) -> usize {
    at.saturating_sub(rows.saturating_sub(1) / 2)
        .min(len.saturating_sub(rows))
}

/// `at` moved by `d`, kept within `0..=last`.
pub fn moved(at: usize, d: isize, last: usize) -> usize {
    at.saturating_add_signed(d).min(last)
}

/// How many characters `n` takes written out, its sign included.
pub fn digits(n: i64) -> usize {
    let width = n
        .unsigned_abs()
        .checked_ilog10()
        .map_or(1, |d| d as usize + 1);
    width + usize::from(n < 0)
}

/// How many of a scrolled run are hidden before and after, as `‹2 5›`;
/// empty when none are.
pub fn hidden_marks(left: usize, right: usize) -> String {
    let left = (left > 0).then(|| format!("‹{left}"));
    let right = (right > 0).then(|| format!("{right}›"));
    [left, right]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .join(" ")
}

/// The home folder, if known.
pub fn home() -> Option<PathBuf> {
    std::env::var_os("HOME").map(PathBuf::from)
}

/// `p` with the home folder as `~`.
pub fn tilde(p: &Path) -> String {
    if let Some(rest) = home().and_then(|h| p.strip_prefix(h).ok().map(Path::to_path_buf)) {
        return if rest.as_os_str().is_empty() {
            "~".into()
        } else {
            format!("~/{}", rest.display())
        };
    }
    p.display().to_string()
}

/// A path as typed into a line: `~` is home, and the rest is made absolute
/// with `.` and `..` resolved (see [`browser::normalize`]).
pub fn typed_path(text: &str) -> PathBuf {
    let path = match (text.strip_prefix('~'), home()) {
        (Some(rest), Some(h)) if rest.is_empty() || rest.starts_with('/') => {
            h.join(rest.trim_start_matches('/'))
        }
        _ => PathBuf::from(text),
    };
    browser::normalize(&path)
}

/// `base/stem`, or `base/stem2`, `base/stem3`, ... : the first that does
/// not exist yet.
pub(crate) fn next_free(base: &Path, stem: &str) -> PathBuf {
    next_free_file(base, stem, "")
}

/// `base/{stem}{ext}`, or `base/{stem}2{ext}`, `base/{stem}3{ext}`, ... :
/// the first that does not exist yet, numbered before the extension.
pub(crate) fn next_free_file(base: &Path, stem: &str, ext: &str) -> PathBuf {
    let mut out = base.join(format!("{stem}{ext}"));
    let mut n = 2;
    while out.exists() {
        out = base.join(format!("{stem}{n}{ext}"));
        n += 1;
    }
    out
}

/// Why `out` cannot take a run's outputs, if it cannot: it is unnamed, a
/// file, or a folder with files in it.
pub(crate) fn output_problem(out: &str) -> Option<String> {
    let path = Path::new(out);
    if out.is_empty() {
        return Some("name a directory".into());
    }
    if path.is_file() {
        return Some(format!("{out} is a file"));
    }
    let non_empty = std::fs::read_dir(path).is_ok_and(|mut d| d.next().is_some());
    non_empty.then(|| format!("{out} already contains files; choose an empty one"))
}

/// The apply key, as footers name it: Ctrl+R, a control byte every
/// terminal sends as it is, so it is never taken for Enter.
pub const APPLY_KEYS: &str = "Ctrl+R";

/// What a view says when Enter (plain or held with Shift, Alt or Ctrl) did
/// nothing: the key that `does` it.
pub fn enter_hint(does: &str) -> String {
    format!("Enter does nothing here; {APPLY_KEYS} {does}")
}

/// A full-screen view fed by the shared key rules (see [`feed`]); the mouse,
/// notices and which keys the terminal reports are data_beans' `Screen`'s.
pub trait View: Screen {
    /// Enter held with Shift, Alt or Ctrl came, which must not act as
    /// Enter: say what the apply key does here.
    fn stray_enter(&mut self) {}

    /// Whether the view has something for the apply key to do; elsewhere it
    /// is dropped, so it never acts as a letter key.
    fn takes_apply(&self) -> bool {
        false
    }

    /// Whether the view takes `key`, a chord (Ctrl, Alt) that is dropped
    /// otherwise; by default none, so no chord ever acts as its letter.
    fn takes_chord(&self, _key: &KeyEvent) -> bool {
        false
    }

    /// `faba kill` asked the view to stop: leave as its quit key does,
    /// without asking, unsaved edits dropped; work that must finish (a
    /// write under way) is let finish, [`run_view`] waiting until the
    /// view is done. By default as Ctrl-C ([`Screen::interrupt`]).
    fn quit(&mut self) {
        self.interrupt();
    }
}

/// The ↑/↓ a turn of the wheel stands for.
pub fn wheel_key(m: &MouseEvent) -> Option<KeyEvent> {
    use ratatui::crossterm::event::MouseEventKind;
    let code = match m.kind {
        MouseEventKind::ScrollUp => KeyCode::Up,
        MouseEventKind::ScrollDown => KeyCode::Down,
        _ => return None,
    };
    Some(KeyEvent::new(code, KeyModifiers::NONE))
}

/// Where a view drew what, for the mouse: rectangles, each with what a
/// click there means. Kept behind a `RefCell` so drawing, which only
/// reads the view, can write it.
pub struct Hits<T> {
    spots: std::cell::RefCell<Vec<(Rect, T)>>,
}

impl<T> Default for Hits<T> {
    fn default() -> Self {
        Self {
            spots: std::cell::RefCell::default(),
        }
    }
}

impl<T: Clone> Hits<T> {
    /// Forget the last draw's; call at the top of a draw.
    pub fn clear(&self) {
        self.spots.borrow_mut().clear();
    }

    pub fn add(&self, area: Rect, what: T) {
        self.spots.borrow_mut().push((area, what));
    }

    /// A list drawn in `area` from item `first`: `what(i)` for item `i` on
    /// each of its rows, up to `len` items.
    pub fn rows(&self, area: Rect, first: usize, len: usize, what: impl Fn(usize) -> T) {
        for (k, i) in (first..len).take(area.height as usize).enumerate() {
            let row = Rect::new(area.x, area.y + k as u16, area.width, 1);
            self.add(row, what(i));
        }
    }

    /// A cell where `what` was drawn, to click in tests.
    #[cfg(test)]
    pub fn spot(&self, what: &T) -> Option<(u16, u16)>
    where
        T: PartialEq,
    {
        let spots = self.spots.borrow();
        let (r, _) = spots.iter().rev().find(|(_, t)| t == what)?;
        Some((r.x, r.y))
    }

    /// What is at a cell: the last drawn there, which is on top.
    pub fn at(&self, column: u16, row: u16) -> Option<T> {
        let pos = ratatui::layout::Position::new(column, row);
        self.spots
            .borrow()
            .iter()
            .rev()
            .find(|(r, _)| r.contains(pos))
            .map(|(_, t)| t.clone())
    }
}

/// Hand `key` to `view` by the rules every view shares: Enter held with
/// Shift, Alt or Ctrl goes to [`View::stray_enter`], other Ctrl and Alt
/// chords are dropped unless the view names them ([`View::takes_chord`]),
/// so Ctrl+S is never `s`, and so is the apply key in a view with nothing for
/// it to do. Views and the widgets in them can take what comes as meant.
pub fn feed(view: &mut impl View, key: KeyEvent) {
    if is_stray(&key) && !view.takes_chord(&key) {
        if key.code == KeyCode::Enter {
            view.stray_enter();
        }
    } else if !is_apply(&key) || view.takes_apply() {
        view.handle_key(key);
    }
}

/// The tool faba's views belong to: their sessions are named
/// `faba-pileup-1`, `faba-qc-2`, and `faba attach` finds them.
pub const TOOL: legume_tui::Tool = legume_tui::Tool::new("faba");

/// A command's views were stopped by `faba kill`: its quit path was taken
/// (work under way finished, edits dropped) and nothing after the view
/// runs. `main` ends quietly on it.
#[derive(Debug)]
pub struct Stopped;

impl std::fmt::Display for Stopped {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("stopped by `faba kill`")
    }
}

impl std::error::Error for Stopped {}

/// Whether `error` is (or wraps) [`Stopped`].
pub fn is_stopped(error: &anyhow::Error) -> bool {
    error.chain().any(|e| e.is::<Stopped>())
}

/// The paragraph a view's long help ends with about its session: how it
/// runs, what goes on while detached (`what_goes_on`, with its verb), the
/// keys, how to return (started again with nothing named, the live
/// sessions are offered), and `more` at the end.
pub fn session_help(what_goes_on: &str, more: &str) -> String {
    let help = legume_tui::cli::SessionHelp::new()
        .what_goes_on(what_goes_on)
        .more(more);
    legume_tui::cli::session_help(&TOOL, help)
}

/// How the views' session starts a child that runs on its own (`faba
/// run`'s pipeline, a download): the user's opt-out passed on; none
/// outside a session (in tests).
pub fn child_env() -> legume_tui::ChildEnv {
    legume_tui::current()
        .map(legume_tui::Session::child_env)
        .unwrap_or_default()
}

thread_local! {
    /// The screen the command's views share while its session is up (see
    /// [`run_view`]), from the first view until [`end_views`].
    static VIEWS: std::cell::RefCell<Option<Tui<'static>>> = const { std::cell::RefCell::new(None) };

    /// Columns the session keys take at the right end of the title row of
    /// the view being drawn ([`Framed`]); 0 outside such a draw.
    static TITLE_KEYS: std::cell::Cell<u16> = const { std::cell::Cell::new(0) };
}

/// Run `view` full screen until it is done, in the command's session
/// (`legume_tui::current`, opened as the command started; a plain one
/// where there is none), its keys fed by [`feed`]. In a session the
/// screen stays up from the first view to [`end_views`], so the command's
/// views follow one another on it. The title row ends with the session
/// keys. A `faba kill` takes the view's quit path ([`View::quit`]): work
/// it must finish (a write under way) is waited for, off the screen, and
/// then the view ends with [`Stopped`].
pub fn run_view<V: View>(view: &mut V) -> anyhow::Result<()> {
    let Some(session) = legume_tui::current() else {
        let plain = legume_tui::Session::plain(&TOOL);
        return run_framed(&mut Tui::new(&plain), view);
    };
    // Taken out while the view runs, so a view started from inside another
    // gets a screen of its own rather than a second borrow.
    let mut tui = VIEWS
        .with_borrow_mut(Option::take)
        .unwrap_or_else(|| Tui::new(session));
    let shown = run_framed(&mut tui, view);
    VIEWS.with_borrow_mut(|views| {
        views.get_or_insert(tui);
    });
    shown
}

/// [`run_view`] on `tui`.
fn run_framed<V: View>(tui: &mut Tui<'_>, view: &mut V) -> anyhow::Result<()> {
    let session = tui.session();
    let mut framed = Framed {
        view,
        session,
        keys: data_beans::interactive::session_keys(),
    };
    // A stop interrupts the view, which takes its quit path (see
    // `Framed::interrupt`), before the loop ends.
    tui.run_with(&mut framed, feed)?;
    if !session.stop_requested() {
        return Ok(());
    }
    let view = framed.view;
    tui.work(|| {
        while !view.done() {
            view.tick();
            std::thread::sleep(data_beans::interactive::ui::TICK);
        }
    })?;
    Err(Stopped.into())
}

/// Bring the views' screen down, with the log records held while it was
/// up: before the command prints or writes on the normal screen, and as it
/// ends. A later view brings it up again.
pub fn end_views() {
    drop(VIEWS.with_borrow_mut(Option::take));
}

/// [`end_views`] as it is dropped: as a command ends, however it ends.
pub struct ViewsEnd;

impl Drop for ViewsEnd {
    fn drop(&mut self) {
        end_views();
    }
}

/// Run `f`, the command's own work between or after its views (what it
/// prints, say), on the normal screen, the views' screen coming back for
/// the next view; where no view has run, as it is. A `faba kill` meanwhile
/// lets `f` finish.
pub fn work<T>(f: impl FnOnce() -> T) -> anyhow::Result<T> {
    let mut tui = VIEWS.with_borrow_mut(Option::take);
    let out = match tui.as_mut() {
        Some(t) => t.work(f),
        None => Ok(f()),
    };
    if let Some(t) = tui {
        VIEWS.with_borrow_mut(|views| {
            views.get_or_insert(t);
        });
    }
    out
}

/// Columns the session keys take at the right end of the title row of the
/// view being drawn, so a mark there can stand clear of them; 0 when none.
pub fn title_keys_width() -> u16 {
    TITLE_KEYS.get()
}

/// A view as [`run_view`] shows it: its own screen with the session keys
/// (`keys`) at the right end of its title row, and a stop taken as its
/// quit path.
pub(crate) struct Framed<'a, V> {
    pub(crate) view: &'a mut V,
    pub(crate) session: &'a legume_tui::Session,
    pub(crate) keys: Vec<(String, String)>,
}

/// The session keys as drawn on a title row: the line, and the columns it
/// takes there, one more for a space before it.
fn keys_line(keys: &[(String, String)]) -> Option<Line<'static>> {
    if keys.is_empty() {
        return None;
    }
    let pairs: Vec<_> = keys.iter().map(|(k, w)| (k.as_str(), w.as_str())).collect();
    Some(data_beans::interactive::ui::help_line(&pairs))
}

/// The session keys drawn over the right end of `area`'s top row, the
/// title there giving way, so the detach key is always on screen.
pub fn draw_session_keys(frame: &mut Frame, area: Rect, keys: &[(String, String)]) {
    let Some(line) = keys_line(keys).filter(|_| area.height > 0) else {
        return;
    };
    let width = (line.width() as u16).min(area.width);
    let at = Rect::new(area.right() - width, area.y, width, 1);
    frame.render_widget(Clear, at);
    frame.render_widget(Paragraph::new(line), at);
}

impl<V: View> View for Framed<'_, V> {
    fn stray_enter(&mut self) {
        self.view.stray_enter();
    }

    fn takes_apply(&self) -> bool {
        self.view.takes_apply()
    }

    fn takes_chord(&self, key: &KeyEvent) -> bool {
        self.view.takes_chord(key)
    }

    fn quit(&mut self) {
        self.view.quit();
    }
}

impl<V: View> Screen for Framed<'_, V> {
    fn render(&mut self, frame: &mut Frame) {
        let keys = keys_line(&self.keys).map_or(0, |l| l.width() as u16 + 1);
        TITLE_KEYS.set(keys);
        self.view.render(frame);
        TITLE_KEYS.set(0);
        draw_session_keys(frame, frame.area(), &self.keys);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.view.handle_key(key);
    }

    /// Ctrl-C as the view takes it; a `kill` as its quit path.
    fn interrupt(&mut self) {
        if self.session.stop_requested() {
            self.view.quit();
        } else {
            self.view.interrupt();
        }
    }

    fn done(&self) -> bool {
        self.view.done()
    }

    fn pending_work(&self) -> Option<String> {
        self.view.pending_work()
    }

    fn do_work(&mut self) {
        self.view.do_work();
    }

    fn tick(&mut self) -> bool {
        self.view.tick()
    }

    fn takes_mouse(&self) -> bool {
        self.view.takes_mouse()
    }

    fn mouse(&mut self, event: MouseEvent) -> bool {
        self.view.mouse(event)
    }

    fn reports_chords(&self) -> bool {
        self.view.reports_chords()
    }

    fn take_notice(&mut self) -> Option<String> {
        self.view.take_notice()
    }

    fn status(&self) -> Option<String> {
        self.view.status()
    }
}

/// Whether a mouse event is one the views act on: a left click or a turn of
/// the wheel.
pub fn is_click_or_wheel(m: &MouseEvent) -> bool {
    use ratatui::crossterm::event::{MouseButton, MouseEventKind};
    matches!(
        m.kind,
        MouseEventKind::Down(MouseButton::Left)
            | MouseEventKind::ScrollUp
            | MouseEventKind::ScrollDown
    )
}

/// Whether `key` is the apply key: Ctrl+R (with Shift or not, which
/// terminals do not tell apart), and nothing else, so a stray key cannot
/// write anything.
pub fn is_apply(key: &KeyEvent) -> bool {
    key.modifiers.contains(KeyModifiers::CONTROL) && matches!(key.code, KeyCode::Char('r' | 'R'))
}

/// The apply key, as a button sends it.
pub fn apply_key() -> KeyEvent {
    KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL)
}

/// Whether `key`, as the terminal sent it, is held with a modifier that
/// makes it another key (Ctrl, Alt, or Shift on Enter) and is not the apply
/// key. Views never see such keys, so Shift+Enter, Alt+Enter or Ctrl+Enter
/// never acts as a plain Enter, nor Ctrl+S as `s`.
fn is_stray(key: &KeyEvent) -> bool {
    let chord = KeyModifiers::CONTROL
        | KeyModifiers::ALT
        | KeyModifiers::SUPER
        | KeyModifiers::META
        | KeyModifiers::HYPER;
    let held = key.modifiers.intersects(chord)
        || (key.code == KeyCode::Enter && key.modifiers.contains(KeyModifiers::SHIFT));
    held && !is_apply(key)
}

/// What has been typed into a file list to narrow it: letters, digits and
/// punctuation type, so a file pane takes no letter keys of its own.
#[derive(Default, Debug, Clone)]
pub struct Find {
    pub text: String,
}

impl Find {
    /// Take `key` if it edits the text: a typed character (not Space, which
    /// panes keep for choosing), or Backspace and Esc while there is text.
    pub fn key(&mut self, key: &KeyEvent) -> bool {
        match key.code {
            KeyCode::Char(c) if c != ' ' && !is_apply(key) => self.text.push(c),
            KeyCode::Backspace if !self.text.is_empty() => {
                self.text.pop();
            }
            KeyCode::Esc if !self.text.is_empty() => self.text.clear(),
            _ => return false,
        }
        true
    }

    /// Where `key` jumps when typing a path: `/` with nothing typed to the
    /// root, `~` to home, and `/` after a name to that folder in `cwd` (the
    /// name as typed, else `dir`, the highlighted folder), so typing
    /// `/data/runs/` walks there.
    pub fn jump(&self, key: &KeyEvent, cwd: &Path, dir: Option<PathBuf>) -> Option<PathBuf> {
        let KeyCode::Char(c @ ('/' | '~')) = key.code else {
            return None;
        };
        match (c, self.text.as_str()) {
            ('/', "") => cwd.ancestors().last().map(Path::to_path_buf),
            ('~', "") => home(),
            ('/', "..") => cwd.parent().map(Path::to_path_buf),
            ('/', name) => Some(cwd.join(name)).filter(|p| p.is_dir()).or(dir),
            _ => None,
        }
    }

    /// A pane title's tail naming the text, if any.
    pub fn tag(&self) -> String {
        if self.text.is_empty() {
            String::new()
        } else {
            format!("· find: {} ", self.text)
        }
    }
}

/// What a narrowed list says when nothing matches the typed text.
pub const NOTHING_MATCHES: &str = "nothing matches; Backspace or Esc";

/// A list of named items narrowed by what has been typed, with a cursor
/// over those shown: a folder's listing, a site's releases or species.
pub struct FindList<T> {
    all: Vec<T>,
    /// `all`'s names in lower case, to match against.
    lower: Vec<String>,
    /// Indices into `all` of the items shown.
    shown: Vec<usize>,
    /// The cursor, among the items shown.
    pub at: usize,
    pub find: Find,
}

impl<T> Default for FindList<T> {
    fn default() -> Self {
        Self {
            all: Vec::new(),
            lower: Vec::new(),
            shown: Vec::new(),
            at: 0,
            find: Find::default(),
        }
    }
}

impl<T: AsRef<str>> FindList<T> {
    /// Show all of `all`, with nothing typed and the cursor at `at`.
    pub fn set(&mut self, all: Vec<T>, at: usize) {
        self.lower = all.iter().map(|x| x.as_ref().to_lowercase()).collect();
        self.shown = (0..all.len()).collect();
        self.all = all;
        self.at = at.min(self.shown.len().saturating_sub(1));
        self.find = Find::default();
    }

    /// The items shown, in order.
    pub fn shown(&self) -> impl Iterator<Item = &T> {
        self.shown.iter().map(|&i| &self.all[i])
    }

    pub fn len(&self) -> usize {
        self.shown.len()
    }

    pub fn is_empty(&self) -> bool {
        self.shown.is_empty()
    }

    /// The item under the cursor.
    pub fn highlighted(&self) -> Option<&T> {
        self.shown.get(self.at).map(|&i| &self.all[i])
    }

    pub fn step(&mut self, d: isize) {
        self.at = moved(self.at, d, self.len().saturating_sub(1));
    }

    /// Put the cursor on the item shown named `name`; whether there is one.
    pub fn select(&mut self, name: &str) -> bool {
        let at = self.shown().position(|x| x.as_ref() == name);
        self.at = at.unwrap_or(self.at);
        at.is_some()
    }

    /// Take `key` if it edits the typed text (see [`Find::key`]) or moves
    /// the cursor (↑/↓, PgUp/PgDn; letters type, so not `j`/`k`).
    pub fn key(&mut self, key: &KeyEvent) -> bool {
        if self.find.key(key) {
            self.narrow();
            return true;
        }
        let d = match key.code {
            KeyCode::Up => -1,
            KeyCode::Down => 1,
            KeyCode::PageUp => -10,
            KeyCode::PageDown => 10,
            _ => return false,
        };
        self.step(d);
        true
    }

    /// Show the items holding the typed text, any case, never `..`; the
    /// cursor on the first that starts with it, else the first.
    fn narrow(&mut self) {
        let text = self.find.text.to_lowercase();
        let lower = &self.lower;
        self.shown = (0..lower.len())
            .filter(|&i| text.is_empty() || (lower[i] != ".." && lower[i].contains(&text)))
            .collect();
        self.at = self
            .shown
            .iter()
            .position(|&i| lower[i].starts_with(&text))
            .unwrap_or(0);
    }

    /// What an empty list says: `empty`, or that nothing matches.
    pub fn empty_note<'a>(&self, empty: &'a str) -> &'a str {
        if self.find.text.is_empty() {
            empty
        } else {
            NOTHING_MATCHES
        }
    }
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
