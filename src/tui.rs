//! Pieces the terminal views share, most of them legume-tui's drawn in
//! faba's styles ([`LOOK`]): pop-ups, buttons and tabs, lists that keep
//! their selection in view and narrow as you type, work on another thread,
//! the screen loop the views run in; and faba's own key rules (the apply
//! key above all, see [`feed`]).

pub mod browser;
pub mod busy;
pub mod help;

pub use busy::{busy, busy_counting};
pub use help::{keys, Helped, KeyList};
pub use legume_tui::browse::{first_visible, home, moved, tilde, typed_path, window, FindList};
pub use legume_tui::hist::plain_axis;
/// Whether a key is the apply key, Ctrl+R (with Shift or not), and nothing
/// else, so a stray key cannot write anything.
pub use legume_tui::is_run_key as is_apply;
pub use legume_tui::keys::{is_click_or_wheel, is_stray, wheel_key};
pub use legume_tui::run_key_event as apply_key;
pub use legume_tui::save::{next_free, next_free_file, output_problem};
#[cfg(test)]
pub use legume_tui::testing::{missing_keys, screen_text};
pub use legume_tui::view::{end_views, title_keys_width, Screen, ViewsEnd};
pub use legume_tui::widgets::{clock, digits, filled, hidden_marks, Hits};
pub use legume_tui::work::{on_thread, Pending};

use std::time::Duration;

use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use ratatui::layout::Rect;
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::Frame;

use ui::{ACCENTED, DIM, HIGHLIGHT, PLAIN};

/// The styles legume-tui draws faba's views in: the cursor, keys and marks
/// highlighted, folders, notes and borders dim, a complaint accented.
pub const LOOK: legume_tui::Look =
    legume_tui::Look::new(PLAIN, HIGHLIGHT, DIM, HIGHLIGHT, DIM, ACCENTED);

/// What faba's views draw with: data-beans' colours, and legume-tui's
/// chrome (in [`LOOK`]) and histograms.
pub mod ui {
    use super::LOOK;
    use ratatui::text::Line;
    use ratatui::widgets::Block;

    use ratatui::style::Style;

    /// The terminal's own text.
    pub const PLAIN: Style = data_beans::interactive::ui::LOOK.page;
    /// What is secondary: folders, notes, borders.
    pub const DIM: Style = data_beans::interactive::ui::LOOK.hint;
    /// data-beans' accent colour: a complaint, a warning.
    pub const ACCENTED: Style = data_beans::interactive::ui::LOOK.accent;
    /// The accent in bold: keys, the cursor, marks.
    pub const HIGHLIGHT: Style = data_beans::interactive::ui::LOOK.highlight();
    pub use legume_tui::hist::{compact, Binned, Binning, HistPlot, MirrorPlot, MirrorSide, Scale};
    pub use legume_tui::view::Screen;

    /// Title row: a badge naming the view, then `title`, then `extra` dim.
    pub fn header(badge: &str, title: &str, extra: &str) -> Line<'static> {
        legume_tui::view::header(badge, title, extra, &LOOK)
    }

    /// A rounded panel titled `title`, lit when `focused`.
    pub fn panel(title: String, focused: bool) -> Block<'static> {
        legume_tui::view::panel(title, focused, &LOOK)
    }

    /// A help line from `(key, what it does)` pairs.
    pub fn help_line(pairs: &[(&str, &str)]) -> Line<'static> {
        legume_tui::view::help_line(pairs, &LOOK)
    }

    /// A footer while typing: `prompt`, the text so far, then `keys`.
    pub fn input_line(prompt: &str, text: &str, keys: &[(&str, &str)]) -> Line<'static> {
        legume_tui::view::input_line(prompt, text, keys, &LOOK)
    }
}

/// The marker of a list's highlighted row.
pub fn marker(on: bool) -> Span<'static> {
    legume_tui::form::marker(on, &LOOK)
}

/// Whether to write and stop rather than open the view: when `asked`
/// (`--batch-process`), or when there is no terminal to show it in.
pub fn batch(asked: bool, command: &str) -> bool {
    if asked {
        return true;
    }
    let terminal = legume_tui::view::terminal_available();
    if !terminal {
        log::info!("no terminal for the view: {command} runs as with --batch-process");
    }
    !terminal
}

/// Clear a `w` x `h` pop-up centred in `area`, frame it in a panel titled
/// `title`, and return the inside.
pub fn popup_frame(frame: &mut Frame, area: Rect, w: u16, h: u16, title: String) -> Rect {
    legume_tui::widgets::popup_frame(frame, area, w, h, &title, &LOOK)
}

/// `lines` in a pop-up sized to them, centred in `area`; the inside.
pub fn popup(frame: &mut Frame, area: Rect, title: &str, lines: Vec<Line<'static>>) -> Rect {
    legume_tui::widgets::popup(frame, area, title, lines, &LOOK)
}

/// A pop-up with `buttons` on its first line, then `body`; each button's
/// spot goes into `hits`.
pub fn button_popup<T: Clone>(
    frame: &mut Frame,
    area: Rect,
    title: &str,
    buttons: &[(&str, T)],
    body: Vec<Line<'static>>,
    hits: &Hits<T>,
) -> Rect {
    legume_tui::widgets::button_popup(frame, area, title, buttons, body, hits, &LOOK)
}

/// A tab line drawn in `area`: the names, `active` lit, then `button` if
/// any; each spot goes into `hits`.
pub fn tab_bar<T: Clone>(
    area: Rect,
    names: &[&str],
    active: usize,
    button: Option<(&str, T)>,
    hits: &Hits<T>,
    tab: impl Fn(usize) -> T,
) -> Line<'static> {
    legume_tui::widgets::tab_bar(area, names, active, button, hits, tab, &LOOK)
}

/// A footer `line` as lines no wider than `width`, a key never parted
/// from what it does.
pub fn wrap_help(line: Line<'static>, width: u16) -> Vec<Line<'static>> {
    legume_tui::text::wrap_help(line, width, HIGHLIGHT)
}

/// A `width`-cell bar for `done` of `total`: the filled part in `on`.
pub fn bar_spans(done: u64, total: u64, width: usize, on: Style) -> [Span<'static>; 2] {
    legume_tui::widgets::bar_spans(done, total, width, on, &LOOK)
}

/// The spinner's frame `elapsed` into a spin: one every 200 ms.
pub fn spinner_frame(elapsed: Duration) -> char {
    legume_tui::work::spinner_frame_every(elapsed, Duration::from_millis(200))
}

/// The apply key, as footers name it: Ctrl+R, a control byte every
/// terminal sends as it is, so it is never taken for Enter: legume-tui's
/// [`legume_tui::RUN_KEY`], capitalised as faba's footers write keys.
pub const APPLY_KEYS: &str = "Ctrl+R";

/// The key that writes a view's HTML page, as footers name it.
pub const REPORT_KEY: &str = "Ctrl+E";

/// Whether `key` is the page's, [`REPORT_KEY`]: Ctrl+E and no other
/// modifier, but the Shift a terminal may add for a capital E. A view
/// names it in [`View::takes_chord`] where it writes a page.
pub fn is_report(key: &KeyEvent) -> bool {
    key.modifiers.difference(KeyModifiers::SHIFT) == KeyModifiers::CONTROL
        && matches!(key.code, KeyCode::Char('e' | 'E'))
}

/// What a view says when Enter (plain or held with Shift, Alt or Ctrl) did
/// nothing: the key that `does` it. [`legume_tui::enter_hint`]'s words,
/// with the key capitalised as faba's footers write keys.
pub fn enter_hint(does: &str) -> String {
    format!("Enter does nothing here; {APPLY_KEYS} {does}")
}

/// A full-screen view fed by the shared key rules (see [`feed`]); the mouse,
/// notices, which keys the terminal reports and the quit a `faba kill`
/// takes are legume-tui's [`Screen`]'s.
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

    /// Every key the view takes as it is now, for the `?` pop-up
    /// ([`Helped`]): the footer's and the ones it leaves out. None while
    /// it is typing, where `?` is typed; and by default.
    fn help(&self) -> Option<KeyList> {
        None
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

/// Whether `error` is (or wraps) a stop by `faba kill`
/// ([`legume_tui::view::Stopped`]).
pub fn is_stopped(error: &anyhow::Error) -> bool {
    legume_tui::view::is_stopped(error.as_ref())
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

/// Run `view` full screen until it is done, in the command's session
/// (a plain one where there is none), its keys fed by [`feed`], the
/// session keys at the right end of its title row
/// ([`legume_tui::view::run_view`]), `?` showing its keys ([`Helped`]). A
/// `faba kill` takes the view's quit
/// path ([`Screen::quit`]) and then ends it with
/// [`legume_tui::view::Stopped`].
pub fn run_view<V: View>(view: &mut V) -> anyhow::Result<()> {
    let mut helped = Helped::new(view);
    Ok(legume_tui::view::run_view(&TOOL, &LOOK, &mut helped, feed)?)
}

/// Run `f`, the command's own work between or after its views, on the
/// normal screen ([`legume_tui::view::work`]).
pub fn work<T>(f: impl FnOnce() -> T) -> anyhow::Result<T> {
    Ok(legume_tui::view::work(f)?)
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
