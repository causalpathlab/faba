//! Drawing the setup view: the header and tabs, each screen's body, the
//! preview and file pop-ups, and the footer of keys.

use std::path::Path;
use std::time::Duration;

use data_beans::interactive::ui::{
    header, help_line, input_line, panel, ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Wrap};
use ratatui::Frame;

use super::child::{ended_word, Failed};
use super::fetch::{Catalogue, Fetch, Source, SOURCES};
use super::form::{self, Kind};
use super::inputs::{InputsFocus, OutState, Role, Row};
use super::steps::Step;
use super::{script, App, Hit, Page, Target};
use crate::tui::{
    button_popup, filled, first_visible, marker, popup, popup_frame, tab_bar, tilde, window,
    APPLY_KEYS,
};

/// Width of the inputs panel on the Inputs screen.
const INPUTS_WIDTH: u16 = 48;
/// Width of the progress bar on the Run screen.
const BAR: usize = 40;

pub(super) fn draw(app: &App, frame: &mut Frame) {
    app.hits.clear();
    let area = frame.area();
    let [top, tabs, body, footer] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(1),
        Constraint::Fill(1),
        Constraint::Length(1),
    ])
    .areas(area);

    let (title, extra) = match (&app.job, app.page) {
        (Some(job), Page::Run) => (format!("faba run · {}", job_state(app)), tilde(&job.out)),
        _ => {
            let (fg, bg) = (app.inputs.fg().len(), app.inputs.bg().len());
            (
                format!("faba run · {} BAMs: {fg} fg, {bg} bg", fg + bg),
                tilde(&app.inputs.bams.cwd),
            )
        }
    };
    frame.render_widget(header("run", &title, &extra), top);

    let pages = [
        (Page::Inputs, "Inputs"),
        (Page::Steps, "Steps"),
        (Page::Flags, "Flags"),
        (Page::Run, "Run"),
    ];
    // A button for what the apply key does, or to stop a run.
    let button = if app.running() {
        ("■ stop", Hit::Stop)
    } else {
        ("▶ preview", Hit::Preview)
    };
    let names: Vec<&str> = pages.iter().map(|(_, n)| *n).collect();
    let active = pages.iter().position(|(p, _)| *p == app.page).unwrap_or(0);
    let tabs_line = tab_bar(tabs, &names, active, Some(button), &app.hits, |i| {
        Hit::Tab(pages[i].0)
    });
    app.hits.add(body, Hit::Body);
    let reserve = (title.chars().count() + extra.chars().count() + 12).max(tabs_line.width());
    frame.render_widget(tabs_line, tabs);
    if area.height >= 2 {
        let corner = Rect::new(area.x, area.y, area.width, 2);
        let reserve = u16::try_from(reserve).unwrap_or(u16::MAX);
        crate::figure::logo::draw_mini_logo(frame.buffer_mut(), corner, reserve);
    }

    match app.page {
        Page::Inputs => draw_inputs(app, frame, body),
        Page::Steps => draw_steps(app, frame, body),
        Page::Flags => draw_flags(app, frame, body),
        Page::Run => draw_run(app, frame, body),
    }
    if app.picking.is_some() {
        draw_picking(app, frame, body);
    }
    if let Some(c) = &app.catalogue {
        draw_catalogue(c, frame, body);
    }
    if let Some(f) = app.fetch.as_ref().filter(|_| app.fetch_shown) {
        draw_fetch(f, frame, body);
    }
    if app.preview {
        draw_preview(app, frame, body);
    }
    if app.ended_shown {
        draw_ended(app, frame, body);
    }
    frame.render_widget(footer_line(app), footer);
}

/// `hh:mm:ss`.
fn clock(d: Duration) -> String {
    let s = d.as_secs();
    format!("{:02}:{:02}:{:02}", s / 3600, s / 60 % 60, s % 60)
}

/// How the run stands, in a word or few.
fn job_state(app: &App) -> String {
    let Some(job) = &app.job else {
        return String::new();
    };
    let Ok(log) = job.log.lock() else {
        return "running".into();
    };
    match &log.ended {
        None => "running".into(),
        Some(Ok(())) => "done".into(),
        Some(Err(Failed::Stopped)) => "stopped".into(),
        Some(Err(_)) => "failed".into(),
    }
}

/// `area` framed in a panel; the inside, a column clear of each side.
fn framed(frame: &mut Frame, area: Rect, title: String, focused: bool) -> Rect {
    let block = panel(title, focused);
    let inner = block.inner(area);
    frame.render_widget(block, area);
    padded(inner)
}

/// `r` less a column each side.
fn padded(r: Rect) -> Rect {
    if r.width < 3 {
        return r;
    }
    Rect::new(r.x + 1, r.y, r.width.saturating_sub(2), r.height)
}

/// Help text as paragraphs: clap's single line breaks are soft.
fn reflow(text: &str) -> Vec<Line<'static>> {
    text.split("\n\n")
        .map(|para| Line::raw(para.split_whitespace().collect::<Vec<_>>().join(" ")))
        .collect()
}

fn draw_inputs(app: &App, frame: &mut Frame, body: Rect) {
    let [left, right] =
        Layout::horizontal([Constraint::Fill(1), Constraint::Length(INPUTS_WIDTH)]).areas(body);
    let inputs = &app.inputs;
    let bams = &inputs.bams;
    let focus = inputs.focus;

    let title = format!(" BAMs · {} {}", tilde(&bams.cwd), bams.list.find.tag());
    let inner = framed(frame, left, title, focus == InputsFocus::Bams);
    app.hits.add(left, Hit::Bams(None));
    let (at, len) = (bams.list.at, bams.list.len());
    let first = first_visible(at, len, inner.height as usize);
    app.hits.rows(inner, first, len, |i| Hit::Bams(Some(i)));
    let batches: std::collections::HashMap<_, _> = inputs.batch_names().into_iter().collect();
    let name_w = bams
        .list
        .shown()
        .filter(|e| !e.dir)
        .map(|e| e.name.chars().count())
        .max()
        .unwrap_or(0)
        .max(16);
    let batch_w = batches
        .values()
        .map(|b| b.chars().count())
        .max()
        .unwrap_or(0)
        .max(10);
    let lines: Vec<Line> = bams
        .list
        .shown()
        .enumerate()
        .map(|(i, e)| {
            let mark = marker(i == bams.list.at && focus == InputsFocus::Bams);
            if e.dir {
                return Line::from(vec![mark, Span::raw(format!("{}/", e.name))]);
            }
            let tag = match inputs.role_of(&e.path) {
                Some(Role::Fg) => Span::styled("[fg] ", HIGHLIGHT),
                Some(Role::Bg) => Span::styled("[bg] ", HIGHLIGHT),
                None => Span::styled("[  ] ", DIM),
            };
            let batch = batches
                .get(&e.path)
                .map_or_else(String::new, |b| format!("batch {b}"));
            let index = if e.tagged {
                Span::styled("indexed", DIM)
            } else {
                Span::styled("no index", ACCENTED)
            };
            Line::from(vec![
                mark,
                tag,
                Span::raw(format!("{:<name_w$}  ", e.name)),
                Span::styled(format!("{batch:<w$}  ", w = batch_w + 6), DIM),
                index,
            ])
        })
        .collect();
    let rows = inner.height as usize;
    let lines = if lines.is_empty() {
        vec![Line::styled(bams.list.empty_note(""), DIM)]
    } else {
        window(lines, bams.list.at, rows)
    };
    frame.render_widget(Paragraph::new(lines), inner);

    let inner = framed(frame, right, " inputs ".into(), focus == InputsFocus::Rows);
    app.hits.add(right, Hit::Rows(None));
    app.hits
        .rows(inner, 0, Row::ALL.len(), |i| Hit::Rows(Some(i)));
    let out = inputs.output();
    let out_path = Path::new(&out);
    let state = if out.is_empty() {
        OutState::Unnamed
    } else {
        app.checked.output
    };
    let out_note = match state {
        OutState::Unnamed => Span::styled("(none: Enter to name it)", DIM),
        OutState::File => Span::styled("  a file", HIGHLIGHT),
        OutState::NotEmpty => Span::styled("  not empty", HIGHLIGHT),
        OutState::Empty => Span::styled("  empty", DIM),
        OutState::New => Span::styled("  new", DIM),
    };
    let threads = app
        .form
        .get("max-threads")
        .map_or_else(String::new, |f| f.shown());
    let mut lines: Vec<Line> = Row::ALL
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut spans = vec![
                marker(i == inputs.row && focus == InputsFocus::Rows),
                Span::raw(format!("{:<12}", row.label())),
            ];
            match row {
                Row::File(f) => spans.push(match f.get(inputs) {
                    Some(p) => Span::raw(tilde(p)),
                    None if f.downloadable() => Span::styled("(none; d downloads)", DIM),
                    None => Span::styled("(none)", DIM),
                }),
                Row::Output => spans.extend([Span::raw(tilde(out_path)), out_note.clone()]),
                Row::Threads => spans.push(Span::raw(threads.clone())),
            }
            Line::from(spans)
        })
        .collect();
    if app.fetch.is_some() {
        lines.push(Line::styled("downloading… (d shows it)", HIGHLIGHT));
    }
    lines.extend([
        Line::raw(""),
        Line::styled("fg: signal, passed as BAMs", DIM),
        Line::styled("bg: control, passed as --control-bam", DIM),
    ]);
    frame.render_widget(Paragraph::new(lines), inner);
}

fn draw_steps(app: &App, frame: &mut Frame, body: Rect) {
    let inner = framed(frame, body, " steps, in order ".into(), true);
    let steps = &app.steps;
    let bg = app.inputs.bg().len();
    let effective = steps.effective(bg > 0);
    let lines: Vec<Line> = Step::ALL
        .iter()
        .enumerate()
        .map(|(i, &s)| {
            let blocked = s.available(bg > 0).is_err();
            let on = effective[i];
            let style = if on { PLAIN } else { DIM };
            let mut spans = vec![
                marker(i == steps.at),
                Span::styled(if on { "[x] " } else { "[ ] " }, style),
                Span::styled(format!("{:<8}", s.label()), style),
                Span::styled(s.about(), style),
            ];
            let extra = match s {
                Step::Depth if steps.is_on(s) => {
                    let kb = steps.depth_kb.trim();
                    if kb.is_empty() {
                        Span::styled("   resolution: ? (Enter to give kb)", HIGHLIGHT)
                    } else {
                        Span::styled(format!("   resolution: {kb} kb"), ACCENTED)
                    }
                }
                Step::Depth => Span::styled("   resolution: off (Space asks for kb)", DIM),
                _ if blocked => Span::styled("   needs a bg BAM", DIM),
                Step::M6a => Span::styled(
                    format!("   {bg} bg BAM{}", if bg == 1 { "" } else { "s" }),
                    DIM,
                ),
                _ => Span::raw(""),
            };
            spans.push(extra);
            Line::from(spans)
        })
        .collect();
    let rows = inner.height as usize;
    let first = first_visible(steps.at, lines.len(), rows);
    app.hits.rows(inner, first, lines.len(), Hit::Step);
    frame.render_widget(Paragraph::new(window(lines, steps.at, rows)), inner);
}

fn draw_flags(app: &App, frame: &mut Frame, body: Rect) {
    let form = &app.form;
    let mut title = format!(" flags · {} changed ", form.changed());
    if app.advanced {
        title.push_str("· advanced ");
    }
    if !app.find.is_empty() {
        title.push_str(&format!("· find {} ", app.find));
    }
    let inner = framed(frame, body, title, true);
    let [list, help] = Layout::vertical([Constraint::Fill(1), Constraint::Length(6)]).areas(inner);

    let visible = app.visible_flags();
    let sel = app.flag_row();
    let complaint = &app.checked.complaint;
    let blamed = complaint
        .as_deref()
        .and_then(|c| form::blamed(c, &form.fields));
    let has_bg = !app.inputs.bg().is_empty();

    let mut lines: Vec<Line> = Vec::new();
    // Each line's flag, by its place among those shown; none for headings.
    let mut flag_of: Vec<Option<usize>> = Vec::new();
    let mut sel_line = 0;
    for h in form.headings() {
        let rows: Vec<usize> = visible
            .iter()
            .copied()
            .filter(|&i| form.fields.get(i).is_some_and(|f| f.heading == h))
            .collect();
        if rows.is_empty() {
            continue;
        }
        let dim = !app.steps.heading_on(&h, has_bg);
        let name = if h.is_empty() { "other" } else { h.as_str() };
        lines.push(if dim {
            Line::styled(format!("{name} (step off)"), DIM)
        } else {
            Line::styled(name.to_string(), HIGHLIGHT)
        });
        flag_of.push(None);
        for i in rows {
            let Some(f) = form.fields.get(i) else {
                continue;
            };
            if Some(i) == sel {
                sel_line = lines.len();
            }
            let changed = !f.is_default();
            let style = if dim { DIM } else { PLAIN };
            let shown = match f.kind {
                Kind::Choice(_) => format!("‹ {} ›", f.shown()),
                _ => f.shown(),
            };
            let mut spans = vec![
                marker(Some(i) == sel),
                Span::styled(if changed { "* " } else { "  " }, ACCENTED),
                Span::styled(format!("{:<24}", f.long), style),
                Span::styled(format!("{shown:<16}"), style),
            ];
            if changed {
                let d = if f.default.is_empty() {
                    "(unset)"
                } else {
                    &f.default
                };
                spans.push(Span::styled(format!("  (default {d})"), DIM));
            }
            if blamed == Some(f.long.as_str()) {
                let c = complaint.as_deref().unwrap_or_default();
                spans.push(Span::styled(format!("  {c}"), ACCENTED));
            }
            lines.push(Line::from(spans));
            flag_of.push(visible.iter().position(|&v| v == i));
        }
    }
    let first = first_visible(sel_line, lines.len(), list.height as usize);
    for (k, line) in (first..lines.len()).take(list.height as usize).enumerate() {
        if let Some(pos) = flag_of[line] {
            let row = Rect::new(list.x, list.y + k as u16, list.width, 1);
            app.hits.add(row, Hit::Flag(pos));
        }
    }
    if lines.is_empty() {
        lines.push(Line::styled(
            format!("no flag matches \"{}\"", app.find),
            DIM,
        ));
    }
    let rows = list.height as usize;
    frame.render_widget(Paragraph::new(window(lines, sel_line, rows)), list);

    let Some(f) = sel.and_then(|i| form.fields.get(i)) else {
        return;
    };
    let block = Block::new()
        .borders(Borders::TOP)
        .border_style(DIM)
        .title(Line::styled(format!(" --{} ", f.long), HIGHLIGHT));
    let text = block.inner(help);
    frame.render_widget(block, help);
    let d = if f.default.is_empty() {
        "(none)"
    } else {
        &f.default
    };
    let mut lines = reflow(&f.long_help);
    lines.push(Line::styled(format!("default {d}"), DIM));
    frame.render_widget(Paragraph::new(lines).wrap(Wrap { trim: true }), text);
}

fn draw_run(app: &App, frame: &mut Frame, body: Rect) {
    let Some(job) = &app.job else {
        let inner = framed(frame, body, " run ".into(), true);
        let line = Line::styled("no run yet: preview it with ".to_string() + APPLY_KEYS, DIM);
        frame.render_widget(line, inner);
        return;
    };
    let Ok(log) = job.log.lock() else {
        return;
    };
    let title = match &log.ended {
        None => {
            let step = log
                .step
                .strip_prefix("Step ")
                .map_or_else(|| "starting".into(), |s| format!("step {s}"));
            format!(" {step} · {} ", clock(job.started.elapsed()))
        }
        Some(Ok(())) => " done ".into(),
        Some(Err(Failed::Stopped)) => " stopped ".into(),
        Some(Err(Failed::Start(r) | Failed::Exit(r))) => format!(" exit: {r} "),
    };
    let inner = framed(frame, body, title, true);
    let bar = log
        .ended
        .is_none()
        .then_some(log.progress.as_ref())
        .flatten();
    let rows = (inner.height as usize).saturating_sub(usize::from(bar.is_some()));
    job.rows.set(rows);
    let len = log.lines.len();
    let scroll = job.scroll.min(len.saturating_sub(rows));
    let end = len.saturating_sub(scroll);
    let start = end.saturating_sub(rows);
    let mut lines: Vec<Line> = log
        .lines
        .iter()
        .skip(start)
        .take(end - start)
        .map(|l| Line::raw(l.clone()))
        .collect();
    if let Some(p) = bar {
        let filled = (p.len > 0).then(|| filled(p.pos, p.len, BAR));
        lines.push(match filled {
            Some(filled) => Line::from(vec![
                Span::styled("█".repeat(filled), ACCENTED),
                Span::styled("░".repeat(BAR.saturating_sub(filled)), DIM),
                Span::raw(format!(" {}/{} {}", p.pos, p.len, p.what)),
            ]),
            // A spinner: no length to fill.
            None => Line::from(vec![
                Span::styled("⠿ ", ACCENTED),
                Span::raw(p.what.clone()),
            ]),
        });
    }
    frame.render_widget(Paragraph::new(lines), inner);
}

fn draw_picking(app: &App, frame: &mut Frame, body: Rect) {
    let Some((row, b)) = &app.picking else {
        return;
    };
    let label = row.label();
    let h = u16::try_from(b.list.len() + 2).unwrap_or(u16::MAX).max(4);
    let h = h.min(body.height.saturating_sub(2).max(3));
    let title = format!(" {label} · {} {}", tilde(&b.cwd), b.list.find.tag());
    let inner = popup_frame(frame, body, 64, h, title);
    let lines = b.lines(
        inner.height as usize,
        inner.width as usize,
        "",
        "nothing to pick here",
    );
    frame.render_widget(Paragraph::new(lines), inner);
}

/// The pop-up telling how the run ended, how long it took, and where its
/// output is.
fn draw_ended(app: &App, frame: &mut Frame, body: Rect) {
    let Some(job) = &app.job else {
        return;
    };
    let ended = job.log.lock().ok().and_then(|l| l.ended.clone());
    let Some(ended) = ended else {
        return;
    };
    let title = format!(" faba run {} ", ended_word(&ended));
    let why = match ended {
        Err(Failed::Start(e) | Failed::Exit(e)) => Some(e),
        _ => None,
    };
    let took = job.took.map_or_else(String::new, clock);
    let mut lines = vec![
        Line::from(vec![Span::styled(" took  ", DIM), Span::raw(took)]),
        Line::from(vec![
            Span::styled(" into  ", DIM),
            Span::raw(tilde(&job.out)),
        ]),
    ];
    if let Some(why) = why {
        lines.push(Line::raw(""));
        lines.push(Line::styled(format!(" {why}"), HIGHLIGHT));
    }
    lines.push(Line::raw(""));
    lines.push(Line::styled(" its log stays on the Run screen", DIM));
    button_popup(frame, body, &title, &[("ok", Hit::Ok)], lines, &app.hits);
}

/// A download's pop-up: what, where, the step and how far through it.
fn draw_fetch(f: &Fetch, frame: &mut Frame, body: Rect) {
    let (step, bytes) = f
        .status
        .lock()
        .map_or_else(|_| (String::new(), None), |s| (s.step.clone(), s.bytes));
    let mb = |b: u64| b as f64 / f64::from(1 << 20);
    let progress = match bytes {
        Some((got, Some(total))) if total > 0 => {
            let w = 36;
            let full = filled(got, total, w);
            format!(
                "{}{} {:>3.0}%  {:.0} / {:.0} MB",
                "█".repeat(full),
                "░".repeat(w - full),
                got.min(total) as f64 / total as f64 * 100.0,
                mb(got),
                mb(total)
            )
        }
        Some((got, _)) => format!("{:.0} MB", mb(got)),
        None => String::new(),
    };
    let lines = vec![
        Line::styled(f.label.clone(), PLAIN),
        Line::styled(format!("into {}", tilde(&f.dir)), DIM),
        Line::raw(""),
        Line::styled(step, HIGHLIGHT),
        Line::styled(progress, HIGHLIGHT),
        Line::raw(""),
        Line::styled("the GTF, then the genome: fetched, unpacked, indexed", DIM),
    ];
    popup(frame, body, " downloading a reference ", lines);
}

/// The download pop-up: the sources, or one source's releases or species.
fn draw_catalogue(c: &Catalogue, frame: &mut Frame, body: Rect) {
    let rows = usize::from(body.height.saturating_sub(10)).clamp(3, 16);
    let mut lines: Vec<Line> = Vec::new();
    let title = match c.source {
        None => {
            lines.extend(SOURCES.iter().enumerate().map(|(i, (name, s))| {
                let fits = match c.bam_chr {
                    Some(chr) if chr == s.chr_named() => "  matches your BAMs",
                    _ => "",
                };
                Line::from(vec![
                    marker(i == c.at),
                    Span::raw(*name),
                    Span::styled(fits, DIM),
                ])
            }));
            " download a reference ".to_string()
        }
        Some(s) => {
            let what = match s {
                Source::Gencode(_) => "release",
                Source::Ensembl => "species",
            };
            if c.loading() {
                lines.push(Line::styled("reading the FTP listing…", DIM));
            } else if let Some(e) = &c.error {
                lines.push(Line::styled(format!("cannot list: {e}"), HIGHLIGHT));
            } else if c.items.is_empty() {
                lines.push(Line::styled(c.items.empty_note(""), DIM));
            } else {
                let items =
                    c.items.shown().enumerate().map(|(i, n)| {
                        Line::from(vec![marker(i == c.items.at), Span::raw(n.clone())])
                    });
                lines.extend(window(items.collect(), c.items.at, rows));
            }
            format!(" choose a {what} {}", c.items.find.tag())
        }
    };
    lines.push(Line::raw(""));
    lines.push(Line::styled(
        "GENCODE names chromosomes chr1, Ensembl 1: take the one your BAMs use",
        DIM,
    ));
    if let Some(chr) = c.bam_chr {
        let named = if chr { "chr1" } else { "1" };
        lines.push(Line::styled(
            format!("your first BAM names them {named}"),
            DIM,
        ));
    }
    lines.push(Line::styled(
        "the GTF and genome come to about 1 GB, 3 GB unpacked",
        DIM,
    ));
    popup(frame, body, &title, lines);
}

fn draw_preview(app: &App, frame: &mut Frame, body: Rect) {
    // What blocks the run and where the script goes come first: a long list
    // of BAMs must not push them out of the pop-up.
    let problems = &app.checked.problems;
    let mut lines: Vec<Line> = if problems.is_empty() {
        vec![Line::styled(" no problems", DIM)]
    } else {
        problems
            .iter()
            .map(|p| Line::styled(format!(" {p}"), HIGHLIGHT))
            .collect()
    };
    let out = app.inputs.output();
    if !out.is_empty() {
        lines.push(Line::styled(
            format!(" saved as {}/{}", tilde(Path::new(&out)), script::SCRIPT),
            DIM,
        ));
    }
    lines.push(Line::raw(""));
    // The folder as named: the script itself runs there, as `-o .`.
    let argv = app.argv_out();
    lines.extend(
        script::command_lines(&argv)
            .into_iter()
            .enumerate()
            .map(|(k, l)| {
                let l = if k == 0 {
                    l.replacen("\"${FABA:-faba}\"", "faba", 1)
                } else {
                    format!("  {l}")
                };
                Line::styled(format!(" {l}"), PLAIN)
            }),
    );
    let buttons = [("▶ start", Hit::Start), ("cancel", Hit::Cancel)];
    button_popup(frame, body, " start this run? ", &buttons, lines, &app.hits);
}

fn footer_line(app: &App) -> Line<'static> {
    if let Some((target, line)) = &app.editing {
        let prompt = match target {
            Target::Output => match &app.note {
                Some(why) => format!("{why}; output folder: "),
                None => "output folder (new or empty): ".to_string(),
            },
            Target::RefDir => "download into: ".to_string(),
            Target::DepthKb => "depth resolution (kb): ".to_string(),
            Target::Flag(i) => app
                .form
                .fields
                .get(*i)
                .map_or_else(String::new, |f| format!("--{}: ", f.long)),
            Target::Find => "find flag: ".to_string(),
        };
        let keys: &[(&str, &str)] = match target {
            Target::Flag(_) => &[("Enter", "set (empty: default)"), ("Esc", "back")],
            _ => &[("Enter", "set"), ("Esc", "back")],
        };
        return input_line(&prompt, line.text().unwrap_or(""), keys);
    }
    if let Some(note) = &app.note {
        return Line::styled(format!(" {note}"), HIGHLIGHT);
    }
    if app.ended_shown {
        return help_line(&[("Enter/Esc", "close, onto the log")]);
    }
    if app.fetch_shown && app.fetch.is_some() {
        return help_line(&[
            ("Esc", "hide; it goes on (d shows it)"),
            ("s", "stop, keeping what came"),
        ]);
    }
    if let Some(c) = &app.catalogue {
        return if c.source.is_none() {
            help_line(&[("↑/↓", "move"), ("Enter", "list it"), ("Esc", "back")])
        } else {
            help_line(&[
                ("type", "find"),
                ("↑/↓", "move"),
                ("Enter", "download"),
                ("Esc", "back"),
            ])
        };
    }
    if let Some((_, b)) = &app.picking {
        let esc = if b.list.find.text.is_empty() {
            "back"
        } else {
            "clear find"
        };
        return help_line(&[
            ("type", "find"),
            ("↑/↓", "move"),
            ("Enter", "pick/open"),
            ("←", "up"),
            ("Esc", esc),
        ]);
    }
    if app.preview {
        return help_line(&[
            (APPLY_KEYS, "save the script and start"),
            ("c", "copy the command"),
            ("p", "print it and leave"),
            ("Esc", "back"),
        ]);
    }
    let mut common = vec![
        ("Tab", "screen"),
        (APPLY_KEYS, "preview"),
        ("c", "copy the command"),
    ];
    if !app.running() {
        common.extend([("p", "print the command and leave"), ("q", "quit")]);
    }
    if app.page == Page::Inputs && app.inputs.focus == InputsFocus::Bams {
        // Letters type into the find here, so only these keys act.
        let mut keys = vec![
            ("type", "find, / and ~ go to a path"),
            ("Space", "fg → bg → off"),
            ("Enter/→", "open"),
            ("←", "up"),
        ];
        if !app.inputs.bams.list.find.text.is_empty() {
            keys.push(("Esc", "clear find"));
        }
        keys.extend([
            ("Tab", "inputs"),
            (APPLY_KEYS, "preview"),
            ("Ctrl-C", "quit"),
        ]);
        return help_line(&keys);
    }
    let mut keys: Vec<(&str, &str)> = match app.page {
        Page::Inputs => vec![
            ("↑/↓", "row"),
            ("Enter", "pick/type"),
            ("d", "download GFF + genome"),
            ("h/Esc", "BAMs"),
        ],
        Page::Steps => vec![("Space", "toggle"), ("↑/↓", "move")],
        Page::Flags => vec![
            ("↑/↓", "move"),
            ("Space", "on/off"),
            ("←/→", "choice"),
            ("Enter", "type"),
            ("/", "find"),
            ("a", "advanced"),
            ("r/R", "reset"),
        ],
        Page::Run => {
            let Some(job) = &app.job else {
                return help_line(&common);
            };
            if job.running() {
                if job.asking {
                    return help_line(&[("s", "again stops the run"), ("Esc", "keep it going")]);
                }
                return help_line(&[
                    ("↑/↓", "scroll"),
                    ("End", "follow"),
                    ("s", "stop (asks; again to kill)"),
                    ("q", "is refused while running"),
                ]);
            }
            vec![("↑/↓", "scroll"), ("End", "follow")]
        }
    };
    keys.extend(common);
    help_line(&keys)
}
