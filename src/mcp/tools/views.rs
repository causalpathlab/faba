//! `view_open`, `view_list`, `view_state` and `view_export`: a view put in
//! front of the user as a detached session they attach to, what they
//! arranged in it read back, and figures drawn from exactly that (see
//! [`crate::tui::live`]).

use crate::mcp::args::{known_keys, push_flag, with_typed, Arguments};
use crate::mcp::run::{absolute, exists, svg_png};
use crate::tui::live::{self, Look, Request};
use crate::tui::TOOL;
use anyhow::{bail, ensure, Context, Result};
use legume_mcp::{Call, Content};
use legume_tui::Listed;
use serde_json::{json, Value};
use std::path::Path;
use std::time::{Duration, Instant};

/// The views `view_open` opens: each by its tool's name, its command and
/// its session's screen, and whether it answers `view_state` and
/// `view_export` yet.
pub(in crate::mcp) const VIEWS: [(&str, &str, bool); 5] = [
    ("qc_view", "qc-view", true),
    ("cluster_view", "cluster-view", false),
    ("pileup", "pileup", false),
    ("metagene", "metagene", false),
    ("pwm", "pwm", false),
];

/// The arguments of a view's tool that write or print in place of the
/// view, refused by `view_open`.
fn batch_only(view: &str) -> &'static [&'static str] {
    match view {
        "qc_view" => &["html", "plots", "out", "figure", "format"],
        "cluster_view" => &["html", "figure", "plot", "plots"],
        "pileup" => &["figure", "out"],
        "metagene" => &["out", "dist_measures", "figure"],
        "pwm" => &["out", "figure"],
        _ => &[],
    }
}

/// How long `view_open` waits for a view that answers to come up.
const OPENING: Duration = Duration::from_secs(60);

/// `view`'s command line from its tool's `args`, for the view itself: no
/// `--batch-process`, nothing written in its place.
fn view_argv(view: &str, a: &Arguments) -> Result<Vec<String>> {
    let (_, command, _) = (VIEWS.iter()).find(|v| v.0 == view).with_context(|| {
        format!("no view {view}: one of qc_view, cluster_view, pileup, metagene, pwm")
    })?;
    for key in batch_only(view) {
        ensure!(
            a.0.get(*key).is_none_or(Value::is_null),
            "`{key}` writes in place of the view: leave it out of `args` (view_export draws what the view shows)"
        );
    }
    let mut argv = vec![command.to_string()];
    match view {
        "qc_view" => argv.push(a.path("input")?),
        "cluster_view" => argv.push(a.path("data")?),
        _ => argv.extend(a.paths_req("data")?),
    }
    match view {
        "pileup" => {
            let genes = a.strings("genes")?;
            if !genes.is_empty() {
                push_flag(&mut argv, "-q", &genes.join(","));
            }
            let regions = a.strings("regions")?;
            if !regions.is_empty() {
                push_flag(&mut argv, "--regions", &regions.join(","));
            }
            if let Some(gtf) = a.path_opt("gtf")? {
                push_flag(&mut argv, "--gtf", &gtf);
            }
            let bam = a.paths("bam")?;
            if !bam.is_empty() {
                argv.push("--bam".into());
                argv.extend(bam);
            }
        }
        "metagene" => {
            if let Some(gff) = a.path_opt("gff")? {
                push_flag(&mut argv, "-g", &gff);
            }
        }
        "pwm" => {
            if let Some(genome) = a.path_opt("genome")? {
                push_flag(&mut argv, "-f", &genome);
            }
        }
        _ => {}
    }
    argv.extend(a.options(view)?);
    Ok(argv)
}

/// The command a user types to attach to `id`: this faba, with the
/// session folder it was started under when that is not the default.
fn attach_line(id: &str) -> String {
    let exe = legume_mcp::own_command(["attach", id]).map_or_else(
        |_| "faba".to_string(),
        |c| c.get_program().to_string_lossy().into_owned(),
    );
    let mut line = format!("{} attach {id}", quoted(&exe));
    if let Some(dir) = std::env::var_os(legume_tui::DIR_VAR).filter(|d| !d.is_empty()) {
        let dir = quoted(&dir.to_string_lossy());
        line = format!("{}={dir} {line}", legume_tui::DIR_VAR);
    }
    line
}

/// `s` as one shell word.
fn quoted(s: &str) -> String {
    let plain = |c: char| c.is_ascii_alphanumeric() || "/._-+=:,@".contains(c);
    match s.chars().all(plain) {
        true => s.to_string(),
        false => format!("'{}'", s.replace('\'', r"'\''")),
    }
}

/// Open a Terminal window on macOS running `line`.
fn open_window(line: &str) -> Result<()> {
    ensure!(
        cfg!(target_os = "macos"),
        "`open_window` opens a Terminal window on macOS only: attach with {line}"
    );
    let script = line.replace('\\', r"\\").replace('"', "\\\"");
    let status = std::process::Command::new("osascript")
        .args([
            "-e",
            &format!("tell application \"Terminal\" to do script \"{script}\""),
            "-e",
            "tell application \"Terminal\" to activate",
        ])
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .status()
        .context("running osascript")?;
    ensure!(
        status.success(),
        "osascript failed to open a Terminal window"
    );
    Ok(())
}

/// faba's live sessions.
fn live_sessions() -> Result<Vec<Listed>> {
    let dir = legume_tui::runtime_dir()?;
    Ok(legume_tui::dir::live(&dir, &TOOL.prefix())?)
}

/// The live session `wanted` names: by its label, its id, its id after
/// `faba-`, or the end of its id (`2` for `faba-qc-view-2`).
fn pick(wanted: &str) -> Result<Listed> {
    let live = live_sessions()?;
    let full = format!("{}{wanted}", TOOL.prefix());
    if let Some(l) = (live.iter())
        .find(|l| l.label.as_deref() == Some(wanted))
        .or_else(|| live.iter().find(|l| l.name == wanted || l.name == full))
    {
        return Ok(l.clone());
    }
    let tail = format!("-{wanted}");
    let ends: Vec<&Listed> = live.iter().filter(|l| l.name.ends_with(&tail)).collect();
    match ends[..] {
        [one] => Ok(one.clone()),
        [] => {
            let names: Vec<&str> = live.iter().map(|l| l.name.as_str()).collect();
            bail!("no live session {wanted}; live: {}", names.join(", "))
        }
        _ => bail!(
            "{wanted} could be any of: {}",
            ends.iter()
                .map(|l| l.name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// What a listing says of session `l`.
fn listed(l: &Listed) -> Value {
    let answers = live::control_path(&l.name).is_ok_and(|p| p.exists());
    json!({
        "session": l.name,
        "label": l.label,
        "attached": l.attached,
        "started": l.age(),
        "status": l.status,
        "command": l.command,
        "answers": answers,
        "attach": attach_line(&l.name),
    })
}

/// What ended session `id` left to say, if anything.
fn left(id: &str) -> Option<String> {
    let path = legume_tui::runtime_dir().ok()?.join(format!("{id}.leave"));
    let text = std::fs::read_to_string(path).ok()?;
    Some(text.trim().to_string()).filter(|t| !t.is_empty())
}

pub(in crate::mcp) fn view_open(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let view = a.str("view").context("`view` is required")?.to_string();
    let args = a.0.get("args").cloned().unwrap_or_else(|| json!({}));
    ensure!(
        args.is_object(),
        "`args` is an object: the arguments of the tool `view` names"
    );
    known_keys(&view, &args)?;
    let va = Arguments(with_typed(&view, args)?);
    let argv = view_argv(&view, &va)?;
    let (_, _, answers) = *VIEWS.iter().find(|v| v.0 == view).expect("checked");
    let screen = argv[0].clone();
    let mut program = legume_mcp::own_command(&argv)?;
    if std::env::var_os("TERM").is_none_or(|t| t.is_empty() || t == "dumb") {
        program.env("TERM", "xterm-256color");
    }
    let label = a.str("session");
    let id = TOOL
        .start_detached(&screen, label, program)
        .with_context(|| format!("starting the {view} view"))?;
    let attach = attach_line(&id);
    let window = a.flag("open_window")?;
    if window {
        open_window(&attach)?;
    }
    let mut v = json!({
        "session": id,
        "attach": attach,
        "window": window,
        "command": argv,
        "answers": answers,
        "note": format!(
            "the view runs detached until the user attaches: ask them to run `{attach}` in a terminal{}; \
             Ctrl+Alt+D detaches again, q leaves the view (and ends the session)",
            if window { " (a Terminal window is open on it)" } else { "" }
        ),
    });
    if !answers {
        v["note"] = json!(format!(
            "{}; this view does not answer view_state or view_export yet",
            v["note"].as_str().unwrap_or_default()
        ));
        return Ok(vec![Content::Json(v)]);
    }
    // Up once it answers: the view has read its inputs and is on screen.
    let start = Instant::now();
    loop {
        match live::ask(&id, &Request::State) {
            Ok(state) => {
                v["state"] = state;
                break;
            }
            Err(e) => {
                // An error it ended with is held on its screen for whoever
                // attaches; said here instead, and the session stopped.
                if let Some(why) = left(&id) {
                    let _ = TOOL.stop(&id, false, Duration::from_secs(5));
                    bail!("the {view} view ended as it opened: {why}");
                }
                if TOOL.find(&id).ok().flatten().is_none() {
                    bail!("the {view} view ended as it opened: {e:#}");
                }
                if start.elapsed() > OPENING || cx.is_cancelled() {
                    v["note"] = json!(format!(
                        "{}; still opening (reading its inputs): view_state once it is up",
                        v["note"].as_str().unwrap_or_default()
                    ));
                    break;
                }
                std::thread::sleep(Duration::from_millis(200));
            }
        }
    }
    Ok(vec![Content::Json(v)])
}

pub(in crate::mcp) fn view_list(_a: &Arguments, _cx: &Call) -> Result<Vec<Content>> {
    let live = live_sessions()?;
    let sessions: Vec<Value> = live.iter().map(listed).collect();
    Ok(vec![Content::Json(json!({ "sessions": sessions }))])
}

pub(in crate::mcp) fn view_state(a: &Arguments, _cx: &Call) -> Result<Vec<Content>> {
    let l = pick(a.str("session").context("`session` is required")?)?;
    let request = match a.0.get("set").filter(|s| !s.is_null()) {
        Some(state) => {
            ensure!(
                state.is_object(),
                "`set` is an object: a state as view_state returns it"
            );
            Request::Set {
                state: state.clone(),
            }
        }
        None => Request::State,
    };
    let state = live::ask(&l.name, &request)?;
    let mut v = listed(&l);
    v["state"] = state;
    Ok(vec![Content::Json(v)])
}

/// The flags of a figure's size and style, left out of a replay's
/// command so the export's own are the ones it takes.
const LOOK_FLAGS: [&str; 5] = [
    "--fig-width",
    "--fig-height",
    "--font-size",
    "--panel-only",
    "--palette",
];

/// The flags of a view's session and of what it writes in place of the
/// view, left out of a replay's command.
const SESSION_FLAGS: [&str; 7] = [
    "--session",
    "--no-detach",
    "--view-state",
    "--figure",
    "--format",
    "--html",
    "--plots",
];

/// The batch command that draws the figure `view_export` wrote, without
/// the view: the view's own command line (`command`), less its session
/// and look flags, with `--view-state` and the look the export drew at
/// (`sheet`).
fn replay(
    command: &[String],
    state_file: &str,
    prefix: &str,
    formats: &[String],
    sheet: &Value,
) -> Vec<String> {
    let takes_value = |f: &str| !matches!(f, "--no-detach" | "--panel-only" | "--plots");
    let mut argv = vec!["faba".to_string()];
    let mut skip = false;
    for w in command {
        if std::mem::take(&mut skip) {
            continue;
        }
        let name = w.split('=').next().unwrap_or(w);
        if LOOK_FLAGS.contains(&name) || SESSION_FLAGS.contains(&name) {
            skip = !w.contains('=') && takes_value(name);
            continue;
        }
        argv.push(w.clone());
    }
    argv.extend([
        format!("--view-state={state_file}"),
        format!("--figure={prefix}"),
    ]);
    argv.push(format!("--format={}", formats.join(",")));
    if let Some(w) = sheet["width"].as_f64() {
        argv.push(format!("--fig-width={w}"));
    }
    if let Some(h) = sheet["height"].as_f64() {
        argv.push(format!("--fig-height={h}"));
    }
    if let Some(pt) = sheet["font_size"].as_f64() {
        argv.push(format!("--font-size={pt}"));
    }
    if sheet["panel_only"] == json!(true) {
        argv.push("--panel-only".into());
    }
    if let Some(p) = sheet["palette"].as_str() {
        argv.push(format!("--palette={p}"));
    }
    argv
}

pub(in crate::mcp) fn view_export(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let l = pick(a.str("session").context("`session` is required")?)?;
    let formats: Vec<crate::figure::Format> = a
        .parse("formats")?
        .unwrap_or_else(|| vec![crate::figure::Format::Pdf, crate::figure::Format::Png]);
    ensure!(
        !formats.is_empty(),
        "`formats` names at least one of svg, pdf, png"
    );
    let taken = |p: &str| {
        let state = live::state_file(p);
        exists(&state) || (formats.iter()).any(|f| exists(&format!("{p}.{}", f.extension())))
    };
    let prefix = absolute(&a.out_or(cx, (&format!("{}.figure", l.name), ""), taken)?);
    ensure!(
        !taken(&prefix),
        "{prefix}: a figure is there; not replacing it: choose another `out`"
    );
    if let Some(dir) = Path::new(&prefix).parent() {
        ensure!(dir.is_dir(), "{}: no such folder for `out`", dir.display());
    }
    let look = Look {
        width: a.parse("width")?,
        height: a.parse("height")?,
        font_size: a.parse("font_size")?,
        panel_only: a.parse("panel_only")?,
        palette: a.parse("palette")?,
    };
    let request = Request::Export {
        prefix: prefix.clone(),
        formats: formats.clone(),
        look,
    };
    let mut v = live::ask(&l.name, &request)?;
    let command: Vec<String> =
        serde_json::from_value(v["state"]["command"].clone()).unwrap_or_default();
    let state_file = v["state_file"].as_str().unwrap_or_default().to_string();
    let names: Vec<String> = formats.iter().map(|f| f.extension().to_string()).collect();
    v["replay"] = json!(replay(&command, &state_file, &prefix, &names, &v["sheet"]));
    v["session"] = json!(l.name);
    let mut content = vec![Content::Json(v.clone())];
    let png = format!("{prefix}.png");
    let svg = format!("{prefix}.svg");
    if exists(&png) {
        content.push(Content::Png(std::fs::read(&png)?));
    } else if exists(&svg) {
        content.push(Content::Png(svg_png(&svg)?));
    }
    Ok(content)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A view's command line is its tool's, for the view: the typed
    /// arguments as flags, nothing written in its place, no batch flag.
    #[test]
    fn a_view_opens_on_its_tools_arguments() {
        let argv = |view: &str, args: Value| {
            let a = Arguments(with_typed(view, args).unwrap());
            view_argv(view, &a)
        };
        let qc = argv(
            "qc_view",
            json!({ "input": "runs/qc", "scale": ["pv=log"], "front": "mut" }),
        )
        .unwrap();
        assert_eq!(qc[..2], ["qc-view", "runs/qc"]);
        for w in [&["--scale", "pv=log"][..], &["--front", "mut"]] {
            assert!(qc.windows(w.len()).any(|x| x == w), "{w:?} in {qc:?}");
        }
        assert!(!qc.iter().any(|w| w == "--batch-process"), "{qc:?}");
        let line = std::iter::once("faba".to_string()).chain(qc);
        assert!(crate::matches(&mut crate::faba_command(), line).is_ok());
        let e = argv("qc_view", json!({ "input": "runs/qc", "html": "p.html" })).unwrap_err();
        assert!(e.to_string().contains("in place of the view"), "{e}");
        let pileup = argv(
            "pileup",
            json!({ "data": ["out/"], "genes": ["GENE1", "GENE2"] }),
        )
        .unwrap();
        assert_eq!(pileup, ["pileup", "out/", "-q", "GENE1,GENE2"]);
        assert!(argv("browser", json!({})).is_err());
    }

    /// A replay draws the export without the view: the view's command less
    /// its session and look, the state and the look the export drew at.
    #[test]
    fn a_replay_is_the_views_command_with_the_state() {
        let command: Vec<String> = [
            "qc-view",
            "runs/qc",
            "--session",
            "S1",
            "--font-size",
            "6",
            "--front=mut",
            "--panel-only",
        ]
        .map(String::from)
        .to_vec();
        let sheet = json!({ "width": 3.5, "height": null, "font_size": 7.0, "panel_only": false, "palette": "gray" });
        let r = replay(&command, "f.view.json", "f", &["svg".into()], &sheet);
        assert_eq!(
            r,
            [
                "faba",
                "qc-view",
                "runs/qc",
                "--front=mut",
                "--view-state=f.view.json",
                "--figure=f",
                "--format=svg",
                "--fig-width=3.5",
                "--font-size=7",
                "--palette=gray"
            ]
        );
        let line = r.into_iter();
        assert!(crate::matches(&mut crate::faba_command(), line).is_ok());
    }

    #[test]
    fn a_shell_word_is_quoted_only_when_it_needs_it() {
        assert_eq!(quoted("/usr/bin/faba"), "/usr/bin/faba");
        assert_eq!(quoted("/my dir/it's"), r"'/my dir/it'\''s'");
    }
}
