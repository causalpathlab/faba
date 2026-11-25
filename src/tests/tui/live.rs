use super::*;

/// A view of one number, its figure a square that says it.
struct Counter {
    n: i64,
}

impl Live for Counter {
    fn view_name(&self) -> &'static str {
        "counter"
    }

    fn view_state(&self) -> Value {
        json!({ "n": self.n })
    }

    fn set_view_state(&mut self, state: &Value) -> anyhow::Result<()> {
        if let Some(n) = state.get("n") {
            self.n = n.as_i64().context("n is a whole number")?;
        }
        Ok(())
    }

    fn own_sheet(&self) -> (Sheet, Figure) {
        (Sheet::OWN, Figure::Qc)
    }

    fn saved_figure(&self, sheet: Sheet) -> anyhow::Result<String> {
        let w = sheet.page.0.unwrap_or(72.0);
        Ok(crate::figure::svg(w, w, |c| {
            c.text(
                4.0,
                12.0,
                &self.n.to_string(),
                8.0,
                crate::figure::Anchor::Start,
                "#000",
            )
        }))
    }
}

/// `request` asked on `path` from another thread while `view`'s loop
/// serves it, as a tick would.
fn asked(control: &Control, view: &mut Counter, request: Request) -> anyhow::Result<Value> {
    let path = control.path.clone();
    let asking = std::thread::spawn(move || ask_at(&path, "faba-test-1", &request));
    while !asking.is_finished() {
        control.serve(view);
        std::thread::sleep(Duration::from_millis(5));
    }
    asking.join().unwrap()
}

/// The view's loop answers what is asked on its socket: its state, with
/// its name and command line; a state set, refused when it is another
/// view's; a figure written with its state beside it, never over a file;
/// and the socket gone with the control.
#[test]
fn a_live_view_answers_on_its_socket() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("faba-test-1.ctl");
    let control = Control::at(path.clone()).unwrap();
    let mut view = Counter { n: 3 };

    let state = asked(&control, &mut view, Request::State).unwrap();
    assert_eq!(state["n"], 3);
    assert_eq!(state["view"], "counter");
    assert!(state["command"].is_array(), "{state}");

    let set = Request::Set {
        state: json!({ "n": 7, "unknown": true }),
    };
    assert_eq!(asked(&control, &mut view, set).unwrap()["n"], 7);
    assert_eq!(view.n, 7);
    let other = Request::Set {
        state: json!({ "view": "qc_view", "n": 1 }),
    };
    let e = asked(&control, &mut view, other).unwrap_err().to_string();
    assert!(e.contains("qc_view's, not counter's"), "{e}");
    assert_eq!(view.n, 7);

    let prefix = dir.path().join("fig").to_string_lossy().into_owned();
    let export = Request::Export {
        prefix: prefix.clone(),
        formats: vec![Format::Svg, Format::Png],
        look: Look {
            width: Some(2.0),
            ..Look::default()
        },
    };
    let done = asked(&control, &mut view, export.clone()).unwrap();
    let written: Vec<String> = serde_json::from_value(done["written"].clone()).unwrap();
    assert_eq!(written, [format!("{prefix}.svg"), format!("{prefix}.png")]);
    let svg = std::fs::read_to_string(&written[0]).unwrap();
    let sheet = sheet_of(&view, &export_look(&export)).unwrap();
    assert_eq!(svg, view.saved_figure(sheet).unwrap());
    assert_eq!(done["sheet"]["width"], 2.0);
    assert_eq!(done["sheet"]["height"], Value::Null);
    let kept = read_state(&state_file(&prefix)).unwrap();
    assert_eq!(kept["n"], 7);
    assert_eq!(kept["view"], "counter");
    let again = asked(&control, &mut view, export).unwrap_err().to_string();
    assert!(again.contains("not replacing it"), "{again}");

    let line = "not json\n";
    let mut s = UnixStream::connect(&path).unwrap();
    s.write_all(line.as_bytes()).unwrap();
    let mut answer = String::new();
    BufReader::new(&s).read_line(&mut answer).unwrap();
    assert!(answer.contains("not a request"), "{answer}");

    drop(control);
    assert!(!path.exists());
    let e = ask_at(&path, "faba-test-1", &Request::State).unwrap_err();
    assert!(e.to_string().contains("has no view answering"), "{e}");
}

fn export_look(r: &Request) -> Look {
    match r {
        Request::Export { look, .. } => look.clone(),
        _ => Look::default(),
    }
}

/// What an export asks for goes over the view's own sheet, a side or a
/// style at a time, checked as the command line checks it.
#[test]
fn a_look_goes_over_the_views_own() {
    let own = Sheet {
        page: (Some(300.0), None),
        look: crate::figure::style::StyleArgs {
            font_size: Some(7.0),
            panel_only: false,
            palette: Palette::Colour,
        },
    };
    let same = Look::default().over(own, Figure::Qc).unwrap();
    assert_eq!((same.page, same.look), (own.page, own.look));
    let asked = Look {
        height: Some(2.0),
        panel_only: Some(true),
        palette: Some(Palette::Gray),
        ..Look::default()
    };
    let s = asked.over(own, Figure::Qc).unwrap();
    assert_eq!(s.page, (Some(300.0), Some(144.0)));
    assert_eq!(s.look.font_size, Some(7.0));
    assert!(s.look.panel_only);
    assert_eq!(s.look.palette, Palette::Gray);
    let small = Look {
        width: Some(0.5),
        ..Look::default()
    };
    assert!(small.over(own, Figure::Qc).is_err());
    let bad_font = Look {
        font_size: Some(-1.0),
        ..Look::default()
    };
    assert!(bad_font.over(own, Figure::Qc).is_err());
}
