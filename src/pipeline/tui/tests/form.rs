use super::*;
use crate::pipeline::args::PipelineArgs;
use crate::pipeline::tui::run_cmd as run;
use clap::FromArgMatches;

#[test]
fn every_flag_but_the_owned_ones_is_a_row() {
    let cmd = run();
    let form = Form::new(&cmd, OWN);
    let longs: Vec<&str> = form.fields.iter().map(|f| f.long.as_str()).collect();
    for a in cmd.get_arguments() {
        let Some(l) = a.get_long() else { continue };
        let countable = !matches!(
            a.get_action(),
            clap::ArgAction::Help | clap::ArgAction::Version | clap::ArgAction::Count
        );
        let owned = OWN.contains(&l)
            || a.get_help_heading()
                .is_some_and(|h| OWNED_HEADINGS.contains(&h));
        assert_eq!(longs.contains(&l), countable && !owned, "{l}");
    }
    assert_eq!(form.headings().first().map(String::as_str), Some("Common"));
}

#[test]
fn changed_values_round_trip_through_clap() {
    let cmd = run();
    let mut form = Form::new(&cmd, OWN);
    // Every switch flipped and every whole-number text bumped by one: the
    // parsed arguments differ from the defaults in exactly those fields.
    for f in &mut form.fields {
        match f.kind {
            Kind::Flag { .. } => f.toggle(),
            Kind::Text => {
                if let Ok(n) = f.default.parse::<u64>() {
                    f.set_text(&(n + 1).to_string());
                }
            }
            Kind::Choice(_) => {}
        }
    }
    let mut argv = vec!["run".to_string()];
    argv.extend(form.argv());
    check(&cmd, &argv).unwrap();
    let args = |argv: &[String]| {
        let m = cmd.clone().try_get_matches_from(argv).unwrap();
        serde_json::to_value(PipelineArgs::from_arg_matches(&m).unwrap()).unwrap()
    };
    let parsed = args(&argv);
    let base = args(&["run".to_string()]);
    // `max_threads` stands for the numbers, `zip` (`--no-zip`) for the switches.
    assert_eq!(
        parsed["max_threads"],
        base["max_threads"].as_u64().unwrap() + 1
    );
    assert_eq!(parsed["zip"], serde_json::json!(false));
    // Nothing else moved: the inputs stay unset.
    assert_eq!(parsed["gff_file"], base["gff_file"]);
}

#[test]
fn a_rejected_value_is_blamed_on_its_row() {
    let cmd = run();
    let mut form = Form::new(&cmd, OWN);
    let i = form
        .fields
        .iter()
        .position(|f| f.kind == Kind::Text && f.default.parse::<u64>().is_ok())
        .unwrap();
    form.fields[i].set_text("not-a-number");
    let mut argv = vec!["run".to_string()];
    argv.extend(form.argv());
    let complaint = check(&cmd, &argv).unwrap_err();
    assert_eq!(
        blamed(&complaint, &form.fields),
        Some(form.fields[i].long.as_str())
    );
}

#[test]
fn prefill_takes_what_the_command_line_gave() {
    let cmd = run();
    let mut form = Form::new(&cmd, OWN);
    let m = cmd
        .clone()
        .try_get_matches_from(["run", "--max-threads", "3", "--no-zip"])
        .unwrap();
    form.prefill(&m);
    assert_eq!(form.get("max-threads").unwrap().value, "3");
    assert!(form.argv().contains(&"--no-zip".to_string()));
    assert_eq!(form.changed(), 2);
}

#[test]
fn headings_follow_the_pipeline_order() {
    let cmd = run();
    let form = Form::new(&cmd, OWN);
    let headings = form.headings();
    let steps = ["Common", "SNP", "count", "ATOI", "m6A", "APA"];
    assert_eq!(headings[..steps.len()], steps, "{headings:?}");
    assert!(
        !headings.iter().any(|h| h.is_empty()),
        "every flag has a heading"
    );
    // Arguments after a flattened group must not inherit its heading.
    for a in cmd.get_arguments() {
        let Some(l) = a.get_long() else { continue };
        if l.starts_with("skip-") || l == "known-snps" || l == "depth-resolution-kb" {
            assert_ne!(a.get_help_heading(), Some("m6A"), "--{l}");
        }
    }
    for l in ["gff", "genome", "output", "control-bam", "known-snps"] {
        assert_eq!(
            cmd.get_arguments()
                .find(|a| a.get_long() == Some(l))
                .and_then(|a| a.get_help_heading()),
            Some("Inputs"),
            "--{l}"
        );
    }
}

#[test]
fn every_step_that_dims_names_a_heading() {
    use crate::pipeline::tui::steps::Step;
    let cmd = run();
    let headings: Vec<&str> = cmd
        .get_arguments()
        .filter_map(|a| a.get_help_heading())
        .collect();
    for s in Step::ALL {
        // Depth's one flag, its resolution, is on the Steps screen.
        if s != Step::Depth {
            assert!(headings.contains(&s.label()), "{}", s.label());
        }
    }
}
