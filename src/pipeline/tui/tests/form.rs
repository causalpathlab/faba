use super::*;
use crate::pipeline::args::PipelineArgs;
use clap::{CommandFactory, Parser};

#[derive(Parser)]
#[command(name = "run")]
struct Run {
    #[command(flatten)]
    args: PipelineArgs,
}

fn run() -> clap::Command {
    let mut c = Run::command();
    c.build();
    c
}

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
        assert_eq!(longs.contains(&l), countable && !OWN.contains(&l), "{l}");
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
    let parsed = serde_json::to_value(Run::try_parse_from(&argv).unwrap().args).unwrap();
    let base = serde_json::to_value(Run::try_parse_from(["run"]).unwrap().args).unwrap();
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
