use super::*;

#[test]
fn defaults_and_skip_flags() {
    let mut s = Steps::default();
    assert!(s.argv(true).is_empty(), "all on but depth: nothing to pass");
    s.at = Step::ALL.iter().position(|x| *x == Step::Atoi).unwrap();
    s.toggle(true);
    assert_eq!(s.argv(true), ["--skip-atoi"]);
    assert!(!s.heading_on("ATOI", true) && s.heading_on("Common", true));
    assert!(s.heading_on("m6A", true) && !s.heading_on("m6A", false));
    assert!(s.heading_on("mixture", false) && s.heading_on("Steps", false));
    s.at = Step::ALL.iter().position(|x| *x == Step::M6a).unwrap();
    s.toggle(true);
    assert_eq!(s.argv(true), ["--skip-atoi", "--skip-m6a"]);
    assert_eq!(
        s.argv(false),
        ["--skip-atoi"],
        "no bg: the pipeline skips m6A itself"
    );
}

#[test]
fn m6a_needs_a_bg_bam() {
    let mut s = Steps::default();
    // Without bg, m6A is simply not run: the pipeline skips it, no flag needed.
    assert!(s.problems().is_empty() && s.argv(false).is_empty());
    s.at = Step::ALL.iter().position(|x| *x == Step::M6a).unwrap();
    s.toggle(false); // turns it off
    assert_eq!(s.toggle(false), Some("m6A needs a bg BAM"));
}

#[test]
fn depth_takes_its_resolution() {
    let mut s = Steps {
        at: Step::ALL.iter().position(|x| *x == Step::Depth).unwrap(),
        ..Steps::default()
    };
    s.toggle(true);
    assert!(s.problems().iter().any(|p| p.contains("resolution")));
    for bad in ["inf", "0", "-1", "NaN"] {
        s.depth_kb = bad.into();
        assert!(!s.problems().is_empty(), "{bad} is no resolution");
    }
    s.depth_kb = "50".into();
    assert!(s.problems().is_empty());
    assert_eq!(s.argv(true), ["--depth-resolution-kb", "50"]);
}
