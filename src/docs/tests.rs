use super::*;

const TEXT: &str = "# Title\n\nintro\n\n\
## 1. Shared\n\nshared body\n\n\
### 1.1 The gene model\n\nmodel body\n\n\
### 1.2 Run records\n\nrecord body\n\n---\n\n\
## 2. `cmd1` and `cmd1-report` — first\n\nfirst body\n\n---\n\n\
## 3. `cmd2` — second, with a gene model\n\nsecond body\n";

fn picked<'a>(all: &'a [Section<'a>], word: &str) -> Vec<&'a str> {
    select(all, word).iter().map(|s| s.heading).collect()
}

#[test]
fn sections_nest_under_their_parent() {
    let all = sections(TEXT);
    let top = &all[0];
    assert_eq!(top.number, Some("1"));
    assert!(TEXT[top.range.clone()].contains("record body"));
    assert!(!TEXT[top.range.clone()].contains("first body"));
    let sub = &all[1];
    assert_eq!(sub.number, Some("1.1"));
    assert!(!TEXT[sub.range.clone()].contains("record body"));
}

#[test]
fn a_command_picks_its_section_only() {
    let all = sections(TEXT);
    assert_eq!(
        picked(&all, "cmd1"),
        ["2. `cmd1` and `cmd1-report` — first"]
    );
    assert_eq!(
        picked(&all, "CMD1-report"),
        ["2. `cmd1` and `cmd1-report` — first"]
    );
}

#[test]
fn a_number_or_a_word_picks_sections() {
    let all = sections(TEXT);
    assert_eq!(picked(&all, "1.2"), ["1.2 Run records"]);
    assert_eq!(picked(&all, "records"), ["1.2 Run records"]);
    // A word in two headings picks both, in document order.
    assert_eq!(
        picked(&all, "gene model"),
        [
            "1.1 The gene model",
            "3. `cmd2` — second, with a gene model"
        ]
    );
    // A subsection inside a picked section is not printed twice.
    assert_eq!(picked(&all, "1"), ["1. Shared"]);
    assert!(picked(&all, "nothing-here").is_empty());
}

#[test]
fn every_subcommand_has_a_section() {
    let (_, _, text) = DOCS[0];
    let all = sections(text);
    for cmd in [
        "dartseq", "atoi", "apa", "count", "snp", "depth", "pileup", "qc", "all",
    ] {
        assert!(!select(&all, cmd).is_empty(), "no section for {cmd}");
    }
}
