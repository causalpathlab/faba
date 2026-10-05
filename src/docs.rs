//! `faba docs` — the method write-ups, compiled into the binary.
//!
//! `include_str!`, not paths read at runtime. The binary is often the only thing on the machine
//! that ran the analysis (installed with `cargo install`, or copied to a cluster with no checkout
//! beside it), and a doc you cannot reach from there is a doc nobody reads. It also means the
//! build breaks if the file is moved or deleted — which enforces that it *exists*, though not
//! that it is *current*.
//!
//! A write-up is split at its `##` and `###` headings, so `faba docs qc` prints the section on
//! `qc` rather than the whole text. The annotation and lineage write-ups live in `senna docs`,
//! with the subcommands they describe.

use anyhow::Result;
use clap::Args;

/// Every write-up, in one place: the topic, a one-line blurb, and the text.
///
/// The listing `faba docs` prints and the text `faba docs <TOPIC>` prints are both read from
/// here, so the index can never advertise a topic the command cannot serve — which is exactly
/// what happens when the two are maintained separately.
const DOCS: &[(&str, &str, &str)] = &[(
    "profiling",
    "METHOD  BAM to per-cell features: m6A, A-to-I, APA, gene counts, SNPs",
    include_str!("../docs/profiling-methods.md"),
)];

#[derive(Args, Debug)]
pub struct DocsArgs {
    #[arg(
        value_name = "TOPIC|SECTION",
        help = "A write-up, or sections of one by keyword (omit to list what there is)",
        long_help = "A whole write-up by its topic, or sections of one by keyword:\n\
                     a subcommand named in a heading (`qc`, `pileup`, `run`),\n\
                     a section number (`8`, `1.3`), or any word of a heading (`records`).\n\
                     Several keywords print several sections, in document order.\n\
                     Omit to list the topics and their sections."
    )]
    pub query: Vec<String>,
}

/// One `##` or `###` section of a write-up.
struct Section<'a> {
    /// `2` for `##`, `3` for `###`.
    level: usize,
    /// The heading line without its `#`s.
    heading: &'a str,
    /// The section number, when the heading starts with one (`8`, `1.3`).
    number: Option<&'a str>,
    /// Byte range of the section in the text, heading included, up to the next heading of the
    /// same or a higher level.
    range: std::ops::Range<usize>,
}

impl Section<'_> {
    /// The subcommands the heading names in backticks.
    fn commands(&self) -> impl Iterator<Item = &str> {
        self.heading.split('`').skip(1).step_by(2)
    }
}

fn sections(text: &str) -> Vec<Section<'_>> {
    let mut starts: Vec<(usize, usize, &str)> = Vec::new();
    let mut at = 0;
    for line in text.split_inclusive('\n') {
        let level = line.bytes().take_while(|&b| b == b'#').count();
        if (level == 2 || level == 3) && line.as_bytes().get(level) == Some(&b' ') {
            starts.push((at, level, line[level..].trim()));
        }
        at += line.len();
    }
    starts
        .iter()
        .enumerate()
        .map(|(i, &(start, level, heading))| {
            let end = starts[i + 1..]
                .iter()
                .find(|s| s.1 <= level)
                .map_or(text.len(), |s| s.0);
            let number = heading
                .split_whitespace()
                .next()
                .map(|w| w.trim_end_matches('.'))
                .filter(|w| !w.is_empty() && w.chars().all(|c| c.is_ascii_digit() || c == '.'));
            Section {
                level,
                heading,
                number,
                range: start..end,
            }
        })
        .collect()
}

/// The sections `word` selects: a subcommand named in a `##` heading, else a section
/// number, else every heading containing the word. A section inside one already chosen is
/// left out.
fn select<'a>(all: &'a [Section<'a>], word: &str) -> Vec<&'a Section<'a>> {
    let w = word.to_lowercase();
    let by_command: Vec<_> = all
        .iter()
        .filter(|s| s.level == 2 && s.commands().any(|c| c.eq_ignore_ascii_case(&w)))
        .collect();
    if !by_command.is_empty() {
        return by_command;
    }
    if let Some(s) = all.iter().find(|s| s.number == Some(w.as_str())) {
        return vec![s];
    }
    outermost(
        all.iter()
            .filter(|s| s.heading.to_lowercase().contains(&w))
            .collect(),
    )
}

/// `chosen` in document order, without repeats or sections inside another one chosen.
fn outermost<'a>(mut chosen: Vec<&'a Section<'a>>) -> Vec<&'a Section<'a>> {
    chosen.sort_by_key(|s| (s.range.start, std::cmp::Reverse(s.range.end)));
    let mut out: Vec<&Section> = Vec::new();
    for s in chosen {
        if out.last().is_none_or(|o| s.range.start >= o.range.end) {
            out.push(s);
        }
    }
    out
}

pub fn run_docs(args: &DocsArgs) -> Result<()> {
    if args.query.is_empty() {
        println!("faba method write-ups: `faba docs <TOPIC>` prints one whole,");
        println!("`faba docs <KEYWORD>...` the sections whose heading names it.\n");
        for (topic, blurb, text) in DOCS {
            println!("  {topic:<14} {blurb}");
            for s in sections(text) {
                let indent = if s.level == 2 { "    " } else { "      " };
                println!("{indent}{}", s.heading);
            }
            println!();
        }
        println!("Annotation and lineage methods moved with their subcommands: `senna docs`.\n");
        return Ok(());
    }
    for word in &args.query {
        if let Some((_, _, text)) = DOCS.iter().find(|(t, _, _)| t.eq_ignore_ascii_case(word)) {
            println!("{text}");
        }
    }
    let words: Vec<&String> = args
        .query
        .iter()
        .filter(|w| !DOCS.iter().any(|(t, _, _)| t.eq_ignore_ascii_case(w)))
        .collect();
    if words.is_empty() {
        return Ok(());
    }
    for (_, _, text) in DOCS {
        let all = sections(text);
        let mut chosen: Vec<&Section> = Vec::new();
        for w in &words {
            let hits = select(&all, w);
            anyhow::ensure!(
                !hits.is_empty(),
                "no section of the write-ups matches `{w}`; run `faba docs` for the list"
            );
            chosen.extend(hits);
        }
        for (i, s) in outermost(chosen).iter().enumerate() {
            if i > 0 {
                println!();
            }
            print!(
                "{}",
                text[s.range.clone()]
                    .trim_end_matches(['\n', '-'])
                    .trim_end()
            );
            println!();
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
