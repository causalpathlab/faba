//! Which of the count fit's batches each of the sites' batches is. A batch
//! met in several folders is named with its folder, as typed in each round
//! (see [`args::qualify_batches`]), so the two rounds may name one batch
//! otherwise: `b1` in one, `runs/a/b1` in the other. Pairs named with
//! `--batch-map` come first; then a batch of the same name; then the one
//! batch that is the other named with its folder. What is still ambiguous
//! is asked on a terminal (see [`Question`]) and refused off one, naming
//! each choice and how to make it.
//!
//! [`args::qualify_batches`]: super::args::qualify_batches

use super::{qualifies, CountCells};
use crate::cluster::progress;
use crate::tui::ui::{header, Screen};
use crate::tui::View;
use legume_tui::widgets::{numbered, Menu, Outcome};
use ratatui::crossterm::event::KeyEvent;
use ratatui::Frame;
use rustc_hash::FxHashSet;

/// One `--batch-map` pair: a site batch and the count fit's batch it is,
/// none to leave it unmatched (`SITES=`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(in crate::cluster) struct Named {
    pub sites: Box<str>,
    pub count: Option<Box<str>>,
}

impl Named {
    /// `SITES=COUNT`, or `SITES=` for none, split at the last `=` (a site
    /// batch may have one in its name; the count fit's batches are named
    /// without).
    pub(in crate::cluster) fn parse(s: &str) -> anyhow::Result<Self> {
        let (sites, count) = s.rsplit_once('=').ok_or_else(|| {
            anyhow::anyhow!(
                "--batch-map {s}: write SITES=COUNT (the sites' batch, then the count fit's), \
                 or SITES= to leave it unmatched"
            )
        })?;
        anyhow::ensure!(!sites.is_empty(), "--batch-map {s}: name the sites' batch");
        Ok(Self {
            sites: sites.into(),
            count: (!count.is_empty()).then(|| count.into()),
        })
    }

    /// Refused when it names a batch the count fit `count` does not have:
    /// what can be checked before the sites are read.
    pub(in crate::cluster) fn check(&self, count: &CountCells) -> anyhow::Result<()> {
        let Some(c) = self.count.as_deref() else {
            return Ok(());
        };
        if count.cells.contains_key(c) {
            return Ok(());
        }
        let mut theirs: Vec<&str> = count.cells.keys().map(|k| &**k).collect();
        theirs.sort_unstable();
        anyhow::bail!(
            "--batch-map {}: the count fit has no batch {c} (its batches: {})",
            self.said(),
            theirs.join(", ")
        )
    }

    /// As `--batch-map` takes it.
    fn said(&self) -> String {
        format!("{}={}", self.sites, self.count.as_deref().unwrap_or(""))
    }
}

/// How a pair was made.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub(in crate::cluster) enum By {
    /// `--batch-map`, or chosen when asked.
    Named,
    /// The same name.
    Same,
    /// The one batch that is the other named with its folder.
    Folder,
    /// A count fit of one batch, sites of one: the same batch.
    Only,
}

/// A site batch and the count fit's batch it is.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub(in crate::cluster) struct Pair {
    pub sites: Box<str>,
    pub count: Box<str>,
    pub by: By,
}

/// Site batches and the count fit's that pair in more than one way: at
/// most one of `choices` (each a site batch and a count batch) is taken.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(in crate::cluster) struct Ambiguity {
    pub choices: Vec<(Box<str>, Box<str>)>,
}

impl Ambiguity {
    /// In words: `b` could be the count fit's `a/b` or `c/b`; `/x/b` and
    /// `/y/b` could each be the count fit's `b`; or the choices listed.
    pub(in crate::cluster) fn said(&self) -> String {
        let sites = distinct(self.choices.iter().map(|c| &*c.0));
        let count = distinct(self.choices.iter().map(|c| &*c.1));
        let quoted = |v: &[&str], and: &str| {
            let q: Vec<String> = v.iter().map(|s| format!("`{s}`")).collect();
            match q.split_last() {
                Some((last, rest)) if !rest.is_empty() => {
                    format!("{} {and} {last}", rest.join(", "))
                }
                _ => q.join(""),
            }
        };
        match (sites.as_slice(), count.as_slice()) {
            ([one], _) => format!(
                "the sites' batch `{one}` could be the count fit's {}",
                quoted(&count, "or")
            ),
            (_, [one]) => format!(
                "the sites' batches {} could each be the count fit's `{one}`",
                quoted(&sites, "and")
            ),
            _ => format!(
                "the sites' batches {} could be the count fit's {} as any of {}",
                quoted(&sites, "and"),
                quoted(&count, "and"),
                (self.choices.iter())
                    .map(|(s, c)| format!("{s}={c}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        }
    }
}

/// `v` without repeats, in order.
fn distinct<'a>(v: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let mut out: Vec<&str> = Vec::new();
    for s in v {
        if !out.contains(&s) {
            out.push(s);
        }
    }
    out
}

/// How the sites' batches pair with the count fit's.
#[derive(Debug, Default, PartialEq)]
pub(in crate::cluster) struct Pairing {
    pub pairs: Vec<Pair>,
    /// Site batches left unmatched though they could be a count batch:
    /// each, and why.
    pub left: Vec<(Box<str>, String)>,
    /// What is still to be chosen.
    pub ambiguous: Vec<Ambiguity>,
    /// What was chosen when asked, as `--batch-map` would name it.
    pub answered: Vec<String>,
}

impl Pairing {
    /// For the within record: each pair and how it was made, and what was
    /// left unmatched and why.
    pub(in crate::cluster) fn recorded(&self) -> serde_json::Value {
        let left: Vec<serde_json::Value> = (self.left.iter())
            .map(|(sites, why)| serde_json::json!({"sites": sites, "why": why}))
            .collect();
        serde_json::json!({"pairs": self.pairs, "unmatched": left})
    }

    /// Off a terminal: what is ambiguous, each choice, and how to make it.
    pub(in crate::cluster) fn refusal(&self) -> String {
        let mut lines =
            vec!["the sites' batches pair with the count fit's in more than one way:".to_string()];
        lines.extend(self.ambiguous.iter().map(|a| format!("  - {}", a.said())));
        let example = (self.ambiguous.first())
            .and_then(|a| a.choices.first())
            .map(|(s, c)| (s.to_string(), c.to_string()))
            .unwrap_or_default();
        lines.push(format!(
            "say which with --batch-map SITES=COUNT, once a pair (SITES= leaves a batch \
             unmatched), e.g. --batch-map '{}={}'; through `faba mcp`, the `cluster` tool's \
             `batch_map`, e.g. {{\"{}\": \"{}\"}}",
            example.0, example.1, example.0, example.1
        ));
        lines.join("\n")
    }
}

impl CountCells {
    /// How `ours` (the sites' batches, each once and in order) pair with
    /// the count fit's: `named` first (refused when it names a batch
    /// neither has), then the same
    /// name, then the one batch that is the other named with its folder
    /// (see [`qualifies`]); a count fit of one batch and sites of one are
    /// the same batch, whatever their names. A batch of the same name wins
    /// over one named with its folder; any other choice between pairings
    /// is left in [`Pairing::ambiguous`].
    pub(in crate::cluster) fn pairing(
        &self,
        ours: &[&str],
        named: &[Named],
    ) -> anyhow::Result<Pairing> {
        let mut theirs: Vec<&str> = self.cells.keys().map(|k| &**k).collect();
        theirs.sort_unstable();
        let mut out = Pairing::default();
        if let ([one], [only], []) = (ours, theirs.as_slice(), named) {
            out.pairs.push(pair(one, only, By::Only));
            return Ok(out);
        }
        let mut done: FxHashSet<&str> = FxHashSet::default();
        let mut taken: FxHashSet<&str> = FxHashSet::default();
        for n in named {
            anyhow::ensure!(
                ours.contains(&&*n.sites),
                "--batch-map {}: the sites have no batch {} (their batches: {})",
                n.said(),
                n.sites,
                ours.join(", ")
            );
            anyhow::ensure!(
                done.insert(&n.sites),
                "--batch-map names the sites' batch {} twice",
                n.sites
            );
            match &n.count {
                Some(c) => {
                    n.check(self)?;
                    anyhow::ensure!(
                        taken.insert(c),
                        "--batch-map names the count fit's batch {c} twice"
                    );
                    out.pairs.push(pair(&n.sites, c, By::Named));
                }
                None => (out.left).push((n.sites.clone(), "--batch-map left it unmatched".into())),
            }
        }
        let rest: Vec<&str> = ours.iter().copied().filter(|o| !done.contains(o)).collect();
        // The same name first: it wins over a folder's.
        for &o in &rest {
            if !self.cells.contains_key(o) {
                continue;
            }
            if taken.insert(o) {
                out.pairs.push(pair(o, o, By::Same));
            } else {
                out.left.push((
                    o.into(),
                    format!("the count fit's {o} was given to another batch by --batch-map"),
                ));
            }
            done.insert(o);
        }
        // Then a batch named with its folder, where only one pairing is
        // possible.
        let alike = |o: &str, t: &str| qualifies(o, t) || qualifies(t, o);
        let mut edges: Vec<Edge> = Vec::new();
        for &o in rest.iter().filter(|o| !done.contains(*o)) {
            let free: Vec<&str> = (theirs.iter().copied())
                .filter(|t| alike(o, t) && !taken.contains(t))
                .collect();
            if free.is_empty() {
                if let Some(t) = theirs.iter().find(|t| alike(o, t)) {
                    out.left.push((
                        o.into(),
                        format!("it could be the count fit's {t}, which another batch has"),
                    ));
                }
            }
            edges.extend(free.into_iter().map(|t| (o, t)));
        }
        // A pairing is sure when neither batch in it could pair otherwise.
        let alone = |&(o, t): &Edge| {
            edges.iter().filter(|e| e.0 == o).count() == 1
                && edges.iter().filter(|e| e.1 == t).count() == 1
        };
        let (sure, open): (Vec<Edge>, Vec<Edge>) = edges.iter().partition(|e| alone(e));
        out.pairs
            .extend(sure.iter().map(|(o, t)| pair(o, t, By::Folder)));
        out.ambiguous = components(&open);
        Ok(out)
    }
}

fn pair(sites: &str, count: &str, by: By) -> Pair {
    Pair {
        sites: sites.into(),
        count: count.into(),
        by,
    }
}

/// A site batch and a count batch it could be.
type Edge<'a> = (&'a str, &'a str);

/// `edges` in groups that share a batch, each an ambiguity, in order.
fn components(edges: &[Edge]) -> Vec<Ambiguity> {
    let mut groups: Vec<Vec<Edge>> = Vec::new();
    for &(o, t) in edges {
        let touching: Vec<usize> = (0..groups.len())
            .filter(|&g| groups[g].iter().any(|e| e.0 == o || e.1 == t))
            .collect();
        let mut merged = vec![(o, t)];
        for &g in touching.iter().rev() {
            merged.extend(groups.remove(g));
        }
        groups.push(merged);
    }
    (groups.into_iter())
        .map(|mut g| {
            g.sort_unstable();
            Ambiguity {
                choices: g.into_iter().map(|(o, t)| (o.into(), t.into())).collect(),
            }
        })
        .collect()
}

/// How the sites' batches `ours` pair with `count`'s, `named` taken first.
/// What is still ambiguous is asked, one at a time, when `ask` says how
/// (on a terminal: see [`ask`]): each answer is taken as if named. Off a
/// terminal it is refused, saying how to choose (see
/// [`Pairing::refusal`]); a question cancelled is [`progress::Cancelled`].
pub(in crate::cluster) fn resolve(
    count: &CountCells,
    ours: &[&str],
    named: &[Named],
    mut ask: Option<&mut Ask<'_>>,
) -> anyhow::Result<Pairing> {
    let given = named.len();
    // Each batch once (few), in order, before any question.
    let mut ours: Vec<&str> = ours.to_vec();
    ours.sort_unstable();
    ours.dedup();
    let mut named = named.to_vec();
    loop {
        let mut pairing = count.pairing(&ours, &named)?;
        let Some(first) = pairing.ambiguous.first() else {
            pairing.answered = named[given..].iter().map(Named::said).collect();
            return Ok(pairing);
        };
        let Some(ask) = ask.as_mut() else {
            anyhow::bail!(pairing.refusal());
        };
        match ask(first)? {
            None => return Err(progress::Cancelled.into()),
            Some(Answer::Take(i)) => {
                let (s, c) = &first.choices[i];
                named.push(Named {
                    sites: s.clone(),
                    count: Some(c.clone()),
                });
            }
            Some(Answer::Neither) => {
                for s in distinct(first.choices.iter().map(|c| &*c.0)) {
                    named.push(Named {
                        sites: s.into(),
                        count: None,
                    });
                }
            }
        }
    }
}

/// Asking which of an ambiguity's pairings is meant: the answer, none
/// when cancelled.
pub(in crate::cluster) type Ask<'a> = dyn FnMut(&Ambiguity) -> anyhow::Result<Option<Answer>> + 'a;

/// What was chosen when asked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::cluster) enum Answer {
    /// The choice at this index.
    Take(usize),
    /// None: the site batches in it left unmatched.
    Neither,
}

/// The question a terminal is asked: which of an ambiguity's pairings is
/// meant, or none, in a [`Menu`]; Esc cancels the run.
pub(in crate::cluster) struct Question {
    menu: Menu<Answer>,
    /// How many rows the menu has, each with its key: `1` on.
    rows: usize,
    /// The answer; `Some(None)` when cancelled.
    pub decision: Option<Option<Answer>>,
}

impl Question {
    pub(in crate::cluster) fn new(a: &Ambiguity) -> Self {
        let mut items: Vec<(String, String, Answer)> = (a.choices.iter().enumerate())
            .map(|(i, (s, c))| {
                let detail = format!("--batch-map '{s}={c}'");
                (
                    format!("sites' {s}  =  count fit's {c}"),
                    detail,
                    Answer::Take(i),
                )
            })
            .collect();
        items.push((
            "neither: leave them unmatched".into(),
            String::new(),
            Answer::Neither,
        ));
        Self {
            rows: items.len(),
            menu: Menu::new(a.said(), numbered(items)),
            decision: None,
        }
    }
}

impl View for Question {
    /// The rows' keys as there are rows: `1-3` for three; past nine,
    /// letters.
    fn help(&self) -> Option<crate::tui::KeyList> {
        let digits = match self.rows.min(9) {
            1 => "1".to_string(),
            n => format!("1-{n}"),
        };
        let what = match self.rows > 9 {
            true => "take that one (then letters)",
            false => "take that one",
        };
        Some(crate::tui::keys(&[
            ("↑/↓", "move"),
            ("PgUp/PgDn", "a page up/down"),
            ("Home/End", "first/last"),
            ("Enter", "take it"),
            (&digits, what),
            ("Esc/q", "cancel the run"),
            ("?", "these keys"),
        ]))
    }
}

impl Screen for Question {
    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        self.decision = Some(None);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        match self.menu.key(key.code) {
            Outcome::Open => {}
            Outcome::Cancel => self.decision = Some(None),
            Outcome::Chosen(answer) => self.decision = Some(Some(answer)),
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        frame.render_widget(header("cluster", "which batch is which?", "--within"), area);
        self.menu.draw(frame, area, &crate::tui::LOOK);
    }
}

/// Ask `a` on the terminal (see [`Question`]); none when cancelled.
pub(in crate::cluster) fn ask(a: &Ambiguity) -> anyhow::Result<Option<Answer>> {
    let mut q = Question::new(a);
    crate::tui::run_view(&mut q)?;
    Ok(q.decision.flatten())
}
