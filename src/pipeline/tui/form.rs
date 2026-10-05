//! The flags of `faba run` as rows of a form, read from its clap
//! definition, and the command line the filled form stands for.
//!
//! Nothing here names a flag of the pipeline apart from [`OWN`]: every row
//! comes from the command's `Arg`s, so a flag added later shows up on its own.

use clap::{Arg, ArgAction, Command};

/// Flags the Inputs and Steps screens own, or that are never passed.
pub const OWN: &[&str] = &[
    "gff",
    "genome",
    "output",
    "control-bam",
    "known-snps",
    "batch-process",
    "skip-snp",
    "skip-count",
    "skip-atoi",
    "skip-apa",
    "depth-resolution-kb",
    "help",
    "version",
    "verbose",
];

/// What a row holds and how it is changed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A switch, given or not. `on` is the value it sets when given.
    Flag { on: bool },
    /// One of a fixed set of values; `""` stands for unset.
    Choice(Vec<String>),
    /// Free text.
    Text,
}

/// How several values reach the command line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Many {
    One,
    /// Joined into one argument: `--flag a,b`.
    Joined(char),
    /// Listed after one flag: `--flag a b`.
    Listed,
    /// The flag repeated: `--flag a --flag b`.
    Repeated,
}

/// One flag of the command.
#[derive(Clone, Debug)]
pub struct Field {
    /// clap's id for the argument.
    pub id: String,
    /// The long name, without `--`.
    pub long: String,
    pub help: String,
    pub long_help: String,
    /// The help heading the flag is listed under; `""` for none.
    pub heading: String,
    pub kind: Kind,
    /// The value clap would use when the flag is not given; `""` for none.
    pub default: String,
    pub value: String,
    /// Hidden from `--help`: shown only with the advanced rows.
    pub advanced: bool,
    many: Many,
}

impl Field {
    fn from_arg(a: &Arg, own: &[&str]) -> Option<Self> {
        let long = a.get_long()?.to_string();
        if own.contains(&long.as_str()) {
            return None;
        }
        let kind = match a.get_action() {
            ArgAction::SetTrue => Kind::Flag { on: true },
            ArgAction::SetFalse => Kind::Flag { on: false },
            ArgAction::Set | ArgAction::Append => {
                let values: Vec<String> = a
                    .get_possible_values()
                    .iter()
                    .filter(|v| !v.is_hide_set())
                    .map(|v| v.get_name().to_string())
                    .collect();
                if values.is_empty() {
                    Kind::Text
                } else {
                    Kind::Choice(values)
                }
            }
            // Counters, help and version: not something to fill in.
            _ => return None,
        };
        let many = many_of(a);
        let sep = match many {
            Many::Joined(c) => c.to_string(),
            _ => ",".to_string(),
        };
        let default = a
            .get_default_values()
            .iter()
            .map(|v| v.to_string_lossy().into_owned())
            .collect::<Vec<_>>()
            .join(&sep);
        let mut kind = kind;
        if let Kind::Choice(values) = &mut kind {
            // A choice with no default can stay unset.
            if default.is_empty() {
                values.insert(0, String::new());
            }
        }
        Some(Field {
            id: a.get_id().to_string(),
            help: a.get_help().map(ToString::to_string).unwrap_or_default(),
            long_help: a
                .get_long_help()
                .or_else(|| a.get_help())
                .map(ToString::to_string)
                .unwrap_or_default(),
            heading: a.get_help_heading().unwrap_or("").to_string(),
            value: default.clone(),
            default,
            advanced: a.is_hide_set(),
            kind,
            many,
            long,
        })
    }

    /// Whether the value is the one clap would use anyway.
    #[must_use]
    pub fn is_default(&self) -> bool {
        self.value.trim() == self.default
    }

    /// Back to what clap would use.
    pub fn reset(&mut self) {
        self.value.clone_from(&self.default);
    }

    /// Typed text; empty means back to the default.
    pub fn set_text(&mut self, text: &str) {
        let t = text.trim();
        if t.is_empty() {
            self.reset();
        } else {
            self.value = t.to_string();
        }
    }

    /// Flip a switch.
    pub fn toggle(&mut self) {
        if let Kind::Flag { .. } = self.kind {
            self.value = if self.value == "true" {
                "false"
            } else {
                "true"
            }
            .into();
        }
    }

    /// Step a choice by `d`, wrapping around.
    pub fn cycle(&mut self, d: isize) {
        if let Kind::Choice(values) = &self.kind {
            let n = values.len() as isize;
            let at = values.iter().position(|v| *v == self.value).unwrap_or(0) as isize;
            self.value
                .clone_from(&values[(at + d).rem_euclid(n) as usize]);
        }
    }

    /// The value as shown: `on` / `off` for a switch, `(unset)` for nothing.
    #[must_use]
    pub fn shown(&self) -> String {
        match self.kind {
            Kind::Flag { .. } => if self.value == "true" { "on" } else { "off" }.into(),
            _ if self.value.trim().is_empty() => "(unset)".into(),
            _ => self.value.clone(),
        }
    }

    /// What this row adds to the command line: nothing when it is the
    /// default.
    #[must_use]
    pub fn argv(&self) -> Vec<String> {
        if self.is_default() {
            return Vec::new();
        }
        let flag = format!("--{}", self.long);
        match self.kind {
            Kind::Flag { on } => {
                if (self.value == "true") == on {
                    vec![flag]
                } else {
                    Vec::new()
                }
            }
            Kind::Choice(_) | Kind::Text => emit(&flag, self.many, values(self.many, &self.value)),
        }
    }
}

fn many_of(a: &Arg) -> Many {
    if let Some(c) = a.get_value_delimiter() {
        Many::Joined(c)
    } else if a.get_num_args().is_some_and(|r| r.max_values() > 1) {
        Many::Listed
    } else if matches!(a.get_action(), ArgAction::Append) {
        Many::Repeated
    } else {
        Many::One
    }
}

/// The values typed into a row: split on commas and spaces when the flag
/// takes several, else the whole text.
fn values(many: Many, text: &str) -> Vec<String> {
    let text = text.trim();
    if text.is_empty() {
        return Vec::new();
    }
    match many {
        Many::One => vec![text.to_string()],
        _ => text
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|v| !v.is_empty())
            .map(str::to_string)
            .collect(),
    }
}

/// `flag` with values `vs`, laid out the way the flag takes them.
fn emit(flag: &str, many: Many, vs: Vec<String>) -> Vec<String> {
    if vs.is_empty() {
        return Vec::new();
    }
    match many {
        Many::One => vec![flag.to_string(), vs[0].clone()],
        Many::Joined(c) => vec![flag.to_string(), vs.join(&c.to_string())],
        Many::Listed => std::iter::once(flag.to_string()).chain(vs).collect(),
        Many::Repeated => vs.into_iter().flat_map(|v| [flag.to_string(), v]).collect(),
    }
}

/// Every flag of `faba run` but `own`, as rows.
#[derive(Clone, Debug)]
pub struct Form {
    pub fields: Vec<Field>,
}

impl Form {
    /// The form for `run`, which must be built.
    pub fn new(run: &Command, own: &[&str]) -> Self {
        let mut form = Form {
            fields: run
                .get_arguments()
                .filter_map(|a| Field::from_arg(a, own))
                .collect(),
        };
        // Stable: clap's order inside a heading, each heading contiguous.
        let order = form.headings();
        form.fields
            .sort_by_key(|f| order.iter().position(|h| *h == f.heading));
        form
    }

    pub fn get(&self, long: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.long == long)
    }

    pub fn get_mut(&mut self, long: &str) -> Option<&mut Field> {
        self.fields.iter_mut().find(|f| f.long == long)
    }

    /// Flags changed from their defaults.
    #[must_use]
    pub fn changed(&self) -> usize {
        self.fields.iter().filter(|f| !f.is_default()).count()
    }

    /// Every changed flag, in form order.
    #[must_use]
    pub fn argv(&self) -> Vec<String> {
        self.fields.iter().flat_map(Field::argv).collect()
    }

    /// The headings in form order, `Common` first.
    #[must_use]
    pub fn headings(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for f in &self.fields {
            if !out.contains(&f.heading) {
                out.push(f.heading.clone());
            }
        }
        out.sort_by_key(|h| h != "Common");
        out
    }

    /// Take the values the command line gave.
    pub fn prefill(&mut self, m: &clap::ArgMatches) {
        for f in &mut self.fields {
            if m.value_source(&f.id) != Some(clap::parser::ValueSource::CommandLine) {
                continue;
            }
            match f.kind {
                // A switch given on the line is the non-default side.
                Kind::Flag { .. } => {
                    f.value = if f.default == "true" { "false" } else { "true" }.into();
                }
                _ => {
                    let raw: Vec<String> = m
                        .get_raw(&f.id)
                        .into_iter()
                        .flatten()
                        .map(|v| v.to_string_lossy().into_owned())
                        .collect();
                    f.value = raw.join(",");
                }
            }
        }
    }
}

/// Whether clap takes `argv` (program name first); its complaint if not.
pub fn check(run: &Command, argv: &[String]) -> Result<(), String> {
    match run.clone().try_get_matches_from(argv) {
        Ok(_) => Ok(()),
        Err(e) => Err(complaint(&e.render().to_string())),
    }
}

/// The first line of a clap error, without its `error: ` lead, with the
/// indented lines under it (such as the missing arguments).
fn complaint(rendered: &str) -> String {
    let mut lines = rendered.lines();
    let first = lines.next().unwrap_or_default().trim();
    let first = first.strip_prefix("error: ").unwrap_or(first);
    std::iter::once(first)
        .chain(
            lines
                .take_while(|l| l.starts_with(char::is_whitespace) && !l.trim().is_empty())
                .map(str::trim),
        )
        .collect::<Vec<_>>()
        .join(" ")
}

/// The flag a clap complaint is about, when it names one of `fields`.
#[must_use]
pub fn blamed<'a>(complaint: &str, fields: &'a [Field]) -> Option<&'a str> {
    fields
        .iter()
        .map(|f| f.long.as_str())
        .filter(|l| {
            complaint.contains(&format!("'--{l}")) || complaint.contains(&format!("--{l} "))
        })
        .max_by_key(|l| l.len())
}

#[cfg(test)]
#[path = "tests/form.rs"]
mod tests;
