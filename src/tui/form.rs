//! A command's flags as rows of a form, read from its clap definition, and
//! the command line the filled form stands for.
//!
//! Nothing here names a flag: every row comes from the command's `Arg`s, so
//! a flag added later shows up on its own.

use crate::script::{emit, many_of, Many};
use crate::tui::marker;
use clap::{Arg, ArgAction, ArgMatches, Command};
use data_beans::interactive::ui::{ACCENTED, DIM};
use ratatui::style::Style;
use ratatui::text::{Line, Span};

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

/// One flag of the command.
#[derive(Clone, Debug)]
pub struct Field {
    /// clap's id for the argument.
    pub id: String,
    /// The long name, without `--`.
    pub long: String,
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
    fn from_arg(a: &Arg) -> Option<Self> {
        let long = a.get_long()?.to_string();
        let heading = a.get_help_heading().unwrap_or("");
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
            long_help: a
                .get_long_help()
                .or_else(|| a.get_help())
                .map(ToString::to_string)
                .unwrap_or_default(),
            heading: heading.to_string(),
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
        self.line()
    }

    /// The flag with the value the row holds, default or not: a switch
    /// only when on its given side.
    #[must_use]
    pub fn line(&self) -> Vec<String> {
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

/// Row `f` as a line of the form: marked when `at`, starred and its
/// default beside it when changed, its name `width` wide in `style` with
/// `shown` (the value as it reads), and `complaint` after it when one blames
/// it.
#[must_use]
pub fn row_line(
    f: &Field,
    at: bool,
    width: usize,
    style: Style,
    shown: String,
    complaint: Option<&str>,
) -> Line<'static> {
    let changed = !f.is_default();
    let mut spans = vec![
        marker(at),
        Span::styled(if changed { "* " } else { "  " }, ACCENTED),
        Span::styled(format!("{:<width$}", f.long), style),
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
    if let Some(c) = complaint {
        spans.push(Span::styled(format!("  {c}"), ACCENTED));
    }
    Line::from(spans)
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

/// The flags of a command that a form shows, as rows.
#[derive(Clone, Debug)]
pub struct Form {
    pub fields: Vec<Field>,
}

impl Form {
    /// The form of `cmd`'s flags that `keep` keeps; `cmd` must be built.
    pub fn new(cmd: &Command, keep: impl Fn(&Arg) -> bool) -> Self {
        let mut form = Form {
            fields: cmd
                .get_arguments()
                .filter(|a| keep(a))
                .filter_map(Field::from_arg)
                .collect(),
        };
        // Stable: clap's order inside a heading, each heading contiguous.
        let order = form.headings();
        form.fields
            .sort_by_key(|f| order.iter().position(|h| *h == f.heading));
        form
    }

    /// The form of `cmd`'s flags named in `longs`, in that order.
    pub fn only(cmd: &Command, longs: &[&str]) -> Self {
        let mut form = Self::new(cmd, |a| a.get_long().is_some_and(|l| longs.contains(&l)));
        form.fields
            .sort_by_key(|f| longs.iter().position(|l| *l == f.long));
        form
    }

    pub fn get(&self, long: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.long == long)
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

    /// Take the values the command line (or the environment) gave.
    pub fn prefill(&mut self, m: &clap::ArgMatches) {
        for f in &mut self.fields {
            let Some(raw) = crate::script::given_words(m, &f.id) else {
                continue;
            };
            f.value = match f.kind {
                // A switch given on the line is the non-default side.
                Kind::Flag { .. } => if f.default == "true" { "false" } else { "true" }.into(),
                _ => raw.join(","),
            };
        }
    }
}

/// Whether clap takes `argv` (program name first); its complaint if not.
pub fn check(cmd: &Command, argv: &[String]) -> Result<(), String> {
    parse(cmd, argv).map(drop)
}

/// What clap makes of `argv` (program name first); its complaint if it
/// refuses it.
pub fn parse(cmd: &Command, argv: &[String]) -> Result<ArgMatches, String> {
    cmd.clone()
        .try_get_matches_from(argv)
        .map_err(|e| complaint(&e.render().to_string()))
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
