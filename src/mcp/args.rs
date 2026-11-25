//! A call's arguments, read and checked, and the flags they become on the
//! command line.

use crate::mcp::schema::{generated, hidden, properties, refused, wanted, Generated, Take, Ty};
use anyhow::{bail, ensure, Context, Result};
use legume_mcp::Call;
use serde_json::{json, Value};

/// A tool's arguments.
pub(super) struct Arguments(pub(super) Value);

impl Arguments {
    pub(super) fn str(&self, key: &str) -> Option<&str> {
        self.0.get(key).and_then(Value::as_str)
    }

    fn req(&self, key: &str) -> Result<String> {
        self.str(key)
            .map(str::to_string)
            .with_context(|| format!("`{key}` is required"))
    }

    /// A path, which goes on the command line where a leading `-` would
    /// read as an option.
    pub(super) fn path(&self, key: &str) -> Result<String> {
        let p = self.req(key)?;
        ensure!(!p.starts_with('-'), "`{key}` {p}: write it as ./{p}");
        Ok(p)
    }

    pub(super) fn path_opt(&self, key: &str) -> Result<Option<String>> {
        match self.str(key) {
            None => Ok(None),
            Some(_) => self.path(key).map(Some),
        }
    }

    pub(super) fn parse<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.0.get(key) {
            None | Some(Value::Null) => Ok(None),
            Some(v) => serde_json::from_value(v.clone())
                .map(Some)
                .with_context(|| format!("reading `{key}`")),
        }
    }

    pub(super) fn flag(&self, key: &str) -> Result<bool> {
        Ok(self.parse(key)?.unwrap_or(false))
    }

    pub(super) fn strings(&self, key: &str) -> Result<Vec<String>> {
        Ok(self.parse(key)?.unwrap_or_default())
    }

    pub(super) fn paths(&self, key: &str) -> Result<Vec<String>> {
        let v = self.strings(key)?;
        if let Some(p) = v.iter().find(|p| p.starts_with('-')) {
            bail!("`{key}` {p}: write it as ./{p}");
        }
        Ok(v)
    }

    pub(super) fn paths_req(&self, key: &str) -> Result<Vec<String>> {
        let v = self.paths(key)?;
        ensure!(!v.is_empty(), "`{key}` is required");
        Ok(v)
    }

    /// `out`, or, when the agent names none, `{name}{ext}` in today's
    /// results folder, `-2`, `-3`… added to `name` (before `ext`) while
    /// `taken` says a name is in use.
    pub(super) fn out_or(
        &self,
        cx: &Call,
        (name, ext): (&str, &str),
        taken: impl Fn(&str) -> bool,
    ) -> Result<String> {
        if self.0.get("out").is_some_and(|v| !v.is_null()) {
            return self.path("out");
        }
        let dir = cx
            .results_dir()
            .context("`out` is required: this server keeps no results folder")?;
        (1..10_000)
            .map(|k| match k {
                1 => dir.join(format!("{name}{ext}")),
                k => dir.join(format!("{name}-{k}{ext}")),
            })
            .map(|p| p.to_string_lossy().into_owned())
            .find(|p| !taken(p))
            .context("no free name in the results folder")
    }

    /// `options`, refusing every name of a flag of `tool`'s command it
    /// does not offer (see [`refused`]), `--`, and a short option hidden
    /// inside a cluster (`-rX`).
    pub(super) fn options(&self, tool: &str) -> Result<Vec<String>> {
        let v = self.strings("options")?;
        for o in &v {
            let name = o.split('=').next().unwrap_or(o);
            let cluster = o.starts_with('-')
                && !o.starts_with("--")
                && o.len() > 2
                && o.chars().nth(1).is_some_and(|c| c.is_ascii_alphabetic());
            ensure!(
                !cluster,
                "{o}: give a short option apart from its value (`-g X`)"
            );
            ensure!(
                name != "--",
                "{o} is not taken in `options`: the tool names the command's inputs itself"
            );
            if let Some(why) = refused(tool, name) {
                bail!("{o} is not taken in `options`: {why}");
            }
        }
        Ok(v)
    }
}

/// Refuses a key of `args` that `tool` does not take (a null one is not
/// given): a flag of the MCMC fits faba no longer has, as `faba cluster`
/// refuses it, but for the values asking for the model left; a flag the
/// tool does not offer, saying why; any other, with the nearest key it
/// takes.
pub(super) fn known_keys(tool: &str, args: &Value) -> Result<()> {
    let Some(props) = properties(tool) else {
        return Ok(());
    };
    let (mut removed, mut unknown) = (Vec::new(), Vec::new());
    for (k, v) in (args.as_object().into_iter().flatten())
        .filter(|(k, v)| !v.is_null() && !props.contains_key(k.as_str()))
    {
        let flag = format!("--{}", k.replace('_', "-"));
        if hidden(tool, k) && crate::cluster::is_removed(&flag) {
            let value = v.as_str().map_or_else(|| v.to_string(), str::to_string);
            if crate::cluster::refused(&flag, &value) {
                removed.push(k.as_str());
            }
            continue;
        }
        unknown.push(match (refused(tool, &flag), nearest(k, props.keys())) {
            (Some(why), _) => format!("`{k}` is not an argument of {tool}: {why}"),
            (None, Some(near)) => {
                format!("`{k}` is not an argument of {tool}: did you mean `{near}`?")
            }
            (None, None) => format!("`{k}` is not an argument of {tool}"),
        });
    }
    removed.sort_unstable();
    ensure!(
        removed.is_empty(),
        "`{}` {}",
        removed.join("`, `"),
        crate::cluster::REMOVED_SAID
    );
    ensure!(unknown.is_empty(), "{}", unknown.join("; "));
    Ok(())
}

/// The key among `keys` nearest `k` by edits, if near enough to be a
/// slip: at most half `k`'s length away.
fn nearest<'k>(k: &str, keys: impl Iterator<Item = &'k String>) -> Option<&'k str> {
    (keys.map(|c| (edits(k, c), c)))
        .filter(|(d, _)| *d <= k.len().max(3) / 2)
        .min_by_key(|(d, _)| *d)
        .map(|(_, c)| c.as_str())
}

/// The fewest one-character insertions, deletions and substitutions
/// that make `a` into `b`.
fn edits(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let next = (diagonal + usize::from(ca != *cb))
                .min(row[j] + 1)
                .min(row[j + 1] + 1);
            diagonal = row[j + 1];
            row[j + 1] = next;
        }
    }
    row[b.len()]
}

impl Generated {
    /// One value given as `x`, as a word of the command line.
    fn word(&self, ty: &Ty, x: &Value) -> Result<String> {
        let key = &self.key;
        match (ty, x) {
            (Ty::Integer { .. }, Value::Number(n)) if n.is_i64() || n.is_u64() => Ok(n.to_string()),
            // 2.0 from a client that has only one kind of number.
            (Ty::Integer { .. }, Value::Number(n))
                if n.as_f64()
                    .is_some_and(|f| f.fract() == 0.0 && f.abs() < 9e15) =>
            {
                Ok(format!("{}", n.as_f64().unwrap_or_default() as i64))
            }
            (Ty::Number, Value::Number(n)) => Ok(n.to_string()),
            (Ty::Text, Value::String(s)) => {
                ensure!(
                    !(self.path && s.starts_with('-')),
                    "`{key}` {s}: write it as ./{s}"
                );
                Ok(s.clone())
            }
            (Ty::Choice(of), Value::String(s)) if of.contains(s) => Ok(s.clone()),
            (Ty::Choice(of), x) => bail!("`{key}` is one of {}, not {x}", of.join(", ")),
            (ty, x) => bail!(
                "`{key}` is {}, not {x}",
                wanted(&Self::property(
                    &Take::Value(ty.clone()),
                    legume_tui::form::Many::One,
                    ""
                ))
            ),
        }
    }

    /// The words the flag is given as for the argument's value `v`.
    fn words(&self, v: &Value) -> Result<Vec<String>> {
        use legume_tui::form::{emit, Many};
        let key = &self.key;
        match (&self.take, v) {
            (Take::Switch, Value::Bool(on)) => Ok(match on {
                true => vec![self.flag.clone()],
                false => vec![],
            }),
            (Take::Count, Value::Number(n)) if n.as_u64().is_some() => Ok(vec![
                self.flag.clone();
                n.as_u64().unwrap_or_default()
                    as usize
            ]),
            (Take::Value(ty), Value::Array(items)) if self.many != Many::One => {
                let vs = (items.iter())
                    .map(|x| self.word(ty, x))
                    .collect::<Result<Vec<_>>>()?;
                if let Many::Joined(c) = self.many {
                    if let Some(x) = vs.iter().find(|x| x.contains(c)) {
                        bail!("`{key}` {x}: give each value apart, none with `{c}`");
                    }
                }
                if let (Many::Listed, Some(x)) = (self.many, vs.iter().find(|x| x.starts_with('-')))
                {
                    bail!("`{key}` {x}: a value starting with `-` reads as a flag here");
                }
                Ok(emit(&self.flag, self.many, vs))
            }
            (Take::Value(ty), x) if self.many == Many::One => {
                Ok(emit(&self.flag, Many::One, vec![self.word(ty, x)?]))
            }
            (take, v) => bail!(
                "`{key}` is {}, not {v}",
                wanted(&Self::property(take, self.many, ""))
            ),
        }
    }
}

/// `tool`'s arguments made from flags, given in `a`, each with the words
/// it becomes, in the command's order; none that becomes no word.
fn typed_words(tool: &str, a: &Arguments) -> Result<Vec<(&'static Generated, Vec<String>)>> {
    let mut made = Vec::new();
    for g in generated(tool) {
        if let Some(v) = a.0.get(&g.key).filter(|v| !v.is_null()) {
            let words = g.words(v)?;
            if !words.is_empty() {
                made.push((g, words));
            }
        }
    }
    Ok(made)
}

/// `tool`'s arguments made from flags, given in `a`, as their flags, in
/// the command's order: the bare flag for a switch, the flag once per
/// value for a list (or the values joined, as the flag takes them).
#[cfg(test)]
pub(super) fn typed_flags(tool: &str, a: &Arguments) -> Result<Vec<String>> {
    Ok(typed_words(tool, a)?
        .into_iter()
        .flat_map(|(_, w)| w)
        .collect())
}

/// `args` with `tool`'s arguments made from flags put ahead of its
/// `options`, which still checks them; a flag given both ways is refused,
/// and so is one in `options` whose argument has another key than the
/// flag (`--fig-width` is `width`).
pub(super) fn with_typed(tool: &str, args: Value) -> Result<Value> {
    let a = Arguments(args);
    let typed = typed_words(tool, &a)?;
    let given = a.strings("options")?;
    let renamed = (generated(tool).iter())
        .filter(|g| g.flag.trim_start_matches("--").replace('-', "_") != g.key);
    for g in renamed {
        let named = |o: &&String| g.names.iter().any(|n| o.split('=').next() == Some(n));
        if let (true, Some(o)) = (
            a.0.get(&g.key).is_none_or(Value::is_null),
            given.iter().find(named),
        ) {
            bail!(
                "{o} is not taken in `options`: this tool has its own argument for it, `{}`",
                g.key
            );
        }
    }
    if typed.is_empty() {
        return Ok(a.0);
    }
    for (g, _) in &typed {
        let named = |o: &&String| g.names.iter().any(|n| o.split('=').next() == Some(n));
        if let Some(o) = given.iter().find(named) {
            bail!(
                "`{}` and `{o}` in `options` both give {}: give it once",
                g.key,
                g.flag
            );
        }
    }
    let typed: Vec<String> = typed.into_iter().flat_map(|(_, w)| w).collect();
    let Arguments(mut args) = a;
    if let Some(m) = args.as_object_mut() {
        m.insert("options".into(), json!([typed, given].concat()));
    }
    Ok(args)
}

pub(super) fn push_flag(argv: &mut Vec<String>, flag: &str, value: &str) {
    argv.extend([flag.to_string(), value.to_string()]);
}

pub(super) fn push_opt(argv: &mut Vec<String>, flag: &str, value: Option<&str>) {
    if let Some(value) = value {
        push_flag(argv, flag, value);
    }
}

/// One comma-separated argument. Clap splits that argument on commas, so a
/// value that already contains one would be read as several.
pub(super) fn comma(key: &str, paths: &[String]) -> Result<String> {
    no_comma(key, paths)?;
    Ok(paths.join(","))
}

pub(super) fn no_comma(key: &str, paths: &[String]) -> Result<()> {
    if let Some(p) = paths.iter().find(|p| p.contains(',')) {
        bail!("`{key}` {p}: a comma here is read as another item; rename it");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::mcp::schema::commands;
    use clap::FromArgMatches;

    use crate::mcp::schema::{generated, tools, value_type, Ty};

    use super::{known_keys, typed_flags, with_typed, Arguments};
    use serde_json::json;
    /// A value of the type `schema` says, as an agent would send it.
    fn sample(schema: &serde_json::Value) -> serde_json::Value {
        let one = |s: &serde_json::Value| match s["type"].as_str() {
            Some("boolean") => json!(true),
            Some("integer") => json!(2),
            Some("number") => json!(0.5),
            _ => s["enum"].get(0).cloned().unwrap_or(json!("x1")),
        };
        match schema["type"].as_str() {
            Some("array") => json!([one(&schema["items"]), one(&schema["items"])]),
            _ => one(schema),
        }
    }

    /// `line` with what clap says is missing from it, until nothing is: a
    /// flag required without another (`--gff` unless `--utr-bed`), or by
    /// one given (`--fig-width` needs `--plots`).
    fn completed(
        cli: &mut clap::Command,
        sub: &clap::Command,
        mut line: Vec<String>,
    ) -> Vec<String> {
        use clap::error::{ContextKind, ContextValue, ErrorKind};
        while let Err(e) = cli.try_get_matches_from_mut(&line) {
            if e.kind() != ErrorKind::MissingRequiredArgument {
                break;
            }
            let Some(ContextValue::Strings(missing)) = e.get(ContextKind::InvalidArg) else {
                break;
            };
            for m in missing.clone() {
                let flag = m.split(' ').next().unwrap_or(&m).to_string();
                let arg = (sub.get_arguments())
                    .find(|a| a.get_long().is_some_and(|l| flag == format!("--{l}")));
                line.push(flag);
                // A value of the flag's type: a number, or a path.
                match arg.filter(|a| a.get_action().takes_values()) {
                    Some(a) if value_type(a) != Ty::Text => line.push("1".into()),
                    Some(_) => line.push("p1".into()),
                    None => {}
                }
            }
        }
        line
    }

    /// Every argument made from a flag, given a value of its schema's
    /// type, becomes words `faba <command>` parses, the flag among what
    /// it was given on the command line.
    #[test]
    fn every_generated_argument_parses_as_its_flag() {
        use clap::parser::ValueSource;
        let specs = tools();
        let mut cli = crate::faba_command();
        let mut checked = 0;
        for (tool, sub) in commands() {
            let spec = specs.iter().find(|s| s.name == tool).unwrap();
            // What the command needs whatever else is given.
            let mut needed = vec!["faba".to_string(), sub.get_name().to_string()];
            for a in sub.get_arguments().filter(|a| a.is_required_set()) {
                if let Some(long) = a.get_long() {
                    needed.push(format!("--{long}"));
                }
                let numeric = a.get_default_values().is_empty() && value_type(a) != Ty::Text;
                needed.push(if numeric { "1" } else { "p1" }.into());
            }
            let needed = completed(&mut cli, &sub, needed);
            let id = |flag: &str| {
                (sub.get_arguments())
                    .find(|a| a.get_long().is_some_and(|l| flag == format!("--{l}")))
                    .map(|a| a.get_id().to_string())
                    .unwrap_or_default()
            };
            for g in generated(&tool) {
                let value = sample(&spec.input_schema["properties"][g.key.as_str()]);
                let a = Arguments(json!({ g.key.as_str(): value }));
                let words =
                    typed_flags(&tool, &a).unwrap_or_else(|e| panic!("{tool} {}: {e}", g.key));
                assert!(!words.is_empty(), "{tool} {}: no words for {value}", g.key);
                // `options`, where the words go, takes them.
                let passed = Arguments(json!({ "options": words })).options(&tool);
                assert!(passed.is_ok(), "{tool} {}: {passed:?}", g.key);
                // A required flag is given once: as the argument under test.
                let mut line = needed.clone();
                if let Some(i) = line.iter().position(|w| *w == g.flag) {
                    line.drain(i..(i + 2).min(line.len()));
                }
                line.extend(words.clone());
                let line = completed(&mut cli, &sub, line);
                let all = cli
                    .try_get_matches_from_mut(&line)
                    .unwrap_or_else(|e| panic!("{tool} {}: {line:?}: {e}", g.key));
                let (_, m) = all.subcommand().unwrap();
                assert_eq!(
                    m.value_source(&id(&g.flag)),
                    Some(ValueSource::CommandLine),
                    "{tool} {}: {line:?}",
                    g.key
                );
                assert!(
                    crate::Cli::from_arg_matches(&all).is_ok(),
                    "{tool} {}: {line:?}",
                    g.key
                );
                checked += 1;
            }
        }
        assert!(checked > 300, "only {checked} arguments made from flags");
    }

    /// A key a tool does not take is refused, naming it and the nearest
    /// key it does, or why its flag is not offered; a call of keys it
    /// takes, or of null ones, passes.
    #[test]
    fn a_key_the_tool_does_not_take_is_refused() {
        let refused = |tool, args| known_keys(tool, &args).unwrap_err().to_string();
        let call = json!({ "bam": ["a.bam"], "gff": "g.gff", "genome": "g.fa" });
        let mut slip = call.clone();
        slip["min_coverge"] = json!(3);
        let e = refused("atoi", slip);
        assert!(
            e.contains("`min_coverge`") && e.contains("`min_coverage`"),
            "{e}"
        );
        let mut output = call.clone();
        output["output"] = json!("o");
        let e = refused("atoi", output);
        assert!(e.contains("`output`") && e.contains("`out`"), "{e}");
        let e = refused("count", json!({ "zzzzzzzz": 1 }));
        assert!(
            e.contains("`zzzzzzzz`") && !e.contains("did you mean"),
            "{e}"
        );
        let mut given = call.clone();
        given["min_coverage"] = json!(3);
        given["out"] = json!("o");
        given["typo"] = json!(null);
        assert!(known_keys("atoi", &given).is_ok());
        // Every key a tool's schema offers is taken.
        for spec in tools() {
            let props = spec.input_schema["properties"].as_object().unwrap();
            let all: serde_json::Map<_, _> = props.keys().map(|k| (k.clone(), json!(1))).collect();
            assert!(known_keys(&spec.name, &all.into()).is_ok(), "{}", spec.name);
        }
    }

    /// `options` refuses every name of a flag the tool does not offer, by
    /// its long name, an alias or its short one: those excluded, and those
    /// the tool has its own argument for.
    #[test]
    fn options_refuse_every_name_of_a_flag_not_offered() {
        let options = |tool, o: &[&str]| Arguments(json!({ "options": o })).options(tool);
        for (tool, o) in [
            ("count", &["--output=z"][..]),
            ("count", &["-o", "z"]),
            ("count", &["-v"]),
            ("count", &["-h"]),
            ("count", &["--"]),
            ("count", &["-g", "g.gff"]),
            ("apa", &["-u", "u.bed"]),
            ("apa", &["--utr-bed", "u.bed"]),
            ("pileup", &["--gtf", "g.gtf"]),
            ("pileup", &["--bam", "a.bam"]),
            ("pileup", &["--svg"]),
            ("pileup", &["--png"]),
            ("pileup", &["--format", "pdf"]),
            ("pileup", &["--no-pdf"]),
            ("pileup", &["--batch-mode"]),
            ("metagene", &["--batch-process"]),
            ("metagene", &["--figure", "f"]),
            ("metagene", &["--dist-measures", "d"]),
            ("metagene", &["--print"]),
            ("cluster_view", &["--session", "s"]),
            ("cluster_view", &["--no-detach"]),
            ("cluster_view", &["--html", "p.html"]),
            ("qc", &["--html=p.html"]),
            ("qc", &["--plots"]),
            ("qc", &["--dry-run"]),
            ("qc_view", &["--html", "p.html"]),
        ] {
            let e = options(tool, o).unwrap_err().to_string();
            assert!(e.contains("not taken in `options`"), "{tool} {o:?}: {e}");
        }
        assert!(options("count", &["--expected-cells", "10"]).is_ok());
        assert!(options("pileup", &["--membership", "m.tsv"]).is_ok());
        // A flag offered under another key: in `options`, named by it.
        let refused = |tool, o: &[&str]| {
            with_typed(tool, json!({ "options": o }))
                .unwrap_err()
                .to_string()
        };
        let e = refused("pileup", &["--height", "5"]);
        assert!(e.contains("own argument for it, `ascii_height`"), "{e}");
        let e = refused("pwm", &["--fig-width=3"]);
        assert!(e.contains("own argument for it, `width`"), "{e}");
    }

    /// A flag given both as its argument and in `options`, under any of
    /// its names, is refused, naming it.
    #[test]
    fn a_flag_given_twice_is_refused() {
        let twice = |tool, args| with_typed(tool, args).unwrap_err().to_string();
        let e = twice(
            "count",
            json!({ "min_mapping_quality": 30, "options": ["--min-mapping-quality", "10"] }),
        );
        assert!(e.contains("--min-mapping-quality"), "{e}");
        let e = twice(
            "pwm",
            json!({ "threads": 2, "options": ["--max-threads=4"] }),
        );
        assert!(e.contains("--threads"), "{e}");
        let e = twice("qc", json!({ "width": 3, "options": ["--fig-width", "4"] }));
        assert!(e.contains("--fig-width"), "{e}");
        // Each once is taken.
        assert!(with_typed(
            "count",
            json!({ "min_mapping_quality": 30, "options": ["--expected-cells", "10"] })
        )
        .is_ok());
    }

    #[test]
    fn typed_arguments_become_flags_ahead_of_options() {
        let args = json!({
            "data": ["d"], "seed": 7, "within_min_cells": 20,
            "options": ["--gene-alpha", "1"],
        });
        let a = Arguments(with_typed("cluster", args).unwrap());
        let argv = a.options("cluster").unwrap();
        assert_eq!(
            argv,
            [
                "--seed",
                "7",
                "--within-min-cells",
                "20",
                "--gene-alpha",
                "1"
            ]
        );
        let flags = |tool, args| typed_flags(tool, &Arguments(args));
        let e = flags("metagene", json!({ "isoforms": "sometimes" })).unwrap_err();
        assert!(e.to_string().contains("one of weighted, all"), "{e}");
        assert!(flags("cluster", json!({ "seed": "7" })).is_err());
        assert!(flags("pileup", json!({ "cell_membership": "-m.tsv" })).is_err());
        assert!(flags("pileup", json!({ "depth": ["-d.zarr"] })).is_err());
        // The ASCII histogram's rows and the figure's height, apart.
        assert_eq!(
            flags("pileup", json!({ "ascii_height": 5, "height": 2.5 })).unwrap(),
            ["--height", "5", "--fig-height", "2.5"]
        );
        let pileup = tools().into_iter().find(|t| t.name == "pileup").unwrap();
        let props = &pileup.input_schema["properties"];
        assert_eq!(props["ascii_height"]["type"], "integer");
        assert_eq!(props["height"]["type"], "number");
        assert_ne!(props["ascii_height"], props["height"]);
        assert_eq!(
            flags("pileup", json!({ "track": ["wt=*_wt_*", "mut=*_mut_*"] })).unwrap(),
            ["--track", "wt=*_wt_*", "--track", "mut=*_mut_*"]
        );
        assert_eq!(
            flags(
                "run",
                json!({ "skip_apa": true, "skip_snp": false, "mixture": true })
            )
            .unwrap(),
            ["--mixture", "--skip-apa"]
        );
        // Not this tool's: left alone.
        assert!(flags("depth", json!({ "gene_type": "lncRNA" }))
            .unwrap()
            .is_empty());
    }

    /// A paper panel's look: the font, the panel alone and the palette, for
    /// every tool drawing a site figure, as its command's flags.
    #[test]
    fn figure_tools_take_the_panels_look() {
        let a = Arguments(json!({ "font_size": 6, "panel_only": true, "palette": "gray" }));
        for tool in ["metagene", "pwm", "qc"] {
            assert_eq!(
                typed_flags(tool, &a).unwrap(),
                ["--font-size", "6", "--panel-only", "--palette", "gray"],
                "{tool}"
            );
        }
        let e = typed_flags("qc", &Arguments(json!({ "palette": "blue" }))).unwrap_err();
        assert!(e.to_string().contains("one of colour, gray"), "{e}");
    }

    /// Axis titles: a text each on the logo and the profile, a `q=TEXT`
    /// each quantity on the site panels, each value whole however it reads
    /// (one led by `-` joined to its flag, lest it read as one).
    #[test]
    fn figure_tools_take_axis_titles() {
        let a = Arguments(json!({ "xlab": "-log10 p (wt vs. mut)", "ylab": "TEXT1, TEXT2" }));
        for tool in ["metagene", "pwm"] {
            assert_eq!(
                typed_flags(tool, &a).unwrap(),
                ["--xlab=-log10 p (wt vs. mut)", "--ylab", "TEXT1, TEXT2"],
                "{tool}"
            );
        }
        let a = Arguments(json!({
            "xlab": ["pv=-log10 p (wt vs. mut)", "fold=TEXT1, TEXT2"],
            "ylab": ["cells=cells with a converted read"],
        }));
        assert_eq!(
            typed_flags("site_plot", &a).unwrap(),
            [
                "--xlab",
                "pv=-log10 p (wt vs. mut)",
                "--xlab",
                "fold=TEXT1, TEXT2",
                "--ylab",
                "cells=cells with a converted read"
            ]
        );
    }
}
