//! The qc view's state as JSON (see [`crate::tui::live`]): what the user
//! arranged, read back by an agent, kept beside an exported figure, and
//! set again from a `--view-state` file; and the figure `s` saves, drawn
//! from it at a size and style asked for.
//!
//! The keys: `modality` on screen (`modalities` those there are),
//! `panel` (thresholds, genes, pairs) Tab is on, `focus` the knob (a
//! threshold's key, its flag without dashes: `site_max_pv`),
//! `thresholds` each knob's value (null when off; `flags` the same as
//! `faba qc` flags), `x_scale` each knob's histogram scale (log, sqrt,
//! linear) and `y_scale` the counts', `two_d` the 2D panels shown,
//! `pairs` the panels chosen per modality (each `[x, y]`, quantities as
//! `--pairs` names them; none, Auto's), `pair_at` the panel in focus,
//! `front` (wt, mut) and `p_axis` (p, neglog10) the arms' look, `titles`
//! the axis titles (`x` and `y` by quantity, `sites` the counts'),
//! `weight` (coverage, converted, sites) the gene and metagene bars',
//! `scales` each quantity's 2D scale (linear, log) and `bins` the
//! rulers' bins. Read only: `saves`, the panel `s` saves now, and `gene`,
//! the gene selected.

use super::*;
use crate::figure::size::Figure as Floor;
use crate::qc::pairs::Axis as Quantity;
use crate::qc::scale::Scales;
use crate::tui::live::Live;
use anyhow::Context;
use serde_json::{json, Map, Value};

/// A knob's key in a state: its flag without dashes (`site_max_pv`).
fn key(c: Criterion) -> String {
    c.flag().trim_start_matches("--").replace('-', "_")
}

/// The knob `k` names.
fn knob(k: &str) -> anyhow::Result<Criterion> {
    (SHOWN.into_iter()).find(|&c| key(c) == k).with_context(|| {
        let known: Vec<String> = SHOWN.into_iter().map(key).collect();
        format!("no knob {k}: one of {}", known.join(", "))
    })
}

fn panel_name(p: Panel) -> &'static str {
    match p {
        Panel::Thresholds => "thresholds",
        Panel::Genes => "genes",
        Panel::Pairs => "pairs",
    }
}

fn weight_name(w: Weight) -> &'static str {
    match w {
        Weight::Coverage => "coverage",
        Weight::Converted => "converted",
        Weight::Sites => "sites",
    }
}

/// A histogram scale by its name.
fn hist_scale(s: &str) -> anyhow::Result<Scale> {
    [Scale::Log, Scale::Sqrt, Scale::Linear]
        .into_iter()
        .find(|x| x.name() == s)
        .with_context(|| format!("no scale {s}: log, sqrt or linear"))
}

/// `v` as a string, or why not.
fn text<'v>(v: &'v Value, what: &str) -> anyhow::Result<&'v str> {
    v.as_str()
        .with_context(|| format!("{what} is {v}, not a string"))
}

/// `v` as an object, or why not.
fn object<'v>(v: &'v Value, what: &str) -> anyhow::Result<&'v Map<String, Value>> {
    v.as_object()
        .with_context(|| format!("{what} is {v}, not an object"))
}

/// An enum serde names, from `v`.
fn named<T: serde::de::DeserializeOwned>(v: &Value, what: &str) -> anyhow::Result<T> {
    serde_json::from_value(v.clone()).with_context(|| format!("{what} {v}"))
}

impl SitePicker<'_> {
    /// The axis titles as JSON: those asked for, by quantity.
    fn titles_json(&self) -> Value {
        let by = |y: bool| -> Map<String, Value> {
            (Quantity::ALL.into_iter())
                .filter_map(|q| Some((q.key().to_string(), json!(self.titles.given(q, y)?))))
                .collect()
        };
        json!({ "x": by(false), "y": by(true), "sites": self.titles.given_sites() })
    }

    /// The titles `v` asks for, in place of all of them.
    fn set_titles(&mut self, v: &Value) -> anyhow::Result<()> {
        let v = object(v, "titles")?;
        let mut titles = crate::qc::site_plot::Titles::default();
        for (axis, y) in [("x", false), ("y", true)] {
            let Some(given) = v.get(axis).filter(|g| !g.is_null()) else {
                continue;
            };
            for (q, t) in object(given, "a title list")? {
                let q = Quantity::parse(q)?;
                let t = (!t.is_null()).then(|| text(t, "a title")).transpose()?;
                titles.set(q, y, t.map(str::to_string));
            }
        }
        if let Some(t) = v.get("sites").filter(|t| !t.is_null()) {
            titles.set_sites(Some(text(t, "the sites' title")?.to_string()));
        }
        self.titles = titles;
        Ok(())
    }

    /// The panel `s` saves now: the focused pair in 2D, else the focused
    /// knob's histogram.
    fn saves(&self) -> Value {
        match self.titled() {
            (x, Some(y)) => json!({ "pair": [x.key(), y.key()] }),
            (q, None) => json!({ "hist": q.key() }),
        }
    }

    /// Every view's 2D axes read on `scales`, their rulers in about
    /// `bins` bins.
    fn rescale(&mut self, scales: Scales, bins: usize) {
        for v in &mut self.views {
            v.cache = pairs::Cache::new(scales, bins);
            v.generation += 1;
        }
    }
}

impl Live for SitePicker<'_> {
    fn view_name(&self) -> &'static str {
        match self.read_only() {
            true => "qc_view",
            false => "qc",
        }
    }

    fn view_state(&self) -> Value {
        let v = self.view();
        let thresholds: Map<String, Value> = (SHOWN.into_iter())
            .map(|c| {
                let at = (!c.is_off(&self.filter)).then(|| c.get(&self.filter));
                (key(c), json!(at))
            })
            .collect();
        let x_scale: Map<String, Value> = (SHOWN.into_iter())
            .map(|c| (key(c), json!(c.scales()[self.x_scale[c as usize]].name())))
            .collect();
        let pairs: Map<String, Value> = (self.views.iter())
            .filter(|v| !v.pairs.is_empty())
            .map(|v| {
                let p: Vec<[&str; 2]> = v.pairs.iter().map(|p| [p.0.key(), p.1.key()]).collect();
                (v.table.modality.to_string(), json!(p))
            })
            .collect();
        let scales: Map<String, Value> = (Quantity::ALL.into_iter())
            .map(|q| (q.key().to_string(), json!(v.cache.scales.of(q))))
            .collect();
        let gene =
            (self.list.selected()).and_then(|g| Some(v.genes.as_ref()?.symbol(g).to_string()));
        json!({
            "title": self.title,
            "modality": v.table.modality,
            "modalities": self.views.iter().map(|v| &v.table.modality).collect::<Vec<_>>(),
            "panel": panel_name(self.panel),
            "focus": key(self.criterion()),
            "thresholds": thresholds,
            "flags": qc_flags(&self.filter),
            "x_scale": x_scale,
            "y_scale": self.y_scale().name(),
            "two_d": self.two_d,
            "pairs": pairs,
            "pair_at": self.pair_at,
            "front": self.arm_look.front,
            "p_axis": self.arm_look.p_axis,
            "titles": self.titles_json(),
            "weight": weight_name(self.weight),
            "scales": scales,
            "bins": v.cache.bins,
            "saves": self.saves(),
            "gene": gene,
        })
    }

    fn set_view_state(&mut self, state: &Value) -> anyhow::Result<()> {
        let s = object(state, "the state")?;
        let get = |k: &str| s.get(k).filter(|v| !v.is_null());
        if let Some(m) = get("modality") {
            let m = text(m, "modality")?;
            let at = (self.views.iter()).position(|v| &*v.table.modality == m);
            let at = at.with_context(|| format!("no {m} sites in this view"))?;
            if at != self.modality {
                // The knob in focus stays, where the modality has it.
                let c = self.criterion();
                self.modality = at;
                self.refresh_genes(true);
                self.focus = (self.view().criteria.iter())
                    .position(|&x| x == c)
                    .unwrap_or(0);
            }
        }
        if let Some(t) = get("thresholds") {
            for (k, v) in object(t, "thresholds")? {
                let c = knob(k)?;
                let raw = match v {
                    Value::Null => c.permissive(),
                    v => v
                        .as_f64()
                        .with_context(|| format!("{k} is {v}, not a number"))?,
                };
                c.set(&mut self.filter, raw);
            }
            for c in SHOWN {
                for view in &mut self.views {
                    view.update(c, &self.filter);
                }
            }
        }
        if get("scales").is_some() || get("bins").is_some() {
            let v = self.view();
            let scales = match get("scales") {
                Some(given) => {
                    let items: Vec<String> = (object(given, "scales")?.iter())
                        .map(|(q, s)| Ok(format!("{q}={}", text(s, "a scale")?)))
                        .collect::<anyhow::Result<_>>()?;
                    Scales::parse(&items)?
                }
                None => v.cache.scales,
            };
            let bins = match get("bins") {
                Some(b) => b.as_u64().context("bins is a whole number")? as usize,
                None => v.cache.bins,
            };
            self.rescale(scales, bins.max(1));
        }
        if let Some(x) = get("x_scale") {
            for (k, v) in object(x, "x_scale")? {
                let c = knob(k)?;
                let want = hist_scale(text(v, "a scale")?)?;
                let at = (c.scales().iter()).position(|&s| s == want);
                self.x_scale[c as usize] =
                    at.with_context(|| format!("{k} has no {} scale", want.name()))?;
            }
        }
        if let Some(y) = get("y_scale") {
            use clap::ValueEnum;
            let y = text(y, "y_scale")?;
            let s = crate::qc::arms::CountScale::from_str(y, false)
                .map_err(|_| anyhow::anyhow!("no y_scale {y}: linear, sqrt or log1p"))?;
            let q = Quantity::of(self.criterion());
            self.arm_look.y.set(q, s);
        }
        if let Some(f) = get("focus") {
            let c = knob(text(f, "focus")?)?;
            let at = self.view().criteria.iter().position(|&x| x == c);
            self.focus =
                at.with_context(|| format!("{} has no {}", self.view().table.modality, key(c)))?;
        }
        if let Some(p) = get("pairs") {
            for (m, list) in object(p, "pairs")? {
                let at = (self.views.iter()).position(|v| &*v.table.modality == m);
                let at = at.with_context(|| format!("pairs: no {m} sites in this view"))?;
                let list: Vec<[String; 2]> = named(list, "pairs")?;
                let chosen = (list.iter())
                    .map(|[x, y]| {
                        let (x, y) = (Quantity::parse(x)?, Quantity::parse(y)?);
                        Ok((x, y, "chosen".to_string()))
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                self.views[at].pairs = chosen;
            }
        }
        if let Some(t) = get("two_d") {
            self.two_d = t.as_bool().context("two_d is true or false")?;
            self.ensure_pairs();
        }
        if let Some(k) = get("pair_at") {
            let k = k.as_u64().context("pair_at is a whole number")? as usize;
            self.pair_at = k.min(self.view().shown_pairs().len().saturating_sub(1));
        }
        if let Some(p) = get("panel") {
            self.panel = match text(p, "panel")? {
                "thresholds" => Panel::Thresholds,
                "genes" => Panel::Genes,
                "pairs" if self.two_d => Panel::Pairs,
                "pairs" => anyhow::bail!("panel pairs is shown only with two_d"),
                other => anyhow::bail!("no panel {other}: thresholds, genes or pairs"),
            };
        }
        if let Some(f) = get("front") {
            self.arm_look.front = named(f, "front")?;
        }
        if let Some(p) = get("p_axis") {
            self.arm_look.p_axis = named(p, "p_axis")?;
        }
        if let Some(t) = get("titles") {
            self.set_titles(t)?;
        }
        if let Some(w) = get("weight") {
            let w = text(w, "weight")?;
            self.weight = [Weight::Coverage, Weight::Converted, Weight::Sites]
                .into_iter()
                .find(|&x| weight_name(x) == w)
                .with_context(|| format!("no weight {w}: coverage, converted or sites"))?;
        }
        // Drawn again whole: the knob's column, the counts, every image.
        self.rebuild_column();
        self.pair_drawn = Default::default();
        self.meta_plot.invalidate();
        self.plot.invalidate();
        Ok(())
    }

    fn own_sheet(&self) -> (Sheet, Floor) {
        let own = self.viewing.as_ref().map_or(Sheet::OWN, |v| v.1);
        (own, Floor::Qc)
    }

    fn saved_figure(&self, sheet: Sheet) -> anyhow::Result<String> {
        match &self.viewing {
            Some((plan, _)) => self.screen_panel(plan, sheet),
            None => self.figure(sheet),
        }
    }
}

/// The figure `s` saves in `faba qc-view` over the site tables, with
/// `state` (a `--view-state` file's) set first, written under `prefix` as
/// `formats` ask: what a live view's export writes, without the view.
/// The paths written.
pub fn write_view_figure(
    input_dir: &str,
    (tables, site_cells): (
        &FxHashMap<Box<str>, SiteTable>,
        &crate::qc::sites::SiteCells,
    ),
    filter: &SiteFilterArgs,
    genes: (Option<&str>, &[Box<str>]),
    (plan, sheet): (crate::qc::site_plot::Plan, Sheet),
    state: Option<&Value>,
    (prefix, formats): (&str, &[crate::figure::Format]),
) -> anyhow::Result<Vec<String>> {
    let scales = (plan.scales, plan.bins_on(&sheet));
    let mut picker = picker(input_dir, tables, site_cells, filter, genes, Some(scales))?;
    picker.view_only(plan, sheet);
    if let Some(state) = state {
        set_from_file(&mut picker, state)?;
    }
    let svg = picker.saved_figure(sheet)?;
    crate::figure::save_as(&svg, prefix, formats)
}

/// `picker` arranged as a `--view-state` file's `state` says, refused
/// when the file is another view's.
pub(super) fn set_from_file(picker: &mut SitePicker, state: &Value) -> anyhow::Result<()> {
    crate::tui::live::check_view(state, picker.view_name())?;
    (picker.set_view_state(state)).context("--view-state")
}
