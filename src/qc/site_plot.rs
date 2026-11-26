//! `faba site-plot`, and the panels `faba qc --plots` writes beside its own
//! figures: one modality's site quantities as paper panels, a histogram
//! each and chosen pairs two at a time, with marginals on the pair's own
//! axes. Every figure of a quantity reads it on one [`Ruler`] (its
//! [`Scale`], range and bins), so a histogram, a 2D axis and its marginal
//! line up tick for tick. Thresholds are optional: given ones are drawn as
//! dashed lines, the sites they remove apart from those kept.
//!
//! It reads the site table and the `_site` matrices (cells per site over
//! every cell: no cell QC), and writes figures only.

use super::args::SiteFilterArgs;
use super::arms::{ArmLook, CountScale, Front, KeyAt, PAxis, YScales};
use super::layout::{scan_input_dir, InputLayout, SITE_MODALITIES};
use super::matrix::{open_matrix, row_nnz_over_columns, select_columns};
use super::pairs::{self, Axis as Quantity, Cache, Sites};
use super::scale::{Ruler, Scale, Scales};
use super::site_tui::Figure;
use super::sites::{accumulate_site_cells, read_site_table, Criterion, SiteTable};
use crate::figure::style::{Sheet, StyleArgs};
use crate::figure::SizeArgs;
use clap::Args;
use data_beans::figure::panel::{DensityPanel, HistPanel, Look, Marginals};
use rustc_hash::FxHashMap;
use std::borrow::Cow;

/// A histogram's own size in points, when none is asked for.
pub const HIST_OWN: (f64, f64) = (252.0, 180.0);
/// A pair's own size in points, when none is asked for.
pub const PAIR_OWN: (f64, f64) = (252.0, 252.0);

/// Which panels, on which scales: `--hist`, `--pairs`, `--marginals`,
/// `--scale` and `--bins`, as `faba site-plot` and `faba qc --plots` take
/// them.
#[derive(Args, Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct PanelArgs {
    #[arg(
        long = "hist",
        value_name = "LIST",
        value_delimiter = ',',
        help = "Quantities to draw a histogram of: pv, fold, log_odds, coverage, converted, edit_ratio, cells (default: every one that applies; `none` for none)"
    )]
    pub hist: Vec<String>,

    #[arg(
        long = "pairs",
        value_name = "X:Y,...",
        value_delimiter = ',',
        help = "Pairs to draw two at a time, x:y, e.g. pv:cells,fold:pv (default: the pairs `faba qc` picks; `none` for none)"
    )]
    pub pairs: Vec<String>,

    #[arg(
        long = "marginals",
        value_name = "LIST",
        value_delimiter = ',',
        help = "Histograms beside each pair on its own axes: top (x), right (y), or both (default: none)"
    )]
    pub marginals: Vec<String>,

    #[arg(
        long = "style",
        value_enum,
        default_value_t = PanelStyle::Panel,
        help = "How the panels are drawn: panel (the paper panel, kept and removed stacked) or view (as the qc view draws its histogram and its 2D panel)"
    )]
    pub style: PanelStyle,

    #[arg(
        long = "marginal-size",
        value_name = "FRACTION",
        help = "With --style view, each marginal's depth as a fraction of the plot's shorter side, above 0 and below 1 (default: 0.2); the panel style's marginals keep their own fifth"
    )]
    pub marginal_size: Option<f64>,

    #[arg(
        long = "p-axis",
        value_enum,
        value_name = "p|neglog10",
        help = "With --style view, how the p-value's ticks read: p itself (1, 0.2, 1e-3) or -log10 p; the bars sit on -log10 p either way (default: p)"
    )]
    pub p_axis: Option<PAxis>,

    #[arg(
        long = "front",
        value_enum,
        value_name = "wt|mut",
        help = "With --style view, which arm is drawn in front: wt (the signal arm) or mut (the control arm), in histograms, squares and marginals alike (default: wt)"
    )]
    pub front: Option<Front>,

    #[arg(
        long = "key",
        default_value_t = false,
        help = "With --style view, keep each histogram's key (its series and their colours) under --panel-only, which drops it otherwise; and draw each 2D panel's key to its kept squares' shades (a bar from 1 site to the busiest square's count) in the corner between its marginals, or in the plot's top right corner"
    )]
    pub key: bool,

    #[arg(
        long = "key-at",
        value_enum,
        value_name = "right|top|inside|alone",
        help = "With --style view, where each histogram's key goes: right (beside the plot, which narrows by its width; a key that would leave the plot under 60% of the panel's width is drawn smaller, else at the top), top (in a row above the plot, two rows or shorter names where one does not fit; the plot keeps its width, losing that row's height), inside (over the plot's top right corner, the plot keeping its width) or alone (none in the panel; {stem}_key, one key for a row of panels) (default: right)"
    )]
    pub key_at: Option<KeyAt>,

    #[arg(
        long = "yscale",
        value_name = "Q=linear|sqrt|log1p,...",
        value_delimiter = ',',
        help = "With --style view, a histogram's count axis, per quantity (pv=sqrt) or one for all (log1p), e.g. pv=sqrt,cells=log1p (default: sqrt)"
    )]
    pub yscale: Vec<String>,

    #[arg(
        long = "marginal-scale",
        value_enum,
        value_name = "linear|sqrt|log1p",
        help = "With --style view, how deep a marginal's bar is for its sites, its axis ticked at a round count over its busiest bar (default: log1p)"
    )]
    pub marginal_scale: Option<CountScale>,

    #[arg(
        long = "count-scale",
        value_enum,
        value_name = "linear|sqrt|log1p",
        help = "With --style view, how dark a 2D square is for its sites (default: log1p)"
    )]
    pub count_scale: Option<CountScale>,

    #[arg(
        long = "ramp",
        value_name = "LOW,HIGH|NAME",
        help = "With --style view, the colours of a 2D panel's kept squares, from the fewest sites to the most (on --count-scale), and of its marginals' kept bars (the most's): two to nine colours, #rrggbb or a name (white, black, gray, red, blue, green, orange, purple, yellow), e.g. white,#d62728; or a named ramp: grays, reds, blues, ylgnbu, viridis (default: the palette's kept ink); the removed squares and the MUT arm keep theirs"
    )]
    pub ramp: Option<String>,

    #[arg(
        long = "threshold-colour",
        visible_alias = "threshold-color",
        value_name = "COLOUR",
        help = "With --style view, the colour of the thresholds' dashed lines on the histograms and of a 2D panel's dashed box, #rrggbb or a name, e.g. black or #d62728 (default: red lines, a black box)"
    )]
    pub threshold_colour: Option<String>,

    #[arg(
        long = "marginal-kept",
        value_enum,
        value_name = "all|panel",
        help = "With --style view, what a 2D panel's marginals count as kept: all (passing every threshold, as its squares are; default) or panel (passing the two thresholds on the panel's own axes, so the kept part of each bar ends at the box's edge)"
    )]
    pub marginal_kept: Option<super::arms::KeptBy>,

    #[arg(
        long = "scale",
        value_name = "Q=linear|log,...",
        value_delimiter = ',',
        help = "A quantity's scale, in every figure of it, e.g. pv=log,fold=linear (default: -log10 p, edit ratio and log odds linear; fold and counts log); log odds is linear only"
    )]
    pub scale: Vec<String>,

    #[arg(
        long = "bins",
        value_name = "N",
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=1000),
        help = "Bins along each quantity's axis (default: about one per 6 points of the figure's width, 20 to 60); with --style view, about this many in each histogram and a pair's squares a side (default: the view's own)"
    )]
    pub bins: Option<usize>,

    #[arg(
        long = "xlab",
        value_name = "Q=TEXT",
        allow_hyphen_values = true,
        help = "A quantity's axis title on x, e.g. pv=\"-log10 p (wt vs. mut)\"; repeat for more (default: the quantity's name)"
    )]
    pub xlab: Vec<String>,

    #[arg(
        long = "ylab",
        value_name = "Q=TEXT",
        allow_hyphen_values = true,
        help = "A quantity's axis title on y, e.g. cells=\"cells with a converted read\"; sites=TEXT for the histograms' counts; repeat for more"
    )]
    pub ylab: Vec<String>,
}

/// Axis titles asked for (`--xlab`, `--ylab`), by quantity: a quantity's
/// title on the axis it was given for, else on the other, else its name;
/// the histograms' counts their own.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Titles {
    x: [Option<String>; 7],
    y: [Option<String>; 7],
    sites: Option<String>,
}

impl Titles {
    /// `--xlab` and `--ylab`, each `q=TEXT`.
    fn parse(xlab: &[String], ylab: &[String]) -> anyhow::Result<Self> {
        let mut t = Self::default();
        for (flag, given, y) in [("--xlab", xlab, false), ("--ylab", ylab, true)] {
            for g in given {
                let Some((q, text)) = g.split_once('=') else {
                    anyhow::bail!("{flag} {g}: write a title as q=TEXT, e.g. pv=TEXT");
                };
                let q = q.trim();
                if q.eq_ignore_ascii_case("sites") {
                    anyhow::ensure!(y, "{flag} {g}: sites are counted on y; give it with --ylab");
                    t.sites = Some(text.into());
                    continue;
                }
                let q = Quantity::parse(q).map_err(|e| anyhow::anyhow!("{flag} {g}: {e}"))?;
                let at = if y { &mut t.y } else { &mut t.x };
                at[q as usize] = Some(text.into());
            }
        }
        Ok(t)
    }

    /// `q`'s title on x.
    pub fn x(&self, q: Quantity) -> String {
        let at = q as usize;
        (self.x[at].as_ref().or(self.y[at].as_ref())).map_or_else(|| q.title().into(), Clone::clone)
    }

    /// `q`'s title on y.
    pub fn y(&self, q: Quantity) -> String {
        let at = q as usize;
        (self.y[at].as_ref().or(self.x[at].as_ref())).map_or_else(|| q.title().into(), Clone::clone)
    }

    /// The histograms' counts' title.
    pub fn sites(&self) -> String {
        self.sites_or("sites")
    }

    /// `q`'s title on x as asked (on either axis), else `own`.
    pub fn x_or(&self, q: Quantity, own: &str) -> String {
        let at = q as usize;
        (self.x[at].as_ref().or(self.y[at].as_ref())).map_or_else(|| own.into(), Clone::clone)
    }

    /// `q`'s title on y as asked (on either axis), else `own`.
    pub fn y_or(&self, q: Quantity, own: &str) -> String {
        let at = q as usize;
        (self.y[at].as_ref().or(self.x[at].as_ref())).map_or_else(|| own.into(), Clone::clone)
    }

    /// The histograms' counts' title as asked, else `own`.
    pub fn sites_or(&self, own: &str) -> String {
        self.sites.clone().unwrap_or_else(|| own.into())
    }

    /// The title asked for `q` on x (`y` false) or on y, none for its own.
    pub fn given(&self, q: Quantity, y: bool) -> Option<&str> {
        let at = if y { &self.y } else { &self.x };
        at[q as usize].as_deref()
    }

    /// The histograms' counts' title asked for, none for their own.
    pub fn given_sites(&self) -> Option<&str> {
        self.sites.as_deref()
    }

    /// Title `q` on x (`y` false) or on y `text`, or its own with none.
    pub fn set(&mut self, q: Quantity, y: bool, text: Option<String>) {
        let at = if y { &mut self.y } else { &mut self.x };
        at[q as usize] = text;
    }

    /// Title the histograms' counts `text`, or their own with none.
    pub fn set_sites(&mut self, text: Option<String>) {
        self.sites = text;
    }
}

/// How the site panels are drawn (`--style`).
#[derive(
    clap::ValueEnum,
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "snake_case")]
pub enum PanelStyle {
    /// The paper panel: kept and removed stacked, the thresholds as edges.
    #[default]
    Panel,
    /// As the qc view draws its knob's histogram and its 2D panel.
    View,
}

/// A marginal's depth, as a fraction of the plot's shorter side, unless
/// `--marginal-size` says otherwise: the panel style's own.
pub const MARGINAL_SIZE: f64 = 0.2;

/// What was asked for, or Auto's.
#[derive(Debug, Clone, PartialEq)]
pub enum Choice<T> {
    Auto,
    Given(Vec<T>),
}

/// [`PanelArgs`] read and checked.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    pub hist: Choice<Quantity>,
    pub pairs: Choice<(Quantity, Quantity)>,
    pub marginals: Marginals,
    pub style: PanelStyle,
    /// The view style's marginals' depth, a fraction of the plot's
    /// shorter side.
    pub marginal_size: f64,
    /// The view style's arms: which in front, how the p-value reads,
    /// whether the key stays alone.
    pub look: ArmLook,
    pub scales: Scales,
    pub bins: Option<usize>,
    pub titles: Titles,
}

impl Plan {
    /// Bins along each quantity's axis on `sheet`: as asked, else about
    /// one per 6 points of a histogram's width.
    pub fn bins_on(&self, sheet: &Sheet) -> usize {
        (self.bins).unwrap_or_else(|| data_beans::figure::bins_for(sheet.fill(|| HIST_OWN).0))
    }
}

/// `items` read by `one`, Auto when there are none, none for `none`.
fn choice<T>(
    items: &[String],
    one: impl Fn(&str) -> anyhow::Result<T>,
) -> anyhow::Result<Choice<T>> {
    match items {
        [] => Ok(Choice::Auto),
        [none] if none.trim().eq_ignore_ascii_case("none") => Ok(Choice::Given(Vec::new())),
        _ => items
            .iter()
            .map(|s| one(s))
            .collect::<anyhow::Result<_>>()
            .map(Choice::Given),
    }
}

impl PanelArgs {
    /// The panels asked for, checked before anything is read.
    pub fn plan(&self) -> anyhow::Result<Plan> {
        let hist = choice(&self.hist, |s| {
            Quantity::parse(s).map_err(|e| anyhow::anyhow!("--hist: {e}"))
        })?;
        let pairs = choice(&self.pairs, |s| {
            let Some((x, y)) = s.split_once(':') else {
                anyhow::bail!("--pairs {s}: write a pair as x:y, e.g. pv:cells");
            };
            let q = |v: &str| Quantity::parse(v).map_err(|e| anyhow::anyhow!("--pairs {s}: {e}"));
            let (x, y) = (q(x)?, q(y)?);
            anyhow::ensure!(x != y, "--pairs {s}: two different quantities");
            Ok((x, y))
        })?;
        let mut marginals = Marginals::default();
        for m in &self.marginals {
            match m.trim().to_ascii_lowercase().as_str() {
                "top" => marginals.top = true,
                "right" => marginals.right = true,
                "both" => (marginals.top, marginals.right) = (true, true),
                "none" => {}
                other => anyhow::bail!("--marginals {other}: top, right or both"),
            }
        }
        let marginal_size = match self.marginal_size {
            None => MARGINAL_SIZE,
            Some(f) => {
                anyhow::ensure!(
                    self.style == PanelStyle::View,
                    "--marginal-size sizes the marginals of --style view only: the panel style's are a fifth of the plot"
                );
                anyhow::ensure!(
                    f > 0.0 && f < 1.0,
                    "--marginal-size {f}: a fraction of the plot's shorter side, above 0 and below 1"
                );
                f
            }
        };
        let view_only = |given: bool, flag: &str| {
            let view = self.style == PanelStyle::View;
            anyhow::ensure!(
                !given || view,
                "{flag} sets how --style view draws the panels; the panel style has no arms to draw apart"
            );
            Ok(())
        };
        view_only(self.p_axis.is_some(), "--p-axis")?;
        view_only(self.front.is_some(), "--front")?;
        view_only(self.key, "--key")?;
        view_only(self.key_at.is_some(), "--key-at")?;
        view_only(!self.yscale.is_empty(), "--yscale")?;
        view_only(self.marginal_scale.is_some(), "--marginal-scale")?;
        view_only(self.count_scale.is_some(), "--count-scale")?;
        view_only(self.ramp.is_some(), "--ramp")?;
        view_only(self.threshold_colour.is_some(), "--threshold-colour")?;
        view_only(self.marginal_kept.is_some(), "--marginal-kept")?;
        let own = ArmLook::default();
        let look = ArmLook {
            front: self.front.unwrap_or_default(),
            p_axis: self.p_axis.unwrap_or_default(),
            key: self.key,
            key_at: self.key_at.unwrap_or_default(),
            y: YScales::parse(&self.yscale)?,
            marginal: self.marginal_scale.unwrap_or(own.marginal),
            count: self.count_scale.unwrap_or(own.count),
            ramp: self
                .ramp
                .as_deref()
                .map(super::arms::Ramp::parse)
                .transpose()?,
            threshold: (self.threshold_colour.as_deref())
                .map(|c| {
                    super::arms::colour(c).map_err(|e| anyhow::anyhow!("--threshold-colour: {e}"))
                })
                .transpose()?,
            kept_by: self.marginal_kept.unwrap_or_default(),
        };
        Ok(Plan {
            hist,
            pairs,
            marginals,
            style: self.style,
            marginal_size,
            look,
            scales: Scales::parse(&self.scale)?,
            bins: self.bins,
            titles: Titles::parse(&self.xlab, &self.ylab)?,
        })
    }
}

/// `faba qc`'s site thresholds, each optional: given ones are drawn and
/// cut the sites into kept and removed; with none, every site is kept.
#[derive(Args, Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct Thresholds {
    #[arg(
        long = "site-max-pv",
        help = "Draw and apply `faba qc`'s p-value threshold (marginal, no FDR)"
    )]
    pub site_max_pv: Option<f32>,
    #[arg(
        long = "site-min-log-odds",
        allow_hyphen_values = true,
        help = "m6A: draw and apply the raw log odds threshold"
    )]
    pub site_min_log_odds: Option<f32>,
    #[arg(
        long = "site-min-coverage",
        help = "Draw and apply the total coverage threshold"
    )]
    pub site_min_coverage: Option<u64>,
    #[arg(
        long = "site-min-converted",
        help = "Draw and apply the converted reads threshold"
    )]
    pub site_min_converted: Option<u64>,
    #[arg(
        long = "site-min-edit-ratio",
        help = "Draw and apply the lower edit ratio threshold"
    )]
    pub site_min_edit_ratio: Option<f32>,
    #[arg(
        long = "site-max-edit-ratio",
        help = "Draw and apply the upper edit ratio threshold"
    )]
    pub site_max_edit_ratio: Option<f32>,
    #[arg(
        long = "site-min-fold",
        help = "m6A: draw and apply the fold over control threshold"
    )]
    pub site_min_fold: Option<f32>,
    #[arg(
        long = "site-min-cells",
        help = "Draw and apply the cells threshold (cells with a converted read, of every cell)"
    )]
    pub site_min_cells: Option<usize>,
}

/// How sites are cut: the thresholds, those drawn and applied, and those
/// Auto ranks its pairs at.
#[derive(Debug, Clone)]
pub struct Cut {
    pub filter: SiteFilterArgs,
    pub on: Vec<Criterion>,
    pub auto: SiteFilterArgs,
}

impl Cut {
    /// Nothing cut; Auto at `faba qc`'s defaults.
    #[cfg(test)]
    pub fn none() -> Self {
        Thresholds::default().cut()
    }

    /// `faba qc`'s cut at `f`: every threshold that is not off.
    #[cfg(test)]
    pub fn of_qc(f: &SiteFilterArgs) -> Self {
        Self {
            filter: f.clone(),
            on: (Criterion::ALL.into_iter())
                .filter(|c| !c.is_off(f))
                .collect(),
            auto: f.clone(),
        }
    }
}

impl Thresholds {
    /// The thresholds given over ones that keep everything, applied; Auto
    /// ranks at them over `faba qc`'s defaults.
    pub fn cut(&self) -> Cut {
        let given = [
            (
                Criterion::MinCoverage,
                self.site_min_coverage.map(|v| v as f64),
            ),
            (
                Criterion::MinConverted,
                self.site_min_converted.map(|v| v as f64),
            ),
            (
                Criterion::MinEditRatio,
                self.site_min_edit_ratio.map(f64::from),
            ),
            (
                Criterion::MaxEditRatio,
                self.site_max_edit_ratio.map(f64::from),
            ),
            (Criterion::MinFold, self.site_min_fold.map(f64::from)),
            (Criterion::MinLogOdds, self.site_min_log_odds.map(f64::from)),
            (Criterion::MaxPv, self.site_max_pv.map(f64::from)),
            (Criterion::MinCells, self.site_min_cells.map(|v| v as f64)),
        ];
        let mut filter = SiteFilterArgs::permissive();
        let mut auto = SiteFilterArgs::default_values();
        let mut on = Vec::new();
        for (c, v) in given {
            if let Some(v) = v {
                c.set(&mut filter, v);
                c.set(&mut auto, v);
                on.push(c);
            }
        }
        Cut { filter, on, auto }
    }
}

/// One quantity as every panel of it draws it: its ruler, each site's
/// value on it (infinities at its ends), and its thresholds' edges.
struct Drawn<'a> {
    ruler: Ruler,
    values: Cow<'a, [f32]>,
    edges: Vec<f64>,
}

impl<'a> Drawn<'a> {
    /// Quantity `q` on `scale` over its shown values `shown`, spread as
    /// `spread` (see [`super::scale::spread`]), reaching the thresholds of `f` in `on`, in
    /// about `bins` bins; `shown` copied only when it has an infinity to
    /// put at an end.
    fn new(
        (q, scale): (Quantity, Scale),
        mut shown: Cow<'a, [f32]>,
        spread: Option<(f64, f64)>,
        (f, on): (&SiteFilterArgs, &[Criterion]),
        bins: usize,
    ) -> Self {
        let (ruler, edges) = Ruler::cut((q, scale), spread, f, |c| on.contains(&c), bins);
        if shown.iter().any(|v| v.is_infinite()) {
            ruler.clamp(shown.to_mut());
        }
        Self {
            ruler,
            values: shown,
            edges,
        }
    }
}

/// One modality's sites, laid out as panels.
pub struct Panels<'a> {
    modality: &'a str,
    hists: Vec<Quantity>,
    pairs: Vec<(Quantity, Quantity)>,
    marginals: Marginals,
    kept: Option<Vec<bool>>,
    drawn: [Option<Drawn<'a>>; 7],
    titles: Titles,
    sheet: Sheet,
    /// With `--style view`, every panel as the view draws it, drawn while
    /// the sites were at hand.
    view: Option<Vec<Figure>>,
}

/// A pair of quantities drawn two at a time, `(x, y)`.
type Pair = (Quantity, Quantity);

/// A histogram's file stem.
fn hist_stem(modality: &str, q: Quantity) -> String {
    format!("{modality}_hist_{}", q.key())
}

/// A pair's file stem.
fn pair_stem(modality: &str, (x, y): Pair) -> String {
    format!("{modality}_{}_vs_{}", y.key(), x.key())
}

/// `hists` and `pairs` of `s` at the thresholds `f` as the qc view draws
/// them (`--style view`), each the size `sheet` asks (else its own): a
/// histogram as the view draws its knob's
/// ([`crate::qc::site_tui::view_hist`]), a pair as its 2D panel
/// ([`pairs::view_pair`]), in `plan`'s bins (else the view's) with its
/// marginals at their size. The cells are named as counted: over the
/// cells the cell QC kept when `kept_cells`, else over every cell.
fn view_figures(
    s: &Sites,
    (f, kept_cells): (&SiteFilterArgs, bool),
    (hists, pairs): (&[Quantity], &[Pair]),
    plan: &Plan,
    sheet: Sheet,
) -> Vec<Figure> {
    use rayon::prelude::*;
    let modality = &s.table.modality;
    let (hist, pair) = (sheet.fill(|| HIST_OWN), sheet.fill(|| PAIR_OWN));
    let (hist_style, pair_style) = (sheet.style(hist), sheet.style(pair));
    let bins = plan.bins.unwrap_or(pairs::BINS);
    let hists = (hists.par_iter()).flat_map_iter(|&q| {
        let stem = hist_stem(modality, q);
        let (svg, key) = super::site_tui::view_hist(
            s,
            (f, q),
            (plan.bins, kept_cells),
            (hist, &hist_style),
            (plan.look, &plan.titles),
        );
        let key = key.map(|svg| Figure {
            stem: format!("{stem}_key"),
            svg,
        });
        std::iter::once(Figure { stem, svg }).chain(key)
    });
    let marginals = (plan.marginals, plan.marginal_size);
    let how = pairs::PairLook {
        arms: plan.look,
        titles: Some(&plan.titles),
    };
    let pairs = (pairs.par_iter()).map(|&p| Figure {
        stem: pair_stem(modality, p),
        svg: pairs::view_pair(s, p, f, bins, (pair, &pair_style), (marginals, &how)),
    });
    hists.chain(pairs).collect()
}

/// The histograms and pairs `plan` asks for of `modality`'s sites, which
/// `criteria` apply to, Auto's pairs from `auto`. A quantity asked for that
/// the sites lack (fold without a control, cells without `_site` matrices)
/// is refused when `strict`, else left out with a warning.
fn chosen(
    modality: &str,
    criteria: &[Criterion],
    plan: &Plan,
    strict: bool,
    auto: impl FnOnce() -> Vec<Pair>,
) -> anyhow::Result<(Vec<Quantity>, Vec<Pair>)> {
    let applies = |q: Quantity| q.criteria().iter().any(|c| criteria.contains(c));
    let lacking = |q: Quantity| -> anyhow::Result<bool> {
        if applies(q) {
            return Ok(false);
        }
        let why = match q {
            Quantity::Cells => "no `_site` matrix to count cells from",
            _ => "no control arm",
        };
        let said = format!("{modality} sites have no {}: {why}", q.key());
        anyhow::ensure!(!strict, "{said}");
        log::warn!("{said}; not drawn");
        Ok(true)
    };
    let mut hists = Vec::new();
    for q in match &plan.hist {
        Choice::Auto => Quantity::ALL.into_iter().filter(|&q| applies(q)).collect(),
        Choice::Given(v) => v.clone(),
    } {
        if !lacking(q)? && !hists.contains(&q) {
            hists.push(q);
        }
    }
    let mut pairs = Vec::new();
    for (x, y) in match &plan.pairs {
        Choice::Auto => auto(),
        Choice::Given(v) => v.clone(),
    } {
        if !lacking(x)? && !lacking(y)? && !pairs.contains(&(x, y)) {
            pairs.push((x, y));
        }
    }
    Ok((hists, pairs))
}

/// The quantities `hists` and `pairs` draw, each once.
fn used(hists: &[Quantity], pairs: &[(Quantity, Quantity)]) -> Vec<Quantity> {
    let mut out: Vec<Quantity> = Vec::new();
    for q in (hists.iter().copied()).chain(pairs.iter().flat_map(|&(x, y)| [x, y])) {
        if !out.contains(&q) {
            out.push(q);
        }
    }
    out
}

impl<'a> Panels<'a> {
    /// The panels `plan` asks for over `t` (with its cells per site, when
    /// `_site` matrices gave them), cut by `cut`, laid out for `sheet`. A
    /// quantity asked for that the sites lack (fold without a control,
    /// cells without `_site` matrices) is refused when `strict`, else left
    /// out with a warning.
    pub fn new(
        t: &'a SiteTable,
        n_cells: Option<&[usize]>,
        cut: &Cut,
        plan: &Plan,
        sheet: Sheet,
        strict: bool,
    ) -> anyhow::Result<Self> {
        use rayon::prelude::*;
        let criteria: Vec<Criterion> = (Criterion::ALL.into_iter())
            .filter(|c| c.applies(t, n_cells.is_some()))
            .collect();
        let bins = plan.bins_on(&sheet);
        let cache = Cache::new(plan.scales, bins);
        // Auto ranks its pairs on the sites failing its thresholds.
        let fails = match plan.pairs {
            Choice::Auto => fails_at(t, n_cells, &criteria, &cut.auto),
            Choice::Given(_) => vec![0; t.len()],
        };
        let sites = Sites {
            table: t,
            n_cells,
            fails: &fails,
            criteria: &criteria,
            cache: &cache,
        };
        let (hists, pairs) = chosen(&t.modality, &criteria, plan, strict, || auto_ranked(&sites))?;
        let used = used(&hists, &pairs);
        used.par_iter().for_each(|&q| {
            sites.column(q);
        });
        let on: Vec<Criterion> = (cut.on.iter().copied())
            .filter(|c| criteria.contains(c))
            .collect();
        let kept = (!on.is_empty()).then(|| {
            (0..t.len())
                .into_par_iter()
                .map(|i| {
                    let cells = n_cells.map(|n| n[i]);
                    !on.iter().any(|c| c.fails(&cut.filter, t, i, cells))
                })
                .collect()
        });
        // The view reads the sites failing the cut's thresholds.
        let view = (plan.style == PanelStyle::View).then(|| {
            let fails = fails_at(t, n_cells, &criteria, &cut.filter);
            let s = Sites {
                fails: &fails,
                ..sites
            };
            // Over every cell: `faba site-plot` runs no cell QC.
            view_figures(&s, (&cut.filter, false), (&hists, &pairs), plan, sheet)
        });
        let mut columns = cache.into_columns();
        let mut drawn: [Option<Drawn>; 7] = Default::default();
        for q in used {
            let column = columns[q as usize]
                .take()
                .expect("a quantity drawn is read");
            drawn[q as usize] = Some(Drawn::new(
                (q, plan.scales.of(q)),
                Cow::Owned(column.shown),
                column.spread,
                (&cut.filter, &on),
                bins,
            ));
        }
        Ok(Self {
            modality: &t.modality,
            hists,
            pairs,
            marginals: plan.marginals,
            kept,
            drawn,
            titles: plan.titles.clone(),
            sheet,
            view,
        })
    }

    /// The panels `plan` asks for over the sites a view holds, at its
    /// thresholds `f`, as [`Panels::new`] draws them at `faba qc`'s cut of
    /// `f` ([`Cut::of_qc`]), a quantity the sites lack left out: the
    /// values read off the view's columns, the sites kept off its failed
    /// checks, Auto's pairs ranked on them.
    pub fn of_view(
        s: &Sites<'a>,
        f: &SiteFilterArgs,
        plan: &Plan,
        sheet: Sheet,
    ) -> anyhow::Result<Self> {
        use rayon::prelude::*;
        let (hists, pairs) = chosen(&s.table.modality, s.criteria, plan, false, || {
            auto_ranked(s)
        })?;
        let on: Vec<Criterion> = (s.criteria.iter().copied())
            .filter(|c| !c.is_off(f))
            .collect();
        let bits = on.iter().fold(0, |b, c| b | c.bit());
        let kept = (!on.is_empty()).then(|| s.fails.par_iter().map(|&m| m & bits == 0).collect());
        let bins = plan.bins_on(&sheet);
        let mut drawn: [Option<Drawn>; 7] = Default::default();
        for q in used(&hists, &pairs) {
            let column = s.column(q);
            drawn[q as usize] = Some(Drawn::new(
                (q, plan.scales.of(q)),
                Cow::Borrowed(&column.shown),
                column.spread,
                (f, &on),
                bins,
            ));
        }
        let view = (plan.style == PanelStyle::View)
            .then(|| view_figures(s, (f, true), (&hists, &pairs), plan, sheet));
        Ok(Self {
            modality: &s.table.modality,
            hists,
            pairs,
            marginals: plan.marginals,
            kept,
            drawn,
            titles: plan.titles.clone(),
            sheet,
            view,
        })
    }

    /// The histograms drawn, in order.
    #[cfg(test)]
    pub fn hists(&self) -> &[Quantity] {
        &self.hists
    }

    /// The pairs drawn, `(x, y)`, in order.
    #[cfg(test)]
    pub fn pairs(&self) -> &[(Quantity, Quantity)] {
        &self.pairs
    }

    fn drawn(&self, q: Quantity) -> &Drawn<'a> {
        self.drawn[q as usize]
            .as_ref()
            .expect("a quantity drawn is laid out")
    }

    /// How a figure `size` points looks: the font of the sheet, in its
    /// palette's kept and removed colours.
    fn look(&self, size: (f64, f64)) -> Look {
        let style = self.sheet.style(size);
        let roles = style.palette.roles();
        Look {
            font: style.font,
            panel_only: style.panel_only,
            kept: roles.kept.into(),
            removed: roles.removed.into(),
        }
    }

    /// The title of a figure, none with `--panel-only`.
    fn title(&self, what: String) -> String {
        match self.sheet.look.panel_only {
            true => String::new(),
            false => format!("{} {what}", self.modality),
        }
    }

    /// The histogram of `q`.
    pub fn hist(&self, q: Quantity) -> HistPanel<'_> {
        let d = self.drawn(q);
        HistPanel {
            axis: d.ruler.axis,
            bins: d.ruler.bins,
            values: &d.values,
            kept: self.kept.as_deref(),
            edges: d.edges.clone(),
            x_title: self.titles.x(q),
            y_title: self.titles.sites(),
            title: self.title(q.title().into()),
        }
    }

    /// The histogram of `q` as a figure `w` by `h` points.
    pub fn hist_figure(&self, q: Quantity, (w, h): (f64, f64)) -> String {
        self.hist(q).svg(&self.look((w, h)), w, h)
    }

    /// `y` against `x`, with the marginals asked for.
    pub fn pair(&self, x: Quantity, y: Quantity) -> DensityPanel<'_> {
        let (dx, dy) = (self.drawn(x), self.drawn(y));
        DensityPanel {
            x: dx.ruler.axis,
            x_bins: dx.ruler.bins,
            y: dy.ruler.axis,
            y_bins: dy.ruler.bins,
            xs: &dx.values,
            ys: &dy.values,
            kept: self.kept.as_deref(),
            x_edges: dx.edges.clone(),
            y_edges: dy.edges.clone(),
            x_title: self.titles.x(x),
            y_title: self.titles.y(y),
            title: self.title(format!("{} by {}", y.title(), x.title())),
            marginals: self.marginals,
        }
    }

    /// Every panel as a figure: `{modality}_hist_{q}`, then
    /// `{modality}_{y}_vs_{x}`, each the size asked (else its own), in
    /// the style asked.
    pub fn figures(&self) -> Vec<Figure> {
        use rayon::prelude::*;
        if let Some(view) = &self.view {
            return view.clone();
        }
        let (hw, hh) = self.sheet.fill(|| HIST_OWN);
        let (pw, ph) = self.sheet.fill(|| PAIR_OWN);
        let pl = self.look((pw, ph));
        let hists = (self.hists.par_iter()).map(|&q| Figure {
            stem: hist_stem(self.modality, q),
            svg: self.hist_figure(q, (hw, hh)),
        });
        let pairs = (self.pairs.par_iter()).map(|&(x, y)| Figure {
            stem: pair_stem(self.modality, (x, y)),
            svg: self.pair(x, y).svg(&pl, pw, ph),
        });
        hists.chain(pairs).collect()
    }
}

/// Per site of `t`, a [`Criterion::bit`] for each of `criteria` it fails
/// at `f`.
fn fails_at(
    t: &SiteTable,
    n_cells: Option<&[usize]>,
    criteria: &[Criterion],
    f: &SiteFilterArgs,
) -> Vec<u8> {
    use rayon::prelude::*;
    (0..t.len())
        .into_par_iter()
        .map(|i| {
            (criteria.iter())
                .filter(|c| c.fails(f, t, i, n_cells.map(|n| n[i])))
                .fold(0, |b, &c| b | c.bit())
        })
        .collect()
}

/// The pairs `faba qc` picks over `s`, at the thresholds its failed checks
/// were read at.
fn auto_ranked(s: &Sites) -> Vec<(Quantity, Quantity)> {
    (pairs::auto_pairs(s, pairs::PANELS).into_iter())
        .map(|p| (p.x, p.y))
        .collect()
}

/// Cells per site of `t`, aligned to it: every cell with a converted read
/// at the site, over every cell of every batch's `_site` matrix (no cell
/// QC); none without such a matrix. The control arm's batches of `dir`
/// ([`super::sites::control_batches`]) counted alone are set on `t`, to
/// draw that arm beside the signal's.
pub fn site_cells(
    layout: &InputLayout,
    t: &mut SiteTable,
    dir: &str,
) -> anyhow::Result<Option<Vec<usize>>> {
    let controls = super::sites::control_batches(std::path::Path::new(dir));
    let (mut acc, mut control) = (FxHashMap::default(), FxHashMap::default());
    let (mut any, mut any_control) = (false, false);
    for m in (layout.matrices.iter()).filter(|m| m.site_modality() == Some(&*t.modality)) {
        any = true;
        let data = open_matrix(&m.path)?;
        let cols = select_columns(data.as_ref(), &data.column_names()?, None);
        let nnz = row_nnz_over_columns(data.as_ref(), &cols)?;
        let rows = data.row_names()?;
        accumulate_site_cells(&rows, &nnz, t.converted_channel(), &mut acc);
        if controls.contains(&m.batch) {
            any_control = true;
            accumulate_site_cells(&rows, &nnz, t.converted_channel(), &mut control);
        }
    }
    super::sites::set_control_cells(t, any_control.then_some(&control));
    Ok(any.then(|| t.cells_per_site(&acc)))
}

#[derive(Args, Debug, serde::Serialize)]
#[command(
    mut_arg("fig_width", |a| a.visible_alias("width")),
    mut_arg("fig_height", |a| a.visible_alias("height"))
)]
pub struct SitePlotArgs {
    #[arg(
        help = "A faba output directory: its {modality}_sites.parquet and {batch}_{modality}_site matrices"
    )]
    pub input_dir: Box<str>,

    #[arg(
        short = 'o',
        long = "output",
        help = "Directory the figures go to, new or empty: nothing is written over"
    )]
    pub output: Box<str>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi",
        help = "Which site table to draw (default: m6a, else atoi)"
    )]
    pub modality: Option<Box<str>>,

    #[command(flatten)]
    pub panels: PanelArgs,

    #[command(flatten)]
    pub site: Thresholds,

    #[command(flatten)]
    pub size: SizeArgs,

    #[command(flatten)]
    pub style: StyleArgs,
}

/// `--modality`, else the first site modality with a table in `layout`.
fn modality_of(
    layout: &InputLayout,
    asked: Option<&str>,
    dir: &str,
) -> anyhow::Result<&'static str> {
    let present: Vec<&'static str> = (SITE_MODALITIES.iter().copied())
        .filter(|m| layout.site_tables.contains_key(*m))
        .collect();
    match asked {
        Some(a) => {
            let m = (SITE_MODALITIES.iter().copied())
                .find(|m| m.eq_ignore_ascii_case(a))
                .ok_or_else(|| {
                    anyhow::anyhow!("--modality {a}: one of {}", SITE_MODALITIES.join(", "))
                })?;
            anyhow::ensure!(present.contains(&m), "{dir}: no {m}_sites.parquet");
            Ok(m)
        }
        None => present.first().copied().ok_or_else(|| {
            anyhow::anyhow!(
                "{dir}: no site table ({}_sites.parquet)",
                SITE_MODALITIES.join("_sites.parquet or ")
            )
        }),
    }
}

pub fn run_site_plot(args: &SitePlotArgs) -> anyhow::Result<()> {
    let plan = args.panels.plan()?;
    let sheet = Sheet::new(
        crate::figure::size::Figure::SitePanels,
        &args.size,
        &args.style,
    )?;
    let (input, out) = (&*args.input_dir, &*args.output);
    if let Some(problem) = crate::tui::output_problem(out) {
        anyhow::bail!("-o {problem}: faba site-plot never writes over files; name another");
    }
    let layout = scan_input_dir(input)?;
    let modality = modality_of(&layout, args.modality.as_deref(), input)?;
    let record = crate::run_record::RunRecord::start("site-plot", out)
        .input("fileset", Some(input))
        .options(args);
    let outcome = draw(
        (&layout, input),
        modality,
        &args.site.cut(),
        &plan,
        sheet,
        out,
    );
    record.finish(&outcome);
    for path in outcome? {
        println!("{path}");
    }
    Ok(())
}

/// The figures of `modality` in `layout`, saved under `out` in parallel;
/// the paths, in the figures' order.
fn draw(
    (layout, dir): (&InputLayout, &str),
    modality: &str,
    cut: &Cut,
    plan: &Plan,
    sheet: Sheet,
    out: &str,
) -> anyhow::Result<Vec<String>> {
    use rayon::prelude::*;
    let mut t = read_site_table(&layout.site_tables[modality], modality)?;
    anyhow::ensure!(
        !t.is_empty(),
        "{modality}_sites.parquet has no sites to draw"
    );
    let n_cells = site_cells(layout, &mut t, dir)?;
    let figures = Panels::new(&t, n_cells.as_deref(), cut, plan, sheet, true)?.figures();
    anyhow::ensure!(
        !figures.is_empty(),
        "nothing to draw: --hist none and --pairs none"
    );
    std::fs::create_dir_all(out)?;
    let saved: Vec<Vec<String>> = (figures.par_iter())
        .map(|f| crate::figure::save_with_png(&f.svg, &format!("{out}/{}", f.stem)))
        .collect::<anyhow::Result<_>>()?;
    let paths: Vec<String> = saved.into_iter().flatten().collect();
    log::info!(
        "wrote {} figures of {modality} sites under {out}",
        figures.len()
    );
    Ok(paths)
}

#[cfg(test)]
#[path = "tests/site_plot.rs"]
mod tests;
