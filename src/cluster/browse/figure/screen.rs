//! The browser's screen as a figure: its three regions laid out as the
//! browser lays them out, each drawn from the same state the browser
//! draws it from. Top left the gene modules against the cell modules
//! (the grid figure, as its diagonal Hinton look or its blocks' bars, each
//! row's genes and shift range beside it); top right one gene module's
//! sites as the sites pane's bars, or each gene module's GO terms aligned
//! with its row, or nothing; below them the pinned genes (see
//! [`super::genes`]). The browser's `s` saves it, and `--screen` writes
//! it from a [`State`] the command line (or a saved view) gives.
//!
//! A [`State`] is plain data, serialisable, naming everything by the
//! names the browser shows (C3, M2, a gene's key or symbol), so a saved
//! view or another tool can draw the same screen again.

use super::{grid, GridArgs, Ink, Look, Refused};
use crate::cluster::browse::go::Table;
use crate::cluster::browse::view::Pin;
use crate::cluster::browse::Modules;
use crate::figure::size::Page;
use crate::figure::style::Sheet;
use crate::figure::{esc, Anchor, Canvas, FAINT, INK, MUTED};
use crate::site_analysis::miami::genemodel::GeneModel;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

/// How the grid orders and draws its rows: the browser's two.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Preset {
    /// Rows and columns in diagonal order, each block a Hinton square of
    /// its z, red above the row's mean and blue below.
    #[default]
    Diagonal,
    /// Rows by shift range, each block's reads as mirror bars.
    Blocks,
}

impl Preset {
    /// The preset that draws the grid in `look`.
    fn of(look: Look) -> Self {
        match look {
            Look::Bars => Preset::Blocks,
            _ => Preset::Diagonal,
        }
    }

    /// The grid's region's title, naming the GO terms when `go` are on it.
    fn title(self, go: bool) -> String {
        match go {
            true => format!("gene modules · {} · GO terms", self.said()),
            false => format!("gene modules · {}", self.said()),
        }
    }

    /// The grid look it draws with.
    pub fn look(self) -> Look {
        match self {
            Preset::Diagonal => Look::Hinton,
            Preset::Blocks => Look::Bars,
        }
    }

    fn said(self) -> &'static str {
        match self {
            Preset::Diagonal => "diagonal · Hinton, z",
            Preset::Blocks => "by shift range · bars",
        }
    }
}

/// What the region right of the grid shows.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Beside {
    /// The gene module in focus's sites, as the sites pane draws them.
    #[default]
    Sites,
    /// Each gene module's top GO terms on its row, q beside each.
    Go,
    /// Nothing: the grid takes the width.
    None,
}

/// The scale the pinned genes' mirror bars are drawn on: each side's
/// reads (the posterior-mean reads at a site, summed over the cells of a
/// track) as they are, their square root, or `log10(1 + reads)`, each
/// track's peak filling a side (the "all" track's its own, the cell
/// modules' one shared): the pileup's.
pub use crate::site_analysis::pileup::tui::Window;
pub use crate::site_analysis::pileup::YScale;

/// Where a pinned gene's model goes: over its tracks (its name under
/// them), or under them (its name over them).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum ModelAt {
    #[default]
    Top,
    Bottom,
}

/// How the pinned genes are drawn on `--screen`, besides their scale.
#[derive(clap::Args, Debug, Clone, Copy, Default, PartialEq)]
pub struct GenesLook {
    #[arg(
        long = "gene-model",
        value_enum,
        default_value_t = ModelAt::Top,
        requires = "screen",
        help = "Where each pinned gene's model goes on --screen: top (over its tracks, its name and module under them) or bottom (under its tracks, its name and module over them)"
    )]
    pub model_at: ModelAt,

    #[arg(
        long = "no-all",
        requires = "screen",
        help = "Leave the pinned genes' \"all\" track (every fitted cell) out on --screen: a track a cell module alone"
    )]
    pub no_all: bool,

    #[arg(
        long = "no-track-labels",
        requires = "screen",
        help = "Leave the pinned genes' track labels (each track's name and cells, left of the genes) out on --screen"
    )]
    pub no_track_labels: bool,

    #[arg(
        long = "two-line-names",
        requires = "screen",
        help = "Name each pinned gene on --screen in two lines, its symbol over its module (\"HOXA10-AS\" / \"M2\"), so a long symbol is not cut to fit beside its module"
    )]
    pub two_line_names: bool,

    #[arg(
        long = "genes-window",
        value_enum,
        requires = "screen",
        help = "What each pinned gene's strip spans on --screen: gene (its model, TSS to TES, as the browser has it; the default but beside a transposed grid) or sites (at its sites, padded by --genes-pad; the default beside a transposed grid, where the genes are drawn equally wide), as faba pileup's --window"
    )]
    pub window: Option<Window>,

    #[arg(
        long = "genes-pad",
        value_name = "BP",
        requires = "screen",
        help = "Bases either side of a pinned gene's sites where its strip is at them (--genes-window sites), from the first to the last; by default centred on them, their span and a quarter of it either side, 1 kb at least, within the gene (as faba pileup's --pad)"
    )]
    pub pad: Option<u32>,

    #[arg(
        long = "genes-xscale",
        value_enum,
        default_value = "linear",
        requires = "screen",
        help = "The pinned genes' strips' x scale on --screen: linear, sqrt or log10 (of 1 + bases) of the distance from the strip's left-most base, the left end spread and the right pressed, the models drawn on it too (as faba pileup's --xscale)"
    )]
    pub x_scale: YScale,
}

/// Which of the screen's regions a figure draws.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Region {
    /// The whole screen.
    #[default]
    All,
    /// The gene modules against the cell modules.
    Grid,
    /// What is right of the grid (--side).
    Side,
    /// The pinned genes.
    Genes,
}

/// What the screen shows, by the names the browser shows: everything a
/// screen is drawn from but the fit itself. Left out, a field is the
/// browser's own as it opens.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct State {
    /// The grid's order and look.
    pub preset: Preset,
    /// The cell modules first in the panes' order (C6, C2), the rest after
    /// them as they were.
    pub cell_order: Vec<String>,
    /// The gene module in focus (M3), whose sites are beside the grid;
    /// the grid's first row when none is named.
    pub module: Option<String>,
    /// What is beside the grid.
    pub side: Beside,
    /// The GO terms a row shows at most (beside the grid as GO).
    pub go_terms: Option<usize>,
    /// The pinned genes, by key or symbol, in order.
    pub pins: Vec<String>,
    /// The scale the pinned genes' bars are drawn on.
    pub genes_yscale: YScale,
}

/// The GO terms a row shows when [`State::go_terms`] names none.
pub const GO_TERMS: usize = 3;

/// A [`State`] with its names found in a fit.
#[derive(Debug, Clone, PartialEq)]
pub struct Placed {
    pub look: Look,
    /// The cell roles in the order the panes show them.
    pub columns: Vec<usize>,
    /// The gene role in focus; none when no module has sites.
    pub module: Option<usize>,
    pub side: Beside,
    pub go_terms: usize,
    pub pins: Vec<Pin>,
    pub genes_yscale: YScale,
}

/// A cell module named as C3 or 3, its role.
pub fn cell_module(m: &Modules, name: &str) -> anyhow::Result<usize> {
    let n: usize = (name.trim().trim_start_matches(['C', 'c']).parse())
        .map_err(|_| anyhow::anyhow!("a cell module is named as C3 or 3, not {name}"))?;
    anyhow::ensure!(
        (1..=m.cell_roles).contains(&n),
        "no cell module C{n}: the fit has C1 to C{}",
        m.cell_roles
    );
    Ok(n - 1)
}

/// Gene `name`'s pin (by key or symbol): its main gene module and its
/// place among that module's genes; refused for a gene with no sites.
pub fn pin_of(m: &Modules, name: &str) -> anyhow::Result<Pin> {
    let (key, role) = m
        .find_gene(name)
        .ok_or_else(|| anyhow::anyhow!("no gene {name} with sites in this fit"))?;
    let gene = (m.genes(role).iter())
        .position(|s| m.gene_of(role, s.start) == key)
        .ok_or_else(|| anyhow::anyhow!("{name} has no sites in its main module"))?;
    Ok(Pin { role, gene })
}

/// The first `n` genes of each of `roles` (each module's genes with the
/// most first, as the browser lists them), by key.
pub fn top_genes(m: &Modules, roles: &[usize], n: usize) -> Vec<String> {
    (roles.iter())
        .flat_map(|&g| (m.genes(g).iter().take(n)).map(move |s| m.gene_of(g, s.start).to_string()))
        .collect()
}

impl State {
    /// Its names found in `m`; refused for one `m` does not have.
    pub fn place(&self, m: &Modules) -> anyhow::Result<Placed> {
        let first = (self.cell_order.iter())
            .map(|c| cell_module(m, c))
            .collect::<anyhow::Result<Vec<usize>>>()?;
        let rest = m.columns.iter().filter(|c| !first.contains(c));
        let columns = first.iter().chain(rest).copied().collect();
        let look = self.preset.look();
        let rows = rows(m, look);
        let module = match &self.module {
            Some(name) => Some(crate::cluster::browse::text::gene_module(m, name)?),
            None => rows.first().copied(),
        };
        let mut pins: Vec<Pin> = Vec::new();
        for name in &self.pins {
            let pin = pin_of(m, name)?;
            if !pins.contains(&pin) {
                pins.push(pin);
            }
        }
        Ok(Placed {
            look,
            columns,
            module,
            side: self.side,
            go_terms: self.go_terms.unwrap_or(GO_TERMS),
            pins,
            genes_yscale: self.genes_yscale,
        })
    }

    /// `m`'s cell modules put in this state's order (see
    /// [`Modules::reorder_columns`]), so the screen drawn of it shows them so.
    pub fn order(&self, m: &mut Modules) -> anyhow::Result<()> {
        let first = (self.cell_order.iter())
            .map(|c| cell_module(m, c))
            .collect::<anyhow::Result<Vec<usize>>>()?;
        m.reorder_columns(&first)
    }
}

/// The grid's rows under `look`, as the browser orders them.
fn rows(m: &Modules, look: Look) -> Vec<usize> {
    match look {
        Look::Bars => m.order.clone(),
        _ => m.diagonal.clone(),
    }
}

/// What the screen is drawn with besides the fit and its state: the GO
/// terms (beside the grid as GO), the genes' models, the grid's extras,
/// the region, the title, and the page and look.
#[derive(Clone, Copy)]
pub struct Ask<'a> {
    pub go: Option<&'a Table>,
    pub models: &'a FxHashMap<Box<str>, GeneModel>,
    pub grid: &'a GridArgs,
    pub region: Region,
    pub genes: GenesLook,
    pub title: &'a str,
    pub sheet: Sheet,
}

/// The screen's own size when no side is asked for, in points: about a
/// page across, the browser's proportions.
const OWN: (f64, f64) = (612.0, 432.0);

/// The screen of `m` as `state` places it (its cell modules already in
/// its order: see [`State::order`]), as `ask` asks; why not, when a
/// region has nothing to draw or the page is too small for the grid.
pub fn draw(m: &Modules, state: &Placed, ask: &Ask) -> Result<String, String> {
    if state.columns != m.columns {
        return Err("the cell modules are not in the state's order: order them first".into());
    }
    let (w, h) = page(ask.sheet.page);
    let pinned = !state.pins.is_empty();
    let beside = match (state.side, state.module) {
        (Beside::Sites, None) => Beside::None,
        (side, _) => side,
    };
    // The text as the browser's: the screen's columns across the page,
    // unless a size is asked for.
    let (modules, sites) = columns(m, state);
    let across = match (ask.grid.transpose && ask.region == Region::All, beside) {
        // Transposed: the grid, a gene module a column, the genes beside
        // it (the sites under them, the GO terms under its columns).
        (true, _) => transposed_columns(m, state),
        (false, Beside::Sites) => modules + sites,
        (false, Beside::Go) => modules + GO_COLUMNS,
        (false, Beside::None) => modules,
    };
    let mut sheet = ask.sheet;
    sheet.page = (Some(w), Some(h));
    sheet.look.font_size = Some(
        (sheet.look.font_size).unwrap_or_else(|| (w / (across * COLUMN)).clamp(FONT.0, FONT.1)),
    );
    let ask = &Ask { sheet, ..*ask };
    let ink = Ink::of(&sheet);
    let f = ink.f;
    let framed = !ink.panel_only;
    match ask.region {
        Region::Side if beside == Beside::None => {
            return Err("nothing is beside the grid: show the sites or GO".into())
        }
        Region::Side if beside == Beside::Go && ask.go.is_none() => {
            return Err("no GO terms to show beside the grid".into())
        }
        Region::Genes if !pinned => return Err("no gene is pinned".into()),
        Region::Side if beside == Beside::Go && ask.grid.transpose => {
            return Err(
                "with --transpose the GO terms are under the grid's columns: draw --region grid"
                    .into(),
            )
        }
        _ => {}
    }
    let mut c = Canvas::new(w, h);
    // The parts drawn as figures of their own, placed over the canvas.
    let mut parts: Vec<String> = Vec::new();
    let head = match framed && ask.region == Region::All {
        true => 2.2 * f,
        false => 0.0,
    };
    if head > 0.0 {
        let room = w - f;
        let title = super::title_text(ask.title, f, room);
        c.text(f / 2.0, 1.5 * f, &esc(&title), f + 1.0, Anchor::Start, INK);
    }
    let body = (0.0, head, w, h - head);
    // The browser's split: the grids over the genes, half the screen
    // each; the grid and its side across as their needs are.
    let (top, genes) = match (ask.region, pinned) {
        (Region::All, true) => {
            let top_h = body.3 * 0.5;
            (
                Some((body.0, body.1, body.2, top_h)),
                Some((body.0, body.1 + top_h, body.2, body.3 - top_h)),
            )
        }
        (Region::All, false) | (Region::Grid | Region::Side, _) => (Some(body), None),
        (Region::Genes, _) => (None, Some(body)),
    };
    let module_name = state
        .module
        .map_or(String::new(), |g| format!("M{}", g + 1));
    let preset = Preset::of(state.look);
    if ask.grid.transpose && ask.region == Region::All {
        let parts = transposed(&mut c, m, state, ask, (body, beside), ink)?;
        return Ok(placed_over(c, &[parts]));
    }
    if let Some(top) = top {
        // The GO terms beside the grid, but for the grid alone; a
        // transposed grid's are its own, under its columns.
        let go_here = beside == Beside::Go
            && ask.go.is_some()
            && (ask.region != Region::Grid || ask.grid.transpose);
        let grid_args = GridArgs {
            ranges: true,
            go_labels: go_here.then_some(state.go_terms),
            ..ask.grid.clone()
        };
        match (ask.region, beside) {
            (Region::Side, Beside::Sites) | (Region::All, Beside::Sites) => {
                let (grid_area, sites_area) = match ask.region {
                    Region::Side => (None, top),
                    _ => {
                        let gw = top.2 * modules / (modules + sites);
                        let left = (top.0, top.1, gw, top.3);
                        (Some(left), (top.0 + gw, top.1, top.2 - gw, top.3))
                    }
                };
                if let Some(area) = grid_area {
                    let title = format!("gene modules · {}", preset.said());
                    let part = grid_region(
                        &mut c,
                        m,
                        (state.look, &grid_args),
                        ask,
                        area,
                        (&title, framed),
                        ink,
                    )?;
                    parts.push(part);
                }
                let n = m.members[state.module.unwrap_or(0)].len();
                let title = format!("sites · {module_name} · {n} sites · module order · bars");
                let inner = frame(&mut c, sites_area, (&title, framed), f);
                let g = state.module.unwrap_or(0);
                let pinned: Vec<usize> = (state.pins.iter())
                    .filter(|p| p.role == g)
                    .map(|p| p.gene)
                    .collect();
                super::sites::draw(&mut c, m, g, &pinned, inner, ink);
            }
            (Region::Side, _) => {
                // The GO terms alone, where they are beside the whole grid.
                let title = "GO terms";
                let inner = frame(&mut c, top, (title, framed), f);
                let sheet = sub_sheet(ask.sheet, inner);
                let (svg, placed) = grid::placed(
                    m,
                    state.look,
                    &grid_args,
                    ask.go,
                    "",
                    sheet.page,
                    Ink::of(&sheet),
                )
                .map_err(|e| refused(e, sheet))?;
                let go_x = placed.go_x.unwrap_or(0.0);
                parts.push(cropped(&svg, inner, go_x));
            }
            (_, _) => {
                let title = preset.title(go_here);
                let part = grid_region(
                    &mut c,
                    m,
                    (state.look, &grid_args),
                    ask,
                    top,
                    (&title, framed),
                    ink,
                )?;
                parts.push(part);
            }
        }
    }
    if let Some(area) = genes {
        let n = state.pins.len();
        let title =
            format!("genes · {n} pinned · posterior-mean reads in the main module, fitted cells");
        let inner = frame(&mut c, area, (&title, framed), f);
        let how = super::genes::How::of(state, ask);
        super::genes::draw(&mut c, m, &state.pins, ask.models, inner, ink, how);
    }
    Ok(placed_over(c, &parts))
}

/// `c` drawn, `parts` placed over it.
fn placed_over(c: Canvas, parts: &[String]) -> String {
    let drawn = c.finish();
    let end = drawn.rfind("</svg>").expect("a drawing ends its svg");
    format!("{}{}{}", &drawn[..end], parts.concat(), &drawn[end..])
}

/// The whole screen of a transposed grid in `body`: the grid (as narrow
/// as it is drawn) with the pinned genes beside it, in the gene modules'
/// order, each cell module's track on its grid row; the sites under them
/// when they are `beside` the grid. The grid's part, to place over `c`.
fn transposed(
    c: &mut Canvas,
    m: &Modules,
    state: &Placed,
    ask: &Ask,
    (body, beside): ((f64, f64, f64, f64), Beside),
    ink: Ink,
) -> Result<String, String> {
    let f = ink.f;
    let framed = !ink.panel_only;
    let (top, sites) = match beside {
        Beside::Sites => {
            let top_h = body.3 * 0.6;
            let low = (body.0, body.1 + top_h, body.2, body.3 - top_h);
            ((body.0, body.1, body.2, top_h), Some(low))
        }
        _ => (body, None),
    };
    let go_here = beside == Beside::Go && ask.go.is_some();
    let grid_args = GridArgs {
        ranges: true,
        go_labels: go_here.then_some(state.go_terms),
        ..ask.grid.clone()
    };
    let title = Preset::of(state.look).title(go_here);
    // The pins in the grid's columns' order, the genes beside it taking
    // what they want of the width, the grid at least a third.
    let order = rows(m, state.look);
    let mut pins = state.pins.clone();
    pins.sort_by_key(|p| {
        order
            .iter()
            .position(|&g| g == p.role)
            .unwrap_or(usize::MAX)
    });
    let how = super::genes::How::of(state, ask);
    let wanted = match pins.is_empty() {
        true => 0.0,
        false => super::genes::wanted(m, &pins, f, how) + 2.0 * f,
    };
    let mut area = (top.0, top.1, (top.2 - wanted).max(top.2 / 3.0), top.3);
    let draw_grid = |inner: (f64, f64, f64, f64)| {
        let sheet = sub_sheet(ask.sheet, inner);
        grid::placed(
            m,
            state.look,
            &grid_args,
            ask.go,
            "",
            sheet.page,
            Ink::of(&sheet),
        )
        .map_err(|e| refused(e, sheet))
    };
    let genes_area = |area: (f64, f64, f64, f64)| {
        let x = area.0 + area.2;
        (x, top.1, top.0 + top.2 - x, top.3)
    };
    // Drawn once to learn how wide it is, then framed that wide (the
    // whole width when squeezed too narrow for its cells: the genes then
    // narrow to what is left); lowered when the genes' "all" track and
    // their models or names over the first row would not fit over it.
    if !pins.is_empty() {
        let inner = inside(area, framed, f);
        let (inner, placed) = match draw_grid(inner) {
            Ok((_, placed)) => (inner, placed),
            Err(_) => {
                area.2 = top.2;
                let inner = inside(area, framed, f);
                (inner, draw_grid(inner)?.1)
            }
        };
        area.2 = (placed.width + area.2 - inner.2).min(area.2);
        let genes_top = inside(genes_area(area), framed, f).1;
        let over = super::genes::over(how.look, !ask.models.is_empty(), f, placed.pitch);
        let short = genes_top + over - (inner.1 + placed.rows_top);
        if short > 0.0 {
            (area.1, area.3) = (area.1 + short, area.3 - short);
        }
    }
    let inner = frame(c, area, (&title, framed), f);
    let (svg, placed) = draw_grid(inner)?;
    if !pins.is_empty() {
        let genes_area = genes_area(area);
        let n = pins.len();
        let said = format!("genes · {n} pinned · posterior-mean reads by cell module");
        let inner_genes = frame(c, genes_area, (&said, framed), f);
        let aligned = super::genes::Aligned {
            top: inner.1 + placed.rows_top,
            pitch: placed.pitch,
        };
        let how = super::genes::How {
            aligned: Some(aligned),
            ..how
        };
        super::genes::draw(c, m, &pins, ask.models, inner_genes, ink, how);
    }
    if let (Some(area), Some(g)) = (sites, state.module) {
        let n = m.members[g].len();
        let title = format!("sites · M{} · {n} sites · module order · bars", g + 1);
        let inner = frame(c, area, (&title, framed), f);
        let pinned: Vec<usize> = (state.pins.iter())
            .filter(|p| p.role == g)
            .map(|p| p.gene)
            .collect();
        super::sites::draw(c, m, g, &pinned, inner, ink);
    }
    Ok(nested(&svg, inner))
}

/// The page asked for, a side not asked for in the screen's own
/// proportions.
fn page(asked: Page) -> (f64, f64) {
    match asked {
        (Some(w), Some(h)) => (w, h),
        (Some(w), None) => (w, w * OWN.1 / OWN.0),
        (None, Some(h)) => (h * OWN.0 / OWN.1, h),
        (None, None) => OWN,
    }
}

/// A terminal column's width, in font sizes, and the text's least and
/// most size when none is asked for, in points.
const COLUMN: f64 = 0.6;
const FONT: (f64, f64) = (5.0, 9.0);
/// The columns the GO terms take beside the grid, about.
const GO_COLUMNS: f64 = 40.0;

/// The columns the grid and the sites need, as the browser splits its
/// panes in proportion to them: every cell module six wide.
fn columns(m: &Modules, state: &Placed) -> (f64, f64) {
    let cells = 6.0 * m.columns.len() as f64;
    let modules = 17.0 + cells + 8.0 + 2.0;
    let label = state.module.map_or(0, |g| {
        (m.members[g].iter())
            .map(|s| {
                let site = &m.sites[s.site];
                crate::site_analysis::pileup::query_symbol(&site.gene)
                    .chars()
                    .count()
                    + 1
                    + crate::tui::digits(site.pos)
                    + 1
            })
            .max()
            .unwrap_or(0)
            .min(24)
    });
    let sites = label as f64 + cells + 27.0 + 4.0;
    (modules, sites)
}

/// The columns a transposed screen needs across: the grid, every gene
/// module six wide past the cell modules' names, and the pinned genes in
/// one row past their tracks' labels.
fn transposed_columns(m: &Modules, state: &Placed) -> f64 {
    let grid = 8.0 + 6.0 * m.gene_roles as f64 + 2.0;
    let genes = match state.pins.is_empty() {
        true => 0.0,
        false => 12.0 + super::genes::one_row(&super::genes::beside(m, &state.pins)),
    };
    grid + genes
}

/// A region's frame and title, as the browser's panes have them, when
/// `framed`: the area inside it.
fn frame(
    c: &mut Canvas,
    (x, y, w, h): (f64, f64, f64, f64),
    (title, framed): (&str, bool),
    f: f64,
) -> (f64, f64, f64, f64) {
    let inner = inside((x, y, w, h), framed, f);
    if !framed {
        return inner;
    }
    let (fx, fy, fw, fh) = (x + 0.3 * f, y + 0.6 * f, w - 0.6 * f, h - 0.9 * f);
    c.line(fx, fy, fx + fw, fy, FAINT, 0.6);
    c.line(fx, fy + fh, fx + fw, fy + fh, FAINT, 0.6);
    c.line(fx, fy, fx, fy + fh, FAINT, 0.6);
    c.line(fx + fw, fy, fx + fw, fy + fh, FAINT, 0.6);
    let said = super::title_text(title, f - 1.0, fw - 2.0 * f);
    let tw = super::width(&said, f - 1.0) + 0.8 * f;
    c.rect(fx + f - 0.4 * f, fy - 0.6 * f, tw, 1.2 * f, "#ffffff");
    c.text(
        fx + f,
        fy + 0.3 * f,
        &esc(&said),
        f - 1.0,
        Anchor::Start,
        MUTED,
    );
    inner
}

/// The area inside a region's frame (see [`frame`]), or inside its
/// margin when it is not framed.
fn inside((x, y, w, h): (f64, f64, f64, f64), framed: bool, f: f64) -> (f64, f64, f64, f64) {
    if !framed {
        return (x + 0.3 * f, y + 0.3 * f, w - 0.6 * f, h - 0.6 * f);
    }
    let pad = 0.5 * f;
    let (fx, fy, fw, fh) = (x + 0.3 * f, y + 0.6 * f, w - 0.6 * f, h - 0.9 * f);
    (fx + pad, fy + 0.8 * f, fw - 2.0 * pad, fh - 0.8 * f - pad)
}

/// A sheet for a part drawn in `area`: its page, the screen's look but
/// the plot alone (the screen frames and titles its regions).
fn sub_sheet(sheet: Sheet, (_, _, w, h): (f64, f64, f64, f64)) -> Sheet {
    let mut look = sheet.look;
    look.panel_only = true;
    Sheet {
        page: (Some(w), Some(h)),
        look,
    }
}

/// Why a part was not drawn, said.
fn refused(e: Refused, sheet: Sheet) -> String {
    match e {
        Refused::Said(why) => why,
        Refused::Small(_) => {
            let inches = |p: Option<f64>| p.unwrap_or(0.0) / crate::figure::size::PT;
            format!(
                "the screen's grid has {:.2} by {:.2} inches, too little for its labels and \
                 cells: give a larger --fig-width and --fig-height, or a smaller --font-size",
                inches(sheet.page.0),
                inches(sheet.page.1)
            )
        }
    }
}

/// The grid in `area` under `look`, framed and titled, drawn as the grid
/// figure draws it (the plot alone) at that size, to place over `c`.
fn grid_region(
    c: &mut Canvas,
    m: &Modules,
    (look, grid_args): (Look, &GridArgs),
    ask: &Ask,
    area: (f64, f64, f64, f64),
    title: (&str, bool),
    ink: Ink,
) -> Result<String, String> {
    let inner = frame(c, area, title, ink.f);
    let sheet = sub_sheet(ask.sheet, inner);
    let (svg, _) = grid::placed(m, look, grid_args, ask.go, "", sheet.page, Ink::of(&sheet))
        .map_err(|e| refused(e, sheet))?;
    Ok(nested(&svg, inner))
}

/// A figure from [`Canvas::finish`] placed at `(x, y)`.
fn nested(svg: &str, (x, y, _, _): (f64, f64, f64, f64)) -> String {
    svg.replacen("<svg ", &format!("<svg x=\"{x}\" y=\"{y}\" "), 1)
}

/// The part of a figure right of `from` placed in `area`, at its own scale.
fn cropped(svg: &str, (x, y, w, h): (f64, f64, f64, f64), from: f64) -> String {
    let inner = svg.replacen("<svg ", &format!("<svg x=\"{}\" y=\"0\" ", -from), 1);
    format!(r#"<svg x="{x}" y="{y}" width="{w}" height="{h}" viewBox="0 0 {w} {h}">{inner}</svg>"#)
}
