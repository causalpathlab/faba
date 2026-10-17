//! The report's Model section: the blockmodel the tables come from, its
//! settings and its plate diagram, drawn as the run set its switches (see
//! [`crate::cluster::unified`]). The page styles the diagram's classes
//! from the shell's theme tokens.

use crate::cluster::unified::{strength, Recorded, Settings};
use serde::Serialize;
use std::fmt::Write;

/// What the page reads of the model: the run's settings and diagram, or
/// why there are none.
#[derive(Serialize)]
#[serde(untagged)]
pub(super) enum Section {
    Drawn(Box<Drawn>),
    Omitted { why: String },
}

/// The run's settings (`shiftSd` none with the shift off, `gamma` none
/// with symmetric weights, `strength` a number or `"inf"`, `cellSd` none at
/// finite strength, where the cells have no level) and its diagram.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Drawn {
    preset: Option<String>,
    cell_roles: usize,
    gene_roles: usize,
    #[serde(serialize_with = "strength::serialize")]
    strength: f64,
    alpha: f64,
    gene_alpha: f64,
    shift_sd: Option<f64>,
    gamma: Option<f64>,
    cell_sd: Option<f64>,
    svg: String,
}

/// The section for the run `run` (see [`super::super::Modules::run`]).
pub(super) fn section(run: &Result<Recorded, String>) -> Section {
    match run {
        Ok(r) => {
            let s = &r.settings;
            Section::Drawn(Box::new(Drawn {
                preset: r.preset.clone(),
                cell_roles: s.cell_roles,
                gene_roles: s.gene_roles,
                strength: s.strength,
                alpha: s.alpha,
                gene_alpha: s.gene_alpha,
                shift_sd: s.shift_sd,
                gamma: s.weights.gamma(),
                cell_sd: s.strength.is_infinite().then_some(s.cell_sd),
                svg: svg(s),
            }))
        }
        Err(why) => Section::Omitted { why: why.clone() },
    }
}

/// The diagram's parts, each a group named by `data-node`, greyed (and
/// marked `data-fixed`) when the run fixes it.
#[derive(Default)]
struct Parts(String);

impl Parts {
    /// Part `name`, drawn by `body`.
    fn add(&mut self, name: &str, fixed: bool, body: &str) {
        let grey = if fixed {
            r#" data-fixed="" opacity="0.45""#
        } else {
            ""
        };
        let _ = write!(self.0, r#"<g data-node="{name}"{grey}>{body}</g>"#);
    }

    /// An edge `name` along the path `d`.
    fn edge(&mut self, name: &str, fixed: bool, d: &str) {
        self.add(name, fixed, &format!(r#"<path class="ed" d="{d}"/>"#));
    }

    /// A fixed setting `name`: a dot at `(x, y)`, `label` to its right
    /// (or its left, at the diagram's edge).
    fn dot(&mut self, name: &str, fixed: bool, (x, y): (u32, u32), label: &str) {
        let (lx, anchor) = if x > 1000 {
            (x - 12, "end")
        } else {
            (x + 12, "start")
        };
        let body = format!(
            r#"<circle class="dot" cx="{x}" cy="{y}" r="5"/><text class="hv" x="{lx}" y="{}" text-anchor="{anchor}">{label}</text>"#,
            y - 4
        );
        self.add(name, fixed, &body);
    }

    /// A variable `name` at `(x, y)` of radius `r`: sampled, or `collapsed`
    /// (integrated out, dashed), `label` in it.
    fn node(
        &mut self,
        name: &str,
        fixed: bool,
        (x, y, r): (u32, u32, u32),
        collapsed: bool,
        label: &str,
    ) {
        let class = if collapsed { "nd cl" } else { "nd" };
        let body = format!(
            r#"<circle class="{class}" cx="{x}" cy="{y}" r="{r}"/><text class="m" x="{x}" y="{}">{label}</text>"#,
            y + 7
        );
        self.add(name, fixed, &body);
    }
}

/// `base` with subscript `sub`, in a node's label.
fn sub(base: &str, sub: &str) -> String {
    format!(r#"{base}<tspan class="sub" dy="5">{sub}</tspan>"#)
}

/// `base` with subscript `sub`, then `after` back on the line, in running
/// text (a setting's label, a formula).
fn lowered(base: &str, sub: &str, after: &str) -> String {
    format!(r#"{base}<tspan class="sub" dy="4">{sub}</tspan><tspan dy="-4">{after}</tspan>"#)
}

/// A setting's label: `base` with subscript `sub`, ` = x`.
fn setting(base: &str, sub: &str, x: f64) -> String {
    lowered(base, sub, &format!(" = {x}"))
}

/// A plate at `(x, y)` of `w` by `h`, `label` in its lower right corner.
fn plate(out: &mut String, (x, y, w, h): (u32, u32, u32, u32), label: &str) {
    let _ = write!(
        out,
        r#"<rect class="pl" x="{x}" y="{y}" width="{w}" height="{h}" rx="6"/><text class="pt" x="{}" y="{}" text-anchor="end">{label}</text>"#,
        x + w - 10,
        y + h - 8
    );
}

/// A switch's flag at `(x, y)`, `w` wide, a line of text each in `lines`.
fn flag(parts: &mut Parts, name: &str, (x, y, w): (u32, u32, u32), lines: &[&str]) {
    let mut body = format!(
        r#"<rect class="fl" x="{x}" y="{y}" width="{w}" height="{}" rx="4"/>"#,
        6 + 16 * lines.len()
    );
    for (i, line) in lines.iter().enumerate() {
        let _ = write!(
            body,
            r#"<text class="ft" x="{}" y="{}">{line}</text>"#,
            x + 10,
            y + 16 + 16 * i as u32
        );
    }
    parts.add(name, false, &body);
}

/// The plate diagram of the blockmodel `s` sets: the shift off draws no
/// `delta`, one gene role greys `psi`, symmetric weights grey `beta`
/// (fixed, no `gamma`), and an infinite strength draws no `p`, the reads
/// drawn from the shifted baseline itself and each cell's own level `nu`.
pub(super) fn svg(s: &Settings) -> String {
    let mut p = Parts::default();
    // Edges first, under the nodes.
    edges(&mut p, s);
    nodes(&mut p, s);
    flags(&mut p, s);
    format!(
        r#"<svg class="plate-svg" viewBox="0 0 1060 640" role="img" aria-label="Plate diagram of this run's blockmodel"><defs><marker id="plate-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="dot" d="M0 0L10 5L0 10z"/></marker></defs>{}{}{}</svg>"#,
        plates(s),
        p.0,
        notes(s)
    )
}

/// The diagram's plates: the blocks' over the shifts with the shift on,
/// and over the profiles at finite strength.
fn plates(s: &Settings) -> String {
    let mut plates = String::new();
    plate(&mut plates, (160, 112, 210, 104), "cells c");
    plate(&mut plates, (500, 112, 210, 104), "genes g");
    plate(&mut plates, (160, 252, 550, 330), "links ℓ = (c, g)");
    plate(&mut plates, (320, 378, 260, 180), "sites s of gene g");
    plate(&mut plates, (770, 404, 270, 200), "sites s");
    let blocks = format!("blocks (a, b): {} × {}", s.cell_roles, s.gene_roles);
    if s.shift_sd.is_some() {
        plate(&mut plates, (770, 252, 250, 110), &blocks);
    }
    if s.strength.is_finite() {
        // Its label clear of `sites s` below it and of `p`'s note above.
        plate(&mut plates, (850, 424, 180, 150), "blocks (a, b)");
    }
    plates
}

/// The diagram's edges, as `s` sets the switches.
fn edges(p: &mut Parts, s: &Settings) {
    let (shift, finite) = (s.shift_sd.is_some(), s.strength.is_finite());
    let (hdp, one_role) = (s.weights.gamma().is_some(), s.gene_roles == 1);
    if hdp {
        p.edge("gamma-beta", false, "M96 66L96 128");
    }
    p.edge("beta-theta", !hdp, "M118 160L236 162");
    p.edge("alpha-theta", false, "M268 70L266 132");
    p.edge("alpha-g-psi", one_role, "M606 70L606 132");
    p.edge("theta-a", false, "M266 190L270 292");
    p.edge("psi-b", one_role, "M606 190L602 292");
    p.edge("a-y", false, "M292 336L428 452");
    p.edge("b-y", false, "M582 336L468 452");
    if shift {
        p.edge("sigma-delta-delta", false, "M900 78L904 282");
    }
    if finite {
        if shift {
            p.edge("delta-p", false, "M926 336L934 470");
        }
        p.edge("mu-p", false, "M836 500L908 500");
        p.edge("strength-p", false, "M1040 216L958 476");
        // Over the baseline, clear of its label and the links' plate label.
        p.edge("p-y", false, "M914 486C800 400 560 420 478 466");
    } else {
        if shift {
            p.edge("delta-y", false, "M898 326C760 380 560 400 472 452");
        }
        p.edge("mu-y", false, "M796 504C700 540 560 500 478 476");
        p.edge("sigma-nu-nu", false, "M404 70L350 144");
        // Right of the cells' plate label and `theta`'s note.
        p.edge("nu-y", false, "M352 178C400 230 420 360 436 442");
    }
}

/// The diagram's settings, variables, and the observed reads.
fn nodes(p: &mut Parts, s: &Settings) {
    let finite = s.strength.is_finite();
    let (hdp, one_role) = (s.weights.gamma(), s.gene_roles == 1);
    if let Some(g) = hdp {
        p.dot("gamma", false, (96, 58), &format!("γ = {g}"));
    }
    p.node("beta", hdp.is_none(), (96, 160, 22), false, "β");
    p.dot("alpha", false, (268, 62), &format!("α = {}", s.alpha));
    p.node("theta", false, (264, 162, 26), true, &sub("θ", "c"));
    if !finite {
        p.dot("sigma-nu", false, (404, 62), &setting("σ", "ν", s.cell_sd));
        p.node("nu", false, (336, 162, 20), false, &sub("ν", "c"));
    }
    p.dot(
        "alpha-g",
        one_role,
        (606, 62),
        &setting("α", "g", s.gene_alpha),
    );
    p.node("psi", one_role, (606, 162, 26), true, &sub("ψ", "g"));
    p.node("a", false, (272, 318, 24), false, &sub("a", "ℓ"));
    p.node("b", false, (600, 318, 24), false, &sub("b", "ℓ"));
    let mut obs = String::from(r#"<circle class="nd obs" cx="448" cy="470" r="30"/>"#);
    let _ = write!(
        obs,
        r#"<text class="m" x="448" y="477">{}</text>"#,
        sub("y", "ℓs")
    );
    p.add("y", false, &obs);
    p.add(
        "mu",
        false,
        &format!(r#"<rect class="nd" x="796" y="480" width="40" height="40" rx="3"/><text class="m" x="816" y="507">{}</text>"#, sub("μ", "s")),
    );
    if let Some(sd) = s.shift_sd {
        p.dot("sigma-delta", false, (900, 70), &setting("σ", "δ", sd));
        p.node("delta", false, (922, 310, 26), false, &sub("δ", "ab"));
    }
    if finite {
        p.dot(
            "strength",
            false,
            (1040, 208),
            &format!("S = {}", s.strength),
        );
        p.node("p", false, (938, 500, 28), true, &sub("p", "s,ab"));
    }
}

/// The switches' flags, as the run set them.
fn flags(p: &mut Parts, s: &Settings) {
    let weights = match s.weights.gamma() {
        Some(_) => "HDP: β drawn",
        None => "symmetric: β = 1/A",
    };
    flag(p, "flag-weights", (16, 224, 140), &[weights]);
    let roles = format!("gene roles B = {}", s.gene_roles);
    flag(p, "flag-roles", (500, 222, 210), &[&roles]);
    let shifted = if s.shift_sd.is_some() {
        "shift on"
    } else {
        "shift off: no block shift"
    };
    flag(p, "flag-shift", (770, 178, 196), &[shifted]);
    let strength = if s.strength.is_finite() {
        format!("strength S = {}", s.strength)
    } else {
        "strength S = ∞: reads Binomial".to_string()
    };
    flag(p, "flag-strength", (770, 140, 230), &[&strength]);
}

/// What each part is, beside it, and in one line under the diagram how
/// the reads are drawn.
fn notes(s: &Settings) -> String {
    let (shift, finite) = (s.shift_sd.is_some(), s.strength.is_finite());
    let mut notes = String::new();
    let mut note = |x: u32, y: u32, anchor: &str, text: &str| {
        let _ = write!(
            notes,
            r#"<text class="an" x="{x}" y="{y}" text-anchor="{anchor}">{text}</text>"#
        );
    };
    note(280, 236, "start", "integrated out");
    if s.weights.gamma().is_none() {
        note(96, 200, "middle", "fixed");
    }
    note(
        616,
        106,
        "start",
        if s.gene_roles == 1 {
            "one role"
        } else {
            "integrated out"
        },
    );
    note(272, 360, "middle", &format!("a ∈ 1…{}", s.cell_roles));
    note(600, 360, "middle", &format!("b ∈ 1…{}", s.gene_roles));
    note(816, 540, "middle", "pooled");
    note(816, 554, "middle", "baseline");
    note(448, 522, "middle", "m converted,");
    note(448, 536, "middle", "u unconverted");
    if finite {
        note(938, 548, "middle", "integrated out");
    }
    let centre = if !finite {
        format!(
            "σ({}{}{}",
            lowered("μ", "s", " + ν"),
            lowered("", "c", " + δ"),
            lowered("", "ab", ")")
        )
    } else if shift {
        format!("σ({}{}", lowered("μ", "s", " + δ"), lowered("", "ab", ")"))
    } else {
        format!("σ({}", lowered("μ", "s", ")"))
    };
    // Infinite strength comes with the shift (see `unified::Model::new`).
    let said = if finite {
        format!("Each block's reads at a site: Beta-binomial, S pseudo-reads around {centre}.")
    } else {
        format!("S = ∞: each link's reads at a site Binomial at {centre}, no p to integrate.")
    };
    let _ = write!(notes, r#"<text class="say" x="160" y="618">{said}</text>"#);
    notes
}
