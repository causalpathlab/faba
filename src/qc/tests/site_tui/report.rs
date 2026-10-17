//! The HTML report: its data against the view's, and `Ctrl+E`.

use super::*;
use crate::test_util::wait_for;

/// The report's chord.
fn ctrl_e() -> KeyEvent {
    KeyEvent::new(KeyCode::Char('e'), KeyModifiers::CONTROL)
}

/// A picker over `t` writing its reports to a folder of its own, which
/// outlives it.
fn report_picker(t: &SiteTable) -> (tempfile::TempDir, SitePicker<'_>) {
    let dir = tempfile::tempdir().unwrap();
    let mut p = picker(t, None, SiteFilterArgs::default_values());
    p.report_dir = dir.path().to_path_buf();
    (dir, p)
}

/// The report's data as the page reads it.
fn report_json(p: &SitePicker) -> serde_json::Value {
    serde_json::to_value(p.report_source().data()).unwrap()
}

/// A value as the page decodes it: whole counts as they are; a real as
/// `Math.fround` reads it, so a number beyond f32's range is an infinity;
/// null is NaN.
fn decoded(v: &serde_json::Value, integer: bool) -> f64 {
    match v.as_f64() {
        None => f64::NAN,
        Some(x) if integer => x,
        Some(x) => x as f32 as f64,
    }
}

/// The report's knob with criterion id `id`.
fn knob(data: &serde_json::Value, id: usize) -> &serde_json::Value {
    data["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|k| k["id"] == id)
        .expect("a knob per criterion")
}

/// What the page recounts from modality `m` of `data`: per site, the bits
/// of the criteria it fails under the thresholds as written.
fn recount(data: &serde_json::Value, m: usize) -> Vec<u8> {
    let md = &data["modalities"][m];
    let n = md["sites"].as_u64().unwrap() as usize;
    let mut fails = vec![0u8; n];
    for (j, id) in md["criteria"].as_array().unwrap().iter().enumerate() {
        let k = knob(data, id.as_u64().unwrap() as usize);
        let integer = k["integer"].as_bool().unwrap();
        let t = decoded(&k["now"], integer);
        let raw = md["raw"][j].as_array().unwrap();
        assert_eq!(raw.len(), n, "a value per site for {k}");
        let max = k["max"].as_bool().unwrap();
        for (i, v) in raw.iter().enumerate() {
            let v = decoded(v, integer);
            if (max && v > t) || (!max && v < t) {
                fails[i] |= 1 << id.as_u64().unwrap();
            }
        }
    }
    fails
}

#[test]
fn the_report_recounts_what_the_view_tallies() {
    let (t, a) = (with_genes(table(M6A, 900)), with_genes(table(ATOI, 500)));
    let start = SiteFilterArgs::default_values();
    let views = vec![
        SiteView::new(&t, Some(cells(t.len())), &start),
        SiteView::new(&a, None, &start),
    ];
    let mut p = SitePicker::new("x", views, start, Meta::Unavailable("none".into()));
    for c in [
        Criterion::MinCoverage,
        Criterion::MaxPv,
        Criterion::MinCells,
    ] {
        focus_on(&mut p, c);
        press(&mut p, KeyCode::Right);
        press(&mut p, KeyCode::Right);
    }
    for m in 0..2 {
        let data = report_json(&p);
        assert_eq!(data["modality"], m);
        let md = &data["modalities"][m];
        // Every criterion the view applies, in its order, with a column.
        let ids: Vec<u64> = p.view().criteria.iter().map(|&c| c as u64).collect();
        let shown: Vec<u64> = md["criteria"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap())
            .collect();
        assert_eq!(shown, ids);
        let fails = recount(&data, m);
        assert_eq!(fails, p.view().fails, "per site");
        let kept = fails.iter().filter(|&&f| f == 0).count();
        assert_eq!(kept, p.tally.kept);
        for &c in &p.view().criteria {
            let n = fails.iter().filter(|&&f| f & c.bit() != 0).count();
            assert_eq!(n, p.tally.fail[c as usize], "fail {c:?}");
        }
        let gene: Vec<usize> = md["gene"]
            .as_array()
            .unwrap()
            .iter()
            .map(|g| g.as_u64().unwrap() as usize)
            .collect();
        let mut kept_genes: Vec<usize> = (0..gene.len())
            .filter(|&i| fails[i] == 0)
            .map(|i| gene[i])
            .collect();
        kept_genes.sort();
        kept_genes.dedup();
        assert_eq!(kept_genes.len(), p.tally.genes);
        // The gene list as the view orders it: most kept first, then most
        // sites, then by name.
        let genes = md["genes"].as_array().unwrap();
        let tallied = |g: usize| {
            let rows = (0..gene.len()).filter(|&i| gene[i] == g);
            let all = rows.clone().count();
            (rows.filter(|&i| fails[i] == 0).count(), all)
        };
        let mut order: Vec<usize> = (0..genes.len()).collect();
        order.sort_by_key(|&g| {
            let (kept, all) = tallied(g);
            let name = genes[g]["name"].as_str().unwrap().to_lowercase();
            (std::cmp::Reverse(kept), std::cmp::Reverse(all), name)
        });
        let listed: Vec<usize> = p.gene_list().iter().map(|&g| g as usize).collect();
        assert_eq!(order, listed);
        for g in 0..genes.len() {
            assert_eq!(tallied(g), p.gene_kept(g), "gene {g}");
        }
        press(&mut p, KeyCode::Char('m'));
    }
}

#[test]
fn infinities_survive_the_reports_json() {
    let t = table(M6A, 600);
    let p = picker(&t, None, SiteFilterArgs::default_values());
    let text = serde_json::to_string(&p.report_source().data()).unwrap();
    let data: serde_json::Value = serde_json::from_str(&text).unwrap();
    let md = &data["modalities"][0];
    let at = |c: Criterion| {
        let ids = md["criteria"].as_array().unwrap();
        md["raw"][ids.iter().position(|v| v == c as usize).unwrap()]
            .as_array()
            .unwrap()
            .clone()
    };
    let fold = at(Criterion::MinFold);
    let infinite = (0..t.len()).filter(|&i| t.fold[i].is_infinite()).count();
    assert!(infinite > 0, "the fixture has clean controls");
    for (i, v) in fold.iter().enumerate() {
        assert_eq!(decoded(v, false), t.fold[i] as f64, "fold at {i}");
    }
    let odds = at(Criterion::MinLogOdds);
    let mut signed = [false; 2];
    for (i, v) in odds.iter().enumerate() {
        let (got, raw) = (decoded(v, false), t.raw_log_odds[i]);
        if raw.is_infinite() {
            assert_eq!(got, raw, "log odds at {i}");
            signed[usize::from(raw > 0.0)] = true;
        } else {
            // Rounded down to an f32, as a lower bound rounds: the same
            // verdict at every f32 threshold.
            assert!(
                got <= raw && (raw - got).abs() <= raw.abs() * 1e-6,
                "at {i}"
            );
        }
    }
    assert_eq!(signed, [true, true], "both infinities in the fixture");
    // An off threshold at an infinity reads back as one.
    let off = knob(&data, Criterion::MinLogOdds as usize);
    assert_eq!(decoded(&off["off"], false), f64::NEG_INFINITY);
}

/// Save `p`'s figure as `prefix` through the save prompt.
fn save_figure(p: &mut SitePicker, prefix: &std::path::Path) {
    press(p, KeyCode::Char('s'));
    for _ in 0..20 {
        press(p, KeyCode::Backspace);
    }
    for ch in prefix.to_str().unwrap().chars() {
        press(p, KeyCode::Char(ch));
    }
    press(p, KeyCode::Enter);
}

#[test]
fn ctrl_e_writes_a_report_past_a_saved_figures_outcome_but_not_its_prompt() {
    let t = with_genes(table(M6A, 300));
    let (dir, mut p) = report_picker(&t);
    save_figure(&mut p, &dir.path().join("view1"));
    assert!(p.controls.footer().is_some(), "the save's outcome shows");
    crate::tui::feed(&mut p, ctrl_e());
    wait_for(|| p.tick());
    let path = dir.path().join("qc_report.html");
    assert!(path.exists(), "written past the outcome");
    let said = p.footer_line().to_string();
    assert!(
        said.contains(&format!("wrote {}", crate::tui::tilde(&path))),
        "{said}"
    );
    // A report finishing behind a save's outcome is said in its place.
    crate::tui::feed(&mut p, ctrl_e());
    save_figure(&mut p, &dir.path().join("view2"));
    wait_for(|| p.report.running().is_none());
    assert!(p.tick(), "said at once");
    let said = p.footer_line().to_string();
    let second = dir.path().join("qc_report2.html");
    assert!(
        said.contains(&format!("wrote {}", crate::tui::tilde(&second))),
        "{said}"
    );
    // Held back while the save prompt is open, typing nothing into it.
    press(&mut p, KeyCode::Char('s'));
    let typed = p.controls.footer().unwrap().to_string();
    crate::tui::feed(&mut p, ctrl_e());
    assert!(p.controls.prompting(), "the prompt stays open");
    assert_eq!(p.controls.footer().unwrap().to_string(), typed);
    assert!(p.report.poll().is_none(), "nothing written");
}

#[test]
fn only_ctrl_e_writes_the_report() {
    let t = with_genes(table(M6A, 300));
    let (dir, mut p) = report_picker(&t);
    let e = |m: KeyModifiers| KeyEvent::new(KeyCode::Char('e'), m);
    for chord in [
        KeyModifiers::CONTROL | KeyModifiers::ALT,
        KeyModifiers::CONTROL | KeyModifiers::SUPER,
        KeyModifiers::CONTROL | KeyModifiers::META,
    ] {
        crate::tui::feed(&mut p, e(chord));
        assert!(matches!(p.mode, Mode::Browse), "{chord:?} types nothing");
        assert!(p.report.poll().is_none() && p.status.is_none(), "{chord:?}");
    }
    // Shift, as a terminal may add for a capital E, is no matter.
    crate::tui::feed(
        &mut p,
        KeyEvent::new(
            KeyCode::Char('E'),
            KeyModifiers::CONTROL | KeyModifiers::SHIFT,
        ),
    );
    wait_for(|| p.tick());
    assert!(dir.path().join("qc_report.html").exists());
}

#[test]
fn a_genes_span_is_its_models_widened_to_its_sites() {
    let model = GeneModel {
        chr: "chr1".into(),
        lo: 1500,
        hi: 2500,
        forward: true,
        exons: vec![(1500, 2500)],
        symbol: "GENE1".into(),
        key: "1_GENE1".into(),
    };
    let sites = [1000, 2000, 2999];
    assert_eq!(gene_span(None, sites.into_iter()), Some((1000, 3000)));
    assert_eq!(
        gene_span(Some(&model), sites.into_iter()),
        Some((1000, 3000))
    );
    assert_eq!(
        gene_span(Some(&model), [2000].into_iter()),
        Some((1500, 2500))
    );
    // Without sites, its model's; with neither, none.
    assert_eq!(
        gene_span(Some(&model), std::iter::empty()),
        Some((1500, 2500))
    );
    assert_eq!(gene_span(None, std::iter::empty()), None);
}

#[test]
fn a_gene_without_sites_spans_its_model_or_nothing_in_view_and_report() {
    let mut t = with_genes(table(M6A, 300));
    // An eighth gene id no site carries.
    t.n_genes = 8;
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let g = 7;
    assert!(p.view().genes.as_ref().unwrap().rows[g].is_empty());
    let select = |p: &mut SitePicker| {
        let at = p.gene_list().iter().position(|&s| s as usize == g).unwrap();
        p.step_gene(-100);
        p.step_gene(at as isize);
        assert_eq!(p.gene(), Some(g));
    };
    select(&mut p);
    assert!(p.gene_profile(40).is_none(), "no span without a model");
    let gene = &report_json(&p)["modalities"][0]["genes"][g];
    assert!(gene["lo"].is_null() && gene["hi"].is_null(), "{gene}");
    // On a model, its span, with no bars.
    let key = p.view().genes.as_ref().unwrap().keys[g].clone();
    p.set_meta(Meta::ready(
        vec![None],
        vec![GeneModel {
            chr: "chr1".into(),
            lo: 5000,
            hi: 6000,
            forward: true,
            exons: vec![(5000, 6000)],
            symbol: "GENE7".into(),
            key,
        }],
    ));
    select(&mut p);
    let profile = p.gene_profile(40).unwrap();
    assert_eq!((profile.edges.min_pos, profile.edges.max_pos), (5000, 6000));
    assert!(profile.all.iter().all(|&n| n == 0));
    let gene = &report_json(&p)["modalities"][0]["genes"][g];
    assert_eq!(
        (gene["lo"].as_i64(), gene["hi"].as_i64()),
        (Some(5000), Some(6000))
    );
}

#[test]
fn a_report_under_way_is_finished_however_the_view_ends() {
    let t = with_genes(table(M6A, 3000));
    let (dir, mut p) = report_picker(&t);
    crate::tui::feed(&mut p, ctrl_e());
    // As an error out of the run drops it, unpolled.
    drop(p);
    let html = std::fs::read_to_string(dir.path().join("qc_report.html")).unwrap();
    assert!(html.ends_with("</html>\n"));
}

#[test]
fn ctrl_e_writes_numbered_reports_where_it_is_told() {
    let t = with_genes(table(M6A, 300));
    let (dir, mut p) = report_picker(&t);
    for name in ["qc_report.html", "qc_report2.html"] {
        crate::tui::feed(&mut p, ctrl_e());
        assert!(matches!(p.mode, Mode::Browse), "not typing");
        wait_for(|| p.tick());
        let path = dir.path().join(name);
        let html = std::fs::read_to_string(&path).expect(name);
        assert!(html.contains(r#"<script id="data" type="application/json">"#));
        let said = p.footer_line().to_string();
        assert!(
            said.contains(&format!("wrote {}", crate::tui::tilde(&path))),
            "{said}"
        );
    }
    // Plain e still types a threshold, and Ctrl+E types nothing into it.
    press(&mut p, KeyCode::Char('e'));
    crate::tui::feed(&mut p, ctrl_e());
    assert!(matches!(&p.mode, Mode::Edit(buf) if buf.is_empty()));
}

#[test]
fn a_report_that_cannot_be_written_says_why() {
    let t = table(M6A, 100);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.report_dir = "/nowhere/at/all".into();
    crate::tui::feed(&mut p, ctrl_e());
    wait_for(|| p.tick());
    let said = p.footer_line().to_string();
    assert!(said.contains("/nowhere/at/all/qc_report.html"), "{said}");
    assert!(!said.contains("wrote"), "{said}");
}

#[test]
fn the_report_page_loads_nothing() {
    let page = crate::qc::site_tui::report::FRAGMENT.to_lowercase();
    for net in [
        "http:",
        "https:",
        "//cdn",
        "src=",
        "<link",
        "url(",
        "fetch(",
        "xmlhttprequest",
        "import ",
    ] {
        assert!(!page.contains(net), "{net}");
    }
}

#[test]
fn the_report_draws_genes_and_the_metagene_from_what_the_view_draws() {
    let t = with_genes(table(M6A, 1500));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let keys = p.view().genes.as_ref().unwrap().keys.clone();
    // A model wider than its gene's sites and one narrower; the rest have
    // none.
    let model = |g: usize, lo, hi, forward| GeneModel {
        chr: "chr1".into(),
        lo,
        hi,
        forward,
        exons: vec![(lo, lo + 100), (hi - 100, hi)],
        symbol: format!("GENE{g}").into(),
        key: keys[g].clone(),
    };
    p.set_meta(Meta::ready(
        vec![Some(meta_layout(t.len()))],
        vec![model(0, 900, 3100, true), model(1, 1500, 2500, false)],
    ));
    focus_on(&mut p, Criterion::MinCoverage);
    press(&mut p, KeyCode::Right);
    press(&mut p, KeyCode::Right);
    for weight in [Weight::Coverage, Weight::Converted, Weight::Sites] {
        p.weight = weight;
        p.refresh_meta();
        let data = report_json(&p);
        let md = &data["modalities"][0];
        // The bar weight per site, absent when the bars count sites.
        match weight {
            Weight::Sites => assert!(md["weight"].is_null()),
            _ => {
                let w = md["weight"].as_array().unwrap();
                assert_eq!(w.len(), t.len());
                for (i, v) in w.iter().enumerate() {
                    assert_eq!(
                        v.as_u64().unwrap(),
                        weight.of(t.raw_columns(), i) as u64,
                        "{i}"
                    );
                }
            }
        }
        // Each gene's span and model as its profile has them.
        let genes = md["genes"].as_array().unwrap();
        for g in 0..genes.len() {
            let at = p.gene_list().iter().position(|&s| s as usize == g).unwrap();
            p.step_gene(-(genes.len() as isize));
            p.step_gene(at as isize);
            assert_eq!(p.gene(), Some(g));
            let profile = p.gene_profile(40).unwrap();
            let gene = &genes[g];
            let span = (gene["lo"].as_i64().unwrap(), gene["hi"].as_i64().unwrap());
            assert_eq!(span, (profile.edges.min_pos, profile.edges.max_pos), "{g}");
            match &profile.model {
                Some(m) => {
                    assert_eq!(gene["chr"], m.chr.to_string());
                    assert_eq!(gene["forward"], m.forward);
                    let exons: Vec<(i64, i64)> =
                        serde_json::from_value(gene["exons"].clone()).unwrap();
                    assert_eq!(exons, m.exons);
                }
                None => assert!(gene["chr"].is_null() && gene["exons"].is_null()),
            }
        }
        // The metagene's bars as the panel draws them.
        let shown = p.meta_counts().unwrap();
        let bars =
            |key: &str| -> Vec<f64> { serde_json::from_value(md["meta"][key].clone()).unwrap() };
        assert_eq!(bars("all"), shown.all);
        assert_eq!(bars("kept"), shown.kept);
        let regions: [usize; 3] = serde_json::from_value(md["meta"]["regions"].clone()).unwrap();
        assert_eq!(regions, shown.regions);
        assert_eq!(md["meta"]["unassigned"], shown.unassigned);
    }
}

#[test]
fn the_reports_values_bin_where_the_view_bins_them() {
    let t = table(M6A, 900);
    let p = picker(&t, Some(cells(t.len())), SiteFilterArgs::default_values());
    let data = report_json(&p);
    let md = &data["modalities"][0];
    let v = p.view();
    let mut nudged = 0;
    for (j, &c) in v.criteria.iter().enumerate() {
        if c.is_integer() {
            continue;
        }
        // The page's value per site: as stored, one f32 on toward the
        // nearest where `near` lists the site.
        let near: Vec<usize> = serde_json::from_value(md["near"][j].clone()).unwrap();
        nudged += near.len();
        for (i, raw) in md["raw"][j].as_array().unwrap().iter().enumerate() {
            let mut got = decoded(raw, false) as f32;
            if near.contains(&i) {
                got = if c.is_max() {
                    got.next_down()
                } else {
                    got.next_up()
                };
            }
            let want = c.raw(&t, i, v.cells(i)) as f32;
            assert!(
                got == want || (got.is_nan() && want.is_nan()),
                "{c:?} at {i}: {got} for {want}"
            );
        }
    }
    assert!(nudged > 0, "the fixture's log odds round both ways");
}

#[test]
fn a_report_written_while_typing_says_so_back_in_the_view() {
    let t = with_genes(table(M6A, 300));
    let (dir, mut p) = report_picker(&t);
    crate::tui::feed(&mut p, ctrl_e());
    press(&mut p, KeyCode::Char('e'));
    wait_for(|| p.report.running().is_none());
    p.tick();
    press(&mut p, KeyCode::Esc);
    assert!(matches!(p.mode, Mode::Browse));
    assert!(p.tick(), "redraws once it says");
    let said = p.footer_line().to_string();
    let path = dir.path().join("qc_report.html");
    assert!(
        said.contains(&format!("wrote {}", crate::tui::tilde(&path))),
        "{said}"
    );
}

#[test]
fn ctrl_e_again_while_writing_keeps_the_first_reports_outcome_in_view() {
    let t = with_genes(table(M6A, 300));
    let (dir, mut p) = report_picker(&t);
    crate::tui::feed(&mut p, ctrl_e());
    crate::tui::feed(&mut p, ctrl_e());
    assert!(p.footer_line().to_string().contains("then again"));
    wait_for(|| p.tick());
    let (first, second) = (
        dir.path().join("qc_report.html"),
        dir.path().join("qc_report2.html"),
    );
    let said = p.footer_line().to_string();
    let want = format!(
        "wrote {}; writing {}",
        crate::tui::tilde(&first),
        crate::tui::tilde(&second)
    );
    assert!(said.contains(&want), "{said}");
    wait_for(|| p.tick());
    let said = p.footer_line().to_string();
    assert!(
        said.contains(&format!("wrote {}", crate::tui::tilde(&second))),
        "{said}"
    );
}
