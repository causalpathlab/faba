//! Loading the tables as modules, and the numbers the grids draw from.

use super::{count, counted, loaded};
use crate::cluster::browse::data::{Count, Measure};
use crate::cluster::browse::text::modules;
use crate::cluster::browse::view::{box_size, eighths, hanging, Ramp, Stretch};
use crate::cluster::browse::{diagonal, standardized};
use ratatui::style::Color;

#[test]
fn the_tables_load_as_modules() {
    let m = loaded();
    assert_eq!((m.cell_roles, m.gene_roles), (2, 3));
    // Fitted cells' memberships summed: B2 was placed.
    let sizes = [0.8 + 0.7 + 0.1, 0.2 + 0.3 + 0.9];
    for (size, expect) in m.cell_sizes.iter().zip(sizes) {
        assert!((size - expect).abs() < 1e-6, "{:?}", m.cell_sizes);
    }
    assert_eq!(m.gene_sizes, vec![2, 1, 0]);
    let names: Vec<(&str, i64)> = m.sites.iter().map(|s| (s.gene.as_ref(), s.pos)).collect();
    let expect = [("GENE1", 100), ("GENE1", 300), ("GENE2", 50), ("GENE3", 10)];
    assert_eq!(names, expect);
    // The widest shift range first; the gene role without sites left out.
    assert_eq!(m.order, vec![1, 0]);
    assert!((m.shift_range(1) - 1.5).abs() < 1e-6);
    let stats = m.stats.as_ref().unwrap();
    assert_eq!((stats.on[3], stats.off[3], stats.cells[3]), (40, 40, 12));
}

/// UMIs over cells: a rate per cell, bars up alone, printed per cell.
#[test]
fn a_per_cell_count_is_a_rate_with_nothing_below() {
    let per_cell = |on, total| Count {
        on,
        total,
        measure: Measure::PerCell,
    };
    let r = per_cell(7.0, 2.0);
    assert_eq!((r.up(), r.off(), r.fraction()), (3.5, 0.0, Some(3.5)));
    assert_eq!(
        (r.pct(), r.number()),
        ("3.5/cell".into(), Some("3.5".into()))
    );
    assert_eq!(per_cell(30.0, 2.0).number().unwrap(), "15");
    assert_eq!(per_cell(1.0, 4.0).pct(), "0.25/cell");
    let none = per_cell(0.0, 0.0);
    assert_eq!(
        (none.up(), none.pct(), none.number()),
        (0.0, "-".into(), None)
    );
    // Read counts keep their unconverted reads and their %.
    let reads = count(1.0, 4.0);
    assert_eq!(
        (reads.up(), reads.off(), reads.pct()),
        (1.0, 3.0, "25%".into())
    );
    // A sum with UMIs per cell in it stays per cell, whichever comes first.
    let mut sum = Count::default();
    sum.add(per_cell(7.0, 2.0));
    sum.add(per_cell(1.0, 2.0));
    assert_eq!(sum, per_cell(8.0, 4.0));
    let total: Count = [Count::default(), per_cell(3.0, 1.0)].into_iter().sum();
    assert_eq!(total.measure, Measure::PerCell);
}

/// A count fit's tables, as `write::count_genes` leaves them, load as
/// modules of genes, each one site, their reads UMIs per cell.
#[test]
fn a_count_fits_tables_load_as_umis_per_cell() {
    let m = counted();
    assert_eq!(m.modality, "count");
    assert_eq!((m.cell_roles, m.gene_roles), (2, 2));
    assert_eq!(m.cell_sizes, vec![2.0, 2.0]);
    assert_eq!(m.gene_sizes, vec![2, 1]);
    let names: Vec<(&str, i64)> = m.sites.iter().map(|s| (s.gene.as_ref(), s.pos)).collect();
    assert_eq!(names, [("G1", 0), ("G2", 0), ("G3", 0)]);
    assert!(m.stats.is_none() && !m.approximate.pileup);
    // G1's and G2's UMIs over their cells, every gene counted in every
    // cell role: G2 has none in cell role 2 and still counts its 2 cells,
    // so module M1 there is 2 UMIs over 4 gene-cells.
    let m1 = m.block(0);
    assert_eq!(m1[0].on, 14.0);
    assert_eq!(m1[0].total, 4.0);
    assert_eq!((m1[1].on, m1[1].total), (2.0, 4.0));
    assert!(m1.iter().all(|r| r.measure == Measure::PerCell));
    assert_eq!(m1[0].number().unwrap(), "3.5");
    let m2 = m.block(1);
    assert_eq!((m2[0].total, m2[0].number()), (2.0, Some("0.00".into())));
    assert_eq!((m2[1].on, m2[1].up()), (6.0, 3.0));
    let g1 = &m.members[0]
        .iter()
        .find(|s| m.sites[s.site].gene.as_ref() == "G1");
    let g1 = g1.expect("G1 in module 1");
    assert!(g1
        .reads
        .iter()
        .all(|r| r.measure == Measure::PerCell && r.off() == 0.0));
}

/// The modules printed, a row each, tab-separated: a count fit's columns
/// are UMIs per cell and cells, a site fit's % and reads.
#[test]
fn a_count_fit_prints_umis_per_cell() {
    let m = counted();
    let lines = modules(&m);
    let head: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(&head[..3], ["module", "main genes", "range"]);
    assert_eq!(head.len(), 3 + 2 * m.cell_roles, "{head:?}");
    for c in 1..=2 {
        assert!(
            head.contains(&format!("C{c} UMI/cell").as_str()),
            "{head:?}"
        );
        assert!(head.contains(&format!("C{c} cells").as_str()), "{head:?}");
    }
    let row = |name: &str| {
        let line = lines
            .iter()
            .find(|l| l.starts_with(&format!("{name}\t")))
            .unwrap();
        line.split('\t').map(str::to_string).collect::<Vec<_>>()
    };
    assert!(lines[1..]
        .iter()
        .all(|l| l.split('\t').count() == head.len()));
    assert_eq!([&row("M1")[3], &row("M1")[5]], ["3.5", "0.50"]);
    assert_eq!([&row("M2")[3], &row("M2")[5]], ["0.00", "3.0"]);
    // A site fit's are % and reads.
    let head = modules(&loaded()).remove(0);
    assert!(
        head.contains("\tC1 %\t") && head.contains("\tC1 reads"),
        "{head}"
    );
}

#[test]
fn a_missing_or_out_of_range_role_is_an_error_not_a_panic() {
    let roles = arrow::array::Int32Array::from(vec![Some(1), Some(3), None, Some(0), Some(4)]);
    let at = |i| crate::cluster::browse::data::role_at(&roles, i, 3);
    assert_eq!((at(0).unwrap(), at(1).unwrap()), (0, 2));
    for i in 2..5 {
        assert!(at(i).is_err(), "row {i}");
    }
}

#[test]
fn a_module_holds_its_genes_reads_per_cell_module() {
    let m = loaded();
    assert_eq!(m.block(0), vec![count(14.0, 22.0), count(1.0, 10.0)]);
    // GENE2's reads in gene role 2 are not its main role's: left out.
    assert_eq!(m.block(1), vec![count(30.0, 40.0), count(10.0, 40.0)]);
}

#[test]
fn a_module_lists_its_deepest_genes_first_each_by_position() {
    let m = loaded();
    let listed = |g: usize| -> Vec<usize> { m.members[g].iter().map(|s| s.site).collect() };
    // Gene role 1: GENE1 (30 reads) before GENE2 (2).
    assert_eq!(listed(0), vec![0, 1, 2]);
    assert_eq!(listed(1), vec![3]);
    assert_eq!(m.genes(0), &[0..2, 2..3]);
    assert_eq!(m.gene_of(0, 2), "GENE2");
}

#[test]
fn rows_standardize_across_cell_modules() {
    let z = standardized(&[
        count(1.0, 10.0),
        count(3.0, 10.0),
        count(5.0, 10.0),
        count(0.0, 0.0),
    ]);
    let s = 1.5f64.sqrt();
    assert!((z[0].unwrap() + s).abs() < 1e-9 && z[1].unwrap().abs() < 1e-9);
    assert!((z[2].unwrap() - s).abs() < 1e-9);
    assert_eq!(z[3], None, "no reads, no z");
}

#[test]
fn the_diagonal_order_puts_each_rows_peak_after_the_last() {
    // Rows peak at columns 2, 0, 3, 1.
    let peak = |c: usize| {
        (0..4)
            .map(|k| Some(if k == c { 2.0 } else { -0.5 }))
            .collect()
    };
    let z: Vec<Vec<Option<f64>>> = [2, 0, 3, 1].into_iter().map(peak).collect();
    let (rows, cols) = diagonal(&z);
    let place: Vec<usize> = rows
        .iter()
        .map(|&r| cols.iter().position(|&c| c == [2, 0, 3, 1][r]).unwrap())
        .collect();
    assert!(place.windows(2).all(|w| w[0] < w[1]), "{rows:?} {cols:?}");
    let mut sorted = rows.clone();
    sorted.sort();
    assert_eq!(sorted, vec![0, 1, 2, 3], "a permutation");
}

#[test]
fn heatmaps_span_their_whole_colour_range() {
    let rows = [
        vec![count(1.0, 100.0), count(50.0, 100.0)],
        vec![count(10.0, 10.0), count(5.0, 10.0)],
    ];
    let s = Stretch::of(rows.iter().map(Vec::as_slice), Ramp::Linear);
    assert_eq!(s.gray(Ramp::Linear, 0.01), 0, "the lowest is black");
    assert_eq!(s.gray(Ramp::Linear, 1.0), 255, "the highest is white");
    assert_eq!(s.deepest, 100.0);
    assert_eq!(s.diverging(s.z_top).0, Color::Rgb(255, 0, 0));
    assert_eq!(s.diverging(-s.z_top).0, Color::Rgb(0, 0, 255));
    assert_eq!(s.diverging(0.0).0, Color::Rgb(255, 255, 255));
    // log10 runs from 0.1% to 100%.
    assert_eq!(Ramp::Log10.of(1e-3), 0.0);
    assert_eq!(Ramp::Log10.of(1.0), 1.0);
}

#[test]
fn boxes_take_area_in_proportion_to_their_reads() {
    assert_eq!(box_size(100.0, 100.0, 6), (5, 2));
    assert_eq!(box_size(1.0, 100.0, 6), (1, 1));
    let (w, _) = box_size(25.0, 100.0, 6);
    assert!((2..=3).contains(&w), "{w}");
    // A compact cell holds a box two columns wide at most.
    assert_eq!(box_size(100.0, 100.0, 3), (2, 2));
}

#[test]
fn a_bar_fills_eighths_and_never_its_last() {
    assert_eq!(eighths(0.0, 10.0, 1), 0);
    assert_eq!(eighths(10.0, 10.0, 1), 7);
    assert_eq!(eighths(10.0, 10.0, 2), 15);
    assert_eq!(eighths(0.001, 10.0, 1), 1, "any reads show");
    assert_eq!(eighths(5.0, 10.0, 1), 4);
    // Hanging bars: the rising block whose blank part is the bar, reversed.
    assert_eq!(hanging(8).0, "█");
    assert_eq!(hanging(2).0, "▆");
    let reversed = ratatui::style::Modifier::REVERSED;
    assert!(hanging(2).1.add_modifier.contains(reversed));
}
