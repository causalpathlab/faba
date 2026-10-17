//! One gene's reads, from the site matrices, by cell module.

use super::loaded;
use crate::cluster::browse::gene;
use crate::cluster::read::Input;
use data_beans::sparse_io::{create_sparse_from_triplets, SparseIoBackend};

#[test]
fn a_gene_is_read_from_the_matrices_by_cell_module() {
    let m = loaded();
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("b_m6a_site.zarr");
    let rows: Vec<Box<str>> = [
        "GENE1/m6a/chr1:100/methylated",
        "GENE1/m6a/chr1:100/unmethylated",
        "GENE1/m6a/chr1:300/methylated",
        "GENE2/m6a/chr1:50/methylated",
    ]
    .map(Box::from)
    .to_vec();
    let cols: Vec<Box<str>> = (0..5).map(|i| format!("B{i}").into()).collect();
    // B0 and B2 in cell role 1, B3 in role 2, B4 without a role.
    let triplets = [
        (0, 0, 2.0),
        (1, 0, 1.0),
        (0, 2, 3.0),
        (1, 3, 4.0),
        (0, 4, 5.0),
        (2, 3, 6.0),
        (3, 0, 9.0),
    ];
    let shape = (rows.len(), cols.len(), triplets.len());
    let file = path.to_string_lossy();
    let mut data =
        create_sparse_from_triplets(&triplets, shape, Some(&file), Some(&SparseIoBackend::Zarr))
            .unwrap();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&cols);
    drop(data);
    let inputs = vec![Input {
        path: file.into(),
        folder: "".into(),
        batch: "b".into(),
        fit: true,
    }];
    let matrices = gene::Matrices::start(inputs, &m);
    let reads = matrices.read(&m, "GENE1").unwrap();
    assert_eq!(reads.chr.as_ref(), "chr1");
    let group = |c: usize| 1 + m.columns.iter().position(|&x| x == c).unwrap();
    // All cells with a cell module: B4 has none, so it is left out, and all
    // cells is the sum of the cell modules.
    assert_eq!(reads.groups[0].0, vec![(100, 5.0), (300, 6.0)]);
    assert_eq!(reads.groups[0].1, vec![(100, 5.0), (300, 0.0)]);
    assert_eq!(reads.groups[group(0)].0, vec![(100, 5.0)]);
    assert_eq!(reads.groups[group(0)].1, vec![(100, 1.0)]);
    // A position with only unconverted reads shows 0 converted there.
    assert_eq!(reads.groups[group(1)].0, vec![(100, 0.0), (300, 6.0)]);
    assert_eq!(reads.groups[group(1)].1, vec![(100, 4.0), (300, 0.0)]);
    assert!(matrices.read(&m, "GENE9").is_err());
}

#[test]
fn a_gene_opened_again_in_place_keeps_its_scroll() {
    use crate::site_analysis::pileup::tui::Exit;
    let mut place = gene::Place::default();
    let exits = [
        Exit::Locus {
            lo: 1,
            hi: 9,
            at: 5,
        },
        Exit::Search("GENE2".into()),
        Exit::Modality { at: 5 },
    ];
    for (k, exit) in exits.into_iter().enumerate() {
        assert_eq!(place.left(exit, 2 * k + 4), None);
        assert_eq!(place.track_from, 2 * k + 4, "scrolled as it was");
    }
    assert_eq!(place.extent, (1, 9));
    assert_eq!(place.left(Exit::Genes, 0), Some(false));
    assert_eq!(place.left(Exit::Quit, 0), Some(true));
}
