//! The README's picture of the run view.

use super::inputs::InputsFocus;
use super::*;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use ratatui::Terminal;

/// The run view on made-up files: three samples as signal and two
/// controls, their indexes beside them, an annotation and a genome; drawn
/// as if in `~/demo`, so no real folder shows. Written where
/// `FABA_SCREENSHOTS_DIR` names.
#[test]
#[ignore = "writes the README's screenshot; run by hand"]
fn readme_screenshots_run() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    for sub in ["reference", "logs"] {
        std::fs::create_dir(dir.join(sub)).unwrap();
    }
    let bams = ["ctrl1", "ctrl2", "s1", "s2", "s3"];
    for name in bams {
        for ext in ["bam", "bam.bai"] {
            std::fs::write(dir.join(format!("{name}.{ext}")), b"").unwrap();
        }
    }
    std::fs::write(dir.join("unsorted.bam"), b"").unwrap();
    let mut a = App::new(run_cmd(), dir.to_path_buf());
    // Space walks a BAM out, to signal, to control.
    for (name, presses) in [("s1", 1), ("s2", 1), ("s3", 1), ("ctrl1", 2), ("ctrl2", 2)] {
        assert!(a.inputs.bams.select(&format!("{name}.bam")));
        for _ in 0..presses {
            crate::tui::feed(
                &mut a,
                KeyEvent::new(KeyCode::Char(' '), KeyModifiers::NONE),
            );
        }
    }
    assert!(a.inputs.bams.select("s1.bam"));
    a.inputs.output = dir.join("out").to_string_lossy().into_owned();
    a.refresh();
    // The same files as if under `~/demo`: only the names drawn change.
    let home = crate::tui::home().expect("a home folder");
    let shown = |p: &std::path::Path| home.join("demo").join(p.strip_prefix(dir).unwrap());
    a.inputs.bams.dir = shown(&a.inputs.bams.dir);
    for p in &mut a.inputs.picked {
        p.path = shown(&p.path);
    }
    a.inputs.gff = Some(home.join("demo/reference/genes.gtf"));
    a.inputs.genome = Some(home.join("demo/reference/genome.fa"));
    a.inputs.output = shown(std::path::Path::new(&a.inputs.output))
        .to_string_lossy()
        .into_owned();
    a.inputs.focus = InputsFocus::Bams;
    let mut term = Terminal::new(TestBackend::new(160, 32)).unwrap();
    term.draw(|f| a.render(f)).unwrap();
    crate::tui::screenshot::save("run", term.backend().buffer());
}
