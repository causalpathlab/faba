//! The Inputs screen: BAMs picked in a browser and tagged fg (signal) or bg
//! (control), the annotation, genome and known SNPs, and the output.

use std::path::{Path, PathBuf};

pub use crate::tui::browser::{is_annotation, normalize, Browser};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Role {
    Fg,
    Bg,
}

#[derive(Clone, Debug)]
pub struct Picked {
    /// Normalized: one spelling per file.
    pub path: PathBuf,
    pub role: Role,
}

impl Picked {
    pub fn new(path: &Path, role: Role) -> Self {
        Picked {
            path: normalize(path),
            role,
        }
    }
}

/// A file row of the Inputs screen, picked in a pop-up browser.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FileRow {
    Gff,
    Genome,
    KnownSnps,
}

impl FileRow {
    pub fn label(self) -> &'static str {
        match self {
            FileRow::Gff => "GFF",
            FileRow::Genome => "genome",
            FileRow::KnownSnps => "known SNPs",
        }
    }

    /// A browser at `start` offering the files this row takes.
    pub fn browser(self, start: PathBuf) -> Browser {
        let keep: fn(&str) -> bool = match self {
            FileRow::Gff => is_annotation,
            FileRow::Genome => |n| {
                [".fa", ".fasta", ".fa.gz", ".fasta.gz"]
                    .iter()
                    .any(|x| n.ends_with(x))
            },
            FileRow::KnownSnps => |n| {
                [".vcf", ".vcf.gz", ".bcf", ".parquet"]
                    .iter()
                    .any(|x| n.ends_with(x))
            },
        };
        Browser::new(start, keep, |_| false).opened()
    }

    /// Whether `d` downloads this row's file (see [`super::fetch`]).
    pub fn downloadable(self) -> bool {
        matches!(self, FileRow::Gff | FileRow::Genome)
    }

    pub fn slot(self, inputs: &mut Inputs) -> &mut Option<PathBuf> {
        match self {
            FileRow::Gff => &mut inputs.gff,
            FileRow::Genome => &mut inputs.genome,
            FileRow::KnownSnps => &mut inputs.known_snps,
        }
    }

    pub fn get(self, inputs: &Inputs) -> Option<&Path> {
        match self {
            FileRow::Gff => inputs.gff.as_deref(),
            FileRow::Genome => inputs.genome.as_deref(),
            FileRow::KnownSnps => inputs.known_snps.as_deref(),
        }
    }
}

/// The rows of the Inputs screen's right panel, in order.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Row {
    File(FileRow),
    Output,
    Threads,
}

impl Row {
    pub const ALL: [Row; 5] = [
        Row::File(FileRow::Gff),
        Row::File(FileRow::Genome),
        Row::File(FileRow::KnownSnps),
        Row::Output,
        Row::Threads,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Row::File(f) => f.label(),
            Row::Output => "output",
            Row::Threads => "threads",
        }
    }
}

pub fn has_index(bam: &Path) -> bool {
    let s = bam.to_string_lossy();
    Path::new(&format!("{s}.bai")).exists() || bam.with_extension("bai").exists()
}

/// What the output row says of the folder named, worked out once per key
/// rather than at each draw.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub enum OutState {
    #[default]
    Unnamed,
    New,
    Empty,
    NotEmpty,
    File,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum InputsFocus {
    Bams,
    Rows,
}

pub struct Inputs {
    pub bams: Browser,
    pub picked: Vec<Picked>,
    pub gff: Option<PathBuf>,
    pub genome: Option<PathBuf>,
    pub known_snps: Option<PathBuf>,
    /// As named; empty until it is.
    pub output: String,
    pub row: usize,
    pub focus: InputsFocus,
}

impl Inputs {
    pub fn new(cwd: PathBuf) -> Self {
        Inputs {
            bams: Browser::new(cwd, |n| n.ends_with(".bam"), has_index).opened(),
            picked: Vec::new(),
            gff: None,
            genome: None,
            known_snps: None,
            output: String::new(),
            row: 0,
            focus: InputsFocus::Bams,
        }
    }

    fn highlighted_bam(&self) -> Option<PathBuf> {
        let e = self.bams.highlighted()?;
        (!e.dir).then(|| e.path.clone())
    }

    /// Set a file row, normalized as picks are.
    pub fn set(&mut self, row: FileRow, path: Option<&Path>) {
        *row.slot(self) = path.map(normalize);
    }

    /// Space: the highlighted BAM from out of the run to fg, fg to bg, and
    /// bg back out.
    pub fn toggle(&mut self) {
        let Some(p) = self.highlighted_bam() else {
            return;
        };
        match self.picked.iter().position(|x| x.path == p) {
            Some(i) if self.picked[i].role == Role::Fg => self.picked[i].role = Role::Bg,
            Some(i) => {
                self.picked.remove(i);
            }
            None => self.picked.push(Picked {
                path: p,
                role: Role::Fg,
            }),
        }
    }

    /// The role of `path`, spelled as [`normalize`] spells it.
    pub fn role_of(&self, path: &Path) -> Option<Role> {
        self.picked.iter().find(|x| x.path == path).map(|x| x.role)
    }

    fn of(&self, role: Role) -> Vec<&Path> {
        self.picked
            .iter()
            .filter(|x| x.role == role)
            .map(|x| x.path.as_path())
            .collect()
    }

    pub fn fg(&self) -> Vec<&Path> {
        self.of(Role::Fg)
    }

    pub fn bg(&self) -> Vec<&Path> {
        self.of(Role::Bg)
    }

    /// The batch name the pipeline gives each picked BAM (fg then bg, in
    /// pick order), named together so they are unique. Falls back to the
    /// file stem if the pipeline cannot name them.
    pub fn batch_names(&self) -> Vec<(PathBuf, String)> {
        let paths: Vec<&Path> = self.fg().into_iter().chain(self.bg()).collect();
        let boxed: Vec<Box<str>> = paths
            .iter()
            .map(|p| p.to_string_lossy().into_owned().into_boxed_str())
            .collect();
        let names: Vec<String> = match crate::common::uniq_batch_names(&boxed) {
            Ok(names) if names.len() == paths.len() => {
                names.iter().map(ToString::to_string).collect()
            }
            _ => paths
                .iter()
                .map(|p| {
                    p.file_stem()
                        .map(|s| s.to_string_lossy().into_owned())
                        .unwrap_or_default()
                })
                .collect(),
        };
        paths.iter().map(|p| p.to_path_buf()).zip(names).collect()
    }

    pub fn argv(&self) -> Vec<String> {
        let s = |p: &Path| p.to_string_lossy().into_owned();
        let mut v: Vec<String> = self.fg().into_iter().map(s).collect();
        for b in self.bg() {
            v.extend(["--control-bam".into(), s(b)]);
        }
        if let Some(g) = &self.gff {
            v.extend(["-g".into(), s(g)]);
        }
        if let Some(f) = &self.genome {
            v.extend(["-f".into(), s(f)]);
        }
        if let Some(k) = &self.known_snps {
            v.extend(["--known-snps".into(), s(k)]);
        }
        v
    }

    /// `faba_out` beside the first picked BAM (else in the browser's
    /// folder), numbered when that exists: what the output line opens with
    /// until a folder is named.
    pub fn suggested_output(&self) -> String {
        let base = self
            .picked
            .first()
            .and_then(|p| p.path.parent())
            .unwrap_or(&self.bams.cwd);
        crate::tui::next_free(base, "faba_out")
            .to_string_lossy()
            .into_owned()
    }

    /// The output folder as named; empty until it is.
    pub fn output(&self) -> String {
        self.output.trim().to_string()
    }

    pub fn out_state(&self) -> OutState {
        let out = self.output();
        let path = Path::new(&out);
        if out.is_empty() {
            OutState::Unnamed
        } else if path.is_file() {
            OutState::File
        } else if std::fs::read_dir(path).is_ok_and(|mut d| d.next().is_some()) {
            OutState::NotEmpty
        } else if path.exists() {
            OutState::Empty
        } else {
            OutState::New
        }
    }

    pub fn output_problem(&self) -> Option<String> {
        let out = self.output();
        if out.is_empty() {
            return Some("no output folder: name one on the output row".into());
        }
        crate::tui::output_problem(&out)
    }

    pub fn problems(&self) -> Vec<String> {
        let mut v = Vec::new();
        if self.fg().is_empty() {
            v.push("no fg BAM: pick one with Space".into());
        }
        if self.gff.is_none() {
            v.push("no GFF".into());
        }
        if self.genome.is_none() {
            v.push("no genome".into());
        }
        v.extend(self.output_problem());
        v
    }
}

#[cfg(test)]
#[path = "tests/inputs.rs"]
mod tests;
