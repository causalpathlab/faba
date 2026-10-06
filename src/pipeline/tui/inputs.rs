//! The Inputs screen: BAMs picked in a browser and tagged fg (signal) or bg
//! (control), the annotation, genome and known SNPs, and the output.

use std::path::{Path, PathBuf};

/// One spelling for a path: absolute, with `.` and `..` resolved by name.
/// Symlinks are kept as given, so a linked BAM keeps the name it was given
/// and its batch is named after it.
pub fn normalize(path: &Path) -> PathBuf {
    let abs = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    let mut out = PathBuf::new();
    for c in abs.components() {
        match c {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                out.pop();
            }
            c => out.push(c),
        }
    }
    out
}

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

#[derive(Clone, Debug)]
pub struct Entry {
    pub name: String,
    pub dir: bool,
    /// The folder's own path, or the file's.
    pub path: PathBuf,
    /// A file with a BAM index beside it.
    pub indexed: bool,
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

    /// The file names the browser offers.
    pub fn ext(self) -> &'static [&'static str] {
        match self {
            FileRow::Gff => &[".gff", ".gtf", ".gff3", ".gff.gz", ".gtf.gz", ".gff3.gz"],
            FileRow::Genome => &[".fa", ".fasta", ".fa.gz", ".fasta.gz"],
            FileRow::KnownSnps => &[".vcf", ".vcf.gz", ".bcf", ".parquet"],
        }
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

/// A folder's subfolders and files with one of `ext`, `..` first.
pub struct Browser {
    pub cwd: PathBuf,
    pub entries: Vec<Entry>,
    pub at: usize,
    ext: &'static [&'static str],
}

impl Browser {
    pub fn new(cwd: PathBuf, ext: &'static [&'static str]) -> Self {
        let mut b = Browser {
            cwd: cwd.clone(),
            entries: Vec::new(),
            at: 0,
            ext,
        };
        b.open(cwd);
        b
    }

    pub fn open(&mut self, dir: PathBuf) {
        let dir = normalize(&dir);
        let mut dirs = Vec::new();
        let mut files = Vec::new();
        if let Ok(rd) = std::fs::read_dir(&dir) {
            for e in rd.flatten() {
                let name = e.file_name().to_string_lossy().into_owned();
                if name.starts_with('.') {
                    continue;
                }
                let is_dir = e.file_type().is_ok_and(|t| t.is_dir()) || e.path().is_dir();
                if is_dir {
                    dirs.push(name);
                } else if self.ext.iter().any(|x| name.ends_with(x)) {
                    files.push(name);
                }
            }
        }
        dirs.sort();
        files.sort();
        let up = dir.parent().map(|p| Entry {
            name: "..".into(),
            dir: true,
            path: p.to_path_buf(),
            indexed: false,
        });
        let dirs = dirs.into_iter().map(|name| Entry {
            path: dir.join(&name),
            name,
            dir: true,
            indexed: false,
        });
        let files = files.into_iter().map(|name| {
            let path = dir.join(&name);
            Entry {
                indexed: has_index(&path),
                path,
                name,
                dir: false,
            }
        });
        self.entries = up.into_iter().chain(dirs).chain(files).collect();
        self.at = 0;
        self.cwd = dir;
    }

    pub fn up(&mut self) {
        if let Some(p) = self.cwd.parent().map(Path::to_path_buf) {
            self.open(p)
        }
    }

    pub fn step(&mut self, d: isize) {
        let n = self.entries.len() as isize;
        if n > 0 {
            self.at = (self.at as isize + d).clamp(0, n - 1) as usize;
        }
    }

    /// Enter: into the highlighted folder; `Some(file)` when it is a file.
    pub fn enter(&mut self) -> Option<PathBuf> {
        let e = self.entries.get(self.at)?;
        let p = e.path.clone();
        if e.dir {
            self.open(p);
            None
        } else {
            Some(p)
        }
    }
}

pub fn has_index(bam: &Path) -> bool {
    let s = bam.to_string_lossy();
    Path::new(&format!("{s}.bai")).exists() || bam.with_extension("bai").exists()
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
    /// As typed; empty means the suggestion.
    pub output: String,
    pub row: usize,
    pub focus: InputsFocus,
}

impl Inputs {
    pub fn new(cwd: PathBuf) -> Self {
        Inputs {
            bams: Browser::new(cwd, &[".bam"]),
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
        let e = self.bams.entries.get(self.bams.at)?;
        (!e.dir).then(|| e.path.clone())
    }

    /// Set a file row, normalized as picks are.
    pub fn set(&mut self, row: FileRow, path: Option<&Path>) {
        *row.slot(self) = path.map(normalize);
    }

    /// Space: the highlighted BAM in or out of the run; new picks are fg.
    pub fn toggle(&mut self) {
        let Some(p) = self.highlighted_bam() else {
            return;
        };
        match self.picked.iter().position(|x| x.path == p) {
            Some(i) => {
                self.picked.remove(i);
            }
            None => self.picked.push(Picked {
                path: p,
                role: Role::Fg,
            }),
        }
    }

    /// `b`: the highlighted picked BAM between fg and bg.
    pub fn flip(&mut self) {
        let Some(p) = self.highlighted_bam() else {
            return;
        };
        if let Some(x) = self.picked.iter_mut().find(|x| x.path == p) {
            x.role = if x.role == Role::Fg {
                Role::Bg
            } else {
                Role::Fg
            };
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
