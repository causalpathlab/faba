//! The Inputs screen: BAMs picked in a browser and tagged fg (signal) or bg
//! (control), the annotation, genome and known SNPs, and the output.

use std::path::{Path, PathBuf};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Role {
    Fg,
    Bg,
}

#[derive(Clone, Debug)]
pub struct Picked {
    pub path: PathBuf,
    pub role: Role,
}

#[derive(Clone, Debug)]
pub struct Entry {
    pub name: String,
    pub dir: bool,
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
        let up = dir.parent().is_some().then(|| Entry {
            name: "..".into(),
            dir: true,
        });
        self.entries = up
            .into_iter()
            .chain(dirs.into_iter().map(|name| Entry { name, dir: true }))
            .chain(files.into_iter().map(|name| Entry { name, dir: false }))
            .collect();
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

    pub fn highlighted(&self) -> Option<PathBuf> {
        let e = self.entries.get(self.at)?;
        Some(if e.name == ".." {
            self.cwd.parent()?.to_path_buf()
        } else {
            self.cwd.join(&e.name)
        })
    }

    /// Enter: into the highlighted folder; `Some(file)` when it is a file.
    pub fn enter(&mut self) -> Option<PathBuf> {
        let e = self.entries.get(self.at)?.clone();
        let p = self.highlighted()?;
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

pub const ROWS: [&str; 5] = ["GFF", "genome", "known SNPs", "output", "threads"];

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
        (!e.dir).then(|| self.bams.cwd.join(&e.name))
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
        match crate::common::uniq_batch_names(&boxed) {
            Ok(names) if names.len() == paths.len() => paths
                .iter()
                .zip(names)
                .map(|(p, n)| (p.to_path_buf(), n.to_string()))
                .collect(),
            _ => paths
                .iter()
                .map(|p| {
                    let stem = p
                        .file_stem()
                        .map(|s| s.to_string_lossy().into_owned())
                        .unwrap_or_default();
                    (p.to_path_buf(), stem)
                })
                .collect(),
        }
    }

    pub fn argv(&self) -> Vec<String> {
        let s = |p: &Path| {
            std::path::absolute(p)
                .unwrap_or_else(|_| p.to_path_buf())
                .to_string_lossy()
                .into_owned()
        };
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
    /// folder), numbered when that exists.
    pub fn suggested_output(&self) -> String {
        let base = self
            .picked
            .first()
            .and_then(|p| p.path.parent())
            .unwrap_or(&self.bams.cwd);
        let mut out = base.join("faba_out");
        let mut n = 2;
        while out.exists() {
            out = base.join(format!("faba_out{n}"));
            n += 1;
        }
        out.to_string_lossy().into_owned()
    }

    pub fn output(&self) -> String {
        if self.output.trim().is_empty() {
            self.suggested_output()
        } else {
            self.output.trim().to_string()
        }
    }

    pub fn output_problem(&self) -> Option<String> {
        let out = self.output();
        let p = Path::new(&out);
        if p.is_file() {
            return Some(format!("{out} is a file"));
        }
        std::fs::read_dir(p)
            .is_ok_and(|mut d| d.next().is_some())
            .then(|| format!("{out} already contains files; choose an empty one"))
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
