//! References to download when no annotation or genome is at hand, chosen
//! from what the GENCODE and Ensembl FTP sites list: a source, then a
//! release or species, whose primary-assembly GTF and genome are found in
//! its folder's listing. Listing and download run on threads; files
//! already in the folder are kept, so an interrupted download picks up
//! where it left.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use ratatui::crossterm::event::{KeyCode, KeyEvent};

use super::child::Stopper;
use crate::tui::{moved, FindList};

const GENCODE: &str = "https://ftp.ebi.ac.uk/pub/databases/gencode";
const ENSEMBL: &str = "https://ftp.ensembl.org/pub/current";

/// Where references come from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Source {
    /// GENCODE, by its folder for one species; chromosomes named `chr1`.
    Gencode(&'static str),
    /// Ensembl's current release, any species; chromosomes named `1`.
    Ensembl,
}

pub const SOURCES: [(&str, Source); 3] = [
    ("human · GENCODE", Source::Gencode("Gencode_human")),
    ("mouse · GENCODE", Source::Gencode("Gencode_mouse")),
    ("any species · Ensembl, current release", Source::Ensembl),
];

impl Source {
    /// The listing to choose from: a GENCODE species' releases, or
    /// Ensembl's species.
    fn catalogue_url(self) -> String {
        match self {
            Source::Gencode(d) => format!("{GENCODE}/{d}/"),
            Source::Ensembl => format!("{ENSEMBL}/gtf/"),
        }
    }

    /// Whether the source names chromosomes `chr1` rather than `1`.
    pub fn chr_named(self) -> bool {
        matches!(self, Source::Gencode(_))
    }

    /// What the listing at [`Source::catalogue_url`] offers, best first:
    /// GENCODE releases newest first, Ensembl species by name.
    pub fn items(self, names: Vec<String>) -> Vec<String> {
        match self {
            Source::Gencode(_) => {
                let mut v: Vec<(u32, String)> = names
                    .into_iter()
                    .filter_map(|n| {
                        let num = n.strip_prefix("release_")?.trim_start_matches('M');
                        Some((num.parse().ok()?, n))
                    })
                    .collect();
                v.sort_by_key(|x| std::cmp::Reverse(x.0));
                v.into_iter().map(|(_, n)| n).collect()
            }
            Source::Ensembl => {
                let mut v: Vec<String> = names
                    .into_iter()
                    .filter(|n| n.chars().all(|c| c.is_ascii_lowercase() || c == '_'))
                    .collect();
                v.sort();
                v
            }
        }
    }
}

/// The names a directory listing links to, folders without their `/`;
/// sort links and links up or away are left out.
pub fn hrefs(html: &str) -> Vec<String> {
    html.split("href=\"")
        .skip(1)
        .filter_map(|s| s.split('"').next())
        .filter(|h| !h.is_empty() && !h.starts_with(['?', '/']) && !h.contains("://"))
        .map(|h| h.trim_end_matches('/').to_string())
        .collect()
}

/// Run curl with `args`, where `stopper` reaches it, calling `each` as it
/// waits; what it wrote to stdout, or why it failed.
fn curl(
    args: &[&std::ffi::OsStr],
    stopper: &Stopper,
    each: impl FnMut(),
) -> anyhow::Result<String> {
    let mut child = Command::new("curl")
        .args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| anyhow::anyhow!("cannot run curl: {e}"))?;
    // Read both pipes on threads, so neither can fill and stall curl.
    let out = child.stdout.take().map(read_on_thread);
    let err = child.stderr.take().map(read_on_thread);
    let status = stopper.wait(child, each)?;
    let text = |h: Option<JoinHandle<String>>| h.and_then(|h| h.join().ok()).unwrap_or_default();
    let (out, err) = (text(out), text(err));
    anyhow::ensure!(!stopper.is_stopped(), "stopped");
    anyhow::ensure!(status.success(), "{}", err.trim());
    Ok(out)
}

/// A listing's names, read with curl; `stopper` stops it.
pub fn list(url: &str, stopper: &Stopper) -> anyhow::Result<Vec<String>> {
    let args = ["-fsSL", "--retry", "2", "-m", "60", url].map(std::ffi::OsStr::new);
    let out = curl(&args, stopper, || {}).map_err(|e| anyhow::anyhow!("{url}: {e}"))?;
    Ok(hrefs(&out))
}

fn read_on_thread(mut r: impl Read + Send + 'static) -> JoinHandle<String> {
    std::thread::spawn(move || {
        let mut s = String::new();
        let _ = r.read_to_string(&mut s);
        s
    })
}

/// A listing being read on a thread: its names, or why not.
type Listing = JoinHandle<Result<Vec<String>, String>>;

/// The download pop-up: a source, then one of its releases or species,
/// narrowed by typing.
pub struct Catalogue {
    /// The source chosen, once one is.
    pub source: Option<Source>,
    /// The sources' cursor, until one is chosen.
    pub at: usize,
    /// The source's releases or species.
    pub items: FindList<String>,
    /// The listing on its way, and what stops it.
    pending: Option<(Listing, Arc<Stopper>)>,
    pub error: Option<String>,
    /// Whether the first BAM picked names chromosomes `chr1`, if known.
    pub bam_chr: Option<bool>,
}

impl Catalogue {
    pub fn new(bam_chr: Option<bool>) -> Self {
        // Start on the source whose chromosome names match the BAMs.
        let at = bam_chr.map_or(0, |chr| {
            SOURCES
                .iter()
                .position(|(_, s)| s.chr_named() == chr)
                .unwrap_or(0)
        });
        Catalogue {
            source: None,
            at,
            items: FindList::default(),
            pending: None,
            error: None,
            bam_chr,
        }
    }

    pub fn loading(&self) -> bool {
        self.pending.is_some()
    }

    /// A listing still on its way is not wanted: stop its curl.
    fn stop_listing(&mut self) {
        if let Some((_, stopper)) = self.pending.take() {
            stopper.stop();
        }
    }

    /// Take `key` if it moves the cursor or, among the items, narrows them.
    pub fn key(&mut self, key: &KeyEvent) -> bool {
        if self.source.is_some() {
            return self.items.key(key);
        }
        let d = match key.code {
            KeyCode::Up => -1,
            KeyCode::Down => 1,
            _ => return false,
        };
        self.at = moved(self.at, d, SOURCES.len() - 1);
        true
    }

    /// Enter: a source starts its listing; a release or species is the
    /// pick.
    pub fn enter(&mut self) -> Option<Pick> {
        let Some(source) = self.source else {
            let (_, s) = *SOURCES.get(self.at)?;
            self.source = Some(s);
            let stopper = Arc::new(Stopper::default());
            let stops = stopper.clone();
            let handle = std::thread::spawn(move || {
                list(&s.catalogue_url(), &stops)
                    .map(|names| s.items(names))
                    .map_err(|e| format!("{e:#}"))
            });
            self.pending = Some((handle, stopper));
            return None;
        };
        let item = self.items.highlighted()?.clone();
        Some(Pick { source, item })
    }

    /// Esc: back to the sources; `false` when there already.
    pub fn back(&mut self) -> bool {
        let Some(s) = self.source.take() else {
            return false;
        };
        self.at = SOURCES.iter().position(|(_, x)| *x == s).unwrap_or(0);
        self.items = FindList::default();
        self.stop_listing();
        self.error = None;
        true
    }

    /// Take the listing once it has arrived; whether it just did.
    pub fn poll(&mut self) -> bool {
        if !self.pending.as_ref().is_some_and(|(h, _)| h.is_finished()) {
            return false;
        }
        let got = self.pending.take().map(|(h, _)| h.join());
        match got {
            Some(Ok(Ok(items))) => self.items.set(items, 0),
            Some(Ok(Err(e))) => self.error = Some(e),
            _ => self.error = Some("the listing failed".into()),
        }
        true
    }
}

impl Drop for Catalogue {
    fn drop(&mut self) {
        self.stop_listing();
    }
}

/// One reference chosen: a source and one of its releases or species.
#[derive(Clone, Debug)]
pub struct Pick {
    pub source: Source,
    pub item: String,
}

impl Pick {
    pub fn label(&self) -> String {
        match self.source {
            Source::Gencode(d) => format!("{} {}", d.replace('_', " "), self.item),
            Source::Ensembl => format!("Ensembl {}", self.item),
        }
    }

    /// The folder name suggested for it.
    pub fn dir_name(&self) -> String {
        match self.source {
            Source::Gencode(d) => format!("{}_{}", d.to_lowercase(), self.item),
            Source::Ensembl => format!("ensembl_{}", self.item),
        }
    }

    /// The URLs of the annotation and the genome, found in the listings.
    fn resolve(&self, stopper: &Stopper) -> anyhow::Result<(String, String)> {
        match self.source {
            Source::Gencode(d) => {
                let url = format!("{GENCODE}/{d}/{}", self.item);
                let names = list(&format!("{url}/"), stopper)?;
                let (gtf, genome) = gencode_files(&names).ok_or_else(|| {
                    anyhow::anyhow!("no primary-assembly GTF and genome in {url}")
                })?;
                Ok((format!("{url}/{gtf}"), format!("{url}/{genome}")))
            }
            Source::Ensembl => {
                let gtf_dir = format!("{ENSEMBL}/gtf/{}", self.item);
                let gtf = ensembl_gtf(&list(&format!("{gtf_dir}/"), stopper)?)
                    .ok_or_else(|| anyhow::anyhow!("no GTF in {gtf_dir}"))?;
                let dna = format!("{ENSEMBL}/fasta/{}/dna", self.item);
                let genome = ensembl_genome(&list(&format!("{dna}/"), stopper)?)
                    .ok_or_else(|| anyhow::anyhow!("no genome in {dna}"))?;
                Ok((format!("{gtf_dir}/{gtf}"), format!("{dna}/{genome}")))
            }
        }
    }
}

/// A GENCODE release's primary-assembly annotation and genome.
pub fn gencode_files(names: &[String]) -> Option<(String, String)> {
    let gtf = names
        .iter()
        .find(|n| n.ends_with(".primary_assembly.annotation.gtf.gz"))?;
    let genome = names
        .iter()
        .find(|n| n.ends_with(".primary_assembly.genome.fa.gz"))?;
    Some((gtf.clone(), genome.clone()))
}

/// An Ensembl species' whole annotation, not one cut to chromosomes,
/// patches or ab initio models.
pub fn ensembl_gtf(names: &[String]) -> Option<String> {
    names
        .iter()
        .find(|n| {
            n.ends_with(".gtf.gz")
                && !n.contains(".chr")
                && !n.contains("abinitio")
                && !n.contains("patch")
        })
        .cloned()
}

/// An Ensembl species' unmasked genome: the primary assembly where there is
/// one, else the top level.
pub fn ensembl_genome(names: &[String]) -> Option<String> {
    ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
        .iter()
        .find_map(|tail| names.iter().find(|n| n.ends_with(&format!(".{tail}"))))
        .cloned()
}

fn file_name(url: &str) -> &str {
    url.rsplit('/').next().unwrap_or(url)
}

/// Where a download has got to.
#[derive(Default)]
pub struct Status {
    pub step: String,
    /// The step's bytes done, and of how many where known.
    pub bytes: Option<(u64, Option<u64>)>,
    /// Bumped with each change, so a view redraws only on one.
    pub seq: u64,
    /// The annotation and genome once ready, or why not.
    pub done: Option<Result<(PathBuf, PathBuf), String>>,
}

/// A download under way.
pub struct Fetch {
    pub label: String,
    /// The folder it goes into.
    pub dir: PathBuf,
    pub status: Arc<Mutex<Status>>,
    stopper: Arc<Stopper>,
}

impl Fetch {
    /// Find `pick`'s files and download them into `dir`, on a thread.
    pub fn start(pick: Pick, dir: PathBuf) -> Fetch {
        let status = Arc::new(Mutex::new(Status::default()));
        let stopper = Arc::new(Stopper::default());
        let label = pick.label();
        let into = dir.clone();
        {
            let (status, stopper) = (status.clone(), stopper.clone());
            std::thread::spawn(move || {
                let said = |step: &str, bytes| {
                    if let Ok(mut st) = status.lock() {
                        st.step = step.into();
                        st.bytes = bytes;
                        st.seq += 1;
                    }
                };
                said("finding the files", None);
                let out = pick
                    .resolve(&stopper)
                    .and_then(|(gtf, genome)| fetch(&gtf, &genome, &dir, &stopper, &said))
                    .map_err(|e| format!("{e:#}"));
                if let Ok(mut st) = status.lock() {
                    st.done = Some(out);
                    st.seq += 1;
                }
            });
        }
        Fetch {
            label,
            dir: into,
            status,
            stopper,
        }
    }

    /// A download with no thread behind it, its status set by hand.
    #[cfg(test)]
    pub fn idle(label: &str, dir: PathBuf) -> Fetch {
        Fetch {
            label: label.into(),
            dir,
            status: Arc::default(),
            stopper: Arc::default(),
        }
    }

    pub fn running(&self) -> bool {
        self.status.lock().is_ok_and(|s| s.done.is_none())
    }

    /// Stop the download; what is half fetched stays as `.part` for the
    /// next try. The thread is not waited for: unpacking and indexing end
    /// at their next check, and the rest goes with the process.
    pub fn stop(&self) {
        self.stopper.stop();
    }
}

/// A download given up (the view closed, or another took its place) stops.
impl Drop for Fetch {
    fn drop(&mut self) {
        self.stop();
    }
}

/// The annotation from `gtf_url`, kept gzipped as faba reads it, and the
/// genome from `genome_url`, unpacked and indexed, both in `dir`.
fn fetch(
    gtf_url: &str,
    genome_url: &str,
    dir: &Path,
    stopper: &Stopper,
    said: &Said<'_>,
) -> anyhow::Result<(PathBuf, PathBuf)> {
    std::fs::create_dir_all(dir).map_err(|e| anyhow::anyhow!("{}: {e}", dir.display()))?;
    let gtf = dir.join(file_name(gtf_url));
    if !gtf.exists() {
        download(gtf_url, &gtf, "annotation", stopper, said)?;
    }
    let gz = dir.join(file_name(genome_url));
    let genome = dir.join(file_name(genome_url).trim_end_matches(".gz"));
    if !genome.exists() {
        if gz == genome {
            // Not gzipped: fetched as it is.
            download(genome_url, &genome, "genome", stopper, said)?;
        } else {
            if !gz.exists() {
                download(genome_url, &gz, "genome", stopper, said)?;
            }
            gunzip(&gz, &genome, stopper, said)?;
            let _ = std::fs::remove_file(&gz);
        }
    }
    let fai = PathBuf::from(format!("{}.fai", genome.display()));
    if !fai.exists() {
        said("genome: indexing", None);
        crate::data::util_htslib::load_fasta_index(&genome.to_string_lossy())?;
    }
    Ok((gtf, genome))
}

/// Fetch `url` into `dest` with curl, by way of `dest.part`, saying how
/// far it has got as it grows.
fn download(
    url: &str,
    dest: &Path,
    what: &str,
    stopper: &Stopper,
    said: &Said<'_>,
) -> anyhow::Result<()> {
    let part = PathBuf::from(format!("{}.part", dest.display()));
    let total = content_length(url);
    let mut last = None;
    let args = [
        "-fsSL".as_ref(),
        "--retry".as_ref(),
        "3".as_ref(),
        "-C".as_ref(),
        "-".as_ref(),
        "-o".as_ref(),
        part.as_os_str(),
        url.as_ref(),
    ];
    curl(&args, stopper, || {
        let got = std::fs::metadata(&part).map_or(0, |m| m.len());
        if last.replace(got >> 20) != Some(got >> 20) {
            said(what, Some((got, total)));
        }
    })
    .map_err(|e| anyhow::anyhow!("{url}: {e}"))?;
    std::fs::rename(&part, dest)?;
    Ok(())
}

/// Say a step and how far through its bytes it is.
type Said<'a> = dyn Fn(&str, Option<(u64, Option<u64>)>) + 'a;

/// The size `url` says it has, asked for in a HEAD request.
fn content_length(url: &str) -> Option<u64> {
    let out = Command::new("curl")
        .args(["-fsSIL", "-m", "30", url])
        .stdin(Stdio::null())
        .output()
        .ok()?;
    // The last answer counts: redirects come first.
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .rev()
        .filter_map(|l| {
            let (k, v) = l.split_once(':')?;
            k.trim()
                .eq_ignore_ascii_case("content-length")
                .then(|| v.trim().parse().ok())?
        })
        .next()
}

/// A reader counting the bytes read through it.
struct Counted<R> {
    inner: R,
    n: std::rc::Rc<std::cell::Cell<u64>>,
}

impl<R: Read> Read for Counted<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.n.set(self.n.get() + n as u64);
        Ok(n)
    }
}

/// Unpack `gz` into `dest` by way of `dest.part`, saying how far through
/// `gz` it is.
fn gunzip(gz: &Path, dest: &Path, stopper: &Stopper, said: &Said<'_>) -> anyhow::Result<()> {
    let part = PathBuf::from(format!("{}.part", dest.display()));
    let total = std::fs::metadata(gz)?.len();
    let read = std::rc::Rc::new(std::cell::Cell::new(0));
    let counted = Counted {
        inner: File::open(gz)?,
        n: read.clone(),
    };
    let mut from = flate2::read::MultiGzDecoder::new(BufReader::new(counted));
    let mut to = BufWriter::new(File::create(&part)?);
    let mut buf = vec![0u8; 1 << 20];
    let mut last = None;
    loop {
        anyhow::ensure!(!stopper.is_stopped(), "stopped");
        let n = from.read(&mut buf)?;
        if n == 0 {
            break;
        }
        to.write_all(&buf[..n])?;
        if last.replace(read.get() >> 20) != Some(read.get() >> 20) {
            said("genome: unpacking", Some((read.get(), Some(total))));
        }
    }
    to.flush()?;
    drop(to);
    std::fs::rename(&part, dest)?;
    Ok(())
}

/// Whether `bam` names its chromosomes `chr1` rather than `1`, from its
/// header; `None` when it cannot be read or names none.
pub fn bam_chr_named(bam: &Path) -> Option<bool> {
    use rust_htslib::bam::Read as _;
    let reader = rust_htslib::bam::Reader::from_path(bam).ok()?;
    let names = reader.header().target_names();
    let first = names.first()?;
    Some(first.starts_with(b"chr"))
}

#[cfg(test)]
#[path = "tests/fetch.rs"]
mod tests;
