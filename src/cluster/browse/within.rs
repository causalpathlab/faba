//! A within record opened: the count fit, and the site fit within each of
//! its cell modules (see [`crate::cluster::within`]).
//!
//! `faba cluster --within` writes `{m}_within.run.json`, which names the
//! count fit's run record and each group's fit. Opened, it is a list: the
//! count fit first, then a line a group (its cells, what became of its
//! fit). Enter opens the fit on the line in the browser every fit opens in,
//! and leaving that comes back to the list. `--group C3` opens one group's
//! fit as any other fit is opened (its levels as text, its page); without
//! it, the levels asked for are the count fit's, the printout is every
//! fit's modules, and `--html` writes a page a fit, each linking to the
//! others.

use super::report::{self, Link};
use super::{gene_models, text, FitRecord, Modules, ModulesArgs, Opened};
use crate::cluster::command::is_count;
use crate::cluster::within::{self, GroupFit, Status};
use crate::site_analysis::pileup::start_gene_models;
use anyhow::Context;
use data_beans::interactive::ui::{header, help_line, panel, Screen, DIM, HIGHLIGHT, PLAIN};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// A within record: the count fit, and the fit within each of its groups.
#[derive(Debug)]
pub(super) struct Nested {
    /// The record's directory, which the groups' are under.
    dir: PathBuf,
    /// Its file name.
    name: String,
    /// The modality fitted within the groups.
    modality: String,
    /// The count fit's run record.
    count: PathBuf,
    groups: Vec<GroupFit>,
}

impl Nested {
    /// The within record `inputs` name: one named alone, else the one in
    /// the directory they open (`-o`, else the first input's) when it
    /// holds no block tables of that modality on their own, which a
    /// directory opens first. Any other input is handed back as the fit's
    /// run record it names, if it names one: that opens the fit, the count
    /// fit's among them.
    pub(super) fn find(
        inputs: &[Box<str>],
        modality: Option<&str>,
        output: Option<&str>,
    ) -> anyhow::Result<Result<Self, Option<FitRecord>>> {
        let mut named = Vec::new();
        for p in inputs {
            named.extend(Self::read(Path::new(p.as_ref()))?);
        }
        anyhow::ensure!(
            named.is_empty() || inputs.len() == 1,
            "name one within record alone, not with other inputs"
        );
        if let Some(nested) = named.pop() {
            anyhow::ensure!(
                output.is_none(),
                "a run record opens its fits as they are (no -o): to fit again, run `faba cluster`"
            );
            if let Some(m) = modality {
                anyhow::ensure!(
                    m.eq_ignore_ascii_case(&nested.modality),
                    "the within record is of {} fits, not --modality {m}",
                    nested.modality
                );
            }
            return Ok(Ok(nested));
        }
        let fit = super::fit_record(inputs)?;
        if fit.is_some() || inputs.is_empty() || is_count(modality) {
            return Ok(Err(fit));
        }
        let dir = crate::cluster::output_dir(output, inputs);
        if super::tables(&dir, modality).is_some() {
            return Ok(Err(None));
        }
        let asked = modality.map(str::to_ascii_lowercase);
        let modalities = match &asked {
            Some(m) => vec![m.as_str()],
            None => crate::qc::layout::SITE_MODALITIES.to_vec(),
        };
        for m in modalities {
            if let Some(nested) = Self::read(&dir.join(within::record_name(m)))? {
                return Ok(Ok(nested));
            }
        }
        Ok(Err(None))
    }

    /// The within record at `path`; none when it is no run record, or one
    /// that lists no groups. The count fit is the one it names, else (the
    /// folder moved) the one beside it.
    fn read(path: &Path) -> anyhow::Result<Option<Self>> {
        let Some((dir, name, json)) = super::record_at(path)? else {
            return Ok(None);
        };
        let Some((modality, count, groups)) = within::listed(&json) else {
            return Ok(None);
        };
        let beside = dir.join(within::COUNT_RECORD);
        let count = if !count.is_file() && beside.is_file() {
            beside
        } else {
            count
        };
        Ok(Some(Self {
            dir,
            name,
            modality,
            count,
            groups,
        }))
    }

    /// Where the record is: a GFF recorded there serves every fit.
    pub(super) fn dir(&self) -> &Path {
        &self.dir
    }

    /// The groups that were fitted.
    pub(super) fn fitted(&self) -> impl Iterator<Item = &GroupFit> {
        self.groups.iter().filter(|f| f.status == Status::Fitted)
    }

    /// The group named `name`: `C3`, `c3` or `3`.
    pub(super) fn group(&self, name: &str) -> anyhow::Result<&GroupFit> {
        let wanted = match name.parse::<usize>() {
            Ok(k) => format!("C{k}"),
            Err(_) => name.to_string(),
        };
        (self.groups.iter())
            .find(|f| f.group.eq_ignore_ascii_case(&wanted))
            .with_context(|| {
                let has: Vec<&str> = self.groups.iter().map(|f| f.group.as_str()).collect();
                format!(
                    "no group {name} in {}: it has {}",
                    self.name,
                    has.join(", ")
                )
            })
    }

    /// The run record of the fit to open: the count fit's, or that within
    /// `group`, refused when the group has none.
    fn record_of(&self, group: Option<&GroupFit>) -> anyhow::Result<PathBuf> {
        match group {
            None => {
                anyhow::ensure!(
                    self.count.is_file(),
                    "the count fit's run record, {}, is not there",
                    self.count.display()
                );
                Ok(self.count.clone())
            }
            Some(fit) => match (&fit.dir, &fit.record) {
                (Some(dir), Some(record)) => Ok(self.dir.join(dir).join(record)),
                _ => anyhow::bail!("{} has no fit to open ({})", fit.group, fit.what()),
            },
        }
    }

    /// The fit to open, as [`super::resolve`] gives one (see
    /// [`Self::record_of`]).
    pub(super) fn open(&self, group: Option<&GroupFit>) -> anyhow::Result<Opened> {
        let record = self.record_of(group)?;
        super::resolve(&[record.to_string_lossy().into()], None, None)
    }

    /// The pages `--html` writes, the count fit's being `html`: that one,
    /// then each fitted group's beside it (see [`group_page`]).
    pub(super) fn pages(&self, html: &Path) -> Vec<(Option<&GroupFit>, PathBuf)> {
        std::iter::once((None, html.to_path_buf()))
            .chain(self.fitted().map(|f| (Some(f), group_page(html, &f.group))))
            .collect()
    }

    /// The modules of the fit [`Self::open`] opens, and where its tables are.
    fn load(&self, group: Option<&GroupFit>) -> anyhow::Result<(PathBuf, Modules)> {
        let (dir, modality, record) = self.open(group)?;
        let modules = Modules::load(&dir, &modality, record)?;
        Ok((dir, modules))
    }

    /// What `args` ask of the fits together (no `--group`): the levels
    /// named, of the count fit; the pages; the printout; else the list.
    pub(super) fn show(&self, args: &ModulesArgs, gtf: Option<Box<str>>) -> anyhow::Result<()> {
        if let Some(path) = &args.html {
            return self.write_html(gtf.as_deref(), Path::new(path.as_ref()));
        }
        if args.shown.any() {
            return super::show(args, self.open(None)?, gtf);
        }
        if crate::tui::batch(args.batch_process, "cluster-view") {
            self.printed().iter().for_each(|line| println!("{line}"));
            return Ok(());
        }
        self.browse(gtf)
    }

    /// The printout: the count fit's modules, a line a group, then the
    /// modules of each group's fit under a heading. A fit that cannot be
    /// read says why in its place.
    pub(super) fn printed(&self) -> Vec<String> {
        let modules = |group: Option<&GroupFit>| match self.load(group) {
            Ok((_, m)) => text::modules(&m),
            Err(e) => vec![format!("cannot be read: {e:#}")],
        };
        let mut out = vec![format!("count fit: {}", self.count.display())];
        out.extend(modules(None));
        out.push(String::new());
        out.push(format!(
            "{} within its cell modules: {} of {} groups fitted ({})",
            self.modality,
            self.fitted().count(),
            self.groups.len(),
            self.dir.join(&self.name).display()
        ));
        out.push("group\tcells\tfit".into());
        out.extend(self.groups.iter().map(GroupFit::line));
        for fit in self.fitted() {
            out.push(String::new());
            out.push(format!("{} within {}", self.modality, fit.group));
            out.extend(modules(Some(fit)));
        }
        out
    }

    /// `--html`: the pages (see [`Self::pages`]), each linking to the
    /// others; gene models from `gtf`.
    fn write_html(&self, gtf: Option<&str>, html: &Path) -> anyhow::Result<()> {
        let all = super::all_models(gtf)?;
        // Every fit read before any page is written: the pages link to
        // each other, so it is all of them or none.
        let mut fits = Vec::new();
        for (group, page) in self.pages(html) {
            let place = match group {
                None => "every cell".to_string(),
                Some(fit) => format!("within {}", fit.group),
            };
            let (_, modules) = self
                .load(group)
                .with_context(|| format!("the fit {place}, for {}", page.display()))?;
            fits.push((group, page, place, modules));
        }
        for (group, page, place, modules) in fits {
            let data = report::data(&modules, &gene_models(&modules, &all))
                .among(&place, self.links(html, group));
            report::save(&data, &page)?;
            log::info!("wrote {}", page.display());
        }
        Ok(())
    }

    /// The links on the page of `here` (none: the count fit's), the count
    /// fit's page being `html`: the count fit, then every group, those
    /// without a fit named and not linked.
    fn links(&self, html: &Path, here: Option<&GroupFit>) -> Vec<Link> {
        let mut out = vec![Link {
            label: "count".into(),
            href: here.is_some().then(|| href(html)),
            here: here.is_none(),
            note: "the count fit: every cell".into(),
        }];
        out.extend(self.groups.iter().map(|fit| {
            let this = here.is_some_and(|h| h.group == fit.group);
            let linked = fit.status == Status::Fitted && !this;
            Link {
                label: fit.group.clone(),
                href: linked.then(|| href(&group_page(html, &fit.group))),
                here: this,
                note: format!("{} cells, {}", fit.cells, fit.what()),
            }
        }));
        out
    }

    /// The list, until it is quit: Enter opens the fit on the line in the
    /// browser, and leaving that comes back to the list as it was. The
    /// gene annotation is asked for once, when the first fit opens, and
    /// read once for every fit opened.
    fn browse(&self, gtf: Option<Box<str>>) -> anyhow::Result<()> {
        let mut list = List::new(self);
        let mut gtf = gtf;
        let mut models = None;
        loop {
            crate::tui::run_view(&mut list)?;
            let Some(row) = list.take_open() else {
                return Ok(());
            };
            let (dir, modules) = match self.load(self.row(row)) {
                Ok(loaded) => loaded,
                Err(e) => {
                    list.status = Some(format!("{e:#}"));
                    continue;
                }
            };
            if models.is_none() {
                let gtf = super::or_asked(gtf.take())?;
                models = Some(start_gene_models(gtf.as_deref()));
            }
            super::browse(&dir, &Arc::new(modules), models.clone().flatten())?;
        }
    }

    /// The group on the list's line `row`; none on the first, the count fit's.
    fn row(&self, row: usize) -> Option<&GroupFit> {
        row.checked_sub(1).and_then(|g| self.groups.get(g))
    }
}

/// The page of the fit within `group`, beside the count fit's page `html`:
/// `{its stem}.{group}.html`.
fn group_page(html: &Path, group: &str) -> PathBuf {
    let stem = html.file_stem().unwrap_or_default().to_string_lossy();
    html.with_file_name(format!("{stem}.{group}.html"))
}

/// The link to the page `path` from one beside it: its file name, the
/// characters a URL reads otherwise escaped.
fn href(path: &Path) -> String {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let mut out = String::with_capacity(name.len());
    for c in name.chars() {
        match c {
            '%' | '#' | '?' | ' ' => out.push_str(&format!("%{:02X}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// The list of a within record's fits: the count fit, then each group.
struct List<'a> {
    nested: &'a Nested,
    /// The line in focus: 0 the count fit, then the groups.
    at: usize,
    /// The line Enter opened.
    open: Option<usize>,
    done: bool,
    status: Option<String>,
}

impl<'a> List<'a> {
    fn new(nested: &'a Nested) -> Self {
        Self {
            nested,
            at: 0,
            open: None,
            done: false,
            status: None,
        }
    }

    /// The line to open, when Enter left the list; none when it was quit.
    /// The list is then ready to be shown again, as it was.
    fn take_open(&mut self) -> Option<usize> {
        self.done = false;
        self.open.take()
    }

    /// Open the fit on the line in focus; a group without one says why.
    fn enter(&mut self) {
        match self.nested.record_of(self.nested.row(self.at)) {
            Err(e) => self.status = Some(format!("{e:#}")),
            Ok(_) => {
                self.open = Some(self.at);
                self.done = true;
            }
        }
    }

    /// The list's lines: the count fit, then each group with its cells
    /// and what became of its fit, those without a fit dim.
    fn lines(&self) -> Vec<Line<'static>> {
        let line = |row: usize, name: &str, what: String, fitted: bool| {
            let (name_style, what_style) = if fitted {
                (HIGHLIGHT, PLAIN)
            } else {
                (DIM, DIM)
            };
            Line::from(vec![
                crate::tui::marker(row == self.at),
                Span::styled(format!("{name:<7}"), name_style),
                Span::styled(what, what_style),
            ])
        };
        let mut out = vec![line(
            0,
            "count",
            "the count fit: gene modules × cell modules over every cell".into(),
            true,
        )];
        out.extend(self.nested.groups.iter().enumerate().map(|(g, fit)| {
            let what = format!("{:>7} cells   {}", fit.cells, fit.what());
            line(g + 1, &fit.group, what, fit.status == Status::Fitted)
        }));
        out
    }
}

impl crate::tui::View for List<'_> {
    /// The apply key opens the fit, as Enter does.
    fn takes_apply(&self) -> bool {
        true
    }
}

impl Screen for List<'_> {
    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.status = None;
        let (at, last) = (self.at, self.nested.groups.len());
        let moved = |d: isize| crate::tui::moved(at, d, last);
        match key.code {
            _ if crate::tui::is_apply(&key) => self.enter(),
            KeyCode::Up | KeyCode::Char('k') => self.at = moved(-1),
            KeyCode::Down | KeyCode::Char('j') => self.at = moved(1),
            KeyCode::PageUp => self.at = moved(-10),
            KeyCode::PageDown => self.at = moved(10),
            KeyCode::Home => self.at = 0,
            KeyCode::End => self.at = last,
            KeyCode::Enter | KeyCode::Right => self.enter(),
            KeyCode::Char('q') | KeyCode::Esc => self.done = true,
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let [top, body, foot] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(1),
            Constraint::Length(1),
        ])
        .areas(frame.area());
        let n = self.nested;
        let title = n.dir.join(&n.name).display().to_string();
        frame.render_widget(header("cluster-view", &title, ""), top);
        let block = panel(
            format!(
                " the count fit, and {} within its cell modules: {} of {} groups fitted ",
                n.modality,
                n.fitted().count(),
                n.groups.len()
            ),
            true,
        );
        let inner = block.inner(body);
        frame.render_widget(block, body);
        let lines = crate::tui::window(self.lines(), self.at, usize::from(inner.height));
        frame.render_widget(Paragraph::new(lines), inner);
        let footer = match &self.status {
            Some(msg) => crate::figure::status_line(msg),
            None => help_line(&[("↑/↓", "move"), ("Enter", "open the fit"), ("q", "quit")]),
        };
        frame.render_widget(Paragraph::new(footer), foot);
    }
}

#[cfg(test)]
mod tests;
