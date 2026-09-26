//! The gene list `faba pileup --interactive` opens when no gene or region
//! was given: every gene with site rows in the inputs, read off the row
//! names alone, filtered as you type.

use data_beans::hdf5_io::resolve_backend_file;
use data_beans::interactive::ui::{
    header, input_line, panel, run_screen, Screen, DIM, HIGHLIGHT, PLAIN,
};
use data_beans::sparse_io::open_sparse_matrix;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::{FxHashMap, FxHashSet};

use super::{fmt_thousands, parse_query, parse_row_name_full, Query};

/// One gene of the inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneEntry {
    /// The row key, e.g. `ENSG..._SYMBOL`.
    pub gene: Box<str>,
    pub chr: Box<str>,
    pub lo: i64,
    pub hi: i64,
    /// Distinct site positions over all inputs.
    pub sites: usize,
}

/// Fold row names into per-gene entries, most sites first.
pub fn catalog_from_rows<'a>(rows: impl IntoIterator<Item = &'a str>) -> Vec<GeneEntry> {
    let mut genes: FxHashMap<&str, (&str, FxHashSet<i64>)> = FxHashMap::default();
    for name in rows {
        if let Some((gene, _, chr, pos)) = parse_row_name_full(name) {
            if !chr.is_empty() {
                genes
                    .entry(gene)
                    .or_insert((chr, FxHashSet::default()))
                    .1
                    .insert(pos);
            }
        }
    }
    let mut out: Vec<GeneEntry> = genes
        .into_iter()
        .map(|(gene, (chr, pos))| GeneEntry {
            gene: gene.into(),
            chr: chr.into(),
            lo: pos.iter().copied().min().unwrap_or(0),
            hi: pos.iter().copied().max().unwrap_or(0),
            sites: pos.len(),
        })
        .collect();
    out.sort_by(|a, b| b.sites.cmp(&a.sites).then_with(|| a.gene.cmp(&b.gene)));
    out
}

/// Every gene with site rows in `files`, from their row names.
pub fn gene_catalog(files: &[Box<str>]) -> anyhow::Result<Vec<GeneEntry>> {
    let mut names: Vec<Box<str>> = Vec::new();
    for file in files {
        let (backend, path) = resolve_backend_file(file, None)?;
        names.extend(open_sparse_matrix(&path, &backend)?.row_names()?);
    }
    Ok(catalog_from_rows(names.iter().map(|n| n.as_ref())))
}

/// What the user chose in the list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Choice {
    /// An entry, by index.
    Gene(usize),
    /// A typed `chr:start-end` or `chr:pos`.
    Locus(String),
    Quit,
}

/// State of the list, independent of the terminal so it can be tested.
pub struct GenePicker<'a> {
    entries: &'a [GeneEntry],
    filter: String,
    /// Indices into `entries` passing the filter.
    shown: Vec<usize>,
    /// Position in `shown`.
    selected: usize,
    offset: usize,
    /// A message for the footer, until the next key.
    status: Option<String>,
    decision: Option<Choice>,
}

impl<'a> GenePicker<'a> {
    pub fn new(entries: &'a [GeneEntry]) -> Self {
        Self {
            entries,
            filter: String::new(),
            shown: (0..entries.len()).collect(),
            selected: 0,
            offset: 0,
            status: None,
            decision: None,
        }
    }

    /// Case-insensitive substring match on the gene key or chromosome.
    fn refilter(&mut self) {
        let f = self.filter.to_lowercase();
        self.shown = (0..self.entries.len())
            .filter(|&i| {
                let e = &self.entries[i];
                f.is_empty() || e.gene.to_lowercase().contains(&f) || e.chr.to_lowercase() == f
            })
            .collect();
        self.selected = 0;
        self.offset = 0;
    }

    fn move_by(&mut self, delta: isize) {
        let last = self.shown.len().saturating_sub(1) as isize;
        self.selected = (self.selected as isize + delta).clamp(0, last) as usize;
    }

    pub fn filter(&self) -> &str {
        &self.filter
    }

    pub fn set_filter(&mut self, filter: &str) {
        self.filter = filter.to_string();
        self.refilter();
    }

    pub fn set_status(&mut self, msg: Option<String>) {
        self.status = msg;
    }

    /// Show the list until the user picks a gene, types a locus, or quits.
    pub fn pick(&mut self) -> anyhow::Result<Choice> {
        self.decision = None;
        run_screen(self)?;
        Ok(self.decision.take().unwrap_or(Choice::Quit))
    }
}

impl Screen for GenePicker<'_> {
    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        self.decision = Some(Choice::Quit);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.status = None;
        match key.code {
            KeyCode::Up => self.move_by(-1),
            KeyCode::Down => self.move_by(1),
            KeyCode::PageUp => self.move_by(-20),
            KeyCode::PageDown => self.move_by(20),
            KeyCode::Home => self.selected = 0,
            KeyCode::End => self.move_by(isize::MAX / 2),
            KeyCode::Enter => {
                if let Some(Query::Locus(..)) = parse_query(&self.filter) {
                    self.decision = Some(Choice::Locus(self.filter.clone()));
                } else if let Some(&i) = self.shown.get(self.selected) {
                    self.decision = Some(Choice::Gene(i));
                }
            }
            KeyCode::Esc if !self.filter.is_empty() => {
                self.filter.clear();
                self.refilter();
            }
            KeyCode::Esc => self.decision = Some(Choice::Quit),
            KeyCode::Backspace => {
                self.filter.pop();
                self.refilter();
            }
            KeyCode::Char(c) if !c.is_control() => {
                self.filter.push(c);
                self.refilter();
            }
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let [top, body, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
        ])
        .areas(frame.area());
        frame.render_widget(
            header(
                "pileup",
                "pick a gene",
                &format!("{} of {} genes", self.shown.len(), self.entries.len()),
            ),
            top,
        );
        let block = panel(" genes by sites ".into(), true);
        let inner = block.inner(body);
        frame.render_widget(block, body);

        let rows = inner.height.saturating_sub(1) as usize;
        if self.selected < self.offset {
            self.offset = self.selected;
        } else if rows > 0 && self.selected >= self.offset + rows {
            self.offset = self.selected + 1 - rows;
        }
        let dim = |t: String| Span::styled(t, DIM);
        let mut lines = vec![Line::from(dim(format!(
            "  {:<32} {:>7}  {}",
            "gene", "sites", "location"
        )))];
        for (k, &i) in self.shown.iter().enumerate().skip(self.offset).take(rows) {
            let e = &self.entries[i];
            let focused = k == self.selected;
            let style = if focused { HIGHLIGHT } else { PLAIN };
            lines.push(Line::from(vec![
                Span::styled(if focused { "▸ " } else { "  " }, HIGHLIGHT),
                Span::styled(format!("{:<32}", e.gene), style),
                Span::raw(format!(" {:>7}  ", e.sites)),
                dim(format!(
                    "{}:{}-{}",
                    e.chr,
                    fmt_thousands(e.lo),
                    fmt_thousands(e.hi)
                )),
            ]));
        }
        frame.render_widget(Paragraph::new(lines), inner);
        let line = match &self.status {
            Some(msg) => Line::from(Span::styled(format!(" {msg}"), HIGHLIGHT)),
            None => input_line(
                "gene or chr:start-end: ",
                &self.filter,
                &[("↑/↓", "move"), ("Enter", "open"), ("Esc", "clear/quit")],
            ),
        };
        frame.render_widget(line, footer);
    }
}

#[cfg(test)]
#[path = "tests/picker.rs"]
mod tests;
