//! The Steps screen: which pipeline steps run.

use crate::pipeline::args::PipelineArgs;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Step {
    Snp,
    Count,
    Depth,
    Atoi,
    M6a,
    Apa,
}

impl Step {
    pub const ALL: [Step; 6] = [
        Step::Snp,
        Step::Count,
        Step::Depth,
        Step::Atoi,
        Step::M6a,
        Step::Apa,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Step::Snp => "SNP",
            Step::Count => "count",
            Step::Depth => "depth",
            Step::Atoi => "ATOI",
            Step::M6a => "m6A",
            Step::Apa => "APA",
        }
    }

    pub fn about(self) -> &'static str {
        match self {
            Step::Snp => "genotype de novo over fg + bg, with known SNPs if given",
            Step::Count => "genes and called cells; these restrict every step below",
            Step::Depth => "per-cell read depth in bins",
            Step::Atoi => "every putative A-to-I site, tested against the error null",
            Step::M6a => "fg vs bg contrast",
            Step::Apa => "poly(A) site usage, runs last",
        }
    }

    /// The flag that turns the step off; depth is off by leaving out its
    /// resolution instead.
    fn skip_flag(self) -> Option<&'static str> {
        match self {
            Step::Snp => Some("--skip-snp"),
            Step::Count => Some("--skip-count"),
            Step::Atoi => Some("--skip-atoi"),
            Step::M6a => Some("--skip-m6a"),
            Step::Apa => Some("--skip-apa"),
            Step::Depth => None,
        }
    }
}

pub struct Steps {
    pub on: [bool; 6],
    pub depth_kb: String,
    pub at: usize,
}

impl Default for Steps {
    fn default() -> Self {
        Steps {
            on: [true, true, false, true, true, true],
            depth_kb: String::new(),
            at: 0,
        }
    }
}

impl Steps {
    pub fn is_on(&self, s: Step) -> bool {
        self.on[s as usize]
    }

    /// Space on the highlighted step; a note when it cannot be turned on.
    pub fn toggle(&mut self, has_bg: bool) -> Option<&'static str> {
        let s = Step::ALL[self.at];
        if s == Step::M6a && !self.is_on(s) && !has_bg {
            return Some("m6A needs a bg BAM");
        }
        self.on[s as usize] ^= true;
        None
    }

    /// The `--skip-*` flags of the steps that are off, and depth's
    /// resolution when it is on. m6A with no bg BAM is skipped by the
    /// pipeline itself, so it needs no flag.
    pub fn argv(&self, has_bg: bool) -> Vec<String> {
        let mut v: Vec<String> = Step::ALL
            .iter()
            .filter(|s| !self.is_on(**s) && !(**s == Step::M6a && !has_bg))
            .filter_map(|s| s.skip_flag())
            .map(String::from)
            .collect();
        if self.is_on(Step::Depth) && !self.depth_kb.trim().is_empty() {
            v.extend(["--depth-resolution-kb".into(), self.depth_kb.trim().into()]);
        }
        v
    }

    pub fn problems(&self) -> Vec<String> {
        let kb_ok = self
            .depth_kb
            .trim()
            .parse::<f32>()
            .is_ok_and(|x| x.is_finite() && x > 0.0);
        if self.is_on(Step::Depth) && !kb_ok {
            vec!["depth is on: give its resolution in kb".into()]
        } else {
            Vec::new()
        }
    }

    /// Whether a flag group's step runs: off when its step is off, and m6A
    /// also without a bg BAM, since the pipeline then skips it. `Common` and
    /// headings that name no step always run.
    pub fn heading_on(&self, heading: &str, has_bg: bool) -> bool {
        Step::ALL
            .iter()
            .find(|s| s.label() == heading)
            .is_none_or(|s| self.is_on(*s) && (*s != Step::M6a || has_bg))
    }

    pub fn prefill(&mut self, a: &PipelineArgs) {
        self.on[Step::Snp as usize] = !a.skip_snp;
        self.on[Step::Count as usize] = !a.skip_count;
        self.on[Step::Atoi as usize] = !a.skip_atoi;
        self.on[Step::M6a as usize] = !a.skip_m6a;
        self.on[Step::Apa as usize] = !a.skip_apa;
        if let Some(kb) = a.depth_resolution_kb {
            self.on[Step::Depth as usize] = true;
            self.depth_kb = kb.to_string();
        }
    }
}

#[cfg(test)]
#[path = "tests/steps.rs"]
mod tests;
