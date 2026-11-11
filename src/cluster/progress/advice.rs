//! What the chains' traces and split R-hat say of a fit, wherever it is
//! looked at: whether it has converged, and one line of advice naming the
//! flags to run it again with.

/// Split R-hat above this, 1.05, counted as not converged (the MCMC fits
/// of faba before 0.19, whose records the browser still reads).
const RHAT_OK: f64 = 1.05;

/// The advice on an MCMC fit (faba before 0.19) short of convergence: fit
/// again by variational inference, not run longer, as the sampler is gone.
pub(in crate::cluster) fn advice_again() -> String {
    "an MCMC fit of faba before 0.19: fit again with `faba cluster` (variational \
     inference) rather than read it"
        .to_string()
}

/// How far below the best chain a chain's mean log joint lags, in median
/// within-chain sds.
const FAR: f64 = 3.0;

/// The median of `v` (taken in any order).
fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    (v[(n - 1) / 2] + v[n / 2]) / 2.0
}

/// The chains whose mean log joint over their traced points after the
/// `warmup` burn-in sweeps is more than [`FAR`] median within-chain sds
/// below the best chain's (the lower of two, or most stuck below one);
/// none to say with fewer than two chains, or a chain with fewer than two
/// points past its burn-in. A chain's level, not its mixing: that is the
/// R-hat's.
pub fn lagging(traces: &[Vec<(usize, f64)>], warmup: usize) -> Option<Vec<usize>> {
    if traces.len() < 2 {
        return None;
    }
    let (means, sds): (Vec<f64>, Vec<f64>) = traces
        .iter()
        .map(|t| {
            let v: Vec<f64> = t.iter().filter(|p| p.0 > warmup).map(|p| p.1).collect();
            (v.len() >= 2).then(|| {
                let n = v.len() as f64;
                let mean = v.iter().sum::<f64>() / n;
                let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
                (mean, var.sqrt())
            })
        })
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .unzip();
    let best = means.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let bar = best - FAR * median(sds);
    Some((0..means.len()).filter(|&i| means[i] < bar).collect())
}

/// What the chains say together.
#[derive(Debug, Clone, PartialEq)]
pub enum Status {
    /// The R-hat is within [`RHAT_OK`] and no chain lags.
    Agree,
    /// One chain, its R-hat (its halves') within bounds: nothing to hold it
    /// against.
    OneChain,
    /// The R-hat is within bounds, but these chains' log joint lags the
    /// best's.
    Lags(Vec<usize>),
    /// The R-hat is past [`RHAT_OK`]; these chains lag too.
    Disagree(Vec<usize>),
    /// No R-hat to go by: too few draws.
    Unknown,
    /// A fit by variational inference: whether its best restart's ELBO
    /// settled within the tolerance, and how far below it the worst
    /// restart ended, per `unit` (count or read).
    Optimised {
        converged: bool,
        stopped: bool,
        asked: usize,
        spread: f64,
        unit: &'static str,
    },
}

/// A fit's convergence: its split R-hat and what the chains say, the
/// `--sweeps` it ran, which the advice doubles, its chains, and what else
/// is short of convergence (see [`Verdict::with_unsettled`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Verdict {
    pub rhat: f64,
    pub status: Status,
    sweeps: usize,
    chains: usize,
    unsettled: Vec<String>,
}

impl Verdict {
    /// The verdict on chains traced as `traces` (each empty for a fit that
    /// kept no trace: the R-hat alone), `warmup` sweeps of burn-in out of
    /// `sweeps`, split R-hat `rhat` (not finite when there were too few
    /// draws). Split R-hat above [`RHAT_OK`] counts as not
    /// converged, one chain's halves' too, as does a chain whose log joint
    /// lags (see [`lagging`]) with the R-hat within it.
    pub fn of(traces: &[Vec<(usize, f64)>], warmup: usize, sweeps: usize, rhat: f64) -> Self {
        let lag = lagging(traces, warmup).unwrap_or_default();
        let status = if !rhat.is_finite() {
            Status::Unknown
        } else if rhat > RHAT_OK {
            Status::Disagree(lag)
        } else if !lag.is_empty() {
            Status::Lags(lag)
        } else if traces.len() == 1 {
            Status::OneChain
        } else {
            Status::Agree
        };
        Self {
            rhat,
            status,
            sweeps,
            chains: traces.len(),
            unsettled: Vec::new(),
        }
    }

    /// The verdict on a fit by variational inference: its best restart
    /// `converged` (its ELBO settled) in `epochs` of the `asked` epochs at
    /// most, or was `stopped` early; `restarts` in all, the worst ending
    /// `spread` below the best per `unit` (count or read).
    pub fn vi(
        converged: bool,
        stopped: bool,
        (epochs, asked): (usize, usize),
        restarts: usize,
        (spread, unit): (f64, &'static str),
    ) -> Self {
        Self {
            rhat: f64::NAN,
            status: Status::Optimised {
                converged,
                stopped,
                asked,
                spread,
                unit,
            },
            sweeps: epochs,
            chains: restarts,
            unsettled: Vec::new(),
        }
    }

    /// The verdict, `unsettled` held against the fit too: what pooling
    /// the chains of an MCMC fit (faba before 0.19) found short of
    /// convergence, as its record names each in a few words.
    pub fn with_unsettled(self, unsettled: Vec<String>) -> Self {
        Self { unsettled, ..self }
    }

    /// Whether the fit can be read as it is.
    pub fn converged(&self) -> bool {
        self.converged_on_traces() && self.unsettled.is_empty()
    }

    /// The status in a few words, the R-hat in them; then what else is
    /// short of convergence.
    pub fn said(&self) -> String {
        let said = self.status_said();
        match self.unsettled.is_empty() {
            true => said,
            false => format!("{said}; not converged: {}", self.unsettled.join(", ")),
        }
    }

    /// The status in a few words, the R-hat in them.
    fn status_said(&self) -> String {
        let rhat = self.rhat;
        match &self.status {
            Status::Agree => format!("chains agree (R-hat {rhat:.2})"),
            Status::OneChain => format!("one chain (R-hat {rhat:.2})"),
            Status::Lags(lag) => format!("{} (R-hat {rhat:.2})", lags(lag)),
            Status::Disagree(_) if self.chains == 1 => {
                format!("R-hat {rhat:.2}: one chain, its halves disagree")
            }
            Status::Disagree(lag) if lag.is_empty() => {
                format!("R-hat {rhat:.2}: chains disagree")
            }
            Status::Disagree(lag) => format!("R-hat {rhat:.2}: chains disagree, {}", lags(lag)),
            Status::Unknown => "R-hat not known".to_string(),
            Status::Optimised {
                converged,
                stopped,
                spread,
                unit,
                ..
            } => {
                let settled = match (converged, stopped) {
                    (true, _) => format!("ELBO settled in {} epochs", self.sweeps),
                    (false, true) => format!("stopped early after {} epochs", self.sweeps),
                    (false, false) => format!("ELBO still rising after {} epochs", self.sweeps),
                };
                match self.chains {
                    0 | 1 => settled,
                    n => format!(
                        "{settled}; {n} restarts, the worst {spread:.4} per {unit} below the best"
                    ),
                }
            }
        }
    }

    /// One line on what to do, naming the flags.
    pub fn advice(&self) -> String {
        if !self.unsettled.is_empty() && self.converged_on_traces() {
            return advice_again();
        }
        match self.status {
            Status::Agree => "nothing to rerun: the fit can be read as it is".to_string(),
            Status::OneChain => {
                "one chain, nothing to hold it against: fit again with `faba cluster` to compare restarts"
                    .to_string()
            }
            Status::Lags(_) | Status::Disagree(_) | Status::Unknown => advice_again(),
            Status::Optimised {
                converged: true, ..
            } => "nothing to rerun: the fit can be read as it is (the best restart's tables)"
                .to_string(),
            Status::Optimised {
                converged: false,
                stopped: true,
                asked,
                ..
            } => {
                format!("stopped before it settled: run again to let it finish (--epochs {asked})")
            }
            Status::Optimised {
                converged: false,
                asked,
                ..
            } => format!("run again with --epochs {}", 2 * asked.max(1)),
        }
    }

    /// Whether the traces and the R-hat alone say the fit can be read.
    fn converged_on_traces(&self) -> bool {
        matches!(
            self.status,
            Status::Agree
                | Status::OneChain
                | Status::Optimised {
                    converged: true,
                    ..
                }
        )
    }

    /// The R-hat and `ok` or `check`, for a header.
    pub fn badge(&self) -> String {
        let word = if self.converged() { "ok" } else { "check" };
        match self.status {
            Status::Optimised { .. } => format!("ELBO {word}"),
            Status::Unknown => format!("R-hat - {word}"),
            _ => format!("R-hat {:.2} {word}", self.rhat),
        }
    }
}

/// `chain 2 lags`, or `chains 1, 2 lag`.
pub(in crate::cluster) fn lags(lag: &[usize]) -> String {
    let names = lag.iter().map(usize::to_string).collect::<Vec<_>>();
    match names.as_slice() {
        [k] => format!("chain {k} lags"),
        _ => format!("chains {} lag", names.join(", ")),
    }
}

#[cfg(test)]
mod tests;
