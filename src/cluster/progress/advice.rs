//! What the chains' traces and split R-hat say of a fit, wherever it is
//! looked at: whether it has converged, and one line of advice naming the
//! flags to run it again with.

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
    /// The R-hat is within [`crate::cluster::RHAT_OK`] and no chain lags.
    Agree,
    /// One chain, its R-hat (its halves') within bounds: nothing to hold it
    /// against.
    OneChain,
    /// The R-hat is within bounds, but these chains' log joint lags the
    /// best's.
    Lags(Vec<usize>),
    /// The R-hat is past [`crate::cluster::RHAT_OK`]; these chains lag too.
    Disagree(Vec<usize>),
    /// No R-hat to go by: too few draws.
    Unknown,
}

/// A fit's convergence: its split R-hat and what the chains say, the
/// `--sweeps` it ran, which the advice doubles, and its chains.
#[derive(Debug, Clone, PartialEq)]
pub struct Verdict {
    pub rhat: f64,
    pub status: Status,
    sweeps: usize,
    chains: usize,
}

impl Verdict {
    /// The verdict on chains traced as `traces` (each empty for a fit that
    /// kept no trace: the R-hat alone), `warmup` sweeps of burn-in out of
    /// `sweeps`, split R-hat `rhat` (not finite when there were too few
    /// draws). Split R-hat above [`crate::cluster::RHAT_OK`] counts as not
    /// converged, one chain's halves' too, as does a chain whose log joint
    /// lags (see [`lagging`]) with the R-hat within it.
    pub fn of(traces: &[Vec<(usize, f64)>], warmup: usize, sweeps: usize, rhat: f64) -> Self {
        let lag = lagging(traces, warmup).unwrap_or_default();
        let status = if !rhat.is_finite() {
            Status::Unknown
        } else if rhat > crate::cluster::RHAT_OK {
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
        }
    }

    /// Whether the fit can be read as it is.
    pub fn converged(&self) -> bool {
        matches!(self.status, Status::Agree | Status::OneChain)
    }

    /// The status in a few words, the R-hat in them.
    pub fn said(&self) -> String {
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
        }
    }

    /// One line on what to do, naming the flags.
    pub fn advice(&self) -> String {
        let more = 2 * self.sweeps.max(1);
        match self.status {
            Status::Agree => "nothing to rerun: the fit can be read as it is".to_string(),
            Status::OneChain => {
                "nothing to hold it against: run again with --chains 4 to see whether others agree"
                    .to_string()
            }
            Status::Lags(_) => format!("run again with --sweeps {more} or more --chains"),
            Status::Disagree(_) if self.chains == 1 => {
                format!("run longer (--sweeps {more}) with --chains 4 and check the trace")
            }
            Status::Disagree(_) => format!("run longer (--sweeps {more}) and check the traces"),
            Status::Unknown => {
                format!("too few draws to compare the chains: run longer (--sweeps {more})")
            }
        }
    }

    /// The R-hat and `ok` or `check`, for a header.
    pub fn badge(&self) -> String {
        let rhat = match self.status {
            Status::Unknown => "-".to_string(),
            _ => format!("{:.2}", self.rhat),
        };
        let word = if self.converged() { "ok" } else { "check" };
        format!("R-hat {rhat} {word}")
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
