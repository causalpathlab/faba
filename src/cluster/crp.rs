//! A Chinese-restaurant-process mixture over cells, every site's converted
//! fraction per cluster integrated out.
//!
//! A cell `c` in cluster `k` has, at each covered site `s`,
//!
//! ```text
//! on_cs | n_cs ~ Binomial(n_cs, p_sk),   p_sk ~ Beta(A p0_s, A (1 - p0_s)),
//! ```
//!
//! with `p0_s` the site's pooled fraction and `A` the prior strength. With
//! `p_sk` integrated out, a cell's move to a cluster is scored by the
//! Beta-binomial predictive of its reads given the cluster's others, so a
//! cluster is only its size and two counts per site. Each sweep:
//! single-cell Gibbs (Neal 2000, Algorithm 3), restricted-Gibbs split-merge
//! (Jain & Neal 2004), and the concentration `alpha` resampled (Escobar &
//! West 1995) under a Gamma(1, 1) prior.
//!
//! Split-merge anchors are two cells drawn at random, and the launch state a
//! random split of their clusters refined by restricted Gibbs scans (Jain &
//! Neal 2004). Chains start from k-means of a fixed picture of the cells
//! with a number of clusters drawn per chain, so they begin over- and
//! under-clustered alike.
//!
//! Reads can come weighted (importance weights that take batches out, see
//! `fit::batch_weights`): the likelihood is then a weighted one, each read's
//! term raised to its weight, which stays conjugate. The counts are no longer
//! whole and the Beta-binomial is taken through Beta functions. The sampler
//! is exact for that weighted target.

use super::read::Weighted;
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::RngExt;
use rand_distr::{Beta, Distribution, Gamma};
use statrs::function::gamma::ln_gamma;

pub(super) use statrs::function::beta::ln_beta;

/// A site's pooled converted fraction from `on` of `total` reads, with half
/// a read of pseudocount: the centre of every cluster's prior there.
pub(super) fn pooled_fraction(on: f64, total: f64) -> f64 {
    (on + 0.5) / (total + 1.0)
}

/// `ln B(x + m, y + u) - ln B(x, y)`: the Beta-binomial predictive of `m`
/// converted and `u` unconverted reads given `Beta(x, y)`, without the
/// binomial coefficient (the same for every cluster). Few reads, the usual
/// case, are summed as logs.
pub(super) fn ln_pred(x: f64, y: f64, m: f64, u: f64) -> f64 {
    if m + u > 16.0 || m.fract() != 0.0 || u.fract() != 0.0 {
        // ln Gamma(z + n) - ln Gamma(z), skipped when n is 0 (one side of
        // most entries).
        let rise = |z: f64, n: f64| {
            if n > 0.0 {
                ln_gamma(z + n) - ln_gamma(z)
            } else {
                0.0
            }
        };
        return rise(x, m) + rise(y, u) - rise(x + y, m + u);
    }
    // Rising factorials as products (at most 16 terms, far from overflow),
    // then two logs.
    let rising = |z: f64, n: f64| (0..n as u32).map(|j| z + j as f64).product::<f64>();
    (rising(x, m) * rising(y, u)).ln() - rising(x + y, m + u).ln()
}

/// The fit's settings.
#[derive(Debug, Clone)]
pub struct Settings {
    /// Prior strength `A`.
    pub strength: f64,
    /// Split-merge proposals per sweep.
    pub split_merge: usize,
    /// Restricted Gibbs scans before each split-merge proposal.
    pub restricted: usize,
    /// Most clusters a chain starts with (drawn from 1 up to this).
    pub start_clusters: usize,
    /// Single-cell Gibbs passes each sweep (off only to test split-merge alone).
    pub gibbs: bool,
    /// Hold `alpha` here instead of resampling it.
    pub fixed_alpha: Option<f64>,
}

/// A cluster: its cells, and its reads per site (dense over sites, one pair
/// of counts for each live cluster only).
struct Cluster {
    n: u32,
    on: Vec<f64>,
    off: Vec<f64>,
}

/// Two groups' reads during a split-merge proposal, reset after each one over
/// the sites it touched only.
struct Scratch {
    on: [Vec<f64>; 2],
    off: [Vec<f64>; 2],
    n: [u32; 2],
    seen: Vec<bool>,
    touched: Vec<u32>,
}

/// One chain's partition and what it carries.
pub struct State {
    z: Vec<u32>,
    /// Cluster slots; `None` is free.
    clusters: Vec<Option<Cluster>>,
    alpha: f64,
    scratch: Scratch,
    /// Emptied clusters kept for reuse, so a cell opening a new cluster does
    /// not allocate (their counts are zero up to rounding).
    spare: Vec<Cluster>,
    /// Split-merge moves accepted so far, as (splits, merges).
    accepted: (usize, usize),
}

/// One kept draw.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Each cell's cluster, numbered 0.. by first appearance.
    pub z: Vec<u16>,
    pub k: usize,
    pub alpha: f64,
    pub log_post: f64,
}

/// One chain's kept draws and its acceptance counts.
#[derive(Debug, Clone)]
pub struct Chain {
    pub samples: Vec<Sample>,
    /// Split-merge moves accepted.
    pub splits: usize,
    pub merges: usize,
}

pub struct Model<'a> {
    reads: &'a Weighted<'a>,
    /// `A p0_s` and `A (1 - p0_s)`, and `ln B` of the two.
    a: Vec<f64>,
    b: Vec<f64>,
    ln_b0: Vec<f64>,
    /// Each cell's predictive under the prior alone (a new cluster).
    lp_new: Vec<f64>,
    /// The picture of the cells chains start from.
    picture: legume_numeric::matrix::dmatrix_io::DMatrix<f32>,
    settings: Settings,
}

impl<'a> Model<'a> {
    pub fn new(reads: &'a Weighted<'a>, sites: usize, settings: Settings) -> Self {
        let (mut pon, mut tot) = (vec![0.0f64; sites], vec![0.0f64; sites]);
        for (i, &s) in reads.reads.site.iter().enumerate() {
            pon[s as usize] += reads.on[i];
            tot[s as usize] += reads.on[i] + reads.off[i];
        }
        let p0: Vec<f64> = (0..sites)
            .map(|s| pooled_fraction(pon[s], tot[s]))
            .collect();
        let a: Vec<f64> = p0.iter().map(|p| settings.strength * p).collect();
        let b: Vec<f64> = p0.iter().map(|p| settings.strength * (1.0 - p)).collect();
        let ln_b0: Vec<f64> = a.iter().zip(&b).map(|(&x, &y)| ln_beta(x, y)).collect();
        let lp_new = (0..reads.cells())
            .map(|c| {
                let (s, m, u) = reads.of(c);
                (0..s.len())
                    .map(|i| ln_pred(a[s[i] as usize], b[s[i] as usize], m[i], u[i]))
                    .sum()
            })
            .collect();
        let picture = super::embed::embed(reads, &p0);
        Self {
            reads,
            a,
            b,
            ln_b0,
            lp_new,
            picture,
            settings,
        }
    }

    fn sites(&self) -> usize {
        self.a.len()
    }

    fn cells(&self) -> usize {
        self.reads.cells()
    }

    /// Cell `c`'s log predictive given `on`/`off` counts per site.
    fn pred(&self, c: usize, on: &[f64], off: &[f64]) -> f64 {
        let (s, m, u) = self.reads.of(c);
        let mut v = 0.0;
        for i in 0..s.len() {
            let j = s[i] as usize;
            v += ln_pred(self.a[j] + on[j], self.b[j] + off[j], m[i], u[i]);
        }
        v
    }

    fn add(&self, cl: &mut Cluster, c: usize, sign: f64) {
        let (s, m, u) = self.reads.of(c);
        for i in 0..s.len() {
            let j = s[i] as usize;
            cl.on[j] = (cl.on[j] + sign * m[i]).max(0.0);
            cl.off[j] = (cl.off[j] + sign * u[i]).max(0.0);
        }
        cl.n = (cl.n as i64 + sign as i64) as u32;
    }

    /// A cluster with no cells, from the chain's spares when it has one.
    fn reuse(&self, st: &mut State) -> Cluster {
        st.spare.pop().unwrap_or_else(|| self.empty_cluster())
    }

    fn empty_cluster(&self) -> Cluster {
        Cluster {
            n: 0,
            on: vec![0.0; self.sites()],
            off: vec![0.0; self.sites()],
        }
    }

    fn free_slot(state: &mut State) -> usize {
        match state.clusters.iter().position(Option::is_none) {
            Some(k) => k,
            None => {
                state.clusters.push(None);
                state.clusters.len() - 1
            }
        }
    }

    /// One single-cell Gibbs pass over every cell, in random order.
    fn gibbs(&self, st: &mut State, rng: &mut SmallRng) {
        let n = self.cells();
        let mut order: Vec<usize> = (0..n).collect();
        order.shuffle(rng);
        let mut lp: Vec<f64> = Vec::new();
        let mut live: Vec<usize> = Vec::new();
        for c in order {
            let k0 = st.z[c] as usize;
            let cl = st.clusters[k0].as_mut().expect("a cell's cluster is live");
            self.add(cl, c, -1.0);
            if cl.n == 0 {
                let mut cl = st.clusters[k0].take().unwrap();
                for &j in self.reads.of(c).0 {
                    cl.on[j as usize] = 0.0;
                    cl.off[j as usize] = 0.0;
                }
                st.spare.push(cl);
            }
            live.clear();
            lp.clear();
            for (k, cl) in st.clusters.iter().enumerate() {
                if let Some(cl) = cl {
                    live.push(k);
                    lp.push((cl.n as f64).ln() + self.pred(c, &cl.on, &cl.off));
                }
            }
            lp.push(st.alpha.ln() + self.lp_new[c]);
            let pick = sample_log(&mut lp, rng);
            let k1 = if pick == live.len() {
                let k = Self::free_slot(st);
                st.clusters[k] = Some(self.reuse(st));
                k
            } else {
                live[pick]
            };
            self.add(st.clusters[k1].as_mut().unwrap(), c, 1.0);
            st.z[c] = k1 as u32;
        }
    }

    /// Put cell `c` on side `side` of the scratch (or take it off).
    fn put(&self, sc: &mut Scratch, c: usize, side: usize, sign: f64) {
        let (s, m, u) = self.reads.of(c);
        for i in 0..s.len() {
            let j = s[i] as usize;
            if !sc.seen[j] {
                sc.seen[j] = true;
                sc.touched.push(j as u32);
            }
            sc.on[side][j] = (sc.on[side][j] + sign * m[i]).max(0.0);
            sc.off[side][j] = (sc.off[side][j] + sign * u[i]).max(0.0);
        }
        sc.n[side] = (sc.n[side] as i64 + sign as i64) as u32;
    }

    /// Log marginal likelihood of a group's reads, over the touched sites
    /// (an untouched site contributes nothing). `which`: side 0, side 1, or
    /// both together (2).
    fn marg(&self, sc: &Scratch, which: usize) -> f64 {
        sc.touched
            .iter()
            .map(|&j| {
                let j = j as usize;
                let (on, off) = match which {
                    2 => (sc.on[0][j] + sc.on[1][j], sc.off[0][j] + sc.off[1][j]),
                    w => (sc.on[w][j], sc.off[w][j]),
                };
                ln_beta(self.a[j] + on, self.b[j] + off) - self.ln_b0[j]
            })
            .sum()
    }

    /// One restricted Gibbs scan of `rest` between the two sides; with
    /// `forced`, the sides are set to it rather than drawn. The log
    /// probability of the assignments made.
    fn scan(
        &self,
        sc: &mut Scratch,
        rest: &[usize],
        side: &mut [usize],
        forced: Option<&[usize]>,
        rng: &mut SmallRng,
    ) -> f64 {
        let mut logq = 0.0;
        for (t, &c) in rest.iter().enumerate() {
            self.put(sc, c, side[t], -1.0);
            let l0 = (sc.n[0] as f64).ln() + self.pred(c, &sc.on[0], &sc.off[0]);
            let l1 = (sc.n[1] as f64).ln() + self.pred(c, &sc.on[1], &sc.off[1]);
            let p0 = 1.0 / (1.0 + (l1 - l0).exp());
            let s1 = match forced {
                Some(f) => f[t],
                None => usize::from(rng.random::<f64>() >= p0),
            };
            logq += if s1 == 0 { p0 } else { 1.0 - p0 }.max(1e-300).ln();
            side[t] = s1;
            self.put(sc, c, s1, 1.0);
        }
        logq
    }

    /// One split-merge proposal; whether it was a split or a merge, if accepted.
    fn split_merge(&self, st: &mut State, rng: &mut SmallRng) -> Option<bool> {
        let n = self.cells();
        if n < 2 {
            return None;
        }
        let i = rng.random_range(0..n);
        let j = rng.random_range(0..n - 1);
        let j = j + usize::from(j >= i);
        let (ki, kj) = (st.z[i], st.z[j]);
        let mut rest: Vec<usize> = (0..n)
            .filter(|&c| c != i && c != j && (st.z[c] == ki || st.z[c] == kj))
            .collect();
        rest.shuffle(rng);
        let split = ki == kj;
        let forced: Option<Vec<usize>> =
            (!split).then(|| rest.iter().map(|&c| usize::from(st.z[c] != ki)).collect());
        let sc = &mut st.scratch;
        self.put(sc, i, 0, 1.0);
        self.put(sc, j, 1, 1.0);
        // Launch: a random split, then restricted scans.
        let mut side: Vec<usize> = rest.iter().map(|_| rng.random_range(0..2)).collect();
        for (t, &c) in rest.iter().enumerate() {
            self.put(sc, c, side[t], 1.0);
        }
        for _ in 0..self.settings.restricted {
            self.scan(sc, &rest, &mut side, None, rng);
        }
        let logq = self.scan(sc, &rest, &mut side, forced.as_deref(), rng);
        let (na, nb) = (sc.n[0] as f64, sc.n[1] as f64);
        let gain = st.alpha.ln() + ln_gamma(na) + ln_gamma(nb) - ln_gamma(na + nb)
            + self.marg(sc, 0)
            + self.marg(sc, 1)
            - self.marg(sc, 2);
        sc.reset();
        let ratio = if split { gain - logq } else { logq - gain };
        let accept = rng.random::<f64>().ln() < ratio;
        let mut outcome = None;
        if accept && split {
            let kb = Self::free_slot(st);
            let mut new = self.reuse(st);
            for (t, &c) in std::iter::once(&j).chain(&rest).enumerate() {
                let s_ = if t == 0 { 1 } else { side[t - 1] };
                if s_ == 1 {
                    self.add(st.clusters[ki as usize].as_mut().unwrap(), c, -1.0);
                    self.add(&mut new, c, 1.0);
                    st.z[c] = kb as u32;
                }
            }
            st.clusters[kb] = Some(new);
            outcome = Some(true);
        } else if accept {
            let mut from = st.clusters[kj as usize].take().unwrap();
            let into = st.clusters[ki as usize].as_mut().unwrap();
            for s in 0..self.sites() {
                into.on[s] += from.on[s];
                into.off[s] += from.off[s];
            }
            into.n += from.n;
            from.n = 0;
            from.on.fill(0.0);
            from.off.fill(0.0);
            st.spare.push(from);
            for z in st.z.iter_mut().filter(|z| **z == kj) {
                *z = ki;
            }
            outcome = Some(false);
        }
        outcome
    }

    /// `alpha` given the number of clusters (Escobar & West 1995), Gamma(1, 1) prior.
    fn resample_alpha(&self, st: &mut State, rng: &mut SmallRng) {
        let k = st.clusters.iter().flatten().count() as f64;
        let n = self.cells() as f64;
        let eta: f64 = Beta::new(st.alpha + 1.0, n).unwrap().sample(rng);
        let rate = 1.0 - eta.ln();
        let odds = k / (n * rate);
        let shape = if rng.random::<f64>() < odds / (1.0 + odds) {
            1.0 + k
        } else {
            k
        };
        st.alpha = Gamma::new(shape, 1.0 / rate).unwrap().sample(rng);
    }

    /// Log collapsed likelihood of a chain's partition.
    fn log_lik(&self, st: &State) -> f64 {
        let mut v = 0.0;
        for cl in st.clusters.iter().flatten() {
            for s in 0..self.sites() {
                if cl.on[s] + cl.off[s] > 1e-9 {
                    v += ln_beta(self.a[s] + cl.on[s], self.b[s] + cl.off[s]) - self.ln_b0[s];
                }
            }
        }
        v
    }

    /// Log CRP prior of a chain's partition given its `alpha`.
    fn log_prior(&self, st: &State) -> f64 {
        let n = self.cells() as f64;
        let mut v = ln_gamma(st.alpha) - ln_gamma(st.alpha + n);
        for cl in st.clusters.iter().flatten() {
            v += st.alpha.ln() + ln_gamma(cl.n as f64);
        }
        v
    }

    /// A chain started from k-means of the picture into a number of
    /// clusters drawn from 1 to `start_clusters`.
    fn start(&self, rng: &mut SmallRng) -> State {
        let k0 = rng
            .random_range(1..=self.settings.start_clusters.max(1))
            .min(self.cells().max(1));
        let z: Vec<u32> = if k0 > 1 {
            super::embed::kmeans(&self.picture, k0, rng.random())
        } else {
            vec![0; self.cells()]
        };
        let k0 = z.iter().max().map_or(1, |&k| k as usize + 1);
        let mut clusters: Vec<Option<Cluster>> =
            (0..k0).map(|_| Some(self.empty_cluster())).collect();
        for (c, &k) in z.iter().enumerate() {
            self.add(clusters[k as usize].as_mut().unwrap(), c, 1.0);
        }
        for cl in clusters.iter_mut() {
            if cl.as_ref().is_some_and(|c| c.n == 0) {
                *cl = None;
            }
        }
        State {
            z,
            clusters,
            alpha: self.settings.fixed_alpha.unwrap_or(1.0),
            scratch: Scratch::new(self.sites()),
            spare: Vec::new(),
            accepted: (0, 0),
        }
    }
}

impl Scratch {
    fn new(sites: usize) -> Self {
        Scratch {
            on: [vec![0.0; sites], vec![0.0; sites]],
            off: [vec![0.0; sites], vec![0.0; sites]],
            n: [0, 0],
            seen: vec![false; sites],
            touched: Vec::new(),
        }
    }

    /// Zero what the last proposal touched.
    fn reset(&mut self) {
        for &j in &self.touched {
            let j = j as usize;
            self.on[0][j] = 0.0;
            self.on[1][j] = 0.0;
            self.off[0][j] = 0.0;
            self.off[1][j] = 0.0;
            self.seen[j] = false;
        }
        self.touched.clear();
        self.n = [0, 0];
    }
}

/// Draw an index with probability proportional to `exp(lp)`; `lp` is
/// overwritten with the unnormalised weights.
fn sample_log(lp: &mut [f64], rng: &mut SmallRng) -> usize {
    let mx = lp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    lp.iter_mut().for_each(|l| *l = (*l - mx).exp());
    let mut r = rng.random::<f64>() * lp.iter().sum::<f64>();
    for (k, &x) in lp.iter().enumerate() {
        r -= x;
        if r <= 0.0 {
            return k;
        }
    }
    lp.len() - 1
}

/// Labels numbered 0.. by first appearance.
pub(super) fn compact(z: &[u32]) -> (Vec<u16>, usize) {
    let mut map: rustc_hash::FxHashMap<u32, u16> = Default::default();
    let out = z
        .iter()
        .map(|&k| {
            let next = map.len() as u16;
            *map.entry(k).or_insert(next)
        })
        .collect();
    (out, map.len())
}

/// One kept draw, with the chain's split-merge acceptance so far (the
/// engine hands a chain's summary only its draws, so the counts ride along
/// and the last draw's are kept).
pub struct Draw {
    sample: Sample,
    splits: usize,
    merges: usize,
}

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = Draw;
    type Result = Chain;

    fn init(&self, rng: &mut SmallRng) -> State {
        self.start(rng)
    }

    /// Gibbs, split-merge, alpha; each move keeps the posterior in detailed
    /// balance.
    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        if self.settings.gibbs {
            self.gibbs(st, rng);
        }
        for _ in 0..self.settings.split_merge {
            match self.split_merge(st, rng) {
                Some(true) => st.accepted.0 += 1,
                Some(false) => st.accepted.1 += 1,
                None => {}
            }
        }
        match self.settings.fixed_alpha {
            Some(a) => st.alpha = a,
            None => self.resample_alpha(st, rng),
        }
    }

    fn collect(&self, st: &State) -> Draw {
        let (z, k) = compact(&st.z);
        Draw {
            sample: Sample {
                z,
                k,
                alpha: st.alpha,
                log_post: self.log_prior(st) + self.log_lik(st),
            },
            splits: st.accepted.0,
            merges: st.accepted.1,
        }
    }

    fn summarize(&self, draws: Vec<Draw>) -> Chain {
        let (splits, merges) = draws.last().map_or((0, 0), |d| (d.splits, d.merges));
        Chain {
            samples: draws.into_iter().map(|d| d.sample).collect(),
            splits,
            merges,
        }
    }
}
