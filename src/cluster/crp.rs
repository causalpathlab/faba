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
//! Split-merge anchors: the first cell at random; the second, a third of the
//! time each, at random, among the first's nearest cells in a fixed picture
//! of the data (likely merges), or from another region of that picture
//! (likely splits of a mixed cluster). The choice depends on the data alone,
//! never on the partition, so its probability is the same for a move and its
//! reverse and cancels: every move keeps detailed balance. The launch state is allocated one cell at a time by predictive
//! (Dahl's sequential allocation), then refined by restricted Gibbs scans;
//! it too depends on nothing in the current split. Chains start from k-means
//! of that picture with a number of clusters drawn per chain, so chains
//! begin over- and under-clustered alike.
//!
//! Reads can come weighted (importance weights that take batches out, see
//! `fit::batch_weights`): the likelihood is then a weighted one, each read's
//! term raised to its weight, which stays conjugate. The counts are no longer
//! whole and the Beta-binomial is taken through Beta functions. The sampler
//! is exact for that weighted target.

use super::read::Reads;
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::{Beta, Distribution, Gamma};
use statrs::function::gamma::ln_gamma;

/// `ln B(a, b)`.
pub(super) fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// `ln B(x + m, y + u) - ln B(x, y)`: the Beta-binomial predictive of `m`
/// converted and `u` unconverted reads given `Beta(x, y)`, without the
/// binomial coefficient (the same for every cluster). Few reads, the usual
/// case, are summed as logs.
pub(super) fn ln_pred(x: f64, y: f64, m: f64, u: f64) -> f64 {
    if m + u > 16.0 || m.fract() != 0.0 || u.fract() != 0.0 {
        return ln_beta(x + m, y + u) - ln_beta(x, y);
    }
    let (m, u) = (m as u32, u as u32);
    let mut v = 0.0;
    for j in 0..m {
        v += (x + j as f64).ln();
    }
    for j in 0..u {
        v += (y + j as f64).ln();
    }
    for j in 0..(m + u) {
        v -= (x + y + j as f64).ln();
    }
    v
}

/// The fit's settings.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Settings {
    /// Prior strength `A`.
    pub strength: f64,
    /// Split-merge proposals per sweep.
    pub split_merge: usize,
    /// Restricted Gibbs scans before each split-merge proposal.
    pub restricted: usize,
    /// Most clusters a chain starts with (drawn from 1 up to this).
    pub start_clusters: usize,
    /// Nearest cells offered as a split-merge anchor's partner (0: anchors
    /// at random only).
    pub neighbours: usize,
    /// Single-cell Gibbs passes each sweep (off only to test split-merge alone).
    pub gibbs: bool,
    /// Hold `alpha` here instead of resampling it.
    pub fixed_alpha: Option<f64>,
}

/// A cluster: its cells, and its reads per site (dense over sites, one pair
/// of counts for each live cluster only).
#[derive(Clone)]
struct Cluster {
    n: u32,
    on: Vec<f64>,
    off: Vec<f64>,
}

/// Two groups' reads during a split-merge proposal, reset after each one over
/// the sites it touched only.
#[derive(Clone)]
struct Scratch {
    on: [Vec<f64>; 2],
    off: [Vec<f64>; 2],
    n: [u32; 2],
    seen: Vec<bool>,
    touched: Vec<u32>,
}

#[derive(Clone)]
pub struct State {
    z: Vec<u32>,
    /// Cluster slots; `None` is free.
    clusters: Vec<Option<Cluster>>,
    alpha: f64,
    scratch: Scratch,
    pub accepted_splits: usize,
    pub accepted_merges: usize,
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
    pub splits: usize,
    pub merges: usize,
}

pub struct Model<'a> {
    reads: &'a Reads,
    /// Each entry's converted and unconverted reads, reweighted.
    on: Vec<f64>,
    off: Vec<f64>,
    /// `A p0_s` and `A (1 - p0_s)`.
    a: Vec<f64>,
    b: Vec<f64>,
    /// Each cell's predictive under the prior alone (a new cluster).
    lp_new: Vec<f64>,
    /// Each cell's nearest cells, for split-merge anchors.
    near: Vec<Vec<u32>>,
    /// Each cell's region of the picture, and the cells of each region.
    region: Vec<u32>,
    regions: Vec<Vec<u32>>,
    /// The picture of the cells chains start from.
    picture: legume_numeric::matrix::dmatrix_io::DMatrix<f32>,
    settings: Settings,
}

impl<'a> Model<'a> {
    /// `weights`: each entry's converted and unconverted reads' weights.
    pub fn new(
        reads: &'a Reads,
        sites: usize,
        settings: Settings,
        weights: Option<(&[f64], &[f64])>,
    ) -> Self {
        let on: Vec<f64> = (0..reads.site.len())
            .map(|i| reads.on[i] as f64 * weights.map_or(1.0, |w| w.0[i]))
            .collect();
        let off: Vec<f64> = (0..reads.site.len())
            .map(|i| reads.off[i] as f64 * weights.map_or(1.0, |w| w.1[i]))
            .collect();
        let (mut pon, mut tot) = (vec![0.0f64; sites], vec![0.0f64; sites]);
        for i in 0..reads.site.len() {
            let s = reads.site[i] as usize;
            pon[s] += on[i];
            tot[s] += on[i] + off[i];
        }
        let p0: Vec<f64> = (0..sites)
            .map(|s| (pon[s] + 0.5) / (tot[s] + 1.0))
            .collect();
        let a: Vec<f64> = p0.iter().map(|p| settings.strength * p).collect();
        let b: Vec<f64> = p0.iter().map(|p| settings.strength * (1.0 - p)).collect();
        let lp_new = (0..reads.cells())
            .map(|c| {
                let r = reads.ptr[c]..reads.ptr[c + 1];
                r.map(|i| {
                    ln_pred(
                        a[reads.site[i] as usize],
                        b[reads.site[i] as usize],
                        on[i],
                        off[i],
                    )
                })
                .sum()
            })
            .collect();
        let picture = super::embed::embed(&reads.ptr, &reads.site, &on, &off, &p0);
        let near = if settings.neighbours > 0 && reads.cells() > 1 {
            super::embed::neighbours(&picture, settings.neighbours.min(reads.cells() - 1))
        } else {
            vec![Vec::new(); reads.cells()]
        };
        let region = if reads.cells() > REGIONS {
            super::embed::kmeans(&picture, REGIONS, REGION_SEED)
        } else {
            (0..reads.cells() as u32).collect()
        };
        let mut regions = vec![Vec::new(); region.iter().max().map_or(0, |&r| r as usize + 1)];
        for (c, &r) in region.iter().enumerate() {
            regions[r as usize].push(c as u32);
        }
        regions.retain(|r| !r.is_empty());
        let region = {
            let mut out = vec![0u32; reads.cells()];
            for (r, cells) in regions.iter().enumerate() {
                for &c in cells {
                    out[c as usize] = r as u32;
                }
            }
            out
        };
        Self {
            reads,
            on,
            off,
            near,
            region,
            regions,
            picture,
            a,
            b,
            lp_new,
            settings,
        }
    }

    fn sites(&self) -> usize {
        self.a.len()
    }

    fn cells(&self) -> usize {
        self.reads.cells()
    }

    /// Cell `c`'s sites and (reweighted) reads.
    fn entries(&self, c: usize) -> (&[u32], &[f64], &[f64]) {
        let r = self.reads.ptr[c]..self.reads.ptr[c + 1];
        (
            &self.reads.site[r.clone()],
            &self.on[r.clone()],
            &self.off[r],
        )
    }

    /// Cell `c`'s log predictive given `on`/`off` counts per site.
    fn pred(&self, c: usize, on: &[f64], off: &[f64]) -> f64 {
        let (s, m, u) = self.entries(c);
        let mut v = 0.0;
        for i in 0..s.len() {
            let j = s[i] as usize;
            v += ln_pred(self.a[j] + on[j], self.b[j] + off[j], m[i], u[i]);
        }
        v
    }

    fn add(&self, cl: &mut Cluster, c: usize, sign: f64) {
        let (s, m, u) = self.entries(c);
        for i in 0..s.len() {
            let j = s[i] as usize;
            cl.on[j] = (cl.on[j] + sign * m[i]).max(0.0);
            cl.off[j] = (cl.off[j] + sign * u[i]).max(0.0);
        }
        cl.n = (cl.n as i64 + sign as i64) as u32;
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
        shuffle(&mut order, rng);
        let mut lp: Vec<f64> = Vec::new();
        let mut live: Vec<usize> = Vec::new();
        for c in order {
            let k0 = st.z[c] as usize;
            let cl = st.clusters[k0].as_mut().expect("a cell's cluster is live");
            self.add(cl, c, -1.0);
            if cl.n == 0 {
                st.clusters[k0] = None;
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
            let pick = sample_log(&lp, rng);
            let k1 = if pick == live.len() {
                let k = Self::free_slot(st);
                st.clusters[k] = Some(self.empty_cluster());
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
        let (s, m, u) = self.entries(c);
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
                ln_beta(self.a[j] + on, self.b[j] + off) - ln_beta(self.a[j], self.b[j])
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

    /// A split-merge anchor's partner, by a rule that looks at the data only:
    /// at random, a near cell, or a cell of another region, a third each.
    fn partner(&self, i: usize, rng: &mut SmallRng) -> usize {
        let n = self.cells();
        let near = &self.near[i];
        let u = rng.random::<f64>();
        if u < 1.0 / 3.0 && !near.is_empty() {
            return near[rng.random_range(0..near.len())] as usize;
        }
        if u < 2.0 / 3.0 && self.regions.len() > 1 {
            let mut r = rng.random_range(0..self.regions.len() - 1);
            r += usize::from(r >= self.region[i] as usize);
            let cells = &self.regions[r];
            return cells[rng.random_range(0..cells.len())] as usize;
        }
        let j = rng.random_range(0..n - 1);
        j + usize::from(j >= i)
    }

    /// One split-merge proposal; whether it was a split or a merge, if accepted.
    fn split_merge(&self, st: &mut State, rng: &mut SmallRng) -> Option<bool> {
        let n = self.cells();
        if n < 2 {
            return None;
        }
        let i = rng.random_range(0..n);
        let j = self.partner(i, rng);
        let (ki, kj) = (st.z[i], st.z[j]);
        let mut rest: Vec<usize> = (0..n)
            .filter(|&c| c != i && c != j && (st.z[c] == ki || st.z[c] == kj))
            .collect();
        shuffle(&mut rest, rng);
        let mut sc = std::mem::replace(&mut st.scratch, Scratch::empty());
        self.put(&mut sc, i, 0, 1.0);
        self.put(&mut sc, j, 1, 1.0);
        // Launch: each cell in turn to the side that predicts it, as drawn.
        let mut side: Vec<usize> = Vec::with_capacity(rest.len());
        for &c in &rest {
            let l0 = (sc.n[0] as f64).ln() + self.pred(c, &sc.on[0], &sc.off[0]);
            let l1 = (sc.n[1] as f64).ln() + self.pred(c, &sc.on[1], &sc.off[1]);
            let s1 = usize::from(rng.random::<f64>() * (1.0 + (l1 - l0).exp()) >= 1.0);
            self.put(&mut sc, c, s1, 1.0);
            side.push(s1);
        }
        for _ in 0..self.settings.restricted {
            self.scan(&mut sc, &rest, &mut side, None, rng);
        }
        let split = ki == kj;
        let forced: Option<Vec<usize>> =
            (!split).then(|| rest.iter().map(|&c| usize::from(st.z[c] != ki)).collect());
        let logq = self.scan(&mut sc, &rest, &mut side, forced.as_deref(), rng);
        let (na, nb) = (sc.n[0] as f64, sc.n[1] as f64);
        let gain = st.alpha.ln() + ln_gamma(na) + ln_gamma(nb) - ln_gamma(na + nb)
            + self.marg(&sc, 0)
            + self.marg(&sc, 1)
            - self.marg(&sc, 2);
        let ratio = if split { gain - logq } else { logq - gain };
        let accept = rng.random::<f64>().ln() < ratio;
        let mut outcome = None;
        if accept && split {
            let kb = Self::free_slot(st);
            let mut new = self.empty_cluster();
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
            let from = st.clusters[kj as usize].take().unwrap();
            let into = st.clusters[ki as usize].as_mut().unwrap();
            for s in 0..self.sites() {
                into.on[s] += from.on[s];
                into.off[s] += from.off[s];
            }
            into.n += from.n;
            for z in st.z.iter_mut().filter(|z| **z == kj) {
                *z = ki;
            }
            outcome = Some(false);
        }
        sc.reset();
        st.scratch = sc;
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

    /// Log joint of the partition and reads (CRP prior times the collapsed likelihood).
    fn log_post(&self, st: &State) -> f64 {
        let n = self.cells() as f64;
        let mut v = ln_gamma(st.alpha) - ln_gamma(st.alpha + n);
        for cl in st.clusters.iter().flatten() {
            v += st.alpha.ln() + ln_gamma(cl.n as f64);
            for s in 0..self.sites() {
                if cl.on[s] + cl.off[s] > 1e-9 {
                    v += ln_beta(self.a[s] + cl.on[s], self.b[s] + cl.off[s])
                        - ln_beta(self.a[s], self.b[s]);
                }
            }
        }
        v
    }
}

impl Scratch {
    fn empty() -> Self {
        Scratch {
            on: [Vec::new(), Vec::new()],
            off: [Vec::new(), Vec::new()],
            n: [0, 0],
            seen: Vec::new(),
            touched: Vec::new(),
        }
    }

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

/// Regions of the picture for far split-merge partners.
const REGIONS: usize = 16;

/// Fixed, so regions are the same for every chain and seed.
const REGION_SEED: u64 = 0x7e61_0105;

/// Fisher-Yates.
fn shuffle<T>(v: &mut [T], rng: &mut SmallRng) {
    for i in (1..v.len()).rev() {
        v.swap(i, rng.random_range(0..=i));
    }
}

/// Draw an index with probability proportional to `exp(lp)`.
fn sample_log(lp: &[f64], rng: &mut SmallRng) -> usize {
    let mx = lp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = lp.iter().map(|&l| (l - mx).exp()).collect();
    let mut r = rng.random::<f64>() * w.iter().sum::<f64>();
    for (k, &x) in w.iter().enumerate() {
        r -= x;
        if r <= 0.0 {
            return k;
        }
    }
    w.len() - 1
}

/// Labels numbered 0.. by first appearance.
fn compact(z: &[u32]) -> (Vec<u16>, usize) {
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

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = (Sample, usize, usize);
    type Result = Chain;

    fn init(&self, rng: &mut SmallRng) -> State {
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
            accepted_splits: 0,
            accepted_merges: 0,
        }
    }

    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        if self.settings.gibbs {
            self.gibbs(st, rng);
        }
        for _ in 0..self.settings.split_merge {
            match self.split_merge(st, rng) {
                Some(true) => st.accepted_splits += 1,
                Some(false) => st.accepted_merges += 1,
                None => {}
            }
        }
        match self.settings.fixed_alpha {
            Some(a) => st.alpha = a,
            None => self.resample_alpha(st, rng),
        }
    }

    fn collect(&self, st: &State) -> (Sample, usize, usize) {
        let (z, k) = compact(&st.z);
        let sample = Sample {
            z,
            k,
            alpha: st.alpha,
            log_post: self.log_post(st),
        };
        (sample, st.accepted_splits, st.accepted_merges)
    }

    fn summarize(&self, samples: Vec<(Sample, usize, usize)>) -> Chain {
        let (splits, merges) = samples.last().map_or((0, 0), |s| (s.1, s.2));
        Chain {
            samples: samples.into_iter().map(|s| s.0).collect(),
            splits,
            merges,
        }
    }
}
