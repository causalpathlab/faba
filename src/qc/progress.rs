//! How far a `faba qc` write has got, shared between the thread writing the
//! fileset and the view showing it.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;

#[derive(Default)]
pub struct Progress {
    total: AtomicUsize,
    started: AtomicUsize,
    step: Mutex<String>,
    finished: AtomicBool,
}

impl Progress {
    /// Progress with nothing left to do.
    pub const fn finished() -> Self {
        Progress {
            total: AtomicUsize::new(0),
            started: AtomicUsize::new(0),
            step: Mutex::new(String::new()),
            finished: AtomicBool::new(true),
        }
    }

    /// Expect `total` steps.
    pub fn plan(&self, total: usize) {
        self.total.store(total, Ordering::Relaxed);
    }

    /// Start the next step, which finishes the one before.
    pub fn next(&self, what: impl Into<String>) {
        let what = what.into();
        log::debug!("qc: {what}");
        *self.step.lock().unwrap_or_else(|e| e.into_inner()) = what;
        self.started.fetch_add(1, Ordering::Relaxed);
    }

    /// The writer stopped, done or not.
    pub fn finish(&self) {
        self.finished.store(true, Ordering::Release);
    }

    pub fn is_finished(&self) -> bool {
        self.finished.load(Ordering::Acquire)
    }

    /// Steps done.
    pub fn done(&self) -> usize {
        let started = self.started.load(Ordering::Relaxed);
        if self.is_finished() {
            started
        } else {
            started.saturating_sub(1)
        }
    }

    /// Steps started, to check a plan against.
    pub fn started(&self) -> usize {
        self.started.load(Ordering::Relaxed)
    }

    /// Steps done, steps planned, and the step running.
    pub fn snapshot(&self) -> (usize, usize, String) {
        let (total, done) = (self.total.load(Ordering::Relaxed), self.done());
        let step = self.step.lock().unwrap_or_else(|e| e.into_inner()).clone();
        (done.min(total), total, step)
    }
}

/// Marks `progress` finished when dropped, so a writer that panics or
/// returns early still lets the view close.
pub struct FinishOnDrop<'a>(pub &'a Progress);

impl Drop for FinishOnDrop<'_> {
    fn drop(&mut self) {
        self.0.finish();
    }
}
