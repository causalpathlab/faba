//! Helpers the tests of several modules share.

/// Wait, yielding, until `done` says so; a test fails after 30 s of it.
pub(crate) fn wait_for(mut done: impl FnMut() -> bool) {
    let start = std::time::Instant::now();
    while !done() {
        assert!(start.elapsed().as_secs() < 30, "waited 30 s in vain");
        std::thread::yield_now();
    }
}
