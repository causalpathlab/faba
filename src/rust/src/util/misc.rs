/// Other utilities
/// make a vector of intervals
pub fn make_intervals(max_size: i64, block_size: i64) -> Vec<(i64, i64)> {
    let mut jobs = vec![];
    for lb in (0..max_size).step_by(block_size as usize) {
        let ub = std::cmp::min(max_size, lb + block_size);
        jobs.push((lb, ub));
    }
    return jobs;
}
