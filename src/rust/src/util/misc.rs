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

/// paste words in a vector of `Box<str>` into `Box<str>`
///
/// * `words`
/// * `indices`
/// * `sep`
#[allow(dead_code)]
pub fn paste(words: &Vec<Box<str>>, indices: &Vec<usize>, sep: &str) -> Box<str> {
    let mut ret = String::new();
    let n = indices.len();
    for (i, j) in indices.iter().enumerate() {
        if let Some(w) = words.get(*j) {
            ret.push_str(w);
        }
        if n > 1 && i < (n - 1) {
            ret.push_str(sep);
        }
    }
    ret.into_boxed_str()
}
