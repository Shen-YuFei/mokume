use mokume_core::{MokumeError, Result};

pub(crate) fn install<T, F>(threads: Option<usize>, operation: F) -> Result<T>
where
    T: Send,
    F: FnOnce() -> Result<T> + Send,
{
    if threads == Some(0) {
        return Err(MokumeError::InvalidInput {
            message: "thread count must be greater than zero".to_owned(),
        });
    }
    let Some(threads) = threads else {
        return operation();
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|error| MokumeError::InvalidInput {
            message: format!("failed to create {threads}-thread worker pool: {error}"),
        })?;
    pool.install(operation)
}

#[cfg(test)]
mod tests {
    use super::install;

    #[test]
    fn consecutive_calls_can_select_different_thread_counts() {
        let first = install(Some(1), || Ok(rayon::current_num_threads()));
        let second = install(Some(4), || Ok(rayon::current_num_threads()));

        assert_eq!(first.ok(), Some(1));
        assert_eq!(second.ok(), Some(4));
    }
}
