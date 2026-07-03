//! Thin `mokume` CLI binary entry point.
//!
//! All command parsing and dispatch lives in the `mokume_cli` library crate so
//! the same logic is reachable both from this standalone binary and from the
//! `mokume-py` PyO3 bindings (`mokume_cli::run_from_args`). This binary keeps the
//! "no Python runtime needed" path for pipeline use (quantms / nextflow).

use std::process::ExitCode;

fn main() -> ExitCode {
    mokume_cli::run()
}
