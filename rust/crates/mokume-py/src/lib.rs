//! PyO3 bindings exposing the mokume Rust compute kernel to Python.
//!
//! maturin builds this crate into the `mokume._mokume` extension module that the
//! `rust/python/mokume/` package imports. It is the FFI boundary between the Rust
//! compute crates and the Python periphery (plotting / tissue maps / reports).
//!
//! The compute commands are reached through [`mokume_cli::run_from_args`], the
//! same clap parsing + dispatch the standalone CLI binary uses, so the flag ->
//! config mapping stays single-sourced. The Python layer in `rust/python/mokume/`
//! builds the argument vector from ergonomic keyword arguments.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

/// The mokume version string (the crate's compile-time version).
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Run a mokume command in-process from an argument vector (no subprocess).
///
/// `args` is the subcommand and its flags WITHOUT the program name, for example
/// `["features2proteins", "--parquet", "x.parquet", "--output", "y.csv"]`. A
/// parse or runtime error is raised as a Python `RuntimeError` rather than
/// exiting the interpreter.
#[pyfunction]
fn run(args: Vec<String>) -> PyResult<()> {
    // clap's `parse_from` expects argv[0] to be the program name.
    let mut argv = Vec::with_capacity(args.len() + 1);
    argv.push("mokume".to_string());
    argv.extend(args);
    mokume_cli::run_from_args(argv).map_err(|error| PyRuntimeError::new_err(error.to_string()))
}

/// Run a mokume command as a CLI and return the process exit code.
///
/// This backs the `mokume` console script (`mokume.__main__:main`). Unlike
/// [`run`], it never raises: clap help/version are printed to stdout (exit code
/// 0), usage errors to stderr (exit code 2), and a runtime failure prints the
/// error and returns 1 — exactly as the standalone binary would, but without
/// exiting the process, so the caller decides via `sys.exit`.
#[pyfunction]
fn run_cli(args: Vec<String>) -> i32 {
    let mut argv = Vec::with_capacity(args.len() + 1);
    argv.push("mokume".to_string());
    argv.extend(args);
    mokume_cli::run_cli_from_args(argv)
}

/// The `mokume._mokume` extension module.
#[pymodule]
fn _mokume(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(run, module)?)?;
    module.add_function(wrap_pyfunction!(run_cli, module)?)?;
    Ok(())
}
