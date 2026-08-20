use std::{
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use mokume_cli::run_from_args;

#[test]
fn log_file_option_writes_file_and_preserves_command_errors(
) -> Result<(), Box<dyn std::error::Error>> {
    let root = temp_root()?;
    let log_file = root.join("logs").join("mokume.log");
    let error = match run_from_args([
        "mokume",
        "--log-file",
        path_str(&log_file)?,
        "features2proteins",
        "--parquet",
        "definitely-missing.feature.parquet",
        "--output",
        "protein.csv",
    ]) {
        Ok(()) => panic!("command must still fail for a missing input"),
        Err(error) => error,
    };

    let message = error.to_string();
    assert!(
        message.contains("input file does not exist: definitely-missing.feature.parquet"),
        "unexpected error: {message}"
    );
    assert!(
        !message.contains("stage `log-file` is not implemented yet"),
        "log-file option must not block command execution: {message}"
    );
    assert!(
        log_file.exists(),
        "log file was not created: {}",
        log_file.display()
    );
    Ok(())
}

fn temp_root() -> Result<PathBuf, Box<dyn std::error::Error>> {
    let timestamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    Ok(tempfile::Builder::new()
        .prefix(&format!("mokume-command-test-{timestamp}-"))
        .tempdir()?
        .keep())
}

fn path_str(path: &Path) -> Result<&str, Box<dyn std::error::Error>> {
    path.to_str()
        .ok_or_else(|| format!("path is not UTF-8: {}", path.display()).into())
}
