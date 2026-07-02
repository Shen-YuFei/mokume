use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum MokumeError {
    #[error("input file does not exist: {path}")]
    MissingInput { path: PathBuf },

    #[error("I/O error for {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("invalid input: {message}")]
    InvalidInput { message: String },

    #[error("invalid memory value `{value}`")]
    InvalidMemory { value: String },

    #[error("stage `{stage}` is not implemented yet")]
    NotImplemented { stage: &'static str },
}

pub type Result<T> = std::result::Result<T, MokumeError>;
