//! Statistical post-processing for mokume: differential expression
//! ([`de`]) and batch correction ([`batch`]).
//!
//! The numerical method implementations are added per file under `de/` and
//! `batch/` as the Rust port reaches parity with the Python `mokume.analysis`
//! and `mokume.postprocessing` packages. Public APIs are defined when the first
//! method lands, not speculatively ahead of a caller.

pub mod batch;
pub mod de;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixOrientation {
    ProteinsBySamples,
    SamplesByProteins,
}
