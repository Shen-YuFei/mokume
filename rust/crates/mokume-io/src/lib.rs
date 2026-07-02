pub mod peptide_parquet;
pub mod qpx;
pub mod sdrf;

pub use arrow::record_batch::RecordBatch;
pub use peptide_parquet::{
    read_peptide_parquet, write_peptide_parquet, PeptideParquetRow, RawPeptideRow, RawPeptideTable,
};
pub use qpx::{flatten_qpx_batch, QpxFeatureRecord, QpxIntensityEntry, QpxParquetReader};
pub use sdrf::{normalize_file_key, normalize_label_key, SdrfRawTable, SdrfRecord, SdrfTable};

pub const DEFAULT_QPX_BATCH_SIZE: usize = 65_536;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QpxReaderConfig {
    pub memory_bytes: Option<u64>,
    pub threads: Option<usize>,
    pub batch_size: usize,
}

impl QpxReaderConfig {
    pub const fn new(memory_bytes: Option<u64>, threads: Option<usize>) -> Self {
        Self {
            memory_bytes,
            threads,
            batch_size: DEFAULT_QPX_BATCH_SIZE,
        }
    }

    pub const fn with_batch_size(mut self, batch_size: usize) -> Self {
        self.batch_size = batch_size;
        self
    }
}

impl Default for QpxReaderConfig {
    fn default() -> Self {
        Self::new(None, None)
    }
}

pub fn crate_ready() -> bool {
    true
}
