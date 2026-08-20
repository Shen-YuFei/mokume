# Mokume Benchmarks

Benchmarking studies for mokume protein quantification methods.

## Available Benchmarks

| Benchmark | Description | Key Finding |
|-----------|-------------|-------------|
| [quant-hela-method-comparison](./quant-hela-method-comparison/) | Cross-experiment consistency using HeLa datasets | iBAQ (log) best for cross-experiment (r=0.74) |
| [batch-quartet-multilab](./batch-quartet-multilab/) | Batch effect correction with Quartet multi-lab data | DirectLFQ + ComBat most reliable |
| [quant-pxd007683-tmt-vs-lfq](./quant-pxd007683-tmt-vs-lfq/) | TMT vs LFQ comparison from Gygi lab | `median-cov` gives lowest CV |

## Quick Summary

### quant-hela-method-comparison
Evaluates iBAQ, DirectLFQ, TopN, Sum across 20 HeLa/human datasets.
- **Winner:** iBAQ (log-transformed) with 0.74 cross-experiment correlation
- **Use for:** Cross-experiment comparisons, absolute quantification

### batch-quartet-multilab
Benchmarks batch correction on 6-lab, 72-sample Quartet data.
- **Winner:** DirectLFQ + ComBat
- **Note:** ~40% missing values, ~15% lab-specific variance persists

### quant-pxd007683-tmt-vs-lfq
Compares TMT and LFQ on same samples (PXD007683).
- **Winner:** `median-cov` normalization (lowest CV)
- **TMT vs LFQ:** TMT better for small fold-changes; both agree on relative abundances

## Running Benchmarks

```bash
pip install "mokume-py[analysis]"   # compute kernel + analysis extras
cd benchmarks/<benchmark-name>
python scripts/run_benchmark.py
```

Quantification calls go through the Rust kernel via `mokume.peptides2protein()` /
`mokume.features2proteins()` (file-based, in-process through PyO3). Normalization,
imputation, and filter benchmarks that test individual components use numpy / scipy
/ scikit-learn directly (the Rust kernel bundles these internally but does not
expose them as standalone Python APIs).

## Data Policy

Large data files (>50MB) are **NOT** committed to git. Each benchmark README documents data sources and download scripts.

## Benchmark Structure

```
benchmarks/<benchmark-name>/
├── README.md           # Overview, results, methodology (expandable)
├── scripts/            # Benchmark scripts
├── data/               # GIT-IGNORED: large data files
├── results/            # CSV metrics (baselines from the original mokume-py runs)
└── figures/            # PNG visualizations
```
