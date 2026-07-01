---
name: proteomics-benchmarks
version: 0.1.0
description: >
  Reference benchmark datasets for proteomics DEA optimization.
  Contains ground-truth metadata, expected fold-changes, and best-known
  method configurations. Used by the evaluator to score candidate
  workflows against known references.
---

# Proteomics Benchmarks

Reference datasets with ground truth for evaluating differential expression
analysis workflows.

## PXD001819

- **Description**: UPS1 spike-in, 27 runs, 9 conditions, LFQ, Orbitrap Velos
- **Ground-truth proteins**: ups1_48_proteins.txt
- **Expected log2FC**:

| Contrast | log2FC |
|----------|--------|
| A_vs_B | 2.0 |
| A_vs_C | 4.32 |
| A_vs_D | 7.64 |
| B_vs_C | 2.32 |
| B_vs_D | 5.64 |
| C_vs_D | 3.32 |

- **Best-known configuration**:
  - **Method**: deqms
  - **FDR method**: ihw
  - **Imputation**: none
  - **Metrics**: AUC 0.996, TP 47, FP 1
