---
name: proteomics-heuristics
version: 0.1.0
description: >
  Domain heuristics for proteomics differential expression analysis.
  Covers data-type priors (LFQ, TMT, DIA), setting-specific top methods,
  ensemble strategies, and conditional rules derived from a 34576-experiment
  benchmark. Used by proposer and reflector to narrow the search space.
---

# Proteomics Heuristics

Domain knowledge distilled from large-scale benchmarking of proteomics
differential expression analysis workflows.

## Data-Type Priors

### LFQ

| Axis | Recommended |
|------|-------------|
| DE | deqms, limma, rots, proda |
| FDR | ihw, bh |
| Imputation | none, mindet, seqknn, missforest |
| Normalization | none, rlr, median |

DEqMS is best when peptide counts are available; limma is a stable baseline.

### TMT

| Axis | Recommended |
|------|-------------|
| DE | limma, proda, deqms |
| FDR | bh |
| Imputation | none, seqknn |
| Normalization | median, none |

TMT has lower missingness; limma is a strong baseline for TMT.

### DIA

| Axis | Recommended |
|------|-------------|
| DE | limma, rots, proda, deqms |
| FDR | ihw |
| Imputation | none, mindet, knn |
| Normalization | none, quantile, rlr |

DIA missing data is MAR not MNAR; avoid MinProb.

## Setting-Specific Top Methods

### FG_DDA (FragPipe DDA)

- **Best DE**: deqms, rots, limma
- **Best normalization**: none, median
- **Best imputation**: mindet, seqknn
- **Avoid DE**: none
- DEqMS dominates DDA when peptide counts are available.

### MQ_DDA (MaxQuant DDA)

- **Best DE**: deqms, limma
- **Best normalization**: none, median
- **Best imputation**: mindet, seqknn
- **Avoid DE**: none

### DIANN_DIA (DIA-NN DIA)

- **Best DE**: limma, rots
- **Best normalization**: none
- **Best imputation**: mindet
- **Avoid DE**: none
- DIA-NN directLFQ + no norm + MinDet + limma is the top workflow.

### spt_DIA (Spectronaut DIA)

- **Best DE**: rots, limma
- **Best normalization**: none
- **Best imputation**: mindet, knn
- ROTS is dominant on Spectronaut DIA (203/315 top-5 workflows).

### FG_TMT (FragPipe TMT)

- **Best DE**: limma, proda
- **Best normalization**: median, none
- **Best imputation**: none, seqknn
- limma is the top DEA for FragPipe TMT.

### MQ_TMT (MaxQuant TMT)

- **Best DE**: proda, limma
- **Best normalization**: median, none
- **Best imputation**: none, seqknn
- proDA is dominant on MaxQuant TMT (best_rank 28.8).

## Ensemble Strategies

### top_k_2of3

- **Methods**: limma, deqms, proda
- **min_k**: 2
- Conservative: protein must be DE in >=2/3 methods.

### top_k_3of4

- **Methods**: limma, rots, deqms, proda
- **min_k**: 3
- Strict: protein must be DE in >=3/4 methods.

## Condition Rules

| Condition | Action | Confidence | Evidence |
|-----------|--------|------------|----------|
| missing_rate < 0.15 | prefer imputation=none | high | No high-performance workflow uses imputation when missing < 15% |
| has_peptide_counts and n_proteins > 200 | prefer de=deqms | high | DEqMS top-1 for FG_DDA (92 top-5) and MQ_DDA (81 top-5) |
| n_per_group < 3 | prefer de=rots or de=limma | medium | ROTS bootstrap handles small samples; limma EB shrinkage is robust |
| n_proteins > 500 | prefer fdr=ihw | high | IHW gains power over BH when many hypotheses are available |
| data_type == DIA | avoid imputation=minprob | high | DIA missing is MAR; MNAR imputation (MinProb) hurts performance |
| data_type == TMT and missing_rate < 0.05 | prefer imputation=none and normalization=median | high | TMT has low missingness; median norm + no imputation is top |
