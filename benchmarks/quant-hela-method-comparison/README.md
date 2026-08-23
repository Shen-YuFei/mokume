# HeLa and Human Protein Quantification Benchmark

This benchmark recomputes six protein-quantification methods with the current
Rust-backed `mokume` distribution. It uses 20 public human QPX datasets for
cross-experiment comparisons and the paired PXD007683 TMT/LFQ inputs for a
cross-technology check.

## Scope and interpretation

- Methods: piBAQ, MaxLFQ, DirectLFQ, Top3, Top10, and Sum.
- piBAQ uses the current Rust core and runtime pyOpenMS digestion. The three
  Nagaraj inputs use their experimental proteases: Trypsin, Lys-C, and Glu-C.
- Every within-dataset method comparison uses the same protein universe for
  that dataset. PXD013658.2 has only one sample and is retained for
  cross-experiment analysis but excluded from the CV calculation.
- CV here is descriptive sample dispersion. The datasets are heterogeneous
  public studies, so this metric must not be presented as a technical-replicate
  reproducibility ranking.
- No ground-truth abundance is available across the 20 studies. The benchmark
  therefore reports metric-specific behavior and does not declare a global
  winner.

## Results

### Cross-experiment consistency

Mean pairwise Pearson correlations across the 20 datasets were 0.7042 for
DirectLFQ, 0.6796 for MaxLFQ, 0.6472 for Sum, 0.6253 for Top3, 0.6221 for piBAQ,
and 0.6192 for Top10. Mean pairwise Spearman correlations were 0.7018, 0.6543,
0.6246, 0.5918, 0.6172, and 0.5726, respectively.

These values measure agreement of median protein profiles after pairwise
matching. They do not measure absolute accuracy and cannot by themselves select
an optimal method.

![Cross-experiment correlation heatmap](figures/correlation_heatmap.png)

### Within-dataset sample dispersion

The original boxplot chart type is retained. Across the 21 eligible datasets,
mean CV was 0.4287 for MaxLFQ, 0.5384 for DirectLFQ, 0.7905 for Top10, 0.8324
for Top3, 0.8351 for piBAQ, and 0.8920 for Sum.

![CV distribution](figures/cv_distribution.png)

Because study design, sample count, and biological heterogeneity differ among
datasets, these CVs are descriptive and should be read together with the
cross-experiment and coverage results.

### PXD007683 TMT/LFQ agreement

All six methods were compared on the same 5,312 proteins. Pearson/Spearman
correlations were 0.8446/0.8524 for piBAQ, 0.7837/0.7823 for MaxLFQ,
0.8226/0.8235 for DirectLFQ, 0.7717/0.7733 for Top3, 0.7740/0.7654 for Top10,
and 0.8245/0.8259 for Sum. This is a cross-technology agreement check, not a
ground-truth recovery score.

## Reproduction

Download the public inputs described in `scripts/config.py`:

```bash
python scripts/01_download_data.py
```

Then run the current Rust refresh. Keep the work directory outside the
repository if the protein matrices should remain disposable.

```bash
python scripts/refresh_rust.py \
  --raw-dir data/raw \
  --fasta /path/to/Homo-sapiens-uniprot-reviewed.fasta \
  --work-dir /tmp/mokume-hela \
  --threads 24 \
  --force
```

The tracked `results/` directory contains summary metrics rather than the large
protein matrices. The two existing PNG filenames and chart types are preserved.

## Data sources

The 20 QPX inputs come from the public PRIDE `ibaqpy-research` resource. The
PXD007683 MSstats inputs come from the public `quantms-benchmark-old` resource;
their exact URLs and filenames are recorded in `scripts/config.py`.
