# mokume

[![Python application](https://github.com/bigbio/mokume/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/mokume/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/mokume/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/mokume/actions/workflows/python-publish.yml)
[![PyPI version](https://badge.fury.io/py/mokume.svg)](https://badge.fury.io/py/mokume)
![PyPI - Downloads](https://img.shields.io/pypi/dm/mokume)

**mokume** is a proteomics quantification toolkit: it turns peptide-level mass
spectrometry intensities into protein expression matrices, with built-in
normalization, imputation, batch correction, and differential expression. It is
designed for the [quantms](https://github.com/bigbio/quantms) ecosystem and works
equally well as a standalone command-line tool.

It supports iBAQ, TopN, MaxLFQ, and DirectLFQ quantification, and ships both as a
fast standalone CLI binary and as a `pip install mokume-rs` wheel.

## Why mokume?

*Mokume-gane* (木目金, "wood-grain metal") is a Japanese metalworking technique
that fuses layers of different metals into a single piece with a distinctive
flowing pattern. mokume does the same with proteomics data: it melds many noisy,
overlapping peptide intensities into one coherent protein expression profile.

mokume is an evolution of [ibaqpy](https://github.com/bigbio/ibaqpy), extended
well beyond iBAQ to a broader range of quantification, normalization, and
differential-expression methods, with the heavy lifting rewritten in Rust.

## Installation

The Python wheel is the recommended way to install mokume:

```bash
pip install mokume-rs
```

> **Note:** `mokume-rs` is not published on PyPI yet. Until it is, build the wheel
> from source (`pip install .` from this `rust/` directory, or `pip install ./rust`
> from the repo root), or use the pure-Python package: `pip install mokume`.

Optional periphery features (plots, reports, tissue maps, iBAQ helpers) are
available as extras:

```bash
pip install "mokume-rs[plotting]"   # matplotlib / seaborn figures
pip install "mokume-rs[reports]"    # interactive plotly QC reports
pip install "mokume-rs[tissuemap]"  # tissue-map + AnnData export
pip install "mokume-rs[ibaq]"       # iBAQ / ProteomicRuler (needs a FASTA)
pip install "mokume-rs[all]"        # everything above
```

To install the standalone CLI binary (no Python runtime required):

```bash
cargo install --path crates/mokume-cli
```

A conda/mamba environment is also provided in `environment.yaml`. See
[docs/installation.md](docs/installation.md) for details.

## Quick start

Run the full pipeline from a quantms feature table to a protein matrix:

```bash
mokume features2proteins \
  --parquet features.parquet \
  --sdrf samples.sdrf.tsv \
  --quant-method maxlfq \
  --output proteins.csv
```

Add differential expression by passing `--de` with one or more contrasts:

```bash
mokume features2proteins \
  --parquet features.parquet \
  --sdrf samples.sdrf.tsv \
  --quant-method maxlfq \
  --de --de-contrasts "Treatment-Control" \
  --output proteins.csv
```

The same kernel runs in-process from Python. Each keyword argument maps to a CLI
flag, and the command writes its result to `--output`:

```python
import mokume

mokume.features2proteins(
    parquet="features.parquet",
    sdrf="samples.sdrf.tsv",
    quant_method="maxlfq",
    output="proteins.csv",
)
```

The wheel and the CLI binary run the identical compute kernel, so a result
produced either way is the same.

## Commands

| Command | What it does |
| --- | --- |
| `features2proteins` | Full pipeline: feature table to a protein quantification matrix |
| `features2peptides` | Aggregate features to peptide-level intensities |
| `peptides2protein` | Roll peptide intensities up to protein quantities |
| `correct-batches` | Standalone ComBat batch correction (with AnnData export) |

## Quantification and methods

mokume covers the common protein quantification strategies and the analysis
steps around them. The full catalog of methods is documented in the
[user guide](docs/user-guide/).

`features2proteins` runs these stages in order:

<p align="center">
  <img src="docs/assets/pipeline.svg" alt="The mokume features2proteins pipeline: source data through quantify, normalize, impute, batch-correct, and differential expression, with the best-known methods at each stage" width="100%">
</p>

- **Quantification:** `maxlfq`, `directlfq`, `ibaq`, `top3`/`topn`, `sum`
  (also `median`, `ratio`, `abd`, `intensity`, `spectral_count`). iBAQ requires
  a FASTA; TopN, MaxLFQ, and Sum do not.
- **Normalization:** run-level and sample-level options including `median`,
  `quantile`, `rlr`, and `loess`.
- **Imputation:** a focused set of imputers, from simple (`mindet`, `knn`) to
  model-based (`qrilc`, `impseq`).
- **Batch correction:** native ComBat (parametric, non-parametric, and
  covariate-aware) that removes technical batch effects while preserving
  biological signal.
- **Differential expression:** `limma`, `deqms`, `rots`, `limrots`,
  `proda`, and an ensemble, with BH (default) or IHW FDR control. LimROTS gives
  the best sensitivity on MaxLFQ data, DEqMS controls false positives better on
  noisier DirectLFQ data, and proDA models missing values as informative
  dropouts.

## How it works

mokume's compute lives once in a Rust kernel and ships two ways:

- a standalone CLI binary (`mokume`) that needs no Python runtime, and
- a `pip install mokume-rs` wheel whose compiled extension runs the same kernel
  in-process (no subprocess, no recomputation).

A thin Python periphery adds plots and reports on top. You never need to know it
is Rust to use it. For the full design, see
[docs/architecture.md](docs/architecture.md).

## Example

A full run on real data: **PXD030304**, a 949-cell-line proteomic panel
(178 M feature rows). mokume's tissue-proteome command (`mokume.tissuemap`)
quantifies the cell lines, scores AdaTiSS tissue specificity, embeds the
samples, and finds tissue markers. Every figure below is rendered by mokume's
**own** visualization on the Rust build — no external plotting code:

<p align="center">
  <img src="docs/assets/pxd030304_tissue_atlas.png" alt="Tissue atlas: the cell-line proteomes embedded and grouped by tissue of origin" width="100%">
</p>

*Tissue atlas — the 949 cell-line proteomes embedded and grouped by their tissue
of origin (`mokume.tissuemap`).*

<p align="center">
  <img src="docs/assets/pxd030304_marker_tsne.png" alt="t-SNE panels of the cell-line proteomes coloured by top tissue-marker expression" width="100%">
</p>

*t-SNE of the same proteomes, each panel coloured by a top tissue-marker's
expression.*

<p align="center">
  <img src="docs/assets/pxd030304_ts_distribution.png" alt="AdaTiSS tissue-specificity score distribution and tissue-specific protein counts" width="100%">
</p>

*AdaTiSS tissue-specificity score distribution (with the GMM-fitted
specific / enriched / housekeeping thresholds) and the tissue-specific protein
count per tissue. `mokume.tissuemap` also emits marker heatmaps, dotplots, and
dendrograms — see [`docs/periphery/tissuemap.md`](docs/periphery/tissuemap.md).*

## Documentation

- [Documentation home](docs/index.md)
- [Quick start](docs/quickstart.md)
- [Installation](docs/installation.md)
- [CLI vs. wheel](docs/cli-vs-wheel.md)
- [Architecture](docs/architecture.md)
- [Benchmarks](benchmarks/)

## Citation

mokume is part of the [quantms](https://github.com/bigbio/quantms) ecosystem and
evolves [ibaqpy](https://github.com/bigbio/ibaqpy). Until a dedicated mokume
paper is available, please cite ibaqpy and quantms — see those repositories for
their current citation details and DOIs.

## Credits, contributing, and license

mokume is developed by the [bigbio](https://github.com/bigbio) community as part
of the quantms ecosystem. Contributions are welcome; see
[the community guide](docs/community.md) for development setup and guidelines.

Licensed under the [MIT License](LICENSE).
