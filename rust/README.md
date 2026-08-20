# mokume

[![Rust](https://github.com/bigbio/mokume/actions/workflows/rust.yml/badge.svg)](https://github.com/bigbio/mokume/actions/workflows/rust.yml)
[![Wheels](https://github.com/bigbio/mokume/actions/workflows/wheels.yml/badge.svg)](https://github.com/bigbio/mokume/actions/workflows/wheels.yml)
[![PyPI version](https://badge.fury.io/py/mokume.svg)](https://badge.fury.io/py/mokume)
![PyPI - Downloads](https://img.shields.io/pypi/dm/mokume)

**mokume** is a proteomics quantification toolkit: it turns peptide-level mass
spectrometry intensities into protein expression matrices, with built-in
normalization, imputation, batch correction, and differential expression. It is
designed for the [quantms](https://github.com/bigbio/quantms) ecosystem and works
through both Python and an installed console command.

It supports piBAQ, TopN, MaxLFQ, and DirectLFQ quantification and ships as the
`mokume` wheel.

## Why mokume?

*Mokume-gane* (木目金, "wood-grain metal") is a Japanese metalworking technique
that fuses layers of different metals into a single piece with a distinctive
flowing pattern. mokume does the same with proteomics data: it melds many noisy,
overlapping peptide intensities into one coherent protein expression profile.

mokume is an evolution of [ibaqpy](https://github.com/bigbio/ibaqpy), extended
well beyond the original iBAQ workflow to a broader range of quantification,
normalization, and
differential-expression methods, with the heavy lifting rewritten in Rust.

## Installation

The Python wheel is the recommended way to install mokume:

```bash
pip install mokume
```

> **Distribution transition:** `mokume<=0.1.0` was pure Python. Starting with
> 0.2.0, this Rust-backed wheel is the default `mokume` distribution; the
> pure-Python implementation is available as `mokume-py`.

Optional periphery features (plots, reports, tissue maps, piBAQ helpers) are
available as extras:

```bash
pip install "mokume[plotting]"   # matplotlib / seaborn figures
pip install "mokume[reports]"    # interactive plotly QC reports
pip install "mokume[tissuemap]"  # tissue-map + AnnData export
pip install "mokume[pibaq]"      # Python piBAQ fallback for unported enzymes
pip install "mokume[analysis]"   # QC/comparison reports + missforest
pip install "mokume[agentic]"    # local MCP service for the Mokume Plugin
pip install "mokume[all]"        # everything above
```

The optional Plugin/MCP workflow requires Python 3.10 or newer.

The wheel also installs the `mokume` console command. A conda/mamba environment
is provided in `environment.yaml`. See
[Installation](../docs/installation.md) for details.

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

The Python API and console command run the identical compute kernel, so a result
produced through either interface is the same.

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
[user guide](../docs/user-guide/).

`features2proteins` runs these stages in order:

<p align="center">
  <img src="../docs/assets/pipeline.svg" alt="The mokume features2proteins pipeline: source data through quantify, normalize, impute, batch-correct, and differential expression, with the best-known methods at each stage" width="100%">
</p>

- **Quantification:** `maxlfq`, `directlfq`, `pibaq`, `top3`/`topn`, `sum`
  (also `median`, `ratio`, `abd`, `intensity`, `spectral_count`). piBAQ requires
  a FASTA; TopN, MaxLFQ, and Sum do not.
- **Normalization:** run-level and sample-level options including `median`,
  `quantile`, `rlr`, and `loess`.
- **Imputation:** a focused set of imputers, from simple (`mindet`, `knn`) to
  model-based (`qrilc`, `impseq`).
- **Batch correction:** native ComBat (parametric, non-parametric, and
  covariate-aware) that removes technical batch effects while preserving
  biological signal.
- **Differential expression:** `limma`, `deqms`, `rots`, `limrots`,
  `proda`, and an ensemble, with BH (default), IHW, BKY, or Storey FDR control
  and fixed or data-driven effect-size gates. LimROTS gives
  the best sensitivity on MaxLFQ data, DEqMS controls false positives better on
  noisier DirectLFQ data, and proDA models missing values as informative
  dropouts.

## How it works

mokume's compute lives once in a Rust kernel shipped in the `mokume` wheel.
The compiled extension runs in-process, and the installed `mokume` console
command reaches that extension through Python without a subprocess.

A thin Python periphery adds plots and reports on top. You never need to know it
is Rust to use it. The repository also contains the separately maintained
pure-Python `mokume-py` distribution. The Rust kernel leads new native
computation work. For the full design and maintenance policy, see
[Architecture](../docs/architecture.md) and
[Maintenance scope](../docs/maintenance-scope.md).

Agentic method recommendation is an installable plugin over this wheel. Its
bundled local MCP configuration calls deterministic profile, policy, and
evaluation services; the agent host owns the model and credentials. See the
[Mokume Plugin guide](../docs/user-guide/agentic-plugin.md).

## Example

A full run on real data: **PXD030304**, a 949-cell-line proteomic panel
(178 M feature rows). mokume's tissue-proteome command (`mokume.tissuemap`)
quantifies the cell lines, scores AdaTiSS tissue specificity, embeds the
samples, and finds tissue markers. Every figure below is rendered by mokume's
**own** visualization on the Rust build — no external plotting code:

<p align="center">
  <img src="../docs/assets/pxd030304_tissue_atlas.png" alt="Tissue atlas: the cell-line proteomes embedded and grouped by tissue of origin" width="100%">
</p>

*Tissue atlas — the 949 cell-line proteomes embedded and grouped by their tissue
of origin (`mokume.tissuemap`).*

<p align="center">
  <img src="../docs/assets/pxd030304_marker_tsne.png" alt="t-SNE panels of the cell-line proteomes coloured by top tissue-marker expression" width="100%">
</p>

*t-SNE of the same proteomes, each panel coloured by a top tissue-marker's
expression.*

<p align="center">
  <img src="../docs/assets/pxd030304_ts_distribution.png" alt="AdaTiSS tissue-specificity score distribution and tissue-specific protein counts" width="100%">
</p>

*AdaTiSS tissue-specificity score distribution (with the GMM-fitted
specific / enriched / housekeeping thresholds) and the tissue-specific protein
count per tissue. `mokume.tissuemap` also emits marker heatmaps, dotplots, and
dendrograms — see [`docs/periphery/tissuemap.md`](../docs/periphery/tissuemap.md).*

## Documentation

- [Documentation home](../docs/index.md)
- [Quick start](../docs/quickstart.md)
- [Installation](../docs/installation.md)
- [Rust wheel](../docs/rust-wheel.md)
- [Architecture](../docs/architecture.md)
- [Maintenance scope](../docs/maintenance-scope.md)
- [Benchmarks](../benchmarks/)

## Citation

mokume is part of the [quantms](https://github.com/bigbio/quantms) ecosystem and
evolves [ibaqpy](https://github.com/bigbio/ibaqpy). Until a dedicated mokume
paper is available, please cite ibaqpy and quantms — see those repositories for
their current citation details and DOIs.

## Credits, contributing, and license

mokume is developed by the [bigbio](https://github.com/bigbio) community as part
of the quantms ecosystem. Contributions are welcome; see
[the community guide](../docs/community.md) for development setup and guidelines.

Licensed under the [MIT License](../LICENSE).
