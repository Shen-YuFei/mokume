# Installation

## From PyPI

```bash
pip install mokume
```

For standard quantification workflows, the core package is enough.
Install optional extras only when you need specific functionality.

### Optional Extras

mokume uses optional dependencies for specialized features:

=== "DirectLFQ"

    ```bash
    pip install mokume[directlfq]
    ```

    Enables DirectLFQ quantification and the DirectLFQ backend for MaxLFQ.

    The extra also pulls in `polars`, which mokume uses to stream the
    long-form parquet into the wide DirectLFQ matrix without materialising
    an intermediate pandas DataFrame. On large studies (>1000 samples) this
    cuts the load step's wall time roughly in half and lets datasets that
    previously OOM-killed pandas pivots (e.g. PXD030304 at 5798 samples)
    complete on a 125 GB host.

=== "Plotting"

    ```bash
    pip install mokume[plotting]
    ```

    Enables volcano plots, heatmaps, PCA, and box plots (matplotlib + seaborn).

=== "Batch Correction"

    ```bash
    pip install mokume[batch-correction]
    ```

    Enables ComBat-based batch correction via the `combat` dependency.

=== "Interactive Reports"

    ```bash
    pip install mokume[reports]
    ```

    Enables interactive HTML reports with plotly.

=== "TissueMap"

    ```bash
    pip install mokume[tissuemap]
    ```

    Enables the `mokume tissuemap` workflow for per-dataset tissue atlas analysis,
    including AdaTiSS tissue-specificity scoring, embeddings, and atlas plots.

=== "Everything"

    ```bash
    pip install mokume[all]
    ```

    Installs all optional dependencies.

## From Source

```bash
git clone https://github.com/bigbio/mokume
cd mokume
pip install .
```

## Using Conda

```bash
mamba env create -f environment.yaml
conda activate mokume
pip install .
```

## Requirements

- Python >= 3.9
- Core dependencies: numpy, pandas, scipy, scikit-learn, pyopenms, pyarrow, duckdb, click, anndata
- DE methods (limma, DEqMS, proDA, LimROTS, ROTS), imputation (BPCA, impSeq, impSeqRob), and VSN normalization are **pure-Python** — no R or rpy2 required
