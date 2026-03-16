# Installation

## From PyPI

```bash
pip install mokume
```

### Optional Extras

mokume uses optional dependencies for specialized features:

=== "DirectLFQ"

    ```bash
    pip install mokume[directlfq]
    ```

    Enables DirectLFQ quantification and the DirectLFQ backend for MaxLFQ.

=== "Plotting"

    ```bash
    pip install mokume[plotting]
    ```

    Enables volcano plots, heatmaps, PCA, and box plots (matplotlib + seaborn).

=== "Batch Correction"

    ```bash
    pip install mokume[batch-correction]
    ```

    Enables ComBat batch correction via the inmoose package.

=== "Interactive Reports"

    ```bash
    pip install mokume[reports]
    ```

    Enables interactive HTML reports with plotly.

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
