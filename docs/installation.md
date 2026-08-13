# Installation

## From PyPI

```bash
pip install mokume-rs
```

!!! warning "`mokume-rs` is not on PyPI yet"

    The Rust wheel has not been published to PyPI yet, so `pip install mokume-rs`
    does not work today. Until it ships, either install the pure-Python package
    (`pip install mokume` — the same import name, but a separately maintained API)
    or build the wheel from source: `pip install ./rust` (needs the Rust toolchain
    and cmake). The `mokume-rs` instructions on this page describe the wheel for
    when it is released.

!!! warning "Install `mokume-rs` **or** `mokume` — never both"

    Both distributions install the same `mokume` import package, so having both
    in one environment makes pip silently overwrite files (it does not detect the
    collision). Keep only one; to switch, `pip uninstall` the other first. Each
    package warns at import time if it finds its sibling installed.

`pip install mokume-rs` gives you the Rust compute kernel (the compiled
`mokume._mokume` extension) and the thin Python API that drives it in-process.
The kernel needs no third-party Python dependencies, so for standard
quantification workflows the core package is enough. Extras pull in only the
Python periphery libraries; install one when you need a specific periphery
command.

!!! note "Standalone CLI binary"

    The kernel also ships as a standalone CLI binary `mokume`, built from
    `rust/crates/mokume-cli` with cargo and requiring no Python at all. The wheel and
    the binary expose the same four compute subcommands over the same kernel.

### Optional Extras

mokume uses optional dependencies for the periphery commands:

=== "Plotting"

    ```bash
    pip install mokume-rs[plotting]
    ```

    Enables the t-SNE visualization, DE plots, and iBAQ QC report periphery
    commands (numpy, pandas, scipy, scikit-learn, matplotlib, seaborn).

=== "Interactive Reports"

    ```bash
    pip install mokume-rs[reports]
    ```

    Enables the interactive HTML report periphery command (plotly).

=== "TissueMap"

    ```bash
    pip install mokume-rs[tissuemap]
    ```

    Enables the `mokume.tissuemap` periphery command for per-dataset tissue atlas
    analysis, including AdaTiSS tissue-specificity scoring, embeddings, and atlas
    plots (scanpy, anndata, umap-learn, combat, matplotlib, seaborn, pyarrow).

=== "iBAQ"

    ```bash
    pip install mokume-rs[ibaq]
    ```

    Enables the pure-Python `mokume.peptides2protein_ibaq` fallback for enzymes
    the Rust kernel does not digest (pyopenms, pyarrow, PyYAML, numpy, pandas,
    scipy).

=== "Analysis"

    ```bash
    pip install mokume-rs[analysis]
    ```

    Enables the QC / workflow-comparison reports and the pure-Python method
    fallback the Rust kernel does not reproduce — `mokume.impute(method="missforest")`
    — plus `mokume.qc_report`
    and `mokume.workflow_comparison` (numpy, pandas, scipy, scikit-learn).

=== "Everything"

    ```bash
    pip install mokume-rs[all]
    ```

    Installs all optional periphery dependencies.

## From Source

mokume builds from two project roots, one per distribution:

- **`mokume` (pure Python)** — `python/pyproject.toml` uses standard PEP 621
  metadata with the **hatchling** build backend
  (`build-backend = "hatchling.build"`, not Poetry). Install it from source with:

  ```bash
  git clone https://github.com/bigbio/mokume
  cd mokume
  pip install ./python          # builds the mokume package via the hatchling backend
  ```

- **`mokume-rs` (Rust wheel)** — the wheel is built with maturin, which compiles
  the `mokume-py` PyO3 crate into the `mokume._mokume` extension and packages it
  with the Python periphery. The maturin project lives in `rust/`:

  ```bash
  pip install ./rust            # builds the extension via the maturin backend
  ```

For a development checkout of the Rust wheel, build the extension in place with
`maturin develop` run from `rust/`; the periphery is plain Python and needs no
build step. The standalone CLI binary is built from the same workspace:

```bash
cargo build --release --manifest-path rust/Cargo.toml -p mokume-cli
```

## Using Conda

```bash
mamba env create -f rust/environment.yaml
conda activate mokume
pip install ./rust
```

## Requirements

- Python >= 3.9 (for the wheel) — the standalone CLI binary needs no Python
- The compute kernel needs **no third-party Python dependencies**; the periphery
  extras pull in numpy, pandas, scipy, scikit-learn, pyopenms, pyarrow, plotly,
  scanpy, and friends only for the periphery commands you run
- The quantification methods (DirectLFQ, MaxLFQ, iBAQ), ComBat batch correction,
  DE methods (limma, DEqMS, proDA, LimROTS, ROTS), and most imputation methods
  run in the **native Rust kernel** — no R or rpy2 required. ComBat is
  oracle-verified against inmoose
- A method stays in the Python periphery because the kernel cannot reproduce
  it cross-language: the `missforest` imputer (`mokume.impute`), in the
  `analysis` extra
