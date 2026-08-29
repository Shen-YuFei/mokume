# Full Pipeline

The pure-Python package (`pip install mokume-py`) exposes an object-oriented pipeline
API. You build one `PipelineConfig` — grouping input, filtering, normalization,
quantification, IRS, batch correction, imputation, and DE into nested dataclasses
— and run it through `run_pipeline` or `QuantificationPipeline`. The result is a
`QpxDataset`: a hierarchical container with the protein matrix, intermediate
levels (peptides / PSMs), sample metadata, and a provenance log.

This is the surface to reach for when you want to script multi-stage analyses,
inspect intermediate levels, or export AnnData — rather than shelling out to the
CLI once per step.

## Building and running a config

`PipelineConfig` takes an `InputConfig` and a set of optional stage configs; each
unspecified stage uses its defaults. Below we quantify with MaxLFQ, apply
median run-normalization and global-median sample-normalization, and turn on KNN
imputation.

=== "Python (package)"

    ```python
    import warnings
    from mokume.pipeline.runner import run_pipeline
    from mokume.pipeline.config import (
        PipelineConfig,
        InputConfig,
        FilterConfig,
        NormalizationConfig,
        QuantificationConfig,
        ImputationConfig,
        BatchCorrectionConfig,
        IRSConfig,
        DEConfig,
    )

    warnings.filterwarnings("ignore")

    config = PipelineConfig(
        input=InputConfig(
            parquet="python/tests/example/feature_wide.parquet",
            sdrf="python/tests/example/PXD020192.sdrf.tsv",
        ),
        filtering=FilterConfig(min_aa=7, min_unique_peptides=2),
        normalization=NormalizationConfig(
            run_method="median",
            sample_method="globalMedian",
        ),
        quantification=QuantificationConfig(method="maxlfq"),
        # Each of these stages is off by default; enable the ones you need.
        imputation=ImputationConfig(enabled=True, method="knn", n_neighbors=3),
        batch=BatchCorrectionConfig(enabled=False),
        irs=IRSConfig(enabled=False),
        de=DEConfig(enabled=False),
    )

    dataset = run_pipeline(config)
    print(type(dataset).__name__)                 # QpxDataset
    print(dataset.get_level("proteins").shape)    # (6, 11) on the fixture
    ```

### Turning on the optional stages

Each stage config has an `enabled` flag plus its own parameters. Flip them on as
your experiment requires:

| Stage config | Enable with | Key parameters |
|--------------|-------------|----------------|
| `ImputationConfig` | `enabled=True` | `method` (`knn`, `minprob`, `missforest`, ...), `n_neighbors`, `quantile`, `shift`, `scale` |
| `BatchCorrectionConfig` | `enabled=True` | `method` (`sample_prefix`, `run`, `column`), `column`, `covariates`, `parametric` |
| `IRSConfig` | `enabled=True` | `reference_regex`, `sdrf_column`, `sdrf_values`, `stat`, `remove_reference` |
| `DEConfig` | `enabled=True` | `contrasts`, `method`, `log2fc_threshold`, `fdr_threshold`, `output` |

For example, to add IRS normalization and differential expression:

=== "Python (package)"

    ```python
    config.irs = IRSConfig(enabled=True, remove_reference=True, stat="median")
    config.de = DEConfig(
        enabled=True,
        contrasts=[("groupA", "groupB")],
        method="limma",
    )
    ```

    (DE needs at least two conditions in the data — see the
    [Differential Expression](differential-expression.md) page for how the shipped
    single-condition fixture constrains this.)

## `QuantificationPipeline.run_dataset`

`run_pipeline` is a functional wrapper; the class form gives you the same result
plus the assembled intermediate levels. `run()` returns just the wide protein
DataFrame, while `run_dataset()` packages everything into a `QpxDataset` with
provenance.

=== "Python (package)"

    ```python
    from mokume.pipeline.features_to_proteins import QuantificationPipeline

    pipeline = QuantificationPipeline(config)

    proteins_df = pipeline.run()          # -> pandas DataFrame (proteins x samples)
    dataset = pipeline.run_dataset()      # -> QpxDataset (proteins + peptides + ...)
    ```

## Working with the `QpxDataset`

The dataset materializes each level on demand and can pivot or export it.

=== "Python (package)"

    ```python
    dataset = QuantificationPipeline(config).run_dataset()

    # get_level(...) returns a DataFrame for one level.
    proteins = dataset.get_level("proteins")      # wide: ProteinName + samples
    peptides = dataset.get_level("peptides")      # long: one row per peptide/sample

    # to_wide_matrix(...) pivots a *long* level to protein x sample. The peptide
    # level is long, so give it the value/index/column names:
    peptide_matrix = dataset.to_wide_matrix(
        level="peptides",
        value_col="NormIntensity",
        protein_col="ProteinName",
        sample_col="SampleID",
    )
    print(peptide_matrix.shape)                    # (11, 10) on the fixture

    # to_anndata(...) exports a level as an AnnData (samples x proteins). The
    # protein level is already wide, so the defaults work:
    adata = dataset.to_anndata(level="proteins")
    print(adata.shape)                             # (10, 6) — 10 samples, 6 proteins
    ```

!!! note "The protein level is already wide"

    `run_pipeline` / `run_dataset` store `proteins` as a wide matrix
    (`ProteinName` first, one column per sample). Call `get_level("proteins")` to
    use it directly, or `to_anndata(level="proteins")` to export it — you do not
    need `to_wide_matrix` on it. Use `to_wide_matrix` for the *long* levels such as
    `peptides`.

## Runtime resources

`RuntimeConfig` controls the DuckDB resource hints used by the pure-Python
pipeline. It does not select a computation implementation.

=== "Python (package)"

    ```python
    from mokume.pipeline.config import RuntimeConfig

    config.runtime = RuntimeConfig(duckdb_memory="80GB", duckdb_threads=24)
    dataset = run_pipeline(config)
    ```

## Choosing the computation package

`run_pipeline` belongs to the pure-Python `mokume-py` distribution and always uses
its Python implementation. To use the Rust implementation, install `mokume`
in a separate environment and call its Rust-backed API directly:

=== "Python (wheel)"

    ```python
    import mokume

    mokume.features2proteins(
        parquet="features.parquet",
        sdrf="experiment.sdrf.tsv",
        output="proteins.csv",
        quant_method="maxlfq",
        threads=24,
    )
    ```

!!! warning "Install one distribution per environment"

    `mokume` and `mokume-py` both provide the `mokume` import package, so they
    overwrite each other's files when installed together. Select the
    distribution at installation time rather than through `RuntimeConfig`.

## What's next

- [Rust Wheel](../rust-wheel.md) — use the Rust kernel from Python or the
  installed console command.
- [Differential Expression](differential-expression.md) — run DE on the
  `QpxDataset` or hand its protein matrix to the separate Mokume Plugin.
- [Python API (package)](../reference/python-api-package.md) — the full generated
  reference for `PipelineConfig`, `QpxDataset`, and the pipeline classes.
