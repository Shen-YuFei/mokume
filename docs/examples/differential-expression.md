# Differential Expression

Test proteins for abundance changes between two conditions. mokume offers several
DE methods — `limma`, `deqms`, `proda`, `limrots`, `rots`, a top-k `ensemble`, and
`auto` (pick a method from the data shape) — reachable from three surfaces:

- the kernel's `features2proteins --de ...` flags (single-sourced Rust);
- the pure-Python `mokume.analysis.DifferentialExpression` class and the
  standalone `run_deqms` / `run_limma` / `run_limrots` / `run_proda` functions;
- the LLM-free agentic optimizer, which searches preprocessing + DE choices for
  you.

## (a) Kernel DE via `features2proteins --de`

The kernel quantifies, then runs DE per contrast and writes one CSV per contrast
(named from `--de-output`). A contrast is written `GroupA vs GroupB` (or
`GroupA-GroupB`); the group names must be **conditions that appear in your SDRF**.

=== "CLI"

    ```bash
    # One inline contrast, limma, one result CSV per contrast
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --de --de-contrasts "NASH vs HL" \
        --de-method limma \
        --de-fdr-method bh \
        --de-output de_results.csv

    # Many contrasts from a TSV file (columns: group1<TAB>group2)
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv \
        --de --de-contrasts-file contrasts.tsv \
        --de-method ensemble \
        --de-output de_results.csv
    ```

    The contrasts file is a tab-separated table with `group1` and `group2`
    header columns (each `<TAB>` below is a literal tab character):

    ```text
    group1<TAB>group2
    NASH<TAB>HL
    Steatosis<TAB>HL
    ```

=== "Python (wheel)"

    ```python
    import mokume

    # de_contrasts takes a list; each entry is "GroupA vs GroupB".
    mokume.features2proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
        de=True,
        de_contrasts=["NASH vs HL"],
        de_method="limma",
        de_fdr_method="bh",
        de_output="de_results.csv",
    )

    # Or point at a group1/group2 TSV instead of inline contrasts.
    mokume.features2proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        de=True,
        de_contrasts_file="contrasts.tsv",
        de_method="ensemble",
        de_output="de_results.csv",
    )
    ```

!!! note "The shipped fixture has only one condition"

    `PXD020192.sdrf.tsv` labels every sample as the single condition `Brain`, so a
    real two-group contrast is not possible from it — the kernel would stop with
    `No samples for 'NASH'. Available: ['Brain']`. The commands above use the
    same illustrative `NASH`/`HL` labels as the [Quick Start](../quickstart.md);
    swap in condition names from your own SDRF. To run DE end-to-end against the
    fixture, use the Python path below, where you assign the sample groups
    yourself.

The kernel writes one CSV per contrast. The leading columns are shared across
methods (`ProteinName`, `log2FC`, `pvalue`, `adj_pvalue`, `significance`); later
columns are method-specific. See
[Computed Values](../reference/computed-values.md#differential-expression-columns)
for the full column reference.

## (b) Python `DifferentialExpression` and the `run_*` functions

The pure-Python package (`pip install mokume`) exposes DE directly on an
in-memory protein matrix, so you control the sample-to-condition mapping. This
runs fully against the fixture.

=== "Python (package)"

    ```python
    import warnings
    from mokume.pipeline.features_to_proteins import QuantificationPipeline
    from mokume.pipeline.config import (
        PipelineConfig,
        InputConfig,
        QuantificationConfig,
    )
    from mokume.analysis import DifferentialExpression

    warnings.filterwarnings("ignore")

    # 1. Quantify to a wide protein matrix (ProteinName + one column per sample).
    config = PipelineConfig(
        input=InputConfig(
            parquet="python/tests/example/feature_wide.parquet",
            sdrf="python/tests/example/PXD020192.sdrf.tsv",
        ),
        quantification=QuantificationConfig(method="maxlfq"),
    )
    proteins = QuantificationPipeline(config).run()

    # 2. Assign each sample column to a condition. Here we split the 10 fixture
    #    samples into two demo groups; use your real conditions in practice.
    samples = [c for c in proteins.columns if c != proteins.columns[0]]
    half = len(samples) // 2
    sample_to_condition = {
        s: ("groupA" if i < half else "groupB")
        for i, s in enumerate(samples)
    }

    # 3. Run a single contrast with the method you want.
    de = DifferentialExpression(method="limma")
    result = de.run(proteins, sample_to_condition, ("groupA", "groupB"))
    print(result[["ProteinName", "log2FC", "pvalue", "adj_pvalue", "significance"]])
    ```

    On the fixture only one protein is fully observed in both groups, so the
    result has a single row:

    ```text
      ProteinName    log2FC    pvalue  adj_pvalue significance
           P09382   -0.8525    0.4226      0.4226    Unchanged
    ```

Swap `method="limma"` for `"deqms"`, `"proda"`, `"limrots"`, or `"rots"` to use a
different test. Each method is also callable as a standalone function that takes a
log2 protein matrix (proteins in the index, samples in the columns) and explicit
sample lists:

=== "Python (package)"

    ```python
    import numpy as np
    from mokume.analysis import run_limma, run_deqms, run_limrots, run_proda

    protein_col = proteins.columns[0]
    samples = [c for c in proteins.columns if c != protein_col]
    half = len(samples) // 2
    samples_a, samples_b = samples[:half], samples[half:]

    # The standalone functions expect a log2-transformed matrix.
    log2_matrix = np.log2(proteins.set_index(protein_col)[samples])

    # run_limma / run_proda take the condition names positionally:
    limma_res = run_limma(log2_matrix, samples_a, samples_b, "groupA", "groupB")
    proda_res = run_proda(log2_matrix, samples_a, samples_b, "groupA", "groupB")

    # run_deqms / run_limrots take the contrast as a (name_a, name_b) tuple:
    deqms_res = run_deqms(log2_matrix, samples_a, samples_b, ("groupA", "groupB"))
    limrots_res = run_limrots(log2_matrix, samples_a, samples_b, ("groupA", "groupB"))

    for name, res in [
        ("limma", limma_res),
        ("proda", proda_res),
        ("deqms", deqms_res),
        ("limrots", limrots_res),
    ]:
        print(name, res.shape)
    ```

    All four return a DataFrame whose first four columns are `ProteinName`,
    `log2FC`, `pvalue`, `adj_pvalue`; the remaining columns are method-specific
    statistics (`t_stat`/`AveExpr`/`B` for limma, `sca_t` for deqms, `d_stat` for
    rots-family methods).

!!! note "Pure-Python reimplementations"

    `deqms`, `proda`, `limma`, and `limrots` are standalone Python reimplementations
    of their R counterparts, built on numpy/scipy. They never call R or rpy2 — there
    is no runtime R dependency and no R-to-Python fallback; the Python path always
    runs.

## (c) Agentic optimization (no LLM required)

The `mokume.agentic` optimizer searches normalization, imputation, and DE choices
for each contrast and reports the best-scoring configuration. Set
`use_llm=False` to run a purely rule-based search with no API key.

=== "Python (package)"

    ```python
    import warnings
    from mokume.pipeline.features_to_proteins import QuantificationPipeline
    from mokume.pipeline.config import (
        PipelineConfig,
        InputConfig,
        QuantificationConfig,
    )
    from mokume.agentic.optimizer import optimize
    from mokume.agentic.config import AgenticConfig

    warnings.filterwarnings("ignore")

    config = PipelineConfig(
        input=InputConfig(
            parquet="python/tests/example/feature_wide.parquet",
            sdrf="python/tests/example/PXD020192.sdrf.tsv",
        ),
        quantification=QuantificationConfig(method="maxlfq"),
    )
    proteins = QuantificationPipeline(config).run()

    samples = [c for c in proteins.columns if c != proteins.columns[0]]
    half = len(samples) // 2
    sample_to_condition = {
        s: ("groupA" if i < half else "groupB")
        for i, s in enumerate(samples)
    }

    agentic = AgenticConfig(
        use_llm=False,                    # rule-based search, no API key
        contrasts=[("groupA", "groupB")],
        max_rounds=2,
        max_experiments=4,
        output_dir="./optimization",
    )
    states = optimize(proteins, sample_to_condition, agentic)

    # states is keyed by "<A>_vs_<B>"; each value is an AgenticState.
    state = states["groupA_vs_groupB"]
    print("best score:", state.best_score)
    print("experiments run:", state.total_experiments)
    print("converged:", state.converged)
    ```

If you already have a `QpxDataset` (see the [Full Pipeline](pipeline.md) page),
call `optimize_from_dataset` instead — it pulls the wide protein matrix off the
dataset and delegates to `optimize`:

=== "Python (package)"

    ```python
    from mokume.agentic.optimizer import optimize_from_dataset

    dataset = QuantificationPipeline(config).run_dataset()
    states = optimize_from_dataset(dataset, sample_to_condition, agentic)
    ```
