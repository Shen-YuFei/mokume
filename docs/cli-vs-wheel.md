# CLI vs Wheel

mokume gives you two ways to reach the same Rust compute kernel, plus a Python
periphery that only the wheel exposes. This page tells you which front door to
use and where each capability actually lives.

- **CLI binary** `mokume` — built with `cargo`, no Python runtime. Pure compute
  for pipelines and CI.
- **Wheel** `pip install mokume-rs` — the PyO3/maturin extension `mokume._mokume`
  runs the same kernel **in-process**, and additionally ships the Python
  periphery (plots, tissue maps, reports, and the pure-Python method fallbacks).

## Which one do I use?

| Your need | Use |
| --- | --- |
| Pure compute (quantify, normalize, impute, DE, ComBat) | **CLI binary** |
| Pipeline / batch / CI with no Python runtime available | **CLI binary** |
| Scripting the compute from Python | **Wheel** (`mokume.features2proteins(...)`, `mokume.run([...])`) |
| Plots (t-SNE, volcano, heatmap, PCA) | **Wheel** (`mokume.tsne_visualization`, `mokume.de_plots`) |
| Tissue proteome maps / atlas | **Wheel** (`mokume.tissuemap`) |
| Interactive HTML DE report | **Wheel** (`mokume.interactive_report`) |
| iBAQ QC report PDF | **Wheel** (`mokume.peptides2protein_qc`) |
| `missforest` imputation | **Wheel** (`mokume.impute(method=...)`) |
| iBAQ for an enzyme the Rust kernel does not digest | **Wheel** (`mokume.peptides2protein_ibaq`) |

The CLI binary is **pure compute**: it does not draw plots, build tissue maps,
or write HTML reports, and it does not shell out to Python to do so.

!!! warning "The CLI binary has exactly four subcommands"
    `features2proteins`, `features2peptides`, `peptides2protein`,
    `correct-batches`. There is **no** `tissuemap`, **no** `tsne-visualization`,
    and **no** `agentic` subcommand. `features2proteins` has **no** `--plot-*`,
    `--interactive-report`, or `--report-output` flags — those moved to the
    wheel. `agentic` remains in the separately installed pure-Python `mokume`
    package.

## Same compute, two front doors

The wheel's compute wrappers parse their arguments through the **same** `clap`
definition the binary uses, so a result computed through the wheel is the result
the binary produces.

=== "CLI binary"

    ```bash
    mokume features2proteins \
      --parquet features.parquet \
      --sdrf experiment.sdrf.tsv \
      --output proteins.csv \
      --quant-method maxlfq \
      --threads 24
    ```

=== "Wheel (kwargs)"

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

=== "Wheel (explicit argv)"

    ```python
    import mokume

    mokume.run([
        "features2proteins",
        "--parquet", "features.parquet",
        "--sdrf", "experiment.sdrf.tsv",
        "--output", "proteins.csv",
        "--quant-method", "maxlfq",
        "--threads", "24",
    ])
    ```

## The kwargs &rarr; flags rule

The compute wrappers (`features2proteins`, `features2peptides`,
`peptides2protein`, `correct_batches`) turn keyword arguments into CLI flags by a
fixed rule:

| Python keyword argument | CLI flag |
| --- | --- |
| `key="value"` | `--key value` (`_` in the key is rewritten to `-`) |
| `key=True` | `--key` |
| `key=[a, b]` | `--key a --key b` (the flag repeats once per item) |
| `key=None` or `key=False` | *skipped* |

So `quant_method="maxlfq"` becomes `--quant-method maxlfq`, `log2=True` becomes
`--log2`, and `de_contrasts=["A vs B", "C vs D"]` becomes
`--de-contrasts "A vs B" --de-contrasts "C vs D"`. When you need flags that
keyword arguments cannot express (e.g. the per-contrast `--contrast KEY A B CSV`
on the periphery plot commands), call `mokume.run([...])` or the explicit-argv
periphery wrappers with a literal argument list.

## Capability &rarr; location delegation map

| Capability | Lives in | Reach it via |
| --- | --- | --- |
| `features2proteins` (quant · norm · IRS · coverage · impute · ComBat · DE) | Rust kernel | CLI `features2proteins` / `mokume.features2proteins` |
| `features2peptides` (filters · factor normalization · peptide export) | Rust kernel | CLI `features2peptides` / `mokume.features2peptides` |
| `peptides2protein` (`ibaq`, `sum`, `top3`, `topn`, `maxlfq`, `directlfq`) | Rust kernel | CLI `peptides2protein` / `mokume.peptides2protein` |
| `correct-batches` (ComBat, native Rust) | Rust kernel | CLI `correct-batches` / `mokume.correct_batches` |
| t-SNE / DE / iBAQ-QC plots | Python periphery | `mokume.tsne_visualization` · `mokume.de_plots` · `mokume.peptides2protein_qc` |
| Tissue proteome maps | Python periphery | `mokume.tissuemap` |
| Interactive HTML DE report | Python periphery | `mokume.interactive_report` |
| QC report / workflow comparison | Python periphery | `mokume.qc_report` · `mokume.workflow_comparison` |
| `missforest` imputation | Python periphery | `mokume.impute(method=...)` |
| iBAQ for unported enzymes (e.g. `CNBr`) | Python periphery | `mokume.peptides2protein_ibaq` |
| `agentic` workflow optimizer | pure-Python `mokume` package | not in the Rust CLI or wheel |

!!! note "DirectLFQ and ComBat are native Rust"
    DirectLFQ and ComBat are part of the kernel and need no extra and no Python.
    The old `directlfq` and `batch-correction` Python extras are removed; ComBat
    is native Rust, oracle-verified vs `inmoose`.
