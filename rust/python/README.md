# mokume-rs — Python wheel

`pip install mokume-rs` gives you a Rust compute kernel with a Python periphery,
in the PyO3/maturin layout used by projects such as polars and pydantic-core
(Python imports a compiled Rust extension).

- The compute-heavy pipeline (expression matrix, normalization, imputation,
  differential expression, batch correction) lives in the Rust crates under
  `crates/` and is exposed in-process through the compiled `mokume._mokume`
  extension that maturin builds from `crates/mokume-py`.
- The periphery lives here in `mokume/commands/`. Plotting and reporting consume
  kernel tables, TissueMap derives downstream atlas outputs from QPX data, and
  explicit Python-only fallbacks handle operations the kernel does not provide.

The standalone Rust CLI binary (`crates/mokume-cli`, built with `cargo`) is the
no-Python compute path for pipelines; it does **not** shell out to Python. The
periphery is reached only through this wheel.

## Install

```bash
pip install mokume-rs                 # compute kernel + Python API
pip install "mokume-rs[plotting]"     # + t-SNE / DE plots / iBAQ QC report
pip install "mokume-rs[tissuemap]"    # + per-dataset tissue proteome analysis
pip install "mokume-rs[reports]"      # + interactive HTML DE report
pip install "mokume-rs[ibaq]"         # + pure-Python iBAQ for unported enzymes
pip install "mokume-rs[analysis]"     # + QC/comparison reports + missforest
pip install "mokume-rs[all]"          # everything
```

From a source checkout, build the extension with `maturin develop` (or
`maturin build`); the periphery is plain Python and needs no build step.

## Compute API

The compute commands run in-process through the same clap parsing the CLI binary
uses, so flags stay single-sourced in Rust:

```python
import mokume

mokume.features2proteins(parquet="features.parquet", output="proteins.csv")
mokume.peptides2protein(method="ibaq", peptides="peptides.parquet",
                        fasta="proteome.fasta", output="proteins.tsv")
# full control over the argument vector:
mokume.run(["correct-batches", "--input", "ibaq.tsv", "--output", "corrected.tsv"])
```

Each wrapper maps keyword arguments to CLI flags (`key=value` → `--key value`
with `_` rewritten to `-`; `key=True` → `--key`; a list repeats the flag; `None`
/ `False` are skipped).

## Periphery commands

Each lives in `mokume.commands.<name>` with a `main(argv)` entry point, is
runnable as `python -m mokume.commands.<name>`, and most have an ergonomic
wrapper on the top-level package:

```python
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")
mokume.tissuemap(scan_dir="./data", output_dir="./out")
mokume.peptides2protein_qc(protein_table="proteins.tsv", qc_report="QC.pdf")
mokume.peptides2protein_ibaq(peptides="peptides.parquet", fasta="proteome.fasta",
                             enzyme="CNBr", output="proteins.tsv")
# de_plots / interactive_report take an explicit argv (the per-contrast
# --contrast KEY A B CSV flag repeats, which keyword arguments cannot express):
mokume.de_plots(["--protein-matrix", "proteins.csv", "--plot-dir", "plots",
                 "--volcano", "--contrast", "c1", "A", "B", "de.csv"])

# QC / comparison reports + the pure-Python method fallbacks (analysis extra):
mokume.qc_report(protein_matrix="proteins.csv", sdrf="x.sdrf.tsv", output="qc.html")
mokume.impute("proteins.csv", method="missforest", output="imputed.csv")
```

The plotting and reporting functions read the tables the kernel wrote without
recomputing them, so their cells match the kernel output. The explicit
`peptides2protein_ibaq` and `impute(method="missforest")` fallbacks instead
compute operations the Rust kernel does not provide.

| Command                       | Extra        | Third-party libraries                                  |
| ----------------------------- | ------------ | ------------------------------------------------------ |
| `mokume.tsne_visualization`   | `plotting`   | numpy, pandas, scipy, scikit-learn, matplotlib, seaborn |
| `mokume.peptides2protein_qc`  | `plotting`   | numpy, pandas, matplotlib, seaborn                     |
| `mokume.de_plots`             | `plotting`   | numpy, pandas, matplotlib, seaborn, scikit-learn       |
| `mokume.interactive_report`   | `reports`    | numpy, pandas, plotly                                  |
| `mokume.tissuemap`            | `tissuemap`  | scanpy, anndata, umap-learn, combat, matplotlib, seaborn, pyarrow |
| `mokume.peptides2protein_ibaq`| `ibaq`       | pyopenms, pyarrow, PyYAML, numpy, pandas, scipy        |
| `mokume.qc_report` / `mokume.workflow_comparison` | `analysis` | numpy, pandas, scipy, scikit-learn         |
| `mokume.impute`                                   | `analysis` | numpy, pandas, scipy, scikit-learn         |

The exact dependency lists are declared in `pyproject.toml`'s
`[project.optional-dependencies]`. `mokume.impute(method="missforest")` reaches the
pure-Python imputer the Rust kernel does not reproduce (it wraps a scikit-learn
estimator); the kernel's `features2proteins` errors point here.
