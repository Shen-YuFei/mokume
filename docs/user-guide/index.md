# User Guide

This section covers the Rust computation commands. The kernel ships both as a standalone CLI binary (`mokume`, built with cargo and no Python) and as a PyO3/maturin wheel (`pip install mokume-rs`) that imports the compiled `mokume._mokume` extension and runs the same commands in-process — no subprocess delegation. These two Rust entry points share one implementation; the separately maintained pure-Python computation package is documented under [Python API (package)](../reference/python-api-package.md).

## Compute commands

The CLI binary exposes exactly four compute subcommands. Each can also be driven from the wheel through a thin keyword wrapper (`mokume.<command>(**kwargs)`) or with an explicit argument list (`mokume.run([...])`).

### [features2proteins: Unified Pipeline](features2proteins.md)

The recommended entry point. Takes raw feature data and produces a protein quantification matrix in one step, with optional normalization, batch correction, IRS, imputation, and differential expression.

### [features2peptides: Peptide Normalization](features2peptides.md)

Normalizes and filters feature-level data into peptide intensities. Use this when you need fine-grained control over the normalization step before quantification.

### [peptides2protein: Protein Quantification](peptides2protein.md)

Quantifies proteins from normalized peptide data. Supports iBAQ (with TPA and ProteomicRuler), TopN, MaxLFQ, DirectLFQ, and Sum.

### [correct-batches: Batch Correction](batch-correct.md)

Standalone batch correction for already-quantified protein data. Combines multiple files and applies native Rust ComBat (oracle-verified vs inmoose) correction.

## Periphery (wheel-only)

Plotting, tissue maps, and interactive reports are not CLI subcommands; they live in the Python wheel under `mokume.commands` and are reached through periphery functions such as `mokume.tsne_visualization`, `mokume.tissuemap`, `mokume.de_plots`, and `mokume.interactive_report`. They read the tables the kernel produced. Install the relevant extra (`plotting`, `reports`, `tissuemap`, `ibaq`, `analysis`, or `all`) to pull in the periphery libraries. The agentic workflow search lives in the separate `mokume_py` package and is not part of this toolkit.
