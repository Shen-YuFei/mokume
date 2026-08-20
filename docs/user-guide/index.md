# User Guide

This section covers the practical usage of each mokume command. The
`mokume` wheel installs a `mokume` console command backed by the compiled
`mokume._mokume` extension. It runs the leading Rust kernel in-process. The
separately maintained pure-Python distribution is documented under
[Python API (package)](../reference/python-api-package.md).

## Compute commands

The console command exposes exactly four compute subcommands. Each can also be
driven through a thin keyword wrapper (`mokume.<command>(**kwargs)`) or with an
explicit argument list (`mokume.run([...])`).

### [features2proteins: Unified Pipeline](features2proteins.md)

The recommended entry point. Takes raw feature data and produces a protein quantification matrix in one step, with optional normalization, batch correction, IRS, imputation, and differential expression.

### [features2peptides: Peptide Normalization](features2peptides.md)

Normalizes and filters feature-level data into peptide intensities. Use this when you need fine-grained control over the normalization step before quantification.

### [peptides2protein: Protein Quantification](peptides2protein.md)

Quantifies proteins from normalized peptide data. Supports piBAQ (with TPA and ProteomicRuler), TopN, MaxLFQ, DirectLFQ, and Sum.

### [correct-batches: Batch Correction](batch-correct.md)

Standalone batch correction for already-quantified protein data. Combines multiple files and applies native Rust ComBat (oracle-verified vs inmoose) correction.

## Periphery (wheel-only)

Plotting, tissue maps, and interactive reports are not CLI subcommands; they
live in the Python wheel under `mokume.commands` and are reached through
periphery functions such as `mokume.tsne_visualization`, `mokume.tissuemap`,
`mokume.de_plots`, and `mokume.interactive_report`. Plotting and reporting
consume kernel tables, while TissueMap derives its downstream atlas from QPX
data. Install the relevant extra (`plotting`, `reports`, `tissuemap`, `pibaq`,
`analysis`, or `all`) to pull in the periphery libraries. The agentic workflow
is provided by the [Mokume Plugin](agentic-plugin.md), whose local MCP service
uses the default Rust-backed wheel.
