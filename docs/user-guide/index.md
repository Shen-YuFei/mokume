# User Guide

This section covers the practical usage of each mokume command, with complete examples for both CLI and Python API.

## Commands

### [features2proteins: Unified Pipeline](features2proteins.md)

The recommended entry point. Takes raw feature data and produces a protein quantification matrix in one step, with optional normalization, batch correction, IRS, differential expression, and visualization.

### [features2peptides: Peptide Normalization](features2peptides.md)

Normalizes and filters feature-level data into peptide intensities. Use this when you need fine-grained control over the normalization step before quantification.

### [peptides2protein: Protein Quantification](peptides2protein.md)

Quantifies proteins from normalized peptide data. Supports iBAQ (with TPA and ProteomicRuler), TopN, MaxLFQ, DirectLFQ, and Sum.

### [correct-batches: Batch Correction](batch-correct.md)

Standalone batch correction for already-quantified protein data. Combines multiple files and applies ComBat correction.

### [tissuemap: Tissue Proteome Atlas](tissuemap.md)

Builds a per-dataset tissue proteome atlas from QPX outputs, including AdaTiSS tissue-specificity scoring, AnnData exports, and atlas-style plots.

### [agentic optimize: LLM-Assisted Workflow Search](agentic-optimization.md)

Searches normalization × imputation × DE × ensemble configurations with an LLM-driven (or rule-based) loop, scoring each against ground truth or unsupervised QC metrics until convergence.

### [Visualization & Reports](visualization.md)

PCA, t-SNE, volcano plots, heatmaps, and interactive HTML QC reports.
