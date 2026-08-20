# mokume Roadmap

Future improvements and features under consideration, based on recent benchmarking studies.

---

## Batch Correction Enhancements

Based on comprehensive benchmarking comparing 7 algorithms × 3 quantification methods × 3 data levels across balanced and confounded scenarios.

### Key Findings

- **Protein-level correction is most robust** - validates current mokume architecture
- **Quantification method interacts with batch correction** - not all combinations work equally well
- **MaxLFQ + Ratio** showed superior performance in clinical validation (1431 plasma samples)
- **Confounded scenarios** (batch mixed with biology) require different approaches than balanced designs

### High Priority

- [ ] **Ratio-based batch correction**: Scale features relative to reference material or batch median
  - Best performer for confounded scenarios (common in real multi-center studies)
  - Simple implementation: `corrected = sample / reference_median`
  - Validated on Quartet reference materials and clinical cohort

- [ ] **Median centering batch correction**: Per-batch median subtraction
  - Tested in benchmark, simple and effective
  - Different from normalization - applied post-quantification

### Medium Priority

- [ ] **RUV-III-C**: Remove Unwanted Variation using negative control proteins
  - Requires defining control proteins not affected by biology
  - Good for studies with spike-in controls

- [ ] **WaveICA2.0**: Wavelet-based Independent Component Analysis
  - Decomposes batch effects using wavelet transform
  - Good benchmark performance

- [ ] **Balanced vs Confounded detection**: Auto-detect scenario type to recommend appropriate method
  - Balanced: biology and batch are independent → most methods work
  - Confounded: biology correlates with batch → only Ratio method effective

### Low Priority

- [ ] Harmony - operates on embeddings, doesn't output expression matrix
- [ ] NormAE - autoencoder approach, marginal gains over simpler methods

### Validated (Already Implemented)

- [x] **Protein-level batch correction** - benchmark confirms most robust strategy
- [x] **ComBat** - native Rust implementation (oracle-verified against inmoose), with covariate support
- [x] **MaxLFQ, Top3/TopN, piBAQ quantification** - all three tested, MaxLFQ recommended

### Evaluation Metrics to Consider

From Quartet benchmark framework:

- [ ] **SNR (Signal-to-Noise Ratio)**: Biological difference / technical variation
- [ ] **Intraclass Correlation Coefficient (ICC)**: Reproducibility across batches
- [ ] **Silhouette score**: Cluster separation for known groups

### References

1. Chen Q, et al. (2025) "Protein-level batch-effect correction enhances robustness in MS-based proteomics" *Nature Communications* [doi:10.1038/s41467-025-64718-y](https://www.nature.com/articles/s41467-025-64718-y)

2. "Correcting batch effects in large-scale multiomics studies using a reference-material-based ratio method" (2023) *Genome Biology* [doi:10.1186/s13059-023-03047-z](https://link.springer.com/article/10.1186/s13059-023-03047-z)

3. "Quartet protein reference materials and datasets for multi-platform assessment of label-free proteomics" (2023) *Genome Biology* [doi:10.1186/s13059-023-03048-y](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03048-y)

---

## Single-Cell Proteomics (SCP) Support

Based on benchmarking of DIA-based single-cell proteomics workflows covering software tools, searching strategies, and postprocessing methods.

### Key Findings

- **Normalization has the largest impact** on downstream analysis
- **Imputation has minor/negative impact** on high-quality DIA data
- **Missing values are prevalent** - proteins near detection limit
- **Batch effects can be mistaken for cell heterogeneity** - critical to address

### High Priority

- [ ] **SCoPE2-style normalization**:
  - Divide by reference channel (5-cell equivalents)
  - Median center columns (per-cell)
  - Mean center rows (per-protein)
  - Standard for TMT-based SCP

- [ ] **Sparsity filtering**:
  - `min_cells_per_protein`: Require protein detected in X% of cells (e.g., 30%)
  - `min_proteins_per_cell`: QC threshold for cell quality (e.g., 500 proteins)
  - Critical for removing low-quality cells and unreliable proteins

- [ ] **SCP QC metrics**:
  - Missingness per cell
  - Proteins quantified per cell
  - CV distributions across cells
  - Batch-specific statistics

### Medium Priority

- [ ] **HarmonizR integration**: Matrix dissection approach for ComBat/limma with missing values
  - Handles features not present in all batches (common in SCP)
  - No imputation required

- [ ] **Carrier channel support**: For TMT-based SCP with carrier (200-cell) and reference (5-cell) channels
  - Carrier boosts identification
  - Reference used for normalization

- [ ] **DIA-NN output optimization**: Direct integration with DIA-NN library-free search
  - Best software for DIA-SCP per benchmarks

### Low Priority (Questionable Benefit for SCP)

- [ ] Complex imputation (KNN, deep learning)
  - Evidence shows minimal benefit for high-quality DIA data
  - Risk of introducing artifacts in sparse single-cell data
  - Better to report missingness honestly

- [ ] scRNA-seq methods (scran, SCnorm)
  - Designed for count data with different statistical properties
  - Proteomics intensities require different approaches

### Data Quality Considerations

From SCP benchmarks:

- Average missingness ~25% for good data, up to 70% with many batches
- Sparsity reduction before analysis is critical
- Cell-level QC more important than aggressive imputation

### References

1. Wang J, et al. (2025) "Benchmarking informatics workflows for data-independent acquisition single-cell proteomics" *Nature Communications* [doi:10.1038/s41467-025-65174-4](https://www.nature.com/articles/s41467-025-65174-4)

2. "Benchmark of Data Integration in Single-Cell Proteomics" (2025) *Analytical Chemistry* [doi:10.1021/acs.analchem.4c04933](https://pubs.acs.org/doi/10.1021/acs.analchem.4c04933)

3. Vanderaa C & Gatto L. "scp: Mass Spectrometry-Based Single-Cell Proteomics Data Analysis" *Bioconductor* [link](https://bioconductor.org/packages/scp)

4. Specht H, et al. (2021) "Single-cell proteomic and transcriptomic analysis of macrophage heterogeneity using SCoPE2" *Genome Biology* [doi:10.1186/s13059-021-02267-5](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02267-5)

---

## Other Future Improvements

### Normalization

- [x] Quantile normalization — in the Rust kernel (`features2proteins --sample-norm quantile`)
- [x] LOESS normalization for intensity-dependent bias — in the Rust kernel (`features2proteins --sample-norm loess`)

### Quantification

- [ ] Phosphoproteomics-specific quantification
- [ ] PTM site-level rollup options

### Integration

- [x] Export to AnnData format — via tissuemap pipeline (`.h5ad`) and `correct-batches --export_anndata`
- [x] Differential expression (limma, DEqMS, proDA, LimROTS, ROTS) — in the Rust kernel (`features2proteins --de`)
- [ ] Nextflow/Snakemake workflow templates
