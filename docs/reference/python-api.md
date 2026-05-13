# Python API

## Quantification

### Factory Function

```python
from mokume.quantification import get_quantification_method, list_quantification_methods

# Create method by name
method = get_quantification_method("top3")        # TopNQuantification(n=3)
method = get_quantification_method("top10")       # TopNQuantification(n=10)
method = get_quantification_method("maxlfq", min_peptides=2, n_jobs=4)

# List available methods
available = list_quantification_methods()
# {'top3': True, 'topn': True, 'maxlfq': True, 'directlfq': False, 'sum': True}
```

### TopNQuantification

```python
from mokume.quantification import TopNQuantification

quant = TopNQuantification(n=3)  # or n=5, n=10, etc.
result = quant.quantify(
    peptides,
    protein_column="ProteinName",
    peptide_column="PeptideSequence",
    intensity_column="NormIntensity",
    sample_column="SampleID",
)
```

### MaxLFQQuantification

```python
from mokume.quantification import MaxLFQQuantification

quant = MaxLFQQuantification(
    min_peptides=2,      # Min peptides for MaxLFQ (uses median for fewer)
    threads=4,           # Parallel cores (-1 for all)
    force_builtin=False, # Force built-in implementation
)
result = quant.quantify(peptides, protein_column="ProteinName", ...)

# Check backend
quant.using_directlfq  # True/False
quant.name             # "MaxLFQ (DirectLFQ)" or "MaxLFQ (built-in)"
```

### DirectLFQQuantification

```python
from mokume.quantification import is_directlfq_available

if is_directlfq_available():
    from mokume.quantification import DirectLFQQuantification
    quant = DirectLFQQuantification(min_nonan=2)
    result = quant.quantify(peptides, protein_column="ProteinName", ...)
```

### AllPeptidesQuantification (Sum)

```python
from mokume.quantification import AllPeptidesQuantification

quant = AllPeptidesQuantification()
result = quant.quantify(peptides, protein_column="ProteinName", ...)
```

### peptides_to_protein (iBAQ)

```python
from mokume.quantification import peptides_to_protein

peptides_to_protein(
    fasta="proteome.fasta",
    peptides="peptides.csv",
    enzyme="Trypsin",
    normalize=True,
    tpa=True,
    ruler=True,
    ploidy=2,
    cpc=200,
    organism="human",
    output="proteins-ibaq.tsv",
    min_aa=7,
    max_aa=30,
    verbose=True,
    qc_report="QC.pdf",
)
```

---

## Pipeline

### PipelineConfig

```python
from mokume.pipeline import QuantificationPipeline, PipelineConfig
from mokume.pipeline.config import (
    InputConfig,
    FilterConfig,
    NormalizationConfig,
    QuantificationConfig,
    IRSConfig,
    BatchCorrectionConfig,
    DEConfig,
    OutputConfig,
)

config = PipelineConfig(
    input=InputConfig(
        parquet="data.parquet",
        sdrf="experiment.sdrf.tsv",
        fasta_file="proteome.fasta",
    ),
    filtering=FilterConfig(
        min_aa=7,
        min_unique_peptides=2,
        remove_contaminants=True,
    ),
    normalization=NormalizationConfig(
        run_method="median",
        sample_method="globalMedian",
    ),
    quantification=QuantificationConfig(
        method="maxlfq",
    ),
)

pipeline = QuantificationPipeline(config)
proteins = pipeline.run()
```

### features_to_proteins (functional)

```python
from mokume.pipeline import features_to_proteins

proteins = features_to_proteins(
    parquet="data.parquet",
    output="proteins.csv",
    sdrf="experiment.sdrf.tsv",
    quant_method="maxlfq",
    run_normalization="median",
    sample_normalization="globalMedian",
    batch_correction=True,
    batch_method="sample_prefix",
    batch_covariates=["characteristics[sex]"],
)
```

---

## Differential Expression

Supported methods: `limrots`, `deqms`, `proda`, `limma`, `rots`, `msstats`,
and `ensemble` (top-k consensus across multiple methods).
All methods are pure-Python reimplementations — no R or rpy2 required.

```python
from mokume.analysis import DifferentialExpression

sample_to_condition = {
    "S1": "HL",
    "S2": "HL",
    "S3": "NASH",
    "S4": "NASH",
}

de = DifferentialExpression(
    method="limrots",
    fdr_method="ihw",
    log2fc_threshold=0.5,
    fdr_threshold=0.05,
)

result = de.run(protein_df, sample_to_condition, ("NASH", "HL"))
```

```python
from mokume.analysis import DifferentialExpression

deqms = DifferentialExpression(
    method="deqms",
    peptide_counts=peptide_counts,
    fdr_method="bh",
)

results = deqms.run_comparisons(
    protein_df,
    sample_to_condition,
    [("NASH", "HL"), ("NASH", "Control")],
)
```

```python
# limma — stable empirical Bayes baseline
de = DifferentialExpression(method="limma", fdr_method="bh")
result = de.run(protein_df, sample_to_condition, ("NASH", "HL"))
```

```python
# Top-k consensus ensemble across multiple methods
from mokume.analysis.ensemble import run_ensemble

ensemble_df = run_ensemble(
    protein_df,
    sample_to_condition,
    contrast=("NASH", "HL"),
    methods=("limrots", "deqms", "proda"),
    min_k=2,
    fdr_method="bh",
    fdr_threshold=0.05,
    log2fc_threshold=0.5,
)
```

---

## Normalization

### Peptide Normalization Pipeline

```python
from mokume.normalization.peptide import peptide_normalization

peptide_normalization(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    min_aa=7,
    min_unique=2,
    remove_ids=None,
    remove_decoy_contaminants=True,
    remove_low_frequency_peptides=True,
    output="peptides.csv",
    skip_normalization=False,
    nmethod="median",          # Feature-level: median, mean, max, global, max_min, iqr, none
    pnmethod="globalMedian",   # Sample-level: globalMedian, conditionMedian, hierarchical, tmm, none
    log2=True,
    save_parquet=False,
)
```

The Python function keeps the historical parameter names `nmethod` and `pnmethod`, even though the CLI now uses `--run-normalization` and `--sample-normalization`.

### IRS Normalization

```python
from mokume.normalization.irs import (
    IRSNormalizer,
    detect_pooled_from_sdrf,
    detect_plexes_from_sdrf,
)

ref_samples = detect_pooled_from_sdrf("experiment.sdrf.tsv")
sample_to_plex = detect_plexes_from_sdrf("experiment.sdrf.tsv")

normalizer = IRSNormalizer(reference_samples=ref_samples, stat="median")
protein_df = normalizer.fit_transform(protein_df, sample_to_plex)
```

### LOESS Normalization

```python
from mokume.normalization import LOESSNormalizer, loess_normalize

# Apply to a sample x protein matrix on log2 scale
loess_df = loess_normalize(log2_df, frac=0.75, reference="median")

normalizer = LOESSNormalizer(frac=0.75, reference="median")
loess_df = normalizer.fit_transform(log2_df)
```

### MBQN Normalization

```python
from mokume.normalization import MBQNNormalizer, mbqn_normalize

mbqn_df = mbqn_normalize(log2_df)
mbqn_df = MBQNNormalizer().fit_transform(log2_df)
```

### VSN Normalization

```python
from mokume.normalization import vsn_normalize

# Input is linear-scale; output is glog2 (variance-stabilised)
stabilised_df = vsn_normalize(linear_df)
```

!!! note
    VSN is intentionally not exposed via `--sample-normalization` because
    its glog2 output is incompatible with the pipeline's downstream
    linear-scale assumptions. Use the standalone API only.

### Imputation Utilities

Low-level imputers (operate on a wide log2 matrix):

```python
from mokume.imputation import (
    impute_minprob,
    impute_mindet,
    impute_qrilc,
    impute_seqknn,
    impute_missforest,
    impute_mice,
    impute_mle,
    impute_bpca,
    impute_gms,
    impute_nbavg,
    impute_impseq,
    impute_impseqrob,
)

minprob_df = impute_minprob(log2_df, quantile=0.01, shift=1.6, scale=0.3)
mindet_df = impute_mindet(log2_df, quantile=0.01)
bpca_df = impute_bpca(log2_df, n_pcs=2)
```

High-level pipeline stage:

```python
from mokume.pipeline import ImputationConfig, ImputationStage, PipelineConfig

config = PipelineConfig(
    imputation=ImputationConfig(enabled=True, method="minprob"),
    # other stage configs ...
)
imputed_df = ImputationStage(config).impute(protein_df)
```

---

## TissueMap Pipeline

```python
from pathlib import Path

from mokume.tissuemap.config import InputConfig, OutputConfig, TissueMapConfig, load_config
from mokume.tissuemap.pipeline import TissueMapPipeline

config = TissueMapConfig(
    n_jobs=8,
    input=InputConfig(scan_dir=Path("QPX_data/tissues-mq/PXD016999")),
    output=OutputConfig(output_dir=Path("./tissuemap_results")),
)

TissueMapPipeline(config).run()
```

```python
from pathlib import Path

from mokume.tissuemap.config import load_config
from mokume.tissuemap.pipeline import TissueMapPipeline

config = load_config(
    Path("tissuemap.yaml"),
    overrides={
        "input.scan_dir": "QPX_data/tissues-mq",
        "output.output_dir": "./tissuemap_results",
        "n_jobs": 8,
    },
)

TissueMapPipeline(config).run()
```

---

## Preprocessing Filters

```python
from mokume.preprocessing.filters import (
    load_filter_config,
    save_filter_config,
    generate_example_config,
    get_filter_pipeline,
)
from mokume.model.filters import PreprocessingFilterConfig

# Generate example
generate_example_config("filters.yaml")

# Load from file
config = load_filter_config("filters.yaml")

# Create programmatically
config = PreprocessingFilterConfig(name="custom", enabled=True)
config.intensity.min_intensity = 1000.0
config.peptide.allowed_charge_states = [2, 3, 4]
config.protein.min_unique_peptides = 2
config.run_qc.max_missing_rate = 0.5

# Apply CLI-style overrides
config.apply_overrides({
    "min_intensity": 500,
    "charge_states": [2, 3],
    "max_missing_rate": 0.3,
})

# Save
save_filter_config(config, "my_filters.yaml")

# Apply
pipeline = get_filter_pipeline(config)
filtered_df, results = pipeline.apply(df)

for result in results:
    print(f"{result.filter_name}: removed {result.removed_count} ({result.removal_rate:.1%})")

summary = pipeline.summary(results)
print(f"Total removed: {summary['total_removed']} / {summary['total_input']}")
```

---

## Postprocessing

### Batch Correction

```python
from mokume.postprocessing import (
    apply_batch_correction,
    detect_batches,
    extract_covariates_from_sdrf,
)

batch_indices = detect_batches(
    sample_ids=df_wide.columns.tolist(),
    method="sample_prefix",  # "sample_prefix", "run", or "column"
)

covariates = extract_covariates_from_sdrf(
    "experiment.sdrf.tsv",
    sample_ids=df_wide.columns.tolist(),
    covariate_columns=["characteristics[sex]"],
)

df_corrected = apply_batch_correction(
    df=df_wide, batch=batch_indices, covs=covariates,
)
```

### Data Reshaping

```python
from mokume.postprocessing.reshape import (
    pivot_wider,
    pivot_longer,
    remove_samples_low_protein_number,
    remove_missing_values,
    describe_expression_metrics,
)

# Long to wide
df_wide = pivot_wider(df, row_name="ProteinName", col_name="SampleID", values="Ibaq")

# Wide to long
df_long = pivot_longer(df_wide, row_name="ProteinName", col_name="SampleID", values="Ibaq")

# Quality filtering
df = remove_samples_low_protein_number(df, min_protein_num=100)
df = remove_missing_values(df, missingness_percentage=20, expression_column="Ibaq")

# Statistics
metrics = describe_expression_metrics(df)
```

---

## IO

### AnnData

```python
from mokume.io.parquet import create_anndata

adata = create_anndata(
    df,
    obs_col="SampleID",
    var_col="ProteinName",
    value_col="Ibaq",
    layer_cols=["IbaqNorm", "IbaqLog"],
    obs_metadata_cols=["Condition"],
    var_metadata_cols=["GeneName"],
)
adata.write("proteins.h5ad")
```

### FASTA

```python
from mokume.io.fasta import (
    load_fasta,
    digest_protein,
    extract_fasta,
    get_protein_molecular_weights,
)

proteins = load_fasta("proteome.fasta")

peptides = digest_protein(
    sequence="MKWVTFISLLFLFSSAYS...",
    enzyme="Trypsin",
    min_aa=7, max_aa=30,
)

unique_peptide_counts, mw_dict, found_proteins = extract_fasta(
    fasta="proteome.fasta",
    enzyme="Trypsin",
    proteins=["P12345", "P67890"],
    min_aa=7, max_aa=30,
    tpa=True,
)

mw_dict = get_protein_molecular_weights("proteome.fasta", ["P12345", "P67890"])
```

---

## Organism Metadata

```python
from mokume.model.organism import OrganismDescription

organisms = OrganismDescription.registered_organisms()
# ['human', 'mouse', 'yeast', 'drome', 'caeel', 'schpo']

human = OrganismDescription.get("human")
print(human.genome_size)
print(human.histone_entries)
```
