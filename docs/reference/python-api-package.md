# Python API (package)

The pure-Python `mokume` package (`pip install mokume`) exposes a class-based API
for quantification and differential expression, built on an OOP data layer
(`QpxDataset`) with a pluggable method registry and a selectable compute backend.
The reference below is rendered directly from the package source, so it always
matches the installed version rather than a hand-written copy.

For the Rust-accelerated wheel's thin in-process API, see
[Python API (wheel)](python-api.md).

## Quantification

::: mokume.quantification.get_quantification_method

::: mokume.quantification.topn.TopNQuantification

::: mokume.quantification.maxlfq.MaxLFQQuantification

::: mokume.quantification.all_peptides.AllPeptidesQuantification

## Differential expression

::: mokume.analysis.DifferentialExpression

## Data layer and pipeline

A `QpxDataset` is the container the pipeline produces. It holds the data levels
(`psms`, `features`, `peptides`, `proteins`), sample metadata, and a provenance
log of the steps applied. `run_pipeline` resolves the configured quantification
method from the `PluginRegistry`, dispatches to the flow matching the method's
`input_level`, and returns a populated `QpxDataset`.

```python
from mokume.pipeline.config import PipelineConfig, InputConfig, QuantificationConfig
from mokume.pipeline.runner import run_pipeline

config = PipelineConfig(
    input=InputConfig(parquet="features.parquet"),
    quantification=QuantificationConfig(method="maxlfq"),
)
dataset = run_pipeline(config)                 # QpxDataset with .proteins populated
protein_matrix = dataset.get_level("proteins")  # protein x sample DataFrame
```

`QuantificationPipeline(config).run_dataset()` returns the same `QpxDataset` from
the pipeline object directly; the legacy `run()` returns the bare protein matrix.

### QpxDataset

::: mokume.core.dataset.QpxDataset

### PluginRegistry

The registry maps a name to a method class within one of the plugin groups
(`quantification`, `normalization.feature`, `normalization.sample`, `imputation`,
`harmonization`, `filter`). Register a class with the `@PluginRegistry.register`
decorator; resolve one with `get`; list a group with `available`.

```python
from mokume.core.registry import PluginRegistry

PluginRegistry.available("quantification")
# ['directlfq', 'ibaq', 'maxlfq', 'median', 'ratio', 'spectral_count', 'sum', ...]
method = PluginRegistry.get("quantification", "maxlfq")
```

::: mokume.core.registry.PluginRegistry

### run_pipeline

::: mokume.pipeline.runner.run_pipeline

## Compute backend

`RuntimeConfig.backend` selects the compute engine. The default `"python"` runs
the pure-Python pipeline; `"rust"` routes supported loading, filtering,
normalization, and quantification settings through the compiled
`mokume._mokume` kernel (installed with the `mokume-rs` wheel), then returns to
Python for postprocessing. When the kernel is absent, the `"rust"` backend
raises a clear error rather than silently falling back. The hybrid profile also
rejects `duckdb_memory`, `duckdb_threads`, and ion alignment other than `None`
or `"none"` before invoking the extension because it cannot honor those settings
with their documented semantics.

```python
from mokume.pipeline.config import RuntimeConfig

config = PipelineConfig(
    input=InputConfig(parquet="features.parquet"),
    quantification=QuantificationConfig(method="sum"),
    runtime=RuntimeConfig(backend="rust"),  # requires the mokume-rs wheel
)
dataset = run_pipeline(config)
```

::: mokume.pipeline.config.RuntimeConfig

## Agentic optimization

`optimize_from_dataset` runs the LLM-driven differential-expression optimizer on
the protein matrix carried by a `QpxDataset`, mirroring the DataFrame-based
`optimize` entry point.

```python
from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import optimize_from_dataset

states = optimize_from_dataset(
    dataset,                                    # a QpxDataset with .proteins populated
    sample_to_condition={"S1": "A", "S2": "B"},
    config=AgenticConfig(use_llm=False),
)
```

::: mokume.agentic.optimizer.optimize_from_dataset
