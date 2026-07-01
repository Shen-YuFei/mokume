"""Data profiling for agentic analysis."""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.profiler")


@dataclass
class DataProfile:
    """Characterization of a protein intensity matrix."""

    n_proteins: int
    n_samples: int
    n_conditions: int
    samples_per_condition: dict[str, int]
    missing_rate: float
    missing_per_sample: dict[str, float]
    median_cv: float | None
    has_peptide_counts: bool
    data_type: str
    batch_fields: list[str]
    intensity_range: tuple[float, float]
    is_log_transformed: bool
    # Enriched signals (P1-8). Default to None / empty so older fixtures keep
    # working without explicitly populating them.
    pca_explained_variance: list[float] = field(default_factory=list)
    per_condition_median_cv: dict[str, float] = field(default_factory=dict)
    pairwise_median_abs_log2fc: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Serialize to JSON-friendly dictionary."""
        return {
            "n_proteins": self.n_proteins,
            "n_samples": self.n_samples,
            "n_conditions": self.n_conditions,
            "samples_per_condition": self.samples_per_condition,
            "missing_rate": round(self.missing_rate, 4),
            "median_cv": round(self.median_cv, 4) if self.median_cv else None,
            "has_peptide_counts": self.has_peptide_counts,
            "data_type": self.data_type,
            "batch_fields": self.batch_fields,
            "intensity_range": [round(v, 2) for v in self.intensity_range],
            "is_log_transformed": self.is_log_transformed,
            "pca_explained_variance": [
                round(v, 4) for v in self.pca_explained_variance
            ],
            "per_condition_median_cv": {
                k: round(v, 4) for k, v in self.per_condition_median_cv.items()
            },
            "pairwise_median_abs_log2fc": {
                k: round(v, 4) for k, v in self.pairwise_median_abs_log2fc.items()
            },
        }


def _detect_log_transformed(values: np.ndarray) -> bool:
    """Heuristic: if max < 40 and median < 30, likely log2-transformed."""
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return False
    return float(np.nanmax(finite)) < 40 and float(np.nanmedian(finite)) < 30


def _compute_cv(
    matrix: pd.DataFrame,
    sample_to_condition: dict[str, str],
) -> float | None:
    """Compute median within-condition CV across proteins."""
    conditions = set(sample_to_condition.values())
    all_cvs = []
    for cond in conditions:
        cols = [s for s in matrix.columns if sample_to_condition.get(s) == cond]
        if len(cols) < 2:
            continue
        subset = matrix[cols]
        row_mean = subset.mean(axis=1)
        row_std = subset.std(axis=1, ddof=1)
        cv = row_std / row_mean.replace(0, np.nan)
        all_cvs.extend(cv.dropna().tolist())
    return float(np.median(all_cvs)) if all_cvs else None


def _compute_per_condition_cv(
    matrix: pd.DataFrame,
    sample_to_condition: dict[str, str],
) -> dict[str, float]:
    """Return per-condition median CV (proteins with >=2 samples in condition)."""
    conditions = sorted(set(sample_to_condition.values()))
    out: dict[str, float] = {}
    for cond in conditions:
        cols = [s for s in matrix.columns if sample_to_condition.get(s) == cond]
        if len(cols) < 2:
            continue
        subset = matrix[cols]
        row_mean = subset.mean(axis=1)
        row_std = subset.std(axis=1, ddof=1)
        cv = (row_std / row_mean.replace(0, np.nan)).dropna()
        if len(cv) > 0:
            out[cond] = float(cv.median())
    return out


def _compute_pca_evr(matrix: pd.DataFrame, n_components: int = 2) -> list[float]:
    """Return explained-variance ratios for the first N principal components.

    Drops proteins with missing values; falls back to an empty list when too
    few samples / proteins remain for a stable decomposition.
    """
    if matrix.shape[1] < 3 or matrix.shape[0] < 5:
        return []
    complete = matrix.dropna(how="any")
    if complete.shape[0] < 5:
        return []
    data = complete.values.T  # samples x proteins
    n_comp = min(n_components, data.shape[0] - 1, data.shape[1])
    if n_comp < 1:
        return []
    try:
        from sklearn.decomposition import PCA
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for agentic PCA profiling. "
            "Install it with: pip install mokume[agentic]"
        ) from exc
    pca = PCA(n_components=n_comp)
    try:
        pca.fit(data)
    except (ValueError, np.linalg.LinAlgError):
        return []
    return [float(v) for v in pca.explained_variance_ratio_]


def _compute_pairwise_log2fc(
    matrix: pd.DataFrame,
    sample_to_condition: dict[str, str],
    is_log: bool,
) -> dict[str, float]:
    """Median |log2FC| of condition-mean profiles for each unordered pair."""
    conditions = sorted(set(sample_to_condition.values()))
    if len(conditions) < 2:
        return {}
    cond_means: dict[str, pd.Series] = {}
    for cond in conditions:
        cols = [s for s in matrix.columns if sample_to_condition.get(s) == cond]
        if not cols:
            continue
        cond_means[cond] = matrix[cols].mean(axis=1, skipna=True)
    out: dict[str, float] = {}
    keys = sorted(cond_means)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            mean_a, mean_b = cond_means[a], cond_means[b]
            if is_log:
                diff = (mean_a - mean_b).dropna()
            else:
                ratio = (mean_a / mean_b.replace(0, np.nan)).replace(
                    [np.inf, -np.inf], np.nan
                )
                ratio = ratio.where(ratio > 0).dropna()
                diff = np.log2(ratio) if len(ratio) > 0 else pd.Series(dtype=float)
            if len(diff) > 0:
                out[f"{a}_vs_{b}"] = float(np.median(np.abs(diff)))
    return out


def _detect_data_type(
    sample_to_condition: dict[str, str],
    n_samples: int,
) -> str:
    """Heuristic data type detection from sample naming patterns."""
    sample_names = list(sample_to_condition.keys())
    joined = " ".join(sample_names).lower()
    if "tmt" in joined or "plex" in joined:
        return "TMT"
    if "dia" in joined or "diann" in joined:
        return "DIA"
    if n_samples <= 100:
        return "LFQ"
    return "unknown"


def _detect_batch_fields(sample_ids: list[str]) -> list[str]:
    """Detect potential batch structure from sample ID prefixes."""
    prefixes = {s.split("-")[0] for s in sample_ids if "-" in s}
    if len(prefixes) > 1:
        return ["sample_prefix"]
    return []


def profile_data(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    peptide_counts: pd.Series | None = None,
) -> DataProfile:
    """Profile a protein intensity matrix for agentic optimization."""
    protein_col = protein_df.columns[0]
    sample_cols = [c for c in protein_df.columns if c != protein_col]
    matrix = protein_df.set_index(protein_col)[sample_cols]

    values = matrix.values.astype(float)
    finite_vals = values[np.isfinite(values) & ~np.isnan(values)]

    # Condition stats
    cond_counts: dict[str, int] = {}
    for s in sample_cols:
        cond = sample_to_condition.get(s, "unknown")
        cond_counts[cond] = cond_counts.get(cond, 0) + 1

    # Missing rate
    total_cells = matrix.size
    missing_count = int(matrix.isna().sum().sum())
    missing_rate = missing_count / total_cells if total_cells > 0 else 0.0

    missing_per_sample = {col: float(matrix[col].isna().mean()) for col in sample_cols}

    is_log = _detect_log_transformed(values)

    return DataProfile(
        n_proteins=len(matrix),
        n_samples=len(sample_cols),
        n_conditions=len(cond_counts),
        samples_per_condition=cond_counts,
        missing_rate=missing_rate,
        missing_per_sample=missing_per_sample,
        median_cv=_compute_cv(matrix, sample_to_condition),
        has_peptide_counts=peptide_counts is not None,
        data_type=_detect_data_type(sample_to_condition, len(sample_cols)),
        batch_fields=_detect_batch_fields(sample_cols),
        intensity_range=(
            float(np.nanmin(finite_vals)) if len(finite_vals) > 0 else 0.0,
            float(np.nanmax(finite_vals)) if len(finite_vals) > 0 else 0.0,
        ),
        is_log_transformed=is_log,
        pca_explained_variance=_compute_pca_evr(matrix),
        per_condition_median_cv=_compute_per_condition_cv(matrix, sample_to_condition),
        pairwise_median_abs_log2fc=_compute_pairwise_log2fc(
            matrix, sample_to_condition, is_log
        ),
    )
