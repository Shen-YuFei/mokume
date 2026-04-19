"""Data profiling for agentic analysis."""

from dataclasses import dataclass

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
        is_log_transformed=_detect_log_transformed(values),
    )
