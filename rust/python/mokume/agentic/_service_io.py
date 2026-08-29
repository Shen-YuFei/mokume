"""Input validation and atomic artifact I/O for the agentic service."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
import csv
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .context import BoundContext, GeneratedRecommendationBlock
    from .knowledge import KnowledgeGraph


def validated_object(
    payload: dict[str, Any] | None,
    allowed_fields: frozenset[str],
    name: str,
) -> dict[str, Any]:
    """Require an object with no fields outside its public contract."""
    if payload is None:
        return {}
    if not isinstance(payload, dict):
        raise ValueError(f"{name} must be an object or null")
    unknown = set(payload) - allowed_fields
    if unknown:
        raise ValueError(f"Unknown {name} fields: {sorted(unknown)}")
    return payload


def parse_contrast(contrast: list[str]) -> tuple[str, str]:
    """Require exactly two distinct, non-empty condition labels."""
    if not isinstance(contrast, list) or len(contrast) != 2:
        raise ValueError("contrast must be a two-item list")
    if not all(isinstance(item, str) and item.strip() for item in contrast):
        raise ValueError("contrast entries must be non-empty strings")
    if contrast[0] == contrast[1]:
        raise ValueError("contrast conditions must be distinct")
    return contrast[0], contrast[1]


def read_matrix(path: str) -> pd.DataFrame:
    """Read a comma- or tab-delimited protein matrix with numeric sample cells."""
    matrix_path = Path(path).expanduser().resolve()
    with matrix_path.open(encoding="utf-8-sig") as handle:
        header = handle.readline()
    separator = "\t" if header.count("\t") > header.count(",") else ","
    columns = next(csv.reader([header], delimiter=separator))
    if any(not column.strip() for column in columns):
        raise ValueError("Protein matrix column names must be non-empty")
    if len(columns) != len(set(columns)):
        raise ValueError("Protein matrix column names must be unique")
    frame = pd.read_csv(matrix_path, sep=separator)
    if frame.shape[1] < 3:
        raise ValueError(
            "Protein matrix requires an identifier and at least two samples"
        )
    identifiers = frame.iloc[:, 0]
    if identifiers.isna().any() or identifiers.astype(str).str.strip().eq("").any():
        raise ValueError("Protein identifiers must be non-empty")
    identifiers = identifiers.astype(str)
    if identifiers.duplicated().any():
        raise ValueError("Protein identifiers must be unique")
    frame[frame.columns[0]] = identifiers
    numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="raise")
    values = numeric.to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError("Protein matrix contains an infinite intensity")
    if not np.isfinite(values).any():
        raise ValueError("Protein matrix contains no finite intensities")
    frame[numeric.columns] = numeric
    return frame


def canonicalize_missing(frame: pd.DataFrame, input_scale: str) -> pd.DataFrame:
    """Use the Rust kernel's missing-value semantics at the service boundary."""
    if input_scale == "log2":
        return frame
    canonical = frame.copy()
    sample_columns = canonical.columns[1:]
    canonical[sample_columns] = canonical[sample_columns].mask(
        canonical[sample_columns] <= 0.0
    )
    if not np.isfinite(canonical[sample_columns].to_numpy(dtype=float)).any():
        raise ValueError("Linear protein matrix contains no positive intensities")
    return canonical


def load_peptide_counts(
    path: str | None,
    protein_df: pd.DataFrame,
) -> pd.Series | None:
    """Read positive integer peptide counts keyed by protein identifier."""
    if path is None:
        return None
    if not isinstance(path, str):
        raise ValueError("peptide_counts must be a string or null")
    require_file(path, "peptide_counts")
    count_path = Path(path).expanduser().resolve()
    with count_path.open(encoding="utf-8-sig") as handle:
        header = handle.readline()
    separator = "\t" if header.count("\t") > header.count(",") else ","
    columns = next(csv.reader([header], delimiter=separator))
    if columns != ["protein", "peptide_count"]:
        raise ValueError(
            "Peptide-count sidecar columns must be protein and peptide_count"
        )
    frame = pd.read_csv(count_path, sep=separator)
    if frame.empty:
        raise ValueError("Peptide-count sidecar contains no proteins")
    proteins = frame["protein"]
    if proteins.isna().any() or proteins.astype(str).str.strip().eq("").any():
        raise ValueError("Peptide-count protein identifiers must be non-empty")
    proteins = proteins.astype(str)
    if proteins.duplicated().any():
        raise ValueError("Peptide-count protein identifiers must be unique")
    counts = pd.to_numeric(frame["peptide_count"], errors="raise")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 1).any():
        raise ValueError("Peptide counts must be finite positive integers")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Peptide counts must be finite positive integers")
    series = pd.Series(values.astype(int), index=proteins)
    matrix_proteins = set(protein_df.iloc[:, 0].astype(str))
    if matrix_proteins.isdisjoint(series.index):
        raise ValueError(
            "Peptide-count sidecar has no protein identifiers in the matrix"
        )
    return series


def require_file(path: str, name: str) -> None:
    """Require absolute file paths at the MCP boundary."""
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        raise ValueError(f"{name} must be an absolute path")
    if not candidate.is_file():
        raise FileNotFoundError(f"{name} not found: {candidate}")


def output_path(path: str) -> Path:
    """Resolve one required absolute output directory."""
    destination = Path(path).expanduser()
    if not destination.is_absolute():
        raise ValueError("output_dir must be an absolute path")
    destination = destination.resolve()
    if destination.exists():
        raise FileExistsError(f"output_dir already exists: {destination}")
    return destination


@contextmanager
def atomic_output_dir(destination: Path) -> Iterator[Path]:
    """Publish a complete evaluation directory or leave no target behind."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{destination.name}.",
            dir=destination.parent,
        )
    )
    committed = False
    try:
        yield staging
        if destination.exists():
            raise FileExistsError(f"output_dir already exists: {destination}")
        staging.rename(destination)
        committed = True
    finally:
        if not committed and staging.exists():
            shutil.rmtree(staging)


def validate_contrast_samples(
    frame: pd.DataFrame,
    conditions: dict[str, str],
    contrast: tuple[str, str],
) -> None:
    """Require two mapped matrix replicates for each requested condition."""
    counts = {
        condition: sum(
            conditions[str(sample)] == condition for sample in frame.columns[1:]
        )
        for condition in contrast
    }
    if any(count < 2 for count in counts.values()):
        raise ValueError(
            "Contrast requires at least two matrix samples per condition: "
            + ", ".join(f"{condition}={count}" for condition, count in counts.items())
        )


def scope_contrast(
    frame: pd.DataFrame,
    conditions: dict[str, str],
    contrast: tuple[str, str],
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Restrict profiling and evaluation to the requested two conditions."""
    selected = [
        sample for sample in frame.columns[1:] if conditions[str(sample)] in contrast
    ]
    scoped = frame[[frame.columns[0], *selected]].copy()
    scoped_conditions = {str(sample): conditions[str(sample)] for sample in selected}
    return scoped, scoped_conditions


def as_linear(
    frame: pd.DataFrame,
    requested: str,
) -> tuple[pd.DataFrame, str]:
    """Return the linear-intensity representation required by the Rust kernel."""
    if requested == "linear":
        return frame, requested
    converted = frame.copy()
    with np.errstate(over="raise", invalid="raise"):
        converted.iloc[:, 1:] = np.exp2(converted.iloc[:, 1:].astype(float))
    return converted, requested


def load_ground_truth(path: str | None) -> set[str] | None:
    """Read one protein identifier per line when Score A truth is available."""
    if path is None:
        return None
    values = {
        line.strip()
        for line in Path(path).expanduser().read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if not values:
        raise ValueError("Ground-truth file contains no protein identifiers")
    return values


def validate_runtime_options(
    fdr_threshold: float,
    input_scale: Any,
    threads: Any,
    expected_direction: Any,
) -> None:
    """Validate user-controlled execution options before writing output."""
    if not 0.0 < fdr_threshold <= 1.0:
        raise ValueError("options.fdr_threshold must be in (0, 1]")
    validate_input_scale(input_scale, "options.input_scale")
    if isinstance(threads, bool) or not isinstance(threads, int):
        raise ValueError("options.threads must be an integer between 1 and 256")
    if not 1 <= threads <= 256:
        raise ValueError("options.threads must be an integer between 1 and 256")
    if expected_direction is not None and expected_direction not in {"UP", "DOWN"}:
        raise ValueError("options.expected_direction must be UP, DOWN, or null")


def validate_input_scale(value: Any, name: str) -> None:
    """Require an explicit matrix scale at both public MCP boundaries."""
    if not isinstance(value, str) or value not in {"linear", "log2"}:
        raise ValueError(f"{name} must be linear or log2")


def audit_fields(
    graph: KnowledgeGraph,
    generated: GeneratedRecommendationBlock,
    context: BoundContext,
    input_scale: str | None,
) -> dict[str, Any]:
    """Return provenance shared by completed and abstained evaluations."""
    return {
        "knowledge_fingerprint": graph.fingerprint,
        "evidence_refs": list(generated.evidence_refs),
        "confidence": generated.confidence,
        "limitations": list(generated.limitations),
        "input_scale": input_scale,
        "diagnostics": [item.to_dict() for item in context.diagnostics],
    }


def write_evaluation(destination: Path, payload: dict[str, Any]) -> dict[str, Any]:
    """Persist and return one strict JSON-safe evaluation artifact."""
    destination.mkdir(parents=True, exist_ok=True)
    safe_payload = json_safe(payload)
    (destination / "evaluation.json").write_text(
        json.dumps(safe_payload, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    return safe_payload


def slug(value: str) -> str:
    """Create a portable filename for one candidate result."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "candidate"


def json_safe(value: Any) -> Any:
    """Replace non-finite numeric diagnostics before JSON/MCP serialization."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value
