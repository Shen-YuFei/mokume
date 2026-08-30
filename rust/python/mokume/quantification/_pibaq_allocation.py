"""DataFrame helpers used by the Rust-backed piBAQ compatibility APIs."""

from collections.abc import Container, Hashable
from typing import List, Sequence

import pandas as pd
from pandas import DataFrame

from mokume.core.constants import (
    EVIDENCE_LEVEL,
    FAMILY_ID,
    FAMILY_SIZE,
    MOLECULARWEIGHT,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PEPTIDE_SEQUENCE,
    PIBAQ,
    PROTEIN_NAME,
    TPA,
)


def _detect_peptide_column(data: pd.DataFrame) -> str:
    """Locate the canonical peptide-sequence column in ``data``."""
    for candidate in (PEPTIDE_CANONICAL, "sequence", PEPTIDE_SEQUENCE):
        if candidate in data.columns:
            return candidate
    raise ValueError(
        "piBAQ requires a peptide-sequence column "
        f"({PEPTIDE_CANONICAL!r}, 'sequence', or {PEPTIDE_SEQUENCE!r})."
    )


def _membership_mask(values: pd.Series, container: Container[str]) -> pd.Series:
    """Test only observed values without materializing a large key domain."""
    return values.map(
        lambda value: isinstance(value, Hashable) and value in container
    ).astype(bool, copy=False)


def _empty_pibaq_frame(group_cols: Sequence[str], include_tpa: bool) -> DataFrame:
    """Return a typed empty frame matching the piBAQ output schema."""
    columns: List[str] = [
        PROTEIN_NAME,
        *group_cols,
        NORM_INTENSITY,
        PIBAQ,
        FAMILY_ID,
        EVIDENCE_LEVEL,
        FAMILY_SIZE,
    ]
    if include_tpa:
        columns.extend([MOLECULARWEIGHT, TPA])
    return pd.DataFrame(columns=columns)
