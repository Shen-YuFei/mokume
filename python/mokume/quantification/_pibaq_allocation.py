"""Family ownership and shared-peptide allocation for piBAQ."""

from collections import defaultdict
from collections.abc import Container, Hashable
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import pandas as pd
from pandas import DataFrame

from mokume.core.constants import (
    EVIDENCE_FAMILY_ONLY,
    EVIDENCE_HIGH,
    EVIDENCE_LEVEL,
    EVIDENCE_MEDIUM,
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
from mokume.quantification.families import Family, families_by_member


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


def _count_unique_anchors(
    peptide_df: DataFrame,
    pep_to_accs: Mapping[str, Set[str]],
    pep_col: str,
) -> Dict[str, int]:
    """Count distinct FASTA-proteotypic peptides observed per accession."""
    counts: "defaultdict[str, int]" = defaultdict(int)
    observed = set(peptide_df[pep_col].dropna().unique())
    for peptide in observed:
        accs = pep_to_accs.get(peptide)
        if accs and len(accs) == 1:
            counts[next(iter(accs))] += 1
    return dict(counts)


def _classify_evidence(
    min_anchor: int,
    max_anchor: int,
    min_required: int,
    high_threshold: int,
) -> str:
    """Bucket a family's anchor support into a human-readable label."""
    if max_anchor < min_required:
        return EVIDENCE_FAMILY_ONLY
    if min_anchor >= high_threshold:
        return EVIDENCE_HIGH
    return EVIDENCE_MEDIUM


@dataclass(frozen=True)
class _FamilyAllocationInputs:
    family: Family
    peptide_to_accessions: Mapping[str, Set[str]]
    accession_to_peptides: Mapping[str, Set[str]]
    owned_peptides: Set[str]
    force_equal_shared: bool


def _split_family_rows(
    inputs: _FamilyAllocationInputs,
    family_df: DataFrame,
    peptide_column: str,
    group_columns: Sequence[str],
) -> Tuple[DataFrame, DataFrame, Dict[str, Set[str]]]:
    """Deduplicate observations and separate proteotypic from shared rows."""
    deduped = (
        family_df.groupby(
            [peptide_column, *group_columns], dropna=False, observed=True
        )[NORM_INTENSITY]
        .max()
        .reset_index()
    )
    if deduped.empty:
        return pd.DataFrame(), pd.DataFrame(), {}

    member_set = set(inputs.family.members)
    family_accessions = {
        peptide: (inputs.peptide_to_accessions.get(peptide, set()) & member_set)
        for peptide in deduped[peptide_column].unique()
    }
    deduped = deduped[
        deduped[peptide_column].map(lambda peptide: bool(family_accessions[peptide]))
    ]
    member_counts = deduped[peptide_column].map(
        lambda peptide: len(family_accessions[peptide])
    )
    return (
        deduped[member_counts == 1].copy(),
        deduped[member_counts > 1].copy(),
        family_accessions,
    )


def _sum_anchor_intensities(
    anchor_rows: DataFrame,
    family_accessions: Mapping[str, Set[str]],
    peptide_column: str,
    group_columns: Sequence[str],
) -> DataFrame:
    """Sum proteotypic peptide intensity for each protein and group."""
    if anchor_rows.empty:
        return pd.DataFrame(columns=[PROTEIN_NAME, *group_columns, NORM_INTENSITY])

    anchor_rows[PROTEIN_NAME] = anchor_rows[peptide_column].map(
        lambda peptide: next(iter(family_accessions[peptide]))
    )
    return (
        anchor_rows.groupby(
            [PROTEIN_NAME, *group_columns], dropna=False, observed=True
        )[NORM_INTENSITY]
        .sum()
        .reset_index()
    )


def _sum_shared_intensities(
    shared_rows: DataFrame,
    family_accessions: Mapping[str, Set[str]],
    anchor_intensities: DataFrame,
    allocation_columns: Sequence[str],
    force_equal_shared: bool,
) -> DataFrame:
    """Allocate shared signal and sum it for each protein and group."""
    peptide_column = allocation_columns[0]
    group_columns = allocation_columns[1:]
    if shared_rows.empty:
        return pd.DataFrame(columns=[PROTEIN_NAME, *group_columns, NORM_INTENSITY])

    shared_rows = shared_rows[[peptide_column, *group_columns, NORM_INTENSITY]].copy()
    shared_rows["_members"] = shared_rows[peptide_column].map(
        lambda peptide: sorted(family_accessions[peptide])
    )
    exploded = shared_rows.explode("_members", ignore_index=True).rename(
        columns={"_members": PROTEIN_NAME}
    )
    exploded = exploded.merge(
        anchor_intensities.rename(columns={NORM_INTENSITY: "_anchor_int"}),
        on=[PROTEIN_NAME, *group_columns],
        how="left",
    )
    exploded["_w"] = exploded["_anchor_int"].astype(float).fillna(0.0)
    allocation_groups = [peptide_column, *group_columns]
    exploded["_wsum"] = exploded.groupby(
        allocation_groups, dropna=False, observed=True
    )["_w"].transform("sum")
    exploded["_member_count"] = exploded.groupby(
        allocation_groups, dropna=False, observed=True
    )[PROTEIN_NAME].transform("size")
    proportional = (exploded["_wsum"] > 0) & (not force_equal_shared)
    exploded["_alloc"] = exploded[NORM_INTENSITY] / exploded["_member_count"]
    exploded.loc[proportional, "_alloc"] = (
        exploded.loc[proportional, NORM_INTENSITY]
        * exploded.loc[proportional, "_w"]
        / exploded.loc[proportional, "_wsum"]
    )
    return (
        exploded.groupby([PROTEIN_NAME, *group_columns], dropna=False, observed=True)[
            "_alloc"
        ]
        .sum()
        .reset_index()
        .rename(columns={"_alloc": NORM_INTENSITY})
    )


def _allocate_family(
    inputs: _FamilyAllocationInputs,
    family_df: DataFrame,
    pep_col: str,
    group_cols: Sequence[str],
) -> DataFrame:
    """Allocate one family's shared signal and calculate member piBAQ values."""
    anchor_rows, shared_rows, family_accessions = _split_family_rows(
        inputs, family_df, pep_col, group_cols
    )
    if not family_accessions:
        return pd.DataFrame()

    anchor_int = _sum_anchor_intensities(
        anchor_rows, family_accessions, pep_col, group_cols
    )
    shared_contrib = _sum_shared_intensities(
        shared_rows,
        family_accessions,
        anchor_int,
        [pep_col, *group_cols],
        inputs.force_equal_shared,
    )
    parts = [frame for frame in (anchor_int, shared_contrib) if not frame.empty]
    if not parts:
        return pd.DataFrame()
    combined = pd.concat(parts, ignore_index=True)
    grouped = (
        combined.groupby([PROTEIN_NAME, *group_cols], dropna=False, observed=True)[
            NORM_INTENSITY
        ]
        .sum()
        .reset_index()
    )
    grouped = grouped[grouped[NORM_INTENSITY] > 0].copy()
    denominators = {
        member: len(
            inputs.accession_to_peptides.get(member, set()) & inputs.owned_peptides
        )
        for member in inputs.family.members
    }
    denom = grouped[PROTEIN_NAME].map(denominators).astype(float)
    grouped[PIBAQ] = grouped[NORM_INTENSITY] / denom.where(denom > 0)
    return grouped


def _assign_peptides_to_owning_family(
    families: Sequence[Family],
    pep_to_accs: Mapping[str, Set[str]],
    anchor_counts: Mapping[str, int],
) -> Dict[str, str]:
    """Assign every peptide to one family with a deterministic razor rule."""
    member_to_family = families_by_member(families)
    owner: Dict[str, str] = {}
    for peptide, accs in pep_to_accs.items():
        best: Optional[Tuple[int, str]] = None
        for acc in accs:
            family = member_to_family.get(acc)
            if family is None:
                continue
            candidate = (-anchor_counts.get(acc, 0), family.family_id)
            if best is None or candidate < best:
                best = candidate
        if best is not None:
            owner[peptide] = best[1]
    return owner


def _invert_peptide_ownership(
    peptide_owner: Mapping[str, str],
) -> Dict[str, Set[str]]:
    """Group peptides by their owning family in a single pass."""
    inverted: "defaultdict[str, Set[str]]" = defaultdict(set)
    for peptide, owner in peptide_owner.items():
        inverted[owner].add(peptide)
    return dict(inverted)


def _assemble_family_input(
    peptide_df: DataFrame,
    family: Family,
    family_to_peptides: Mapping[str, Set[str]],
    pep_col: str,
) -> DataFrame:
    """Restrict the observed peptide table to rows owned by one family."""
    owned_peptides = family_to_peptides.get(family.family_id, set())
    if not owned_peptides:
        return peptide_df.iloc[0:0]
    return peptide_df[peptide_df[pep_col].isin(owned_peptides)]


def _annotate_family_metadata(
    block: DataFrame,
    family: Family,
    evidence: str,
) -> DataFrame:
    """Attach family identity, evidence, and size columns."""
    block = block.copy()
    block[FAMILY_ID] = family.family_id
    block[EVIDENCE_LEVEL] = evidence
    block[FAMILY_SIZE] = family.size
    return block


def _finalize_tpa(out: DataFrame, mw_map: Mapping[str, float]) -> DataFrame:
    """Attach molecular weights and TPA once to the combined result."""
    molecular_weights = out[PROTEIN_NAME].map(mw_map).fillna(0.0)
    out[MOLECULARWEIGHT] = molecular_weights.replace(0.0, 1.0)
    out[TPA] = out[NORM_INTENSITY] / out[MOLECULARWEIGHT]
    return out


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
