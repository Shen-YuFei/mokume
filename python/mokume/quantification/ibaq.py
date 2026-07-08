"""
iBAQ protein quantification method.

This module provides the iBAQ (intensity-Based Absolute Quantification)
method for protein quantification, including normalization and proteomic
ruler calculations.

The default quantification path is **piBAQ** (paralog-aware iBAQ, the
anchor-gated, family-rollup-fallback algorithm; :func:`compute_pibaq`).
For each protein family discovered via the FASTA-driven peptide-sharing
graph, two branches are possible:

- **Per-protein proportional allocation** when *at least one* family
  member has a detected proteotypic peptide. For each sample/condition
  group independently, every shared peptide's intensity is distributed
  across the members it maps to **in proportion to those members'
  per-group proteotypic (unique-peptide) intensity** -- so a member
  absent from a given sample receives ~none of that sample's shared
  signal (gpGrouper-style area-weighted distribution; see Saltzman et
  al. 2018 *Mol Cell Proteomics*). Each member's (proteotypic +
  allocated-shared) intensity is then divided by the number of
  theoretically observable peptides the family owns for that member
  under the cross-family razor (proteotypic AND shared), keeping iBAQ
  comparable across proteins and symmetric between numerator and
  denominator.
- **Family-level rollup** only when *no* member has a detected
  proteotypic peptide (e.g. the fully-overlapping actin family). The
  whole family receives a single iBAQ obtained from the union of family
  peptide intensities divided by the family's owned theoretical peptide
  count, replicated across members. The fallback ensures every detected
  protein family contributes to the output table even when the data
  cannot resolve members individually, rather than silently dropping
  rows as the upstream ibaqpy baseline does.

Note that the fallback gate is deliberately permissive: a family where
some members have unique anchors and others do not stays on the
proportional branch (the anchorless members fall to ~0 via the per-group
weighting) rather than collapsing every member onto one shared value.

"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from mokume.core.constants import (
    CONCENTRATION_NM,
    CONDITION,
    COPYNUMBER,
    EVIDENCE_FAMILY_ONLY,
    EVIDENCE_HIGH,
    EVIDENCE_LEVEL,
    EVIDENCE_MEDIUM,
    FAMILY_ID,
    FAMILY_SIZE,
    IBAQ,
    IBAQ_LOG,
    IBAQ_NORMALIZED,
    IBAQ_PPB,
    MOLECULARWEIGHT,
    MOLES_NMOL,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PEPTIDE_SEQUENCE,
    PROTEIN_NAME,
    SAMPLE_ID,
    TPA,
    WEIGHT_NG,
    is_parquet,
)
from mokume.core.logger import get_logger, log_execution_time, log_function_call
from mokume.io.fasta import (
    digest_fasta_full,
)
from mokume.io.fasta import extract_fasta as _extract_fasta_io
from mokume.model.organism import OrganismDescription
from mokume.plotting import is_plotting_available
from mokume.quantification.families import (
    Family,
    discover_families,
    families_by_member,
    load_families_yaml,
    merge_overrides,
)

# Proteomic Ruler constants
AVAGADRO: float = 6.02214129e23
AVERAGE_BASE_PAIR_MASS: float = 617.96

# Get a logger for this module
logger = get_logger("mokume.quantification.ibaq")


@log_function_call(logger)
def normalize_ibaq(res: DataFrame) -> DataFrame:
    """
    Normalize the ibaq values using the total ibaq of the sample.

    The resulted ibaq values are then multiplied by 100'000'000
    (PRIDE database normalization) for the ibaq ppb and log10 shifted
    by 10 (ProteomicsDB).

    Parameters
    ----------
    res : pd.DataFrame
        Input DataFrame with iBAQ values.

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized iBAQ values.
    """
    # rIBAQ: divide each protein's iBAQ by its (sample, condition) total.
    # ``groupby(...).transform`` is the vectorised, future-proof equivalent of
    # the per-group ``apply(normalize)`` it replaces (no DeprecationWarning,
    # no implicit grouping-column handling).
    group_totals = res.groupby([SAMPLE_ID, CONDITION], observed=True)[IBAQ].transform(
        "sum"
    )
    res[IBAQ_NORMALIZED] = res[IBAQ] / group_totals
    # Normalization method used by Proteomics DB 10 + log10(ibaq/sum(ibaq))
    res[IBAQ_LOG] = np.where(
        res[IBAQ_NORMALIZED] > 0, np.log10(res[IBAQ_NORMALIZED]) + 10, 0
    )
    # Normalization used by PRIDE Team (no log transformation) (ibaq/total_ibaq) * 100'000'000
    res[IBAQ_PPB] = res[IBAQ_NORMALIZED] * 100_000_000
    return res


@log_function_call(logger, level=logging.DEBUG)
def handle_nonstandard_aa(aa_seq: str):
    """
    Remove any nonstandard amino acid from the sequence.

    Parameters
    ----------
    aa_seq : str
        Protein sequence from database.

    Returns
    -------
    tuple
        A tuple containing a list of nonstandard amino acids and the cleaned sequence.
    """
    standard_aa = "ARNDBCEQZGHILKMFPSTWYV"
    nonstandard_aa_lst = [aa for aa in aa_seq if aa not in standard_aa]
    considered_seq = "".join([aa for aa in aa_seq if aa in standard_aa])
    return nonstandard_aa_lst, considered_seq


@log_function_call(logger)
def extract_fasta(
    fasta: str, enzyme: str, proteins: List, min_aa: int, max_aa: int, tpa: bool
):
    """Forward to :func:`mokume.io.fasta.extract_fasta`.

    Kept as a thin wrapper so existing callers can keep importing
    ``extract_fasta`` from :mod:`mokume.quantification.ibaq`. All the
    digestion, uniqueness, and MW logic — including the cross-protein unique
    peptide counting — lives in a single implementation in
    :mod:`mokume.io.fasta` to avoid drift between the two public entry points.
    See that function's docstring for the full contract.

    .. note::
        This is the **legacy ibaqpy-style** denominator helper. It returns a
        *proteotypic-only* (cross-protein unique) theoretical peptide count,
        which is symmetric with a *proteotypic-only numerator* (shared
        peptides discarded). It is **not** used by the default piBAQ path
        (:func:`compute_pibaq` / :func:`peptides_to_protein` /
        ``features2proteins --quant-method ibaq``), which instead pairs a
        *shared-aware numerator* (proportional shared-peptide reallocation)
        with a *total-potential* denominator (proteotypic + shared) via
        :func:`mokume.io.fasta.digest_fasta_full`. Do not mix the two
        conventions: the proteotypic-only and total-potential denominators
        are each only correct alongside their matching numerator. Retained
        for direct API consumers; prefer ``compute_pibaq`` for new code.
    """
    return _extract_fasta_io(fasta, enzyme, proteins, min_aa, max_aa, tpa)


class ConcentrationWeightByProteomicRuler:
    """
    Calculate protein copy number, moles, weight, and concentration using proteomic ruler.

    This uses a proteomic ruler approach to estimate the copy number, moles,
    and weight of proteins in a dataset based on their normalized intensity and molecular
    weight.
    """

    organism: OrganismDescription
    ploidy: int
    concentration_per_cell: float
    dna_mass: float

    def __init__(
        self, organism: OrganismDescription, ploidy: int, concentration_per_cell: float
    ):
        self.organism = organism
        self.ploidy = ploidy
        self.concentration_per_cell = concentration_per_cell
        self.dna_mass = (
            self.ploidy * self.organism.genome_size * AVERAGE_BASE_PAIR_MASS / AVAGADRO
        )

    def total_histone_intensities(self, protein_intensities: pd.DataFrame) -> float:
        histones = set(self.organism.histone_entries)
        is_histone_mask = protein_intensities[PROTEIN_NAME].isin(histones)
        histone_intensities = max(
            protein_intensities[is_histone_mask][NORM_INTENSITY].sum(), 1.0
        )
        return histone_intensities

    def apply_ruler(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        protein_intensities = protein_intensities.copy()
        histone_intensity = self.total_histone_intensities(protein_intensities)

        protein_intensities[COPYNUMBER] = (
            protein_intensities[NORM_INTENSITY]
            / histone_intensity
            * self.dna_mass
            * AVAGADRO
            / protein_intensities[MOLECULARWEIGHT]
        )

        protein_intensities[MOLES_NMOL] = protein_intensities[COPYNUMBER] * (
            1e9 / AVAGADRO
        )
        protein_intensities[WEIGHT_NG] = (
            protein_intensities[MOLES_NMOL] * protein_intensities[MOLECULARWEIGHT]
        )

        volume = (
            protein_intensities[WEIGHT_NG].sum() / 1e-9 / self.concentration_per_cell
        )
        protein_intensities[CONCENTRATION_NM] = volume * protein_intensities[MOLES_NMOL]
        return protein_intensities

    def __call__(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        return self.apply_ruler(protein_intensities)

    def apply_by_condition(self, protein_intensities: pd.DataFrame):
        protein_intensities = protein_intensities.groupby([CONDITION]).apply(self)
        return protein_intensities


# ---------------------------------------------------------------------------
# piBAQ (paralog-aware iBAQ) -- anchor-gated proportional allocation
#                              + family-rollup fallback
# ---------------------------------------------------------------------------

# Small constant added to per-member allocation weights so a peptide can still
# be split when *all* members of a family have zero unique anchors -- this
# only matters for sanity in the proportional branch boundary case; the
# fallback branch is normally chosen first.
_WEIGHT_FLOOR = 1e-9


def _detect_peptide_column(data: pd.DataFrame) -> str:
    """Locate the canonical peptide-sequence column in ``data``.

    The piBAQ algorithm needs the raw peptide sequence to look up its
    FASTA-derived accession set. The mokume pipeline emits
    :data:`PEPTIDE_CANONICAL`; the raw QPX feature parquet uses
    ``sequence``; some legacy callers use :data:`PEPTIDE_SEQUENCE`. The
    first column found in that priority order wins.
    """
    for candidate in (PEPTIDE_CANONICAL, "sequence", PEPTIDE_SEQUENCE):
        if candidate in data.columns:
            return candidate
    raise ValueError(
        "piBAQ requires a peptide-sequence column "
        f"({PEPTIDE_CANONICAL!r}, 'sequence', or {PEPTIDE_SEQUENCE!r})."
    )


def _count_unique_anchors(
    peptide_df: DataFrame,
    pep_to_accs: Mapping[str, Set[str]],
    pep_col: str,
) -> Dict[str, int]:
    """Count distinct FASTA-proteotypic peptides observed per canonical accession.

    A peptide contributes to a protein's anchor count when (a) it appears
    at least once in ``peptide_df`` and (b) the FASTA digest maps it to
    exactly that one accession. The denominator for the proportional
    branch's eligibility test is the *member's* unique-anchor count.
    """
    counts: "defaultdict[str, int]" = defaultdict(int)
    observed = set(peptide_df[pep_col].dropna().unique())
    for peptide in observed:
        accs = pep_to_accs.get(peptide)
        if accs and len(accs) == 1:
            counts[next(iter(accs))] += 1
    return dict(counts)


def _classify_evidence(min_anchor: int, min_required: int, high_threshold: int) -> str:
    """Bucket a family's minimum anchor count into a human-readable label."""
    if min_anchor < min_required:
        return EVIDENCE_FAMILY_ONLY
    if min_anchor >= high_threshold:
        return EVIDENCE_HIGH
    return EVIDENCE_MEDIUM


def _proportional_branch(
    family: Family,
    family_df: DataFrame,
    pep_to_accs: Mapping[str, Set[str]],
    acc_to_peps: Mapping[str, Set[str]],
    owned_peptides: Set[str],
    pep_col: str,
    group_cols: Sequence[str],
) -> DataFrame:
    """Per-protein iBAQ via per-sample, unique-intensity-weighted allocation.

    The output is one row per ``(member, *group_cols)`` combination. Within
    each group the shared-peptide intensity is split across the members it
    maps to **in proportion to those members' per-group proteotypic
    intensity**, so a member that is absent from a given sample receives ~0
    of that sample's shared signal (only the :data:`_WEIGHT_FLOOR`
    tie-breaker), instead of a globally fixed ratio leaking phantom signal
    into samples where the member was never detected.

    Each member's summed (proteotypic + allocated-shared) intensity is then
    divided by the number of theoretically observable peptides the family
    *owns* for that member under the cross-family razor (``owned_peptides``
    intersected with the member's digest). Restricting the denominator to
    owned peptides keeps it symmetric with the numerator, which only ever
    sums family-owned peptides, and closes the cross-family seam where a
    peptide razor-assigned to another family would otherwise inflate this
    member's denominator without contributing to its numerator.
    """
    group_keys = list(group_cols)
    member_set = set(family.members)

    # De-duplicate by (peptide, *group_cols) BEFORE allocation: an upstream
    # razor table may list the same shared peptide once per assigned protein
    # (e.g. DIA-NN's razor rows mirroring intensity onto every member of a
    # group with identical intensity). We take MAX rather than SUM because
    # razor mirror rows describe the SAME detected signal repeated.
    deduped = (
        family_df.groupby([pep_col, *group_keys], dropna=False, observed=True)[
            NORM_INTENSITY
        ]
        .max()
        .reset_index()
    )
    if deduped.empty:
        return pd.DataFrame()

    # Within-family accession set per observed peptide: size 1 -> proteotypic
    # anchor for that member; size >1 -> shared among those members.
    fam_accs = {
        pep: (pep_to_accs.get(pep, set()) & member_set)
        for pep in deduped[pep_col].unique()
    }
    deduped = deduped[deduped[pep_col].map(lambda p: len(fam_accs[p]) > 0)]
    if deduped.empty:
        return pd.DataFrame()

    n_members = deduped[pep_col].map(lambda p: len(fam_accs[p]))
    anchor_rows = deduped[n_members == 1].copy()
    shared_rows = deduped[n_members > 1].copy()

    # Per-(member, group) proteotypic intensity -- the within-family anchor
    # signal that both drives the shared split and forms the numerator base.
    if not anchor_rows.empty:
        anchor_rows[PROTEIN_NAME] = anchor_rows[pep_col].map(
            lambda p: next(iter(fam_accs[p]))
        )
        anchor_int = (
            anchor_rows.groupby(
                [PROTEIN_NAME, *group_keys], dropna=False, observed=True
            )[NORM_INTENSITY]
            .sum()
            .reset_index()
        )
    else:
        anchor_int = pd.DataFrame(columns=[PROTEIN_NAME, *group_keys, NORM_INTENSITY])

    # Allocate each shared peptide WITHIN EACH GROUP by the members' per-group
    # anchor intensity (absent members weigh ~0 via _WEIGHT_FLOOR).
    if not shared_rows.empty:
        shared_rows = shared_rows[[pep_col, *group_keys, NORM_INTENSITY]].copy()
        shared_rows["_members"] = shared_rows[pep_col].map(
            lambda p: sorted(fam_accs[p])
        )
        exploded = shared_rows.explode("_members", ignore_index=True).rename(
            columns={"_members": PROTEIN_NAME}
        )
        exploded = exploded.merge(
            anchor_int.rename(columns={NORM_INTENSITY: "_anchor_int"}),
            on=[PROTEIN_NAME, *group_keys],
            how="left",
        )
        exploded["_w"] = exploded["_anchor_int"].fillna(0.0) + _WEIGHT_FLOOR
        exploded["_wsum"] = exploded.groupby(
            [pep_col, *group_keys], dropna=False, observed=True
        )["_w"].transform("sum")
        exploded["_alloc"] = (
            exploded[NORM_INTENSITY] * exploded["_w"] / exploded["_wsum"]
        )
        shared_contrib = (
            exploded.groupby([PROTEIN_NAME, *group_keys], dropna=False, observed=True)[
                "_alloc"
            ]
            .sum()
            .reset_index()
            .rename(columns={"_alloc": NORM_INTENSITY})
        )
    else:
        shared_contrib = pd.DataFrame(
            columns=[PROTEIN_NAME, *group_keys, NORM_INTENSITY]
        )

    # Concatenate only the non-empty branch frames. The empty placeholders
    # above are object-dtype (``pd.DataFrame(columns=...)`` has no values to
    # infer from); including one in the concat lets newer pandas (3.0 dropped
    # the legacy "exclude empty/all-NA entries when inferring result dtypes"
    # rule) poison ``NormIntensity`` to object. That object dtype then flows
    # through the groupby-sum into an object-dtype iBAQ and breaks the
    # ``np.log10`` step in :func:`normalize_ibaq`.
    parts = [frame for frame in (anchor_int, shared_contrib) if not frame.empty]
    if not parts:
        return pd.DataFrame()
    combined = pd.concat(parts, ignore_index=True)
    grouped = (
        combined.groupby([PROTEIN_NAME, *group_keys], dropna=False, observed=True)[
            NORM_INTENSITY
        ]
        .sum()
        .reset_index()
    )

    denominators = {
        m: len(acc_to_peps.get(m, set()) & owned_peptides) for m in family.members
    }
    denom = grouped[PROTEIN_NAME].map(denominators).astype(float)
    grouped[IBAQ] = grouped[NORM_INTENSITY] / denom.where(denom > 0)
    return grouped


def _family_rollup_branch(
    family: Family,
    family_df: DataFrame,
    owned_peptides: Set[str],
    pep_col: str,
    group_cols: Sequence[str],
) -> DataFrame:
    """Family-level iBAQ shared identically across every member.

    Each peptide-sample observation contributes once (we de-duplicate by
    ``(peptide, *group_cols)`` before summing), preventing the double
    counting that would arise if a shared peptide appeared under several
    protein rows in the upstream razor table. The denominator is the number
    of theoretical peptides the family owns under the cross-family razor
    (``len(owned_peptides)``), symmetric with the family-summed numerator.
    """
    denominator = len(owned_peptides) if owned_peptides else 1
    # MAX deduplicates razor mirror rows; see ``_proportional_branch`` for
    # the rationale. After per-peptide collapse we sum across all family
    # peptides to obtain the per-sample family intensity total.
    deduped_per_peptide = (
        family_df.groupby([pep_col, *group_cols], dropna=False, observed=True)[
            NORM_INTENSITY
        ]
        .max()
        .reset_index()
    )
    family_total = (
        deduped_per_peptide.groupby(list(group_cols), dropna=False, observed=True)[
            NORM_INTENSITY
        ]
        .sum()
        .reset_index()
    )
    family_total[IBAQ] = family_total[NORM_INTENSITY] / denominator

    # Replicate the family-level value for every member so callers retain a
    # row per protein. ``DataFrame.merge`` on a cross key is the vectorised
    # equivalent of a nested loop and stays fast even for large samples.
    members_df = pd.DataFrame({PROTEIN_NAME: list(family.members)})
    members_df["_key"] = 1
    family_total["_key"] = 1
    replicated = family_total.merge(members_df, on="_key").drop(columns="_key")
    cols = [PROTEIN_NAME, *group_cols, NORM_INTENSITY, IBAQ]
    return replicated[cols]


def _assign_peptides_to_owning_family(
    families: Sequence[Family],
    pep_to_accs: Mapping[str, Set[str]],
    anchor_counts: Mapping[str, int],
) -> Dict[str, str]:
    """Cross-family razor: every observed peptide gets exactly one owner.

    A peptide that maps to proteins in more than one family would otherwise
    be processed (and its intensity claimed) by every touching family,
    inflating the global iBAQ total. We resolve the ambiguity once, up
    front, by assigning the peptide to the family whose strongest member
    (largest unique-anchor count) contains it. Ties on anchor count fall
    back to lexicographic family_id so the assignment is reproducible.
    """
    member_to_family = families_by_member(families)
    owner: Dict[str, str] = {}
    for peptide, accs in pep_to_accs.items():
        best: Optional[Tuple[int, str]] = None
        for acc in accs:
            family = member_to_family.get(acc)
            if family is None:
                continue
            score = anchor_counts.get(acc, 0)
            candidate = (-score, family.family_id)
            if best is None or candidate < best:
                best = candidate
        if best is not None:
            owner[peptide] = best[1]
    return owner


def _invert_peptide_ownership(
    peptide_owner: Mapping[str, str],
) -> Dict[str, Set[str]]:
    """Group peptides by their owning family in a single pass.

    Building this inverted index once -- O(N_peptides) -- avoids the
    quadratic cost of scanning ``peptide_owner`` inside every family's
    :func:`_assemble_family_input` call. On a human FASTA with ~975k
    peptides and ~17k families that loop is the dominant runtime.
    """
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
    """Restrict ``peptide_df`` to rows owned by this family.

    Each peptide is processed by exactly one family thanks to the upstream
    razor assignment in :func:`_assign_peptides_to_owning_family`; this
    function simply selects the subset owned by the requested family,
    using the pre-built ``family_to_peptides`` index for O(1) lookup.
    """
    owned_peptides = family_to_peptides.get(family.family_id, set())
    if not owned_peptides:
        return peptide_df.iloc[0:0]
    return peptide_df[peptide_df[pep_col].isin(owned_peptides)]


def _annotate_family_metadata(
    block: DataFrame,
    family: Family,
    evidence: str,
) -> DataFrame:
    """Attach the (FamilyId, EvidenceLevel, FamilySize) columns in one pass."""
    block = block.copy()
    block[FAMILY_ID] = family.family_id
    block[EVIDENCE_LEVEL] = evidence
    block[FAMILY_SIZE] = family.size
    return block


def _attach_tpa(
    block: DataFrame,
    family: Family,
    branch_is_fallback: bool,
    mw_map: Mapping[str, float],
) -> DataFrame:
    """Compute TPA per the rule in plan section 2.4.

    Per-protein branch: ``TPA = NormIntensity / MW(member)``.
    Fallback branch: ``TPA = family_intensity / sum(MW for member in family)``,
    so every replicated family row reports the same TPA -- mirroring the
    fact that the iBAQ is shared across the family.
    """
    if branch_is_fallback:
        family_mw = sum(mw_map.get(m, 0.0) for m in family.members)
        if family_mw <= 0:
            family_mw = 1.0
        block[MOLECULARWEIGHT] = family_mw
    else:
        mw_series = block[PROTEIN_NAME].map(mw_map).fillna(0.0)
        block[MOLECULARWEIGHT] = mw_series.replace(0.0, 1.0)
    block[TPA] = block[NORM_INTENSITY] / block[MOLECULARWEIGHT]
    return block


def _empty_pibaq_frame(group_cols: Sequence[str], include_tpa: bool) -> DataFrame:
    """Return a typed empty frame matching the piBAQ output schema."""
    columns: List[str] = [
        PROTEIN_NAME,
        *group_cols,
        NORM_INTENSITY,
        IBAQ,
        FAMILY_ID,
        EVIDENCE_LEVEL,
        FAMILY_SIZE,
    ]
    if include_tpa:
        columns.extend([MOLECULARWEIGHT, TPA])
    return pd.DataFrame(columns=columns)


def compute_pibaq(
    peptide_df: DataFrame,
    accession_to_peptides: Mapping[str, Set[str]],
    peptide_to_accessions: Mapping[str, Set[str]],
    families: Sequence[Family],
    *,
    mw_map: Optional[Mapping[str, float]] = None,
    min_anchors: int = 1,
    high_anchor_threshold: int = 3,
    extra_group_cols: Optional[Sequence[str]] = None,
) -> DataFrame:
    """Compute per-protein iBAQ with anchor-gated family fallback.

    Parameters
    ----------
    peptide_df : DataFrame
        Long-format peptide table with at minimum the columns
        :data:`PROTEIN_NAME`, :data:`SAMPLE_ID`, :data:`NORM_INTENSITY`,
        and one of (:data:`PEPTIDE_CANONICAL`, ``'sequence'``,
        :data:`PEPTIDE_SEQUENCE`). :data:`CONDITION` is honoured when
        present; additional grouping columns may be requested via
        ``extra_group_cols``.
    accession_to_peptides, peptide_to_accessions : Mapping
        FASTA digest indices, typically produced by
        :func:`mokume.io.fasta.digest_fasta_full`.
    families : Sequence[Family]
        Output of :func:`mokume.quantification.families.discover_families`
        (optionally post-processed by
        :func:`mokume.quantification.families.merge_overrides`).
    mw_map : Mapping[str, float], optional
        Per-canonical molecular weight (only consumed when callers request
        TPA via :func:`peptides_to_protein` ``tpa=True``).
    min_anchors : int, optional
        Minimum unique-anchor count *per family member* required to stay
        on the proportional branch. Defaults to 1 -- the gpGrouper-style
        threshold (Saltzman et al. 2018 MCP 17:2270).
    high_anchor_threshold : int, optional
        Threshold (in unique anchors of the *weakest* family member) at
        which ``EvidenceLevel`` is reported as ``high``. Defaults to 3.
    extra_group_cols : Sequence[str], optional
        Additional dataframe columns to retain in the groupby keys, on
        top of (SampleID, Condition).

    Returns
    -------
    DataFrame
        Long-format result containing one row per ``(member,
        SampleID, Condition, *extra_group_cols)``. Includes the new
        :data:`FAMILY_ID`, :data:`EVIDENCE_LEVEL`, :data:`FAMILY_SIZE`
        metadata columns, and :data:`MOLECULARWEIGHT` + :data:`TPA` when
        ``mw_map`` is supplied.
    """
    if peptide_df.empty or not families:
        include_tpa = mw_map is not None
        group_cols = _resolve_group_cols(peptide_df, extra_group_cols)
        return _empty_pibaq_frame(group_cols, include_tpa)

    pep_col = _detect_peptide_column(peptide_df)
    group_cols = _resolve_group_cols(peptide_df, extra_group_cols)

    # Keep only peptides present in the FASTA digest; matching is driven by
    # the peptide-sequence column against the digest index, not by the
    # per-row protein-name strings.
    working = peptide_df[peptide_df[pep_col].isin(peptide_to_accessions)]
    if working.empty:
        include_tpa = mw_map is not None
        return _empty_pibaq_frame(group_cols, include_tpa)

    anchor_counts = _count_unique_anchors(working, peptide_to_accessions, pep_col)
    peptide_owner = _assign_peptides_to_owning_family(
        families, peptide_to_accessions, anchor_counts
    )
    family_to_peptides = _invert_peptide_ownership(peptide_owner)
    # Restrict to peptides that any family claims, then further to families
    # with at least one observed peptide. On a 25k-protein human FASTA the
    # vast majority of families are isolated singletons with no detection,
    # so iterating only the subset with data drops the per-call hot loop
    # from O(|families|) to O(|families with data|).
    all_owned_peptides: Set[str] = set()
    for peps in family_to_peptides.values():
        all_owned_peptides |= peps
    working = working[working[pep_col].isin(all_owned_peptides)]
    observed_peptides = set(working[pep_col].dropna().unique())

    results: List[DataFrame] = []
    n_fallback = 0
    for family in families:
        owned = family_to_peptides.get(family.family_id, set())
        if not owned & observed_peptides:
            continue
        family_input = _assemble_family_input(
            working, family, family_to_peptides, pep_col
        )
        if family_input.empty:
            continue
        anchors_in_family = {m: anchor_counts.get(m, 0) for m in family.members}
        min_anchor = min(anchors_in_family.values())
        max_anchor = max(anchors_in_family.values())
        # Fall back to a single family-level iBAQ only when NO member has a
        # proteotypic anchor (e.g. the fully-overlapping actin family). When
        # some members are resolvable, stay on the proportional branch: the
        # anchorless members fall to ~0 via the per-group intensity weighting
        # rather than collapsing every member onto one shared value.
        is_fallback = max_anchor < min_anchors
        if is_fallback:
            block = _family_rollup_branch(
                family,
                family_input,
                owned,
                pep_col,
                group_cols,
            )
            n_fallback += 1
        else:
            block = _proportional_branch(
                family,
                family_input,
                peptide_to_accessions,
                accession_to_peptides,
                owned,
                pep_col,
                group_cols,
            )
        if block.empty:
            continue
        evidence = _classify_evidence(min_anchor, min_anchors, high_anchor_threshold)
        block = _annotate_family_metadata(block, family, evidence)
        if mw_map is not None:
            block = _attach_tpa(block, family, is_fallback, mw_map)
        results.append(block)

    if not results:
        include_tpa = mw_map is not None
        return _empty_pibaq_frame(group_cols, include_tpa)

    out = pd.concat(results, ignore_index=True)
    n_processed = len(results)
    logger.info(
        "piBAQ: %d families processed (%d fallback to family rollup, %d per-protein)",
        n_processed,
        n_fallback,
        n_processed - n_fallback,
    )
    return out


def _resolve_group_cols(
    peptide_df: DataFrame,
    extra_group_cols: Optional[Sequence[str]],
) -> List[str]:
    """Build the canonical groupby key list, dropping absent columns."""
    cols: List[str] = [SAMPLE_ID]
    if CONDITION in peptide_df.columns:
        cols.append(CONDITION)
    if extra_group_cols:
        for col in extra_group_cols:
            if col in peptide_df.columns and col not in cols:
                cols.append(col)
    return cols


def _resolve_families(
    accession_to_peptides: Mapping[str, Set[str]],
    peptide_to_accessions: Mapping[str, Set[str]],
    *,
    families_yaml: Optional[Path],
    min_shared: int,
) -> List[Family]:
    """Combine automatic CC discovery with optional YAML overrides."""
    auto = discover_families(
        accession_to_peptides, peptide_to_accessions, min_shared=min_shared
    )
    if families_yaml is None:
        return auto
    overrides = load_families_yaml(Path(families_yaml))
    return merge_overrides(auto, overrides)


def _compute_pibaq_table(
    data: DataFrame,
    fasta: str,
    enzyme: str,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    *,
    families_yaml: Optional[Path],
    min_shared: int,
    min_anchors: int = 1,
    high_anchor_threshold: int = 3,
) -> DataFrame:
    """piBAQ (paralog-aware iBAQ) driver -- the default code path.

    Digests the FASTA once (with UniProt isoform collapse), discovers
    families on the shared-peptide graph, merges in any YAML overrides,
    then dispatches to :func:`compute_pibaq`. Returns a long-format
    table with the new ``FamilyId`` / ``EvidenceLevel`` / ``FamilySize``
    metadata columns appended.

    Unlike the legacy baseline this driver does NOT pre-strip
    ``unique != 1`` rows: the family-rollup fallback for indistinguishable
    paralogs (e.g. actin) needs the shared-peptide intensity to compute a
    family-level iBAQ, and the proportional branch likewise re-derives
    its own allocation from the FASTA digest. Razor-mirror rows that
    list the same peptide observation under multiple proteins are
    de-duplicated inside the per-family branches via
    ``groupby(peptide, *group_cols).max``.
    """
    accession_to_peptides, peptide_to_accessions, accession_to_mw = digest_fasta_full(
        fasta,
        enzyme,
        min_aa,
        max_aa,
        canonicalize_isoforms=True,
        compute_mw=tpa,
    )

    families = _resolve_families(
        accession_to_peptides,
        peptide_to_accessions,
        families_yaml=families_yaml,
        min_shared=min_shared,
    )

    mw_map: Optional[Mapping[str, float]] = accession_to_mw if tpa else None
    res = compute_pibaq(
        data,
        accession_to_peptides,
        peptide_to_accessions,
        families,
        mw_map=mw_map,
        min_anchors=min_anchors,
        high_anchor_threshold=high_anchor_threshold,
    )
    return res


@log_execution_time(logger)
def peptides_to_protein(
    fasta: str,
    peptides: str,
    enzyme: str,
    normalize: bool,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    ruler: bool,
    ploidy: int,
    cpc: float,
    organism: str,
    output: str,
    verbose: bool,
    qc_report: str,
    *,
    families_yaml: Optional[str] = None,
    min_shared: int = 2,
    min_anchors: int = 1,
    high_anchor_threshold: int = 3,
) -> None:
    """
    Compute iBAQ values for peptides and generate a QC report.

    Parameters
    ----------
    fasta : str
        Fasta file used to perform the peptide identification.
    peptides : str
        Peptide intensity file.
    enzyme : str
        Enzyme used to digest the protein sample.
    normalize : bool
        Use normalization steps.
    min_aa : int
        Minimum number of amino acids to consider a peptide.
    max_aa : int
        Maximum number of amino acids to consider a peptide.
    tpa : bool
        Calculate TPA values.
    ruler : bool
        Calculate protein weight and concentration using a proteomic ruler approach.
    ploidy : int
        Ploidy of the organism.
    cpc : float
        Concentration per cell.
    organism : str
        Organism name.
    output : str
        Output file path.
    verbose : bool
        Print additional information.
    qc_report : str
        PDF file to store multiple QC images.
    families_yaml : str, optional
        Path to a YAML file declaring explicit family overrides used by
        :func:`mokume.quantification.families.load_families_yaml`. When
        ``None`` (default) family discovery is purely data-driven.
    min_shared : int, optional
        Minimum number of distinct peptides two proteins must share to be
        placed in the same automatically discovered family. Defaults to 2.
    min_anchors : int, optional
        Minimum proteotypic ("anchor") peptides a family member needs for
        the family to stay on the per-protein proportional branch. The
        family falls back to a single rolled-up iBAQ only when *no* member
        reaches this threshold. Defaults to 1.
    high_anchor_threshold : int, optional
        Minimum anchor count (of the weakest member) for a family to be
        labelled ``EvidenceLevel == "high"``. Defaults to 3.
    """
    if organism:
        organism_descr = OrganismDescription.get(organism)
        if organism_descr is None:
            raise KeyError(f"Could not resolve organism description for {organism}")
    else:
        organism_descr = None

    if ruler:
        if not ploidy or not cpc or not organism or not tpa:
            raise ValueError(
                "Arguments `ploidy`, `cpc`, `organism` and `tpa` are required for calculate protein weight(ng) and concentration(nM)"
            )

    # load data
    if is_parquet(peptides):
        data = pd.read_parquet(peptides)
    else:
        data = pd.read_csv(peptides)
    data[NORM_INTENSITY] = data[NORM_INTENSITY].astype(float)
    data = data.dropna(subset=[NORM_INTENSITY])
    data = data[data[NORM_INTENSITY] > 0]

    res = _compute_pibaq_table(
        data,
        fasta,
        enzyme,
        min_aa,
        max_aa,
        tpa,
        families_yaml=Path(families_yaml) if families_yaml else None,
        min_shared=min_shared,
        min_anchors=min_anchors,
        high_anchor_threshold=high_anchor_threshold,
    )

    # normalize ibaq
    if normalize:
        res = normalize_ibaq(res)
        res = res.dropna(subset=[IBAQ_NORMALIZED])
        plot_column = IBAQ_PPB
    else:
        res = res.dropna(subset=[IBAQ])
        plot_column = IBAQ

    res = res.reset_index(drop=True)

    # calculate protein weight and concentration
    if ruler:
        concentration_by_ruler = ConcentrationWeightByProteomicRuler(
            organism_descr, ploidy, cpc
        )
        res = concentration_by_ruler.apply_by_condition(res)

    # Print the distribution of the protein IBAQ values
    if verbose:
        if not is_plotting_available():
            logger.warning(
                "QC report skipped: plotting dependencies not installed. "
                "Install with: pip install mokume[plotting]"
            )
        else:
            from mokume.plotting import PdfPages, plot_box_plot, plot_distributions

            plot_width = len(set(res[SAMPLE_ID])) * 0.5 + 10
            pdf = PdfPages(qc_report)
            density1 = plot_distributions(
                res,
                plot_column,
                SAMPLE_ID,
                log2=True,
                width=plot_width,
                title="{} Distribution".format(plot_column),
            )
            box1 = plot_box_plot(
                res,
                plot_column,
                SAMPLE_ID,
                log2=True,
                width=plot_width,
                title="{} Distribution".format(plot_column),
                violin=False,
            )
            pdf.savefig(density1, bbox_inches="tight")
            pdf.savefig(box1, bbox_inches="tight")
            if tpa:
                density2 = plot_distributions(
                    res,
                    TPA,
                    SAMPLE_ID,
                    log2=True,
                    width=plot_width,
                    title="TPA Distribution",
                )
                box2 = plot_box_plot(
                    res,
                    TPA,
                    SAMPLE_ID,
                    log2=True,
                    width=plot_width,
                    title="{} Distribution".format(TPA),
                    violin=False,
                )
                pdf.savefig(density2, bbox_inches="tight")
                pdf.savefig(box2, bbox_inches="tight")
            if ruler:
                density3 = plot_distributions(
                    res,
                    COPYNUMBER,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{} Distribution".format(COPYNUMBER),
                )
                box3 = plot_box_plot(
                    res,
                    COPYNUMBER,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{} Distribution".format(COPYNUMBER),
                    violin=False,
                )
                pdf.savefig(density3, bbox_inches="tight")
                pdf.savefig(box3, bbox_inches="tight")
                density4 = plot_distributions(
                    res,
                    CONCENTRATION_NM,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{} Distribution".format(CONCENTRATION_NM),
                )
                box4 = plot_box_plot(
                    res,
                    CONCENTRATION_NM,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{} Distribution".format(CONCENTRATION_NM),
                    violin=False,
                )
                pdf.savefig(density4, bbox_inches="tight")
                pdf.savefig(box4, bbox_inches="tight")
            pdf.close()

    res.to_csv(output, sep="\t", index=False)
