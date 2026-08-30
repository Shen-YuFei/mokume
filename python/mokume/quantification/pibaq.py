"""Paralog-aware intensity-Based Absolute Quantification (piBAQ).

piBAQ retains the theoretical-peptide scaling of iBAQ while making shared
peptide handling explicit and auditable.

The quantification path is **piBAQ** (paralog-aware iBAQ;
:func:`compute_pibaq`). For each sample/condition group independently,
every shared peptide's intensity is distributed across the family members
it maps to using the gpGrouper area rule (Saltzman et al. 2018 *Mol Cell
Proteomics*):

- when at least one mapped member has proteotypic (unique-peptide)
  intensity, allocation is proportional to those per-group intensities;
  a member with zero proteotypic intensity receives exactly zero;
- when every mapped member has zero proteotypic intensity, the shared
  intensity is split equally among them.

The allocation conserves each observed shared-peptide intensity. Each
member's (proteotypic + allocated-shared) intensity is then divided by the
number of theoretically observable peptides the family owns for that
member under the cross-family razor (proteotypic AND shared), keeping the
numerator and denominator symmetric.

"""

import importlib
import logging
from dataclasses import dataclass
from typing import (
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

import numpy as np
import pandas as pd
from pandas import DataFrame

from mokume.core.constants import (
    CONCENTRATION_NM,
    CONDITION,
    COPYNUMBER,
    EVIDENCE_FAMILY_ONLY,
    PIBAQ,
    PIBAQ_LOG,
    PIBAQ_NORMALIZED,
    PIBAQ_PPB,
    MOLECULARWEIGHT,
    MOLES_NMOL,
    NORM_INTENSITY,
    PROTEIN_NAME,
    SAMPLE_ID,
    WEIGHT_NG,
)
from mokume.core.logger import get_logger, log_function_call
from mokume.io.fasta import extract_fasta as _extract_fasta_io
from mokume.model.organism import OrganismDescription
from mokume.quantification._pibaq_allocation import (
    _allocate_family,
    _annotate_family_metadata,
    _assemble_family_input,
    _assign_peptides_to_owning_family,
    _classify_evidence,
    _count_unique_anchors,
    _detect_peptide_column,
    _empty_pibaq_frame,
    _FamilyAllocationInputs,
    _finalize_tpa,
    _invert_peptide_ownership,
    _membership_mask,
)
from mokume.quantification.families import Family


def peptides_to_protein(*args: object, **kwargs: object) -> None:
    """Run the file-oriented workflow without coupling it to module import."""
    workflow = importlib.import_module("mokume.quantification._pibaq_workflow")
    getattr(workflow, "peptides_to_protein")(*args, **kwargs)


# Proteomic Ruler constants
AVAGADRO: float = 6.02214129e23
AVERAGE_BASE_PAIR_MASS: float = 617.96

# Get a logger for this module
logger = get_logger("mokume.quantification.pibaq")


def _bind_arguments(
    function_name: str,
    args: Sequence[object],
    kwargs: Mapping[str, object],
    positional_names: Sequence[str],
    defaults: Mapping[str, object],
) -> dict[str, object]:
    """Bind a legacy call without changing its accepted argument forms."""
    if len(args) > len(positional_names):
        raise TypeError(
            f"{function_name}() takes {len(positional_names)} positional arguments "
            f"but {len(args)} were given"
        )
    values = dict(zip(positional_names, args))
    duplicate = next((name for name in values if name in kwargs), None)
    if duplicate is not None:
        raise TypeError(f"{function_name}() got multiple values for '{duplicate}'")
    allowed = set(positional_names) | set(defaults)
    unexpected = next((name for name in kwargs if name not in allowed), None)
    if unexpected is not None:
        raise TypeError(
            f"{function_name}() got an unexpected keyword argument '{unexpected}'"
        )
    values.update(kwargs)
    missing = [name for name in positional_names if name not in values]
    if missing:
        raise TypeError(f"{function_name}() missing required argument: '{missing[0]}'")
    for name, default in defaults.items():
        values.setdefault(name, default)
    return values


@log_function_call(logger)
def normalize_pibaq(res: DataFrame) -> DataFrame:
    """Normalize piBAQ values by the per-sample piBAQ total.

    The resulting relative piBAQ values are then multiplied by 100'000'000
    (PRIDE database normalization) for the piBAQ ppb and log10 shifted
    by 10 (ProteomicsDB).

    Parameters
    ----------
    res : pd.DataFrame
        Input DataFrame with a :data:`PIBAQ` column.

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized piBAQ values.
    """
    # Relative piBAQ: divide each protein's value by its group total.
    # ``groupby(...).transform`` is the vectorised, future-proof equivalent of
    # the per-group ``apply(normalize)`` it replaces (no DeprecationWarning,
    # no implicit grouping-column handling).
    group_totals = res.groupby([SAMPLE_ID, CONDITION], observed=True)[PIBAQ].transform(
        "sum"
    )
    res[PIBAQ_NORMALIZED] = res[PIBAQ] / group_totals
    # Normalization method used by ProteomicsDB: 10 + log10(piBAQ / total piBAQ).
    res[PIBAQ_LOG] = np.where(
        res[PIBAQ_NORMALIZED] > 0, np.log10(res[PIBAQ_NORMALIZED]) + 10, 0
    )
    # PRIDE normalization (no log transform): piBAQ / total piBAQ * 100'000'000.
    res[PIBAQ_PPB] = res[PIBAQ_NORMALIZED] * 100_000_000
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


_EXTRACT_FASTA_ARGUMENTS = (
    "fasta",
    "enzyme",
    "proteins",
    "min_aa",
    "max_aa",
    "tpa",
)


@log_function_call(logger)
def extract_fasta(*args: object, **kwargs: object):
    """Forward to :func:`mokume.io.fasta.extract_fasta`.

    Kept as a thin wrapper so existing callers can keep importing
    ``extract_fasta`` from :mod:`mokume.quantification.pibaq`. All the
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
        ``features2proteins --quant-method pibaq``), which instead pairs a
        *shared-aware numerator* (proportional shared-peptide reallocation)
        with a *total-potential* denominator (proteotypic + shared) via
        :func:`mokume.io.fasta.digest_fasta_full`. Do not mix the two
        conventions: the proteotypic-only and total-potential denominators
        are each only correct alongside their matching numerator. Retained
        for direct API consumers; prefer ``compute_pibaq`` for new code.
    """
    bound = _bind_arguments("extract_fasta", args, kwargs, _EXTRACT_FASTA_ARGUMENTS, {})
    return _extract_fasta_io(**bound)


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
        """Return the summed intensity of accession- or entry-matched histones."""
        # ``ProteinName`` carries UniProt accessions (from ``pg_accessions`` in the
        # quantms feature files), so match against ``histone_proteins`` (accessions)
        # in addition to ``histone_entries`` (entry names). Matching only entry
        # names found no histones in an accession-keyed table, so the intensity fell
        # back to 1.0 and inflated CopyNumber. The two id forms never collide, so the
        # union is safe whichever form ``ProteinName`` uses. See issue #19.
        histones = set(self.organism.histone_proteins) | set(
            self.organism.histone_entries
        )
        is_histone_mask = protein_intensities[PROTEIN_NAME].isin(histones)
        histone_intensities = max(
            protein_intensities[is_histone_mask][NORM_INTENSITY].sum(), 1.0
        )
        return histone_intensities

    def apply_ruler(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        """Calculate copy number, amount, weight, and concentration columns."""
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
        """Apply the proteomic ruler to one condition's protein rows."""
        return self.apply_ruler(protein_intensities)

    def apply_by_condition(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
        """Apply the proteomic ruler independently within each condition."""
        # Apply the ruler independently per condition. pandas 3.0 excludes the
        # grouping column from ``apply``'s result, silently dropping
        # ``Condition`` from the output table; older pandas kept it. We keep the
        # original per-condition computation unchanged and only realign the
        # label when it is missing. ``group_keys=False`` preserves the original
        # row index, so the assignment maps each row to its own condition; we
        # then restore the original column order (pandas 3.0 appends the
        # re-added column) so the output schema matches across pandas versions.
        result = protein_intensities.groupby([CONDITION], group_keys=False).apply(self)
        if CONDITION not in result.columns:
            result[CONDITION] = protein_intensities[CONDITION]
            added = [c for c in result.columns if c not in protein_intensities.columns]
            result = result[list(protein_intensities.columns) + added]
        return result


# ---------------------------------------------------------------------------
# piBAQ (paralog-aware iBAQ) -- exact gpGrouper-style shared allocation
# ---------------------------------------------------------------------------


_COMPUTE_PIBAQ_ARGUMENTS = (
    "peptide_df",
    "accession_to_peptides",
    "peptide_to_accessions",
    "families",
)
_COMPUTE_PIBAQ_DEFAULTS = {
    "mw_map": None,
    "min_anchors": 1,
    "high_anchor_threshold": 3,
    "extra_group_cols": None,
}


@dataclass(frozen=True)
class _ComputePibaqOptions:
    mw_map: Optional[Mapping[str, float]]
    min_anchors: int
    high_anchor_threshold: int
    extra_group_cols: Optional[Sequence[str]]


@dataclass(frozen=True)
class _ComputePibaqRequest:
    peptide_df: DataFrame
    accession_to_peptides: Mapping[str, Set[str]]
    peptide_to_accessions: Mapping[str, Set[str]]
    families: Sequence[Family]
    options: _ComputePibaqOptions


@dataclass(frozen=True)
class _PreparedPibaq:
    peptide_column: str
    group_columns: List[str]
    working: DataFrame
    anchor_counts: Mapping[str, int]
    family_to_peptides: Mapping[str, Set[str]]
    observed_peptides: Set[str]


def compute_pibaq(*args: object, **kwargs: object) -> DataFrame:
    """Compute per-protein piBAQ with exact shared-peptide allocation.

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
        Unique-anchor threshold. If no family member reaches it, shared signal
        is split equally and evidence is ``family_only``. Defaults to 1.
    high_anchor_threshold : int, optional
        Threshold (in unique anchors of the *weakest* family member) at
        which ``EvidenceLevel`` is reported as ``high``. Defaults to 3.
    extra_group_cols : Sequence[str], optional
        Additional dataframe columns to retain in the groupby keys, on
        top of (SampleID, Condition).

    Returns
    -------
    DataFrame
        Long-format result containing positive per-member estimates for
        ``(member, SampleID, Condition, *extra_group_cols)``. Includes the
        :data:`FAMILY_ID`, :data:`EVIDENCE_LEVEL`, :data:`FAMILY_SIZE`
        metadata columns, and :data:`MOLECULARWEIGHT` + :data:`TPA` when
        ``mw_map`` is supplied.
    """
    bound = _bind_arguments(
        "compute_pibaq",
        args,
        kwargs,
        _COMPUTE_PIBAQ_ARGUMENTS,
        _COMPUTE_PIBAQ_DEFAULTS,
    )
    options = _ComputePibaqOptions(
        mw_map=bound.pop("mw_map"),
        min_anchors=bound.pop("min_anchors"),
        high_anchor_threshold=bound.pop("high_anchor_threshold"),
        extra_group_cols=bound.pop("extra_group_cols"),
    )
    return _run_compute_pibaq(_ComputePibaqRequest(options=options, **bound))


setattr(
    compute_pibaq,
    "__text_signature__",
    "(peptide_df, accession_to_peptides, peptide_to_accessions, families, *, "
    "mw_map=None, min_anchors=1, high_anchor_threshold=3, extra_group_cols=None)",
)


def _empty_compute_result(
    request: _ComputePibaqRequest, group_columns: Sequence[str]
) -> DataFrame:
    """Return the empty output schema for a computation request."""
    return _empty_pibaq_frame(group_columns, request.options.mw_map is not None)


def _prepare_pibaq(
    request: _ComputePibaqRequest, group_columns: List[str]
) -> Optional[_PreparedPibaq]:
    """Filter observed peptides and build family ownership indices."""
    peptide_column = _detect_peptide_column(request.peptide_df)
    working = request.peptide_df[
        _membership_mask(
            request.peptide_df[peptide_column], request.peptide_to_accessions
        )
    ]
    if working.empty:
        return None
    anchor_counts = _count_unique_anchors(
        working, request.peptide_to_accessions, peptide_column
    )
    peptide_owner = _assign_peptides_to_owning_family(
        request.families, request.peptide_to_accessions, anchor_counts
    )
    family_to_peptides = _invert_peptide_ownership(peptide_owner)
    working = working[_membership_mask(working[peptide_column], peptide_owner)]
    return _PreparedPibaq(
        peptide_column=peptide_column,
        group_columns=group_columns,
        working=working,
        anchor_counts=anchor_counts,
        family_to_peptides=family_to_peptides,
        observed_peptides=set(working[peptide_column].dropna().unique()),
    )


def _family_evidence(
    request: _ComputePibaqRequest,
    prepared: _PreparedPibaq,
    family: Family,
) -> str:
    """Classify one family from its observed unique-anchor counts."""
    anchors = [prepared.anchor_counts.get(member, 0) for member in family.members]
    return _classify_evidence(
        min(anchors),
        max(anchors),
        request.options.min_anchors,
        request.options.high_anchor_threshold,
    )


def _quantify_prepared_family(
    request: _ComputePibaqRequest,
    prepared: _PreparedPibaq,
    family: Family,
) -> Optional[Tuple[DataFrame, str]]:
    """Compute and annotate one observed family's piBAQ block."""
    owned = prepared.family_to_peptides.get(family.family_id, set())
    if not owned & prepared.observed_peptides:
        return None
    family_input = _assemble_family_input(
        prepared.working,
        family,
        prepared.family_to_peptides,
        prepared.peptide_column,
    )
    if family_input.empty:
        return None
    evidence = _family_evidence(request, prepared, family)
    block = _allocate_family(
        _FamilyAllocationInputs(
            family=family,
            peptide_to_accessions=request.peptide_to_accessions,
            accession_to_peptides=request.accession_to_peptides,
            owned_peptides=owned,
            force_equal_shared=evidence == EVIDENCE_FAMILY_ONLY,
        ),
        family_input,
        prepared.peptide_column,
        prepared.group_columns,
    )
    if block.empty:
        return None
    return _annotate_family_metadata(block, family, evidence), evidence


def _collect_family_results(
    request: _ComputePibaqRequest, prepared: _PreparedPibaq
) -> Tuple[List[DataFrame], int]:
    """Collect non-empty family blocks and count family-only results."""
    results: List[DataFrame] = []
    family_only = 0
    for family in request.families:
        quantified = _quantify_prepared_family(request, prepared, family)
        if quantified is None:
            continue
        block, evidence = quantified
        results.append(block)
        family_only += evidence == EVIDENCE_FAMILY_ONLY
    return results, family_only


def _run_compute_pibaq(request: _ComputePibaqRequest) -> DataFrame:
    """Execute the unchanged piBAQ calculation for a bound request."""
    group_columns = _resolve_group_cols(
        request.peptide_df, request.options.extra_group_cols
    )
    if request.peptide_df.empty or not request.families:
        return _empty_compute_result(request, group_columns)
    prepared = _prepare_pibaq(request, group_columns)
    if prepared is None:
        return _empty_compute_result(request, group_columns)
    results, family_only = _collect_family_results(request, prepared)
    if not results:
        return _empty_compute_result(request, group_columns)
    output = pd.concat(results, ignore_index=True)
    if request.options.mw_map is not None:
        output = _finalize_tpa(output, request.options.mw_map)
    processed = len(results)
    logger.info(
        "piBAQ: %d families processed (%d family-only, %d member-resolving)",
        processed,
        family_only,
        processed - family_only,
    )
    return output


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
