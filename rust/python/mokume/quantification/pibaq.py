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

import logging
from pathlib import Path
from typing import (
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
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
from mokume.quantification.families import (
    Family,
    discover_families,
    load_families_yaml,
    merge_overrides,
)

# Proteomic Ruler constants
AVAGADRO: float = 6.02214129e23
AVERAGE_BASE_PAIR_MASS: float = 617.96

# Get a logger for this module
logger = get_logger("mokume.quantification.pibaq")


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


@log_function_call(logger)
def extract_fasta(
    fasta: str, enzyme: str, proteins: List, min_aa: int, max_aa: int, tpa: bool
):
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

    def apply_by_condition(self, protein_intensities: pd.DataFrame) -> pd.DataFrame:
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
    if peptide_df.empty or not families:
        include_tpa = mw_map is not None
        group_cols = _resolve_group_cols(peptide_df, extra_group_cols)
        return _empty_pibaq_frame(group_cols, include_tpa)

    pep_col = _detect_peptide_column(peptide_df)
    group_cols = _resolve_group_cols(peptide_df, extra_group_cols)

    # Keep only peptides present in the FASTA digest; matching is driven by
    # the peptide-sequence column against the digest index, not by the
    # per-row protein-name strings.
    working = peptide_df[_membership_mask(peptide_df[pep_col], peptide_to_accessions)]
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
    working = working[_membership_mask(working[pep_col], peptide_owner)]
    observed_peptides = set(working[pep_col].dropna().unique())

    results: List[DataFrame] = []
    n_family_only = 0
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
        evidence = _classify_evidence(
            min_anchor, max_anchor, min_anchors, high_anchor_threshold
        )
        block = _allocate_family(
            _FamilyAllocationInputs(
                family=family,
                peptide_to_accessions=peptide_to_accessions,
                accession_to_peptides=accession_to_peptides,
                owned_peptides=owned,
                force_equal_shared=evidence == EVIDENCE_FAMILY_ONLY,
            ),
            family_input,
            pep_col,
            group_cols,
        )
        if block.empty:
            continue
        if evidence == EVIDENCE_FAMILY_ONLY:
            n_family_only += 1
        block = _annotate_family_metadata(block, family, evidence)
        results.append(block)

    if not results:
        include_tpa = mw_map is not None
        return _empty_pibaq_frame(group_cols, include_tpa)

    out = pd.concat(results, ignore_index=True)
    if mw_map is not None:
        out = _finalize_tpa(out, mw_map)
    n_processed = len(results)
    logger.info(
        "piBAQ: %d families processed (%d family-only, %d member-resolving)",
        n_processed,
        n_family_only,
        n_processed - n_family_only,
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
    ``unique != 1`` rows: piBAQ allocates shared-peptide intensity from the
    FASTA digest. Razor-mirror rows that list the same peptide observation
    under multiple proteins are de-duplicated inside the family allocator via
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
    Compute piBAQ values for peptides and generate a QC report.

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
        Unique-anchor threshold. If no family member reaches it, shared signal
        is split equally and evidence is ``family_only``. Defaults to 1.
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

    # Normalize piBAQ.
    if normalize:
        res = normalize_pibaq(res)
        res = res.dropna(subset=[PIBAQ_NORMALIZED])
        plot_column = PIBAQ_PPB
    else:
        res = res.dropna(subset=[PIBAQ])
        plot_column = PIBAQ

    res = res.reset_index(drop=True)

    # calculate protein weight and concentration
    if ruler:
        concentration_by_ruler = ConcentrationWeightByProteomicRuler(
            organism_descr, ploidy, cpc
        )
        res = concentration_by_ruler.apply_by_condition(res)

    # Plot the distribution of the protein piBAQ values.
    if verbose:
        if not is_plotting_available():
            logger.warning(
                "QC report skipped: plotting dependencies not installed. "
                "Install: pip install mokume-py[plotting] (pure Python) or "
                "pip install mokume[plotting] (Rust wheel)"
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
