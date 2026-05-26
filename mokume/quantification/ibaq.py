"""
iBAQ protein quantification method.

This module provides the iBAQ (intensity-Based Absolute Quantification)
method for protein quantification, including normalization and proteomic
ruler calculations.
"""

import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from pandas import DataFrame, Series

from mokume.core.constants import (
    CONCENTRATION_NM,
    CONDITION,
    COPYNUMBER,
    IBAQ,
    IBAQ_LOG,
    IBAQ_NORMALIZED,
    IBAQ_PPB,
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
from mokume.io.fasta import extract_fasta as _extract_fasta_io
from mokume.model.organism import OrganismDescription
from mokume.plotting import is_plotting_available

# Proteomic Ruler constants
AVAGADRO: float = 6.02214129e23
AVERAGE_BASE_PAIR_MASS: float = 617.96

# Get a logger for this module
logger = get_logger("mokume.quantification.ibaq")


@log_function_call(logger, level=logging.DEBUG)
def normalize(group):
    """
    Normalize the ibaq values using the total ibaq of the sample.

    This method is called rIBAQ, originally published in
    https://pubs.acs.org/doi/10.1021/pr401017h

    Parameters
    ----------
    group : pd.DataFrame
        Dataframe with all the ibaq values.

    Returns
    -------
    pd.DataFrame
        Dataframe with the normalized ibaq values.
    """
    group[IBAQ_NORMALIZED] = group[IBAQ] / group[IBAQ].sum()
    return group


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
    res = res.groupby([SAMPLE_ID, CONDITION]).apply(normalize)
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
    peptide counting that this PR introduces — lives in a single
    implementation in :mod:`mokume.io.fasta` to avoid drift between the two
    public entry points. See that function's docstring for the full contract.
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


class PeptideProteinMapper:
    """Map peptides to proteins and calculate iBAQ values."""

    _peptide_protein_ratio: dict[str, float]
    unique_peptide_counts: dict[str, int]
    map_size: dict[str, int]
    protein_mass_map: dict[str, float]

    def __init__(
        self,
        unique_peptide_counts: Optional[dict[str, int]] = None,
        map_size: Optional[dict[str, int]] = None,
        protein_mass_map: Optional[dict[str, float]] = None,
    ):
        self.unique_peptide_counts = unique_peptide_counts or {}
        self.map_size = map_size or {}
        self.protein_mass_map = protein_mass_map or {}
        self._peptide_protein_ratio = {}

    def peptide_protein_ratio(self, protein_group: str):
        if protein_group in self._peptide_protein_ratio:
            return self._peptide_protein_ratio[protein_group]

        proteins_list = protein_group.split(";")
        total = 0
        for prot in proteins_list:
            total += self.unique_peptide_counts[prot]

        if not proteins_list:
            val = self._peptide_protein_ratio[protein_group] = 0
        else:
            val = self._peptide_protein_ratio[protein_group] = total / len(
                proteins_list
            )
        return val

    def get_average_nr_peptides_unique_by_group(
        self, pdrow: Series
    ) -> Union[float, Series]:
        """Calculate the average number of unique peptides per protein group."""
        average_peptides_per_protein = self.peptide_protein_ratio(pdrow.name[0])

        if average_peptides_per_protein > 0:
            return (
                pdrow.NormIntensity
                / self.map_size[pdrow.name]
                / average_peptides_per_protein
            )

        return np.nan

    def protein_group_mass(self, protein_group: str):
        """Calculate the molecular weight of a protein group."""
        mw_list = [self.protein_mass_map[i] for i in protein_group.split(";")]
        return sum(mw_list)


def _keep_only_proteotypic_rows(data: pd.DataFrame) -> pd.DataFrame:
    """Drop shared-peptide rows so that both the iBAQ and TPA numerators
    aggregate only proteotypic intensity.

    The mokume pipeline (``features2proteins``) already strips ``unique != 1``
    rows upstream and drops the ``unique`` column, so this filter is a no-op
    on that path. When callers feed a raw QPX feature parquet directly
    (e.g. via the ``peptides2protein`` CLI), this prevents shared-homologue
    signal from inflating the numerator of large homologous families
    (myosin, tubulin, histone, ...) — the same signal that would otherwise
    double-count once it appears under every protein the peptide maps to.
    Together with the cross-protein unique denominator in
    :func:`extract_fasta` this keeps the iBAQ and TPA ratios symmetric.
    """
    if "unique" not in data.columns:
        return data
    n_before = len(data)
    data = data[data["unique"].isin([1, "1", True])]
    n_dropped = n_before - len(data)
    if n_dropped:
        logger.info(
            "Dropped %d/%d shared-peptide rows (unique != 1) before "
            "iBAQ/TPA aggregation; numerators will sum proteotypic "
            "intensity only.",
            n_dropped,
            n_before,
        )
    return data


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

    data = _keep_only_proteotypic_rows(data)

    # get fasta info
    proteins = data[PROTEIN_NAME].unique().tolist()
    proteins = sum([i.split(";") for i in proteins], [])

    unique_peptide_counts, mw_dict, found_proteins = extract_fasta(
        fasta, enzyme, proteins, min_aa, max_aa, tpa
    )

    data = data[data[PROTEIN_NAME].isin(found_proteins)]

    # data processing
    logger.info("Processing data with %d rows", len(data))
    logger.debug("Data sample: \n%s", data.head().to_string())
    map_size = data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION]).size().to_dict()
    res = pd.DataFrame(
        data.groupby([PROTEIN_NAME, SAMPLE_ID, CONDITION])[NORM_INTENSITY].sum()
    )

    protein_mapper = PeptideProteinMapper(
        unique_peptide_counts=unique_peptide_counts,
        map_size=map_size,
        protein_mass_map=mw_dict,
    )

    # ibaq
    res[IBAQ] = res.apply(protein_mapper.get_average_nr_peptides_unique_by_group, 1)
    res = res.reset_index()

    # normalize ibaq
    if normalize:
        res = normalize_ibaq(res)
        res = res.dropna(subset=[IBAQ_NORMALIZED])
        plot_column = IBAQ_PPB
    else:
        res = res.dropna(subset=[IBAQ])
        plot_column = IBAQ

    res = res.reset_index(drop=True)

    # tpa
    if tpa:
        res[MOLECULARWEIGHT] = (
            res[PROTEIN_NAME]
            .apply(protein_mapper.protein_group_mass)
            .fillna(1.0)
            .replace(0.0, 1.0)
        )
        res[TPA] = res[NORM_INTENSITY] / res[MOLECULARWEIGHT]

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
