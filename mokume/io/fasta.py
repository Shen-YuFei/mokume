"""
FASTA file handling utilities.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from pyopenms import AASequence, FASTAFile, ProteaseDigestion

from mokume.core.constants import build_accession_map, get_accession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def load_fasta(fasta_path: str) -> List:
    """
    Load a FASTA file and return the list of protein entries.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    List
        List of FASTA protein entries.
    """
    fasta_proteins = []
    FASTAFile().load(fasta_path, fasta_proteins)
    return fasta_proteins


def digest_protein(
    sequence: str,
    enzyme: str,
    min_aa: int = 7,
    max_aa: int = 30,
) -> List[str]:
    """
    Digest a protein sequence using a specified enzyme.

    Parameters
    ----------
    sequence : str
        Protein amino acid sequence.
    enzyme : str
        Name of the enzyme to use for digestion.
    min_aa : int
        Minimum peptide length to include.
    max_aa : int
        Maximum peptide length to include.

    Returns
    -------
    List[str]
        List of peptide sequences.
    """
    digestor = ProteaseDigestion()
    digestor.setEnzyme(enzyme)

    digest = []
    digestor.digest(AASequence.fromString(sequence), digest, min_aa, max_aa)
    return [str(pep.toString()) for pep in digest]


_NONSTANDARD_AA = {"X", "B", "Z", "J", "U", "O"}


def _strip_nonstandard_aa(sequence: str) -> str:
    """Remove non-standard amino acids from a protein sequence."""
    for aa in _NONSTANDARD_AA:
        sequence = sequence.replace(aa, "")
    return sequence


def extract_fasta(
    fasta: str,
    enzyme: str,
    proteins: List[str],
    min_aa: int,
    max_aa: int,
    tpa: bool = False,
) -> Tuple[Dict[str, int], Dict[str, float], Set[str]]:
    """
    Extract protein information from a FASTA file using a specified enzyme for digestion.

    The number of enzyme-specific theoretical peptides reported per protein
    is the count of peptides that are *unique across the full FASTA database*
    (proteotypic). Counting all theoretical peptides — including those shared
    with homologous proteins — would systematically deflate iBAQ for large
    homologous families (myosin, tubulin, actin, histone, keratin) by 3-20x
    because the denominator would be inflated by shared peptides while the
    numerator only contains proteotypic intensity.

    Sequences with nonstandard amino acids (X/B/Z/J/U/O) are stripped before
    digestion. Callers asking for proteins that are absent from the FASTA are
    silently omitted from ``uniquepepcounts`` and ``mw_dict``; only when
    *no* requested protein is found does the function raise ``ValueError``.

    Parameters
    ----------
    fasta : str
        Path to the FASTA file containing protein sequences.
    enzyme : str
        Name of the enzyme used for protein digestion.
    proteins : List[str]
        List of protein accessions to search for in the FASTA file.
    min_aa : int
        Minimum number of amino acids for peptides to be considered.
    max_aa : int
        Maximum number of amino acids for peptides to be considered.
    tpa : bool
        If True, calculate molecular weights for Total Protein Approach (TPA).

    Returns
    -------
    Tuple[Dict[str, int], Dict[str, float], Set[str]]
        A tuple containing:
        - uniquepepcounts: Dictionary mapping caller-provided protein names to
          cross-protein unique theoretical peptide counts.
        - mw_dict: Dictionary mapping caller-provided protein names to
          molecular weights (empty if ``tpa=False``).
        - found_proteins: Set of caller-provided protein names (i.e. the
          original entries from ``proteins``) that resolved to a FASTA entry.

    Raises
    ------
    ValueError
        If none of the specified proteins are found in the FASTA file.
    """
    acc_to_originals, protein_accessions = build_accession_map(proteins)

    fasta_proteins = load_fasta(fasta)
    digestor = ProteaseDigestion()
    digestor.setEnzyme(enzyme)

    accession_to_peptides, peptide_to_accessions, accession_to_mw = _digest_full_fasta(
        fasta_proteins,
        digestor,
        min_aa,
        max_aa,
        tpa,
        keep_peptides_for=protein_accessions,
    )
    uniquepepcounts, mw_dict, found_proteins = _select_unique_counts_for_callers(
        protein_accessions,
        acc_to_originals,
        accession_to_peptides,
        peptide_to_accessions,
        accession_to_mw,
        tpa,
    )

    if len(found_proteins) == 0:
        raise ValueError(
            "None of the specified proteins were found in the FASTA file. "
            "Please check that the protein accessions match."
        )

    logger.info("Found %d proteins in FASTA file", len(found_proteins))
    return uniquepepcounts, mw_dict, found_proteins


def _digest_one(sequence, digestor: ProteaseDigestion, min_aa: int, max_aa: int):
    """Digest a single sequence and return its peptide set.

    Uses ``pep.toString()`` rather than ``str(pep)`` to keep peptide string
    keys stable across pyOpenMS versions (``__str__`` has varied; ``toString``
    is part of the documented API).
    """
    digest: List = []
    digestor.digest(sequence, digest, min_aa, max_aa)
    return {pep.toString() for pep in digest}


def _digest_full_fasta(
    fasta_proteins,
    digestor: ProteaseDigestion,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    keep_peptides_for: Set[str],
):
    """Digest every protein in the FASTA once and return the inverted index
    needed to compute cross-protein peptide uniqueness.

    The inverted index ``peptide_to_accessions`` is populated for every FASTA
    entry (cross-protein uniqueness requires it). To avoid wasting memory on
    large proteomes, ``accession_to_peptides`` and ``accession_to_mw`` are
    only retained for accessions in ``keep_peptides_for``.

    Per-entry parse, digest, and weight failures are logged and skipped so
    that one bad sequence cannot poison quantification of the rest of the
    proteome.
    """
    accession_to_peptides: Dict[str, Set[str]] = {}
    peptide_to_accessions: "defaultdict[str, Set[str]]" = defaultdict(set)
    accession_to_mw: Dict[str, float] = {}
    for entry in fasta_proteins:
        accession = get_accession(entry.identifier)
        try:
            aa_sequence = AASequence.fromString(_strip_nonstandard_aa(entry.sequence))
            peps = _digest_one(aa_sequence, digestor, min_aa, max_aa)
        except (ValueError, RuntimeError) as exc:
            logger.warning("Skipping %s: %s", accession, exc)
            continue
        for pep in peps:
            peptide_to_accessions[pep].add(accession)
        if accession in keep_peptides_for:
            accession_to_peptides[accession] = peps
            if tpa:
                try:
                    accession_to_mw[accession] = aa_sequence.getMonoWeight()
                except (ValueError, RuntimeError) as exc:
                    logger.warning("Skipping MW for %s: %s", accession, exc)
    return accession_to_peptides, peptide_to_accessions, accession_to_mw


def _emit_mw_for_caller(
    originals: List[str],
    accession: str,
    accession_to_mw: Dict[str, float],
    mw_dict: Dict[str, float],
) -> None:
    """Copy the molecular weight of ``accession`` to each caller name when known."""
    if accession not in accession_to_mw:
        return
    mw = accession_to_mw[accession]
    for orig in originals:
        mw_dict[orig] = mw


def _select_unique_counts_for_callers(
    protein_accessions: Set[str],
    acc_to_originals: Dict[str, List[str]],
    accession_to_peptides: Dict[str, Set[str]],
    peptide_to_accessions: Dict[str, Set[str]],
    accession_to_mw: Dict[str, float],
    tpa: bool,
):
    """Emit per-caller cross-protein unique peptide counts (and MW)."""
    found_proteins: Set[str] = set()
    uniquepepcounts: Dict[str, int] = {}
    mw_dict: Dict[str, float] = {}
    for accession in protein_accessions:
        peps = accession_to_peptides.get(accession)
        if peps is None:
            continue
        originals = acc_to_originals[accession]
        found_proteins.update(originals)
        unique_count = sum(
            1 for pep in peps if peptide_to_accessions[pep] == {accession}
        )
        for orig in originals:
            uniquepepcounts[orig] = unique_count
        if tpa:
            _emit_mw_for_caller(originals, accession, accession_to_mw, mw_dict)
    return uniquepepcounts, mw_dict, found_proteins


def get_protein_molecular_weights(
    fasta: str,
    proteins: List[str],
) -> Dict[str, float]:
    """
    Calculate molecular weights for a list of proteins from a FASTA file.

    Parameters
    ----------
    fasta : str
        Path to the FASTA file.
    proteins : List[str]
        List of protein accessions.

    Returns
    -------
    Dict[str, float]
        Dictionary mapping protein accessions to molecular weights.
    """
    acc_to_originals, protein_accessions = build_accession_map(proteins)

    fasta_proteins = load_fasta(fasta)
    mw_dict: Dict[str, float] = {}

    for entry in fasta_proteins:
        accession = get_accession(entry.identifier)
        if accession in protein_accessions:
            aa_sequence = AASequence.fromString(_strip_nonstandard_aa(entry.sequence))
            mw = aa_sequence.getMonoWeight()
            for orig in acc_to_originals[accession]:
                mw_dict[orig] = mw

    return mw_dict
