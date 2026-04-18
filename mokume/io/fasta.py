"""
FASTA file handling utilities.
"""

import logging
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

    This function processes a FASTA file to extract proteins, performs in-silico digestion
    using a specified enzyme, calculates unique peptide counts, and optionally computes
    molecular weights. It handles sequences with nonstandard amino acids by removing them
    and raises an error if none of the specified proteins are found in the FASTA file.

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
        - uniquepepcounts: Dictionary mapping protein accessions to unique peptide counts.
        - mw_dict: Dictionary mapping protein accessions to molecular weights (empty if tpa=False).
        - found_proteins: Set of protein accessions that were found in the FASTA file.

    Raises
    ------
    ValueError
        If none of the specified proteins are found in the FASTA file.
    """
    acc_to_originals, protein_accessions = build_accession_map(proteins)

    fasta_proteins = load_fasta(fasta)
    found_proteins: Set[str] = set()
    uniquepepcounts: Dict[str, int] = {}
    mw_dict: Dict[str, float] = {}

    digestor = ProteaseDigestion()
    digestor.setEnzyme(enzyme)

    for entry in fasta_proteins:
        accession = get_accession(entry.identifier)
        if accession in protein_accessions:
            originals = acc_to_originals[accession]
            found_proteins.update(originals)
            aa_sequence = AASequence.fromString(
                _strip_nonstandard_aa(entry.sequence)
            )

            # Perform digestion
            digest = []
            digestor.digest(aa_sequence, digest, min_aa, max_aa)

            unique_digest = set(digest)
            for orig in originals:
                uniquepepcounts[orig] = len(unique_digest)

            # Calculate molecular weight if needed for TPA
            if tpa:
                mw = aa_sequence.getMonoWeight()
                for orig in originals:
                    mw_dict[orig] = mw

    if len(found_proteins) == 0:
        raise ValueError(
            "None of the specified proteins were found in the FASTA file. "
            "Please check that the protein accessions match."
        )

    logger.info(f"Found {len(found_proteins)} proteins in FASTA file")
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
            aa_sequence = AASequence.fromString(
                _strip_nonstandard_aa(entry.sequence)
            )
            mw = aa_sequence.getMonoWeight()
            for orig in acc_to_originals[accession]:
                mw_dict[orig] = mw

    return mw_dict
