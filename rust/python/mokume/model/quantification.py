"""
Protein quantification method enumerations for the mokume package.

This module provides an enumeration of protein quantification methods
supported by the package.
"""

from enum import Enum, auto


class QuantificationMethod(Enum):
    """
    Enumeration of protein quantification methods.

    Attributes
    ----------
    IBAQ : auto
        Intensity-Based Absolute Quantification.
        Divides the sum of peptide intensities by the number of
        theoretically observable peptides.
    TOP3 : auto
        Top3 quantification method.
        Uses the average of the three most intense peptides.
    TOPN : auto
        TopN quantification method.
        Uses the average of the N most intense peptides.
    ALL_PEPTIDES : auto
        Sum of all peptide intensities.
    MAXLFQ : auto
        MaxLFQ algorithm for label-free quantification.
        Uses delayed normalization and maximal peptide ratio extraction.
    TMT_ABUNDANCE : auto
        TMT abundance method (``abd``).
        Median of log2-transformed peptide intensities per protein.
    TMT_REPORTER_INTENSITY : auto
        TMT reporter intensity method (``intensity``).
        Sum of raw reporter intensities per protein (linear space).
    RATIO : auto
        TMT ratio quantification (Proteome Sciences protocol).
        log2(sample / reference) per PSM, median-aggregated to protein.
    """

    IBAQ = auto()
    TOP3 = auto()
    TOPN = auto()
    ALL_PEPTIDES = auto()
    MAXLFQ = auto()
    TMT_ABUNDANCE = auto()
    TMT_REPORTER_INTENSITY = auto()
    RATIO = auto()

    @classmethod
    def from_str(cls, name: str) -> "QuantificationMethod":
        """
        Convert a string to a QuantificationMethod.

        Parameters
        ----------
        name : str
            The name of the quantification method.

        Returns
        -------
        QuantificationMethod
            The quantification method.

        Raises
        ------
        KeyError
            If the name does not match any quantification method.
        """
        name_ = name.lower().replace("-", "_").replace(" ", "_")
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(name)

    @property
    def description(self) -> str:
        """
        Get a human-readable description of the quantification method.

        Returns
        -------
        str
            Description of the quantification method.
        """
        descriptions = {
            QuantificationMethod.IBAQ: "Intensity-Based Absolute Quantification",
            QuantificationMethod.TOP3: "Average of top 3 most intense peptides",
            QuantificationMethod.TOPN: "Average of top N most intense peptides",
            QuantificationMethod.ALL_PEPTIDES: "Sum of all peptide intensities",
            QuantificationMethod.MAXLFQ: "MaxLFQ label-free quantification",
            QuantificationMethod.TMT_ABUNDANCE: (
                "TMT abundance: median of log2 peptide intensities"
            ),
            QuantificationMethod.TMT_REPORTER_INTENSITY: (
                "TMT reporter intensity: sum of raw reporter intensities"
            ),
            QuantificationMethod.RATIO: (
                "TMT ratio: log2(sample/reference) per PSM, median-aggregated"
            ),
        }
        return descriptions.get(self, "Unknown method")
