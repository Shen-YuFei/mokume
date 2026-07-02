"""
Enumeration types for preprocessing filters.
"""

from enum import Enum, auto


class RazorPeptideHandling(Enum):
    """How to handle razor (shared) peptides."""

    KEEP = auto()  # Keep all razor peptides
    REMOVE = auto()  # Remove all razor peptides
    ASSIGN_TO_TOP = auto()  # Assign to protein with most peptides

    @classmethod
    def from_str(cls, name: str) -> "RazorPeptideHandling":
        """Convert string to enum value (case-insensitive)."""
        name_ = name.lower().replace("-", "_").replace(" ", "_")
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(f"Unknown razor peptide handling: {name}")


class ProteinGroupingStrategy(Enum):
    """Protein grouping strategies."""

    NONE = auto()  # No grouping
    SUBSUMPTION = auto()  # Group proteins whose peptides are subsets
    PARSIMONY = auto()  # Minimal set of proteins explaining peptides

    @classmethod
    def from_str(cls, name: str) -> "ProteinGroupingStrategy":
        """Convert string to enum value (case-insensitive)."""
        name_ = name.lower().replace("-", "_").replace(" ", "_")
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(f"Unknown protein grouping strategy: {name}")


class FilterLevel(Enum):
    """Levels at which filtering can be applied."""

    FEATURE = auto()
    PEPTIDE = auto()
    PROTEIN = auto()
    RUN = auto()
    SAMPLE = auto()

    @classmethod
    def from_str(cls, name: str) -> "FilterLevel":
        """Convert string to enum value (case-insensitive)."""
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(f"Unknown filter level: {name}")


class FilterAction(Enum):
    """Actions taken when a filter condition is met."""

    REMOVE = auto()  # Remove matching items
    KEEP = auto()  # Keep only matching items
    FLAG = auto()  # Flag but don't remove

    @classmethod
    def from_str(cls, name: str) -> "FilterAction":
        """Convert string to enum value (case-insensitive)."""
        name_ = name.lower()
        for k, v in cls._member_map_.items():
            if k.lower() == name_:
                return v
        raise KeyError(f"Unknown filter action: {name}")
