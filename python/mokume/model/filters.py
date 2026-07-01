"""
Filter configuration models for preprocessing.

This module provides dataclasses for configuring preprocessing filters
at the intensity, peptide, protein, and sample QC levels.
"""

from dataclasses import dataclass, field, asdict
from typing import ClassVar, Optional, List


@dataclass
class IntensityFilterConfig:
    """
    Configuration for intensity-based filters.

    Attributes
    ----------
    min_intensity : float
        Minimum intensity threshold. Features below this are removed.
    cv_threshold : float, optional
        Maximum coefficient of variation (CV) threshold for technical replicates.
    min_replicate_agreement : int
        Minimum number of replicates where a feature must be detected.
    quantile_lower : float
        Lower quantile threshold for filtering (0.0-1.0).
    quantile_upper : float
        Upper quantile threshold for filtering (0.0-1.0).
    remove_zero_intensity : bool
        Whether to remove zero intensity values.
    """

    registry: ClassVar[dict[str, "IntensityFilterConfig"]] = {}

    name: str = "default"
    min_intensity: float = 0.0
    cv_threshold: Optional[float] = None
    min_replicate_agreement: int = 1
    quantile_lower: float = 0.0
    quantile_upper: float = 1.0
    remove_zero_intensity: bool = True

    @classmethod
    def get(cls, name: str, default=None) -> "Optional[IntensityFilterConfig]":
        """Retrieve a configuration from the registry."""
        return cls.registry.get(name.lower(), default)

    def __post_init__(self):
        """Register this configuration in the class registry."""
        if self.name:
            self.registry[self.name.lower()] = self

    def to_dict(self) -> dict:
        """Convert to dictionary, excluding registry."""
        d = asdict(self)
        return d


@dataclass
class PeptideFilterConfig:
    """
    Configuration for peptide-level filters.

    Attributes
    ----------
    min_search_score : float, optional
        Minimum search engine score threshold.
    allowed_charge_states : list[int], optional
        List of allowed charge states (e.g., [2, 3, 4]).
    exclude_modifications : list[str]
        List of modification names to exclude.
    max_missed_cleavages : int, optional
        Maximum number of missed cleavages allowed.
    fdr_threshold : float
        Maximum FDR threshold (e.g., 0.01 for 1%).
    min_peptide_length : int
        Minimum peptide length in amino acids.
    max_peptide_length : int
        Maximum peptide length in amino acids.
    exclude_sequence_patterns : list[str]
        Regex patterns for sequences to exclude.
    require_unique_peptides : bool
        Whether to require peptides unique to one protein.
    """

    registry: ClassVar[dict[str, "PeptideFilterConfig"]] = {}

    name: str = "default"
    min_search_score: Optional[float] = None
    allowed_charge_states: Optional[List[int]] = None
    exclude_modifications: List[str] = field(default_factory=list)
    max_missed_cleavages: Optional[int] = None
    fdr_threshold: float = 0.01
    min_peptide_length: int = 7
    max_peptide_length: int = 50
    exclude_sequence_patterns: List[str] = field(default_factory=list)
    require_unique_peptides: bool = False

    @classmethod
    def get(cls, name: str, default=None) -> "Optional[PeptideFilterConfig]":
        """Retrieve a configuration from the registry."""
        return cls.registry.get(name.lower(), default)

    def __post_init__(self):
        """Register this configuration in the class registry."""
        if self.name:
            self.registry[self.name.lower()] = self

    def to_dict(self) -> dict:
        """Convert to dictionary, excluding registry."""
        d = asdict(self)
        return d


@dataclass
class ProteinFilterConfig:
    """
    Configuration for protein-level filters.

    Attributes
    ----------
    fdr_threshold : float
        Maximum protein-level FDR threshold.
    min_coverage : float
        Minimum sequence coverage (0.0-1.0).
    min_peptides : int
        Minimum number of peptides per protein.
    min_unique_peptides : int
        Minimum number of unique peptides per protein.
    razor_peptide_handling : str
        How to handle razor peptides: 'keep', 'remove', 'assign_to_top'.
    protein_grouping : str
        Protein grouping strategy: 'none', 'subsumption', 'parsimony'.
    remove_contaminants : bool
        Whether to remove contaminant proteins.
    remove_decoys : bool
        Whether to remove decoy proteins.
    contaminant_patterns : list[str]
        Patterns identifying contaminant proteins.
    """

    registry: ClassVar[dict[str, "ProteinFilterConfig"]] = {}

    name: str = "default"
    fdr_threshold: float = 0.01
    min_coverage: float = 0.0
    min_peptides: int = 1
    min_unique_peptides: int = 2
    razor_peptide_handling: str = "keep"
    protein_grouping: str = "none"
    remove_contaminants: bool = True
    remove_decoys: bool = True
    contaminant_patterns: List[str] = field(
        default_factory=lambda: ["CONTAMINANT", "ENTRAP", "DECOY"]
    )

    @classmethod
    def get(cls, name: str, default=None) -> "Optional[ProteinFilterConfig]":
        """Retrieve a configuration from the registry."""
        return cls.registry.get(name.lower(), default)

    def __post_init__(self):
        """Register this configuration in the class registry."""
        if self.name:
            self.registry[self.name.lower()] = self

    def to_dict(self) -> dict:
        """Convert to dictionary, excluding registry."""
        d = asdict(self)
        return d


@dataclass
class RunQCFilterConfig:
    """
    Configuration for run/sample quality control filters.

    Attributes
    ----------
    min_total_intensity : float
        Minimum total intensity for a run to be included.
    min_identified_features : int
        Minimum number of identified features per run.
    min_identified_proteins : int
        Minimum number of identified proteins per run.
    min_sample_correlation : float, optional
        Minimum pairwise correlation between replicates.
    max_missing_rate : float
        Maximum fraction of missing values allowed per run.
    """

    registry: ClassVar[dict[str, "RunQCFilterConfig"]] = {}

    name: str = "default"
    min_total_intensity: float = 0.0
    min_identified_features: int = 0
    min_identified_proteins: int = 0
    min_sample_correlation: Optional[float] = None
    max_missing_rate: float = 1.0

    @classmethod
    def get(cls, name: str, default=None) -> "Optional[RunQCFilterConfig]":
        """Retrieve a configuration from the registry."""
        return cls.registry.get(name.lower(), default)

    def __post_init__(self):
        """Register this configuration in the class registry."""
        if self.name:
            self.registry[self.name.lower()] = self

    def to_dict(self) -> dict:
        """Convert to dictionary, excluding registry."""
        d = asdict(self)
        return d


@dataclass
class PreprocessingFilterConfig:
    """
    Complete preprocessing filter configuration combining all filter types.

    This is the main configuration class that aggregates all filter settings
    and can be loaded from or saved to YAML/JSON files.
    """

    registry: ClassVar[dict[str, "PreprocessingFilterConfig"]] = {}

    name: str = "default"
    intensity: IntensityFilterConfig = field(default_factory=IntensityFilterConfig)
    peptide: PeptideFilterConfig = field(default_factory=PeptideFilterConfig)
    protein: ProteinFilterConfig = field(default_factory=ProteinFilterConfig)
    run_qc: RunQCFilterConfig = field(default_factory=RunQCFilterConfig)

    # Processing options
    enabled: bool = True
    strict_mode: bool = False  # If True, fail on any filter error
    log_filtered_counts: bool = True

    @classmethod
    def get(cls, name: str, default=None) -> "Optional[PreprocessingFilterConfig]":
        """Retrieve a configuration from the registry."""
        return cls.registry.get(name.lower(), default)

    @classmethod
    def from_dict(cls, data: dict) -> "PreprocessingFilterConfig":
        """Create configuration from a dictionary."""
        intensity_data = data.get("intensity", {})
        peptide_data = data.get("peptide", {})
        protein_data = data.get("protein", {})
        run_qc_data = data.get("run_qc", {})

        # Create nested configs with unique names to avoid registry conflicts
        config_name = data.get("name", "custom")

        return cls(
            name=config_name,
            intensity=IntensityFilterConfig(
                name=f"{config_name}_intensity", **intensity_data
            ),
            peptide=PeptideFilterConfig(name=f"{config_name}_peptide", **peptide_data),
            protein=ProteinFilterConfig(name=f"{config_name}_protein", **protein_data),
            run_qc=RunQCFilterConfig(name=f"{config_name}_run_qc", **run_qc_data),
            enabled=data.get("enabled", True),
            strict_mode=data.get("strict_mode", False),
            log_filtered_counts=data.get("log_filtered_counts", True),
        )

    def __post_init__(self):
        """Register this configuration in the class registry."""
        if self.name:
            self.registry[self.name.lower()] = self

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "intensity": self.intensity.to_dict(),
            "peptide": self.peptide.to_dict(),
            "protein": self.protein.to_dict(),
            "run_qc": self.run_qc.to_dict(),
            "enabled": self.enabled,
            "strict_mode": self.strict_mode,
            "log_filtered_counts": self.log_filtered_counts,
        }

    def apply_overrides(self, overrides: dict) -> None:
        """
        Apply CLI overrides to the configuration.

        Parameters
        ----------
        overrides : dict
            Dictionary of override values. Keys should use dot notation
            (e.g., 'intensity.min_intensity') or flat names
            (e.g., 'min_intensity' for intensity filters).
        """
        # Intensity overrides
        if "min_intensity" in overrides and overrides["min_intensity"] is not None:
            self.intensity.min_intensity = overrides["min_intensity"]
        if "cv_threshold" in overrides and overrides["cv_threshold"] is not None:
            self.intensity.cv_threshold = overrides["cv_threshold"]
        if (
            "min_replicate_agreement" in overrides
            and overrides["min_replicate_agreement"] is not None
        ):
            self.intensity.min_replicate_agreement = overrides[
                "min_replicate_agreement"
            ]

        # Peptide overrides
        if "charge_states" in overrides and overrides["charge_states"] is not None:
            self.peptide.allowed_charge_states = overrides["charge_states"]
        if (
            "max_missed_cleavages" in overrides
            and overrides["max_missed_cleavages"] is not None
        ):
            self.peptide.max_missed_cleavages = overrides["max_missed_cleavages"]
        if (
            "min_peptide_length" in overrides
            and overrides["min_peptide_length"] is not None
        ):
            self.peptide.min_peptide_length = overrides["min_peptide_length"]
        if (
            "exclude_modifications" in overrides
            and overrides["exclude_modifications"] is not None
        ):
            self.peptide.exclude_modifications = overrides["exclude_modifications"]

        # Protein overrides
        if "protein_fdr" in overrides and overrides["protein_fdr"] is not None:
            self.protein.fdr_threshold = overrides["protein_fdr"]
        if (
            "min_unique_peptides" in overrides
            and overrides["min_unique_peptides"] is not None
        ):
            self.protein.min_unique_peptides = overrides["min_unique_peptides"]
        if (
            "remove_contaminants" in overrides
            and overrides["remove_contaminants"] is not None
        ):
            self.protein.remove_contaminants = overrides["remove_contaminants"]
        if "remove_decoys" in overrides and overrides["remove_decoys"] is not None:
            self.protein.remove_decoys = overrides["remove_decoys"]

        # Run QC overrides
        if "min_features" in overrides and overrides["min_features"] is not None:
            self.run_qc.min_identified_features = overrides["min_features"]
        if (
            "max_missing_rate" in overrides
            and overrides["max_missing_rate"] is not None
        ):
            self.run_qc.max_missing_rate = overrides["max_missing_rate"]
