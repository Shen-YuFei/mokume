"""
Configuration file I/O for preprocessing filters.
"""

import json
from pathlib import Path
from typing import Union, Optional

from mokume.model.filters import PreprocessingFilterConfig
from mokume.core.logger import get_logger


logger = get_logger("mokume.preprocessing.filters.io")

# Lazy import for optional PyYAML
_YAML_AVAILABLE = None


def _check_yaml_available() -> bool:
    """Check if PyYAML is installed."""
    global _YAML_AVAILABLE
    if _YAML_AVAILABLE is None:
        try:
            import yaml  # noqa: F401

            _YAML_AVAILABLE = True
        except ImportError:
            _YAML_AVAILABLE = False
    return _YAML_AVAILABLE


def _import_yaml():
    """Import yaml with helpful error message."""
    if not _check_yaml_available():
        raise ImportError(
            "YAML configuration support requires the 'pyyaml' package.\n"
            "Install with: pip install mokume-py[directlfq]\n"
            "Or: pip install pyyaml"
        )
    import yaml

    return yaml


def is_filter_config_available() -> bool:
    """
    Check if filter configuration loading is available.

    JSON is always available, YAML requires pyyaml.

    Returns
    -------
    bool
        True (JSON is always available).
    """
    return True


def load_filter_config(config_path: Union[str, Path]) -> PreprocessingFilterConfig:
    """
    Load filter configuration from YAML or JSON file.

    Parameters
    ----------
    config_path : Union[str, Path]
        Path to configuration file (.yaml, .yml, or .json).

    Returns
    -------
    PreprocessingFilterConfig
        Loaded filter configuration.

    Raises
    ------
    ValueError
        If file format is not supported.
    FileNotFoundError
        If config file does not exist.
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    suffix = config_path.suffix.lower()

    if suffix in (".yaml", ".yml"):
        yaml = _import_yaml()
        with open(config_path, "r") as f:
            data = yaml.safe_load(f)
    elif suffix == ".json":
        with open(config_path, "r") as f:
            data = json.load(f)
    else:
        raise ValueError(
            f"Unsupported config format: {suffix}. Use .yaml, .yml, or .json"
        )

    logger.info("Loaded filter configuration from %s", config_path)
    return PreprocessingFilterConfig.from_dict(data)


def save_filter_config(
    config: PreprocessingFilterConfig,
    output_path: Union[str, Path],
    format: Optional[str] = None,
) -> None:
    """
    Save filter configuration to YAML or JSON file.

    Parameters
    ----------
    config : PreprocessingFilterConfig
        Filter configuration to save.
    output_path : Union[str, Path]
        Output file path.
    format : str, optional
        Output format ('yaml' or 'json'). Inferred from extension if not provided.
    """
    output_path = Path(output_path)

    if format is None:
        suffix = output_path.suffix.lower()
        if suffix in (".yaml", ".yml"):
            format = "yaml"
        elif suffix == ".json":
            format = "json"
        else:
            format = "yaml"  # Default to YAML

    data = config.to_dict()

    # Remove internal name fields from nested configs
    for key in ["intensity", "peptide", "protein", "run_qc"]:
        if key in data and "name" in data[key]:
            del data[key]["name"]

    if format == "yaml":
        yaml = _import_yaml()
        with open(output_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    else:
        with open(output_path, "w") as f:
            json.dump(data, f, indent=2)

    logger.info("Saved filter configuration to %s", output_path)


def generate_example_config(
    output_path: Union[str, Path],
    format: Optional[str] = None,
) -> None:
    """
    Generate an example configuration file with default values and comments.

    Parameters
    ----------
    output_path : Union[str, Path]
        Output file path.
    format : str, optional
        Output format ('yaml' or 'json').
    """
    output_path = Path(output_path)

    if format is None:
        suffix = output_path.suffix.lower()
        if suffix in (".yaml", ".yml"):
            format = "yaml"
        elif suffix == ".json":
            format = "json"
        else:
            format = "yaml"

    if format == "yaml":
        # Generate YAML with comments
        yaml_content = """# Mokume Preprocessing Filter Configuration
# This file defines quality filters applied during peptide normalization

name: example_config

# Global options
enabled: true              # Set to false to disable all filtering
strict_mode: false         # If true, fail on any filter error
log_filtered_counts: true  # Log how many items each filter removes

# Intensity-based filters
intensity:
  min_intensity: 0.0           # Minimum intensity threshold (0 = no filter)
  cv_threshold: null           # Maximum CV across replicates (null = no filter)
  min_replicate_agreement: 1   # Min replicates where feature must be detected
  quantile_lower: 0.0          # Lower quantile for outlier removal (0-1)
  quantile_upper: 1.0          # Upper quantile for outlier removal (0-1)
  remove_zero_intensity: true  # Remove features with zero intensity

# Peptide-level filters
peptide:
  min_search_score: null          # Min search engine score (null = no filter)
  allowed_charge_states: null     # e.g., [2, 3, 4] or null for all charges
  exclude_modifications: []       # Modification names to exclude, e.g., ["Oxidation"]
  max_missed_cleavages: null      # Max missed cleavages (null = no filter)
  fdr_threshold: 0.01             # Peptide FDR threshold (requires q_value column)
  min_peptide_length: 7           # Minimum peptide length in amino acids
  max_peptide_length: 50          # Maximum peptide length in amino acids
  exclude_sequence_patterns: []   # Regex patterns to exclude
  require_unique_peptides: false  # Require peptides unique to one protein

# Protein-level filters
protein:
  fdr_threshold: 0.01         # Protein FDR threshold
  min_coverage: 0.0           # Minimum sequence coverage (0-1)
  min_peptides: 1             # Minimum total peptides per protein
  min_unique_peptides: 2      # Minimum unique peptides per protein
  razor_peptide_handling: keep   # How to handle shared peptides: keep, remove, assign_to_top
  protein_grouping: none         # Grouping strategy: none, subsumption, parsimony
  remove_contaminants: true      # Remove contaminant proteins
  remove_decoys: true            # Remove decoy proteins
  contaminant_patterns:          # Patterns identifying contaminants
    - CONTAMINANT
    - ENTRAP
    - DECOY

# Run/Sample QC filters
run_qc:
  min_total_intensity: 0.0      # Min total intensity per run
  min_identified_features: 0    # Min features per run
  min_identified_proteins: 0    # Min proteins per run
  min_sample_correlation: null  # Min correlation between samples (null = no filter)
  max_missing_rate: 1.0         # Max missing value rate (0-1)
"""
        with open(output_path, "w") as f:
            f.write(yaml_content)
    else:
        # Generate JSON
        config = PreprocessingFilterConfig(name="example_config")
        save_filter_config(config, output_path, format="json")

    logger.info("Generated example filter configuration at %s", output_path)
