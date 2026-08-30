pub(crate) const EXAMPLE_FILTER_CONFIG_YAML: &str = r#"# Mokume Preprocessing Filter Configuration
# This file defines quality filters applied during peptide normalization

name: example_config

# Global options
enabled: true              # Set to false to disable all filtering

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
  allowed_charge_states: null     # e.g., [2, 3, 4] or null for all charges
  exclude_modifications: []       # Modification names to exclude, e.g., ["Oxidation"]
  max_missed_cleavages: null      # Max missed cleavages (null = no filter)
  fdr_threshold: null              # Peptide q-value cutoff (null = no filter)
  min_peptide_length: 7           # Minimum peptide length in amino acids
  max_peptide_length: 50          # Maximum peptide length in amino acids
  exclude_sequence_patterns: []   # Regex patterns to exclude

# Protein-level filters
protein:
  fdr_threshold: null           # Protein-group q-value cutoff (null = no filter)
  min_unique_peptides: 2      # Minimum unique peptides per protein
  razor_peptide_handling: keep   # How to handle shared peptides: keep, remove, assign_to_top
  remove_contaminants: true      # Remove contaminant proteins
  remove_decoys: true            # Remove decoy proteins
  contaminant_patterns:          # Patterns identifying contaminants
    - CONTAMINANT
    - CONTAM_
    - ENTRAP
    - DECOY

# Run/Sample QC filters
run_qc:
  min_total_intensity: 0.0      # Min total intensity per run
  min_identified_features: 0    # Min features per run
  min_identified_proteins: 0    # Min proteins per run
  max_missing_rate: 1.0         # Max missing value rate per run (0-1)
"#;

pub(crate) const EXAMPLE_FILTER_CONFIG_JSON: &str = r#"{
  "name": "example_config",
  "intensity": {
    "min_intensity": 0.0,
    "cv_threshold": null,
    "min_replicate_agreement": 1,
    "quantile_lower": 0.0,
    "quantile_upper": 1.0,
    "remove_zero_intensity": true
  },
  "peptide": {
    "allowed_charge_states": null,
    "exclude_modifications": [],
    "max_missed_cleavages": null,
    "fdr_threshold": null,
    "min_peptide_length": 7,
    "max_peptide_length": 50,
    "exclude_sequence_patterns": []
  },
  "protein": {
    "fdr_threshold": null,
    "min_unique_peptides": 2,
    "razor_peptide_handling": "keep",
    "remove_contaminants": true,
    "remove_decoys": true,
    "contaminant_patterns": [
      "CONTAMINANT",
      "CONTAM_",
      "ENTRAP",
      "DECOY"
    ]
  },
  "run_qc": {
    "min_total_intensity": 0.0,
    "min_identified_features": 0,
    "min_identified_proteins": 0,
    "max_missing_rate": 1.0
  },
  "enabled": true
}
"#;
