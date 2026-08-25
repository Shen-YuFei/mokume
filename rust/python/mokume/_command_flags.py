"""Validated keyword-to-CLI translation for the thin Python wrappers."""

from __future__ import annotations

_ALLOWED = {
    "features2peptides": set(
        "parquet sdrf min_aa min_unique min_unique_peptides keep_shared_peptides "
        "remove_ids remove_decoy_contaminants remove_contaminants "
        "remove_low_frequency_peptides output skip_normalization run_normalization "
        "sample_normalization log2 save_parquet irs_channel irs_autodetect_regex "
        "irs_stat irs_scope aggregation_level filter_config generate_filter_config "
        "filter_min_intensity filter_cv_threshold filter_charge_states "
        "filter_max_missed_cleavages filter_peptide_fdr filter_score "
        "filter_exclude_modifications filter_min_unique_peptides filter_protein_fdr "
        "filter_min_features filter_max_missing_rate".split()
    ),
    "features2proteins": set(
        "parquet msstats psm output output_format sdrf quant_method method min_aa "
        "min_unique min_unique_peptides keep_contaminants remove_contaminants "
        "run_normalization sample_normalization normalization_proteins "
        "normalization_proteins_file fasta fasta_file pibaq_enzyme pibaq_max_aa "
        "pibaq_min_shared pibaq_families_yaml pibaq_min_anchors directlfq_cores "
        "directlfq_num_cores directlfq_min_nonan directlfq_num_samples_quadratic "
        "export_peptides export_ions batch_correction batch_method batch_column "
        "batch_covariates batch_nonparametric batch_parametric batch_mean_only "
        "batch_ref irs irs_reference_samples irs_reference_sample irs_sdrf_column "
        "irs_sdrf_values irs_reference_regex irs_stat irs_remove_reference "
        "coverage_threshold min_sample_correlation sample_correlation_threshold "
        "ratio_fraction_merge impute impute_method impute_quantile impute_shift "
        "impute_scale impute_n_neighbors de differential_expression de_contrasts "
        "de_contrasts_file de_method de_ensemble_methods de_ensemble_min_k "
        "de_log2fc de_log2fc_threshold de_effect_size_gate de_fdr de_fdr_threshold "
        "de_fdr_method de_output memory threads".split()
    ),
    "peptides2protein": set(
        "fasta peptides method enzyme normalize min_aa max_aa tpa ruler ploidy "
        "organism cpc output verbose qc_report threads min_nonan families "
        "families_yaml min_shared min_anchors high_anchor_threshold".split()
    ),
    "correct-batches": set(
        "folder pattern comment sep output sample_id_column protein_id_column "
        "pibaq_raw_column pibaq_corrected_column export_anndata".split()
    ),
    "visualize": {"folder", "pattern", "output"},
    "tissuemap": set(
        "scan_dir output_dir config config_path generate_config tmt_dataset "
        "tmt_datasets n_jobs dpi imputation_method embedding_method".split()
    ),
    "peptides2protein_qc": {
        "protein_table",
        "qc_report",
        "plot_column",
        "tpa",
        "ruler",
    },
    "peptides2protein_pibaq": set(
        "peptides fasta enzyme output min_aa max_aa ploidy organism cpc qc_report "
        "families families_yaml min_shared min_anchors high_anchor_threshold "
        "normalize tpa ruler verbose".split()
    ),
}

_BOOLEAN = {
    "features2peptides": set(
        "keep_shared_peptides remove_decoy_contaminants remove_contaminants "
        "remove_low_frequency_peptides skip_normalization log2 save_parquet".split()
    ),
    "features2proteins": set(
        "keep_contaminants remove_contaminants batch_correction batch_nonparametric "
        "batch_parametric batch_mean_only irs irs_remove_reference impute de "
        "differential_expression".split()
    ),
    "peptides2protein": {"normalize", "tpa", "ruler", "verbose"},
    "correct-batches": {"export_anndata"},
    "visualize": set(),
    "tissuemap": set(),
    "peptides2protein_qc": {"tpa", "ruler"},
    "peptides2protein_pibaq": {"normalize", "tpa", "ruler", "verbose"},
}

_ALIASES = {
    "features2peptides": {
        "min_unique_peptides": "min-unique",
        "remove_contaminants": "remove-decoy-contaminants",
    },
    "features2proteins": {
        "method": "quant-method",
        "min_unique_peptides": "min-unique",
        "normalization_proteins_file": "normalization-proteins",
        "fasta_file": "fasta",
        "pibaq_families_yaml": "pibaq-families",
        "directlfq_num_cores": "threads",
        "sample_correlation_threshold": "min-sample-correlation",
        "differential_expression": "de",
        "de_log2fc_threshold": "de-log2fc",
        "de_fdr_threshold": "de-fdr",
    },
    "peptides2protein": {"families_yaml": "families"},
    "tissuemap": {"config_path": "config", "tmt_datasets": "tmt-dataset"},
    "peptides2protein_pibaq": {"families_yaml": "families"},
}

_INVERTED_BOOLEAN = {
    ("features2proteins", "remove_contaminants"): "keep-contaminants",
    ("features2proteins", "batch_parametric"): "batch-nonparametric",
}

_CSV = {
    "features2peptides": {"filter_charge_states", "filter_exclude_modifications"},
    "features2proteins": {
        "batch_covariates",
        "irs_reference_samples",
        "irs_sdrf_values",
        "de_contrasts",
        "de_ensemble_methods",
    },
}

_REPEAT = {
    "features2proteins": {"irs_reference_sample"},
    "tissuemap": {"tmt_dataset", "tmt_datasets"},
}


def flags_for(command: str, kwargs: dict[str, object]) -> list[str]:
    """Translate validated wrapper kwargs into the command's exact argv."""
    allowed = _ALLOWED[command]
    unknown = sorted(set(kwargs) - allowed)
    if unknown:
        names = ", ".join(unknown)
        raise TypeError(f"{command} got unexpected keyword argument(s): {names}")

    args: list[str] = []
    for key, value in kwargs.items():
        if value is None:
            continue
        inverted = _INVERTED_BOOLEAN.get((command, key))
        if inverted is not None:
            _append_inverted_boolean(args, command, key, value, inverted)
            continue
        flag = "--" + _ALIASES.get(command, {}).get(key, key.replace("_", "-"))
        if key in _BOOLEAN[command]:
            _append_boolean(args, command, key, value, flag)
        elif isinstance(value, (list, tuple)):
            _append_sequence(args, command, key, value, flag)
        elif isinstance(value, bool):
            raise TypeError(f"{command}.{key} expects a value, not bool")
        else:
            args.extend((flag, str(value)))
    return args


def _append_boolean(
    args: list[str], command: str, key: str, value: object, flag: str
) -> None:
    if not isinstance(value, bool):
        raise TypeError(f"{command}.{key} expects bool")
    if value:
        args.append(flag)


def _append_inverted_boolean(
    args: list[str], command: str, key: str, value: object, inverse_flag: str
) -> None:
    if not isinstance(value, bool):
        raise TypeError(f"{command}.{key} expects bool")
    if not value:
        args.append("--" + inverse_flag)


def _append_sequence(
    args: list[str], command: str, key: str, values: list | tuple, flag: str
) -> None:
    if key in _CSV.get(command, set()):
        if values:
            args.extend((flag, ",".join(str(value) for value in values)))
        return
    if key in _REPEAT.get(command, set()):
        for value in values:
            args.extend((flag, str(value)))
        return
    raise TypeError(f"{command}.{key} does not accept a sequence")
