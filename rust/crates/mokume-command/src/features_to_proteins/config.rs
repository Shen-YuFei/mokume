use mokume_core::{
    BatchCorrectionConfig, DifferentialExpressionConfig, DirectLfqConfig, FeatureToProteinsConfig,
    FilterConfig, ImputationConfig, InputConfig, IrsConfig, MaxLfqConfig, MokumeError,
    NormalizationConfig, OutputConfig, OutputFormat, PibaqConfig, QuantMethod, RatioConfig,
    RuntimeConfig,
};

mod differential_expression;
mod imputation;

use super::Features2ProteinsArgs;
use crate::parsers::{QuantMethodArg, DEFAULT_TOPN_PEPTIDES};
use differential_expression::resolve_differential_expression;
use imputation::resolve_imputation;

fn invalid_input(message: impl Into<String>) -> MokumeError {
    MokumeError::InvalidInput {
        message: message.into(),
    }
}

struct QuantificationOptions {
    method: QuantMethod,
    topn_peptides: usize,
    run_normalization: String,
    sample_normalization: String,
}

struct ResolvedOptions {
    quantification: QuantificationOptions,
    batch: BatchCorrectionConfig,
    irs: IrsConfig,
    ratio: RatioConfig,
    imputation: ImputationConfig,
    differential_expression: DifferentialExpressionConfig,
}

pub(super) fn into_config(
    args: Features2ProteinsArgs,
) -> mokume_core::Result<FeatureToProteinsConfig> {
    let quantification = resolve_quantification(&args)?;
    let batch = resolve_batch(&args)?;
    let irs = resolve_irs(&args, quantification.method)?;
    let ratio = resolve_ratio(&args, quantification.method)?;
    let imputation = resolve_imputation(&args)?;
    let differential_expression = resolve_differential_expression(&args, quantification.method)?;
    Ok(build_config(
        &args,
        ResolvedOptions {
            quantification,
            batch,
            irs,
            ratio,
            imputation,
            differential_expression,
        },
    ))
}

fn resolve_quantification(
    args: &Features2ProteinsArgs,
) -> mokume_core::Result<QuantificationOptions> {
    let QuantMethodArg { method, topn } = args.quant_method;
    validate_lfq_options(args, method)?;
    if method == QuantMethod::Pibaq && args.min_unique.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "piBAQ defines its denominator independently; do not pass --min-unique"
                .to_owned(),
        });
    }
    validate_input_for_quantification(args, method)?;
    validate_method_requirements(args, method)?;
    let manages_normalization = matches!(
        method,
        QuantMethod::DirectLfq
            | QuantMethod::Ratio
            | QuantMethod::PeptideCount
            | QuantMethod::SpectralCount
    );
    let (run_default, sample_default) = default_normalization(args, method, manages_normalization);
    let run_normalization = args
        .run_normalization
        .as_deref()
        .map_or(run_default, |method| method.replace('-', "_"));
    let sample_normalization = args
        .sample_normalization
        .as_deref()
        .map_or(sample_default, |method| method.replace('-', ""));
    validate_normalization_options(
        args,
        method,
        manages_normalization,
        &run_normalization,
        &sample_normalization,
    )?;
    Ok(QuantificationOptions {
        method,
        topn_peptides: topn.unwrap_or(DEFAULT_TOPN_PEPTIDES),
        run_normalization,
        sample_normalization,
    })
}

/// Default run and sample normalization for `method`. Global-median takes the
/// median of each sample's detected features, which moves with detection depth
/// and spike-in composition, so MaxLFQ defaults to the hierarchical (DirectLFQ)
/// alignment of shared peptide species unless normalization proteins ask for a
/// median over those proteins.
fn default_normalization(
    args: &Features2ProteinsArgs,
    method: QuantMethod,
    manages_normalization: bool,
) -> (String, String) {
    if manages_normalization {
        return ("none".to_owned(), "none".to_owned());
    }
    let defaults = NormalizationConfig::default();
    if method == QuantMethod::MaxLfq && args.normalization_proteins.is_none() {
        (defaults.run_method, "hierarchical".to_owned())
    } else {
        (defaults.run_method, defaults.sample_method)
    }
}

fn validate_input_for_quantification(
    args: &Features2ProteinsArgs,
    method: QuantMethod,
) -> mokume_core::Result<()> {
    if method == QuantMethod::SpectralCount {
        if args.psm.is_none() || args.parquet.is_none() {
            return Err(MokumeError::InvalidInput {
                message: "spectral_count requires matching QPX inputs via --psm and --parquet"
                    .to_owned(),
            });
        }
    } else if args.psm.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "--psm only applies to --quant-method spectral_count".to_owned(),
        });
    }
    if args.sdrf.is_none()
        && args
            .sample_normalization
            .as_deref()
            .is_some_and(|method| method.eq_ignore_ascii_case("condition-median"))
    {
        return Err(MokumeError::InvalidInput {
            message: "conditionmedian sample normalization requires --sdrf option".to_owned(),
        });
    }
    Ok(())
}

/// Inputs, outputs, and FASTA options that only some quantification methods use.
fn validate_method_requirements(
    args: &Features2ProteinsArgs,
    method: QuantMethod,
) -> mokume_core::Result<()> {
    if method == QuantMethod::Ratio && args.msstats.is_some() {
        return Err(invalid_input(
            "Ratio quantification requires PSM-level QPX input; MSstats feature tables do not contain PSM evidence",
        ));
    }
    if method == QuantMethod::Pibaq && args.fasta.is_none() {
        return Err(invalid_input(
            "piBAQ quantification requires --fasta option",
        ));
    }
    if method != QuantMethod::Pibaq && pibaq_options_supplied(args) {
        return Err(invalid_input(
            "piBAQ FASTA/digestion options require --quant-method pibaq",
        ));
    }
    if method == QuantMethod::Ratio && args.sdrf.is_none() {
        return Err(invalid_input("Ratio quantification requires --sdrf option"));
    }
    if args.export_peptides.is_some()
        && matches!(
            method,
            QuantMethod::DirectLfq | QuantMethod::Ratio | QuantMethod::SpectralCount
        )
    {
        return Err(invalid_input(format!(
            "export-peptides is not supported by {method} quantification"
        )));
    }
    if args.coverage_threshold.is_some() && args.sdrf.is_none() {
        return Err(invalid_input("coverage-threshold requires --sdrf option"));
    }
    Ok(())
}

fn pibaq_options_supplied(args: &Features2ProteinsArgs) -> bool {
    args.fasta.is_some()
        || args.pibaq_enzyme.is_some()
        || args.pibaq_max_aa.is_some()
        || args.pibaq_min_shared.is_some()
        || args.pibaq_families_yaml.is_some()
        || args.pibaq_min_anchors.is_some()
}

fn validate_normalization_options(
    args: &Features2ProteinsArgs,
    method: QuantMethod,
    manages_normalization: bool,
    run_normalization: &str,
    sample_normalization: &str,
) -> mokume_core::Result<()> {
    if args.normalization_proteins.is_some()
        && !matches!(
            sample_normalization.to_ascii_lowercase().as_str(),
            "globalmedian" | "conditionmedian"
        )
    {
        return Err(invalid_input(
            "--normalization-proteins requires globalmedian or conditionmedian sample normalization",
        ));
    }
    if manages_normalization
        && (!run_normalization.eq_ignore_ascii_case("none")
            || !sample_normalization.eq_ignore_ascii_case("none")
            || args.normalization_proteins.is_some())
    {
        let reason = if matches!(
            method,
            QuantMethod::PeptideCount | QuantMethod::SpectralCount
        ) {
            "does not use intensity normalization"
        } else {
            "manages normalization internally"
        };
        return Err(invalid_input(format!(
            "{method} {reason}; use --run-normalization none and --sample-normalization \
             none, and do not pass --normalization-proteins"
        )));
    }
    Ok(())
}

fn validate_lfq_options(
    args: &Features2ProteinsArgs,
    method: QuantMethod,
) -> mokume_core::Result<()> {
    if method != QuantMethod::DirectLfq && args.export_ions.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "--export-ions requires --quant-method directlfq".to_owned(),
        });
    }
    if method != QuantMethod::DirectLfq && args.directlfq_min_nonan.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "--directlfq-min-nonan requires --quant-method directlfq".to_owned(),
        });
    }
    if method != QuantMethod::DirectLfq && args.directlfq_num_samples_quadratic.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "--directlfq-num-samples-quadratic only applies to DirectLFQ".to_owned(),
        });
    }
    if method != QuantMethod::DirectLfq && args.directlfq_no_sample_normalization {
        return Err(MokumeError::InvalidInput {
            message: "--directlfq-no-sample-normalization only applies to DirectLFQ".to_owned(),
        });
    }
    if method != QuantMethod::MaxLfq && args.maxlfq_min_ratio_count.is_some() {
        return Err(MokumeError::InvalidInput {
            message: "--maxlfq-min-ratio-count requires --quant-method maxlfq".to_owned(),
        });
    }
    if method != QuantMethod::MaxLfq && args.stabilize {
        return Err(MokumeError::InvalidInput {
            message: "--stabilize requires --quant-method maxlfq".to_owned(),
        });
    }
    Ok(())
}

fn resolve_batch(args: &Features2ProteinsArgs) -> mokume_core::Result<BatchCorrectionConfig> {
    let method_supplied = args.batch_method.is_some();
    let method = args.batch_method.clone().map_or_else(
        || BatchCorrectionConfig::default().method,
        |value| value.replace('-', "_"),
    );
    if !args.batch_correction
        && (method_supplied
            || args.batch_column.is_some()
            || !args.batch_covariate.is_empty()
            || args.batch_nonparametric
            || args.batch_mean_only
            || args.batch_ref.is_some())
    {
        return Err(MokumeError::InvalidInput {
            message: "batch options require --batch-correction".to_owned(),
        });
    }
    if args.batch_correction {
        validate_batch_columns(args, &method)?;
    }
    Ok(BatchCorrectionConfig {
        enabled: args.batch_correction,
        method,
        column: args.batch_column.clone(),
        covariates: (!args.batch_covariate.is_empty()).then(|| args.batch_covariate.clone()),
        parametric: !args.batch_nonparametric,
        mean_only: args.batch_mean_only,
        ref_batch: args.batch_ref.clone(),
    })
}

fn validate_batch_columns(args: &Features2ProteinsArgs, method: &str) -> mokume_core::Result<()> {
    let column_method = method.eq_ignore_ascii_case("column");
    if column_method && args.batch_column.is_none() {
        return Err(invalid_input(
            "Batch correction with method 'column' requires --batch-column option",
        ));
    }
    if !column_method && args.batch_column.is_some() {
        return Err(invalid_input(
            "--batch-column requires --batch-method column",
        ));
    }
    if args.sdrf.is_none() && (args.batch_column.is_some() || !args.batch_covariate.is_empty()) {
        return Err(invalid_input(
            "Batch correction with --batch-column or --batch-covariate requires --sdrf option",
        ));
    }
    Ok(())
}

fn resolve_irs(
    args: &Features2ProteinsArgs,
    quantification: QuantMethod,
) -> mokume_core::Result<IrsConfig> {
    let reference_samples =
        (!args.irs_reference_sample.is_empty()).then(|| args.irs_reference_sample.clone());
    let selector_count = validate_irs_selectors(args, reference_samples.is_some())?;
    validate_irs_mode(args, quantification, selector_count)?;
    let defaults = IrsConfig::default();
    Ok(IrsConfig {
        enabled: args.irs,
        reference_samples,
        sdrf_column: args.irs_sdrf_column.clone(),
        sdrf_values: (!args.irs_sdrf_value.is_empty()).then(|| args.irs_sdrf_value.clone()),
        reference_regex: args
            .irs_reference_regex
            .clone()
            .unwrap_or(defaults.reference_regex),
        stat: args.irs_stat.clone().unwrap_or(defaults.stat),
        remove_reference: args.irs_remove_reference,
    })
}

fn validate_irs_selectors(
    args: &Features2ProteinsArgs,
    has_reference_samples: bool,
) -> mokume_core::Result<usize> {
    if args.irs_sdrf_column.is_some() == args.irs_sdrf_value.is_empty() {
        return Err(MokumeError::InvalidInput {
            message: "--irs-sdrf-column and --irs-sdrf-value must be provided together".to_owned(),
        });
    }
    let selector_count = usize::from(has_reference_samples)
        + usize::from(args.irs_sdrf_column.is_some())
        + usize::from(args.irs_reference_regex.is_some());
    if selector_count > 1 {
        return Err(MokumeError::InvalidInput {
            message: "choose one reference selector: samples, SDRF column+values, or regex"
                .to_owned(),
        });
    }
    Ok(selector_count)
}

fn validate_irs_mode(
    args: &Features2ProteinsArgs,
    quantification: QuantMethod,
    selector_count: usize,
) -> mokume_core::Result<()> {
    if matches!(
        quantification,
        QuantMethod::PeptideCount | QuantMethod::SpectralCount
    ) && args.irs
    {
        return Err(MokumeError::InvalidInput {
            message: format!("{quantification} quantification cannot apply IRS"),
        });
    }
    if quantification == QuantMethod::Ratio {
        validate_ratio_reference_options(args)?;
    } else if !args.irs
        && (selector_count > 0 || args.irs_stat.is_some() || args.irs_remove_reference)
    {
        return Err(MokumeError::InvalidInput {
            message: "IRS options require --irs".to_owned(),
        });
    }
    if args.irs && args.sdrf.is_none() {
        return Err(invalid_input("IRS options require --sdrf option"));
    }
    Ok(())
}

/// Ratio takes its reference samples from the IRS selectors but runs no IRS.
fn validate_ratio_reference_options(args: &Features2ProteinsArgs) -> mokume_core::Result<()> {
    if args.irs {
        return Err(MokumeError::InvalidInput {
            message: "Ratio quantification cannot also apply IRS".to_owned(),
        });
    }
    if args.irs_sdrf_column.is_some()
        || !args.irs_sdrf_value.is_empty()
        || args.irs_stat.is_some()
        || args.irs_remove_reference
    {
        return Err(MokumeError::InvalidInput {
            message: "Ratio accepts --irs-reference-sample or --irs-reference-regex; IRS-only options require --irs"
                .to_owned(),
        });
    }
    Ok(())
}

fn resolve_ratio(
    args: &Features2ProteinsArgs,
    quantification: QuantMethod,
) -> mokume_core::Result<RatioConfig> {
    if args.ratio_fraction_merge.is_some() && quantification != QuantMethod::Ratio {
        return Err(MokumeError::InvalidInput {
            message: "--ratio-fraction-merge only applies to --quant-method ratio".to_owned(),
        });
    }
    Ok(RatioConfig {
        fraction_merge: args
            .ratio_fraction_merge
            .clone()
            .unwrap_or_else(|| RatioConfig::default().fraction_merge),
    })
}

fn build_config(
    args: &Features2ProteinsArgs,
    resolved: ResolvedOptions,
) -> FeatureToProteinsConfig {
    let quantification = resolved.quantification;
    FeatureToProteinsConfig {
        input: input_config(args),
        output: output_config(args),
        filtering: FilterConfig {
            min_aa: args.min_aa,
            min_unique_peptides: args.min_unique.unwrap_or(
                if quantification.method == QuantMethod::Pibaq {
                    0
                } else {
                    FilterConfig::default().min_unique_peptides
                },
            ),
            remove_contaminants: !args.keep_contaminants,
        },
        normalization: NormalizationConfig {
            run_method: quantification.run_normalization,
            sample_method: quantification.sample_normalization,
            normalization_proteins: args.normalization_proteins.clone(),
        },
        quantification: quantification.method,
        topn_peptides: quantification.topn_peptides,
        maxlfq: MaxLfqConfig {
            min_ratio_count: args
                .maxlfq_min_ratio_count
                .unwrap_or(MaxLfqConfig::default().min_ratio_count),
            stabilize: args.stabilize,
        },
        pibaq: pibaq_config(args),
        directlfq: directlfq_config(args),
        batch: resolved.batch,
        irs: resolved.irs,
        coverage_threshold: args.coverage_threshold,
        sample_correlation_threshold: args.min_sample_correlation,
        ratio: resolved.ratio,
        imputation: resolved.imputation,
        differential_expression: resolved.differential_expression,
        runtime: RuntimeConfig {
            memory: args.memory.clone(),
            threads: args.threads,
        },
    }
}

fn input_config(args: &Features2ProteinsArgs) -> InputConfig {
    InputConfig {
        parquet: args.parquet.clone(),
        msstats: args.msstats.clone(),
        psm: args.psm.clone(),
        sdrf: args.sdrf.clone(),
        fasta: args.fasta.clone(),
    }
}

fn output_config(args: &Features2ProteinsArgs) -> OutputConfig {
    OutputConfig {
        protein_matrix: args.output.clone(),
        export_peptides: args.export_peptides.clone(),
        export_ions: args.export_ions.clone(),
        format: OutputFormat::PythonCompatible,
    }
}

fn pibaq_config(args: &Features2ProteinsArgs) -> PibaqConfig {
    let defaults = PibaqConfig::default();
    PibaqConfig {
        enzyme: args.pibaq_enzyme.clone().unwrap_or(defaults.enzyme),
        max_aa: args.pibaq_max_aa.unwrap_or(defaults.max_aa),
        min_shared: args.pibaq_min_shared.unwrap_or(defaults.min_shared),
        families_yaml: args.pibaq_families_yaml.clone(),
        min_anchors: args.pibaq_min_anchors.unwrap_or(defaults.min_anchors),
        high_anchor_threshold: defaults.high_anchor_threshold,
    }
}

fn directlfq_config(args: &Features2ProteinsArgs) -> DirectLfqConfig {
    let defaults = DirectLfqConfig::default();
    DirectLfqConfig {
        min_nonan: args.directlfq_min_nonan.unwrap_or(defaults.min_nonan),
        num_samples_quadratic: args
            .directlfq_num_samples_quadratic
            .unwrap_or(defaults.num_samples_quadratic),
        normalize_samples: !args.directlfq_no_sample_normalization,
    }
}
