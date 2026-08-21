//! Differential-expression step for the features2proteins pipeline.
//!
//! This mirrors `mokume.pipeline.stages.PostprocessingStage
//! .run_differential_expression` for limma, deqms, rots, limrots, proDA, and
//! ensemble. RNG/optimizer-driven methods are faithful ports rather than
//! cell-exact p-value reproductions. The flow per contrast:
//!   1. derive a sample -> condition map from the SDRF `factor value[*]`
//!      column (with the same `_shorten_factor_label` reduction the Python
//!      `detect_condition_from_sdrf` applies);
//!   2. split the protein matrix sample columns into the two contrast groups;
//!   3. reject infinite intensities, then log2-transform positive finite cells;
//!      zero, non-positive, and missing values become NaN;
//!   4. run the configured test, apply BH/IHW/BKY/Storey as appropriate, resolve
//!      the fixed or data-driven effect-size gate, classify, and sort;
//!   5. write one DE result table per contrast (the column set differs by
//!      method: limma carries the moderated-fit `t_stat`/`AveExpr`/`B` extras,
//!      deqms carries the spectraCounteBayes
//!      `sca_t`/`sca_pvalue`/`sca_adj_pvalue` plus a `peptide_count` column, and
//!      rots carries a single `d_stat` column).
//!
//! `--de-method ensemble` is also wired: it runs each configured member method
//! (default members `[limrots, deqms, proda]`, overridable via
//! `--de-ensemble-methods`), keeps the non-empty results, and fuses them with the
//! deterministic top-k consensus combiner (`combine_de_results`), emitting the
//! consensus column set (`n_methods_up`/`n_methods_down`/`methods_significant`).
//! The combine step is cell-exact vs Python; the end-to-end ensemble inherits
//! the RNG/optimizer non-determinism of whichever members it runs.
//!
//! BH, IHW, BKY, and Storey are wired. BKY/Storey use the Python reliability
//! gate and fall back to BH when pi0 is not trustworthy; ROTS/LimROTS preserve
//! their permutation FDR. Ensemble members use the requested correction and the
//! combined Fisher p-values are corrected again, matching Python. A fixed
//! `--de-log2fc` threshold or the Python-compatible `auto` mixture gate is
//! applied to both member and ensemble classifications; `null_quantile` is also
//! available explicitly.
//!
//! Unsupported methods and options are rejected upstream in
//! `validate_de_subset`; this module never silently substitutes a different
//! test. The in-pipeline deqms path feeds the per-protein
//! unique-canonical-peptide counts captured at ingest (mirroring Python's
//! `_load_de_peptide_counts` -> `_build_count_vector`), so the
//! spectraCounteBayes count moderation runs when the data has a real spread of
//! counts and at least 10 testable proteins; with all-equal counts (or fewer
//! than 10 valid points) it falls back to plain eBayes (== limma), exactly as
//! mokume does.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use mokume_core::{DifferentialExpressionConfig, MokumeError, Result, SampleId};
use mokume_io::SdrfTable;
use mokume_stats::de::{
    adaptive_adjust, combine_de_results, deqms_two_group, estimate_effect_size_gate,
    ihw_correction, limma_two_group, limrots_two_group, proda_two_group, rots_two_group,
    AdaptiveFdrMethod, AppliedFdrMethod, DeResult, EffectSizeGateMethod, EnsembleResult,
    Significance,
};

use crate::{create_parent_dir, invalid_input, ProteinMatrix};

/// A parsed contrast: condition A versus condition B.
struct Contrast {
    cond_a: String,
    cond_b: String,
}

impl Contrast {
    /// Key used for the output filename when more than one contrast runs,
    /// matching the Python `f"{contrast[0]}-{contrast[1]}"`.
    fn key(&self) -> String {
        format!("{}-{}", self.cond_a, self.cond_b)
    }
}

/// Run differential expression for every configured contrast and write the
/// result tables. Called from `run_features_to_proteins` after the protein
/// matrix is finalized (post IRS / coverage filter / imputation), mirroring
/// the Python stage order.
pub(crate) fn run_differential_expression(
    matrix: &ProteinMatrix,
    sdrf: Option<&SdrfTable>,
    config: &DifferentialExpressionConfig,
    drop_empty_samples: bool,
) -> Result<()> {
    let Some(sdrf) = sdrf else {
        return Err(invalid_input(
            "differential expression requires an SDRF file (--sdrf)",
        ));
    };
    let contrasts = parse_contrasts(config)?;
    if contrasts.is_empty() {
        return Ok(());
    }

    let condition_by_sample = condition_by_sample(sdrf);
    let samples = matrix.sample_columns(drop_empty_samples);

    for contrast in &contrasts {
        let results = run_one_contrast(matrix, &samples, &condition_by_sample, contrast, config)?;
        if let Some(output) = de_output_path(config, contrast, contrasts.len()) {
            match &results {
                MatrixDifferentialExpressionResults::Standard(rows) => {
                    let DeMethod::Member(method) = de_method(config)? else {
                        return Err(invalid_input(
                            "ensemble results cannot use the single-method writer",
                        ));
                    };
                    write_de_csv(&output, rows, contrast, method, matrix)?;
                }
                MatrixDifferentialExpressionResults::Ensemble(rows) => {
                    write_ensemble_csv(&output, rows, contrast)?;
                }
            }
        }
    }
    Ok(())
}

/// A contrast's DE output: either a single-method table (`Vec<DeResult>`) or the
/// ensemble's combined table (`Vec<EnsembleResult>`, carrying the extra
/// consensus columns). Both carry the requested FDR adjustment (or the
/// documented fallback) and are sorted by adjusted p-value.
#[derive(Debug, Clone)]
pub enum MatrixDifferentialExpressionResults {
    Standard(Vec<DeResult>),
    Ensemble(Vec<EnsembleResult>),
}

/// Run one two-group differential-expression analysis on a row-major linear
/// protein matrix without loading QPX or writing a result file.
///
/// The first `n_a` columns belong to group A and the next `n_b` columns to
/// group B. Infinite intensities are rejected; non-positive values and NaN are
/// treated as missing, and positive finite values are log2-transformed,
/// matching the full pipeline. `auto` is intentionally rejected because a
/// matrix-only call has no quantification method from which to resolve it.
pub fn differential_expression_matrix(
    proteins: &[String],
    values: &[Vec<f64>],
    n_a: usize,
    n_b: usize,
    peptide_counts: Option<&[f64]>,
    config: &DifferentialExpressionConfig,
    threads: Option<usize>,
) -> Result<MatrixDifferentialExpressionResults> {
    validate_config(config, false)?;
    validate_matrix_de_inputs(proteins, values, n_a, n_b, peptide_counts)?;

    crate::threading::install(threads, || {
        let rows = log2_matrix(values);
        let row_refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let default_counts;
        let counts = match peptide_counts {
            Some(counts) => counts,
            None => {
                default_counts = vec![1.0; values.len()];
                &default_counts
            }
        };
        let prepared = Prepared {
            proteins,
            rows: &row_refs,
            peptide_counts: counts,
            n_a,
            n_b,
        };
        match de_method(config)? {
            DeMethod::Ensemble => Ok(MatrixDifferentialExpressionResults::Ensemble(run_ensemble(
                &prepared, config,
            )?)),
            DeMethod::Member(method) => Ok(MatrixDifferentialExpressionResults::Standard(
                run_member(method, &prepared, config),
            )),
        }
    })
}

fn validate_matrix_de_inputs(
    proteins: &[String],
    values: &[Vec<f64>],
    n_a: usize,
    n_b: usize,
    peptide_counts: Option<&[f64]>,
) -> Result<()> {
    if n_a == 0 || n_b == 0 {
        return Err(invalid_input(
            "matrix differential expression requires at least one sample per group",
        ));
    }
    let width = n_a
        .checked_add(n_b)
        .ok_or_else(|| invalid_input("group sample counts overflow matrix width"))?;
    let actual_width = crate::matrix::validate_rectangular(values)?;
    crate::matrix::validate_no_infinite(values, "differential-expression")?;
    if actual_width != width {
        return Err(invalid_input(format!(
            "matrix has {actual_width} columns but group sizes require {width}"
        )));
    }
    if proteins.len() != values.len() {
        return Err(invalid_input(format!(
            "matrix has {} rows but {} protein names were supplied",
            values.len(),
            proteins.len()
        )));
    }
    if let Some(counts) = peptide_counts {
        if counts.len() != values.len() {
            return Err(invalid_input(format!(
                "matrix has {} rows but {} peptide counts were supplied",
                values.len(),
                counts.len()
            )));
        }
    }
    Ok(())
}

fn log2_matrix(values: &[Vec<f64>]) -> Vec<Vec<f64>> {
    values
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| {
                    if value.is_finite() && *value > 0.0 {
                        value.log2()
                    } else {
                        f64::NAN
                    }
                })
                .collect()
        })
        .collect()
}

/// Run a single contrast and return the sorted DE table (single-method or
/// ensemble).
fn run_one_contrast(
    matrix: &ProteinMatrix,
    samples: &[(SampleId, &str)],
    condition_by_sample: &HashMap<String, String>,
    contrast: &Contrast,
    config: &DifferentialExpressionConfig,
) -> Result<MatrixDifferentialExpressionResults> {
    let samples_a = group_samples(samples, condition_by_sample, &contrast.cond_a);
    let samples_b = group_samples(samples, condition_by_sample, &contrast.cond_b);
    if samples_a.is_empty() {
        return Err(invalid_input(format!(
            "No samples for '{}'. Available: {}",
            contrast.cond_a,
            available_conditions(condition_by_sample)
        )));
    }
    if samples_b.is_empty() {
        return Err(invalid_input(format!(
            "No samples for '{}'. Available: {}",
            contrast.cond_b,
            available_conditions(condition_by_sample)
        )));
    }

    // Build the log2 protein x sample matrix with group A columns first, then
    // group B, matching the column order the two-group methods expect.
    let ordered: Vec<SampleId> = samples_a.iter().chain(samples_b.iter()).copied().collect();
    let (proteins, rows) = matrix.log2_rows(&ordered)?;
    let row_refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let peptide_counts = build_count_vector(matrix, &proteins);
    let prepared = Prepared {
        proteins: &proteins,
        rows: &row_refs,
        peptide_counts: &peptide_counts,
        n_a: samples_a.len(),
        n_b: samples_b.len(),
    };

    // Method is validated upstream in `validate_de_subset`; only `limma`,
    // `deqms`, `rots`, `limrots`, `proda`, and `ensemble` reach here.
    match de_method(config)? {
        DeMethod::Ensemble => Ok(MatrixDifferentialExpressionResults::Ensemble(run_ensemble(
            &prepared, config,
        )?)),
        DeMethod::Member(method) => Ok(MatrixDifferentialExpressionResults::Standard(run_member(
            method, &prepared, config,
        ))),
    }
}

/// The per-contrast inputs shared by every single-method run: the matrix (for
/// deqms counts), the contrast-ordered protein names, the log2 rows, and the two
/// group sizes.
struct Prepared<'a> {
    proteins: &'a [String],
    rows: &'a [&'a [f64]],
    peptide_counts: &'a [f64],
    n_a: usize,
    n_b: usize,
}

/// Run one concrete single-method DE test on the prepared contrast data,
/// returning its configured-FDR-adjusted, classified, sorted `Vec<DeResult>`.
fn run_member(
    method: MemberMethod,
    prepared: &Prepared<'_>,
    config: &DifferentialExpressionConfig,
) -> Vec<DeResult> {
    let results = run_method_bh(method, prepared, config);
    finalize_member_results(method, results, config)
}

fn finalize_member_results(
    method: MemberMethod,
    mut results: Vec<DeResult>,
    config: &DifferentialExpressionConfig,
) -> Vec<DeResult> {
    // LimROTS and ROTS preserve their own permutation FDR. All other methods
    // apply the configured correction over the raw p-values.
    if !matches!(method, MemberMethod::Limrots | MemberMethod::Rots) {
        if is_ihw(config) {
            apply_ihw(&mut results, config);
        } else if let Some(adaptive_method) = adaptive_fdr_method(config) {
            apply_adaptive_fdr(&mut results, adaptive_method, config);
        }
    }
    let gate = resolved_effect_size_gate(&results, config);
    classify_and_sort(&mut results, gate, config);
    results
}

/// Run one single-method DE test with its built-in Benjamini-Hochberg
/// correction (the method's own `bh_adjust` + sort). The configured alternative
/// correction is layered on by [`run_member`].
fn run_method_bh(
    method: MemberMethod,
    prepared: &Prepared<'_>,
    config: &DifferentialExpressionConfig,
) -> Vec<DeResult> {
    let (proteins, rows, n_a, n_b) = (prepared.proteins, prepared.rows, prepared.n_a, prepared.n_b);
    let (fdr, log2fc) = (config.fdr_threshold, config.log2fc_threshold);
    match method {
        MemberMethod::Rots => rots_two_group(proteins, rows, n_a, n_b, fdr, log2fc),
        MemberMethod::Limrots => limrots_two_group(proteins, rows, n_a, n_b, fdr, log2fc),
        MemberMethod::Proda => proda_two_group(proteins, rows, n_a, n_b, fdr, log2fc),
        MemberMethod::Deqms => {
            // Per-protein unique-canonical-peptide counts captured at ingest,
            // aligned to the contrast's protein order and defaulting to 1 for any
            // protein absent from the map -- mirroring Python's
            // `_build_count_vector` (deqms.py:147), which reindexes the
            // `_load_de_peptide_counts` Series onto the matrix protein order and
            // fills missing entries with 1. With a real spread of counts the
            // spectraCounteBayes moderation runs; an all-equal/all-ones vector
            // still trips the `ptp(log2 counts) < 1e-10` guard and falls back to
            // plain eBayes (== limma), so behaviour degrades gracefully.
            deqms_two_group(
                proteins,
                rows,
                n_a,
                n_b,
                prepared.peptide_counts,
                fdr,
                log2fc,
            )
        }
        MemberMethod::Limma => limma_two_group(proteins, rows, n_a, n_b, fdr, log2fc),
    }
}

/// IHW uses 5 covariate bins, matching mokume's `_ihw_correction(..., n_bins=5)`
/// default (differential_expression.py:393).
const IHW_N_BINS: usize = 5;

/// Whether the configured FDR method is IHW (case-insensitive, trimmed),
/// matching mokume's `self.fdr_method = config["fdr_method"].lower()` and the
/// `== "ihw"` check (differential_expression.py:61, 245).
fn is_ihw(config: &DifferentialExpressionConfig) -> bool {
    config.fdr_method.trim().eq_ignore_ascii_case("ihw")
}

fn adaptive_fdr_method(config: &DifferentialExpressionConfig) -> Option<AdaptiveFdrMethod> {
    match config.fdr_method.trim().to_ascii_lowercase().as_str() {
        "bky" => Some(AdaptiveFdrMethod::Bky),
        "storey" => Some(AdaptiveFdrMethod::Storey),
        _ => None,
    }
}

/// Replace each row's BH adjusted p-value with the IHW-adjusted value. Shared
/// finalization then re-classifies and re-sorts.
///
/// The covariate is the per-protein row-mean of the two `mean_<cond>` columns
/// (`_ihw_covariate`, differential_expression.py:328-330): mokume's
/// `de_df[mean_cols].mean(axis=1)` skips NaN, so a protein observed in only one
/// group uses that group's mean; both-NaN yields a NaN covariate (excluded from
/// binning by `ihw_correction`). The fallback to `n_a + n_b`
/// (differential_expression.py:331-332) cannot trigger here because the
/// two-group methods always emit both mean columns.
fn apply_ihw(results: &mut [DeResult], config: &DifferentialExpressionConfig) {
    if results.is_empty() {
        return;
    }
    let pvalues: Vec<f64> = results.iter().map(|row| row.p_value).collect();
    let covariate: Vec<f64> = results.iter().map(ihw_covariate).collect();
    let adjusted = ihw_correction(&pvalues, &covariate, config.fdr_threshold, IHW_N_BINS);

    for (row, &adj_p_value) in results.iter_mut().zip(&adjusted) {
        row.adj_p_value = adj_p_value;
    }
}

fn apply_adaptive_fdr(
    results: &mut [DeResult],
    method: AdaptiveFdrMethod,
    config: &DifferentialExpressionConfig,
) {
    let pvalues = results.iter().map(|row| row.p_value).collect::<Vec<_>>();
    let (adjusted, applied) = adaptive_adjust(&pvalues, method, config.fdr_threshold);
    if applied == AppliedFdrMethod::Bh {
        tracing::info!(
            requested = config.fdr_method,
            "adaptive FDR was not reliable; applied BH"
        );
    }
    for (row, adj_p_value) in results.iter_mut().zip(adjusted) {
        row.adj_p_value = adj_p_value;
    }
}

fn resolved_effect_size_gate(results: &[DeResult], config: &DifferentialExpressionConfig) -> f64 {
    let Some(method) = effect_size_gate_method(config) else {
        return config.log2fc_threshold;
    };
    let fold_changes = results
        .iter()
        .map(|row| row.log2_fold_change)
        .collect::<Vec<_>>();
    estimate_effect_size_gate(&fold_changes, method, config.log2fc_threshold)
}

fn effect_size_gate_method(config: &DifferentialExpressionConfig) -> Option<EffectSizeGateMethod> {
    match config.effect_size_gate.as_deref().map(str::trim) {
        Some(method) if method.eq_ignore_ascii_case("mixture") => {
            Some(EffectSizeGateMethod::Mixture)
        }
        Some(method) if method.eq_ignore_ascii_case("null_quantile") => {
            Some(EffectSizeGateMethod::NullQuantile)
        }
        _ => None,
    }
}

fn classify_and_sort(
    results: &mut [DeResult],
    log2fc_threshold: f64,
    config: &DifferentialExpressionConfig,
) {
    for row in results.iter_mut() {
        row.significance = classify_significance(
            row.adj_p_value,
            row.log2_fold_change,
            config.fdr_threshold,
            log2fc_threshold,
        );
    }
    results.sort_by(|left, right| left.adj_p_value.total_cmp(&right.adj_p_value));
}

/// Per-protein IHW covariate: the mean of the two group means, skipping NaN
/// (pandas `mean(axis=1)` semantics). Returns NaN only when both group means are
/// NaN; a single finite mean is returned as-is.
fn ihw_covariate(row: &DeResult) -> f64 {
    match (row.mean_a.is_finite(), row.mean_b.is_finite()) {
        (true, true) => (row.mean_a + row.mean_b) / 2.0,
        (true, false) => row.mean_a,
        (false, true) => row.mean_b,
        (false, false) => f64::NAN,
    }
}

/// Significance call mirroring mokume's `_finalize_results`: a non-finite
/// adjusted p-value or log2FC is NotTested, else UP when adj-p < FDR and log2FC
/// > +threshold, DOWN for the symmetric negative threshold, or Unchanged.
fn classify_significance(
    adj_p_value: f64,
    log2_fold_change: f64,
    fdr_threshold: f64,
    log2fc_threshold: f64,
) -> Significance {
    if !adj_p_value.is_finite() || !log2_fold_change.is_finite() {
        return Significance::NotTested;
    }
    if adj_p_value < fdr_threshold && log2_fold_change > log2fc_threshold {
        Significance::Up
    } else if adj_p_value < fdr_threshold && log2_fold_change < -log2fc_threshold {
        Significance::Down
    } else {
        Significance::Unchanged
    }
}

/// Run the multi-method ensemble on the prepared contrast data, porting
/// `run_ensemble` (ensemble.py:23): run each configured member method, keep only
/// the members that produced a non-empty result (Python's `if not result.empty`),
/// then fuse with `combine_de_results`. The member list defaults to
/// `[limrots, deqms, proda]` (the Python default) and is overridable via
/// `--de-ensemble-methods`. Invalid member configurations are rejected before
/// any member runs.
fn run_ensemble(
    prepared: &Prepared<'_>,
    config: &DifferentialExpressionConfig,
) -> Result<Vec<EnsembleResult>> {
    let mut members: Vec<(String, Vec<DeResult>)> = Vec::new();
    for name in validated_ensemble_member_names(config)? {
        let method = parse_ensemble_method(&name)?;
        let result = run_member(method, prepared, config);
        // Python only adds members whose result is non-empty.
        if !result.is_empty() {
            members.push((name, result));
        }
    }
    // `run_ensemble` returns an empty frame when no member produced results;
    // `combine_de_results` over an empty member set returns no rows.
    let combined = combine_de_results(
        &members,
        config.ensemble_min_k,
        config.log2fc_threshold,
        config.fdr_threshold,
    );
    Ok(finalize_ensemble_results(combined, config))
}

fn finalize_ensemble_results(
    mut combined: Vec<EnsembleResult>,
    config: &DifferentialExpressionConfig,
) -> Vec<EnsembleResult> {
    if let Some(method) = adaptive_fdr_method(config) {
        let pvalues = combined.iter().map(|row| row.p_value).collect::<Vec<_>>();
        let (adjusted, applied) = adaptive_adjust(&pvalues, method, config.fdr_threshold);
        if applied == AppliedFdrMethod::Bh {
            tracing::info!(
                requested = config.fdr_method,
                "ensemble adaptive FDR was not reliable; applied BH"
            );
        }
        for (row, adj_p_value) in combined.iter_mut().zip(adjusted) {
            row.adj_p_value = adj_p_value;
        }
    }
    let gate = resolved_ensemble_effect_size_gate(&combined, config);
    for row in &mut combined {
        row.significance =
            classify_ensemble(row, config.ensemble_min_k, gate, config.fdr_threshold);
    }
    combined.sort_by(|left, right| left.adj_p_value.total_cmp(&right.adj_p_value));
    combined
}

fn resolved_ensemble_effect_size_gate(
    results: &[EnsembleResult],
    config: &DifferentialExpressionConfig,
) -> f64 {
    let Some(method) = effect_size_gate_method(config) else {
        return config.log2fc_threshold;
    };
    let fold_changes = results
        .iter()
        .map(|row| row.log2_fold_change)
        .collect::<Vec<_>>();
    estimate_effect_size_gate(&fold_changes, method, config.log2fc_threshold)
}

fn classify_ensemble(
    row: &EnsembleResult,
    min_k: usize,
    log2fc_threshold: f64,
    fdr_threshold: f64,
) -> Significance {
    if !row.adj_p_value.is_finite() || !row.log2_fold_change.is_finite() {
        Significance::NotTested
    } else if row.n_methods_up >= min_k
        && row.adj_p_value < fdr_threshold
        && row.log2_fold_change > log2fc_threshold
    {
        Significance::Up
    } else if row.n_methods_down >= min_k
        && row.adj_p_value < fdr_threshold
        && row.log2_fold_change < -log2fc_threshold
    {
        Significance::Down
    } else {
        Significance::Unchanged
    }
}

pub(crate) fn validate_config(
    config: &DifferentialExpressionConfig,
    allow_auto: bool,
) -> Result<()> {
    let method = config.method.trim().to_ascii_lowercase();
    let is_ensemble = method == "ensemble";
    match method.as_str() {
        "limma" | "deqms" | "rots" | "limrots" | "proda" | "ensemble" => {}
        "auto" if allow_auto => {}
        _ => {
            return Err(invalid_input(format!(
                "unknown DE method `{}`",
                config.method
            )))
        }
    }
    if !matches!(
        config.fdr_method.trim().to_ascii_lowercase().as_str(),
        "bh" | "ihw" | "bky" | "storey"
    ) {
        return Err(MokumeError::NotImplemented {
            stage: "differential-expression-fdr-method",
        });
    }
    if config.ensemble_methods.is_some() && !is_ensemble {
        return Err(invalid_input(
            "--de-ensemble-methods only applies to --de-method ensemble",
        ));
    }
    if is_ensemble {
        validated_ensemble_member_names(config)?;
    }
    validate_de_thresholds(config)?;
    validate_effect_size_gate(config)
}

fn validate_de_thresholds(config: &DifferentialExpressionConfig) -> Result<()> {
    if !config.log2fc_threshold.is_finite() || config.log2fc_threshold < 0.0 {
        return Err(invalid_input(
            "--de-log2fc must be a finite, non-negative value",
        ));
    }
    if !config.fdr_threshold.is_finite() || !(0.0..=1.0).contains(&config.fdr_threshold) {
        return Err(invalid_input("--de-fdr must be between 0 and 1"));
    }
    Ok(())
}

fn validate_effect_size_gate(config: &DifferentialExpressionConfig) -> Result<()> {
    if let Some(gate) = config.effect_size_gate.as_deref() {
        if !matches!(
            gate.trim().to_ascii_lowercase().as_str(),
            "mixture" | "null_quantile"
        ) {
            return Err(invalid_input(format!(
                "unknown effect-size gate method `{gate}`"
            )));
        }
    }
    Ok(())
}

/// Return canonical ensemble members after validating the member and consensus
/// configuration. An omitted list uses the documented default; an explicitly
/// empty list is invalid.
pub(crate) fn validated_ensemble_member_names(
    config: &DifferentialExpressionConfig,
) -> Result<Vec<String>> {
    let default_methods = ["limrots", "deqms", "proda"];
    let methods: Vec<&str> = match &config.ensemble_methods {
        None => default_methods.to_vec(),
        Some(methods) if methods.is_empty() => {
            return Err(invalid_input("ensemble member list must not be empty"));
        }
        Some(methods) => methods.iter().map(String::as_str).collect(),
    };

    let mut seen = HashSet::new();
    let mut canonical = Vec::with_capacity(methods.len());
    for method in methods {
        let name = method.trim().to_ascii_lowercase();
        if name.is_empty() {
            return Err(invalid_input("ensemble member name must not be empty"));
        }
        if name == "ensemble" {
            return Err(invalid_input("nested ensemble members are not supported"));
        }
        parse_ensemble_method(&name)?;
        if !seen.insert(name.clone()) {
            return Err(invalid_input(format!("duplicate ensemble member `{name}`")));
        }
        canonical.push(name);
    }

    if config.ensemble_min_k == 0 || config.ensemble_min_k > canonical.len() {
        return Err(invalid_input(format!(
            "ensemble min-k must be between 1 and the number of member methods ({}); got {}",
            canonical.len(),
            config.ensemble_min_k
        )));
    }
    Ok(canonical)
}

/// Build the DEqMS per-protein count vector aligned to `proteins` (the contrast
/// row order), defaulting to `1` for any protein absent from the matrix's
/// per-protein count map. This is the Rust port of Python's `_build_count_vector`
/// (deqms.py:147): `counts = np.ones(...); counts[valid] = aligned[valid]`, where
/// `aligned = peptide_counts.reindex(proteins)`.
///
/// Count semantics vs Python (verified by cross-language pairing on a synthetic
/// 15-protein, 3-vs-3 dataset):
///   * KERNEL (counts agree, single-accession proteins): both sides count the
///     distinct canonical-peptide sequences per protein, keyed by the parsed
///     protein name `proteins` carries. Here `peptide_counts_by_name` resolves
///     the same key the matrix uses, so the count vector matches Python.
///   * Niche 1 -- accession format. Python's `_load_de_peptide_counts`
///     (stages.py:1810) groups the RAW parquet by the un-parsed `anchor_protein`
///     column (`sp|P1000|NAME`), but the matrix `ProteinName` is run through
///     `parse_uniprot_accession` (aggregation.py:207) to `P1000`. With the
///     standard `db|ACC|NAME` UniProt format the reindex therefore matches
///     nothing, every count collapses to the all-ones default, and Python's
///     spectraCounteBayes degrades to plain eBayes ("all counts identical").
///     This Rust port keys both the count map and `proteins` by the same parsed
///     name, so it does NOT degrade -- it runs real count-aware moderation.
///   * Niche 2 -- protein groups. For a multi-accession group the matrix name is
///     the `;`-joined parsed group (`P1000;P1001`), while Python's
///     `anchor_protein` is a single accession, so the reindex misses and Python
///     keeps the all-ones default for that group; this port carries the real
///     group-level count.
///   * Niche 3 -- data pipeline. Python counts distinct sequences from the raw
///     parquet (pre-filter), this port from the post-filter accepted stream
///     (`unique_peptides`, lib.rs:`protein_peptide_counts`). A shared / dropped
///     peptide can shift a single-accession count by one.
///
/// In all three niches only the DEqMS-specific `peptide_count` / `sca_t` /
/// `sca_pvalue` columns differ; the limma layer (`log2FC` / `pvalue`) and the
/// protein quantification matrix remain cell-for-cell identical.
fn build_count_vector(matrix: &ProteinMatrix, proteins: &[String]) -> Vec<f64> {
    let counts_by_name = matrix.peptide_counts_by_name();
    proteins
        .iter()
        .map(|protein| {
            counts_by_name
                .get(protein)
                .map_or(1.0, |count| *count as f64)
        })
        .collect()
}

/// The DE methods the pipeline runs. Validation upstream guarantees the
/// configured `--de-method` is one of these by the time a contrast runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeMethod {
    Member(MemberMethod),
    Ensemble,
}

/// A concrete method that can safely execute as one ensemble member.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MemberMethod {
    Limma,
    Deqms,
    Rots,
    Limrots,
    Proda,
}

/// Resolve the top-level `--de-method` string (already validated to one of the
/// supported methods, including `ensemble`) to a [`DeMethod`].
fn de_method(config: &DifferentialExpressionConfig) -> Result<DeMethod> {
    parse_de_method(&config.method)
}

/// Map a validated method name to [`DeMethod`] without a scientific fallback.
fn parse_de_method(name: &str) -> Result<DeMethod> {
    let method = match name.trim().to_ascii_lowercase().as_str() {
        "limma" => DeMethod::Member(MemberMethod::Limma),
        "deqms" => DeMethod::Member(MemberMethod::Deqms),
        "rots" => DeMethod::Member(MemberMethod::Rots),
        "limrots" => DeMethod::Member(MemberMethod::Limrots),
        "proda" => DeMethod::Member(MemberMethod::Proda),
        "ensemble" => DeMethod::Ensemble,
        _ => return Err(invalid_input(format!("unknown DE method `{name}`"))),
    };
    Ok(method)
}

/// Parse one concrete ensemble member, rejecting nested ensembles.
fn parse_ensemble_method(name: &str) -> Result<MemberMethod> {
    match parse_de_method(name) {
        Err(_) => Err(invalid_input(format!(
            "unknown ensemble member `{name}`; supported members: limma, deqms, rots, limrots, proda"
        ))),
        Ok(DeMethod::Ensemble) => {
            Err(invalid_input("nested ensemble members are not supported"))
        }
        Ok(DeMethod::Member(method)) => Ok(method),
    }
}

/// Sample ids whose condition equals `condition`, in matrix column order.
fn group_samples(
    samples: &[(SampleId, &str)],
    condition_by_sample: &HashMap<String, String>,
    condition: &str,
) -> Vec<SampleId> {
    samples
        .iter()
        .filter(|(_, name)| condition_by_sample.get(*name).map(String::as_str) == Some(condition))
        .map(|(sample, _)| *sample)
        .collect()
}

/// Render the sorted unique condition labels from the SDRF-derived map as a
/// Python-style list literal (`['a', 'b']`) for the "No samples" error, mirroring
/// `_split_samples`'s `sorted(set(sample_to_condition.values()))`. Listing every
/// mapped condition (not just those present in the matrix columns) lets a user
/// whose contrast label does not match see the full set of valid labels.
fn available_conditions(condition_by_sample: &HashMap<String, String>) -> String {
    let mut conditions = condition_by_sample.values().cloned().collect::<Vec<_>>();
    conditions.sort();
    conditions.dedup();
    let inner = conditions
        .iter()
        .map(|condition| format!("'{condition}'"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{inner}]")
}

/// Parse the configured contrasts into condition pairs, matching the Python
/// `_parse_de_contrasts`: split on " vs " first, else on the first "-";
/// unparseable entries are skipped with no error.
fn parse_contrasts(config: &DifferentialExpressionConfig) -> Result<Vec<Contrast>> {
    let Some(raw) = &config.contrasts else {
        return Err(invalid_input(
            "differential expression requires explicit contrasts via --de-contrasts",
        ));
    };
    let mut contrasts = Vec::new();
    for text in raw {
        let Some((left, right)) = split_contrast(text) else {
            continue;
        };
        let cond_a = left.trim().to_owned();
        let cond_b = right.trim().to_owned();
        if cond_a.is_empty() || cond_b.is_empty() {
            continue;
        }
        contrasts.push(Contrast { cond_a, cond_b });
    }
    Ok(contrasts)
}

/// Split a single contrast string into its two condition labels. " vs " takes
/// precedence over "-", matching the Python ordering.
fn split_contrast(text: &str) -> Option<(&str, &str)> {
    if let Some(index) = text.find(" vs ") {
        let (left, right) = text.split_at(index);
        return Some((left, &right[" vs ".len()..]));
    }
    text.split_once('-')
}

/// The output path for one contrast, or `None` when no `--de-output` is set.
/// One contrast writes to `--de-output` directly; multiple contrasts append
/// the contrast key (`base_<key>.ext`), matching Python's `_de_output_path`.
fn de_output_path(
    config: &DifferentialExpressionConfig,
    contrast: &Contrast,
    n_contrasts: usize,
) -> Option<PathBuf> {
    let output = config.output.as_ref()?;
    if n_contrasts == 1 {
        return Some(output.clone());
    }

    let key = contrast.key();
    let parent = output.parent();
    let stem = output.file_stem().and_then(|stem| stem.to_str());
    let ext = output.extension().and_then(|ext| ext.to_str());
    let file_name = match (stem, ext) {
        (Some(stem), Some(ext)) => format!("{stem}_{key}.{ext}"),
        (Some(stem), None) => format!("{stem}_{key}.csv"),
        // No usable stem: fall back to a key-named file with a .csv suffix.
        (None, _) => format!("{key}.csv"),
    };
    Some(match parent {
        Some(parent) if !parent.as_os_str().is_empty() => parent.join(file_name),
        _ => PathBuf::from(file_name),
    })
}

/// Derive a sample -> condition map from the SDRF, reproducing
/// `mokume.normalization.irs.detect_condition_from_sdrf`: the first
/// `factor value[*]` column, reduced by `_shorten_factor_label`, keyed by both
/// the `source name` and the `comment[data file]` name. The data-file keys let
/// run-level protein matrices (whose columns are run-file stems, not source
/// names) resolve their conditions.
fn condition_by_sample(sdrf: &SdrfTable) -> HashMap<String, String> {
    let mut map = HashMap::new();
    // Pass 1: source name -> condition (Python's dedup'd source-name loop, first
    // value wins).
    for record in sdrf.records() {
        let Some(condition) = record.condition.as_deref() else {
            continue;
        };
        map.entry(record.sample_accession.clone())
            .or_insert_with(|| shorten_factor_label(condition));
    }
    // Pass 2: data-file name + extension-stripped stem (+ uppercase variants),
    // filling only keys not already mapped so source-name entries win
    // (irs.py:393-401). The stem strips `.raw`/`.mzML`/`.d`/`.wiff`
    // case-insensitively so the key matches the extension-less `run_file_name`
    // the protein matrix uses as its sample columns.
    for record in sdrf.records() {
        let Some(condition) = record.condition.as_deref() else {
            continue;
        };
        let raw = record.data_file.as_str();
        if raw.is_empty() {
            continue;
        }
        let short = shorten_factor_label(condition);
        let stem = strip_run_extension(raw);
        for variant in [
            raw.to_owned(),
            stem.clone(),
            raw.to_uppercase(),
            stem.to_uppercase(),
        ] {
            if !variant.is_empty() {
                map.entry(variant).or_insert_with(|| short.clone());
            }
        }
    }
    map
}

/// Strip a single trailing run extension (case-insensitive). Beyond Python's
/// `re.sub(r"\.(raw|mzML|d)$", ...)` this also strips `.wiff`: the protein-matrix
/// sample columns come from the QPX `run_file_name`, which is already
/// extension-less, so the SDRF `comment[data file]` key must drop the run
/// extension to match it. Without `.wiff`, DE on Sciex (`.wiff`) datasets finds
/// no samples for any contrast even though the conditions are present.
fn strip_run_extension(name: &str) -> String {
    let lower = name.to_ascii_lowercase();
    for ext in [".raw", ".mzml", ".d", ".wiff"] {
        if lower.ends_with(ext) {
            if let Some(stem) = name.get(..name.len() - ext.len()) {
                return stem.to_owned();
            }
        }
    }
    name.to_owned()
}

/// Port of `mokume.normalization.irs._shorten_factor_label`. Verbose
/// ontology-annotated labels like `CT=mixture;QY=12500 amol;CN=UPS1` are
/// reduced to the most informative token (`QY`, then `CN`, then `NT`, else the
/// first value); plain labels pass through unchanged.
fn shorten_factor_label(raw: &str) -> String {
    if !raw.contains(';') || !raw.contains('=') {
        return raw.to_owned();
    }
    let mut parts: Vec<(String, String)> = Vec::new();
    for token in raw.split(';') {
        if let Some((key, value)) = token.split_once('=') {
            parts.push((key.trim().to_owned(), value.trim().to_owned()));
        }
    }
    for key in ["QY", "CN", "NT"] {
        if let Some((_, value)) = parts.iter().find(|(name, _)| name == key) {
            return value.clone();
        }
    }
    parts
        .first()
        .map(|(_, value)| value.clone())
        .unwrap_or_else(|| raw.to_owned())
}

/// Write a DE result table, reproducing the per-method column order of
/// `DifferentialExpression(method=...).run()`.
///
/// limma emits the moderated-fit extras (`extra_cols=("t_stat","AveExpr","B")`):
/// `ProteinName, log2FC, pvalue, adj_pvalue, t_stat, AveExpr, B,
///  mean_<cond_a>, mean_<cond_b>, n_a, n_b, significance`.
///
/// deqms emits the spectraCounteBayes columns (deqms.py:219-232) plus the
/// `significance` column `_finalize_results` appends:
/// `ProteinName, log2FC, pvalue, adj_pvalue, sca_t, sca_pvalue, sca_adj_pvalue,
///  mean_<cond_a>, mean_<cond_b>, n_a, n_b, peptide_count, log_pvalue,
///  significance`. Here
/// `pvalue == sca_pvalue` and `adj_pvalue == sca_adj_pvalue` (both BH of
/// `sca_pvalue`), exactly as Python sets them, and `peptide_count` is the
/// per-protein unique-canonical-peptide count captured at ingest, aligned by
/// protein name (default 1 for any protein with no recorded count -- see
/// `build_count_vector`).
fn write_de_csv(
    path: &Path,
    results: &[DeResult],
    contrast: &Contrast,
    method: MemberMethod,
    matrix: &ProteinMatrix,
) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path).map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);

    let header = match method {
        MemberMethod::Limma => format!(
            "ProteinName,log2FC,pvalue,adj_pvalue,t_stat,AveExpr,B,mean_{},mean_{},n_a,n_b,significance",
            contrast.cond_a, contrast.cond_b
        ),
        MemberMethod::Deqms => format!(
            "ProteinName,log2FC,pvalue,adj_pvalue,sca_t,sca_pvalue,sca_adj_pvalue,mean_{},mean_{},n_a,n_b,peptide_count,log_pvalue,significance",
            contrast.cond_a, contrast.cond_b
        ),
        // rots' finalize_de_result uses extra_cols=("d_stat",) (rots.py:331), so
        // a single d_stat column sits between adj_pvalue and the mean columns.
        MemberMethod::Rots => format!(
            "ProteinName,log2FC,pvalue,adj_pvalue,d_stat,mean_{},mean_{},n_a,n_b,significance",
            contrast.cond_a, contrast.cond_b
        ),
        // limrots emits a single `t_stat` column (== the LimROTS d-statistic)
        // between adj_pvalue and the mean columns (limrots.py:319-329).
        MemberMethod::Limrots => format!(
            "ProteinName,log2FC,pvalue,adj_pvalue,t_stat,mean_{},mean_{},n_a,n_b,significance",
            contrast.cond_a, contrast.cond_b
        ),
        // proda emits the Wald `t_stat` as its single extra column
        // (proda.py:700-714 lists `t_stat` right after adj_pvalue).
        MemberMethod::Proda => format!(
            "ProteinName,log2FC,pvalue,adj_pvalue,t_stat,mean_{},mean_{},n_a,n_b,significance",
            contrast.cond_a, contrast.cond_b
        ),
    };
    write_line(&mut writer, path, &header)?;

    // deqms emits a per-protein `peptide_count` column; resolve the same aligned
    // count map `build_count_vector` uses (default 1 for absent proteins) so the
    // written value matches Python's `counts` column exactly. Other methods never
    // read this map.
    let counts_by_name = match method {
        MemberMethod::Deqms => matrix.peptide_counts_by_name(),
        _ => HashMap::new(),
    };

    for row in results {
        let prefix = format!(
            "{},{},{},{}",
            row.protein,
            format_float(row.log2_fold_change),
            format_float(row.p_value),
            format_float(row.adj_p_value),
        );
        let extras = match method {
            MemberMethod::Limma => format!(
                "{},{},{},",
                format_float(row.t_statistic),
                format_float(row.ave_expr),
                format_float(row.b),
            ),
            // sca_t, sca_pvalue (== pvalue), sca_adj_pvalue (== adj_pvalue).
            MemberMethod::Deqms => format!(
                "{},{},{},",
                format_float(row.t_statistic),
                format_float(row.p_value),
                format_float(row.adj_p_value),
            ),
            // rots' single d_stat extra column (carried in t_statistic).
            MemberMethod::Rots => format!("{},", format_float(row.t_statistic)),
            // limrots' single t_stat extra column (== d_stat, carried in
            // t_statistic).
            MemberMethod::Limrots => format!("{},", format_float(row.t_statistic)),
            // proda's single Wald t_stat extra column (carried in t_statistic).
            MemberMethod::Proda => format!("{},", format_float(row.t_statistic)),
        };
        // deqms appends a `peptide_count` column before `significance`; it carries
        // the per-protein unique-peptide count aligned by name (default 1 for any
        // protein with no recorded count), matching Python's `_build_count_vector`.
        let suffix = match method {
            MemberMethod::Deqms => {
                format!(
                    ",{},{}",
                    counts_by_name.get(&row.protein).copied().unwrap_or(1),
                    format_float(row.log_p_value),
                )
            }
            _ => String::new(),
        };
        let line = format!(
            "{prefix},{extras}{},{},{},{}{suffix},{}",
            format_float(row.mean_a),
            format_float(row.mean_b),
            row.n_a,
            row.n_b,
            significance_label(row.significance),
        );
        write_line(&mut writer, path, &line)?;
    }
    writer.flush().map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })
}

/// Write the ensemble DE table, reproducing the column set and order of
/// `combine_de_results` (ensemble.py:158-199). The Python frame builds the
/// `rows` dict first (`Protein, log2FC, pvalue, n_methods_up, n_methods_down,
/// methods_significant`), then assigns `adj_pvalue` and `significance` as new
/// columns, so the emitted order is:
/// `Protein, log2FC, pvalue, n_methods_up, n_methods_down, methods_significant,
///  adj_pvalue, significance`. Rows arrive already sorted by `adj_pvalue`.
///
/// `methods_significant` joins member names with commas and so is CSV-quoted
/// (pandas `QUOTE_MINIMAL`) via [`csv_quote`]; an empty value (no member called
/// the protein UP/DOWN) needs no quoting.
fn write_ensemble_csv(path: &Path, results: &[EnsembleResult], contrast: &Contrast) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path).map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);

    // The ensemble column set is contrast-independent (no per-group mean columns);
    // `contrast` is kept in the signature for symmetry with `write_de_csv` and in
    // case a future port adds the per-condition means.
    let _ = contrast;
    write_line(
        &mut writer,
        path,
        "Protein,log2FC,pvalue,n_methods_up,n_methods_down,methods_significant,adj_pvalue,significance",
    )?;

    for row in results {
        let line = format!(
            "{},{},{},{},{},{},{},{}",
            row.protein,
            format_float(row.log2_fold_change),
            format_float(row.p_value),
            row.n_methods_up,
            row.n_methods_down,
            csv_quote(&row.methods_significant),
            format_float(row.adj_p_value),
            significance_label(row.significance),
        );
        write_line(&mut writer, path, &line)?;
    }
    writer.flush().map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })
}

fn write_line(writer: &mut impl Write, path: &Path, line: &str) -> Result<()> {
    writeln!(writer, "{line}").map_err(|source| MokumeError::Io {
        path: path.to_path_buf(),
        source,
    })
}

/// Significance label matching pandas' `significance` strings.
fn significance_label(significance: Significance) -> &'static str {
    match significance {
        Significance::Up => "UP",
        Significance::Down => "DOWN",
        Significance::Unchanged => "Unchanged",
        Significance::NotTested => "NotTested",
    }
}

/// Format a float for the DE CSV. Non-finite values (e.g. a NaN log2FC coerced
/// by the P20 guard) are written as an empty cell, matching how pandas renders
/// NaN to CSV.
fn format_float(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        String::new()
    }
}

/// Minimal CSV quoting for a text cell, matching pandas' default
/// `QUOTE_MINIMAL`: a field containing a comma, a double quote, or a newline is
/// wrapped in double quotes with embedded quotes doubled. The ensemble's
/// `methods_significant` column joins member names with commas, so it must be
/// quoted to stay a single valid CSV field.
fn csv_quote(field: &str) -> String {
    if field.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_owned()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};
    use std::error::Error;
    use std::time::{SystemTime, UNIX_EPOCH};

    use mokume_core::{DifferentialExpressionConfig, ProteinId, SampleId};
    use mokume_io::SdrfTable;
    use mokume_stats::de::{DeResult, EnsembleResult, Significance};

    use super::{
        available_conditions, build_count_vector, classify_significance, condition_by_sample,
        differential_expression_matrix, finalize_ensemble_results, finalize_member_results,
        run_differential_expression, shorten_factor_label, split_contrast,
        validated_ensemble_member_names, MemberMethod,
    };
    use crate::{CellKey, ProteinMatrix, ProteinValues};

    // `super`/`crate` pull in the crate's single-parameter `Result` alias, so
    // the test helpers spell out the two-parameter standard `Result`.
    type TestResult<T> = std::result::Result<T, Box<dyn Error>>;

    // Sample column headers (= SDRF source names) in matrix order; A* are
    // groupA, B* are groupB.
    const SAMPLES: [&str; 6] = ["A1", "A2", "A3", "B1", "B2", "B3"];
    const PROTEINS: [&str; 6] = ["P1", "P2", "P3", "P4", "P5", "P6"];

    // Raw (linear) intensities, proteins x samples, identical to the Python
    // oracle fixture in this module's golden test comment below.
    const RAW: [[f64; 6]; 6] = [
        [1000.0, 1149.0, 870.0, 4000.0, 4290.0, 3732.0],
        [32000.0, 34000.0, 30000.0, 8000.0, 9200.0, 6800.0],
        [256.0, 270.0, 244.0, 257.0, 248.0, 256.0],
        [
            1048576.0, 4194304.0, 262144.0, 2097152.0, 524288.0, 8388608.0,
        ],
        [32.0, 32.2, 31.8, 64.0, 64.4, 63.6],
        [512.0, 724.0, 362.0, 480.0, 660.0, 350.0],
    ];

    // Golden table captured from the Python oracle (see exact command below).
    // Tuple: (protein, log2FC, pvalue, adj_pvalue, t_stat, significance).
    const EXPECTED: &[(&str, f64, f64, f64, f64, &str)] = &[
        (
            "P1",
            -2.0004868432854863,
            3.0352542907348282e-05,
            7.47541102122209e-05,
            -16.39236762936801,
            "DOWN",
        ),
        (
            "P2",
            2.0090616097754097,
            3.737705510611045e-05,
            7.47541102122209e-05,
            15.653604075943404,
            "UP",
        ),
        (
            "P3",
            0.015910691677655464,
            0.7393468894849923,
            0.8128513866443309,
            0.353558782338873,
            "Unchanged",
        ),
        (
            "P4",
            -1.0,
            0.5436591577043687,
            0.8128513866443309,
            -0.6554181358332211,
            "Unchanged",
        ),
        (
            "P5",
            -1.0,
            8.847368379958774e-08,
            5.308421027975264e-07,
            -58.99007257281155,
            "DOWN",
        ),
        (
            "P6",
            0.09175595082658106,
            0.8128513866443309,
            0.8128513866443309,
            0.25075791047305773,
            "Unchanged",
        ),
    ];

    // limma's AveExpr is the row mean of the log2 expression matrix over ALL
    // samples (limma's Amean), independent of the contrast -- NOT the log2FC.
    // Computed as mean(log2(RAW[row])); deterministic and cell-exact.
    const EXPECTED_AVE_EXPR: &[(&str, f64)] = &[
        ("P1", 10.96584974099084),
        ("P2", 13.959371292060668),
        ("P3", 7.994562299032412),
        ("P4", 20.5),
        ("P5", 5.499981214541417),
        ("P6", 8.954019282642177),
    ];

    fn temp_dir(tag: &str) -> TestResult<std::path::PathBuf> {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        Ok(tempfile::Builder::new()
            .prefix(&format!("mokume-de-{tag}-{nanos}-"))
            .tempdir()?
            .keep())
    }

    /// Build a `ProteinMatrix` directly from the raw fixture, with every protein
    /// kept and every sample present.
    fn fixture_matrix() -> TestResult<ProteinMatrix> {
        let mut proteins = mokume_core::StringIdRegistry::<ProteinId>::new();
        let mut samples = mokume_core::StringIdRegistry::<SampleId>::new();
        let mut sample_ids = Vec::with_capacity(SAMPLES.len());
        for name in SAMPLES {
            sample_ids.push(samples.get_or_insert(name).ok_or("sample id overflow")?);
        }
        let mut values = HashMap::new();
        let mut allowed = HashSet::new();
        for (row, name) in PROTEINS.iter().enumerate() {
            let protein = proteins.get_or_insert(name).ok_or("protein id overflow")?;
            allowed.insert(protein);
            for (col, sample) in sample_ids.iter().enumerate() {
                values.insert(
                    CellKey {
                        protein,
                        sample: *sample,
                    },
                    RAW[row][col],
                );
            }
        }
        Ok(ProteinMatrix {
            proteins,
            samples,
            allowed_proteins: allowed,
            excluded_samples: HashSet::new(),
            // No recorded peptide counts: the deqms path then defaults every
            // protein to 1 (the all-ones `_build_count_vector` fallback), which is
            // what this 6-protein fixture exercises (< 10 valid points -> eBayes).
            peptide_counts: HashMap::new(),
            values: ProteinValues::Cells(values),
        })
    }

    // KERNEL pin for the DEqMS count vector: when a protein name is present in
    // the matrix's per-protein unique-canonical count map, `build_count_vector`
    // carries that real count, aligned by the parsed name the contrast row order
    // uses. A name absent from the map keeps the all-ones default -- the same
    // fallback Python reaches when its `anchor_protein`-keyed reindex misses (the
    // niche documented on `build_count_vector`). Single-accession proteins with a
    // recorded count are the cell-for-cell-consistent KERNEL; the all-ones
    // fallback is the niche boundary.
    #[test]
    fn deqms_count_vector_uses_real_count_and_defaults_missing_to_one() -> TestResult<()> {
        let mut proteins = mokume_core::StringIdRegistry::<ProteinId>::new();
        let p1 = proteins.get_or_insert("P1").ok_or("protein id overflow")?;
        let p2 = proteins.get_or_insert("P2").ok_or("protein id overflow")?;
        // P3 is a contrast row with no recorded count (e.g. a `;`-joined group or
        // a protein Python's anchor reindex would miss): it must default to 1.
        let _p3 = proteins.get_or_insert("P3").ok_or("protein id overflow")?;

        let mut peptide_counts = HashMap::new();
        peptide_counts.insert(p1, 4usize);
        peptide_counts.insert(p2, 2usize);

        let matrix = ProteinMatrix {
            proteins,
            samples: mokume_core::StringIdRegistry::<SampleId>::new(),
            allowed_proteins: HashSet::new(),
            excluded_samples: HashSet::new(),
            peptide_counts,
            values: ProteinValues::Cells(HashMap::new()),
        };

        let order = vec!["P1".to_owned(), "P2".to_owned(), "P3".to_owned()];
        let counts = build_count_vector(&matrix, &order);

        assert_eq!(
            counts,
            vec![4.0, 2.0, 1.0],
            "real counts for recorded proteins, all-ones default for the missing one"
        );
        Ok(())
    }

    fn fixture_sdrf() -> TestResult<SdrfTable> {
        let mut tsv = String::from("source name\tcomment[data file]\tfactor value[phenotype]\n");
        for name in SAMPLES {
            let condition = if name.starts_with('A') {
                "groupA"
            } else {
                "groupB"
            };
            tsv.push_str(&format!("{name}\t{name}.raw\t{condition}\n"));
        }
        Ok(SdrfTable::from_reader(tsv.as_bytes())?)
    }

    fn de_config(output: &std::path::Path, method: &str) -> DifferentialExpressionConfig {
        DifferentialExpressionConfig {
            enabled: true,
            contrasts: Some(vec!["groupA vs groupB".to_string()]),
            method: method.to_string(),
            output: Some(output.to_path_buf()),
            ..DifferentialExpressionConfig::default()
        }
    }

    #[test]
    fn adaptive_fdr_and_auto_gate_are_wired_into_member_finalization() {
        let mut pvalues = (0..900)
            .map(|index| (index as f64 + 0.5) / 900.0)
            .collect::<Vec<_>>();
        pvalues.extend((0..300).map(|index| 0.05 * ((index as f64 + 0.5) / 300.0).powi(3)));
        let rows = pvalues
            .into_iter()
            .enumerate()
            .map(|(index, p_value)| DeResult {
                protein: format!("P{index:04}"),
                log2_fold_change: if index < 900 {
                    -0.02 + 0.04 * index as f64 / 899.0
                } else {
                    0.25 + 0.10 * (index - 900) as f64 / 299.0
                },
                p_value,
                log_p_value: p_value.ln(),
                adj_p_value: p_value,
                t_statistic: 0.0,
                ave_expr: 0.0,
                b: 0.0,
                mean_a: 0.0,
                mean_b: 0.0,
                n_a: 3,
                n_b: 3,
                significance: Significance::Unchanged,
            })
            .collect::<Vec<_>>();
        let storey_config = DifferentialExpressionConfig {
            fdr_method: "storey".to_string(),
            effect_size_gate: Some("mixture".to_string()),
            ..DifferentialExpressionConfig::default()
        };

        let bky_config = DifferentialExpressionConfig {
            fdr_method: "bky".to_string(),
            effect_size_gate: Some("mixture".to_string()),
            ..DifferentialExpressionConfig::default()
        };
        let bky = finalize_member_results(MemberMethod::Limma, rows.clone(), &bky_config);
        assert_eq!(bky.iter().filter(|row| row.adj_p_value < 0.05).count(), 168);

        let finalized = finalize_member_results(MemberMethod::Limma, rows, &storey_config);
        assert_eq!(
            finalized
                .iter()
                .filter(|row| row.adj_p_value < 0.05)
                .count(),
            187
        );
        assert!(finalized
            .iter()
            .any(|row| row.significance == Significance::Up));
        assert!(finalized
            .iter()
            .all(|row| { row.significance != Significance::Up || row.log2_fold_change < 0.5 }));
    }

    #[test]
    fn adaptive_fdr_and_auto_gate_reach_ensemble_combination() {
        let mut pvalues = (0..900)
            .map(|index| (index as f64 + 0.5) / 900.0)
            .collect::<Vec<_>>();
        pvalues.extend((0..300).map(|index| 0.05 * ((index as f64 + 0.5) / 300.0).powi(3)));
        let rows = pvalues
            .into_iter()
            .enumerate()
            .map(|(index, p_value)| EnsembleResult {
                protein: format!("P{index:04}"),
                log2_fold_change: if index < 900 {
                    -0.02 + 0.04 * index as f64 / 899.0
                } else {
                    0.25 + 0.10 * (index - 900) as f64 / 299.0
                },
                p_value,
                adj_p_value: p_value,
                n_methods_up: usize::from(index >= 900) * 2,
                n_methods_down: 0,
                methods_significant: String::new(),
                significance: Significance::Unchanged,
            })
            .collect::<Vec<_>>();
        let config = DifferentialExpressionConfig {
            ensemble_min_k: 2,
            fdr_method: "storey".to_string(),
            effect_size_gate: Some("mixture".to_string()),
            ..DifferentialExpressionConfig::default()
        };

        let finalized = finalize_ensemble_results(rows, &config);
        assert_eq!(
            finalized
                .iter()
                .filter(|row| row.adj_p_value < 0.05)
                .count(),
            187
        );
        assert!(finalized
            .iter()
            .any(|row| row.significance == Significance::Up));
    }

    #[test]
    fn condition_by_sample_maps_data_file_stems() -> TestResult<()> {
        // source name != run-file stem, so run-level matrix columns (the stems)
        // only resolve via the data-file pass. `.raw`/`.mzML`/`.wiff` are all
        // stripped (case-insens) so the stem matches the extension-less
        // run_file_name the protein matrix uses as columns.
        let tsv = concat!(
            "source name\tcomment[data file]\tfactor value[phenotype]\n",
            "Sample-1\trun_alpha.raw\tgroupA\n",
            "Sample-2\tRUN_BETA.mzML\tgroupB\n",
            "Sample-3\trun_gamma.wiff\tgroupA\n",
        );
        let sdrf = SdrfTable::from_reader(tsv.as_bytes())?;
        let map = condition_by_sample(&sdrf);

        // Source names still resolve.
        assert_eq!(map.get("Sample-1").map(String::as_str), Some("groupA"));
        // `.raw` stem resolves (the run-level fix).
        assert_eq!(map.get("run_alpha").map(String::as_str), Some("groupA"));
        assert_eq!(map.get("run_alpha.raw").map(String::as_str), Some("groupA"));
        // `.mzML` stem strips case-insensitively; uppercase variant present too.
        assert_eq!(map.get("RUN_BETA").map(String::as_str), Some("groupB"));
        // `.wiff` is stripped too, so the extension-less stem (the protein-matrix
        // column form) and the full name both resolve.
        assert_eq!(map.get("run_gamma").map(String::as_str), Some("groupA"));
        assert_eq!(
            map.get("run_gamma.wiff").map(String::as_str),
            Some("groupA")
        );
        Ok(())
    }

    #[test]
    fn available_conditions_lists_all_map_values_python_style() {
        // Mirrors `_split_samples`'s `sorted(set(sample_to_condition.values()))`:
        // every mapped condition appears (sorted, deduped), rendered as a
        // Python list literal with single quotes, regardless of which keys a
        // protein matrix actually carries.
        let mut map = HashMap::new();
        map.insert("Sample-1".to_string(), "groupB".to_string());
        map.insert("run_alpha".to_string(), "groupB".to_string());
        map.insert("Sample-2".to_string(), "groupA".to_string());
        map.insert("Sample-3".to_string(), "groupA".to_string());
        assert_eq!(available_conditions(&map), "['groupA', 'groupB']");

        assert_eq!(available_conditions(&HashMap::new()), "[]");
    }

    #[test]
    fn split_contrast_prefers_vs_over_dash() {
        assert_eq!(split_contrast("A vs B"), Some(("A", "B")));
        assert_eq!(split_contrast("A-B"), Some(("A", "B")));
        // " vs " wins even when dashes are present in the labels.
        assert_eq!(split_contrast("a-1 vs b-2"), Some(("a-1", "b-2")));
    }

    #[test]
    fn matrix_de_rejects_infinite_input() {
        let proteins = vec!["P1".to_owned()];
        let values = vec![vec![1.0, f64::INFINITY, 2.0, 3.0]];
        let config = DifferentialExpressionConfig {
            enabled: true,
            method: "limma".to_owned(),
            ..DifferentialExpressionConfig::default()
        };

        let Err(error) =
            differential_expression_matrix(&proteins, &values, 2, 2, None, &config, Some(2))
        else {
            panic!("Inf must be rejected");
        };
        assert!(error.to_string().contains("contains an infinite value"));
    }

    #[test]
    fn classification_marks_non_finite_results_not_tested() {
        assert_eq!(
            classify_significance(f64::NAN, 1.0, 0.05, 0.5),
            Significance::NotTested
        );
        assert_eq!(
            classify_significance(0.01, f64::NAN, 0.05, 0.5),
            Significance::NotTested
        );
    }

    #[test]
    fn shorten_factor_label_matches_python() {
        assert_eq!(shorten_factor_label("groupA"), "groupA");
        assert_eq!(
            shorten_factor_label("CT=mixture;QY=12500 amol;CN=UPS1"),
            "12500 amol"
        );
        assert_eq!(shorten_factor_label("CT=mixture;CN=UPS1"), "UPS1");
        // No ';' or no '=' -> unchanged.
        assert_eq!(shorten_factor_label("plain-label"), "plain-label");
    }

    // Golden test: the Rust DE dispatcher reproduces the Python oracle's limma
    // table for a 6-protein / 3-vs-3 design (raw intensities, log2 inside the
    // dispatcher). Oracle command:
    //   conda run -n Bigbio python  # against
    //   mokume.analysis.differential_expression
    //       .DifferentialExpression(method="limma", log2fc_threshold=0.5,
    //                               fdr_threshold=0.05)
    //   .run(protein_df, {A1..A3: groupA, B1..B3: groupB}, ("groupA","groupB"))
    // on protein_df with raw intensities matching RAW above. Assertions are by
    // protein name (the adj-p sort can tie); deterministic limma matched to
    // relative 1e-6.
    #[test]
    fn de_dispatcher_matches_python_limma_oracle() -> TestResult<()> {
        let dir = temp_dir("oracle")?;
        let output = dir.join("de_results.csv");
        let matrix = fixture_matrix()?;
        let sdrf = fixture_sdrf()?;

        run_differential_expression(&matrix, Some(&sdrf), &de_config(&output, "limma"), false)?;

        let text = std::fs::read_to_string(&output)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("missing header")?;
        assert_eq!(
            header,
            "ProteinName,log2FC,pvalue,adj_pvalue,t_stat,AveExpr,B,mean_groupA,mean_groupB,n_a,n_b,significance"
        );

        let mut rows: HashMap<String, Vec<String>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        for line in lines {
            let fields = line.split(',').map(ToOwned::to_owned).collect::<Vec<_>>();
            let protein = fields.first().ok_or("empty row")?.clone();
            order.push(protein.clone());
            rows.insert(protein, fields);
        }
        assert_eq!(rows.len(), EXPECTED.len(), "row count must match oracle");

        // The emitted order must be sorted ascending by adj_pvalue: P5 (smallest)
        // first; the two tied 7.475e-05 rows (P1, P2) keep input order.
        assert_eq!(order.first().map(String::as_str), Some("P5"));

        for &(protein, log2fc, pvalue, adj, t_stat, sig) in EXPECTED {
            let fields = rows
                .get(protein)
                .ok_or_else(|| format!("missing {protein}"))?;
            assert_close(parse(fields, 1)?, log2fc, "log2FC", protein);
            assert_close(parse(fields, 2)?, pvalue, "pvalue", protein);
            assert_close(parse(fields, 3)?, adj, "adj_pvalue", protein);
            assert_close(parse(fields, 4)?, t_stat, "t_stat", protein);
            // AveExpr is the row mean of the log2 matrix (limma's Amean), not
            // the log2FC.
            let ave_expr = EXPECTED_AVE_EXPR
                .iter()
                .find(|&&(p, _)| p == protein)
                .map(|&(_, v)| v)
                .ok_or_else(|| format!("missing AveExpr for {protein}"))?;
            assert_close(parse(fields, 5)?, ave_expr, "AveExpr", protein);
            assert_eq!(
                fields.get(6).map(String::as_str),
                Some("0"),
                "B for {protein}"
            );
            assert_eq!(
                fields.get(9).map(String::as_str),
                Some("3"),
                "n_a for {protein}"
            );
            assert_eq!(
                fields.get(10).map(String::as_str),
                Some("3"),
                "n_b for {protein}"
            );
            assert_eq!(
                fields.get(11).map(String::as_str),
                Some(sig),
                "sig for {protein}"
            );
        }
        Ok(())
    }

    // ===== IHW (--de-fdr-method ihw) dispatcher fixtures =====
    //
    // A 12-protein / 3-vs-3 fixture (>= 2 * n_bins = 10, so IHW actually runs
    // rather than falling back to BH) with a spread of effect sizes AND mean
    // intensities, so the IHW covariate (the row-mean of the two mean_ columns)
    // bins informatively. Shared by the LimROTS/ROTS "preserve own FDR under
    // ihw" regression below.
    const IHW_SAMPLES: [&str; 6] = ["A1", "A2", "A3", "B1", "B2", "B3"];
    const IHW_PROTEINS: [&str; 12] = [
        "P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12",
    ];
    const IHW_RAW: [[f64; 6]; 12] = [
        [1000.0, 1149.0, 870.0, 4000.0, 4290.0, 3732.0],
        [32000.0, 34000.0, 30000.0, 8000.0, 9200.0, 6800.0],
        [256.0, 270.0, 244.0, 257.0, 248.0, 256.0],
        [
            1048576.0, 4194304.0, 262144.0, 2097152.0, 524288.0, 8388608.0,
        ],
        [32.0, 32.2, 31.8, 64.0, 64.4, 63.6],
        [512.0, 724.0, 362.0, 480.0, 660.0, 350.0],
        [5000.0, 5100.0, 4900.0, 5050.0, 4950.0, 5000.0],
        [200.0, 205.0, 195.0, 410.0, 400.0, 420.0],
        [90000.0, 92000.0, 88000.0, 45000.0, 46000.0, 44000.0],
        [1500.0, 1520.0, 1480.0, 1505.0, 1495.0, 1500.0],
        [700.0, 690.0, 710.0, 1400.0, 1390.0, 1410.0],
        [60000.0, 61000.0, 59000.0, 60500.0, 59500.0, 60000.0],
    ];

    /// Build a 12-protein `ProteinMatrix` from `IHW_RAW` (every protein kept,
    /// every sample present), mirroring `fixture_matrix` for the wider fixture.
    fn ihw_fixture_matrix() -> TestResult<ProteinMatrix> {
        let mut proteins = mokume_core::StringIdRegistry::<ProteinId>::new();
        let mut samples = mokume_core::StringIdRegistry::<SampleId>::new();
        let mut sample_ids = Vec::with_capacity(IHW_SAMPLES.len());
        for name in IHW_SAMPLES {
            sample_ids.push(samples.get_or_insert(name).ok_or("sample id overflow")?);
        }
        let mut values = HashMap::new();
        let mut allowed = HashSet::new();
        for (row, name) in IHW_PROTEINS.iter().enumerate() {
            let protein = proteins.get_or_insert(name).ok_or("protein id overflow")?;
            allowed.insert(protein);
            for (col, sample) in sample_ids.iter().enumerate() {
                values.insert(
                    CellKey {
                        protein,
                        sample: *sample,
                    },
                    IHW_RAW[row][col],
                );
            }
        }
        Ok(ProteinMatrix {
            proteins,
            samples,
            allowed_proteins: allowed,
            excluded_samples: HashSet::new(),
            peptide_counts: HashMap::new(),
            values: ProteinValues::Cells(values),
        })
    }

    fn ihw_fixture_sdrf() -> TestResult<SdrfTable> {
        let mut tsv = String::from("source name\tcomment[data file]\tfactor value[phenotype]\n");
        for name in IHW_SAMPLES {
            let condition = if name.starts_with('A') {
                "groupA"
            } else {
                "groupB"
            };
            tsv.push_str(&format!("{name}\t{name}.raw\t{condition}\n"));
        }
        Ok(SdrfTable::from_reader(tsv.as_bytes())?)
    }

    // Regression: LimROTS and ROTS carry their own permutation-based FDR (a
    // `bh_adjust` over the permutation p-values), so another FDR request must
    // NOT overwrite it -- mirroring the Python
    // `_run_limrots`/`_run_rots` routing through `_classify_and_sort` instead of
    // `_finalize_results`. The internal PRNG is deterministic run-to-run, so the
    // permutation FDR is stable: running each method under `bh` and every
    // alternative must yield identical adj_pvalues per protein.
    #[test]
    fn de_dispatcher_limrots_rots_preserve_own_fdr() -> TestResult<()> {
        fn adj_pvalues(method: &str, fdr_method: &str) -> TestResult<HashMap<String, f64>> {
            let dir = temp_dir(&format!("preserve-fdr-{method}-{fdr_method}"))?;
            let output = dir.join("de_results.csv");
            let matrix = ihw_fixture_matrix()?;
            let sdrf = ihw_fixture_sdrf()?;
            let config = DifferentialExpressionConfig {
                enabled: true,
                contrasts: Some(vec!["groupA vs groupB".to_string()]),
                method: method.to_string(),
                fdr_method: fdr_method.to_string(),
                output: Some(output.clone()),
                ..DifferentialExpressionConfig::default()
            };
            run_differential_expression(&matrix, Some(&sdrf), &config, false)?;
            let text = std::fs::read_to_string(&output)?;
            let mut lines = text.lines();
            lines.next().ok_or("missing header")?;
            let mut adj = HashMap::new();
            for line in lines {
                let fields = line.split(',').map(ToOwned::to_owned).collect::<Vec<_>>();
                let protein = fields.first().ok_or("empty row")?.clone();
                adj.insert(protein, parse(&fields, 3)?);
            }
            Ok(adj)
        }

        for method in ["limrots", "rots"] {
            let bh = adj_pvalues(method, "bh")?;
            for fdr_method in ["ihw", "bky", "storey"] {
                let alternative = adj_pvalues(method, fdr_method)?;
                assert!(!bh.is_empty(), "{method}: produced no rows");
                assert_eq!(
                    bh.len(),
                    alternative.len(),
                    "{method}: bh/{fdr_method} row count differs"
                );
                for (protein, bh_adj) in &bh {
                    let alternative_adj = alternative
                        .get(protein)
                        .ok_or_else(|| format!("{method}: missing {protein} under {fdr_method}"))?;
                    assert_exact(
                        *alternative_adj,
                        *bh_adj,
                        "permutation_fdr_preserved",
                        protein,
                    );
                }
            }
        }
        Ok(())
    }

    // Pipeline wiring test for `--de-method deqms`. The 6-protein fixture has
    // fewer than 10 valid points, so the in-pipeline deqms (with its all-ones
    // default counts) falls back to plain eBayes and reproduces the limma oracle
    // table (EXPECTED). This checks: (1) the deqms-specific header (sca_t,
    // sca_pvalue, sca_adj_pvalue, peptide_count, log_pvalue); (2) the
    // peptide_count column is the all-ones fallback default (1); (3)
    // sca_pvalue == pvalue and sca_adj_pvalue == adj_pvalue; (4) the fallback
    // values equal the limma oracle. The spectraCounteBayes count path itself is
    // covered cell-by-cell in mokume-stats' deqms unit tests with explicit counts.
    #[test]
    fn de_dispatcher_deqms_falls_back_to_limma_oracle() -> TestResult<()> {
        let dir = temp_dir("oracle-deqms")?;
        let output = dir.join("de_results.csv");
        let matrix = fixture_matrix()?;
        let sdrf = fixture_sdrf()?;

        run_differential_expression(&matrix, Some(&sdrf), &de_config(&output, "deqms"), false)?;

        let text = std::fs::read_to_string(&output)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("missing header")?;
        assert_eq!(
            header,
            "ProteinName,log2FC,pvalue,adj_pvalue,sca_t,sca_pvalue,sca_adj_pvalue,mean_groupA,mean_groupB,n_a,n_b,peptide_count,log_pvalue,significance"
        );

        let mut rows: HashMap<String, Vec<String>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        for line in lines {
            let fields = line.split(',').map(ToOwned::to_owned).collect::<Vec<_>>();
            let protein = fields.first().ok_or("empty row")?.clone();
            order.push(protein.clone());
            rows.insert(protein, fields);
        }
        assert_eq!(rows.len(), EXPECTED.len(), "row count must match oracle");
        // Adjusted-p sort: P5 first (the fallback == limma ordering).
        assert_eq!(order.first().map(String::as_str), Some("P5"));

        for &(protein, log2fc, pvalue, adj, t_stat, sig) in EXPECTED {
            let fields = rows
                .get(protein)
                .ok_or_else(|| format!("missing {protein}"))?;
            assert_deqms_row(fields, protein, log2fc, pvalue, adj, t_stat, sig)?;
        }
        Ok(())
    }

    // ROTS log2FC is `nanmean(A) - nanmean(B)` (rots.py:313) -- the SAME sign
    // convention as limma's design-matrix coefficient and limrots -- so P1 (low
    // in A, high in B) is -2.0, matching the limma table. Captured from numpy on
    // log2(RAW): nanmean(A)-nanmean(B) per protein. Deterministic and cell-exact
    // regardless of the RNG.
    const EXPECTED_ROTS_LOG2FC: &[(&str, f64)] = &[
        ("P1", -2.0004868432854863),
        ("P2", 2.0090616097754097),
        ("P3", 0.015910691677656352),
        ("P4", -1.0),
        ("P5", -1.0),
        ("P6", 0.09175595082658106),
    ];

    // Pipeline wiring test for `--de-method rots`. ROTS is a faithful RNG-based
    // port (fixed internal PRNG), so its permutation p-values are deterministic
    // run-to-run but NOT bit-matched to Python (Python is itself seed-unstable).
    // This test asserts the structural contract that IS deterministic:
    //   (1) the rots-specific header (a single `d_stat` extra column, no
    //       AveExpr/B, no peptide_count);
    //   (2) log2FC cell-exact vs Python (the deterministic nanmean(A)-nanmean(B),
    //       the same sign convention as limma/limrots);
    //   (3) every pvalue/adj_pvalue in [0,1] and the rows sorted by adj_pvalue;
    //   (4) the d_stat column is finite.
    // The deterministic ROTS helpers (compute_d_stat/group_stats/calculate_p/
    // grids/rank_abs) are covered cell-by-cell in mokume-stats' rots unit tests.
    #[test]
    fn de_dispatcher_rots_wires_dstat_column_and_log2fc() -> TestResult<()> {
        let dir = temp_dir("oracle-rots")?;
        let output = dir.join("de_results.csv");
        let matrix = fixture_matrix()?;
        let sdrf = fixture_sdrf()?;

        run_differential_expression(&matrix, Some(&sdrf), &de_config(&output, "rots"), false)?;

        let text = std::fs::read_to_string(&output)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("missing header")?;
        // rots header: a single d_stat column between adj_pvalue and the means.
        assert_eq!(
            header,
            "ProteinName,log2FC,pvalue,adj_pvalue,d_stat,mean_groupA,mean_groupB,n_a,n_b,significance"
        );

        let mut rows: HashMap<String, Vec<String>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        for line in lines {
            let fields = line.split(',').map(ToOwned::to_owned).collect::<Vec<_>>();
            let protein = fields.first().ok_or("empty row")?.clone();
            order.push(protein.clone());
            rows.insert(protein, fields);
        }
        assert_eq!(rows.len(), EXPECTED_ROTS_LOG2FC.len(), "row count");

        for &(protein, log2fc) in EXPECTED_ROTS_LOG2FC {
            let fields = rows
                .get(protein)
                .ok_or_else(|| format!("missing {protein}"))?;
            assert_rots_row(fields, protein, log2fc)?;
        }
        assert_rots_sorted_by_adj_pvalue(&order, &rows)?;
        Ok(())
    }

    // LimROTS log2FC is `beta` = the full-data eBayes contrast coefficient
    // (mean_A - mean_B, the SAME sign convention as limma, the OPPOSITE of rots),
    // so P1 (low in A, high in B) is -2.0 here. Captured from limrots' full-data
    // eBayes on log2(RAW) (see scratchpad limrots oracle). Deterministic and
    // cell-exact regardless of the RNG.
    const EXPECTED_LIMROTS_LOG2FC: &[(&str, f64)] = &[
        ("P1", -2.0004868432854863),
        ("P2", 2.0090616097754097),
        ("P3", 0.015910691677655464),
        ("P4", -1.0),
        ("P5", -1.0),
        ("P6", 0.09175595082658106),
    ];

    // Pipeline wiring test for `--de-method limrots`. LimROTS is a faithful
    // RNG-based port (fixed internal PRNG), so its permutation p-values are
    // deterministic run-to-run but NOT bit-matched to Python (Python is itself
    // seed-unstable). This test asserts the structural contract that IS
    // deterministic:
    //   (1) the limrots-specific header (a single `t_stat` extra column, no
    //       AveExpr/B, no peptide_count);
    //   (2) log2FC cell-exact vs Python (the deterministic full-data eBayes beta,
    //       in limrots' mean_A - mean_B convention);
    //   (3) every pvalue/adj_pvalue in [0,1] and the rows sorted by adj_pvalue;
    //   (4) the t_stat (== d_stat) column is finite.
    // The deterministic LimROTS helpers (d_stat / boot_ebayes reorder /
    // p-value counting / s2_post) are covered cell-by-cell in mokume-stats'
    // limrots unit tests.
    #[test]
    fn de_dispatcher_limrots_wires_tstat_column_and_log2fc() -> TestResult<()> {
        let dir = temp_dir("oracle-limrots")?;
        let output = dir.join("de_results.csv");
        let matrix = fixture_matrix()?;
        let sdrf = fixture_sdrf()?;

        run_differential_expression(&matrix, Some(&sdrf), &de_config(&output, "limrots"), false)?;

        let text = std::fs::read_to_string(&output)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("missing header")?;
        // limrots header: a single t_stat column between adj_pvalue and the means.
        assert_eq!(
            header,
            "ProteinName,log2FC,pvalue,adj_pvalue,t_stat,mean_groupA,mean_groupB,n_a,n_b,significance"
        );

        let mut rows: HashMap<String, Vec<String>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        for line in lines {
            let fields = line.split(',').map(ToOwned::to_owned).collect::<Vec<_>>();
            let protein = fields.first().ok_or("empty row")?.clone();
            order.push(protein.clone());
            rows.insert(protein, fields);
        }
        assert_eq!(rows.len(), EXPECTED_LIMROTS_LOG2FC.len(), "row count");

        // The limrots field layout matches rots (single extra column at index 4),
        // so the rots row assertions apply unchanged.
        for &(protein, log2fc) in EXPECTED_LIMROTS_LOG2FC {
            let fields = rows
                .get(protein)
                .ok_or_else(|| format!("missing {protein}"))?;
            assert_rots_row(fields, protein, log2fc)?;
        }
        assert_rots_sorted_by_adj_pvalue(&order, &rows)?;
        Ok(())
    }

    /// Assert one rots output row. Field layout: ProteinName(0), log2FC(1),
    /// pvalue(2), adj_pvalue(3), d_stat(4), mean_A(5), mean_B(6), n_a(7), n_b(8),
    /// significance(9). log2FC is cell-exact (deterministic); pvalue/adj_pvalue
    /// are checked in-range and d_stat finite (the RNG-driven values are not
    /// bit-matched -- see the rots module docs).
    fn assert_rots_row(fields: &[String], protein: &str, log2fc: f64) -> TestResult<()> {
        assert_exact(parse(fields, 1)?, log2fc, "log2FC", protein);
        let p = parse(fields, 2)?;
        let adj = parse(fields, 3)?;
        assert!(
            (0.0..=1.0).contains(&p),
            "{protein} pvalue {p} out of range"
        );
        assert!(
            (0.0..=1.0).contains(&adj),
            "{protein} adj_pvalue {adj} out of range"
        );
        assert!(parse(fields, 4)?.is_finite(), "{protein} d_stat not finite");
        assert_eq!(field(fields, 7), Some("3"), "n_a {protein}");
        assert_eq!(field(fields, 8), Some("3"), "n_b {protein}");
        Ok(())
    }

    /// Assert the emitted rots rows are in ascending adj_pvalue order (BH sort).
    fn assert_rots_sorted_by_adj_pvalue(
        order: &[String],
        rows: &HashMap<String, Vec<String>>,
    ) -> TestResult<()> {
        let mut prev = f64::NEG_INFINITY;
        for name in order {
            let adj = parse(rows.get(name).ok_or("row")?, 3)?;
            assert!(adj >= prev, "rots rows not sorted by adj_pvalue at {name}");
            prev = adj;
        }
        Ok(())
    }

    /// Assert one deqms output row against the limma-fallback oracle. Field
    /// layout: ProteinName(0), log2FC(1), pvalue(2), adj_pvalue(3), sca_t(4),
    /// sca_pvalue(5), sca_adj_pvalue(6), mean_A(7), mean_B(8), n_a(9), n_b(10),
    /// peptide_count(11), log_pvalue(12), significance(13). In the fallback
    /// `sca_t == limma t`, `sca_pvalue == pvalue`,
    /// `sca_adj_pvalue == adj_pvalue`, and the `peptide_count` is the all-ones
    /// default.
    fn assert_deqms_row(
        fields: &[String],
        protein: &str,
        log2fc: f64,
        pvalue: f64,
        adj: f64,
        t_stat: f64,
        sig: &str,
    ) -> TestResult<()> {
        assert_close(parse(fields, 1)?, log2fc, "log2FC", protein);
        assert_close(parse(fields, 2)?, pvalue, "pvalue", protein);
        assert_close(parse(fields, 3)?, adj, "adj_pvalue", protein);
        assert_close(parse(fields, 4)?, t_stat, "sca_t", protein);
        assert_close(parse(fields, 5)?, pvalue, "sca_pvalue", protein);
        assert_close(parse(fields, 6)?, adj, "sca_adj_pvalue", protein);
        assert_eq!(field(fields, 9), Some("3"), "n_a {protein}");
        assert_eq!(field(fields, 10), Some("3"), "n_b {protein}");
        assert_eq!(field(fields, 11), Some("1"), "peptide_count {protein}");
        assert_close(parse(fields, 12)?, pvalue.ln(), "log_pvalue", protein);
        assert_eq!(field(fields, 13), Some(sig), "sig {protein}");
        Ok(())
    }

    fn field(fields: &[String], index: usize) -> Option<&str> {
        fields.get(index).map(String::as_str)
    }

    fn parse(fields: &[String], index: usize) -> TestResult<f64> {
        Ok(fields
            .get(index)
            .ok_or_else(|| format!("missing field {index}"))?
            .parse::<f64>()?)
    }

    fn assert_close(actual: f64, expected: f64, label: &str, protein: &str) {
        let tolerance = expected.abs().max(1.0) * 1e-6;
        assert!(
            (actual - expected).abs() <= tolerance,
            "{protein}.{label}: got {actual}, expected {expected}"
        );
    }

    // Cell-exact assertion for deterministic paths (absolute 1e-9).
    fn assert_exact(actual: f64, expected: f64, label: &str, protein: &str) {
        assert!(
            (actual - expected).abs() <= 1e-9,
            "{protein}.{label}: got {actual}, expected {expected}"
        );
    }

    /// Split one CSV line into fields, honouring `QUOTE_MINIMAL` quoting (so the
    /// ensemble's quoted `methods_significant` cell stays a single field). Only
    /// the features the ensemble writer emits are handled (double-quoted fields
    /// with doubled embedded quotes); this is a test helper, not a full parser.
    fn split_csv_line(line: &str) -> Vec<String> {
        let mut fields = Vec::new();
        let mut current = String::new();
        let mut in_quotes = false;
        let mut chars = line.chars().peekable();
        while let Some(c) = chars.next() {
            match c {
                '"' if in_quotes => {
                    if chars.peek() == Some(&'"') {
                        current.push('"');
                        chars.next();
                    } else {
                        in_quotes = false;
                    }
                }
                '"' => in_quotes = true,
                ',' if !in_quotes => {
                    fields.push(std::mem::take(&mut current));
                }
                other => current.push(other),
            }
        }
        fields.push(current);
        fields
    }

    /// An ensemble DE config with an explicit member list.
    fn ensemble_config(output: &std::path::Path, members: &[&str]) -> DifferentialExpressionConfig {
        DifferentialExpressionConfig {
            enabled: true,
            contrasts: Some(vec!["groupA vs groupB".to_string()]),
            method: "ensemble".to_string(),
            ensemble_methods: Some(members.iter().map(|m| (*m).to_string()).collect()),
            ensemble_min_k: 2,
            output: Some(output.to_path_buf()),
            ..DifferentialExpressionConfig::default()
        }
    }

    #[test]
    fn ensemble_members_are_canonicalized_once() -> TestResult<()> {
        let output = temp_dir("ensemble-member-canonicalization")?.join("de.csv");
        let mut config = ensemble_config(&output, &[" LimMA ", "DeQMS"]);
        config.ensemble_min_k = 1;

        assert_eq!(
            validated_ensemble_member_names(&config)?,
            vec!["limma".to_string(), "deqms".to_string()]
        );
        Ok(())
    }

    // End-to-end property test of `--de-method ensemble` over a member subset that
    // MIXES a deterministic member (limma) with an RNG member (rots). Bit-exact
    // end-to-end parity is NOT claimed here: rots/limrots are RNG-based and
    // proda is optimizer-driven (prior accepted decision -- see the rots/limrots/
    // proda module docs), so the combined p inherits their non-determinism. The
    // assertions are the contract that IS deterministic regardless of the RNG:
    //   (1) the ensemble column set/order;
    //   (2) every combined pvalue (when present) and adj_pvalue in [0, 1];
    //   (3) rows sorted ascending by adj_pvalue;
    //   (4) consensus UP/DOWN only where at least min_k members agree on the
    //       direction (n_methods_up >= 2 for UP, n_methods_down >= 2 for DOWN),
    //       and never both;
    //   (5) median log2FC is finite for every protein both members reported.
    #[test]
    fn de_dispatcher_ensemble_with_rng_member_is_valid() -> TestResult<()> {
        let dir = temp_dir("oracle-ensemble-rng")?;
        let output = dir.join("de_results.csv");
        let matrix = fixture_matrix()?;
        let sdrf = fixture_sdrf()?;

        run_differential_expression(
            &matrix,
            Some(&sdrf),
            &ensemble_config(&output, &["limma", "rots"]),
            false,
        )?;

        let text = std::fs::read_to_string(&output)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("missing header")?;
        assert_eq!(
            header,
            "Protein,log2FC,pvalue,n_methods_up,n_methods_down,methods_significant,adj_pvalue,significance"
        );

        let mut prev_adj = f64::NEG_INFINITY;
        let mut n_rows = 0usize;
        for line in lines {
            let fields = split_csv_line(line);
            n_rows += 1;
            let adj = assert_ensemble_row_valid(&fields)?;
            // adj_pvalue, when present, is non-decreasing (BH sort).
            if let Some(adj) = adj {
                let protein = fields.first().map(String::as_str).unwrap_or("?");
                assert!(
                    adj >= prev_adj,
                    "{protein} rows not sorted by adj_pvalue ({adj} < {prev_adj})"
                );
                prev_adj = adj;
            }
        }
        assert_eq!(n_rows, PROTEINS.len(), "ensemble emits one row per protein");
        Ok(())
    }

    /// Validate one ensemble output row's RNG-independent invariants (finite
    /// log2FC, pvalue/adj_pvalue in `[0,1]`, consensus direction agreeing with
    /// the member counts) and return the row's adjusted p-value when present, so
    /// the caller can check the global BH sort order.
    fn assert_ensemble_row_valid(fields: &[String]) -> TestResult<Option<f64>> {
        let protein = fields.first().ok_or("empty row")?.clone();
        assert!(parse(fields, 1)?.is_finite(), "{protein} log2FC not finite");
        assert_unit_interval(field(fields, 2), &protein, "pvalue")?;
        let adj = assert_unit_interval(field(fields, 6), &protein, "adj_pvalue")?;

        // Consensus: UP requires n_methods_up >= 2, DOWN requires
        // n_methods_down >= 2, and never both directions at once.
        let n_up = field(fields, 3).ok_or("n_up")?.parse::<usize>()?;
        let n_down = field(fields, 4).ok_or("n_down")?.parse::<usize>()?;
        let sig = field(fields, 7).ok_or("sig")?;
        assert!(
            !(n_up >= 2 && n_down >= 2),
            "{protein} cannot be both UP and DOWN consensus"
        );
        match sig {
            "UP" => assert!(n_up >= 2, "{protein} UP with only {n_up} up members"),
            "DOWN" => assert!(
                n_down >= 2,
                "{protein} DOWN with only {n_down} down members"
            ),
            "Unchanged" => {}
            other => return Err(format!("{protein} bad significance `{other}`").into()),
        }
        Ok(adj)
    }

    /// Parse an optional CSV cell as a probability in `[0,1]`; an empty/absent
    /// cell (a NaN combined or adjusted p) yields `None`.
    fn assert_unit_interval(
        raw: Option<&str>,
        protein: &str,
        label: &str,
    ) -> TestResult<Option<f64>> {
        match raw {
            Some(text) if !text.is_empty() => {
                let value = text.parse::<f64>()?;
                assert!(
                    (0.0..=1.0).contains(&value),
                    "{protein} {label} {value} out of range"
                );
                Ok(Some(value))
            }
            _ => Ok(None),
        }
    }
}
