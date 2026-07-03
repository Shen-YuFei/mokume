use mokume_core::{ProteinId, SampleId};

use crate::support::{fill_missing_by_sample, finite_sample_values};

/// Standard-deviation floor mirroring Python's `if sd_obs < 1e-10: sd_obs = 0.1`.
/// A near-constant column would otherwise place the imputation center exactly on
/// the truncation point; the floor keeps it strictly below.
const SD_FLOOR_THRESHOLD: f64 = 1e-10;
const SD_FLOOR_VALUE: f64 = 0.1;

/// Multiplier applied to the estimated standard deviation before drawing,
/// matching the Python default `tune_sigma = 1.0`. It scales only the draw
/// spread `sd_imp = sd_obs * tune_sigma`, never the distribution center, so it
/// does not influence the deterministic fill below. There is no `ImputationConfig`
/// field for it, so the default is fixed here.
const TUNE_SIGMA: f64 = 1.0;

/// QRILC (Quantile Regression Imputation of Left-Censored data).
///
/// Faithful port of `mokume.imputation.qrilc.impute_qrilc`. For each sample
/// (column) a truncated normal distribution is fitted to the lower tail of the
/// observed log2 abundances and missing entries are replaced from that tail.
/// Per sample the Python reference computes:
///   - `observed` = the finite values in that column (requires at least two)
///   - `sd_obs` = sample standard deviation (`numpy.std(observed, ddof=1)`),
///     floored to `0.1` when below `1e-10`
///   - `trunc` = the minimum observed value (the left-censoring point)
///   - `mu_imp` = `trunc - sd_obs` (the imputation distribution center)
///   - `sd_imp` = `sd_obs * tune_sigma` (the draw spread)
///
/// It then draws `N(mu_imp, sd_imp)` for each missing cell and clamps each draw
/// to be at most `trunc`.
///
/// Determinism note: the Python implementation is STOCHASTIC. It samples each
/// replacement from `numpy.random.normal(mu_imp, sd_imp)`, so values change with
/// the global NumPy seed and cannot match cross-language cell-for-cell.
/// Following the same convention used by `minprob` and `gms` in this crate, the
/// Rust port deterministically fills the distribution center `mu_imp` (the
/// expected value of those draws) instead of sampling, and never calls an RNG.
/// The clamp `min(draw, trunc)` is a no-op for the center: `mu_imp = trunc -
/// sd_obs < trunc` whenever `sd_obs > 0` (it is floored to `0.1`), so the center
/// already lies at or below the truncation point. The fill is therefore exactly
/// `trunc - sd_obs`, which is reproducible for a given matrix.
///
/// Skip rule: a column with no missing cell, or with fewer than two observed
/// values, is left untouched (`continue` in Python). This matches the reference
/// exactly and is deterministic.
///
/// Standard deviation: QRILC needs the *sample* standard deviation (`ddof = 1`),
/// unlike `mokume_core::stats::finite_sd` which is the population deviation
/// (`ddof = 0`); the sample variant is computed locally below.
///
/// Values arrive already in log2 space (the pipeline applies `log2` before and
/// `2**` after), so this routine works directly on the supplied values and
/// returns only the originally-missing cells. No linear algebra is required, so
/// `crate::linalg::solve_linear_system` is not used here.
///
/// Reference: Wei R, et al. Missing Value Imputation Approach for Mass
/// Spectrometry-based Metabolomics Data. Sci Rep. 2018;8:663.
pub(crate) fn qrilc_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let fills = samples
        .iter()
        .filter_map(|sample| {
            let observed = finite_sample_values(proteins, *sample, value_at);
            sample_center(&observed).map(|center| (*sample, center))
        })
        .collect::<std::collections::HashMap<_, _>>();
    fill_missing_by_sample(proteins, samples, &fills, value_at)
}

/// Deterministic QRILC fill for a single sample from its observed log2 values.
///
/// Returns `None` (no fill) when there are fewer than two observed values,
/// matching Python's `if len(observed) < 2: continue`. Otherwise returns the
/// truncated-normal center `mu_imp = trunc - sd_obs`, where `trunc` is the
/// minimum observed value and `sd_obs` is the sample standard deviation
/// (`ddof = 1`) floored to `0.1`.
fn sample_center(observed: &[f64]) -> Option<f64> {
    if observed.len() < 2 {
        return None;
    }
    let sd_obs = sample_std_ddof1(observed)?;
    let sd_obs = if sd_obs < SD_FLOOR_THRESHOLD {
        SD_FLOOR_VALUE
    } else {
        sd_obs
    };
    // sd_imp = sd_obs * TUNE_SIGMA scales only the draw spread, which the
    // deterministic center does not depend on; referenced here for parity.
    let _sd_imp = sd_obs * TUNE_SIGMA;
    let trunc = observed.iter().copied().fold(f64::INFINITY, f64::min);
    let center = trunc - sd_obs;
    center.is_finite().then_some(center)
}

/// Sample standard deviation with `ddof = 1`, matching
/// `numpy.std(values, ddof=1)`. Returns `None` for fewer than two values (the
/// caller already guards this). The variance is `sum((x - mean)^2) / (n - 1)`.
fn sample_std_ddof1(values: &[f64]) -> Option<f64> {
    let n = values.len();
    if n < 2 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / (n - 1) as f64;
    Some(variance.sqrt())
}
