//! ComBat batch correction.
//!
//! Port of `inmoose.pycombat.pycombat_norm` (the empirical-Bayes ComBat of
//! Johnson, Li & Rabinovic 2007) as mokume invokes it through
//! `postprocessing/batch_correction.py`. It covers the `ref_batch` (reference
//! batch left unmodified) and `mean_only` (additive effect only) options, plus:
//!   - optional covariates (`covar_mod`), whose biological signal is preserved;
//!   - both the parametric (`par_prior=true`) and non-parametric
//!     (`par_prior=false`) prior estimations.
//!
//! With no covariates and a one-hot batch design the two least-squares solves
//! in the reference reduce to per-batch means, so [`combat_parametric`] keeps a
//! solver-free fast path. When covariates are supplied the design becomes dense
//! ([batch one-hot | covariate columns]) and standardization regresses out
//! batch + covariates via a Gaussian-elimination solve of the normal equations
//! `B_hat = solve(design design^T, design data^T)` (`calculate_mean_var` in the
//! reference). The covariate columns mirror what mokume feeds inmoose: each
//! categorical-encoded covariate enters as a single numeric design column (the
//! `dmatrix("+".join(cols))` patsy build treats integer codes as continuous and
//! drops the redundant intercept).
//!
//! All reductions use population variance (`ddof = 0`), matching numpy. The
//! parametric `it_sol` fixed-point uses the reference's relative-change stop at
//! `1e-4`; the non-parametric path ports `int_eprior`'s deterministic numerical
//! integration (no RNG, `O(features^2)` per batch).

/// `it_sol` relative-change convergence threshold (`conv` in the reference).
const CONVERGENCE: f64 = 1e-4;
/// Hard iteration cap (`exit_iteration` in the reference).
const MAX_ITERATIONS: usize = 1_000_000;

/// Options for [`combat`], mirroring the `pycombat_norm` kwargs mokume exposes.
/// `ref_batch` is a batch *label* (not an index); `mean_only` skips the
/// multiplicative (variance) batch effect; `par_prior` selects the parametric
/// empirical-Bayes prior (`true`, the inmoose default) versus the non-parametric
/// `int_eprior` integration (`false`).
#[derive(Debug, Clone, Copy)]
pub struct ComBatParams {
    pub ref_batch: Option<usize>,
    pub mean_only: bool,
    pub par_prior: bool,
}

impl Default for ComBatParams {
    fn default() -> Self {
        // `par_prior` defaults to `true`, matching `pycombat_norm`'s signature.
        Self {
            ref_batch: None,
            mean_only: false,
            par_prior: true,
        }
    }
}

/// Parametric ComBat on a feature x sample matrix, no covariates.
///
/// Thin wrapper over [`combat`] preserving the original signature for callers
/// (the `correct-batches` CLI and the no-covariate golden tests). Forces
/// `par_prior = true` regardless of `params`, matching the historical
/// parametric-only contract.
pub fn combat_parametric(
    data: &[Vec<f64>],
    batch: &[usize],
    params: ComBatParams,
) -> Vec<Vec<f64>> {
    combat(
        data,
        batch,
        None,
        ComBatParams {
            par_prior: true,
            ..params
        },
    )
}

/// ComBat on a feature x sample matrix, with optional covariates.
///
/// `data[g]` is feature `g` across all samples; `batch[n]` is the batch label of
/// sample `n`. `covariates`, when present, is sample-major
/// (`covariates[n]` is the covariate row for sample `n`, all rows equal length),
/// matching the `samples x covariates` layout mokume builds from the SDRF before
/// handing it to `pycombat_norm(covar_mod=...)`. Returns the batch-corrected
/// matrix in the feature x sample orientation. Batches are processed in
/// ascending label order. The caller must ensure every batch has at least two
/// samples and `batch.len()` equals the sample count; a degenerate input (fewer
/// than two distinct batches, or an empty matrix) is returned unchanged.
pub fn combat(
    data: &[Vec<f64>],
    batch: &[usize],
    covariates: Option<&[Vec<f64>]>,
    params: ComBatParams,
) -> Vec<Vec<f64>> {
    let n_features = data.len();
    let n_samples = batch.len();
    if n_features == 0 || n_samples == 0 {
        return data.to_vec();
    }

    // Distinct batch labels in ascending order, and each batch's sample indices.
    let mut labels = batch.to_vec();
    labels.sort_unstable();
    labels.dedup();
    let n_batch = labels.len();
    if n_batch < 2 {
        return data.to_vec();
    }
    let batch_index = labels
        .iter()
        .enumerate()
        .map(|(index, label)| (*label, index))
        .collect::<std::collections::HashMap<_, _>>();
    let sample_batch = batch
        .iter()
        .map(|label| batch_index.get(label).copied().unwrap_or(0))
        .collect::<Vec<_>>();
    let reference = params
        .ref_batch
        .and_then(|label| batch_index.get(&label).copied());
    let mut batches_ind = vec![Vec::new(); n_batch];
    for (sample, &k) in sample_batch.iter().enumerate() {
        batches_ind[k].push(sample);
    }
    let batch_sizes = batches_ind.iter().map(Vec::len).collect::<Vec<_>>();

    // A covariate matrix with at least one column makes the design dense and
    // routes through the general normal-equations solve. With no usable
    // covariates the one-hot batch design reduces every solve to per-batch means
    // (the original solver-free path), which we keep for fidelity and speed.
    let covar_columns = covariates
        .and_then(|rows| rows.first().map(Vec::len))
        .filter(|&columns| columns > 0);

    let Standardized {
        s_data,
        grand_mean,
        std_dev,
        cov_mean,
    } = match covar_columns {
        Some(columns) => {
            match standardize_with_covariates(
                data,
                covariates.unwrap_or(&[]),
                columns,
                &sample_batch,
                &batches_ind,
                &batch_sizes,
                n_batch,
                reference,
            ) {
                Some(standardized) => standardized,
                // Singular normal equations (e.g. a covariate confounded with
                // the batches): fall back to returning the input unchanged
                // rather than emitting NaNs.
                None => return data.to_vec(),
            }
        }
        None => standardize_one_hot(
            data,
            &sample_batch,
            &batches_ind,
            &batch_sizes,
            n_batch,
            reference,
        ),
    };

    // Per-batch L/S estimates on standardized data. gamma_hat is the per-batch
    // mean of the standardized data (the batch-design solve with a one-hot
    // batch_design), identical with and without covariates.
    let gamma_hat = (0..n_batch)
        .map(|k| {
            (0..n_features)
                .map(|g| mean(batches_ind[k].iter().map(|&n| s_data[g][n])))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let delta_hat = (0..n_batch)
        .map(|k| {
            (0..n_features)
                .map(|g| {
                    if params.mean_only {
                        1.0
                    } else {
                        variance_pop(
                            batches_ind[k].iter().map(|&n| s_data[g][n]),
                            gamma_hat[k][g],
                        )
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    // Empirical-Bayes shrinkage per batch.
    let mut gamma_star = vec![vec![0.0; n_features]; n_batch];
    let mut delta_star = vec![vec![1.0; n_features]; n_batch];
    for k in 0..n_batch {
        let gamma_bar = mean(gamma_hat[k].iter().copied());
        let t2 = variance_pop(gamma_hat[k].iter().copied(), gamma_bar);
        if params.mean_only {
            // Closed-form additive effect (n = 1); multiplicative effect = 1.
            gamma_star[k] = (0..n_features)
                .map(|g| (t2 * gamma_hat[k][g] + gamma_bar) / (t2 + 1.0))
                .collect();
        } else if params.par_prior {
            let (a_prior, b_prior) = inverse_gamma_prior(&delta_hat[k]);
            let priors = BatchPriors {
                gamma_bar,
                t2,
                a_prior,
                b_prior,
            };
            let (g_star, d_star) = it_sol(
                &s_data,
                &batches_ind[k],
                &gamma_hat[k],
                &delta_hat[k],
                priors,
            );
            gamma_star[k] = g_star;
            delta_star[k] = d_star;
        } else {
            let (g_star, d_star) =
                int_eprior(&s_data, &batches_ind[k], &gamma_hat[k], &delta_hat[k]);
            gamma_star[k] = g_star;
            delta_star[k] = d_star;
        }
    }

    // The reference batch is left unmodified.
    if let Some(r) = reference {
        gamma_star[r] = vec![0.0; n_features];
        delta_star[r] = vec![1.0; n_features];
    }

    // Adjust and map back to the original scale; reference-batch samples are
    // restored to the raw input exactly. `stand_mean[g][n] = grand_mean[g] +
    // cov_mean[g][n]` (the covariate contribution, zero in the one-hot path).
    (0..n_features)
        .map(|g| {
            (0..n_samples)
                .map(|n| {
                    let k = sample_batch[n];
                    if reference == Some(k) {
                        data[g][n]
                    } else {
                        let adjusted = (s_data[g][n] - gamma_star[k][g]) / delta_star[k][g].sqrt();
                        let stand_mean = grand_mean[g] + cov_mean[g][n];
                        adjusted * std_dev[g] + stand_mean
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Standardization outputs shared by both design paths. `s_data` is the
/// standardized matrix; `grand_mean[g]` and `std_dev[g]` are per-feature;
/// `cov_mean[g][n]` is the per-sample covariate contribution to the standardized
/// mean (all zero when there are no covariates).
struct Standardized {
    s_data: Vec<Vec<f64>>,
    grand_mean: Vec<f64>,
    std_dev: Vec<f64>,
    cov_mean: Vec<Vec<f64>>,
}

/// One-hot (no-covariate) standardization: the two reference solves collapse to
/// per-batch means, so `B_hat[k][g]` is the batch-`k` mean of feature `g`.
/// Mirrors `calculate_mean_var` / `standardise_data` for a pure batch design.
fn standardize_one_hot(
    data: &[Vec<f64>],
    sample_batch: &[usize],
    batches_ind: &[Vec<usize>],
    batch_sizes: &[usize],
    n_batch: usize,
    reference: Option<usize>,
) -> Standardized {
    let n_features = data.len();
    let n_samples = sample_batch.len();

    let b_hat = (0..n_batch)
        .map(|k| {
            (0..n_features)
                .map(|g| mean(batches_ind[k].iter().map(|&n| data[g][n])))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    // grand_mean: the reference batch intercept when a reference is set,
    // otherwise the size-weighted average of the batch intercepts.
    let grand_mean = (0..n_features)
        .map(|g| match reference {
            Some(r) => b_hat[r][g],
            None => (0..n_batch)
                .map(|k| batch_sizes[k] as f64 / n_samples as f64 * b_hat[k][g])
                .sum::<f64>(),
        })
        .collect::<Vec<_>>();

    // var_pooled[g]: residual variance (ddof=0) about the batch means, taken
    // over the reference batch only when a reference is set, otherwise all
    // samples.
    let var_pooled = (0..n_features)
        .map(|g| match reference {
            Some(r) => {
                let sum_sq = batches_ind[r]
                    .iter()
                    .map(|&n| {
                        let residual = data[g][n] - b_hat[r][g];
                        residual * residual
                    })
                    .sum::<f64>();
                sum_sq / batch_sizes[r] as f64
            }
            None => {
                let sum_sq = (0..n_samples)
                    .map(|n| {
                        let residual = data[g][n] - b_hat[sample_batch[n]][g];
                        residual * residual
                    })
                    .sum::<f64>();
                sum_sq / n_samples as f64
            }
        })
        .collect::<Vec<_>>();

    let std_dev = var_pooled.iter().map(|v| v.sqrt()).collect::<Vec<_>>();
    let s_data = (0..n_features)
        .map(|g| {
            (0..n_samples)
                .map(|n| (data[g][n] - grand_mean[g]) / std_dev[g])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Standardized {
        s_data,
        grand_mean,
        std_dev,
        cov_mean: vec![vec![0.0; n_samples]; n_features],
    }
}

/// Covariate-aware standardization via the dense normal-equations solve, porting
/// `calculate_mean_var` / `calculate_stand_mean` / `standardise_data`.
///
/// The design (rows = predictors, columns = samples) is
/// `[batch one-hot (n_batch rows) | covariate columns]`. When a reference batch
/// is set its one-hot row is forced to all-ones (matching inmoose's
/// `batch_mod[:, ref] = 1`). `B_hat = solve(design design^T, design data^T)` is
/// the predictor x feature coefficient matrix. Returns `None` if the normal
/// equations are singular.
#[allow(clippy::too_many_arguments)]
fn standardize_with_covariates(
    data: &[Vec<f64>],
    covariates: &[Vec<f64>],
    covar_columns: usize,
    sample_batch: &[usize],
    batches_ind: &[Vec<usize>],
    batch_sizes: &[usize],
    n_batch: usize,
    reference: Option<usize>,
) -> Option<Standardized> {
    let n_features = data.len();
    let n_samples = sample_batch.len();
    if covariates.len() != n_samples {
        return None;
    }

    // design[row][sample]: batch one-hot rows then covariate columns.
    let n_cols = n_batch + covar_columns;
    let mut design = vec![vec![0.0; n_samples]; n_cols];
    for (sample, &k) in sample_batch.iter().enumerate() {
        design[k][sample] = 1.0;
    }
    // Reference batch row is forced all-ones (inmoose `batch_mod[:, ref] = 1`).
    if let Some(r) = reference {
        for value in &mut design[r] {
            *value = 1.0;
        }
    }
    for c in 0..covar_columns {
        for (sample, row) in covariates.iter().enumerate() {
            design[n_batch + c][sample] = *row.get(c)?;
        }
    }

    // Normal equations: gram = design design^T (n_cols x n_cols);
    // rhs[g] = design data[g]^T (length n_cols). B_hat[col][g] solves gram.
    let gram = (0..n_cols)
        .map(|i| {
            (0..n_cols)
                .map(|j| {
                    (0..n_samples)
                        .map(|n| design[i][n] * design[j][n])
                        .sum::<f64>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut b_hat = vec![vec![0.0; n_features]; n_cols];
    for g in 0..n_features {
        let rhs = (0..n_cols)
            .map(|i| {
                (0..n_samples)
                    .map(|n| design[i][n] * data[g][n])
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let solution = solve_linear_system(&gram, &rhs)?;
        for (col, value) in solution.into_iter().enumerate() {
            b_hat[col][g] = value;
        }
    }

    // grand_mean[g]: reference batch intercept, else size-weighted batch
    // intercepts. With a reference, B_hat[ref] is that batch's intercept.
    let grand_mean = (0..n_features)
        .map(|g| match reference {
            Some(r) => b_hat[r][g],
            None => (0..n_batch)
                .map(|k| batch_sizes[k] as f64 / n_samples as f64 * b_hat[k][g])
                .sum::<f64>(),
        })
        .collect::<Vec<_>>();

    // fitted[g][n] = (design^T B_hat)[n][g]: the full regression prediction.
    let fitted = |g: usize, n: usize| -> f64 {
        (0..n_cols)
            .map(|col| design[col][n] * b_hat[col][g])
            .sum::<f64>()
    };

    // var_pooled[g]: residual variance about the full fit, over the reference
    // batch only when set, otherwise all samples (`calculate_mean_var`).
    let var_pooled = (0..n_features)
        .map(|g| match reference {
            Some(r) => {
                let sum_sq = batches_ind[r]
                    .iter()
                    .map(|&n| {
                        let residual = data[g][n] - fitted(g, n);
                        residual * residual
                    })
                    .sum::<f64>();
                sum_sq / batch_sizes[r] as f64
            }
            None => {
                let sum_sq = (0..n_samples)
                    .map(|n| {
                        let residual = data[g][n] - fitted(g, n);
                        residual * residual
                    })
                    .sum::<f64>();
                sum_sq / n_samples as f64
            }
        })
        .collect::<Vec<_>>();
    let std_dev = var_pooled.iter().map(|v| v.sqrt()).collect::<Vec<_>>();

    // cov_mean[g][n]: the covariate-only contribution to the standardized mean
    // (`calculate_stand_mean` zeroes the batch rows before the dot product).
    let cov_mean = (0..n_features)
        .map(|g| {
            (0..n_samples)
                .map(|n| {
                    (n_batch..n_cols)
                        .map(|col| design[col][n] * b_hat[col][g])
                        .sum::<f64>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    // s_data[g][n] = (data - grand_mean - cov_mean) / std_dev.
    let s_data = (0..n_features)
        .map(|g| {
            (0..n_samples)
                .map(|n| (data[g][n] - grand_mean[g] - cov_mean[g][n]) / std_dev[g])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Some(Standardized {
        s_data,
        grand_mean,
        std_dev,
        cov_mean,
    })
}

/// Solve the square system `a x = b` by Gauss-Jordan elimination with partial
/// pivoting (the same scheme as `mokume-imputation`'s solver, inlined here to
/// avoid a cross-crate dependency). Returns `None` when `a` is singular,
/// mirroring `np.linalg.solve` raising `LinAlgError`.
fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut augmented = a
        .iter()
        .zip(b)
        .map(|(row, &rhs)| {
            let mut extended = row.clone();
            extended.push(rhs);
            extended
        })
        .collect::<Vec<Vec<f64>>>();

    for column in 0..n {
        let pivot = (column..n).max_by(|&left, &right| {
            augmented[left][column]
                .abs()
                .total_cmp(&augmented[right][column].abs())
        })?;
        if augmented[pivot][column].abs() < 1e-300 {
            return None;
        }
        augmented.swap(column, pivot);

        let pivot_value = augmented[column][column];
        for value in &mut augmented[column] {
            *value /= pivot_value;
        }

        let (before, rest) = augmented.split_at_mut(column);
        let (pivot_row, after) = rest.split_first_mut()?;
        for row in before.iter_mut().chain(after.iter_mut()) {
            let factor = row[column];
            if factor != 0.0 {
                for (target, &source) in row.iter_mut().zip(pivot_row.iter()) {
                    *target -= factor * source;
                }
            }
        }
    }

    Some(augmented.iter().map(|row| row[n]).collect())
}

/// Method-of-moments inverse-gamma prior `(a_prior, b_prior)` from a batch's
/// per-feature variances (`compute_prior` in the reference).
fn inverse_gamma_prior(delta_hat: &[f64]) -> (f64, f64) {
    let m = mean(delta_hat.iter().copied());
    let s2 = variance_pop(delta_hat.iter().copied(), m);
    let a_prior = (2.0 * s2 + m * m) / s2;
    let b_prior = (m * s2 + m * m * m) / s2;
    (a_prior, b_prior)
}

/// Per-batch empirical-Bayes prior hyperparameters for [`it_sol`].
struct BatchPriors {
    gamma_bar: f64,
    t2: f64,
    a_prior: f64,
    b_prior: f64,
}

/// `it_sol`: parametric empirical-Bayes fixed-point for one batch's
/// `(gamma_star, delta_star)`.
fn it_sol(
    s_data: &[Vec<f64>],
    samples: &[usize],
    gamma_hat: &[f64],
    delta_hat: &[f64],
    priors: BatchPriors,
) -> (Vec<f64>, Vec<f64>) {
    let BatchPriors {
        gamma_bar,
        t2,
        a_prior,
        b_prior,
    } = priors;
    let n = samples.len() as f64;
    let t2_n = t2 * n;
    let mut g_old = gamma_hat.to_vec();
    let mut d_old = delta_hat.to_vec();
    let n_features = gamma_hat.len();

    for _ in 0..MAX_ITERATIONS {
        let g_new = (0..n_features)
            .map(|g| (t2_n * gamma_hat[g] + d_old[g] * gamma_bar) / (t2_n + d_old[g]))
            .collect::<Vec<_>>();
        let d_new = (0..n_features)
            .map(|g| {
                let sum_sq = samples
                    .iter()
                    .map(|&sample| {
                        let residual = s_data[g][sample] - g_new[g];
                        residual * residual
                    })
                    .sum::<f64>();
                (0.5 * sum_sq + b_prior) / (0.5 * n + a_prior - 1.0)
            })
            .collect::<Vec<_>>();

        let change = (0..n_features)
            .map(|g| {
                let gamma_change = (g_new[g] - g_old[g]).abs() / g_old[g].abs();
                let delta_change = (d_new[g] - d_old[g]).abs() / d_old[g].abs();
                gamma_change.max(delta_change)
            })
            .fold(0.0_f64, f64::max);

        g_old = g_new;
        d_old = d_new;
        if change <= CONVERGENCE {
            break;
        }
    }
    (g_old, d_old)
}

/// `int_eprior`: non-parametric empirical-Bayes estimation for one batch's
/// `(gamma_star, delta_star)`, ported from the reference's deterministic
/// numerical integration (no Monte Carlo / RNG despite the docstring's wording).
///
/// For each feature `i` the prior is integrated over all *other* features'
/// `(gamma_hat, delta_hat)` in this batch: with `x` the standardized batch row
/// of feature `i` (length `n` = batch size), `g`/`d` the other features'
/// estimates, and `sum2[k] = sum_j (x_j - g_k)^2`, the likelihood weight is
/// `LH[k] = (1 / (pi * 2 d_k))^(n/2) * exp(-sum2[k] / (2 d_k))`. The posterior
/// means are `sum(g LH) / sum(LH)` and `sum(d LH) / sum(LH)`. NaN weights are
/// zeroed (`np.nan_to_num`); if every weight underflows to zero they are reset
/// to `exp(-745)` (the reference's underflow guard) so the ratio stays defined.
fn int_eprior(
    s_data: &[Vec<f64>],
    samples: &[usize],
    gamma_hat: &[f64],
    delta_hat: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n_features = gamma_hat.len();
    let n = samples.len() as f64;
    let half_n = n / 2.0;
    let mut g_star = vec![0.0; n_features];
    let mut d_star = vec![0.0; n_features];

    for i in 0..n_features {
        // Likelihood weight of feature i's batch row under every OTHER feature's
        // (gamma_hat, delta_hat) prior.
        let mut weights = Vec::with_capacity(n_features.saturating_sub(1));
        for k in 0..n_features {
            if k == i {
                continue;
            }
            let g_k = gamma_hat[k];
            let two_d_k = 2.0 * delta_hat[k];
            let sum_sq = samples
                .iter()
                .map(|&sample| {
                    let residual = s_data[i][sample] - g_k;
                    residual * residual
                })
                .sum::<f64>();
            let mut weight =
                (1.0 / (std::f64::consts::PI * two_d_k)).powf(half_n) * (-sum_sq / two_d_k).exp();
            if weight.is_nan() {
                // np.nan_to_num: NaN -> 0.0.
                weight = 0.0;
            }
            weights.push((gamma_hat[k], delta_hat[k], weight));
        }

        let total: f64 = weights.iter().map(|&(_, _, w)| w).sum();
        if total == 0.0 {
            // Reference underflow guard: every zero weight becomes exp(-745).
            let floor = (-745.0_f64).exp();
            for entry in &mut weights {
                if entry.2 == 0.0 {
                    entry.2 = floor;
                }
            }
        }
        let denom: f64 = weights.iter().map(|&(_, _, w)| w).sum();
        let g_num: f64 = weights.iter().map(|&(g, _, w)| g * w).sum();
        let d_num: f64 = weights.iter().map(|&(_, d, w)| d * w).sum();
        g_star[i] = g_num / denom;
        d_star[i] = d_num / denom;
    }

    (g_star, d_star)
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut count = 0usize;
    for value in values {
        sum += value;
        count += 1;
    }
    if count == 0 {
        0.0
    } else {
        sum / count as f64
    }
}

/// Population variance (`ddof = 0`) about a precomputed mean.
fn variance_pop(values: impl Iterator<Item = f64>, mean: f64) -> f64 {
    let mut sum_sq = 0.0;
    let mut count = 0usize;
    for value in values {
        let deviation = value - mean;
        sum_sq += deviation * deviation;
        count += 1;
    }
    if count == 0 {
        0.0
    } else {
        sum_sq / count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::{combat, combat_parametric, ComBatParams};

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "actual={actual} expected={expected}"
        );
    }

    // Oracle from inmoose `pycombat_norm` (parametric, no covariates): 4
    // features x 6 samples, 2 batches of 3. Reference output captured by
    // `conda run -n Bigbio python /tmp/combat_oracle.py`.
    #[test]
    fn matches_inmoose_pycombat_oracle() {
        let data = vec![
            vec![10.0, 11.0, 9.5, 20.0, 21.0, 19.0],
            vec![5.0, 6.0, 4.0, 8.0, 7.5, 9.0],
            vec![1.0, 2.0, 1.5, 3.0, 2.5, 4.0],
            vec![50.0, 52.0, 48.0, 30.0, 31.0, 29.0],
        ];
        let batch = [0, 0, 0, 1, 1, 1];
        let corrected = combat_parametric(&data, &batch, ComBatParams::default());

        let expected = [
            [
                14.864380299149872,
                15.864215988028887,
                14.364462454710365,
                15.135873445153031,
                16.131298305630697,
                14.140448584675365,
            ],
            [
                6.570_036_097_877_598,
                7.530_315_584_683_231,
                5.6097566110719645,
                6.420_867_458_192_403,
                5.896_843_353_383_13,
                7.468_915_667_810_949,
            ],
            [
                1.8221230342951276,
                2.8335983703721137,
                2.3278607023336204,
                2.1750301068609867,
                1.6832478972769493,
                3.158594526029062,
            ],
            [
                40.137_786_918_348_62,
                42.01029911676347,
                38.265_274_719_933_77,
                39.880718344542714,
                40.965961295069526,
                38.795475394015895,
            ],
        ];
        for (feature, expected_row) in expected.iter().enumerate() {
            for (sample, &want) in expected_row.iter().enumerate() {
                assert_close(corrected[feature][sample], want, 1e-6);
            }
        }
    }

    fn oracle_input() -> Vec<Vec<f64>> {
        vec![
            vec![10.0, 11.0, 9.5, 20.0, 21.0, 19.0],
            vec![5.0, 6.0, 4.0, 8.0, 7.5, 9.0],
            vec![1.0, 2.0, 1.5, 3.0, 2.5, 4.0],
            vec![50.0, 52.0, 48.0, 30.0, 31.0, 29.0],
        ]
    }

    // Reference batch 0: batch-0 columns are returned as the raw input, batch-1
    // shifted onto batch 0. Oracle from `pycombat_norm(..., ref_batch=0)`.
    #[test]
    fn matches_inmoose_ref_batch_oracle() {
        let data = oracle_input();
        let batch = [0, 0, 0, 1, 1, 1];
        let params = ComBatParams {
            ref_batch: Some(0),
            mean_only: false,
            ..ComBatParams::default()
        };
        let corrected = combat_parametric(&data, &batch, params);
        let expected = [
            [
                10.0,
                11.0,
                9.5,
                10.19824941563275,
                11.048785837051422,
                9.347712994214078,
            ],
            [
                5.0,
                6.0,
                4.0,
                4.837155840740258,
                4.340315404980416,
                5.830836712259942,
            ],
            [
                1.0,
                2.0,
                1.5,
                1.3696047042858863,
                0.9721857919446014,
                2.164442528968456,
            ],
            [
                50.0,
                52.0,
                48.0,
                49.92144658970713,
                50.972711406794325,
                48.87018177261993,
            ],
        ];
        for (feature, expected_row) in expected.iter().enumerate() {
            for (sample, &want) in expected_row.iter().enumerate() {
                assert_close(corrected[feature][sample], want, 1e-6);
            }
        }
    }

    // mean_only: multiplicative batch effect fixed at 1, additive effect from
    // the closed-form shrinkage. Oracle from `pycombat_norm(..., mean_only=True)`.
    #[test]
    fn matches_inmoose_mean_only_oracle() {
        let data = oracle_input();
        let batch = [0, 0, 0, 1, 1, 1];
        let params = ComBatParams {
            ref_batch: None,
            mean_only: true,
            ..ComBatParams::default()
        };
        let corrected = combat_parametric(&data, &batch, params);
        let expected = [
            [
                14.763385418186564,
                15.763385418186564,
                14.263385418186564,
                15.236614581813432,
                16.236614581813434,
                14.236614581813432,
            ],
            [
                6.545876060049075,
                7.545876060049075,
                5.545876060049075,
                6.454123939950924,
                5.954123939950924,
                7.454123939950925,
            ],
            [
                1.8171160603917724,
                2.8171160603917724,
                2.3171160603917724,
                2.182883939608227,
                1.6828839396082271,
                3.182883939608227,
            ],
            [
                40.3786752916639,
                42.3786752916639,
                38.3786752916639,
                39.621_324_708_336_1,
                40.621_324_708_336_1,
                38.621_324_708_336_1,
            ],
        ];
        for (feature, expected_row) in expected.iter().enumerate() {
            for (sample, &want) in expected_row.iter().enumerate() {
                assert_close(corrected[feature][sample], want, 1e-6);
            }
        }
    }

    // Covariate-aware standardization (parametric prior). One categorical
    // covariate, sample-major (`samples x covariates`), values [0,1,0,1,0,1],
    // exactly as mokume feeds `pycombat_norm(covar_mod=...)`. The covariate
    // enters the patsy design as a single numeric column (integer codes are
    // treated as continuous, the redundant intercept is dropped). Oracle from
    // `pycombat_norm(..., covar_mod=[[0],[1],[0],[1],[0],[1]])`.
    #[test]
    fn matches_inmoose_covariate_oracle() {
        let data = oracle_input();
        let batch = [0, 0, 0, 1, 1, 1];
        let covariates = vec![
            vec![0.0],
            vec![1.0],
            vec![0.0],
            vec![1.0],
            vec![0.0],
            vec![1.0],
        ];
        let corrected = combat(&data, &batch, Some(&covariates), ComBatParams::default());
        let expected = [
            [
                14.880324738490788,
                15.96298968871664,
                14.343584760612632,
                15.107684002440037,
                16.054497868898974,
                14.168468155058392,
            ],
            [
                6.392477045131027,
                7.376279455981547,
                5.327686688533107,
                6.668850019948931,
                6.126252158939517,
                7.6120528719363785,
            ],
            [
                1.6610698901054048,
                2.674853261167651,
                2.2162033743543903,
                2.351060726841906,
                1.8228943013580798,
                3.275950258885035,
            ],
            [
                40.00325131995268,
                42.06893034764898,
                37.8981648756386,
                40.01045716540908,
                40.919729685609504,
                39.06230143958027,
            ],
        ];
        for (feature, expected_row) in expected.iter().enumerate() {
            for (sample, &want) in expected_row.iter().enumerate() {
                assert_close(corrected[feature][sample], want, 1e-9);
            }
        }
    }

    // Non-parametric prior (`par_prior=false`), no covariates: replaces the
    // `it_sol` fixed-point with the deterministic `int_eprior` integration.
    // Oracle from `pycombat_norm(..., par_prior=False)`.
    #[test]
    fn matches_inmoose_nonparametric_oracle() {
        let data = oracle_input();
        let batch = [0, 0, 0, 1, 1, 1];
        let params = ComBatParams {
            ref_batch: None,
            mean_only: false,
            par_prior: false,
        };
        let corrected = combat(&data, &batch, None, params);
        let expected = [
            [
                11.9691855098242,
                12.85894203082681,
                11.524307249322895,
                18.267871372322563,
                19.113025751618693,
                17.42271699302643,
            ],
            [
                6.02218632450599,
                7.313180773241795,
                4.731191875770184,
                6.809831239944923,
                6.387254112580213,
                7.6549854946743405,
            ],
            [
                2.1690285772457,
                3.0587850982483094,
                2.6139068377470047,
                1.7718155933606998,
                1.1893332208499905,
                2.936780338382118,
            ],
            [
                51.40103502432606,
                53.180548066331276,
                49.62152198232084,
                29.823293554358948,
                30.668447809087468,
                28.978139299630435,
            ],
        ];
        for (feature, expected_row) in expected.iter().enumerate() {
            for (sample, &want) in expected_row.iter().enumerate() {
                assert_close(corrected[feature][sample], want, 1e-9);
            }
        }
    }
}
