// proDA EM driver and Wald test. Included from `proda.rs`.
//
// Ports proda.py:514 `run_proda`: seed-dependent imputation, OLS init, prior
// estimation, the 20-iteration EM loop (proda.py:569) with its convergence
// check, and the final per-protein Wald test (proda.py:662-691).

/// Per-protein proDA statistics returned by [`run_proda`], one entry per input
/// row (in input order; the public wrapper applies the testability filter and
/// BH adjustment). `t_stat` is the Wald statistic (the proDA extra column).
pub(crate) struct ProdaStat {
    pub log2_fold_change: f64,
    pub t_stat: f64,
    pub p_value: f64,
}

/// proda.py:520 seed for the initial imputation (only seeds the EM trajectory).
const DEFAULT_SEED: u64 = 42;
/// proda.py:569 EM iteration cap.
const EM_ITERS: usize = 20;
/// proda.py:657 EM convergence threshold on the aggregate parameter change.
const EM_TOL: f64 = 1e-3;

/// The current EM parameter state (priors + per-sample dropout curves).
struct EmState {
    mu0: f64,
    sigma20: f64,
    tau20: f64,
    df0: f64,
    rho: Vec<f64>,
    zeta: Vec<f64>,
    pred: Vec<Vec<f64>>,
    pred_var: Vec<Vec<f64>>,
}

/// Predictions / variances / s2 / df extracted from one round of fits
/// (proda.py:603 `_extract`).
struct Extracted {
    pred: Vec<Vec<f64>>,
    pred_var: Vec<Vec<f64>>,
    s2: Vec<f64>,
    df: Vec<f64>,
}

/// proda.py:514 `run_proda`. `mat` is the protein x sample log2 matrix (group A
/// columns first, then B) with NaN for missing; returns one [`ProdaStat`] per
/// protein in input order.
fn run_proda(mat: &[Vec<f64>], n_a: usize, n_b: usize) -> Vec<ProdaStat> {
    let n_samples = n_a + n_b;
    let x = Design::two_group(n_a, n_b);

    // proda.py:545-558 initial imputation (SplitMix64 instead of MT19937).
    let y_init = impute_initial(mat, n_samples, DEFAULT_SEED);

    // proda.py:560 OLS init.
    let ols = ols_per_protein(&y_init, &x);
    // proda.py:562-564 prior init.
    let (mu0, sigma20) = location_prior(&ols.pred, &ols.pred_var, None, None);
    let (rho, zeta) = fit_dropout_curves(mat, &x, &ols.pred, &ols.pred_var);
    let (tau20, df0) = variance_prior(&ols.s2, &ols.df);

    let mut state = EmState {
        mu0,
        sigma20,
        tau20,
        df0,
        rho,
        zeta,
        pred: ols.pred,
        pred_var: ols.pred_var,
    };

    let mut fits_reg: Vec<PdLmFit> = Vec::new();
    for _ in 0..EM_ITERS {
        let step = em_step(mat, &x, &state);
        let converged = step.err < EM_TOL;
        fits_reg = step.fits_reg;
        state = step.next_state;
        if converged {
            break;
        }
    }
    wald_test(mat, n_samples, x.n_cols, &fits_reg)
}

/// proda.py:545-558 — per-column imputation: missing entries are drawn from
/// `Normal(q10, sd/5)` over the observed values of that column.
fn impute_initial(mat: &[Vec<f64>], n_samples: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = SplitMix64::new(seed);
    let mut y: Vec<Vec<f64>> = mat.to_vec();
    for col in 0..n_samples {
        let obs: Vec<f64> = mat
            .iter()
            .filter_map(|row| {
                let v = row[col];
                if v.is_nan() {
                    None
                } else {
                    Some(v)
                }
            })
            .collect();
        let n_miss = mat.iter().filter(|row| row[col].is_nan()).count();
        if n_miss == 0 || obs.is_empty() {
            continue;
        }
        let q10 = percentile10(&obs);
        let sd5 = if obs.len() > 1 {
            sample_std(&obs) / 5.0
        } else {
            1.0
        };
        for row in y.iter_mut() {
            if row[col].is_nan() {
                row[col] = q10 + sd5 * rng.next_gaussian();
            }
        }
    }
    y
}

/// One EM iteration outcome (proda.py:569-658).
struct EmStep {
    err: f64,
    fits_reg: Vec<PdLmFit>,
    next_state: EmState,
}

/// proda.py:569-658 — a single EM iteration: unregularized + regularized fits,
/// re-estimate the priors / dropout curves, and compute the change `err`.
fn em_step(mat: &[Vec<f64>], x: &Design, state: &EmState) -> EmStep {
    let n_prot = mat.len();
    let priors = Priors {
        mu0: state.mu0,
        sigma20: state.sigma20,
        tau20: state.tau20,
        df0: state.df0,
    };

    let fits_unreg: Vec<PdLmFit> = (0..n_prot)
        .map(|i| pd_lm_fit(&mat[i], x, &state.rho, &state.zeta, None))
        .collect();
    let fits_reg: Vec<PdLmFit> = (0..n_prot)
        .map(|i| pd_lm_fit(&mat[i], x, &state.rho, &state.zeta, Some(priors)))
        .collect();

    let unreg = extract(&fits_unreg, x, &state.pred, &state.pred_var);
    let reg = extract(&fits_reg, x, &state.pred, &state.pred_var);

    let (mu0_new, sigma20_new) =
        location_prior(&reg.pred, &reg.pred_var, Some(&unreg.pred), Some(&unreg.pred_var));
    let (rho_new, zeta_new) = fit_dropout_curves(mat, x, &reg.pred, &reg.pred_var);

    // v_mask = (s2_unreg > 0) & (df_unreg > 0) (proda.py:629).
    let mut s2_v = Vec::new();
    let mut df_v = Vec::new();
    for (s2, df) in unreg.s2.iter().zip(&unreg.df) {
        if *s2 > 0.0 && *df > 0.0 {
            s2_v.push(*s2);
            df_v.push(*df);
        }
    }
    let (tau20_new, df0_new) = variance_prior(&s2_v, &df_v);

    let err = em_error(state, mu0_new, sigma20_new, &rho_new, &zeta_new, tau20_new, df0_new);

    EmStep {
        err,
        fits_reg,
        next_state: EmState {
            mu0: mu0_new,
            sigma20: sigma20_new,
            tau20: tau20_new,
            df0: df0_new,
            rho: rho_new,
            zeta: zeta_new,
            pred: reg.pred,
            pred_var: reg.pred_var,
        },
    }
}

/// proda.py:603-620 `_extract`: build Pred / Pred_var from a fit list, falling
/// back to the previous values for proteins whose coef is NaN.
fn extract(
    fits: &[PdLmFit],
    x: &Design,
    prev_pred: &[Vec<f64>],
    prev_pred_var: &[Vec<f64>],
) -> Extracted {
    let n_samples = x.n_rows;
    let mut pred = Vec::with_capacity(fits.len());
    let mut pred_var = Vec::with_capacity(fits.len());
    let mut s2 = vec![0.0; fits.len()];
    let mut df = vec![0.0; fits.len()];
    for (idx, fit) in fits.iter().enumerate() {
        if fit.coef.iter().any(|b| b.is_nan()) {
            pred.push(prev_pred[idx].clone());
            pred_var.push(prev_pred_var[idx].clone());
            continue;
        }
        pred.push(design_times(x, &fit.coef));
        let pv: Vec<f64> = (0..n_samples)
            .map(|j| {
                let xj = x.row(j);
                let v = mat_vec(&fit.coef_var, xj);
                dot(xj, &v)
            })
            .collect();
        pred_var.push(pv);
        s2[idx] = fit.s2;
        df[idx] = fit.df;
    }
    Extracted {
        pred,
        pred_var,
        s2,
        df,
    }
}

/// proda.py:632-651 — the aggregate parameter-change error driving convergence.
/// Each term is the squared difference of nan-mean-reduced parameter vectors;
/// scalar parameters reduce to a plain squared difference.
fn em_error(
    state: &EmState,
    mu0_new: f64,
    sigma20_new: f64,
    rho_new: &[f64],
    zeta_new: &[f64],
    tau20_new: f64,
    df0_new: f64,
) -> f64 {
    let zi_new = zeta_to_inv(zeta_new);
    let zi_old = zeta_to_inv(&state.zeta);
    let df0_term_new = if df0_new > 0.0 { 1.0 / df0_new } else { 0.0 };
    let df0_term_old = if state.df0 > 0.0 { 1.0 / state.df0 } else { 0.0 };

    sq_diff(mu0_new, state.mu0)
        + sq_diff(sigma20_new, state.sigma20)
        + sq_diff(nanmean_over_len(rho_new), nanmean_over_len(&state.rho))
        + sq_diff(nanmean_over_len(&zi_new), nanmean_over_len(&zi_old))
        + sq_diff(tau20_new, state.tau20)
        + sq_diff(df0_term_new, df0_term_old)
}

/// proda.py:632-633 `where(abs(zeta) > 1e-100, 1/zeta, nan)`.
fn zeta_to_inv(zeta: &[f64]) -> Vec<f64> {
    zeta.iter()
        .map(|z| {
            if z.abs() > 1e-100 {
                1.0 / z
            } else {
                f64::NAN
            }
        })
        .collect()
}

/// proda.py:636-637 `nansum(v) / max(len(v), 1)` — note this divides by the FULL
/// length (NaN entries contribute 0 to the sum but still count in the length).
fn nanmean_over_len(v: &[f64]) -> f64 {
    let sum: f64 = v.iter().filter(|x| x.is_finite()).sum();
    sum / (v.len().max(1)) as f64
}

fn sq_diff(a: f64, b: f64) -> f64 {
    (a - b) * (a - b)
}

/// proda.py:662-696 — per-protein Wald test on the final regularized fits.
/// `contrast = [1, -1]`, two-sided p from the Student-t with the per-protein
/// empirical-Bayes-moderated df (`fit.df`), not the naive `n_samples - p`.
fn wald_test(mat: &[Vec<f64>], n_samples: usize, p: usize, fits: &[PdLmFit]) -> Vec<ProdaStat> {
    let fallback_df = n_samples.saturating_sub(p).max(1) as f64;
    let contrast = [1.0, -1.0];
    let mut out = Vec::with_capacity(fits.len());
    for fit in fits {
        let beta = &fit.coef;
        let has_nan = beta.iter().any(|b| b.is_nan());
        let log2fc = if has_nan {
            f64::NAN
        } else {
            contrast[0] * beta[0] + contrast[1] * beta[1]
        };
        // se_fc = sqrt(contrast @ diag(se^2) @ contrast) (proda.py:671-672).
        let se2: Vec<f64> = fit.se_coef.iter().map(|s| s * s).collect();
        let se_fc = (contrast[0] * contrast[0] * se2.first().copied().unwrap_or(f64::NAN)
            + contrast[1] * contrast[1] * se2.get(1).copied().unwrap_or(f64::NAN))
        .sqrt();

        // Moderated df the fit already computed; df -> inf is proDA's normal limit.
        let moderated_df = if !fit.df.is_finite() {
            1e6
        } else if fit.df <= 0.0 {
            fallback_df
        } else {
            fit.df
        };

        let (t_stat, p_value) = if se_fc < 1e-10 || se_fc.is_nan() || log2fc.is_nan() {
            (0.0, 1.0)
        } else {
            let t = log2fc / se_fc;
            (t, 2.0 * student_t_sf(t.abs(), moderated_df))
        };
        out.push(ProdaStat {
            log2_fold_change: log2fc,
            t_stat,
            p_value,
        });
    }
    let _ = mat;
    out
}

// --- statistics helpers for imputation ---

/// `np.percentile(x, 10)` (linear interpolation, the numpy default).
fn percentile10(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n == 1 {
        return v[0];
    }
    // numpy: rank = (n - 1) * q, with q = 0.10.
    let rank = (n as f64 - 1.0) * 0.10;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        v[lo]
    } else {
        let frac = rank - lo as f64;
        v[lo] * (1.0 - frac) + v[hi] * frac
    }
}

/// `np.std(x, ddof=1)` (sample standard deviation).
fn sample_std(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
    var.sqrt()
}

/// Two-group proDA differential expression on a log2 protein x sample matrix.
///
/// `rows[i]` labels `proteins[i]`; the first `n_a` columns are condition A, the
/// next `n_b` condition B. Returns per-protein [`ProdaStat`] in INPUT order
/// (the public DE wrapper applies the testability filter and BH adjustment).
pub(crate) fn proda_run(rows: &[&[f64]], n_a: usize, n_b: usize) -> Vec<ProdaStat> {
    let mat: Vec<Vec<f64>> = rows.iter().map(|r| r.to_vec()).collect();
    run_proda(&mat, n_a, n_b)
}
