"""
Probabilistic dropout differential expression analysis (pure Python).

Reimplements the proDA model: sigmoid dropout curve estimation, EM-based
parameter estimation with empirical Bayes priors, and Wald test for
differential abundance, without rpy2 or R dependencies.

Reference
---------
Ahlmann-Eltze C, Anders S. proDA: Probabilistic Dropout Analysis for
Identifying Differentially Abundant Proteins in Label-Free Quantitative Mass
Spectrometry. *bioRxiv*. 2019. doi:10.1101/661496
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import f as f_dist
from scipy.stats import norm as sp_norm
from scipy.stats import t as t_dist

from mokume.analysis._helpers import (
    bh_adjust,
    filter_testable,
    per_group_summary,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.proda")


def _invprobit(x, rho, zeta, log=False, oneminus=False):
    """Dropout / observation probability — matches R proDA:::invprobit."""
    x, rho, zeta = (np.asarray(v, float) for v in (x, rho, zeta))
    if np.nansum(np.sign(zeta)) < 0:
        z = (x - rho) / np.maximum(np.abs(zeta), 1e-100)
        return (
            (sp_norm.logcdf if log else sp_norm.cdf)(z)
            if oneminus
            else (sp_norm.logsf if log else sp_norm.sf)(z)
        )
    z = (x - rho) / np.maximum(np.abs(zeta), 1e-100)
    return (
        (sp_norm.logsf if log else sp_norm.sf)(z)
        if oneminus
        else (sp_norm.logcdf if log else sp_norm.cdf)(z)
    )


def _trimmed_mean(a, trim=0.2):
    """Trimmed mean matching R mean(x, trim=0.2)."""
    a = np.asarray(a, float).ravel()
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return 0.0
    n = len(a)
    lo = int(np.floor(n * trim))
    hi = n - lo
    return float(np.mean(np.sort(a)[lo:hi]))


def _fit_dropout_curves(Y, X, Pred, Pred_var):
    """Per-sample dropout curves — matches R proDA:::dropout_curves."""
    n_samples = X.shape[0]
    mu0 = _trimmed_mean(Pred, 0.2)
    diff2 = (mu0 - Pred) ** 2
    sigma20 = _trimmed_mean(diff2, 0.2)
    if sigma20 == 0:
        sigma20 = 5.0

    rho = np.full(n_samples, np.nan)
    zeta = np.full(n_samples, np.nan)

    for col in range(n_samples):
        y_col = Y[:, col]
        obs_mask = ~np.isnan(y_col)
        if obs_mask.all():
            continue
        yo = y_col[obs_mask]
        predm = Pred[~obs_mask, col]
        pvm = np.maximum(Pred_var[~obs_mask, col], 1e-10)

        def _neg_ll(par, _yo=yo, _predm=predm, _pvm=pvm):
            r, zinv = par
            if zinv >= 0:
                return 1e10
            z = 1.0 / zinv
            ll = sp_norm.logpdf(r, mu0, np.sqrt(sigma20))
            ll += min(np.log(abs(zinv)), np.log(10000))
            ll += np.sum(_invprobit(_yo, r, z, log=True, oneminus=True))
            zs = np.sign(zinv) * np.sqrt(1.0 / zinv**2 + _pvm)
            ll += np.sum(_invprobit(_predm, r, zs, log=True))
            return -ll

        x0 = np.array([mu0, -1.0 / np.sqrt(sigma20)])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = minimize(
                _neg_ll,
                x0,
                method="Nelder-Mead",
                options={"maxiter": 1000, "xatol": 1e-8, "fatol": 1e-8},
            )
        rho[col] = res.x[0]
        zinv = res.x[1] if abs(res.x[1]) <= 10000 else -10000
        zeta[col] = 1.0 / zinv
    return rho, zeta


def _location_prior(X, Pred_reg, Pred_var_reg, Pred_unreg=None, Pred_var_unreg=None):
    """Estimate location prior — matches R proDA:::location_prior."""
    if Pred_unreg is None:
        Pred_unreg = Pred_reg
    if Pred_var_unreg is None:
        Pred_var_unreg = Pred_var_reg
    mu0 = _trimmed_mean(Pred_reg, 0.2)
    pr = Pred_reg.ravel()
    pu = Pred_unreg.ravel()
    pvr = Pred_var_reg.ravel()
    pvu = Pred_var_unreg.ravel()
    pred_flat = np.where(np.isnan(pu), pr, pu)
    pv_flat = np.where(np.isnan(pvu), pvr, pvu)
    above = pr > mu0
    p = pred_flat[above] - mu0
    pv = pv_flat[above]

    def _obj(A):
        w = 1.0 / (2.0 * (A + pv) ** 2)
        return float(np.sum((p**2 - pv) * w) / np.sum(w) - A)

    lo, hi = 0.0, 1000.0
    try:
        if np.sign(_obj(lo)) != np.sign(_obj(hi)):
            from scipy.optimize import brentq

            sigma20 = brentq(_obj, lo, hi)
        else:
            sigma20 = float(np.sum(p**2) / max(len(p), 1))
    except Exception:
        sigma20 = float(np.sum(p**2) / max(len(p), 1))
    return mu0, max(sigma20, 1e-10)


def _variance_prior(s2_arr, df_arr):
    """Estimate empirical Bayes variance prior — matches R proDA:::variance_prior."""
    s2 = np.asarray(s2_arr, dtype=float)
    df = np.asarray(df_arr, dtype=float)
    ok = np.isfinite(s2) & np.isfinite(df) & (s2 > 0) & (df > 0)
    s2, df = s2[ok], df[ok]
    if len(s2) < 2:
        return 1.0, 1.0

    def _neg_ll(par):
        tau, dfinv = par
        if tau <= 0 or dfinv <= 0:
            return 1e10
        df0 = 1.0 / dfinv
        ll = np.sum(f_dist.logpdf(s2 / tau, df, df0) - np.log(tau))
        return -ll

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = minimize(
            _neg_ll,
            [1.0, 1.0],
            method="Nelder-Mead",
            options={"maxiter": 10000, "xatol": 1e-8, "fatol": 1e-8},
        )
    tau20 = max(res.x[0], 1e-10)
    df0 = 1.0 / max(res.x[1], 1e-10)
    return tau20, df0


def _pd_lm_fit(y, X, rho, zeta, mu0=None, sigma20=None, tau20=None, df0=None, loc_df=3):
    """Per-protein MLE — matches R proDA:::pd_lm.fit with priors."""
    obs = ~np.isnan(y)
    n_obs = int(obs.sum())
    n, p = len(y), X.shape[1]
    mod_loc = mu0 is not None and sigma20 is not None
    mod_var = tau20 is not None and df0 is not None

    if n_obs == 0:
        return {
            "coef": np.full(p, np.nan),
            "s2": np.nan,
            "se_coef": np.full(p, np.nan),
            "coef_var": np.full((p, p), np.nan),
            "df": np.nan,
            "n_obs": 0,
        }

    Xo, yo = X[obs], y[obs]
    Xm = X[~obs]
    rho_m = np.asarray(rho)[~obs] if not obs.all() else np.array([])
    zeta_m = np.asarray(zeta)[~obs] if not obs.all() else np.array([])

    nan_dc = np.isnan(rho_m)
    no_dropout = obs.all() or nan_dc.all()

    if no_dropout and not mod_var and not mod_loc:
        beta = np.linalg.lstsq(Xo, yo, rcond=None)[0]
        resid = yo - Xo @ beta
        s2 = float(np.sum(resid**2)) / max(n - p, 1)
        try:
            XtXi = np.linalg.inv(Xo.T @ Xo)
        except np.linalg.LinAlgError:
            XtXi = np.linalg.pinv(Xo.T @ Xo)
        Vcov = XtXi * s2
        se = np.sqrt(np.maximum(np.diag(Vcov), 0))
        return {
            "coef": beta,
            "s2": s2,
            "se_coef": se,
            "coef_var": Vcov,
            "df": float(max(n - p, 1)),
            "n_obs": n_obs,
        }

    has_dropout = len(Xm) > 0 and not nan_dc.all()
    rm = rho_m[~nan_dc] if has_dropout else np.array([])
    zm = zeta_m[~nan_dc] if has_dropout else np.array([])
    Xm_d = Xm[~nan_dc] if has_dropout else np.empty((0, p))

    if mod_loc:
        beta_init = np.zeros(p)
        beta_init[0] = mu0
    else:
        beta_init = np.full(p, float(np.mean(yo)))
    s2_init = (df0 * tau20 / (df0 + 2)) if mod_var else 1.0

    def _neg_ll(par):
        beta = par[:p]
        sigma2 = max(par[p], 1e-100)
        ll = 0.0
        if mod_loc:
            ll += np.sum(
                t_dist.logpdf(X @ beta, loc_df, loc=mu0, scale=np.sqrt(sigma20))
            )
        if mod_var:
            ll += (
                -(df0 / 2 + 1) * np.log(sigma2)
                - df0 * tau20 / (2 * sigma2)
                + np.log(sigma2)
            )
        ll += float(np.sum(sp_norm.logpdf(yo, Xo @ beta, np.sqrt(sigma2))))
        if has_dropout:
            mm = Xm_d @ beta
            zs = zm * np.sqrt(1.0 + sigma2 / zm**2)
            ll += float(np.sum(_invprobit(mm, rm, zs, log=True)))
        return -ll

    def _neg_ll_log(par_log):
        par = par_log.copy()
        par[p] = np.exp(np.clip(par_log[p], -50, 50))
        return _neg_ll(par)

    def _grad_log(par_log):
        beta = par_log[:p]
        sigma2 = np.exp(np.clip(par_log[p], -50, 50))
        g_b = np.zeros(p)
        g_s = 0.0
        resid = yo - Xo @ beta
        g_b -= (1.0 / sigma2) * (Xo.T @ resid)
        g_s += len(yo) / (2 * sigma2) - np.sum(resid**2) / (2 * sigma2**2)
        if mod_loc:
            diff = X @ beta - mu0
            g_b += X.T @ ((loc_df + 1) * diff / (loc_df * sigma20 + diff**2))
        if mod_var:
            g_s += df0 / (2 * sigma2) - df0 * tau20 / (2 * sigma2**2)
        if has_dropout:
            mm = Xm_d @ beta
            zs2 = zm**2 + sigma2
            zs_abs = np.sqrt(zs2)
            z = (mm - rm) / zs_abs
            ratio = sp_norm.pdf(z) / np.maximum(sp_norm.sf(z), 1e-300)
            g_b += Xm_d.T @ (ratio / zs_abs)
            g_s -= float(np.sum(ratio * z / (2 * zs2)))
        g_s *= sigma2
        return np.concatenate([g_b, [g_s]])

    _nan_result = {
        "coef": np.full(p, np.nan),
        "s2": np.nan,
        "se_coef": np.full(p, np.nan),
        "coef_var": np.full((p, p), np.nan),
        "df": np.nan,
        "n_obs": n_obs,
    }
    x0_log = np.concatenate([beta_init, [np.log(max(s2_init, 1e-100))]])
    bounds = [(None, None)] * p + [(-50, 50)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res_log = minimize(
                _neg_ll_log,
                x0_log,
                jac=_grad_log,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 1000, "ftol": 1e-15, "gtol": 1e-12},
            )
            x_bfgs = res_log.x.copy()
            x_bfgs[p] = np.clip(x_bfgs[p], -50, 50)
            res_log = minimize(
                _neg_ll_log,
                x_bfgs,
                jac=_grad_log,
                method="BFGS",
                options={"maxiter": 200, "gtol": 1e-14},
            )
        except Exception:
            return _nan_result

    xopt = res_log.x.copy()
    xopt[p] = np.clip(xopt[p], -50, 50)
    if not np.all(np.isfinite(xopt)):
        return _nan_result
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _nr in range(5):
            g = _grad_log(xopt)
            if not np.all(np.isfinite(g)) or np.max(np.abs(g)) < 1e-14:
                break
            eps_h = 1e-5
            dim = len(xopt)
            H = np.zeros((dim, dim))
            for j in range(dim):
                e = np.zeros(dim)
                e[j] = eps_h
                H[:, j] = (_grad_log(xopt + e) - _grad_log(xopt - e)) / (2 * eps_h)
            H = 0.5 * (H + H.T)
            try:
                step = np.linalg.solve(H, -g)
                if np.max(np.abs(step)) < 1e-14:
                    break
                xopt += step
                xopt[p] = np.clip(xopt[p], -50, 50)
            except np.linalg.LinAlgError:
                break

    class _Res:
        pass

    res = _Res()
    res.x = xopt.copy()
    res.x[p] = np.exp(xopt[p])

    beta = res.x[:p]
    fit_sigma2 = max(res.x[p], 1e-100)

    try:
        resid_h = yo - Xo @ beta
        dss_o = float(np.sum((fit_sigma2 - 2.0 * resid_h**2) / (2.0 * fit_sigma2**3)))
        dss_m = 0.0
        if has_dropout:
            mm_h = Xm_d @ beta
            zs2_h = zm**2 + fit_sigma2
            zs_abs_h = np.sqrt(zs2_h)
            z_h = (mm_h - rm) / zs_abs_h
            ratio_h = sp_norm.pdf(z_h) / np.maximum(sp_norm.sf(z_h), 1e-300)
            imr_h = -ratio_h / zs_abs_h
            diff_h = mm_h - rm
            dss_m = float(
                np.sum(
                    diff_h
                    / (4.0 * zs2_h**2)
                    * imr_h
                    * (3.0 - diff_h * imr_h - diff_h**2 / zs2_h)
                )
            )
        dss_p = 0.0
        if mod_var:
            dss_p = (
                (1.0 + df0 / 2.0) / fit_sigma2**2
                - df0 * tau20 / fit_sigma2**3
                - 1.0 / fit_sigma2**2
            )
        neg_hess_s2 = -(dss_p + dss_o + dss_m)
        fit_sigma2_var = 1.0 / max(neg_hess_s2, 1e-100)
        if fit_sigma2_var < 0:
            raise ValueError("negative sigma2 variance")
    except Exception:
        fit_sigma2_var = 2.0 * fit_sigma2**2 / max(n_obs, 1)

    n_approx_raw = 2.0 * fit_sigma2**2 / max(fit_sigma2_var, 1e-100)
    rss_approx = 2.0 * fit_sigma2**3 / max(fit_sigma2_var, 1e-100)
    s2_approx = rss_approx / max(n_approx_raw - p, 0.001)
    if s2_approx < 0 or n_approx_raw <= p:
        df_approx = 0.001
        n_approx_raw = df_approx + p
        s2_approx = np.sqrt(fit_sigma2_var * n_approx_raw**3 / (2.0 * df_approx**2))
    elif mod_var:
        n_approx_raw -= df0
        df_approx = n_approx_raw - p + df0
        if df_approx > 30 * n and df_approx > 100:
            df_approx = np.inf
            s2_approx = tau20
    else:
        df_approx = n_approx_raw - p

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            zetastar2 = zm**2 + fit_sigma2

            def _analytic_hess_bb(sigma2_val):
                h_o = -(Xo.T @ Xo) / sigma2_val
                h_p = np.zeros((p, p))
                if mod_loc:
                    d_loc = X @ beta - mu0
                    t_pf = (loc_df * sigma20 - d_loc**2) / (
                        loc_df * sigma20 + d_loc**2
                    ) ** 2
                    h_p = -(loc_df + 1) * (X.T @ np.diag(t_pf) @ X)
                h_m = np.zeros((p, p))
                if has_dropout:
                    mm_c = Xm_d @ beta
                    zs_abs_c = np.sqrt(zetastar2)
                    z_c = (mm_c - rm) / zs_abs_c
                    rat = sp_norm.pdf(z_c) / np.maximum(sp_norm.sf(z_c), 1e-300)
                    imr_c = -rat / zs_abs_c
                    w_m = imr_c**2 + (mm_c - rm) / zetastar2 * imr_c
                    h_m = -(Xm_d.T @ np.diag(w_m) @ Xm_d)
                return h_p + h_o + h_m

            hess_fit = _analytic_hess_bb(fit_sigma2)
            Var_coef_fit = np.linalg.inv(-hess_fit)

            par_fit = res.x.copy()
            out_factor = 8.0
            offset = _neg_ll(par_fit)
            cf_diag = np.ones(p)
            for idx in range(p):
                if np.any(np.isnan(beta)):
                    continue
                others = [j for j in range(p) if j != idx]
                V_sub = Var_coef_fit[np.ix_(others, others)]
                V_cross = Var_coef_fit[idx, others]
                try:
                    V_sub_inv = np.linalg.inv(V_sub)
                    vc_corr = Var_coef_fit[idx, idx] - V_cross @ V_sub_inv @ V_cross
                except np.linalg.LinAlgError:
                    vc_corr = Var_coef_fit[idx, idx]
                if vc_corr <= 0 or np.isnan(vc_corr):
                    continue
                shift_vec = np.zeros(p + 1)
                shift_vec[idx] = np.sqrt(out_factor * vc_corr)
                par_shift = par_fit.copy()
                par_shift[:p] += shift_vec[:p]
                diff_ll = abs(_neg_ll(par_shift) - offset)
                if diff_ll > 0:
                    cf_diag[idx] = (diff_ll / (out_factor / 2.0)) ** (-0.5)
            CF = np.diag(cf_diag)

            hess_s2a = _analytic_hess_bb(s2_approx)
            neg_hess_s2a = -hess_s2a
            if np.any(np.diag(neg_hess_s2a) < 0):
                neg_hess_s2a = np.zeros((p, p))
            Var_coef_s2a = (
                np.linalg.inv(neg_hess_s2a)
                if np.all(np.diag(neg_hess_s2a) > 0)
                else np.full((p, p), np.inf)
            )
            Var_coef_unbiased = CF @ Var_coef_s2a @ CF

            large_var = np.diag(Var_coef_fit) > 1e6
            for idx in np.where(large_var)[0]:
                Var_coef_unbiased[idx, idx] = np.inf
            se = np.sqrt(np.maximum(np.diag(Var_coef_unbiased), 0))
            coef_var = Var_coef_unbiased
        except Exception:
            se = np.full(p, np.nan)
            coef_var = np.full((p, p), np.nan)

    return {
        "coef": beta,
        "s2": float(s2_approx),
        "se_coef": se,
        "coef_var": coef_var,
        "df": float(max(df_approx, 0.001)),
        "n_obs": n_obs,
    }


def _ols_per_protein(mat, X):
    """OLS fit per protein, return Pred, Pred_var, s2, df arrays."""
    n_prot, n_samples = mat.shape
    p = X.shape[1]
    try:
        XtXi = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        XtXi = np.linalg.pinv(X.T @ X)
    Pred = np.zeros_like(mat)
    Pred_var = np.zeros_like(mat)
    s2_arr = np.zeros(n_prot)
    df_arr = np.zeros(n_prot)
    for i in range(n_prot):
        beta_i = np.linalg.lstsq(X, mat[i], rcond=None)[0]
        pred_i = X @ beta_i
        Pred[i] = pred_i
        resid = mat[i] - pred_i
        s2_i = float(np.sum(resid**2)) / max(n_samples - p, 1)
        s2_arr[i] = s2_i
        df_arr[i] = max(n_samples - p, 1)
        for j in range(n_samples):
            Pred_var[i, j] = float(X[j] @ XtXi @ X[j]) * s2_i
    return Pred, Pred_var, s2_arr, df_arr


def run_proda(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    cond_a: str,
    cond_b: str,
    seed: int = 42,
) -> pd.DataFrame:
    """Run proDA differential expression (pure Python implementation).

    Parameters
    ----------
    seed : int
        RNG seed for the initial missing-value imputation step before EM.
        Different seeds yield slightly different EM trajectories.
    """
    sample_order = list(samples_a) + list(samples_b)
    full_matrix = log2_matrix[sample_order].dropna(how="all")
    if full_matrix.empty:
        logger.warning("No proteins passed DE filtering for proDA")
        return pd.DataFrame()

    mat = full_matrix.values
    gene_names = list(full_matrix.index)
    n_prot, n_samples = mat.shape

    n_a = len(samples_a)
    X = np.zeros((n_samples, 2))
    X[:n_a, 0] = 1.0
    X[n_a:, 1] = 1.0

    rng = np.random.RandomState(seed)
    Y_init = mat.copy()
    for col in range(n_samples):
        col_vals = mat[:, col]
        col_miss = np.isnan(col_vals)
        n_miss = int(col_miss.sum())
        if n_miss == 0:
            continue
        col_obs = col_vals[~col_miss]
        if len(col_obs) == 0:
            continue
        q10_col = float(np.percentile(col_obs, 10))
        sd5_col = float(np.std(col_obs, ddof=1)) / 5.0 if len(col_obs) > 1 else 1.0
        Y_init[col_miss, col] = rng.normal(q10_col, sd5_col, size=n_miss)

    Pred, Pred_var, s2_init, df_init = _ols_per_protein(Y_init, X)

    mu0, sigma20 = _location_prior(X, Pred, Pred_var)
    rho, zeta = _fit_dropout_curves(mat, X, Pred, Pred_var)
    tau20, df0 = _variance_prior(s2_init, df_init)

    fits = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _iter in range(20):
            _nan_fit = {
                "coef": np.full(X.shape[1], np.nan),
                "s2": np.nan,
                "se_coef": np.full(X.shape[1], np.nan),
                "coef_var": np.full((X.shape[1], X.shape[1]), np.nan),
                "df": np.nan,
                "n_obs": 0,
            }
            fits_unreg = []
            for i in range(n_prot):
                try:
                    fit_u = _pd_lm_fit(mat[i], X, rho, zeta)
                except Exception:
                    fit_u = _nan_fit.copy()
                fits_unreg.append(fit_u)

            fits_reg = []
            for i in range(n_prot):
                try:
                    fit_r = _pd_lm_fit(
                        mat[i],
                        X,
                        rho,
                        zeta,
                        mu0=mu0,
                        sigma20=sigma20,
                        tau20=tau20,
                        df0=df0,
                    )
                except Exception:
                    fit_r = _nan_fit.copy()
                fits_reg.append(fit_r)

            def _extract(fit_list):
                P = np.zeros_like(mat)
                PV = np.zeros_like(mat)
                s2a = np.zeros(n_prot)
                dfa = np.zeros(n_prot)
                for idx, f in enumerate(fit_list):
                    b = f["coef"]
                    if np.any(np.isnan(b)):
                        P[idx] = Pred[idx]
                        PV[idx] = Pred_var[idx]
                        continue
                    P[idx] = X @ b
                    Vcov = f["coef_var"]
                    for jj in range(n_samples):
                        PV[idx, jj] = float(X[jj] @ Vcov @ X[jj])
                    s2a[idx] = f["s2"]
                    dfa[idx] = f["df"]
                return P, PV, s2a, dfa

            Pred_unreg, PV_unreg, s2_unreg, df_unreg = _extract(fits_unreg)
            Pred_reg, PV_reg, _, _ = _extract(fits_reg)

            mu0_new, sigma20_new = _location_prior(
                X, Pred_reg, PV_reg, Pred_unreg, PV_unreg
            )
            rho_new, zeta_new = _fit_dropout_curves(mat, X, Pred_reg, PV_reg)
            v_mask = (s2_unreg > 0) & (df_unreg > 0)
            tau20_new, df0_new = _variance_prior(s2_unreg[v_mask], df_unreg[v_mask])

            zi_new = np.where(np.abs(zeta_new) > 1e-100, 1.0 / zeta_new, np.nan)
            zi_old = np.where(np.abs(zeta) > 1e-100, 1.0 / zeta, np.nan)
            err = sum(
                (
                    np.nansum(v_new) / max(len(np.atleast_1d(v_new)), 1)
                    - np.nansum(v_old) / max(len(np.atleast_1d(v_old)), 1)
                )
                ** 2
                for v_new, v_old in [
                    ([mu0_new], [mu0]),
                    ([sigma20_new], [sigma20]),
                    (rho_new, rho),
                    (zi_new, zi_old),
                    ([tau20_new], [tau20]),
                    (
                        [1.0 / df0_new if df0_new > 0 else 0],
                        [1.0 / df0 if df0 > 0 else 0],
                    ),
                ]
            )
            Pred, Pred_var = Pred_reg, PV_reg
            rho, zeta = rho_new, zeta_new
            mu0, sigma20 = mu0_new, sigma20_new
            tau20, df0 = tau20_new, df0_new
            fits = fits_reg
            if err < 1e-3:
                break

    logger.info("proDA dropout: rho=%s, zeta=%s", np.round(rho, 3), np.round(zeta, 3))

    contrast = np.array([1.0, -1.0])
    wald_df = n_samples - X.shape[1]
    results = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, gene in enumerate(gene_names):
            fit = fits[i]
            beta = fit["coef"]
            log2fc = float(contrast @ beta) if not np.any(np.isnan(beta)) else np.nan
            se2 = fit["se_coef"] ** 2
            se_fc = float(np.sqrt(contrast @ np.diag(se2) @ contrast))

            if se_fc < 1e-10 or np.isnan(se_fc) or np.isnan(log2fc):
                t_stat, pvalue = 0.0, 1.0
            else:
                t_stat = log2fc / se_fc
                pvalue = 2.0 * float(t_dist.sf(abs(t_stat), wald_df))

            results.append(
                {
                    "ProteinName": gene,
                    "log2FC": log2fc,
                    "se": se_fc,
                    "t_stat": t_stat,
                    "pvalue": pvalue,
                    "n_obs": fit["n_obs"],
                    "avg_abundance": float(np.nanmean(mat[i])),
                    "df": wald_df,
                }
            )

    raw = pd.DataFrame(results)
    testable = filter_testable(full_matrix, samples_a, samples_b)
    raw = raw[raw["ProteinName"].isin(testable.index)]
    summary = per_group_summary(full_matrix, samples_a, samples_b, cond_a, cond_b)
    merged = raw.merge(summary, on="ProteinName", how="left")
    merged["adj_pvalue"] = bh_adjust(merged["pvalue"].values)

    columns = [
        "ProteinName",
        "log2FC",
        "pvalue",
        "adj_pvalue",
        "t_stat",
        f"mean_{cond_a}",
        f"mean_{cond_b}",
        "n_a",
        "n_b",
        "df",
        "se",
        "avg_abundance",
        "n_obs",
    ]
    available = [c for c in columns if c in merged.columns]
    return merged[available].sort_values("adj_pvalue").reset_index(drop=True)
