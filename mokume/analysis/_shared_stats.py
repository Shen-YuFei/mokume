import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.shared_stats")


def _trigamma(x):
    return polygamma(1, x)


def _tetragamma(x):
    return polygamma(2, x)


def _collect_protein_stats(
    log2_matrix: pd.DataFrame,
    sample_groups: tuple[list[str], list[str]],
    contrast: tuple[str, str],
) -> pd.DataFrame:
    samples_a, samples_b = sample_groups
    cond_a, cond_b = contrast
    rows = []
    for protein in log2_matrix.index:
        va = log2_matrix.loc[protein, samples_a].dropna().values
        vb = log2_matrix.loc[protein, samples_b].dropna().values
        if len(va) < 2 or len(vb) < 2:
            continue
        na, nb = len(va), len(vb)
        ma, mb = np.mean(va), np.mean(vb)
        ss = np.sum((va - ma) ** 2) + np.sum((vb - mb) ** 2)
        df_res = (na - 1) + (nb - 1)
        rows.append(
            {
                "ProteinName": protein,
                "log2FC": ma - mb,
                f"mean_{cond_a}": ma,
                f"mean_{cond_b}": mb,
                "n_a": na,
                "n_b": nb,
                "s2": ss / df_res,
                "df_residual": df_res,
            }
        )
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _fit_f_prior(
    s2: np.ndarray, df: np.ndarray,
) -> tuple[float, float]:
    valid = np.isfinite(s2) & (s2 > 0) & np.isfinite(df) & (df > 0)
    s2, df = s2[valid], df[valid]

    if len(s2) < 3:
        logger.warning("Too few proteins for eBayes estimation, using defaults")
        return 3.0, float(np.median(s2)) if len(s2) else 1.0

    unique_df = np.unique(df)
    if len(unique_df) == 1:
        d0, s0_sq = _fit_f_prior_balanced(s2, unique_df[0])
    else:
        d0, s0_sq = _fit_f_prior_unbalanced(s2, df)

    return max(d0, 0.01), max(s0_sq, 1e-15)


def _fit_f_prior_balanced(
    s2: np.ndarray, d: float,
) -> tuple[float, float]:
    z = np.log(s2)
    z_mean, z_var = np.mean(z), np.var(z, ddof=1)
    target = z_var - _trigamma(d / 2.0)
    d0 = 1e6 if target <= 0 else _solve_trigamma(target) * 2.0
    log_s0 = (
        z_mean
        - digamma(d / 2.0) + np.log(d / 2.0)
        + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    )
    return d0, np.exp(log_s0)


def _fit_f_prior_unbalanced(
    s2: np.ndarray, df: np.ndarray,
) -> tuple[float, float]:
    z = np.log(s2)
    e_z = digamma(df / 2.0) - np.log(df / 2.0)
    z_res = z - e_z
    target = np.var(z_res, ddof=1) - np.mean(_trigamma(df / 2.0))
    d0 = 1e6 if target <= 0 else _solve_trigamma(target) * 2.0
    log_s0 = np.mean(z_res) + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    return d0, np.exp(log_s0)


def _solve_trigamma(target: float) -> float:
    if target <= 0:
        return 1e10
    if target > 1e6:
        return 1.0 / target

    x = 0.5 + 1.0 / target
    for _ in range(50):
        fx = _trigamma(x) - target
        dfx = _tetragamma(x)
        if abs(dfx) < 1e-20:
            break
        step = fx / dfx
        while x - step <= 0:
            step /= 2.0
        x = x - step
        if abs(fx) < 1e-12:
            break

    return max(x, 1e-10)
