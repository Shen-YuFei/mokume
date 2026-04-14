"""
LimROTS: limma + ROTS bootstrap-optimized test statistic.

Combines limma's empirical Bayes moderated statistics with ROTS-style
bootstrap optimization of the test statistic d/(a1 + a2*s), yielding
superior sensitivity while controlling false positives.

Reference
---------
Anwar AM, Jeba A, Lahti L, Coffey E. LimROTS: a hybrid method
integrating empirical Bayes and reproducibility-optimized statistics
for robust differential expression analysis. Bioinformatics. 2025.
"""

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from mokume.analysis.differential_expression import (
    _collect_protein_stats,
    _fit_f_prior,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.limrots")


def _limma_moderated_ds(stats_df: pd.DataFrame):
    """Compute limma moderated d (|logFC|) and s (posterior SE)."""
    d0, s0_sq = _fit_f_prior(
        stats_df["s2"].values, stats_df["df_residual"].values,
    )
    df_post = stats_df["df_residual"].values + d0
    s2_post = (
        d0 * s0_sq
        + stats_df["df_residual"].values * stats_df["s2"].values
    ) / df_post
    stdev_unscaled = np.sqrt(
        1.0 / stats_df["n_a"].values + 1.0 / stats_df["n_b"].values,
    )
    d = np.abs(stats_df["log2FC"].values)
    s = np.sqrt(s2_post) * stdev_unscaled
    return d, s, d0, df_post


_SHARED: dict = {}


def _init_worker(shared_dict: dict) -> None:
    """Initialize worker process with shared data."""
    _SHARED.update(shared_dict)


def _worker_bootstrap(b: int):
    """Thin wrapper that reads shared data for ProcessPoolExecutor."""
    return _bootstrap_one(
        b,
        _SHARED["log2_mat"],
        _SHARED["samples_a"],
        _SHARED["samples_b"],
        _SHARED["all_samples"],
        _SHARED["boot_sa"],
        _SHARED["boot_sb"],
        _SHARED["perm_orders"],
        _SHARED["cond_a"],
        _SHARED["cond_b"],
        _SHARED["n_prot"],
        _SHARED["prot_to_idx"],
        _SHARED["n_a"],
    )


def _bootstrap_one(
    b, log2_mat, samples_a, samples_b, all_samples,
    boot_sa, boot_sb, perm_orders,
    cond_a, cond_b, n_prot, prot_to_idx, n_a,
):
    """Run one bootstrap + permutation iteration."""
    d_col = np.full(n_prot, np.nan)
    s_col = np.full(n_prot, np.nan)
    pd_col = np.full(n_prot, np.nan)
    ps_col = np.full(n_prot, np.nan)

    # Bootstrap
    try:
        sa_b = [samples_a[i] for i in boot_sa[b]]
        sb_b = [samples_b[i] for i in boot_sb[b]]
        bst = _collect_protein_stats(log2_mat, sa_b, sb_b, cond_a, cond_b)
        if not bst.empty:
            bd0, bs0 = _fit_f_prior(
                bst["s2"].values, bst["df_residual"].values,
            )
            dfp = bst["df_residual"].values + bd0
            s2p = (
                bd0 * bs0 + bst["df_residual"].values * bst["s2"].values
            ) / dfp
            su = np.sqrt(1.0 / bst["n_a"].values + 1.0 / bst["n_b"].values)
            for j, p in enumerate(bst["ProteinName"].values):
                if p in prot_to_idx:
                    idx = prot_to_idx[p]
                    d_col[idx] = np.abs(bst["log2FC"].values[j])
                    s_col[idx] = np.sqrt(s2p[j]) * su[j]
    except Exception as exc:
        logger.warning("Bootstrap %d failed: %s", b, exc)

    # Permutation
    try:
        po = perm_orders[b]
        pa = [all_samples[i] for i in po[:n_a]]
        pb = [all_samples[i] for i in po[n_a:]]
        pst = _collect_protein_stats(log2_mat, pa, pb, cond_a, cond_b)
        if not pst.empty:
            pd0, ps0 = _fit_f_prior(
                pst["s2"].values, pst["df_residual"].values,
            )
            pdfp = pst["df_residual"].values + pd0
            ps2p = (
                pd0 * ps0 + pst["df_residual"].values * pst["s2"].values
            ) / pdfp
            psu = np.sqrt(
                1.0 / pst["n_a"].values + 1.0 / pst["n_b"].values,
            )
            for j, p in enumerate(pst["ProteinName"].values):
                if p in prot_to_idx:
                    idx = prot_to_idx[p]
                    pd_col[idx] = np.abs(pst["log2FC"].values[j])
                    ps_col[idx] = np.sqrt(ps2p[j]) * psu[j]
    except Exception as exc:
        logger.warning("Permutation %d failed: %s", b, exc)

    return b, d_col, s_col, pd_col, ps_col


def _optimize_ssq(D, S, pD, pS, n_boot, n_prot):
    """Find optimal fudge factor via z-score maximization."""
    s_pos = S[S > 0]
    if len(s_pos) == 0:
        return 0.01
    ssq_grid = np.quantile(s_pos, np.linspace(0, 0.95, 20))
    n_tops = [10, 25, 50, 100, min(200, n_prot // 4)]
    best_z = -np.inf
    best_ssq = ssq_grid[0]

    for ssq in ssq_grid:
        for n_top in n_tops:
            overlaps = []
            perm_overlaps = []
            for i in range(0, n_boot - 1, 2):
                t1 = D[:, i] / (ssq + S[:, i])
                t2 = D[:, i + 1] / (ssq + S[:, i + 1])
                top1 = set(np.argsort(-t1)[:n_top])
                top2 = set(np.argsort(-t2)[:n_top])
                overlaps.append(len(top1 & top2) / n_top)

                pt1 = pD[:, i] / (ssq + pS[:, i])
                pt2 = pD[:, i + 1] / (ssq + pS[:, i + 1])
                ptop1 = set(np.argsort(-pt1)[:n_top])
                ptop2 = set(np.argsort(-pt2)[:n_top])
                perm_overlaps.append(len(ptop1 & ptop2) / n_top)

            mean_ol = np.mean(overlaps)
            mean_pol = np.mean(perm_overlaps)
            sd_ol = max(np.std(overlaps, ddof=1), 1e-10)
            z = (mean_ol - mean_pol) / sd_ol
            if z > best_z:
                best_z = z
                best_ssq = ssq

    return best_ssq


def _permutation_fdr(t_final, perm_t, n_boot):
    """Compute FDR from permutation distribution."""
    n_prot = len(t_final)
    abs_t = np.abs(t_final)
    abs_perm = np.abs(perm_t)
    order = np.argsort(-abs_t)

    fdr = np.ones(n_prot)
    for rank_i, idx in enumerate(order):
        thresh = abs_t[idx]
        n_above = rank_i + 1
        perm_counts = np.sum(abs_perm >= thresh, axis=0)
        fdr[idx] = np.median(perm_counts) / n_above
    fdr = np.minimum(fdr, 1.0)

    # Monotonize
    for i in range(len(order) - 2, -1, -1):
        if i + 1 < len(order):
            fdr[order[i]] = min(fdr[order[i]], fdr[order[i + 1]])
    return fdr


def run_limrots(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    cond_a: str,
    cond_b: str,
    n_boot: int = 100,
    n_threads: int | None = None,
    seed: int = 42,
) -> pd.DataFrame:
    """LimROTS: limma + ROTS bootstrap-optimized DE test.

    Parameters
    ----------
    log2_matrix : pd.DataFrame
        Log2-transformed protein intensity matrix (proteins × samples).
    samples_a, samples_b : list[str]
        Sample column names for the two conditions.
    cond_a, cond_b : str
        Condition labels for the contrast.
    n_boot : int
        Number of bootstrap/permutation iterations.
    n_threads : int or None
        Number of parallel workers. ``None`` uses all available CPU cores.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        DE results with columns: ProteinName, log2FC, pvalue, adj_pvalue.
    """
    rng = np.random.default_rng(seed)

    # Step 1: Original limma stats
    stats_df = _collect_protein_stats(
        log2_matrix, samples_a, samples_b, cond_a, cond_b,
    )
    if stats_df.empty:
        return pd.DataFrame()

    d_orig, s_orig, d0, df_post = _limma_moderated_ds(stats_df)
    n_prot = len(d_orig)
    all_samples = samples_a + samples_b
    n_a = len(samples_a)
    prot_index = stats_df["ProteinName"].values
    prot_to_idx = {p: i for i, p in enumerate(prot_index)}

    # Pre-generate random indices
    boot_sa = [
        rng.choice(len(samples_a), len(samples_a), replace=True)
        for _ in range(n_boot)
    ]
    boot_sb = [
        rng.choice(len(samples_b), len(samples_b), replace=True)
        for _ in range(n_boot)
    ]
    perm_orders = [
        rng.permutation(len(all_samples)) for _ in range(n_boot)
    ]

    # Step 2-3: Parallel bootstrap + permutation
    D = np.full((n_prot, n_boot), np.nan)
    S = np.full((n_prot, n_boot), np.nan)
    pD = np.full((n_prot, n_boot), np.nan)
    pS = np.full((n_prot, n_boot), np.nan)

    if n_threads is None:
        n_threads = os.cpu_count() or 4
    logger.info("LimROTS: %d bootstrap iterations (%d workers)", n_boot, n_threads)

    shared_dict = {
        "log2_mat": log2_matrix,
        "samples_a": samples_a,
        "samples_b": samples_b,
        "all_samples": all_samples,
        "boot_sa": boot_sa,
        "boot_sb": boot_sb,
        "perm_orders": perm_orders,
        "cond_a": cond_a,
        "cond_b": cond_b,
        "n_prot": n_prot,
        "prot_to_idx": prot_to_idx,
        "n_a": n_a,
    }

    with ProcessPoolExecutor(
        max_workers=n_threads,
        initializer=_init_worker,
        initargs=(shared_dict,),
    ) as pool:
        for b, dc, sc, pdc, psc in pool.map(_worker_bootstrap, range(n_boot)):
            D[:, b] = dc
            S[:, b] = sc
            pD[:, b] = pdc
            pS[:, b] = psc

    D = np.nan_to_num(D)
    S = np.nan_to_num(S, nan=1.0)
    pD = np.nan_to_num(pD)
    pS = np.nan_to_num(pS, nan=1.0)

    # Step 4: Optimize fudge factor
    best_ssq = _optimize_ssq(D, S, pD, pS, n_boot, n_prot)
    logger.info("LimROTS: optimal s0 = %.4f", best_ssq)

    # Step 5: Final optimized test statistic → permutation FDR
    t_final = d_orig / (best_ssq + s_orig)
    perm_t = pD / (best_ssq + pS)
    fdr = _permutation_fdr(t_final, perm_t, n_boot)

    # Parametric p-values for ranking / AUC
    pvalues = 2.0 * stats.t.sf(np.abs(t_final), df=df_post)

    result = pd.DataFrame({
        "ProteinName": prot_index,
        "log2FC": stats_df["log2FC"].values,
        "pvalue": pvalues,
        "adj_pvalue": fdr,
    })
    return result
