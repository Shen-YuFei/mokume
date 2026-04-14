"""LimROTS: limma + ROTS bootstrap-optimized test statistic.

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
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats

from mokume.analysis._shared_stats import (
    _collect_protein_stats,
    _fit_f_prior,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.limrots")


_LIMROTS_OPTION_DEFAULTS = {
    "n_boot": 100,
    "n_threads": None,
    "seed": 42,
}


class _ResampledStats(NamedTuple):
    d: np.ndarray
    s: np.ndarray
    pd: np.ndarray
    ps: np.ndarray


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


def _empty_iteration_columns(n_prot: int) -> tuple[np.ndarray, np.ndarray]:
    """Return empty d and s columns for one resampling iteration."""
    return np.full(n_prot, np.nan), np.full(n_prot, np.nan)


def _moderated_columns(
    stats_df: pd.DataFrame,
    prot_to_idx: dict[str, int],
    n_prot: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Map moderated LimROTS statistics back to the full protein index."""
    d_col, s_col = _empty_iteration_columns(n_prot)
    if stats_df.empty:
        return d_col, s_col

    d_vals, s_vals, _, _ = _limma_moderated_ds(stats_df)
    for protein, d_val, s_val in zip(stats_df["ProteinName"].values, d_vals, s_vals):
        idx = prot_to_idx.get(protein)
        if idx is None:
            continue
        d_col[idx] = d_val
        s_col[idx] = s_val
    return d_col, s_col


def _safe_iteration_columns(
    label: str,
    b: int,
    left_samples: list[str],
    right_samples: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Collect moderated columns for one resampling iteration."""
    try:
        stats_df = _collect_protein_stats(
            _SHARED["log2_mat"],
            (left_samples, right_samples),
            (_SHARED["cond_a"], _SHARED["cond_b"]),
        )
        return _moderated_columns(
            stats_df,
            _SHARED["prot_to_idx"],
            _SHARED["n_prot"],
        )
    except Exception as exc:
        logger.warning("%s %d failed: %s", label, b, exc)
        return _empty_iteration_columns(_SHARED["n_prot"])


def _worker_bootstrap(b: int):
    """Thin wrapper that reads shared data for ProcessPoolExecutor."""
    return _bootstrap_one(b)


def _bootstrap_one(b: int):
    """Run one bootstrap and permutation iteration from shared worker state."""
    samples_a = _SHARED["samples_a"]
    samples_b = _SHARED["samples_b"]

    # Bootstrap
    sa_b = [samples_a[i] for i in _SHARED["boot_sa"][b]]
    sb_b = [samples_b[i] for i in _SHARED["boot_sb"][b]]
    d_col, s_col = _safe_iteration_columns("Bootstrap", b, sa_b, sb_b)

    # Permutation
    all_samples = _SHARED["all_samples"]
    split = _SHARED["n_a"]
    permuted = _SHARED["perm_orders"][b]
    pa = [all_samples[i] for i in permuted[:split]]
    pb = [all_samples[i] for i in permuted[split:]]
    pd_col, ps_col = _safe_iteration_columns("Permutation", b, pa, pb)

    return b, d_col, s_col, pd_col, ps_col


def _generate_sampling_plan(
    rng: np.random.Generator,
    samples_a: list[str],
    samples_b: list[str],
    n_boot: int,
) -> dict:
    """Pre-generate bootstrap and permutation indices."""
    all_samples = samples_a + samples_b
    return {
        "all_samples": all_samples,
        "boot_sa": [
            rng.choice(len(samples_a), len(samples_a), replace=True)
            for _ in range(n_boot)
        ],
        "boot_sb": [
            rng.choice(len(samples_b), len(samples_b), replace=True)
            for _ in range(n_boot)
        ],
        "perm_orders": [rng.permutation(len(all_samples)) for _ in range(n_boot)],
        "n_a": len(samples_a),
    }


def _run_parallel_iterations(
    n_prot: int,
    n_boot: int,
    n_threads: int,
    shared_dict: dict,
) -> _ResampledStats:
    """Run LimROTS bootstrap and permutation iterations in parallel."""
    d_mat = np.full((n_prot, n_boot), np.nan)
    s_mat = np.full((n_prot, n_boot), np.nan)
    pd_mat = np.full((n_prot, n_boot), np.nan)
    ps_mat = np.full((n_prot, n_boot), np.nan)

    with ProcessPoolExecutor(
        max_workers=n_threads,
        initializer=_init_worker,
        initargs=(shared_dict,),
    ) as pool:
        for b, d_col, s_col, pd_col, ps_col in pool.map(_worker_bootstrap, range(n_boot)):
            d_mat[:, b] = d_col
            s_mat[:, b] = s_col
            pd_mat[:, b] = pd_col
            ps_mat[:, b] = ps_col

    return _ResampledStats(
        np.nan_to_num(d_mat),
        np.nan_to_num(s_mat, nan=1.0),
        np.nan_to_num(pd_mat),
        np.nan_to_num(ps_mat, nan=1.0),
    )


def _build_limrots_result(
    stats_df: pd.DataFrame,
    prot_index: np.ndarray,
    pvalues: np.ndarray,
    fdr: np.ndarray,
) -> pd.DataFrame:
    """Assemble the final LimROTS result table."""
    return pd.DataFrame(
        {
            "ProteinName": prot_index,
            "log2FC": stats_df["log2FC"].values,
            "pvalue": pvalues,
            "adj_pvalue": fdr,
        }
    )


def _build_shared_state(
    log2_matrix: pd.DataFrame,
    sample_groups: tuple[list[str], list[str]],
    contrast: tuple[str, str],
    prot_index: np.ndarray,
    sampling_plan: dict,
) -> dict:
    """Build the shared worker state for LimROTS resampling."""
    samples_a, samples_b = sample_groups
    cond_a, cond_b = contrast
    return {
        "log2_mat": log2_matrix,
        "samples_a": samples_a,
        "samples_b": samples_b,
        "cond_a": cond_a,
        "cond_b": cond_b,
        "n_prot": len(prot_index),
        "prot_to_idx": {protein: idx for idx, protein in enumerate(prot_index)},
        **sampling_plan,
    }


def _top_overlap(values_a: np.ndarray, values_b: np.ndarray, n_top: int) -> float:
    """Compute overlap fraction between top-ranked proteins."""
    top_a = set(np.argsort(-values_a)[:n_top])
    top_b = set(np.argsort(-values_b)[:n_top])
    return len(top_a & top_b) / n_top


def _separation_zscore(resampled: _ResampledStats, ssq: float, n_top: int) -> float:
    """Measure bootstrap-versus-permutation separation for one ssq value."""
    overlaps = []
    perm_overlaps = []
    for idx in range(0, resampled.d.shape[1] - 1, 2):
        overlaps.append(
            _top_overlap(
                resampled.d[:, idx] / (ssq + resampled.s[:, idx]),
                resampled.d[:, idx + 1] / (ssq + resampled.s[:, idx + 1]),
                n_top,
            )
        )
        perm_overlaps.append(
            _top_overlap(
                resampled.pd[:, idx] / (ssq + resampled.ps[:, idx]),
                resampled.pd[:, idx + 1] / (ssq + resampled.ps[:, idx + 1]),
                n_top,
            )
        )

    if not overlaps:
        return -np.inf
    sd_ol = np.std(overlaps, ddof=1) if len(overlaps) > 1 else 1e-10
    return (np.mean(overlaps) - np.mean(perm_overlaps)) / max(sd_ol, 1e-10)


def _optimize_ssq(resampled: _ResampledStats, n_prot: int):
    """Find optimal fudge factor via z-score maximization."""
    s_pos = resampled.s[resampled.s > 0]
    if len(s_pos) == 0:
        return 0.01
    ssq_grid = np.quantile(s_pos, np.linspace(0, 0.95, 20))
    n_tops = [10, 25, 50, 100, min(200, n_prot // 4)]
    best_z = -np.inf
    best_ssq = ssq_grid[0]

    for ssq in ssq_grid:
        for n_top in n_tops:
            z = _separation_zscore(resampled, ssq, n_top)
            if z > best_z:
                best_z = z
                best_ssq = ssq

    return best_ssq


def _permutation_fdr(t_final, perm_t):
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


def _prepare_limrots_run(
    log2_matrix: pd.DataFrame,
    sample_groups: tuple[list[str], list[str]],
    contrast: tuple[str, str],
    options: dict,
):
    """Prepare shared inputs and original statistics for LimROTS."""
    stats_df = _collect_protein_stats(log2_matrix, sample_groups, contrast)
    if stats_df.empty:
        return None

    d_orig, s_orig, _, df_post = _limma_moderated_ds(stats_df)
    prot_index = stats_df["ProteinName"].values
    sampling_plan = _generate_sampling_plan(
        np.random.default_rng(options["seed"]),
        sample_groups[0],
        sample_groups[1],
        options["n_boot"],
    )
    return (
        stats_df,
        prot_index,
        d_orig,
        s_orig,
        df_post,
        options["n_threads"] or os.cpu_count() or 4,
        _build_shared_state(log2_matrix, sample_groups, contrast, prot_index, sampling_plan),
    )


def _final_limrots_outputs(
    d_orig: np.ndarray,
    s_orig: np.ndarray,
    df_post: np.ndarray,
    resampled: _ResampledStats,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the final parametric p-values and permutation FDR."""
    best_ssq = _optimize_ssq(resampled, len(d_orig))
    logger.info("LimROTS: optimal s0 = %.4f", best_ssq)
    t_final = d_orig / (best_ssq + s_orig)
    perm_t = resampled.pd / (best_ssq + resampled.ps)
    return 2.0 * stats.t.sf(np.abs(t_final), df=df_post), _permutation_fdr(t_final, perm_t)


def run_limrots(
    log2_matrix: pd.DataFrame,
    samples_a: list[str],
    samples_b: list[str],
    contrast: tuple[str, str],
    **options,
) -> pd.DataFrame:
    """Run LimROTS differential expression with bootstrap-optimized moderation."""
    unknown = set(options) - set(_LIMROTS_OPTION_DEFAULTS)
    if unknown:
        unknown_args = ", ".join(sorted(unknown))
        raise TypeError(f"Unexpected LimROTS options: {unknown_args}")

    settings = {**_LIMROTS_OPTION_DEFAULTS, **options}
    prepared = _prepare_limrots_run(log2_matrix, (samples_a, samples_b), contrast, settings)
    if prepared is None:
        return pd.DataFrame()

    stats_df, prot_index, d_orig, s_orig, df_post, n_threads, shared_dict = prepared

    logger.info("LimROTS: %d bootstrap iterations (%d workers)", settings["n_boot"], n_threads)
    resampled = _run_parallel_iterations(shared_dict["n_prot"], settings["n_boot"], n_threads, shared_dict)
    pvalues, fdr = _final_limrots_outputs(d_orig, s_orig, df_post, resampled)
    return _build_limrots_result(stats_df, prot_index, pvalues, fdr)
