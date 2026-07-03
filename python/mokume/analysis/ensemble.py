"""
Multi-method ensemble for differential expression analysis.

Runs multiple DE methods on the same data and combines results
using a top-k consensus strategy: a protein is deemed significant
only when at least *k* methods agree on its direction.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

from mokume.analysis.differential_expression import DifferentialExpression
from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.ensemble")


def run_ensemble(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    contrast: tuple[str, str],
    methods: Sequence[str] = ("limrots", "deqms", "proda"),
    min_k: int = 2,
    fdr_method: str = "bh",
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 0.5,
    peptide_counts: pd.Series | None = None,
) -> pd.DataFrame:
    """Run multiple DE methods and combine with top-k consensus.

    Parameters
    ----------
    protein_df : pd.DataFrame
        Protein intensity matrix (first column = protein IDs).
    sample_to_condition : dict[str, str]
        Sample-to-condition mapping.
    contrast : tuple[str, str]
        Condition pair (A, B) for comparison.
    methods : sequence of str
        DE methods to run (default: limrots, deqms, proda).
    min_k : int
        Minimum number of methods that must agree on significance
        for a protein to be called DE (default 2).
    fdr_method : str
        FDR correction method passed to each sub-method.
    fdr_threshold : float
        FDR threshold for each sub-method.
    log2fc_threshold : float
        Log2 fold-change threshold.
    peptide_counts : pd.Series or None
        Peptide counts (passed to DEqMS if included).

    Returns
    -------
    pd.DataFrame
        Ensemble DE results with columns: ProteinName, log2FC (median),
        pvalue (Fisher combined), adj_pvalue, significance,
        n_methods_up, n_methods_down, methods_significant.
    """
    individual_results: dict[str, pd.DataFrame] = {}
    for method in methods:
        try:
            de = DifferentialExpression(
                method=method,
                fdr_method=fdr_method,
                fdr_threshold=fdr_threshold,
                log2fc_threshold=log2fc_threshold,
                peptide_counts=peptide_counts,
            )
            result = de.run(protein_df, sample_to_condition, contrast)
            if not result.empty:
                individual_results[method] = result
                logger.info("Ensemble member %s: %d proteins", method, len(result))
        except (ValueError, RuntimeError, KeyError, ArithmeticError) as exc:
            logger.warning("Ensemble member %s failed: %s", method, exc)

    if not individual_results:
        logger.warning("No ensemble members produced results")
        return pd.DataFrame()

    return combine_de_results(
        individual_results,
        min_k=min_k,
        log2fc_threshold=log2fc_threshold,
        fdr_threshold=fdr_threshold,
    )


def combine_de_results(
    results: dict[str, pd.DataFrame],
    min_k: int = 2,
    log2fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Combine DE results from multiple methods using top-k consensus.

    Parameters
    ----------
    results : dict[str, pd.DataFrame]
        Mapping of method name to its DE result DataFrame.
        Each must have columns: ProteinName, log2FC, pvalue, adj_pvalue,
        significance.
    min_k : int
        Minimum agreement count.
    log2fc_threshold : float
        Log2 fold-change threshold for re-classification.
    fdr_threshold : float
        FDR threshold for re-classification.

    Returns
    -------
    pd.DataFrame
        Combined results.
    """
    all_proteins: set[str] = set()
    for df in results.values():
        all_proteins.update(df["ProteinName"].tolist())

    rows = []
    for protein in sorted(all_proteins):
        fc_values = []
        p_values = []
        n_up = 0
        n_down = 0
        contributing: list[str] = []

        for method, df in results.items():
            match = df[df["ProteinName"] == protein]
            if match.empty:
                continue
            row = match.iloc[0]
            fc_values.append(row["log2FC"])
            pval = row["pvalue"]
            if np.isfinite(pval) and pval > 0:
                p_values.append(pval)
            sig = row.get("significance", "Unchanged")
            if sig == "UP":
                n_up += 1
                contributing.append(method)
            elif sig == "DOWN":
                n_down += 1
                contributing.append(method)

        median_fc = float(np.median(fc_values)) if fc_values else 0.0

        if len(p_values) >= 2:
            _, combined_p = combine_pvalues(p_values, method="fisher")
        elif len(p_values) == 1:
            combined_p = p_values[0]
        else:
            combined_p = np.nan

        rows.append(
            {
                "ProteinName": protein,
                "log2FC": median_fc,
                "pvalue": combined_p,
                "n_methods_up": n_up,
                "n_methods_down": n_down,
                "methods_significant": ",".join(contributing),
            }
        )

    ensemble = pd.DataFrame(rows)
    if ensemble.empty:
        return ensemble

    # BH correction on combined p-values
    valid = np.isfinite(ensemble["pvalue"].values)
    adj = np.full(len(ensemble), np.nan)
    if valid.any():
        from statsmodels.stats.multitest import multipletests

        adj[valid] = multipletests(
            ensemble.loc[valid, "pvalue"].values, method="fdr_bh"
        )[1]
    ensemble["adj_pvalue"] = adj

    # Consensus significance
    ensemble["significance"] = "Unchanged"
    ensemble.loc[
        (ensemble["n_methods_up"] >= min_k)
        & (ensemble["adj_pvalue"] < fdr_threshold)
        & (ensemble["log2FC"] > log2fc_threshold),
        "significance",
    ] = "UP"
    ensemble.loc[
        (ensemble["n_methods_down"] >= min_k)
        & (ensemble["adj_pvalue"] < fdr_threshold)
        & (ensemble["log2FC"] < -log2fc_threshold),
        "significance",
    ] = "DOWN"

    ensemble = ensemble.sort_values("adj_pvalue").reset_index(drop=True)

    n_up = (ensemble["significance"] == "UP").sum()
    n_down = (ensemble["significance"] == "DOWN").sum()
    logger.info(
        "Ensemble results: %d proteins, %d UP, %d DOWN (k>=%d)",
        len(ensemble),
        n_up,
        n_down,
        min_k,
    )
    return ensemble
