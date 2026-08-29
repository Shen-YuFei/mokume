"""
Censored-aware missing value imputation for proteomics data.

Implements multiple imputation strategies that account for the
non-random nature of missingness in mass spectrometry data:

- **MinProb**: Draws replacement values from the low tail of a
  fitted normal distribution (Perseus-style).
- **MinDet**: Replaces missing values with a fixed quantile of the
  per-sample observed values.

References
----------
- Lazar C, et al. Accounting for the Multiple Natures of Missing Values
  in Label-Free Quantitative Proteomics. J Proteome Res. 2016;15(4):1116-25.

"""

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.censored")


def impute_minprob(
    data: pd.DataFrame,
    quantile: float = 0.01,
    shift: float = 1.6,
    scale: float = 0.3,
    random_state: int = 42,
) -> pd.DataFrame:
    """Perform MinProb imputation using a low-tail normal draw per sample.

    Each missing entry is drawn from ``N(q_low - shift*sd, (scale*sd)**2)`` where
    ``q_low`` is the ``quantile``-th observed value. The center is anchored at the
    low observed quantile, which is more aggressive (lower) than Perseus's
    mean-anchored MinProb; adjust ``quantile`` / ``shift`` / ``scale`` to match a
    specific convention.

    ``random_state`` seeds the draws so repeated runs on the same matrix agree.
    Without it a benchmark cannot rank this method against a deterministic one:
    the run-to-run spread gets read as a difference between methods.
    """
    result = data.copy()
    rng = np.random.default_rng(random_state)
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) == 0:
            continue

        q_low = np.quantile(observed, quantile)
        sd = np.std(observed) * scale
        if sd < 1e-10:
            sd = 0.1
        mu = q_low - shift * sd

        imputed = rng.normal(mu, sd, size=missing.sum())
        result.loc[missing, col] = imputed
        n_imputed += missing.sum()

    logger.info("MinProb imputation: %d values imputed", n_imputed)
    return result


def impute_mindet(
    data: pd.DataFrame,
    quantile: float = 0.01,
) -> pd.DataFrame:
    """Apply MinDet imputation using a per-sample observed quantile."""
    result = data.copy()
    n_imputed = 0

    for col in data.columns:
        missing = data[col].isna()
        if not missing.any():
            continue
        observed = data[col].dropna().values
        if len(observed) == 0:
            continue
        fill_val = np.quantile(observed, quantile)
        result.loc[missing, col] = fill_val
        n_imputed += missing.sum()

    logger.info("MinDet imputation: %d values imputed (q=%.3f)", n_imputed, quantile)
    return result


def _knn_impute(data: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """KNN imputation wrapper for the dispatcher."""
    try:
        from sklearn.impute import KNNImputer
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for KNN imputation. "
            "Install it with: pip install mokume-py[imputation]"
        ) from exc
    imputer = KNNImputer(
        n_neighbors=kwargs.get("n_neighbors", 5),
        weights=kwargs.get("weights", "uniform"),
    )
    result = pd.DataFrame(
        imputer.fit_transform(data),
        index=data.index,
        columns=data.columns,
    )
    n_imputed = data.isna().sum().sum()
    logger.info("KNN imputation: %d values imputed", n_imputed)
    return result


def _dispatch_seqknn(data, **kwargs):
    from mokume.imputation.seqknn import impute_seqknn

    return impute_seqknn(data, **kwargs)


def _dispatch_missforest(data, **kwargs):
    from mokume.imputation.missforest import impute_missforest

    return impute_missforest(data, **kwargs)


def _dispatch_impseq(data, **kwargs):
    from mokume.imputation.impseq import impute_impseq

    return impute_impseq(data, **kwargs)


def _dispatch_impseqrob(data, **kwargs):
    from mokume.imputation.impseqrob import impute_impseqrob

    return impute_impseqrob(data, **kwargs)


def _dispatch_bpca(data, **kwargs):
    from mokume.imputation.bpca import impute_bpca

    return impute_bpca(data, **kwargs)


def _dispatch_qrilc(data, **kwargs):
    from mokume.imputation.qrilc import impute_qrilc

    return impute_qrilc(data, **kwargs)


def _dispatch_gms(data, **kwargs):
    from mokume.imputation.gms import impute_gms

    return impute_gms(data, **kwargs)


_DISPATCH = {
    "minprob": impute_minprob,
    "mindet": impute_mindet,
    "knn": _knn_impute,
    "seqknn": _dispatch_seqknn,
    "missforest": _dispatch_missforest,
    "impseq": _dispatch_impseq,
    "impseqrob": _dispatch_impseqrob,
    "bpca": _dispatch_bpca,
    "qrilc": _dispatch_qrilc,
    "gms": _dispatch_gms,
}


def impute_censored(
    data: pd.DataFrame,
    method: str = "minprob",
    **kwargs,
) -> pd.DataFrame:
    """Dispatch censored-aware imputation by method name."""
    method = method.lower()
    if method == "none":
        logger.info("No imputation applied (method='none')")
        return data.copy()
    handler = _DISPATCH.get(method)
    if handler is None:
        valid = ", ".join(sorted([*_DISPATCH, "none"]))
        raise ValueError(f"Unknown imputation method '{method}'. Choose from: {valid}")
    return handler(data, **kwargs)
