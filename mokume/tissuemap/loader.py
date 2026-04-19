"""
QPX parquet data loading for TissueMap.

Handles both LFQ and TMT datasets with automatic GIS detection
for internal-standard normalization.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

# Patterns that identify GIS / reference / pool channels in organism_part
_GIS_PATTERNS = re.compile(
    r"global\s*internal\s*standard|reference|pool", re.IGNORECASE
)


def _parse_samples_col(samples_val) -> list[dict]:
    """Parse the ``run.samples`` column (list/ndarray of dicts) into a list."""
    if samples_val is None:
        return []
    if isinstance(samples_val, (list, np.ndarray)):
        return list(samples_val)
    return []


def build_run_to_sample_map(ds_dir: Path, ds_id: str) -> pd.DataFrame:
    """Return DataFrame: run_file_name, tmt_label, sample_accession, fraction, batch."""
    run = pq.read_table(ds_dir / "qpx_output" / f"{ds_id}.run.parquet").to_pandas()
    rows: list[dict] = []
    for _, r in run.iterrows():
        for s in _parse_samples_col(r["samples"]):
            raw_label = s.get("label", "LFQ")
            label = raw_label.replace("AC=MS:1002038;NT=", "").strip()
            # Normalize common LFQ label variants
            if label.lower() in ("label free sample", "label-free", "lfq"):
                label = "LFQ"
            rows.append(
                {
                    "run_file_name": r["run_file_name"],
                    "fraction": r.get("fraction", 1),
                    "sample_accession": s["sample_accession"],
                    "tmt_label": label,
                    "batch": r["run_file_name"],
                }
            )
    return pd.DataFrame(rows)


def _detect_gis_accessions(
    ds_dir: Path,
    ds_id: str,
) -> set[str]:
    """Auto-detect GIS / reference sample accessions from sample.parquet.

    Reads ``organism_part`` and flags samples whose tissue annotation matches
    known GIS patterns (e.g. "global internal standard", "reference", "pool").
    Returns the set of sample_accession strings that are GIS channels.
    """
    smp = pq.read_table(ds_dir / "qpx_output" / f"{ds_id}.sample.parquet").to_pandas()
    if "organism_part" not in smp.columns:
        return set()

    gis_accessions = set(
        smp.loc[
            smp["organism_part"].fillna("").str.contains(_GIS_PATTERNS),
            "sample_accession",
        ]
    )
    if gis_accessions:
        logger.info(
            "[%s] GIS auto-detect: %d reference accessions found",
            ds_id,
            len(gis_accessions),
        )
    return gis_accessions


def _read_and_explode_features(
    ds_dir: Path,
    ds_id: str,
    feature_prefix: str | None,
) -> pd.DataFrame:
    """Read feature parquet, filter decoys, explode intensities."""
    prefix = feature_prefix or ds_id
    parquet_path = ds_dir / "qpx_output" / f"{prefix}.feature.parquet"
    logger.info("[%s] Reading %s", ds_id, parquet_path.name)

    feat = pq.read_table(
        parquet_path,
        columns=["anchor_protein", "run_file_name", "intensities", "is_decoy"],
    ).to_pandas()
    feat = feat[~feat["is_decoy"].fillna(False)]
    feat = feat.dropna(subset=["intensities"])
    feat = feat[feat["intensities"].map(lambda x: x is not None and len(x) > 0)]
    long = feat[["anchor_protein", "run_file_name", "intensities"]].explode(
        "intensities"
    )
    long = long.rename(columns={"anchor_protein": "protein"})

    # Extract label and intensity from dicts — single-pass list comprehension
    raw = long["intensities"].tolist()
    labels = [d.get("label", "LFQ") if isinstance(d, dict) else str(d) for d in raw]
    values = [d.get("intensity", 0.0) if isinstance(d, dict) else 0.0 for d in raw]
    long["tmt_label"] = labels
    long["intensity"] = values
    long = long.drop(columns=["intensities"])
    long["intensity"] = pd.to_numeric(long["intensity"], errors="coerce").fillna(0.0)
    long = long[long["intensity"] > 0]

    if long.empty:
        raise ValueError(f"{ds_id}: no valid intensity values")
    return long


def _normalize_tmt(
    long: pd.DataFrame,
    run_map: pd.DataFrame,
    ds_dir: Path,
    ds_id: str,
) -> pd.DataFrame:
    """Merge with run map and apply GIS normalization for TMT data."""
    long = long.merge(run_map, on=["run_file_name", "tmt_label"], how="left")
    long = long.dropna(subset=["sample_accession"])

    gis_accessions = _detect_gis_accessions(ds_dir, ds_id)
    if not gis_accessions:
        fallback_labels = {"TMT126", "TMT131"}
        gis_rows = run_map[run_map["tmt_label"].isin(fallback_labels)]
        gis_accessions = set(gis_rows["sample_accession"])
        logger.warning(
            "[%s] GIS auto-detect found nothing, falling back to "
            "labels %s (%d accessions)",
            ds_id,
            fallback_labels,
            len(gis_accessions),
        )

    is_gis = long["sample_accession"].isin(gis_accessions)
    gis_data = long[is_gis]
    bio_data = long[~is_gis].copy()

    if gis_data.empty:
        logger.warning("[%s] No GIS data found, skipping normalization", ds_id)
        return bio_data

    gis_mean = (
        gis_data.groupby(["protein", "run_file_name"])["intensity"]
        .mean()
        .reset_index()
        .rename(columns={"intensity": "gis_intensity"})
    )
    bio_data = bio_data.merge(gis_mean, on=["protein", "run_file_name"], how="inner")
    bio_data = bio_data[bio_data["gis_intensity"] > 0].copy()
    bio_data["intensity"] = bio_data["intensity"] / bio_data["gis_intensity"]
    bio_data = bio_data.drop(columns=["gis_intensity"])
    logger.info(
        "[%s] TMT GIS norm: %d GIS accessions, %d bio rows",
        ds_id,
        len(gis_accessions),
        len(bio_data),
    )
    return bio_data


def _pivot_and_annotate(
    merged: pd.DataFrame,
    ds_dir: Path,
    ds_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Annotate tissues, pivot to matrix, and build sample metadata."""
    prot_samp = (
        merged.groupby(["protein", "sample_accession"])["intensity"].sum().reset_index()
    )

    smp = pq.read_table(ds_dir / "qpx_output" / f"{ds_id}.sample.parquet").to_pandas()
    smp = smp[["sample_accession", "organism_part"]].rename(
        columns={"organism_part": "tissue"}
    )
    prot_samp = prot_samp.merge(smp, on="sample_accession", how="left")
    prot_samp = prot_samp[
        prot_samp["tissue"].notna()
        & ~prot_samp["tissue"]
        .str.lower()
        .isin(["not available", "n/a", "", "global internal standard"])
    ]

    mat = prot_samp.pivot_table(
        index="protein",
        columns="sample_accession",
        values="intensity",
        aggfunc="sum",
    )
    tissue_meta = (
        prot_samp[["sample_accession", "tissue"]]
        .drop_duplicates("sample_accession")
        .set_index("sample_accession")
    )
    return mat, tissue_meta


def _filter_low_sample_tissues(
    mat: pd.DataFrame,
    meta: pd.DataFrame,
    ds_id: str,
    min_tissue_samples: int,
    low_sample_warning_threshold: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop tissues with too few samples and warn about borderline ones."""
    tissue_counts = meta["tissue"].value_counts()
    if min_tissue_samples > 1:
        valid_tissues = tissue_counts[tissue_counts >= min_tissue_samples].index
        dropped = tissue_counts[tissue_counts < min_tissue_samples]
        if len(dropped) > 0:
            logger.warning(
                "[%s] Dropped %d tissues with < %d samples: %s",
                ds_id,
                len(dropped),
                min_tissue_samples,
                ", ".join(dropped.index[:5]),
            )
        valid_samples = meta[meta["tissue"].isin(valid_tissues)].index
        mat = mat[mat.columns.intersection(valid_samples)]
        meta = meta.loc[meta.index.intersection(valid_samples)]

    if low_sample_warning_threshold > 0:
        low_sample = tissue_counts[tissue_counts < low_sample_warning_threshold]
        if len(low_sample) > 0:
            logger.warning(
                "[%s] %d tissues have < %d samples (TS scores may be unreliable): %s",
                ds_id,
                len(low_sample),
                low_sample_warning_threshold,
                ", ".join(f"{t}({c})" for t, c in low_sample.items()),
            )
    return mat, meta


def _detect_quant_type(
    long: pd.DataFrame,
    is_tmt: bool | None,
) -> tuple[bool, str]:
    """Auto-detect quantification type from TMT labels."""
    labels_unique = long["tmt_label"].unique()
    detected_tmt = any("TMT" in label or "iTRAQ" in label for label in labels_unique)
    if is_tmt is None:
        is_tmt = detected_tmt
    return is_tmt, "TMT" if is_tmt else "LFQ"


def _merge_quant_data(
    long: pd.DataFrame,
    run_map: pd.DataFrame,
    ds_dir: Path,
    ds_id: str,
    is_tmt: bool,
) -> pd.DataFrame:
    """Merge feature data with sample mapping (TMT or LFQ path)."""
    if is_tmt:
        merged = _normalize_tmt(long, run_map, ds_dir, ds_id)
    else:
        run_map_lfq = run_map.drop_duplicates("run_file_name")[
            ["run_file_name", "sample_accession", "fraction"]
        ]
        merged = long.merge(run_map_lfq, on="run_file_name", how="left")
    return merged.dropna(subset=["sample_accession"])


def _attach_metadata(
    tissue_meta: pd.DataFrame,
    run_map: pd.DataFrame,
    quant_type: str,
    ds_id: str,
) -> None:
    """Attach batch, quant_type, and dataset columns to tissue metadata."""
    batch_map = run_map.drop_duplicates("sample_accession").set_index(
        "sample_accession"
    )["batch"]
    tissue_meta["batch"] = tissue_meta.index.map(batch_map).fillna("unknown")
    tissue_meta["quant_type"] = quant_type
    tissue_meta["dataset"] = ds_id


def load_dataset(
    ds_dir: Path,
    ds_id: str,
    *,
    is_tmt: bool | None = None,
    feature_prefix: str | None = None,
    min_tissue_samples: int = 1,
    low_sample_warning_threshold: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load a single QPX dataset and return (mat, meta).

    Parameters
    ----------
    ds_dir : Path
        Directory containing ``qpx_output/`` with parquet files.
    ds_id : str
        Dataset identifier (e.g. "PXD016999").
    is_tmt : bool | None
        Force TMT mode.  ``None`` = auto-detect from labels.
    feature_prefix : str | None
        Prefix for feature parquet (defaults to *ds_id*).
    min_tissue_samples : int
        Tissues with fewer samples are dropped.

    Returns
    -------
    mat : pd.DataFrame
        Proteins (rows) x samples (columns), raw intensities.
    meta : pd.DataFrame
        Index = sample_accession, columns = tissue, batch, quant_type.
    """
    long = _read_and_explode_features(ds_dir, ds_id, feature_prefix)
    run_map = build_run_to_sample_map(ds_dir, ds_id)

    is_tmt, quant_type = _detect_quant_type(long, is_tmt)
    merged = _merge_quant_data(long, run_map, ds_dir, ds_id, is_tmt)
    mat, tissue_meta = _pivot_and_annotate(merged, ds_dir, ds_id)
    _attach_metadata(tissue_meta, run_map, quant_type, ds_id)

    mat, tissue_meta = _filter_low_sample_tissues(
        mat,
        tissue_meta,
        ds_id,
        min_tissue_samples,
        low_sample_warning_threshold,
    )

    logger.info(
        "[%s] Loaded: %d proteins x %d samples, %d tissues (%s)",
        ds_id,
        mat.shape[0],
        mat.shape[1],
        tissue_meta["tissue"].nunique(),
        quant_type,
    )
    return mat, tissue_meta
