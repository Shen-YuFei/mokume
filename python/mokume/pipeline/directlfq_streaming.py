"""Streaming DirectLFQ protein estimation helpers."""

from __future__ import annotations

import logging
import multiprocessing as mp
from typing import Iterable, Optional

import numpy as np
import pandas as pd

LOGGER = logging.getLogger(__name__)


def iter_protein_groups(normed_df: pd.DataFrame) -> Iterable[tuple[int, pd.DataFrame]]:
    """Yield DirectLFQ protein groups without materialising all groups at once."""
    protein_names = normed_df.index.get_level_values(0).to_numpy()
    seen_proteins = set()
    last_protein = object()
    groups_are_contiguous = True
    for protein_name in protein_names:
        if protein_name != last_protein:
            if protein_name in seen_proteins:
                groups_are_contiguous = False
                break
            seen_proteins.add(protein_name)
            last_protein = protein_name

    if not groups_are_contiguous:
        LOGGER.info("Sorting DirectLFQ frame by protein/ion before streaming")
        normed_df = normed_df.sort_index()
        protein_names = normed_df.index.get_level_values(0).to_numpy()

    ion_names = normed_df.index.get_level_values(1).to_numpy()
    normed_array = normed_df.to_numpy()

    switch_indices = np.where(protein_names[:-1] != protein_names[1:])[0] + 1
    starts = np.insert(switch_indices, 0, 0)
    stops = np.append(switch_indices, len(protein_names))

    for group_idx, (start, stop) in enumerate(zip(starts, stops)):
        index = pd.MultiIndex.from_arrays(
            [protein_names[start:stop], ion_names[start:stop]],
            names=normed_df.index.names,
        )
        yield group_idx, pd.DataFrame(normed_array[start:stop], index=index)


def _estimate_protein_group(args):
    """Worker: estimate one protein and return only its name + profile.

    Returning just the protein name and profile (not the large
    ``shifted_peptides`` ion frame) keeps the data sent back from worker
    processes small, which is what lets the parallel path stay memory-bounded.
    """
    group_idx, protein_df, num_samples_quadratic, min_nonan = args
    from directlfq.protein_intensity_estimation import (
        calculate_peptide_and_protein_intensities,
    )

    protein_profile, shifted_peptides = calculate_peptide_and_protein_intensities(
        group_idx,
        protein_df,
        num_samples_quadratic,
        min_nonan,
    )
    if protein_profile is None:
        return None
    return shifted_peptides.index.get_level_values(0)[0], protein_profile


def estimate_protein_intensities_streamed(
    normed_df: pd.DataFrame,
    min_nonan: int,
    num_samples_quadratic: int,
    num_cores: Optional[int] = None,
) -> pd.DataFrame:
    """Estimate DirectLFQ protein intensities without retaining ion-level results.

    This mirrors the external directlfq package's protein estimation path, but
    keeps only the final protein profiles. It deliberately does not support
    normalized ion-table export, because that table is the large intermediate
    this path is meant to avoid.

    When ``num_cores > 1`` the per-protein estimation runs across a process pool.
    Unlike directlfq's own multiprocessing, proteins are streamed from
    :func:`iter_protein_groups` and dispatched in bounded batches (so at most a
    few cores' worth of protein frames are materialised at once) and workers
    return only the protein name and profile, keeping peak memory close to the
    sequential path. Output row order is preserved, so results are identical to
    the sequential path regardless of core count.
    """
    from directlfq import config as lfq_config

    protein_ids: list[str] = []
    profile_rows: list[np.ndarray] = []

    total_groups = normed_df.index.get_level_values(0).nunique()
    LOGGER.info("%d lfq-groups total", total_groups)

    if num_cores is not None and num_cores > 1:
        LOGGER.info("Estimating DirectLFQ proteins in parallel (%d cores)", num_cores)
        batch_size = num_cores * 4
        with mp.Pool(processes=num_cores) as pool:
            batch: list[tuple] = []
            for group_idx, protein_df in iter_protein_groups(normed_df):
                batch.append((group_idx, protein_df, num_samples_quadratic, min_nonan))
                if len(batch) >= batch_size:
                    _append_results(
                        pool.map(_estimate_protein_group, batch),
                        protein_ids,
                        profile_rows,
                    )
                    batch = []
            if batch:
                _append_results(
                    pool.map(_estimate_protein_group, batch), protein_ids, profile_rows
                )
    else:
        from directlfq.protein_intensity_estimation import (
            calculate_peptide_and_protein_intensities,
        )

        for group_idx, protein_df in iter_protein_groups(normed_df):
            protein_profile, shifted_peptides = (
                calculate_peptide_and_protein_intensities(
                    group_idx,
                    protein_df,
                    num_samples_quadratic,
                    min_nonan,
                )
            )
            if protein_profile is None:
                continue
            protein_ids.append(shifted_peptides.index.get_level_values(0)[0])
            profile_rows.append(protein_profile)

    if not profile_rows:
        return pd.DataFrame(columns=[lfq_config.PROTEIN_ID, *normed_df.columns])

    protein_index = pd.Index(protein_ids, name=lfq_config.PROTEIN_ID)
    protein_df = 2 ** pd.DataFrame(
        profile_rows,
        index=protein_index,
        columns=normed_df.columns,
    )
    protein_df = protein_df.replace(np.nan, 0).reset_index()
    return protein_df


def _append_results(results, protein_ids: list, profile_rows: list) -> None:
    """Append non-empty (name, profile) worker results, preserving order."""
    for result in results:
        if result is not None:
            protein_ids.append(result[0])
            profile_rows.append(result[1])
