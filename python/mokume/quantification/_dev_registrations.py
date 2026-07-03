"""Registry bindings for quantification methods whose source files are
byte-locked (sidecar-drift CI gate) and therefore cannot carry an
``@PluginRegistry.register`` decorator inline.

At present this covers iBAQ. ``mokume.quantification.ibaq`` is locked and
exposes iBAQ as a *function* core (:func:`compute_pibaq` /
:func:`peptides_to_protein`, plus :class:`ConcentrationWeightByProteomicRuler`)
rather than a :class:`ProteinQuantificationMethod` subclass. To expose it
through :class:`~mokume.core.registry.PluginRegistry` we register a thin
adapter that *delegates* to that locked core -- it adds no quantification
math of its own. The delegation target (``compute_pibaq`` plus the FASTA
digest / family discovery it needs) is exactly the core used by
``FeaturesToProteins`` iBAQ stage and the standalone ``peptides2protein``
CLI, so all three entry points agree on iBAQ values when configured
identically.

Import this module for its side effect of registering the adapter::

    import mokume.quantification._dev_registrations  # noqa: F401
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
    IBAQ,
)
from mokume.core.logger import get_logger
from mokume.core.registry import PluginRegistry
from mokume.quantification.base import ProteinQuantificationMethod

# Delegation targets -- all imported FROM the locked iBAQ / families sources.
# Nothing here reimplements their algorithm; the adapter only marshals inputs.
from mokume.quantification.ibaq import compute_pibaq
from mokume.quantification.families import (
    discover_families,
    load_families_yaml,
    merge_overrides,
)

logger = get_logger("mokume.quantification._dev_registrations")


@PluginRegistry.register("quantification", "ibaq")
class IbaqQuantification(ProteinQuantificationMethod):
    """Registry adapter for paralog-aware iBAQ (piBAQ).

    This is a thin wrapper around the locked
    :func:`mokume.quantification.ibaq.compute_pibaq` core. iBAQ needs
    external context that a bare ``peptide_df`` cannot supply -- namely a
    FASTA digest (to count theoretically observable peptides per protein)
    and the protein families discovered from it. Those are built here with
    the same locked helpers the pipeline uses
    (:func:`mokume.io.fasta.digest_fasta_full`,
    :func:`mokume.quantification.families.discover_families`) and handed to
    ``compute_pibaq``.

    Parameters
    ----------
    fasta : str, optional
        Path to the FASTA used for identification. Required before
        :meth:`quantify` is called.
    enzyme : str
        Digestion enzyme. Defaults to ``"Trypsin"``.
    min_aa, max_aa : int
        Peptide length bounds for the in-silico digest. Default 7 / 30.
    min_shared : int
        Minimum shared peptides for two proteins to join a family.
    min_anchors, high_anchor_threshold : int
        Anchor-gating knobs forwarded verbatim to ``compute_pibaq``.
    families_yaml : str, optional
        Optional YAML family-override file.
    """

    def __init__(
        self,
        fasta: Optional[str] = None,
        enzyme: str = "Trypsin",
        min_aa: int = 7,
        max_aa: int = 30,
        min_shared: int = 2,
        min_anchors: int = 1,
        high_anchor_threshold: int = 3,
        families_yaml: Optional[str] = None,
    ):
        self.fasta = fasta
        self.enzyme = enzyme
        self.min_aa = min_aa
        self.max_aa = max_aa
        self.min_shared = min_shared
        self.min_anchors = min_anchors
        self.high_anchor_threshold = high_anchor_threshold
        self.families_yaml = families_yaml

    @property
    def name(self) -> str:
        return "iBAQ"

    @property
    def input_level(self) -> str:
        # iBAQ consumes assembled, normalized peptide intensities -- the
        # "standard" flow, same level as TopN / sum / median.
        return "peptides"

    def quantify(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str = PROTEIN_NAME,
        peptide_column: str = PEPTIDE_CANONICAL,
        intensity_column: str = NORM_INTENSITY,
        sample_column: str = SAMPLE_ID,
        run_column: Optional[str] = None,
    ) -> pd.DataFrame:
        """Compute iBAQ by delegating to the locked ``compute_pibaq`` core.

        The FASTA is a configuration input supplied at construction time
        (``IbaqQuantification(fasta=...)``); it is mandatory here. The
        method returns a long-format table (protein x sample) with the
        :data:`~mokume.core.constants.IBAQ` intensity renamed to the
        standard ``"Intensity"`` column for parity with the other
        quantification methods.
        """
        fasta_path = self.fasta
        if not fasta_path:
            raise ValueError(
                "iBAQ quantification requires a FASTA file. Provide it via "
                "IbaqQuantification(fasta=...)."
            )

        # Local import keeps io.fasta out of module import time (mirrors the
        # pipeline stage, which imports digest_fasta_full lazily).
        from mokume.io.fasta import digest_fasta_full

        accession_to_peptides, peptide_to_accessions, _ = digest_fasta_full(
            fasta=fasta_path,
            enzyme=self.enzyme,
            min_aa=self.min_aa,
            max_aa=self.max_aa,
            canonicalize_isoforms=True,
            compute_mw=False,
        )
        families = discover_families(
            accession_to_peptides,
            peptide_to_accessions,
            min_shared=self.min_shared,
        )
        if self.families_yaml:
            overrides = load_families_yaml(Path(self.families_yaml))
            families = merge_overrides(families, overrides)

        long_df = compute_pibaq(
            peptide_df,
            accession_to_peptides,
            peptide_to_accessions,
            families,
            mw_map=None,
            min_anchors=self.min_anchors,
            high_anchor_threshold=self.high_anchor_threshold,
        )

        if long_df.empty:
            return pd.DataFrame(columns=[protein_column, sample_column, "Intensity"])

        # Normalize output to the shared "Intensity" contract.
        if IBAQ in long_df.columns:
            long_df = long_df.rename(columns={IBAQ: "Intensity"})
        return long_df
