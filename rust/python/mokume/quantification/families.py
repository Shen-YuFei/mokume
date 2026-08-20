"""Protein-family discovery for piBAQ (paralog-aware iBAQ) quantification.

When two or more proteins share a substantial fraction of their tryptic
peptides, shotgun proteomics cannot reliably attribute observed intensity
to a single member. Examples include the actin family (Actc1/Actb/Actg1
share >99% sequence identity), histone variants, tropomyosins, and the
many UniProt splice isoforms that map onto the same gene.

This module groups such proteins into *families* so that downstream
quantification (see :mod:`mokume.quantification.pibaq`) can either keep
per-protein resolution where evidence allows or fall back to a single
family-level estimate where it does not.

Two layers are applied in sequence, mirroring the gpGrouper protein
grouping conventions (Saltzman et al. 2018 *MCP* 17:2270):

1. **UniProt isoform collapse** -- already handled upstream by
   :func:`mokume.io.fasta.digest_fasta_full` when
   ``canonicalize_isoforms=True``.
2. **Shared-peptide connected components** -- :func:`discover_families`
   below builds the bipartite peptide-protein graph and groups proteins
   that share at least ``min_shared`` distinct digested peptides.

A third optional layer accepts an explicit YAML override
(:func:`load_families_yaml` + :func:`merge_overrides`) so that audit-
critical analyses can pin family definitions independent of the
data-driven CC.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@dataclass(frozen=True)
class Family:
    """A group of proteins that must be quantified jointly.

    Attributes
    ----------
    family_id : str
        Stable identifier for the family. By default the canonical
        accession with the most distinct digested peptides; YAML
        overrides may set a human-readable name (``ACT``, ``HIST_H2A``).
    members : tuple[str, ...]
        Canonical accessions of every member, sorted lexicographically
        for reproducibility.
    """

    family_id: str
    members: Tuple[str, ...]

    @property
    def size(self) -> int:
        return len(self.members)

    def __post_init__(self) -> None:
        if not self.members:
            raise ValueError(f"Family {self.family_id!r} has no members")
        # ``dataclass(frozen=True)`` blocks normal assignment; ``object.__setattr__``
        # is the canonical escape used to coerce inputs in post-init hooks. We
        # sort once here so callers can pass ``members`` in any order and still
        # get the documented lexicographic ordering used for stable diffs and
        # regression assertions.
        sorted_members = tuple(sorted(self.members))
        if sorted_members != self.members:
            object.__setattr__(self, "members", sorted_members)


def _build_shared_peptide_edges(
    accession_to_peptides: Mapping[str, Set[str]],
    peptide_to_accessions: Mapping[str, Set[str]],
    min_shared: int,
) -> Dict[str, Set[str]]:
    """Return the adjacency list of the shared-peptide graph.

    Only edges whose two endpoints share at least ``min_shared`` distinct
    peptides are kept. Iterating peptides (with the inverted index already
    available) is asymptotically cheaper than enumerating all O(N^2)
    protein pairs because the per-peptide accession set is typically
    small (1-5 elements).
    """
    pair_counts: "defaultdict[Tuple[str, str], int]" = defaultdict(int)
    for accessions in peptide_to_accessions.values():
        if len(accessions) < 2:
            continue
        sorted_accs = sorted(accessions)
        for i, acc_a in enumerate(sorted_accs):
            for acc_b in sorted_accs[i + 1 :]:
                pair_counts[(acc_a, acc_b)] += 1

    adjacency: Dict[str, Set[str]] = {acc: set() for acc in accession_to_peptides}
    for (acc_a, acc_b), count in pair_counts.items():
        if count >= min_shared:
            adjacency.setdefault(acc_a, set()).add(acc_b)
            adjacency.setdefault(acc_b, set()).add(acc_a)
    return adjacency


def _connected_components(adjacency: Mapping[str, Set[str]]) -> List[Set[str]]:
    """Iterative connected-component traversal (BFS) for the adjacency map.

    Recursion is intentionally avoided to keep stack usage bounded on
    proteomes with very large connected components (e.g. histone or
    keratin clusters).
    """
    seen: Set[str] = set()
    components: List[Set[str]] = []
    for start in adjacency:
        if start in seen:
            continue
        stack: List[str] = [start]
        component: Set[str] = set()
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            component.add(node)
            stack.extend(adjacency[node] - seen)
        components.append(component)
    return components


def _choose_family_id(
    members: Iterable[str],
    accession_to_peptides: Mapping[str, Set[str]],
) -> str:
    """Pick the canonical accession with the most peptides as the family ID.

    Ties are broken by lexicographic order so that the choice is stable
    run-to-run.
    """
    return max(
        members,
        key=lambda acc: (len(accession_to_peptides.get(acc, ())), acc),
    )


def discover_families(
    accession_to_peptides: Mapping[str, Set[str]],
    peptide_to_accessions: Mapping[str, Set[str]],
    *,
    min_shared: int = 2,
) -> List[Family]:
    """Group proteins into families using shared-peptide connected components.

    Parameters
    ----------
    accession_to_peptides : Mapping[str, Set[str]]
        Per-accession digested peptide sets, typically the first return
        value of :func:`mokume.io.fasta.digest_fasta_full`.
    peptide_to_accessions : Mapping[str, Set[str]]
        Inverted peptide-to-accessions index from the same source.
    min_shared : int, optional
        Minimum number of distinct peptides two proteins must share to be
        placed in the same family. Defaults to ``2``: a single shared
        peptide is rarely sufficient evidence of homology and may instead
        reflect chance digestion collisions or sequence-search ambiguity.

    Returns
    -------
    list[Family]
        One :class:`Family` per connected component. Every accession in
        ``accession_to_peptides`` appears in exactly one family, including
        isolated proteins which become singleton families.
    """
    if min_shared < 1:
        raise ValueError(f"min_shared must be >= 1 (got {min_shared})")
    adjacency = _build_shared_peptide_edges(
        accession_to_peptides, peptide_to_accessions, min_shared
    )
    components = _connected_components(adjacency)

    families: List[Family] = []
    for component in components:
        family_id = _choose_family_id(component, accession_to_peptides)
        members = tuple(sorted(component))
        families.append(Family(family_id=family_id, members=members))

    multi = sum(1 for f in families if f.size > 1)
    logger.info(
        "Discovered %d families (%d singletons, %d multi-member)",
        len(families),
        len(families) - multi,
        multi,
    )
    return families


def load_families_yaml(path: Path) -> List[Family]:
    """Parse a user-supplied YAML file declaring explicit family overrides.

    Expected schema::

        families:
          - name: ACT
            members: [P60709, P63261, P68133]
          - name: HIST_H2A
            members: [P0C0S5, Q96QV6, P04908]

    Parameters
    ----------
    path : Path
        Path to the YAML file.

    Returns
    -------
    list[Family]
        One :class:`Family` per ``families[*]`` entry.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the YAML file is malformed (missing ``families`` key, missing
        ``name`` or ``members`` on an entry, empty ``members`` list, or
        duplicate ``name`` values).
    """
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - import error path
        raise ImportError(
            "Reading --families requires PyYAML. Install with: pip install pyyaml"
        ) from exc

    if not path.exists():
        raise FileNotFoundError(f"Families YAML not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle) or {}

    entries = document.get("families")
    if not isinstance(entries, list):
        raise ValueError(
            f"{path}: top-level key 'families' must be a list of "
            f"{{name, members}} entries"
        )

    families: List[Family] = []
    seen_names: Set[str] = set()
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise ValueError(f"{path}: entry #{index} must be a mapping")
        name = entry.get("name")
        members_raw = entry.get("members")
        if not isinstance(name, str) or not name:
            raise ValueError(f"{path}: entry #{index} is missing a 'name' string")
        if name in seen_names:
            raise ValueError(f"{path}: duplicate family name {name!r}")
        if not isinstance(members_raw, Sequence) or isinstance(members_raw, str):
            raise ValueError(
                f"{path}: entry {name!r} 'members' must be a list of accessions"
            )
        members = tuple(sorted(str(m) for m in members_raw if m))
        if not members:
            raise ValueError(f"{path}: entry {name!r} has no members")
        families.append(Family(family_id=name, members=members))
        seen_names.add(name)

    logger.info("Loaded %d family overrides from %s", len(families), path)
    return families


def merge_overrides(
    auto_families: Sequence[Family],
    overrides: Optional[Sequence[Family]],
) -> List[Family]:
    """Apply user-supplied overrides on top of automatically discovered families.

    Any accession that appears in an override family is removed from the
    automatically discovered families; the override family is then taken
    verbatim. Singleton auto-families whose member is consumed by an
    override are dropped, and remaining auto-families are kept in their
    original form.

    Parameters
    ----------
    auto_families : Sequence[Family]
        Output of :func:`discover_families`.
    overrides : Sequence[Family], optional
        Manual families from :func:`load_families_yaml`. When ``None`` or
        empty, ``auto_families`` is returned unchanged.

    Returns
    -------
    list[Family]
        Merged family list with overrides taking priority.
    """
    if not overrides:
        return list(auto_families)

    pinned: Set[str] = set()
    for override in overrides:
        pinned.update(override.members)

    merged: List[Family] = list(overrides)
    for family in auto_families:
        remaining = tuple(m for m in family.members if m not in pinned)
        if not remaining:
            continue
        if remaining == family.members:
            merged.append(family)
        else:
            # Some -- but not all -- members were claimed by an override.
            # Re-emit the survivors as a smaller auto-family so they are
            # not silently dropped.
            new_id = _smallest_canonical(remaining)
            merged.append(Family(family_id=new_id, members=remaining))
    logger.info(
        "Merged %d overrides with %d auto-families -> %d total",
        len(overrides),
        len(auto_families),
        len(merged),
    )
    return merged


def _smallest_canonical(members: Sequence[str]) -> str:
    """Pick a deterministic representative when no peptide-count info exists."""
    return min(members)


def families_by_member(families: Sequence[Family]) -> Dict[str, Family]:
    """Index families by member accession for O(1) lookup.

    Raises
    ------
    ValueError
        If the same accession appears in more than one family (should be
        impossible when :func:`merge_overrides` is the producer).
    """
    index: Dict[str, Family] = {}
    for family in families:
        for member in family.members:
            if member in index:
                raise ValueError(
                    f"Accession {member!r} appears in both "
                    f"{index[member].family_id!r} and {family.family_id!r}"
                )
            index[member] = family
    return index
