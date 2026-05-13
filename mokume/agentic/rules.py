"""Rule-based configuration generation (fallback when LLM unavailable)."""

from itertools import product
from pathlib import Path

import yaml

from mokume.agentic.profiler import DataProfile
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.rules")

_KNOWLEDGE_DIR = Path(__file__).parent / "knowledge"


def _load_heuristics() -> dict:
    """Load heuristics from YAML knowledge base."""
    path = _KNOWLEDGE_DIR / "heuristics.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _get_data_type_priors(heuristics: dict, data_type: str) -> dict:
    """Get recommended methods for the given data type."""
    priors = heuristics.get("data_type_priors", {})
    dt = priors.get(data_type, priors.get("LFQ", {}))
    return {
        "de_methods": list(dt.get("recommended_de", ["deqms", "limma", "proda"])),
        "fdr_methods": list(dt.get("recommended_fdr", ["bh", "ihw"])),
        "normalizations": list(dt.get("recommended_normalization", ["none"]))[:2],
        "imputations": list(dt.get("recommended_imputation", ["none"])),
    }


def _refine_de_methods(de: list[str], profile: DataProfile) -> list[str]:
    """Add DE methods based on data characteristics."""
    min_grp = min(profile.samples_per_condition.values(), default=3)
    if min_grp < 3 and "rots" not in de:
        de.append("rots")
    if profile.has_peptide_counts and "deqms" not in de:
        de.insert(0, "deqms")
    return de


def _apply_condition_rules(space: dict, profile: DataProfile) -> dict:
    """Refine search space based on data profile conditions."""
    space["de_methods"] = _refine_de_methods(space["de_methods"], profile)
    if profile.n_proteins > 500 and "ihw" not in space["fdr_methods"]:
        space["fdr_methods"].insert(0, "ihw")
    if profile.missing_rate < 0.15:
        others = [m for m in space["imputations"] if m != "none"]
        space["imputations"] = ["none"] + others
    if "none" not in space["normalizations"]:
        space["normalizations"].insert(0, "none")
    return space


def rule_propose(
    profile: DataProfile,
    max_configs: int = 18,
) -> list[CandidateConfig]:
    """Generate candidate configs using domain heuristics."""
    heuristics = _load_heuristics()
    space = _get_data_type_priors(heuristics, profile.data_type)
    space = _apply_condition_rules(space, profile)

    thresholds = [0.5, 1.0]
    combos = list(
        product(
            space["de_methods"],
            space["fdr_methods"],
            space["normalizations"],
            space["imputations"],
            thresholds,
        )
    )

    # Trim to budget
    if len(combos) > max_configs:
        combos = combos[:max_configs]

    configs = []
    for de, fdr, norm, imp, thr in combos:
        name = f"{de}_{fdr}_{norm}_{imp}_fc{thr}"
        configs.append(
            CandidateConfig(
                name=name,
                de_method=de,
                fdr_method=fdr,
                normalization=norm,
                imputation=imp,
                log2fc_threshold=thr,
                reasoning=f"Rule-based: {profile.data_type} prior",
            )
        )

    logger.info("Rule engine generated %d candidate configs", len(configs))
    return configs


def load_prompts() -> dict[str, str]:
    """Load LLM prompt templates from YAML."""
    path = _KNOWLEDGE_DIR / "prompts.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_benchmarks() -> dict:
    """Load benchmark reference data from YAML."""
    path = _KNOWLEDGE_DIR / "benchmarks.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)
