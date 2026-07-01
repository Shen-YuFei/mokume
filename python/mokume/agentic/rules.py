"""Rule-based configuration generation (fallback when LLM unavailable)."""

from itertools import product

from mokume.agentic.profiler import DataProfile
from mokume.agentic.skill_loader import (
    extract_prompts,
    load_skill,
    load_skill_data,
)
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.rules")


def _load_heuristics() -> dict:
    """Load heuristics from the proteomics-heuristics skill data."""
    return load_skill_data("proteomics-heuristics")


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


def _build_ensemble_configs(
    heuristics: dict,
    space: dict,
    profile: DataProfile,
) -> list[CandidateConfig]:
    """Materialise ensemble strategies from heuristics as CandidateConfigs."""
    strategies = heuristics.get("ensemble_strategies", {}) or {}
    if not strategies:
        return []

    baseline_fdr = space["fdr_methods"][0]
    baseline_norm = space["normalizations"][0]
    baseline_imp = space["imputations"][0]

    configs: list[CandidateConfig] = []
    for strat_name, strat in strategies.items():
        methods = strat.get("methods") or []
        min_k = int(strat.get("min_k", 2))
        if not methods:
            continue
        ensemble_str = ",".join(methods)
        configs.append(
            CandidateConfig(
                name=f"ensemble_{strat_name}_{baseline_fdr}_{baseline_norm}",
                de_method="ensemble",
                fdr_method=baseline_fdr,
                normalization=baseline_norm,
                imputation=baseline_imp,
                log2fc_threshold=0.5,
                ensemble=ensemble_str,
                ensemble_k=min_k,
                reasoning=(
                    f"Rule-based ensemble ({strat_name}): "
                    f"top-{min_k} consensus across {ensemble_str}"
                ),
            )
        )
    return configs


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

    # Reserve room for ensemble candidates so they are never dropped on trim
    ensemble_configs = _build_ensemble_configs(heuristics, space, profile)
    budget_for_singles = max(1, max_configs - len(ensemble_configs))
    if len(combos) > budget_for_singles:
        combos = combos[:budget_for_singles]

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

    configs.extend(ensemble_configs)
    logger.info("Rule engine generated %d candidate configs", len(configs))
    return configs


def load_prompts() -> dict[str, str]:
    """Load LLM prompt templates from skill markdown files.

    Returns a dict with keys: proposal_system, proposal_user,
    reflection_system, reflection_user, report_system, report_user.
    """
    result: dict[str, str] = {}
    for skill_name, prefix in [
        ("proteomics-proposal", "proposal"),
        ("proteomics-reflection", "reflection"),
        ("proteomics-reporting", "report"),
    ]:
        skill = load_skill(skill_name)
        prompts = extract_prompts(skill)
        result[f"{prefix}_system"] = prompts.get("system", "")
        result[f"{prefix}_user"] = prompts.get("user", "")
    return result


def load_benchmarks() -> dict:
    """Load benchmark reference data from the proteomics-benchmarks skill."""
    return load_skill_data("proteomics-benchmarks")
