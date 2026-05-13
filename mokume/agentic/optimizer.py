"""Main orchestrator for the agentic optimization loop."""

from dataclasses import dataclass

import pandas as pd

from mokume.agentic.config import AgenticConfig
from mokume.agentic.evaluator import (
    compute_score_ground_truth,
    compute_score_unsupervised,
    evaluate,
)
from mokume.agentic.profiler import DataProfile, profile_data
from mokume.agentic.proposer import propose_configs
from mokume.agentic.reflector import reflect
from mokume.agentic.reporter import save_outputs
from mokume.agentic.runner import PreprocessCache, run_experiment
from mokume.agentic.state import (
    AgenticState,
    AuditEntry,
    CandidateConfig,
    EvaluationResult,
    RoundResult,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.optimizer")


@dataclass
class RoundContext:
    """Bundled data needed for running optimization rounds."""

    protein_df: pd.DataFrame
    sample_to_condition: dict
    contrast: tuple
    ground_truth: set | None
    peptide_counts: pd.Series | None
    config: AgenticConfig
    cache: PreprocessCache | None = None


def _trim_to_budget(
    configs: list[CandidateConfig],
    state: AgenticState,
    max_experiments: int,
) -> list[CandidateConfig]:
    """Trim configs to fit within remaining budget."""
    remaining = max_experiments - state.total_experiments
    if remaining <= 0:
        return []
    return configs[:remaining]


def _score_gt(results: list[EvaluationResult], config: AgenticConfig) -> None:
    """Assign ground-truth scores to results."""
    max_tp = max((r.tp or 0) for r in results) if results else 1
    total = max_tp + max((r.fp or 0) for r in results) if results else 1
    for r in results:
        r.score = compute_score_ground_truth(r, max_tp, total, config.weights)


def _score_unsup(results: list[EvaluationResult], config: AgenticConfig) -> None:
    """Assign unsupervised scores to results."""
    max_de = max(r.n_de_up + r.n_de_down for r in results) if results else 1
    for r in results:
        r.score = compute_score_unsupervised(r, max_de, config.weights)


def _score_results(
    results: list[EvaluationResult],
    has_ground_truth: bool,
    config: AgenticConfig,
) -> list[EvaluationResult]:
    """Compute composite scores for all results."""
    if has_ground_truth:
        _score_gt(results, config)
    else:
        _score_unsup(results, config)
    return results


def _run_round(
    round_num: int,
    configs: list[CandidateConfig],
    ctx: RoundContext,
) -> RoundResult:
    """Execute one round of experiments."""
    results = []
    for cfg in configs:
        try:
            de_df = run_experiment(
                cfg,
                ctx.protein_df,
                ctx.sample_to_condition,
                ctx.contrast,
                ctx.peptide_counts,
                cache=ctx.cache,
            )
            result = evaluate(
                cfg,
                de_df,
                ctx.protein_df,
                ctx.sample_to_condition,
                ctx.ground_truth,
            )
            results.append(result)
        except (ValueError, KeyError, RuntimeError, ArithmeticError) as exc:
            logger.error("Experiment %s failed: %s", cfg.name, exc)
            results.append(
                EvaluationResult(
                    config_name=cfg.name,
                    config=cfg.to_dict(),
                    score=-1.0,
                )
            )

    results = _score_results(
        results,
        ctx.ground_truth is not None,
        ctx.config,
    )
    results.sort(key=lambda r: r.score, reverse=True)

    best_name = results[0].config_name if results else ""
    return RoundResult(
        round_num=round_num,
        configs=configs,
        results=results,
        best_config_name=best_name,
    )


def _get_candidates(
    round_num: int,
    state: AgenticState,
    profile: DataProfile,
    config: AgenticConfig,
    seen: set[str],
) -> list[CandidateConfig] | None:
    """Get candidates for this round, or None to stop."""
    if round_num == 1:
        return propose_configs(profile, config, seen=seen)
    ref = state.rounds[-1].reflection
    if ref and ref.next_configs:
        # Drop reflector candidates already in seen set
        return [c for c in ref.next_configs if c.signature() not in seen]
    return None


def _update_best(state: AgenticState, rnd: RoundResult) -> None:
    """Update global best if this round improved."""
    if not rnd.results or rnd.results[0].score <= state.best_score:
        return
    state.best_score = rnd.results[0].score
    idx = next(i for i, c in enumerate(rnd.configs) if c.name == rnd.best_config_name)
    state.best_config = rnd.configs[idx]


def optimize_contrast(
    ctx: RoundContext,
    profile: DataProfile,
) -> AgenticState:
    """Run the full optimization loop for a single contrast."""
    state = AgenticState()
    seen_signatures: set[str] = set()
    logger.info("Optimizing contrast: %s vs %s", ctx.contrast[0], ctx.contrast[1])

    for round_num in range(1, ctx.config.max_rounds + 1):
        candidates = _get_candidates(
            round_num, state, profile, ctx.config, seen_signatures
        )
        if candidates is None:
            break
        candidates = _trim_to_budget(
            candidates,
            state,
            ctx.config.max_experiments,
        )
        if not candidates:
            logger.info("Round %d: no new candidates after dedup, stopping", round_num)
            break

        state.audit_trail.append(
            AuditEntry(
                step="propose",
                round_num=round_num,
                data=[c.name for c in candidates],
            )
        )

        rnd = _run_round(round_num, candidates, ctx)
        state.rounds.append(rnd)
        state.total_experiments += len(candidates)
        seen_signatures.update(c.signature() for c in candidates)
        _update_best(state, rnd)

        ref_result = reflect(profile, state.rounds, ctx.config)
        rnd.reflection = ref_result
        if ref_result.converged:
            state.converged = True
            break

    return state


def optimize(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    config: AgenticConfig,
    ground_truth: set[str] | None = None,
    peptide_counts: pd.Series | None = None,
) -> dict[str, AgenticState]:
    """Run optimization for all contrasts in config."""
    profile = profile_data(protein_df, sample_to_condition, peptide_counts)
    logger.info(
        "Data profile: %d proteins, %d samples, %.1f%% missing",
        profile.n_proteins,
        profile.n_samples,
        profile.missing_rate * 100,
    )

    all_states = {}
    for contrast in config.contrasts or []:
        ctx = RoundContext(
            protein_df=protein_df,
            sample_to_condition=sample_to_condition,
            contrast=contrast,
            ground_truth=ground_truth,
            peptide_counts=peptide_counts,
            config=config,
            cache=PreprocessCache(),
        )
        key = f"{contrast[0]}_vs_{contrast[1]}"
        state = optimize_contrast(ctx, profile)
        if ctx.cache is not None:
            stats = ctx.cache.stats()
            logger.info(
                "Preprocess cache for %s: %d hits, %d misses, %d unique combos",
                key,
                stats["hits"],
                stats["misses"],
                stats["unique_combos"],
            )
        save_outputs(profile, state, config)
        all_states[key] = state

    return all_states
