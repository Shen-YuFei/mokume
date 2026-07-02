"""Main orchestrator for the agentic optimization loop."""

from dataclasses import dataclass

import pandas as pd

from joblib import Parallel, delayed

from mokume.agentic.config import AgenticConfig
from mokume.agentic.evaluator import (
    compute_score_ground_truth,
    compute_score_unsupervised,
    evaluate,
)
from mokume.agentic.profiler import DataProfile, profile_data
from mokume.agentic.proposer import propose_configs
from mokume.agentic.reflector import reflect
from mokume.agentic.reporter import (
    contrast_slug,
    load_state_snapshot,
    save_outputs,
    save_state_snapshot,
)
from mokume.agentic.runner import (
    PreprocessCache,
    _apply_imputation,
    _apply_normalization,
    run_experiment,
)
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


def _processed_matrix(cfg: CandidateConfig, ctx: RoundContext) -> pd.DataFrame:
    """Return the per-config normalized/imputed matrix used by run_experiment.

    CV / missingness must be measured on this processed matrix, not the raw
    shared ``protein_df``; otherwise they are identical constants across every
    candidate. When a cache is present the (norm, imp) prefix was just computed
    by run_experiment, so this is a cheap hit returning the same matrix.
    """
    if ctx.cache is not None:
        return ctx.cache.get_or_compute(
            cfg.normalization, cfg.imputation, ctx.protein_df
        )
    normed_df = _apply_normalization(ctx.protein_df, cfg.normalization)
    return _apply_imputation(normed_df, cfg.imputation)


def _evaluate_one(cfg: CandidateConfig, ctx: RoundContext) -> EvaluationResult:
    """Run + evaluate a single config, returning a sentinel result on failure."""
    try:
        de_df = run_experiment(
            cfg,
            ctx.protein_df,
            ctx.sample_to_condition,
            ctx.contrast,
            ctx.peptide_counts,
            cache=ctx.cache,
        )
        return evaluate(
            cfg,
            de_df,
            _processed_matrix(cfg, ctx),
            ctx.sample_to_condition,
            ctx.ground_truth,
        )
    except (ValueError, KeyError, RuntimeError, ArithmeticError) as exc:
        logger.error("Experiment %s failed: %s", cfg.name, exc)
        return EvaluationResult(
            config_name=cfg.name,
            config=cfg.to_dict(),
            score=-1.0,
        )


def _run_round(
    round_num: int,
    configs: list[CandidateConfig],
    ctx: RoundContext,
) -> RoundResult:
    """Execute one round of experiments (parallel when ctx.config.n_jobs != 1)."""
    n_jobs = ctx.config.n_jobs
    if n_jobs == 1 or len(configs) <= 1:
        results = [_evaluate_one(cfg, ctx) for cfg in configs]
    else:
        # Threading backend keeps the shared PreprocessCache useful and
        # avoids pickling the (potentially large) protein matrix per worker.
        results = list(
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(_evaluate_one)(cfg, ctx) for cfg in configs
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


def _load_resume_state(ctx: RoundContext) -> tuple[AgenticState, set[str], int]:
    """Return (state, seen_signatures, start_round) honouring ctx.config.resume."""
    if not ctx.config.resume:
        return AgenticState(), set(), 1
    loaded = load_state_snapshot(ctx.config.output_dir, ctx.contrast)
    if loaded is None:
        return AgenticState(), set(), 1
    if loaded.converged:
        logger.info(
            "Resumed state for %s is already converged; skipping contrast",
            ctx.contrast,
        )
    seen = {cfg.signature() for rnd in loaded.rounds for cfg in rnd.configs}
    start_round = len(loaded.rounds) + 1
    return loaded, seen, start_round


def optimize_contrast(
    ctx: RoundContext,
    profile: DataProfile,
) -> AgenticState:
    """Run the full optimization loop for a single contrast."""
    state, seen_signatures, start_round = _load_resume_state(ctx)
    logger.info(
        "Optimizing contrast: %s vs %s (starting at round %d)",
        ctx.contrast[0],
        ctx.contrast[1],
        start_round,
    )
    if state.converged:
        return state

    for round_num in range(start_round, ctx.config.max_rounds + 1):
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

        # Checkpoint after each completed round so --resume can pick up here
        save_state_snapshot(ctx.config.output_dir, ctx.contrast, state)

        if state.converged:
            break

    return state


@dataclass
class ContrastJob:
    """Inputs needed to optimize one contrast (shared across serial / parallel)."""

    protein_df: pd.DataFrame
    sample_to_condition: dict[str, str]
    config: AgenticConfig
    ground_truth: set[str] | None
    peptide_counts: pd.Series | None
    profile: DataProfile
    cache: PreprocessCache


def _run_one_contrast(
    contrast: tuple[str, str],
    job: ContrastJob,
) -> tuple[str, AgenticState]:
    """Run optimization for a single contrast against the shared cache."""
    ctx = RoundContext(
        protein_df=job.protein_df,
        sample_to_condition=job.sample_to_condition,
        contrast=contrast,
        ground_truth=job.ground_truth,
        peptide_counts=job.peptide_counts,
        config=job.config,
        cache=job.cache,
    )
    state = optimize_contrast(ctx, job.profile)
    save_outputs(job.profile, state, contrast, job.config)
    return contrast_slug(contrast), state


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

    contrasts = list(config.contrasts or [])
    # (norm, imp) is contrast-independent, so we share one cache instance
    # across every contrast in this run.
    job = ContrastJob(
        protein_df=protein_df,
        sample_to_condition=sample_to_condition,
        config=config,
        ground_truth=ground_truth,
        peptide_counts=peptide_counts,
        profile=profile,
        cache=PreprocessCache(),
    )
    max_concurrent = config.max_concurrent_contrasts

    if max_concurrent == 1 or len(contrasts) <= 1:
        pairs = [_run_one_contrast(contrast, job) for contrast in contrasts]
    else:
        pairs = list(
            Parallel(n_jobs=max_concurrent, backend="threading")(
                delayed(_run_one_contrast)(contrast, job) for contrast in contrasts
            )
        )

    all_states = dict(pairs)

    stats = job.cache.stats()
    logger.info(
        "Shared preprocess cache: %d hits, %d misses, %d unique combos",
        stats["hits"],
        stats["misses"],
        stats["unique_combos"],
    )

    return all_states
