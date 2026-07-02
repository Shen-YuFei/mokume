"""Report generation for agentic analysis (LLM + template fallback)."""

import json
from pathlib import Path

import yaml

from mokume.agentic.config import AgenticConfig
from mokume.agentic.llm_client import LLMUnavailableError, call_text
from mokume.agentic.profiler import DataProfile
from mokume.agentic.rules import load_prompts
from mokume.agentic.state import AgenticState, EvaluationResult
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.reporter")


def contrast_slug(contrast: tuple[str, str]) -> str:
    """Filesystem-safe slug for a contrast pair (e.g. 'A_vs_B')."""
    return f"{contrast[0]}_vs_{contrast[1]}"


def checkpoint_path(output_dir: str | Path, contrast: tuple[str, str]) -> Path:
    """Return the state snapshot path for a contrast."""
    return Path(output_dir) / contrast_slug(contrast) / "state.json"


def save_state_snapshot(
    output_dir: str | Path,
    contrast: tuple[str, str],
    state: AgenticState,
) -> Path:
    """Persist the per-contrast state for --resume."""
    path = checkpoint_path(output_dir, contrast)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(state.to_dict(), indent=2, default=str),
        encoding="utf-8",
    )
    logger.debug("Saved state snapshot to %s", path)
    return path


def load_state_snapshot(
    output_dir: str | Path,
    contrast: tuple[str, str],
) -> AgenticState | None:
    """Load a per-contrast state snapshot if one exists."""
    path = checkpoint_path(output_dir, contrast)
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Failed to read checkpoint %s (%s); starting fresh", path, exc)
        return None
    state = AgenticState.from_dict(data)
    logger.info(
        "Resumed state from %s: %d rounds, %d experiments",
        path,
        len(state.rounds),
        state.total_experiments,
    )
    return state


def _all_results(state: AgenticState) -> list[EvaluationResult]:
    """Flatten all results across rounds."""
    results = []
    for rnd in state.rounds:
        results.extend(rnd.results)
    return sorted(results, key=lambda r: r.score, reverse=True)


def _results_to_tsv(results: list[EvaluationResult]) -> str:
    """Format results as TSV string."""
    header = "config\tscore\ttp\tfp\tauc\tn_up\tn_down\tmedian_cv\tmissing_rate"
    lines = [header]
    for r in results:
        lines.append(
            f"{r.config_name}\t{r.score:.4f}\t{r.tp}\t{r.fp}\t"
            f"{r.auc}\t{r.n_de_up}\t{r.n_de_down}\t"
            f"{r.median_cv}\t{r.missing_rate:.4f}"
        )
    return "\n".join(lines)


def _template_report(
    profile: DataProfile,
    state: AgenticState,
    results: list[EvaluationResult],
) -> str:
    """Generate a template-based Markdown report."""
    best = results[0] if results else None
    lines = [
        "# Agentic Analysis Report",
        "",
        "## Data Profile",
        f"- Proteins: {profile.n_proteins}",
        f"- Samples: {profile.n_samples}",
        f"- Conditions: {profile.n_conditions}",
        f"- Missing rate: {profile.missing_rate:.1%}",
        f"- Data type: {profile.data_type}",
        "",
        "## Optimization Summary",
        f"- Rounds: {len(state.rounds)}",
        f"- Total experiments: {state.total_experiments}",
        f"- Converged: {state.converged}",
        "",
    ]

    if best:
        lines.extend(
            [
                "## Best Configuration",
                f"- **Name**: {best.config_name}",
                f"- **Score**: {best.score:.4f}",
                f"- **TP**: {best.tp}, **FP**: {best.fp}, **AUC**: {best.auc}",
                f"- **DE UP**: {best.n_de_up}, **DE DOWN**: {best.n_de_down}",
                "",
            ]
        )

    lines.append("## All Results")
    lines.append("")
    lines.append("| Config | Score | TP | FP | AUC | UP | DOWN |")
    lines.append("|--------|-------|----|----|-----|----|----- |")
    for r in results:
        lines.append(
            f"| {r.config_name} | {r.score:.3f} | "
            f"{r.tp} | {r.fp} | {r.auc} | "
            f"{r.n_de_up} | {r.n_de_down} |"
        )

    return "\n".join(lines)


def _llm_report(
    profile: DataProfile,
    state: AgenticState,
    results: list[EvaluationResult],
    config: AgenticConfig,
) -> str:
    """Generate report using LLM."""
    prompts = load_prompts()
    best = results[0] if results else None

    system_msg = prompts["report_system"]
    user_msg = prompts["report_user"].format(
        data_profile_json=json.dumps(profile.to_dict(), indent=2),
        all_results_table=_results_to_tsv(results),
        best_config=json.dumps(best.config if best else {}, indent=2),
    )
    return call_text(system_msg, user_msg, config)


def save_outputs(
    profile: DataProfile,
    state: AgenticState,
    contrast: tuple[str, str],
    config: AgenticConfig,
) -> Path:
    """Save all optimization outputs to the per-contrast output dir."""
    out = Path(config.output_dir) / contrast_slug(contrast)
    out.mkdir(parents=True, exist_ok=True)

    results = _all_results(state)

    # 1. summary.md
    try:
        if config.use_llm:
            report_text = _llm_report(profile, state, results, config)
        else:
            report_text = _template_report(profile, state, results)
    except (LLMUnavailableError, ValueError, ConnectionError):
        report_text = _template_report(profile, state, results)

    (out / "summary.md").write_text(report_text, encoding="utf-8")

    # 2. best_config.yaml
    if state.best_config:
        (out / "best_config.yaml").write_text(
            yaml.dump(state.best_config.to_dict(), default_flow_style=False),
            encoding="utf-8",
        )

    # 3. all_results.tsv
    (out / "all_results.tsv").write_text(
        _results_to_tsv(results),
        encoding="utf-8",
    )

    # 4. data_profile.json
    (out / "data_profile.json").write_text(
        json.dumps(profile.to_dict(), indent=2),
        encoding="utf-8",
    )

    # 5. audit_trail.json
    audit = [
        {"step": e.step, "round": e.round_num, "data": str(e.data)}
        for e in state.audit_trail
    ]
    (out / "audit_trail.json").write_text(
        json.dumps(audit, indent=2),
        encoding="utf-8",
    )

    logger.info("Outputs saved to %s", out)
    return out
