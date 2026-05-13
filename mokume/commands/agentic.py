"""CLI commands for agentic analysis."""

from pathlib import Path

import click
import pandas as pd

from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import optimize
from mokume.normalization.irs import detect_condition_from_sdrf

_ENV_KEY_MAP = {
    "llm_api_key": "DEEPSEEK_API_KEY",
    "llm_base_url": "OPENAI_BASE_URL",
}


def _persist_to_dotenv(kwargs: dict) -> None:
    """Save LLM settings to .env if provided via CLI and not already stored."""
    dotenv_path = Path.cwd() / ".env"
    existing = {}
    if dotenv_path.exists():
        for line in dotenv_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, _, v = line.partition("=")
                existing[k.strip()] = v.strip()

    updates = {}
    for kwarg_key, env_key in _ENV_KEY_MAP.items():
        val = kwargs.get(kwarg_key)
        if val and existing.get(env_key) != val:
            updates[env_key] = val

    if not updates:
        return

    existing.update(updates)
    lines = [f"{k}={v}" for k, v in existing.items()]
    dotenv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    saved = ", ".join(updates.keys())
    click.echo(f"Saved to .env: {saved}")


def _parse_contrasts(raw: str) -> list[tuple[str, str]]:
    """Parse comma-separated contrast string into tuples."""
    parsed = []
    for c in raw.split(","):
        c = c.strip()
        if " vs " in c:
            parts = c.split(" vs ", 1)
        elif "-" in c:
            parts = c.split("-", 1)
        else:
            raise click.BadParameter(f"Invalid contrast: '{c}', use 'A vs B'")
        parsed.append((parts[0].strip(), parts[1].strip()))
    return parsed


def _load_ground_truth(path: str | None) -> set[str] | None:
    """Load ground truth protein set from file."""
    if not path:
        return None
    with open(path, encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip()}


@click.group("agentic")
def agentic_cmd():
    """AI-assisted analysis optimization."""


@agentic_cmd.command("optimize")
@click.option(
    "--protein-matrix",
    required=True,
    type=click.Path(exists=True),
    help="Protein intensity matrix (TSV/CSV, first col = protein ID).",
)
@click.option(
    "--sdrf",
    required=True,
    type=click.Path(exists=True),
    help="SDRF file for condition mapping.",
)
@click.option(
    "--contrasts",
    required=True,
    help="Contrasts to test, comma-separated (e.g. 'A vs B,A vs C').",
)
@click.option(
    "--ground-truth",
    default=None,
    type=click.Path(exists=True),
    help="File with ground truth protein IDs (one per line).",
)
@click.option(
    "--expected-fc",
    default=None,
    type=click.Path(exists=True),
    help="YAML file with expected log2FC per contrast.",
)
@click.option("--max-rounds", default=5, type=int, help="Max optimization rounds.")
@click.option("--max-experiments", default=30, type=int, help="Max total experiments.")
@click.option(
    "--llm-provider",
    default="deepseek",
    type=click.Choice(["deepseek", "custom"]),
    help="LLM provider preset (default: deepseek).",
)
@click.option(
    "--llm-base-url",
    default=None,
    envvar="OPENAI_BASE_URL",
    help="OpenAI-compatible API base URL. Required for custom provider.",
)
@click.option(
    "--llm-api-key",
    default=None,
    envvar=["DEEPSEEK_API_KEY", "OPENAI_API_KEY"],
    help="API key for the LLM service. Reads DEEPSEEK_API_KEY or OPENAI_API_KEY.",
)
@click.option(
    "--llm-model", default=None, help="LLM model name (auto-set from provider)."
)
@click.option(
    "--llm-thinking/--no-llm-thinking",
    "llm_thinking",
    default=True,
    show_default=True,
    help="Enable DeepSeek V4 Pro deep-thinking mode for proposer / reflector.",
)
@click.option(
    "--reasoning-effort",
    "llm_reasoning_effort",
    default="high",
    type=click.Choice(["high", "max"], case_sensitive=False),
    show_default=True,
    help="Reasoning depth when --llm-thinking is enabled (DeepSeek only).",
)
@click.option(
    "--no-llm",
    is_flag=True,
    default=False,
    help="Disable LLM, use rule-based mode only.",
)
@click.option(
    "--output-dir", "-o", default="./optimization", help="Output directory for results."
)
def optimize_cmd(**kwargs):
    """Run agentic optimization for differential expression analysis."""
    _persist_to_dotenv(kwargs)
    pm = kwargs["protein_matrix"]
    sep = "\t" if pm.endswith(".tsv") else ","
    protein_df = pd.read_csv(pm, sep=sep)
    s2c = detect_condition_from_sdrf(kwargs["sdrf"])
    parsed = _parse_contrasts(kwargs["contrasts"])
    gt_set = _load_ground_truth(kwargs["ground_truth"])

    config = AgenticConfig(
        sdrf=kwargs["sdrf"],
        protein_matrix=pm,
        ground_truth=kwargs["ground_truth"],
        expected_fc=kwargs["expected_fc"],
        max_rounds=kwargs["max_rounds"],
        max_experiments=kwargs["max_experiments"],
        use_llm=not kwargs["no_llm"],
        llm_provider=kwargs["llm_provider"],
        llm_base_url=kwargs["llm_base_url"],
        llm_api_key=kwargs["llm_api_key"],
        llm_model=kwargs["llm_model"],
        llm_thinking=kwargs["llm_thinking"],
        llm_reasoning_effort=kwargs["llm_reasoning_effort"].lower(),
        contrasts=parsed,
        output_dir=kwargs["output_dir"],
    )

    results = optimize(protein_df, s2c, config, ground_truth=gt_set)
    for key, state in results.items():
        best = state.best_config
        click.echo(
            f"[{key}] Best: {best.name if best else 'N/A'} "
            f"(score={state.best_score:.4f}, "
            f"experiments={state.total_experiments})"
        )
