"""Configuration dataclasses for the agentic analysis module."""

from dataclasses import dataclass, field
from typing import Literal, Optional

# Provider presets: (base_url, default_model)
_PROVIDER_PRESETS: dict[str, tuple[str, str]] = {
    "deepseek": ("https://api.deepseek.com", "deepseek-v4-pro"),
}


@dataclass
class ScoreWeights:
    """Weights for composite score computation."""

    # Ground truth mode (Mode A)
    w_auc: float = 0.4
    w_tp: float = 0.4
    w_fp: float = 0.2

    # Unsupervised mode (Mode B)
    w_de: float = 0.4
    w_cv: float = 0.3
    w_miss: float = 0.3


@dataclass
class AgenticConfig:
    """Configuration for the agentic optimization loop."""

    # Input paths
    qpx_dir: Optional[str] = None
    sdrf: Optional[str] = None
    protein_matrix: Optional[str] = None

    # Ground truth (optional)
    ground_truth: Optional[str] = None
    expected_fc: Optional[str] = None

    # Optimization budget
    max_rounds: int = 5
    max_experiments: int = 30
    # Parallel execution: number of candidates to run concurrently within a
    # single round. 1 keeps the historical sequential behaviour; -1 uses all
    # cores. Joblib uses a threading backend so the PreprocessCache stays
    # shared across workers.
    n_jobs: int = 1
    # Number of contrasts to optimize concurrently. Default 1 keeps the
    # original behaviour and is also LLM-rate-limit friendly; increase to
    # parallelize multi-contrast runs (PreprocessCache is shared across
    # contrasts since (norm, imp) is contrast-independent).
    max_concurrent_contrasts: int = 1
    # Rule-based convergence: declare converged when the best score has not
    # improved by more than this delta for 2 consecutive rounds.
    convergence_score_gap: float = 0.02
    # If True, optimize() will try to load a previous state.json per contrast
    # from output_dir and resume after the last completed round.
    resume: bool = False

    # LLM settings (OpenAI-compatible API)
    use_llm: bool = True
    llm_provider: Literal["deepseek", "custom"] = "deepseek"
    llm_base_url: Optional[str] = None  # auto-set from provider if None
    llm_api_key: Optional[str] = None
    llm_model: Optional[str] = None  # auto-set from provider if None
    llm_temperature: float = 0.1  # ignored when thinking mode is enabled
    llm_thinking: bool = True  # deepseek thinking mode (V4 Pro)
    llm_reasoning_effort: str = "high"  # "high" or "max" (deepseek only)
    # Default to DeepSeek V4 Pro's documented 384K output cap so thinking
    # mode reasoning_content + tool-call arguments never hit a length cut.
    # Override on the AgenticConfig if your custom provider has a lower cap.
    llm_max_tokens: int = 384000

    def __post_init__(self) -> None:
        """Apply provider presets for base_url and model if not explicitly set."""
        preset = _PROVIDER_PRESETS.get(self.llm_provider)
        if preset:
            if self.llm_base_url is None:
                self.llm_base_url = preset[0]
            if self.llm_model is None:
                self.llm_model = preset[1]
        else:
            if self.llm_base_url is None:
                raise ValueError("llm_base_url is required for custom provider")
            if self.llm_model is None:
                raise ValueError("llm_model is required for custom provider")

    # Evaluation
    weights: ScoreWeights = field(default_factory=ScoreWeights)
    fdr_threshold: float = 0.05

    # DE settings
    contrasts: Optional[list] = None

    # Output
    output_dir: str = "./optimization"
