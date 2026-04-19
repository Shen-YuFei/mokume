"""Configuration dataclasses for the agentic analysis module."""

from dataclasses import dataclass, field
from typing import Optional


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

    # LLM settings (OpenAI-compatible API)
    use_llm: bool = True
    llm_base_url: Optional[str] = None
    llm_api_key: Optional[str] = None
    llm_model: str = "deepseek-chat"
    llm_temperature: float = 0.1

    # Evaluation
    weights: ScoreWeights = field(default_factory=ScoreWeights)
    fdr_threshold: float = 0.05

    # DE settings
    contrasts: Optional[list] = None

    # Output
    output_dir: str = "./optimization"
