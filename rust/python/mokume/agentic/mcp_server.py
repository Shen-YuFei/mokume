"""Local stdio MCP server for the installable Mokume Plugin."""

from __future__ import annotations

import argparse
import importlib
from pathlib import Path
from typing import Any


def create_server(knowledge: str):
    """Build the MCP server around one immutable knowledge snapshot."""
    try:
        fast_mcp = getattr(importlib.import_module("mcp.server.fastmcp"), "FastMCP")
        service_module = importlib.import_module("mokume.agentic.service")
        service_type = getattr(service_module, "RecommendationService")
        request_type = getattr(service_module, "EvaluationRequest")
    except ImportError as exc:
        raise RuntimeError(
            "Mokume MCP dependencies are missing; install 'mokume[agentic]'"
        ) from exc

    server = fast_mcp("mokume")
    service = service_type(knowledge)

    @server.tool()
    def inspect_dataset(
        protein_matrix: str,
        sdrf: str,
        metadata: dict[str, str | None] | None = None,
    ) -> dict[str, Any]:
        """Profile a protein matrix and return policy-filtered Mokume evidence."""
        return service.inspect_dataset(
            protein_matrix,
            sdrf,
            metadata,
        )

    @server.tool()
    def evaluate_recommendation(
        protein_matrix: str,
        sdrf: str,
        contrast: list[str],
        recommendation: dict[str, Any],
        options: dict[str, Any],
    ) -> dict[str, Any]:
        """Validate and run a scoped candidate block with the Rust kernel."""
        runtime_options = dict(options)
        output_dir = runtime_options.pop("output_dir", None)
        if not isinstance(output_dir, str) or not output_dir:
            raise ValueError("options.output_dir must be a non-empty string")
        return service.evaluate_recommendation(
            request_type(
                protein_matrix,
                sdrf,
                contrast,
                recommendation,
                output_dir,
                runtime_options,
            )
        )

    return server


def main(argv: list[str] | None = None) -> int:
    """Run the plugin MCP server over stdio."""
    parser = argparse.ArgumentParser(prog="mokume mcp serve")
    parser.add_argument("--knowledge", required=True)
    args = parser.parse_args(argv)
    knowledge = Path(args.knowledge).expanduser().resolve()
    if not knowledge.is_file():
        parser.error(f"knowledge catalog not found: {knowledge}")
    create_server(str(knowledge)).run(transport="stdio")
    return 0
