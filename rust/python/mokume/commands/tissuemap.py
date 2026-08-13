#!/usr/bin/env python3
"""Per-dataset tissue proteome analysis command (``mokume.tissuemap``).

A periphery command in the mokume wheel: it re-implements no compute kernel; it
drives the tissue-atlas pipeline (PCA / t-SNE / UMAP, AdaTiSS specificity,
markers, atlas plots) over a dataset directory. Exposed in-process as
``mokume.tissuemap`` and runnable as ``python -m mokume.commands.tissuemap``.

Run requirements: the ``tissuemap`` extra (scanpy, umap-learn, combat, anndata,
matplotlib, seaborn, ...): ``pip install mokume-rs[tissuemap]``. The
``mokume.imputation`` subpackage that ``tissuemap`` depends on ships in the
wheel.

argv contract:
    --scan-dir          PATH    dataset directory (required unless --generate-config)
    --output-dir        PATH    output directory (default: tissuemap_output)
    --config            PATH    YAML configuration file
    --generate-config   PATH    write default YAML template and exit
    --tmt-dataset       STR     TMT dataset id (repeatable)
    --n-jobs            INT     threads (default: 8)
    --dpi               INT     plot DPI
    --imputation-method STR
    --embedding-method  {tsne,umap}
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional


_IMPUTATION_CHOICES = (
    "mindet",
    "minprob",
    "qrilc",
    "knn",
    "seqknn",
    "missforest",
    "impseq",
    "median",
)


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="tissuemap.py",
        description="Per-dataset tissue proteome analysis (mokume command).",
    )
    parser.add_argument("--scan-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("tissuemap_output"))
    parser.add_argument("--config", dest="config_path", type=Path, default=None)
    parser.add_argument("--generate-config", type=Path, default=None)
    parser.add_argument(
        "--tmt-dataset", dest="tmt_datasets", action="append", default=[]
    )
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=None)
    parser.add_argument(
        "--imputation-method",
        dest="imputation_method",
        type=str.lower,
        choices=_IMPUTATION_CHOICES,
        default=None,
    )
    parser.add_argument(
        "--embedding-method",
        dest="embedding_method",
        type=str.lower,
        choices=("tsne", "umap"),
        default=None,
    )
    return parser.parse_args(argv)


def _build_config(args: argparse.Namespace):
    """Mirror ``mokume.commands.tissuemap._build_config``."""
    from mokume.tissuemap.config import (
        EmbeddingConfig,
        InputConfig,
        OutputConfig,
        PlottingConfig,
        TissueMapConfig,
        load_config,
    )

    scan_dir: Path = args.scan_dir
    output_dir: Path = args.output_dir
    config_path: Optional[Path] = args.config_path
    tmt_datasets = tuple(args.tmt_datasets)
    n_jobs: int = args.n_jobs
    dpi: Optional[int] = args.dpi
    imputation_method: Optional[str] = args.imputation_method
    embedding_method: Optional[str] = args.embedding_method

    if config_path is not None:
        overrides: dict[str, object] = {
            "input.scan_dir": str(scan_dir),
            "output.output_dir": str(output_dir),
            "n_jobs": n_jobs,
        }
        if tmt_datasets:
            overrides["input.tmt_datasets"] = list(tmt_datasets)
        if dpi is not None:
            overrides["plotting.dpi"] = dpi
        if imputation_method is not None:
            overrides["embedding.imputation_method"] = imputation_method.lower()
        if embedding_method is not None:
            overrides["embedding.embedding_method"] = embedding_method.lower()
        return load_config(config_path, overrides=overrides)

    embedding_kwargs: dict[str, str] = {}
    if imputation_method is not None:
        embedding_kwargs["imputation_method"] = imputation_method.lower()
    if embedding_method is not None:
        embedding_kwargs["embedding_method"] = embedding_method.lower()

    return TissueMapConfig(
        n_jobs=n_jobs,
        input=InputConfig(scan_dir=scan_dir, tmt_datasets=list(tmt_datasets)),
        plotting=PlottingConfig(dpi=dpi if dpi is not None else 250),
        output=OutputConfig(output_dir=output_dir),
        embedding=EmbeddingConfig(**embedding_kwargs),
    )


def main(argv: list[str]) -> int:
    args = _parse_args(argv)

    try:
        from mokume.tissuemap.config import generate_default_yaml
        from mokume.tissuemap.pipeline import TissueMapPipeline
    except ImportError as exc:  # pragma: no cover - environment dependent
        print(
            "error: failed to import mokume tissuemap dependencies: "
            f"{exc}\nInstall them with: pip install mokume-rs[tissuemap]",
            file=sys.stderr,
        )
        return 1

    # `--generate-config` writes a template and exits (mirrors the click body).
    if args.generate_config is not None:
        generate_default_yaml(args.generate_config)
        print(f"Default config written to {args.generate_config}")
        return 0

    if args.scan_dir is None:
        print(
            "error: --scan-dir is required when not using --generate-config.",
            file=sys.stderr,
        )
        return 2

    if not args.scan_dir.exists() or not args.scan_dir.is_dir():
        print(
            f"error: --scan-dir does not exist or is not a directory: {args.scan_dir}",
            file=sys.stderr,
        )
        return 2

    config = _build_config(args)
    TissueMapPipeline(config).run()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
