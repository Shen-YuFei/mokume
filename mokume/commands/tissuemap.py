"""
CLI command for the TissueMap pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click


@click.command("tissuemap", short_help="Per-dataset tissue proteome analysis")
@click.option(
    "--scan-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory containing QPX dataset(s). Can be a single PXD directory "
    "or a parent with PXD* subdirectories. [required unless --generate-config]",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("tissuemap_output"),
    show_default=True,
    help="Output directory for results.",
)
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="YAML configuration file. Use --generate-config to create a template.",
)
@click.option(
    "--generate-config",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write default YAML config template to this path and exit.",
)
@click.option(
    "--tmt-dataset",
    "tmt_datasets",
    multiple=True,
    default=(),
    help="Dataset ID(s) that are TMT (repeat for multiple). "
    "Auto-detected if not specified.",
)
@click.option(
    "--n-jobs",
    type=int,
    default=8,
    show_default=True,
    help="Threads for t-SNE and parallel dataset processing.",
)
@click.option(
    "--dpi",
    type=int,
    default=None,
    help="Plot resolution (DPI).  [default: 250]",
)
def tissuemap_cmd(
    scan_dir: Optional[Path],
    output_dir: Path,
    config_path: Optional[Path],
    generate_config: Optional[Path],
    tmt_datasets: tuple[str, ...],
    n_jobs: int,
    dpi: Optional[int],
) -> None:
    """Per-dataset tissue proteome analysis.

    \b
    Quick start:
      mokume tissuemap --scan-dir ./data --output-dir ./results

    \b
    Advanced (with config file):
      mokume tissuemap --generate-config tissuemap.yaml
      mokume tissuemap --scan-dir ./data --config tissuemap.yaml
    """
    from mokume.tissuemap.config import (
        InputConfig,
        OutputConfig,
        PlottingConfig,
        TissueMapConfig,
        generate_default_yaml,
        load_config,
    )
    from mokume.tissuemap.pipeline import TissueMapPipeline

    # Handle --generate-config: write template and exit
    if generate_config is not None:
        generate_default_yaml(generate_config)
        click.echo(f"Default config written to {generate_config}")
        return

    if scan_dir is None:
        raise click.UsageError("--scan-dir is required when not using --generate-config.")

    # Build config: YAML file → CLI overrides → defaults
    if config_path is not None:
        overrides: dict[str, object] = {}
        overrides["input.scan_dir"] = str(scan_dir)
        overrides["output.output_dir"] = str(output_dir)
        if tmt_datasets:
            overrides["input.tmt_datasets"] = list(tmt_datasets)
        overrides["n_jobs"] = n_jobs
        if dpi is not None:
            overrides["plotting.dpi"] = dpi
        config = load_config(config_path, overrides=overrides)
    else:
        config = TissueMapConfig(
            n_jobs=n_jobs,
            input=InputConfig(
                scan_dir=scan_dir,
                tmt_datasets=list(tmt_datasets),
            ),
            plotting=PlottingConfig(
                dpi=dpi if dpi is not None else 250,
            ),
            output=OutputConfig(
                output_dir=output_dir,
            ),
        )

    pipeline = TissueMapPipeline(config)
    pipeline.run()
