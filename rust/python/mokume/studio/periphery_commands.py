"""Studio contracts for Mokume's Python periphery commands."""

from __future__ import annotations

import contextlib
import importlib
import io
from pathlib import Path
from typing import Any


def _flag(
    name: str,
    *value_names: str,
    short: str | None = None,
    default: tuple[str, ...] = (),
    choices: tuple[str, ...] = (),
) -> dict[str, Any]:
    count = len(value_names)
    return {
        "id": name.replace("-", "_"),
        "long": name,
        "short": short,
        "help": None,
        "default": list(default),
        "possible_values": list(choices),
        "required": False,
        "repeat": False,
        "value_names": list(value_names),
        "value_arity": {"min": count, "max": count},
        "conflicts": [],
        "global": False,
    }


def _required_flag(
    name: str,
    *value_names: str,
    short: str | None = None,
    repeat: bool = False,
) -> dict[str, Any]:
    flag = _flag(name, *value_names, short=short)
    flag.update(required=True, repeat=repeat)
    return flag


def _repeat_flag(name: str, *value_names: str) -> dict[str, Any]:
    flag = _flag(name, *value_names)
    flag["repeat"] = True
    return flag


PERIPHERY_COMMAND_SPECS: tuple[dict[str, Any], ...] = (
    {
        "path": ["tissuemap"],
        "help": "Build a tissue proteome atlas [requires: mokume[tissuemap]]",
        "flags": [
            _required_flag("input", "DIR", short="i"),
            _required_flag("outdir", "DIR", short="o"),
            _flag("config", "FILE"),
            _repeat_flag("tmt-dataset", "ID"),
            _flag("threads", "N", short="t", default=("8",)),
            _flag("dpi", "N", default=("250",)),
            _flag(
                "impute-method",
                "METHOD",
                choices=(
                    "mindet",
                    "minprob",
                    "qrilc",
                    "knn",
                    "seqknn",
                    "missforest",
                    "impseq",
                    "median",
                ),
            ),
            _flag("embedding-method", "METHOD", choices=("tsne", "umap")),
        ],
    },
    {
        "path": ["plot", "pca"],
        "help": "Render a PCA-by-condition plot [requires: mokume[analysis]]",
        "flags": [
            _required_flag("protein-matrix", "FILE", short="p"),
            _required_flag("sdrf", "FILE", short="s"),
            _required_flag("output", "FILE", short="o"),
        ],
    },
    {
        "path": ["plot", "tsne"],
        "help": "Render a t-SNE plot from protein tables [requires: mokume[analysis]]",
        "flags": [
            _required_flag("input", "DIR", short="i"),
            _flag("pattern", "GLOB", short="p", default=("proteins.tsv",)),
            _required_flag("output", "FILE", short="o"),
        ],
    },
    {
        "path": ["plot", "de"],
        "help": "Render volcano plots and DE heatmaps [requires: mokume[analysis]]",
        "flags": [
            _required_flag("protein-matrix", "FILE", short="p"),
            _required_flag("outdir", "DIR", short="o"),
            _flag("sdrf", "FILE", short="s"),
            _flag("volcano", default=("false",)),
            _flag("heatmap", default=("false",)),
            _required_flag(
                "contrast",
                "KEY",
                "GROUP_A",
                "GROUP_B",
                "FILE",
                repeat=True,
            ),
            _flag("log2fc", "VALUE", default=("0.5",)),
            _flag("fdr", "FRACTION", default=("0.05",)),
            _repeat_flag("highlight-protein", "PROTEIN"),
        ],
    },
    {
        "path": ["interactive-report"],
        "help": "Build an interactive DE HTML report [requires: mokume[analysis]]",
        "flags": [
            _required_flag("protein-matrix", "FILE", short="p"),
            _required_flag("sdrf", "FILE", short="s"),
            _required_flag("output", "FILE", short="o"),
            _required_flag(
                "contrast",
                "KEY",
                "GROUP_A",
                "GROUP_B",
                "FILE",
                repeat=True,
            ),
            _flag("log2fc", "VALUE", default=("0.5",)),
            _flag("fdr", "FRACTION", default=("0.05",)),
            _repeat_flag("highlight-protein", "PROTEIN"),
        ],
    },
)

PERIPHERY_COMMAND_PATHS = {
    tuple(command["path"]) for command in PERIPHERY_COMMAND_SPECS
}

_COMMAND_MODULES = {
    ("tissuemap",): ("tissuemap", None),
    ("plot", "pca"): ("de_plots", "pca"),
    ("plot", "tsne"): ("visualize", None),
    ("plot", "de"): ("de_plots", "de"),
    ("interactive-report",): ("interactive_report", None),
}


def validate_periphery_command(argv: list[str]) -> None:
    """Parse and apply non-compute validation for one periphery command."""
    path, module, mode, args = _command_request(argv)
    parsed = _parse_args(path, module, mode, args)
    try:
        validator = _COMMAND_VALIDATORS.get(path)
        if validator is not None:
            validator(module, parsed)
    except SystemExit as exc:
        raise ValueError(_exit_message(exc, "invalid workflow parameters")) from None


def _validate_tissuemap(module, parsed) -> None:
    if parsed.scan_dir is None:
        raise ValueError("--input is required in Mokume Studio")
    if parsed.output_dir is None:
        raise ValueError("--outdir is required in Mokume Studio")
    if parsed.n_jobs is not None and parsed.n_jobs < 1:
        raise ValueError("--threads must be greater than zero")
    if parsed.dpi is not None and parsed.dpi < 1:
        raise ValueError("--dpi must be greater than zero")
    config_error = getattr(module, "_config_error")(
        getattr(module, "_build_config")(parsed)
    )
    if config_error:
        raise ValueError(config_error)


def _validate_tsne(_module, parsed) -> None:
    if "/" in parsed.pattern or "\\" in parsed.pattern:
        raise ValueError("--pattern must match files directly inside --input")


def _validate_de(module, parsed) -> None:
    _validate_contrast_keys(parsed.contrast)
    getattr(module, "_validate_de_args")(parsed)


def _validate_interactive_report(module, parsed) -> None:
    if parsed.output is None:
        raise ValueError("--output is required in Mokume Studio")
    if not parsed.contrast:
        raise ValueError("--contrast is required")
    _validate_contrast_keys(parsed.contrast)
    getattr(module, "_validate_args")(parsed)


_COMMAND_VALIDATORS = {
    ("tissuemap",): _validate_tissuemap,
    ("plot", "tsne"): _validate_tsne,
    ("plot", "de"): _validate_de,
    ("interactive-report",): _validate_interactive_report,
}


def run_periphery_command(argv: list[str]) -> None:
    """Execute one validated periphery command in the Studio worker."""
    path, module, mode, args = _command_request(argv)
    main = getattr(module, "main")
    try:
        result = main(args, mode=mode) if mode is not None else main(args)
    except SystemExit as exc:
        if exc.code in (0, None):
            return
        raise RuntimeError(_exit_message(exc, f"{' '.join(path)} failed")) from None
    if result not in (0, None):
        raise RuntimeError(f"{' '.join(path)} failed with exit code {result}")


def periphery_output_paths(argv: list[str], outputs: list[Path]) -> list[Path]:
    """Expand report output names that depend on repeated contrast keys."""
    if tuple(argv[:1]) != ("interactive-report",):
        return outputs
    output = _option_value(argv, "--output")
    keys = [
        argv[index + 1]
        for index, argument in enumerate(argv[:-1])
        if argument == "--contrast"
    ]
    if output is None or not keys:
        return outputs
    module = importlib.import_module("mokume.commands.interactive_report")
    resolve = getattr(module, "_resolve_output_html")
    return [Path(resolve(output, key, len(keys))) for key in keys]


def _command_request(argv: list[str]):
    for path in sorted(PERIPHERY_COMMAND_PATHS, key=len, reverse=True):
        if tuple(argv[: len(path)]) == path:
            module_name, mode = _COMMAND_MODULES[path]
            module = importlib.import_module(f"mokume.commands.{module_name}")
            return path, module, mode, argv[len(path) :]
    raise ValueError("command is not a Mokume periphery workflow")


def _parse_args(path, module, mode, args):
    stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(stderr):
            parser = getattr(module, "_parse_args")
            return parser(args, mode) if mode is not None else parser(args)
    except SystemExit as exc:
        message = stderr.getvalue().strip().splitlines()
        fallback = f"invalid parameters for {' '.join(path)}"
        raise ValueError(
            message[-1] if message else _exit_message(exc, fallback)
        ) from None


def _validate_contrast_keys(contrasts) -> None:
    for key, *_rest in contrasts:
        if not key or key in {".", ".."} or "/" in key or "\\" in key:
            raise ValueError("--contrast KEY must be a file-safe name")


def _option_value(argv: list[str], option: str) -> str | None:
    try:
        return argv[argv.index(option) + 1]
    except (ValueError, IndexError):
        return None


def _exit_message(exc: SystemExit, fallback: str) -> str:
    return str(exc.code) if isinstance(exc.code, str) else fallback
