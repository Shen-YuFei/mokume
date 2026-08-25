"""Unified entry point for the Rust-backed ``mokume`` wheel."""

import argparse
import importlib
import sys


_PERIPHERY_COMMANDS = {
    "tsne-visualization": "visualize",
    "tissuemap": "tissuemap",
    "de-plots": "de_plots",
    "interactive-report": "interactive_report",
}

_ROOT_HELP = """Mokume proteomics quantification and analysis toolkit

Usage:
  mokume [NATIVE OPTIONS] <COMPUTE COMMAND> [COMMAND OPTIONS]
  mokume <PERIPHERY COMMAND> [COMMAND OPTIONS]

Rust-native compute commands:
  features2proteins  Quantify proteins from a QPX feature parquet file
  features2peptides  Convert features to peptide-level output
  peptides2protein   Compute protein quantities from peptide-level input
  correct-batches    Correct batch effects in protein quantification output

Optional analysis and visualization commands:
  tsne-visualization  Render a t-SNE plot [requires: mokume[plotting]]
  tissuemap           Build a tissue proteome atlas [requires: mokume[tissuemap]]
  de-plots            Render DE plots [requires: mokume[plotting]]
  interactive-report  Build a DE HTML report [requires: mokume[reports]]

Other commands:
  help                Print this message or command-specific help

Native compute options:
  -v, --log-level <LEVEL>  [default: debug] [possible values: debug, info, warn]
      --log-file <PATH>

General options:
  -h, --help     Print help
  -V, --version  Print version

Run `mokume <COMMAND> --help` for command-specific options.
"""


def _print_root_help(stream):
    print(_ROOT_HELP, file=stream, end="")


def _periphery_request(args):
    if args and args[0] in _PERIPHERY_COMMANDS:
        return args[0], args[1:]
    if len(args) == 2 and args[0] == "help" and args[1] in _PERIPHERY_COMMANDS:
        return args[1], ["--help"]
    return None


def _run_periphery(command, args):
    module = importlib.import_module(f"mokume.commands.{_PERIPHERY_COMMANDS[command]}")
    result = getattr(module, "main")(args)
    return 0 if result is None else result


def _render_requested_pibaq_qc(args, package):
    """Render the QC PDF after the native command has written its piBAQ table."""
    try:
        command_index = args.index("peptides2protein")
    except ValueError:
        return
    command_args = args[command_index + 1 :]
    qc_requested = any(
        option == "--verbose"
        or option in {"--qc_report", "--qc-report"}
        or option.startswith(("--qc_report=", "--qc-report="))
        for option in command_args
    )
    if not qc_requested:
        return
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--method", default="pibaq")
    parser.add_argument("--qc_report", "--qc-report", default="QCprofile.pdf")
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--tpa", action="store_true")
    parser.add_argument("--ruler", action="store_true")
    parsed, _ = parser.parse_known_args(command_args)
    if parsed.method.lower() != "pibaq":
        return
    package.peptides2protein_qc(
        protein_table=parsed.output,
        qc_report=parsed.qc_report,
        plot_column="PiBAQPpb" if parsed.normalize else "PiBAQ",
        tpa=parsed.tpa,
        ruler=parsed.ruler,
    )


def main():
    """Dispatch wheel periphery, plugin service, or Rust-native computation."""
    args = sys.argv[1:]
    if not args:
        _print_root_help(sys.stderr)
        raise SystemExit(2)
    if args in (["-h"], ["--help"], ["help"]):
        _print_root_help(sys.stdout)
        raise SystemExit(0)
    periphery = _periphery_request(args)
    if periphery is not None:
        command, command_args = periphery
        raise SystemExit(_run_periphery(command, command_args))
    if args[:2] == ["mcp", "serve"]:
        module = importlib.import_module("mokume.agentic.mcp_server")
        try:
            raise SystemExit(getattr(module, "main")(args[2:]))
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from None
    if args and args[0] == "mcp":
        raise SystemExit("Usage: mokume mcp serve --knowledge PATH")
    package = importlib.import_module("mokume")
    code = getattr(package, "_run_cli")(args)
    if code == 0:
        try:
            _render_requested_pibaq_qc(args, package)
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from None
    raise SystemExit(code)


if __name__ == "__main__":
    main()
