"""Entry point for ``python -m mokume`` and the ``mokume`` console script.

Runs the compute CLI in-process through the Rust extension (no subprocess).
clap handles help/version/usage errors with the usual exit codes.
"""

import argparse
import importlib
import sys


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
    """Dispatch the plugin service or the Rust compute CLI."""
    args = sys.argv[1:]
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
