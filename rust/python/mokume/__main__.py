"""Entry point for ``python -m mokume`` and the ``mokume`` console script.

Runs the compute CLI in-process through the Rust extension (no subprocess).
clap handles help/version/usage errors with the usual exit codes.
"""

import importlib
import sys


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
    raise SystemExit(getattr(package, "_run_cli")(args))


if __name__ == "__main__":
    main()
