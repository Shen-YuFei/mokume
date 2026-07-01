"""Entry point for ``python -m mokume`` and the ``mokume`` console script.

Runs the compute CLI in-process through the Rust extension (no subprocess), so
``pip install mokume`` provides the same ``mokume`` command the standalone Rust
binary does. clap handles help/version/usage errors with the usual exit codes.
"""

import sys

from mokume._mokume import run_cli


def main():
    """Console-script entry: dispatch argv to the Rust CLI and exit with its code."""
    raise SystemExit(run_cli(sys.argv[1:]))


if __name__ == "__main__":
    main()
