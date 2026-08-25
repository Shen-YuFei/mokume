#!/usr/bin/env python3
"""Native-backed piBAQ compatibility command.

The normal Rust-backed ``peptides2protein`` and ``features2proteins`` piBAQ paths
support every protease registered in the installed pyOpenMS ``ProteaseDB``.
This module preserves the established file-oriented Python command contract and
forwards it through ``mokume.quantification.pibaq.peptides_to_protein``. That
compatibility API delegates shared-peptide allocation, theoretical denominators,
evidence classification, TPA, and output generation to the same native Rust core
as ``mokume quantify peptides2protein``; it is neither an enzyme fallback nor a separate
full-Python implementation.

The output uses the native command's TSV schema, so downstream consumers see the
same contract through either entry point.

argv contract (runnable via ``python -m mokume.commands.peptides2protein_pibaq``):

    python -m mokume.commands.peptides2protein_pibaq \
        --peptides <peptides.csv|.tsv|.parquet> \
        --fasta <proteome.fasta> \
        --enzyme <protease name> \
        --output <proteins.tsv> \
        --min-aa 7 --max-aa 30 \
        --ploidy 2 --organism human --cpc 200 \
        --min-shared 2 --min-anchors 1 --high-anchor-threshold 3 \
        [--normalize] [--tpa] [--ruler] [--verbose] \
        [--qc-report <QCprofile.pdf>] [--families <families.yaml>]

The DirectLFQ knobs (``--threads`` / ``--directlfq-min-nonan``) are not forwarded:
``peptides_to_protein`` takes no such parameters, and only the piBAQ method
reaches this command. TopN needs no knob at all -- its N is spelled in the
method name (``top5``).

Exit-code convention (see python/README.md):
    0   success
    1   a runtime/import failure (e.g. mokume not installed, bad enzyme)
    2   bad arguments (argparse)

Exit code 0 on success; non-zero (with a message on stderr) otherwise.
"""

import argparse
import importlib
import sys


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Compute piBAQ through the native-backed compatibility API."
    )
    parser.add_argument(
        "--peptides",
        required=True,
        help="Peptide intensity table (CSV / TSV / parquet).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Protein FASTA used to derive theoretical peptide counts.",
    )
    parser.add_argument(
        "--enzyme",
        required=True,
        help="Digestion enzyme name understood by pyOpenMS ProteaseDigestion.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination protein table (tab-separated, or parquet if .parquet).",
    )
    parser.add_argument("--min-aa", type=int, default=7)
    parser.add_argument("--max-aa", type=int, default=30)
    parser.add_argument("--ploidy", type=int, default=None)
    parser.add_argument("--organism", default=None)
    parser.add_argument("--cpc", type=float, default=None)
    parser.add_argument(
        "--qc-report",
        default=None,
        help="Render a QC PDF at this path; --verbose is not required.",
    )
    parser.add_argument(
        "--families",
        dest="families",
        default=None,
        help="Optional YAML declaring explicit piBAQ family overrides.",
    )
    parser.add_argument("--min-shared", type=int, default=2)
    parser.add_argument("--min-anchors", type=int, default=1)
    parser.add_argument("--high-anchor-threshold", type=int, default=3)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--tpa", action="store_true")
    parser.add_argument("--ruler", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def _validate_options(args):
    if args.ruler and not args.tpa:
        raise SystemExit("piBAQ command aborted: --ruler requires --tpa")
    if not args.ruler and any(
        value is not None for value in (args.ploidy, args.organism, args.cpc)
    ):
        raise SystemExit(
            "piBAQ command aborted: --ploidy/--organism/--cpc require --ruler"
        )
    if args.ploidy is not None and args.ploidy < 1:
        raise SystemExit("piBAQ command aborted: --ploidy must be greater than zero")
    if args.cpc is not None and args.cpc <= 0:
        raise SystemExit("piBAQ command aborted: --cpc must be greater than zero")


def _resolve_options(args):
    _validate_options(args)
    ploidy = 2 if args.ploidy is None else args.ploidy
    organism = "human" if args.organism is None else args.organism
    cpc = 200.0 if args.cpc is None else args.cpc
    qc_report = "QCprofile.pdf" if args.qc_report is None else args.qc_report
    return ploidy, organism, cpc, qc_report


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    ploidy, organism, cpc, qc_report = _resolve_options(args)

    # Import lazily so a missing install yields an actionable message (exit 1)
    # rather than an opaque traceback before argparse even runs.
    try:
        peptides_to_protein = getattr(
            importlib.import_module("mokume.quantification.pibaq"),
            "peptides_to_protein",
        )
    except (ImportError, AttributeError) as exc:
        raise SystemExit(
            "piBAQ command aborted: the mokume package could not be imported "
            f"({exc}). Install the Rust-backed distribution with: pip install mokume"
        ) from exc

    # Mirror the Python CLI ``peptides2protein`` piBAQ branch exactly: it forwards
    # every option straight into ``peptides_to_protein``. The enzyme validity
    # (and any organism / ruler guard) is enforced inside mokume, so a bad enzyme
    # surfaces here as the same error the Python CLI would raise.
    try:
        peptides_to_protein(
            fasta=args.fasta,
            peptides=args.peptides,
            enzyme=args.enzyme,
            normalize=args.normalize,
            min_aa=args.min_aa,
            max_aa=args.max_aa,
            tpa=args.tpa,
            ruler=args.ruler,
            ploidy=ploidy,
            cpc=cpc,
            organism=organism,
            output=args.output,
            verbose=args.verbose or args.qc_report is not None,
            qc_report=qc_report,
            families_yaml=args.families,
            min_shared=args.min_shared,
            min_anchors=args.min_anchors,
            high_anchor_threshold=args.high_anchor_threshold,
        )
    except (ValueError, KeyError, RuntimeError) as exc:
        # Bad enzyme / organism / ruler-guard failures arrive here. Re-raise as a
        # SystemExit so the caller sees a non-zero exit with a clear
        # message rather than a raw traceback.
        raise SystemExit(f"piBAQ command aborted: {exc}") from exc

    print(f"piBAQ protein table written to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
