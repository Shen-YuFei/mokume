#!/usr/bin/env python3
"""Pure-Python iBAQ command for digestion enzymes outside the Rust-ported set.

The Rust ``peptides2protein`` iBAQ path natively computes the protein table for
the enzymes whose pyOpenMS cleavage rules are ported (Trypsin[/P], Lys-C[/P],
Arg-C[/P], Chymotrypsin[/P], Glu-C, Asp-N, Lys-N, PepsinA). pyOpenMS knows many
more proteases (CNBr, V8-DE, unspecific cleavage, ...); for those the Rust kernel
has no cleavage rule, so the mokume wheel computes the whole iBAQ table here in
pure Python, calling ``mokume.quantification.ibaq.peptides_to_protein`` exactly
as the upstream Python CLI ``mokume peptides2protein --method ibaq`` would.

Design (same "copy-py" provenance as ``peptides2protein_qc.py``): no compute is
duplicated in Rust here -- the command calls the first-class
``mokume.quantification.ibaq`` and forwards the original argv. The output is ``res.to_csv(output, sep='\t')``
inside ``peptides_to_protein``, i.e. byte-for-byte the same TSV schema the Rust
native path writes (ProteinName, SampleID, Condition, NormIntensity, Ibaq,
FamilyId, EvidenceLevel, FamilySize, plus the optional TPA / normalize / ruler
columns), so a downstream consumer cannot tell which branch produced it.

argv contract (runnable via ``python -m mokume.commands.peptides2protein_ibaq``):

    python -m mokume.commands.peptides2protein_ibaq \
        --peptides <peptides.csv|.tsv|.parquet> \
        --fasta <proteome.fasta> \
        --enzyme <protease name> \
        --output <proteins.tsv> \
        --min-aa 7 --max-aa 30 \
        --ploidy 2 --organism human --cpc 200 \
        --min-shared 2 --min-anchors 1 --high-anchor-threshold 3 \
        [--normalize] [--tpa] [--ruler] [--verbose] \
        [--qc-report <QCprofile.pdf>] [--families <families.yaml>]

The TopN / DirectLFQ knobs (``--topn_n`` / ``--threads`` / ``--min_nonan``) are
not forwarded: ``peptides_to_protein`` takes no such parameters, and only the
iBAQ method reaches this command.

Exit-code convention (see python/README.md):
    0   success
    1   a runtime/import failure (e.g. mokume not installed, bad enzyme)
    2   bad arguments (argparse)

Exit code 0 on success; non-zero (with a message on stderr) otherwise.
"""

import argparse
import sys


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description=(
            "Compute the peptides2protein iBAQ table in pure Python via "
            "mokume.quantification.ibaq, for enzymes outside the "
            "Rust-ported set."
        )
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
    parser.add_argument("--ploidy", type=int, default=2)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--cpc", type=float, default=200.0)
    parser.add_argument(
        "--qc-report",
        default="QCprofile.pdf",
        help="PDF for the verbose QC images (only used with --verbose).",
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


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    # Import lazily so a missing install yields an actionable message (exit 1)
    # rather than an opaque traceback before argparse even runs.
    try:
        from mokume.quantification.ibaq import peptides_to_protein
    except ImportError as exc:
        raise SystemExit(
            "iBAQ command aborted: the mokume package could not be imported "
            "({0}). Install its third-party dependencies with: "
            "pip install mokume[ibaq]".format(exc)
        )

    # Mirror the Python CLI ``peptides2protein`` iBAQ branch exactly: it forwards
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
            ploidy=args.ploidy,
            cpc=args.cpc,
            organism=args.organism,
            output=args.output,
            verbose=args.verbose,
            qc_report=args.qc_report,
            families_yaml=args.families,
            min_shared=args.min_shared,
            min_anchors=args.min_anchors,
            high_anchor_threshold=args.high_anchor_threshold,
        )
    except (ValueError, KeyError, RuntimeError) as exc:
        # Bad enzyme / organism / ruler-guard failures arrive here. Re-raise as a
        # SystemExit so the caller sees a non-zero exit with a clear
        # message rather than a raw traceback.
        raise SystemExit("iBAQ command aborted: {0}".format(exc))

    print("iBAQ protein table written to {0}".format(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
