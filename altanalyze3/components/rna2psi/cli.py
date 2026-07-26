"""Command-line interface for rna2psi (AML splicing-PSI imputation from RNA).

Subcommands:
    model-info    print bundle metadata
    predict       impute PSI from a gene-by-sample counts matrix (TSV/CSV)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DEFAULT_BUNDLE_PATH, load_bundle


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="rna2psi")
    sub = p.add_subparsers(dest="cmd", required=True)

    mi = sub.add_parser("model-info")
    mi.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))

    pc = sub.add_parser("predict")
    pc.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    pc.add_argument("--input", required=True, help="genes(ENSG) x samples counts matrix")
    pc.add_argument("--output", required=True, help="events x samples imputed PSI (TSV)")
    pc.add_argument("--sep", default="\t")
    pc.add_argument("--gene-axis", default="rows", choices=["rows", "columns"],
                    help="'rows' (default): genes are the index; 'columns': samples are the index")
    pc.add_argument("--normalized", action="store_true",
                    help="input is already cp10k+log1p / cellHarmony-lite scaled")

    args = p.parse_args(argv)
    if args.cmd == "model-info":
        print(json.dumps(load_bundle(args.bundle).model_info(), indent=2, default=str))
        return 0

    bundle = load_bundle(args.bundle)
    res = bundle.predict_from_file(args.input, sep=args.sep, gene_axis=args.gene_axis,
                                   normalized=args.normalized)
    # write events x samples to mirror the AltAnalyze PSI EventAnnotation layout
    out = res.predictions.T
    out.index.name = "UID"
    out.to_csv(args.output, sep="\t", na_rep="", float_format="%.3f")
    print(json.dumps(res.summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
