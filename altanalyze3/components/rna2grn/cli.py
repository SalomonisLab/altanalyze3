"""Command-line interface for rna2grn (GRN imputation from RNA).

Subcommands:
    list-references print the named references shipped with the module
    model-info      print bundle metadata
    predict-csv     predict GRNs from a pseudobulk-by-gene CSV/TSV
    predict-h5ad    predict GRNs from an h5ad (optionally grouped into pseudobulks)
    predict-10x     predict GRNs from a 10x filtered_feature_bc_matrix.h5

Every predict subcommand takes ``--reference {leukemia,lung}`` (default: leukemia)
or an explicit ``--bundle <path>``; passing both is an error.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import (DEFAULT_REFERENCE, REFERENCE_BUNDLES, Rna2GrnBundle,
                  available_references)


def _write(result, output: str, summary_json: str | None):
    result.predictions.to_csv(output)
    summ = {k: v for k, v in result.summary.items() if k != "neighbors"}
    if "neighbors" in result.summary:
        nbr_path = str(Path(output).with_suffix("")) + ".neighbors.csv"
        result.summary["neighbors"].to_csv(nbr_path)
        summ["neighbors_csv"] = nbr_path
    print(json.dumps(summ, indent=2, default=str))
    if summary_json:
        Path(summary_json).write_text(json.dumps(summ, indent=2, default=str))


def _add_bundle_args(sp):
    sp.add_argument("--bundle", default=None,
                    help="explicit bundle path (mutually exclusive with --reference)")
    sp.add_argument("--reference", default=None,
                    choices=sorted(REFERENCE_BUNDLES),
                    help=f"named reference (default: {DEFAULT_REFERENCE})")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="rna2grn")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list-references")

    mi = sub.add_parser("model-info")
    _add_bundle_args(mi)

    pc = sub.add_parser("predict-csv")
    _add_bundle_args(pc)
    pc.add_argument("--input", required=True)
    pc.add_argument("--output", required=True)
    pc.add_argument("--sep", default=",")
    pc.add_argument("--transpose", action="store_true", help="input is gene-by-pseudobulk")
    pc.add_argument("--normalized", action="store_true", help="input is already CP10k+log1p")
    pc.add_argument("--neighbors", action="store_true")
    pc.add_argument("--summary-json", default=None)

    ph = sub.add_parser("predict-h5ad")
    _add_bundle_args(ph)
    ph.add_argument("--input", required=True)
    ph.add_argument("--output", required=True)
    ph.add_argument("--groupby", default=None)
    ph.add_argument("--layer", default=None)
    ph.add_argument("--gene-symbol-col", default=None)
    ph.add_argument("--min-cells", type=int, default=1)
    ph.add_argument("--neighbors", action="store_true")
    ph.add_argument("--summary-json", default=None)

    px = sub.add_parser("predict-10x")
    _add_bundle_args(px)
    px.add_argument("--input", required=True)
    px.add_argument("--output", required=True)
    px.add_argument("--groupby", default=None)
    px.add_argument("--neighbors", action="store_true")
    px.add_argument("--summary-json", default=None)

    args = p.parse_args(argv)
    if args.cmd == "list-references":
        refs = {n: dict(path=pth, present=Path(pth).exists(),
                        default=(n == DEFAULT_REFERENCE))
                for n, pth in available_references().items()}
        print(json.dumps(refs, indent=2))
        return 0
    if args.cmd == "model-info":
        b = Rna2GrnBundle.load(args.bundle, reference=args.reference)
        print(json.dumps(b.model_info(), indent=2, default=str))
        return 0

    bundle = Rna2GrnBundle.load(args.bundle, reference=args.reference)
    if args.cmd == "predict-csv":
        df = pd.read_csv(args.input, sep=args.sep, index_col=0)
        if args.transpose:
            df = df.T
        res = bundle.predict_from_dataframe(df, normalized=args.normalized,
                                            return_neighbors=args.neighbors)
        _write(res, args.output, args.summary_json)
    elif args.cmd == "predict-h5ad":
        res = bundle.predict_from_h5ad(args.input, groupby=args.groupby, layer=args.layer,
                                       gene_symbol_col=args.gene_symbol_col,
                                       min_cells=args.min_cells, return_neighbors=args.neighbors)
        _write(res, args.output, args.summary_json)
    elif args.cmd == "predict-10x":
        res = bundle.predict_from_10x_h5(args.input, groupby=args.groupby,
                                         return_neighbors=args.neighbors)
        _write(res, args.output, args.summary_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
