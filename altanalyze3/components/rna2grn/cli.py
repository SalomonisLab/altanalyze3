"""Command-line interface for rna2grn (GRN imputation from RNA).

Subcommands:
    model-info     print bundle metadata
    predict-csv    predict GRNs from a pseudobulk-by-gene CSV/TSV
    predict-h5ad   predict GRNs from an h5ad (optionally grouped into pseudobulks)
    predict-10x    predict GRNs from a 10x filtered_feature_bc_matrix.h5
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DEFAULT_BUNDLE_PATH, Rna2GrnBundle


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


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="rna2grn")
    sub = p.add_subparsers(dest="cmd", required=True)

    mi = sub.add_parser("model-info")
    mi.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))

    pc = sub.add_parser("predict-csv")
    pc.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    pc.add_argument("--input", required=True)
    pc.add_argument("--output", required=True)
    pc.add_argument("--sep", default=",")
    pc.add_argument("--transpose", action="store_true", help="input is gene-by-pseudobulk")
    pc.add_argument("--normalized", action="store_true", help="input is already CP10k+log1p")
    pc.add_argument("--neighbors", action="store_true")
    pc.add_argument("--summary-json", default=None)

    ph = sub.add_parser("predict-h5ad")
    ph.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    ph.add_argument("--input", required=True)
    ph.add_argument("--output", required=True)
    ph.add_argument("--groupby", default=None)
    ph.add_argument("--layer", default=None)
    ph.add_argument("--gene-symbol-col", default=None)
    ph.add_argument("--min-cells", type=int, default=1)
    ph.add_argument("--neighbors", action="store_true")
    ph.add_argument("--summary-json", default=None)

    px = sub.add_parser("predict-10x")
    px.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    px.add_argument("--input", required=True)
    px.add_argument("--output", required=True)
    px.add_argument("--groupby", default=None)
    px.add_argument("--neighbors", action="store_true")
    px.add_argument("--summary-json", default=None)

    args = p.parse_args(argv)
    if args.cmd == "model-info":
        print(json.dumps(Rna2GrnBundle.load(args.bundle).model_info(), indent=2, default=str))
        return 0

    bundle = Rna2GrnBundle.load(args.bundle)
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
