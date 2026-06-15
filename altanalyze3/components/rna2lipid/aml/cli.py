"""Command-line interface for the AML rna2lipid model (lipid imputation from RNA).

Separate from the lung/IPF rna2lipid CLI. Subcommands:
    model-info      print bundle metadata
    predict-csv     impute lipids from a samples-by-gene CSV/TSV
    predict-h5ad    impute lipids from an h5ad (optionally grouped into pseudobulks)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DEFAULT_BUNDLE_PATH, load_bundle


def _summary(summary):
    s = {k: v for k, v in summary.items() if k != "n_cells_per_group"}
    print(json.dumps(s, indent=2, default=str))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="rna2lipid-aml")
    sub = p.add_subparsers(dest="cmd", required=True)

    mi = sub.add_parser("model-info")
    mi.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))

    pc = sub.add_parser("predict-csv")
    pc.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    pc.add_argument("--input", required=True)
    pc.add_argument("--output", required=True)
    pc.add_argument("--sep", default=",")
    pc.add_argument("--transpose", action="store_true", help="input is gene-by-sample")
    pc.add_argument("--normalized", action="store_true", help="input is already CP10k+log1p")

    ph = sub.add_parser("predict-h5ad")
    ph.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH))
    ph.add_argument("--input", required=True)
    ph.add_argument("--output", required=True, help="CSV of imputed abundance (samples x molecules)")
    ph.add_argument("--h5ad-out", default=None, help="also write an AnnData (X=imputed, obs carried over)")
    ph.add_argument("--groupby", default=None)
    ph.add_argument("--layer", default=None)
    ph.add_argument("--gene-symbol-col", default=None)
    ph.add_argument("--min-cells", type=int, default=1)
    ph.add_argument("--normalized", action="store_true")

    args = p.parse_args(argv)
    if args.cmd == "model-info":
        print(json.dumps(load_bundle(args.bundle).model_info(), indent=2, default=str))
        return 0

    bundle = load_bundle(args.bundle)
    if args.cmd == "predict-csv":
        df = pd.read_csv(args.input, sep=args.sep, index_col=0)
        if args.transpose:
            df = df.T
        res = bundle.predict_from_dataframe(df, normalized=args.normalized)
        res.predictions.to_csv(args.output)
        _summary(res.summary)
    elif args.cmd == "predict-h5ad":
        import anndata as ad
        adata = ad.read_h5ad(args.input)
        res = bundle.predict_from_adata(
            adata, groupby=args.groupby, layer=args.layer,
            gene_symbol_col=args.gene_symbol_col, min_cells=args.min_cells,
            normalized=args.normalized)
        res.predictions.to_csv(args.output)
        if args.h5ad_out:
            out = bundle.impute_anndata(
                adata, groupby=args.groupby, layer=args.layer,
                gene_symbol_col=args.gene_symbol_col, min_cells=args.min_cells,
                normalized=args.normalized)
            out.write_h5ad(args.h5ad_out)
            res.summary["h5ad_out"] = args.h5ad_out
        _summary(res.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
