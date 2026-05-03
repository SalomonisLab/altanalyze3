"""Command-line interface for rna2adt.

Subcommands:
    model-info     Print the bundle metadata.
    predict-h5ad   Run the loaded bundle on an .h5ad and dump predicted ADTs.
    predict-csv    Run the loaded bundle on a sample-by-gene CSV/TSV.

Training and benchmarking live in ``rna2adt.training`` and have their own
``__main__`` entrypoint.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DEFAULT_BUNDLE_PATH, load_bundle


def _write_outputs(predictions: pd.DataFrame, output_path: Path, summary, summary_path: Path | None) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_path)
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n")


def _add_bundle_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bundle", default=str(DEFAULT_BUNDLE_PATH),
                        help="Path to the rna2adt bundle pickle.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RNA -> ADT (CITE-seq) imputation")
    sub = parser.add_subparsers(dest="command", required=True)

    info = sub.add_parser("model-info", help="Print bundle metadata.")
    _add_bundle_arg(info)

    pcsv = sub.add_parser("predict-csv", help="Predict ADTs from a sample-by-gene CSV/TSV.")
    _add_bundle_arg(pcsv)
    pcsv.add_argument("--input", required=True)
    pcsv.add_argument("--output", required=True)
    pcsv.add_argument("--sep", default=",")
    pcsv.add_argument("--transpose", action="store_true")
    pcsv.add_argument("--summary-json")

    ph5 = sub.add_parser("predict-h5ad", help="Predict ADTs from a single-cell h5ad.")
    _add_bundle_arg(ph5)
    ph5.add_argument("--input", required=True)
    ph5.add_argument("--output", required=True)
    ph5.add_argument("--layer")
    ph5.add_argument("--gene-symbol-col")
    ph5.add_argument("--groupby")
    ph5.add_argument("--chunk-size", type=int, default=4096)
    ph5.add_argument("--summary-json")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "model-info":
        bundle = load_bundle(args.bundle)
        print(json.dumps(bundle.model_info(), indent=2, sort_keys=True, default=str))
        return 0

    if args.command == "predict-csv":
        bundle = load_bundle(args.bundle)
        expression = pd.read_csv(args.input, sep=args.sep, index_col=0)
        if args.transpose:
            expression = expression.T
        result = bundle.predict_from_dataframe(expression)
        _write_outputs(result.predictions, Path(args.output), result.summary,
                       Path(args.summary_json) if args.summary_json else None)
        print(json.dumps(result.summary, indent=2, sort_keys=True, default=str))
        return 0

    if args.command == "predict-h5ad":
        bundle = load_bundle(args.bundle)
        result = bundle.predict_from_h5ad(
            args.input,
            layer=args.layer,
            gene_symbol_col=args.gene_symbol_col,
            groupby=args.groupby,
            chunk_size=args.chunk_size,
        )
        _write_outputs(result.predictions, Path(args.output), result.summary,
                       Path(args.summary_json) if args.summary_json else None)
        print(json.dumps(result.summary, indent=2, sort_keys=True, default=str))
        return 0

    parser.error(f"Unhandled command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
