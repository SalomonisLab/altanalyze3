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
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run RNA-to-lipid predictions from tabular expression matrices or h5ad files."
    )
    parser.add_argument(
        "--bundle",
        default=str(DEFAULT_BUNDLE_PATH),
        help="Path to the bundled scikit-learn model pickle.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("model-info", help="Print bundle metadata.")

    predict_csv = subparsers.add_parser("predict-csv", help="Predict lipids from a sample-by-gene CSV/TSV.")
    predict_csv.add_argument("--input", required=True, help="Path to the expression matrix.")
    predict_csv.add_argument("--output", required=True, help="Where to write predicted lipids as CSV.")
    predict_csv.add_argument(
        "--summary-json",
        help="Optional JSON path for run metadata.",
    )
    predict_csv.add_argument(
        "--sep",
        default=",",
        help="Field separator for the input expression matrix. Default: ','.",
    )
    predict_csv.add_argument(
        "--transpose",
        action="store_true",
        help="Transpose the input matrix before prediction.",
    )

    predict_h5ad = subparsers.add_parser("predict-h5ad", help="Predict lipids from a single-cell h5ad file.")
    predict_h5ad.add_argument("--input", required=True, help="Path to the input h5ad.")
    predict_h5ad.add_argument("--output", required=True, help="Where to write predicted lipids as CSV.")
    predict_h5ad.add_argument(
        "--summary-json",
        help="Optional JSON path for run metadata.",
    )
    predict_h5ad.add_argument(
        "--layer",
        help="Optional AnnData layer to use instead of X.",
    )
    predict_h5ad.add_argument(
        "--gene-symbol-col",
        help="Optional adata.var column containing gene symbols. Defaults to var_names.",
    )
    predict_h5ad.add_argument(
        "--groupby",
        help="Optional adata.obs column name used to average predicted lipids by group.",
    )
    predict_h5ad.add_argument(
        "--chunk-size",
        type=int,
        default=2048,
        help="Number of cells to process per chunk. Default: 2048.",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    bundle = load_bundle(args.bundle)

    if args.command == "model-info":
        print(json.dumps(bundle.model_info(), indent=2, sort_keys=True))
        return 0

    if args.command == "predict-csv":
        expression = pd.read_csv(args.input, sep=args.sep, index_col=0)
        if args.transpose:
            expression = expression.T
        result = bundle.predict_from_dataframe(expression)
        _write_outputs(
            result.predictions,
            Path(args.output),
            result.summary,
            Path(args.summary_json) if args.summary_json else None,
        )
        print(json.dumps(result.summary, indent=2, sort_keys=True))
        return 0

    if args.command == "predict-h5ad":
        result = bundle.predict_from_h5ad(
            args.input,
            layer=args.layer,
            gene_symbol_col=args.gene_symbol_col,
            groupby=args.groupby,
            chunk_size=args.chunk_size,
        )
        _write_outputs(
            result.predictions,
            Path(args.output),
            result.summary,
            Path(args.summary_json) if args.summary_json else None,
        )
        print(json.dumps(result.summary, indent=2, sort_keys=True))
        return 0

    parser.error(f"Unhandled command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
