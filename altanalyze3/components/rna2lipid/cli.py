from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DEFAULT_BUNDLE_PATH, load_bundle
from .evaluation import run_evaluation
from .pipeline import DEFAULT_CONFIG_PATH
from .single_cell_benchmark import benchmark_single_cell_assignments
from .training import run_training


def _write_outputs(predictions: pd.DataFrame, output_path: Path, summary, summary_path: Path | None) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_path)
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def _add_bundle_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--bundle",
        default=str(DEFAULT_BUNDLE_PATH),
        help="Path to the scikit-learn model pickle.",
    )


def _add_config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="Path to the JSON training/evaluation configuration.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run RNA-to-lipid predictions from tabular expression matrices or h5ad files."
    )
    parser.add_argument(
        "--bundle",
        default=str(DEFAULT_BUNDLE_PATH),
        help="Path to the bundled scikit-learn model pickle.",
    )
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="Path to the JSON training/evaluation configuration.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    model_info_parser = subparsers.add_parser("model-info", help="Print bundle metadata.")
    _add_bundle_arg(model_info_parser)

    train_parser = subparsers.add_parser("train", help="Train a reproducible rna2lipid bundle from the configured inputs.")
    _add_config_arg(train_parser)
    train_parser.add_argument(
        "--output-bundle",
        help="Optional override for the trained bundle pickle path.",
    )
    train_parser.add_argument(
        "--output-dir",
        help="Optional override for the training report directory.",
    )

    evaluate_parser = subparsers.add_parser("evaluate", help="Run donor-held-out evaluation and optional external structure checks.")
    _add_bundle_arg(evaluate_parser)
    _add_config_arg(evaluate_parser)
    evaluate_parser.add_argument(
        "--eval-bundle",
        help="Bundle path to use for the external pediatric evaluation. Defaults to --output-bundle if provided, otherwise --bundle.",
    )
    evaluate_parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where evaluation outputs should be written.",
    )
    evaluate_parser.add_argument(
        "--max-folds",
        type=int,
        help="Optional cap on the number of donor-held-out folds, useful for smoke testing.",
    )
    evaluate_parser.add_argument(
        "--skip-external",
        action="store_true",
        help="Skip the external pediatric structure-preservation evaluation.",
    )

    benchmark_parser = subparsers.add_parser(
        "benchmark-single-cell",
        help="Benchmark single-cell lipid assignment strategies on a 10x .h5 or .h5ad input.",
    )
    _add_bundle_arg(benchmark_parser)
    benchmark_parser.add_argument("--input", required=True, help="Input single-cell .h5 or .h5ad file.")
    benchmark_parser.add_argument("--output-dir", required=True, help="Directory where benchmark outputs should be written.")
    benchmark_parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed used for PCA/UMAP/kmeans steps.",
    )
    benchmark_parser.add_argument(
        "--knn-alpha",
        type=float,
        default=0.5,
        help="Weight retained from the original per-cell prediction during kNN smoothing.",
    )
    benchmark_parser.add_argument(
        "--coarse-clusters",
        type=int,
        help="Optional override for the number of coarse RNA kmeans clusters.",
    )
    benchmark_parser.add_argument(
        "--metacell-target-size",
        type=int,
        default=50,
        help="Approximate cells per metacell for the metacell assignment benchmark.",
    )

    predict_csv = subparsers.add_parser("predict-csv", help="Predict lipids from a sample-by-gene CSV/TSV.")
    _add_bundle_arg(predict_csv)
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
    _add_bundle_arg(predict_h5ad)
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

    if args.command == "model-info":
        bundle = load_bundle(args.bundle)
        print(json.dumps(bundle.model_info(), indent=2, sort_keys=True))
        return 0

    if args.command == "train":
        summary = run_training(
            config_path=args.config,
            bundle_path=args.output_bundle,
            output_dir=args.output_dir,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if args.command == "evaluate":
        evaluation_bundle = args.eval_bundle or args.bundle
        summary = run_evaluation(
            config_path=args.config,
            bundle_path=evaluation_bundle,
            output_dir=args.output_dir,
            max_folds=args.max_folds,
            skip_external=args.skip_external,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if args.command == "benchmark-single-cell":
        summary = benchmark_single_cell_assignments(
            input_path=args.input,
            output_dir=args.output_dir,
            bundle_path=args.bundle,
            random_state=args.random_state,
            knn_alpha=args.knn_alpha,
            coarse_clusters=args.coarse_clusters,
            metacell_target_size=args.metacell_target_size,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if args.command == "predict-csv":
        bundle = load_bundle(args.bundle)
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
        bundle = load_bundle(args.bundle)
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
