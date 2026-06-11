from __future__ import annotations

import argparse
import json
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_LR_TABLE = PACKAGE_DIR / "resources" / "__auto__ligand_receptor.tsv"
DEFAULT_RESPONSE_MATRIX = PACKAGE_DIR / "resources" / "__auto__response_signatures.tsv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fast receptor-ligand communication scoring with optional receiver-response evidence."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    score = subparsers.add_parser("score", help="Score state-state ligand-receptor interactions.")
    input_group = score.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--h5ad", type=Path, help="Input AnnData file.")
    input_group.add_argument("--expression", type=Path, help="Cell-by-gene CSV/TSV expression matrix.")
    score.add_argument("--metadata", type=Path, help="Cell metadata CSV/TSV. Required with --expression.")
    score.add_argument("--output", type=Path, required=True, help="Output TSV for scored interactions.")
    score.add_argument(
        "--lr-table",
        type=Path,
        default=DEFAULT_LR_TABLE,
        help="Ligand-receptor TSV/CSV. Defaults to the built upstream bundle for the inferred/provided species.",
    )
    score.add_argument(
        "--response-matrix",
        type=Path,
        default=DEFAULT_RESPONSE_MATRIX,
        help="Response matrix, rows are LR/ligand/pathway keys. Defaults to the built upstream bundle for the inferred/provided species when available.",
    )
    score.add_argument("--no-response", action="store_true", help="Disable receiver-response scoring.")
    score.add_argument("--lr-source", action="append", help="Restrict LR rows to this source value, e.g. CellChatDB or NicheNet. Repeat to allow multiple sources.")
    score.add_argument("--state-pair-output", type=Path, help="Optional state-pair summary TSV.")
    score.add_argument("--state-expression-output", type=Path, help="Optional state pseudobulk expression TSV.")
    score.add_argument("--state-key", default="cell_state", help="Metadata/obs column defining cell states.")
    score.add_argument("--layer", help="AnnData layer to use instead of X.")
    score.add_argument("--gene-symbol-col", help="AnnData var column containing gene symbols.")
    score.add_argument("--species", choices=["human", "mouse"], help="Optional species hint used to resolve built upstream bundles.")
    score.add_argument("--min-cells", type=int, default=1, help="Minimum cells per state.")
    score.add_argument("--min-ligand-expr", type=float, default=0.01, help="Minimum sender ligand expression.")
    score.add_argument("--min-receptor-expr", type=float, default=0.01, help="Minimum receiver receptor expression.")
    score.add_argument("--min-lr-expression-score", type=float, default=0.001, help="Minimum raw LR expression score.")
    score.add_argument(
        "--max-lr-candidates-per-state-pair",
        type=int,
        default=25,
        help="Keep the top N raw LR candidates per sender/receiver state pair before response scoring. Use 0 for no cap.",
    )
    score.add_argument("--exclude-self-edges", action="store_true", help="Skip sender=receiver state edges.")
    score.add_argument("--baseline-state", help="Receiver baseline state used for response deltas.")
    score.add_argument("--top-n", type=int, help="Keep only the top N scored edges.")
    score.add_argument("--summary-json", type=Path, help="Optional JSON summary path.")

    registry = subparsers.add_parser("source-registry", help="Print the planned autonomous source registry.")
    registry.add_argument(
        "--path",
        type=Path,
        default=Path(__file__).resolve().parent / "configs" / "source_registry.json",
        help="Path to the source registry JSON.",
    )

    benchmark = subparsers.add_parser("benchmark-h5ad", help="Run donor/sample stability benchmark on one h5ad.")
    benchmark.add_argument("--h5ad", type=Path, required=True, help="Input AnnData file.")
    benchmark.add_argument("--output-dir", type=Path, required=True, help="Directory for benchmark outputs.")
    benchmark.add_argument("--state-key", default="cell_state", help="Metadata/obs column defining cell states.")
    benchmark.add_argument("--split-key", default="Donor", help="Metadata/obs column used for split stability.")
    benchmark.add_argument(
        "--lr-table",
        type=Path,
        default=DEFAULT_LR_TABLE,
        help="Ligand-receptor TSV/CSV. Defaults to the built upstream bundle for the inferred/provided species.",
    )
    benchmark.add_argument(
        "--response-matrix",
        type=Path,
        default=DEFAULT_RESPONSE_MATRIX,
        help="Response matrix used for receiver-response scoring. Defaults to the built upstream bundle for the inferred/provided species when available.",
    )
    benchmark.add_argument("--no-response", action="store_true", help="Disable receiver-response scoring.")
    benchmark.add_argument("--lr-source", action="append", help="Restrict LR rows to this source value, e.g. CellChatDB or NicheNet. Repeat to allow multiple sources.")
    benchmark.add_argument("--layer", help="AnnData layer to use instead of X.")
    benchmark.add_argument("--gene-symbol-col", help="AnnData var column containing gene symbols.")
    benchmark.add_argument("--species", choices=["human", "mouse"], help="Optional species hint used to resolve built upstream bundles.")
    benchmark.add_argument("--min-cells", type=int, default=20, help="Minimum cells per state.")
    benchmark.add_argument("--min-ligand-expr", type=float, default=0.01, help="Minimum sender ligand expression.")
    benchmark.add_argument("--min-receptor-expr", type=float, default=0.01, help="Minimum receiver receptor expression.")
    benchmark.add_argument("--min-lr-expression-score", type=float, default=0.001, help="Minimum raw LR expression score.")
    benchmark.add_argument(
        "--max-lr-candidates-per-state-pair",
        type=int,
        default=25,
        help="Keep the top N raw LR candidates per sender/receiver state pair before response scoring. Use 0 for no cap.",
    )
    benchmark.add_argument("--include-self-edges", action="store_true", help="Include sender=receiver state edges.")
    benchmark.add_argument("--top-n-stability", type=int, default=100, help="Top edge count for split-vs-full Jaccard.")

    build_response = subparsers.add_parser(
        "build-response-matrix",
        help="Build a scorer-ready response matrix from perturbation signatures.",
    )
    build_response.add_argument("--input", type=Path, required=True, help="Input long or wide response table.")
    build_response.add_argument("--output", type=Path, required=True, help="Output response matrix TSV.")
    build_response.add_argument("--manifest", type=Path, help="Optional JSON manifest path.")
    build_response.add_argument("--input-format", choices=["auto", "long", "wide"], default="auto")
    build_response.add_argument("--id-col", default="signature", help="Long-table signature column.")
    build_response.add_argument("--gene-col", default="gene", help="Long-table gene column.")
    build_response.add_argument("--score-col", default="score", help="Long-table response score column.")
    build_response.add_argument("--top-genes", type=int, default=200, help="Keep top absolute genes per signature.")
    build_response.add_argument("--min-abs-score", type=float, default=0.0, help="Zero scores below this absolute value.")
    build_response.add_argument("--no-l2-normalize", action="store_true", help="Do not L2-normalize each signature.")

    upstream = subparsers.add_parser(
        "build-upstream-resources",
        help="Build larger fastComm LR/response resources directly from CellChatDB and NicheNet v2.",
    )
    upstream.add_argument("--species", choices=["human", "mouse"], required=True, help="Species bundle to build.")
    upstream.add_argument("--output-dir", type=Path, help="Processed bundle output directory.")
    upstream.add_argument("--raw-dir", type=Path, help="Raw download cache directory.")
    upstream.add_argument("--top-targets-per-signature", type=int, default=250, help="Keep top ligand/pathway target genes per signature.")
    upstream.add_argument("--min-target-weight", type=float, default=0.001, help="Drop NicheNet target weights below this value.")
    upstream.add_argument("--skip-nichenet-lr-network", action="store_true", help="Use CellChatDB only for LR definitions.")
    upstream.add_argument("--overwrite", action="store_true", help="Re-download and rebuild even if cached files exist.")

    exemplar = subparsers.add_parser("exemplar-report", help="Summarize top scored interactions for inspection.")
    exemplar.add_argument("--scores", type=Path, required=True, help="fastComm scores TSV.")
    exemplar.add_argument("--output-tsv", type=Path, required=True, help="Output exemplar TSV.")
    exemplar.add_argument("--output-md", type=Path, help="Optional Markdown report.")
    exemplar.add_argument("--split-stability", type=Path, help="Optional split_stability.tsv from benchmark-h5ad.")
    exemplar.add_argument("--top-n", type=int, default=25, help="Number of exemplar rows to report.")
    exemplar.add_argument("--min-score", type=float, default=0.0, help="Minimum fastComm score.")
    exemplar.add_argument("--allow-multiple-per-state-pair", action="store_true", help="Do not collapse to the best LR per state pair.")
    exemplar.add_argument("--title", default="fastComm Exemplar Interactions", help="Markdown report title.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "source-registry":
        print(args.path.read_text(encoding="utf-8"))
        return 0

    if args.command == "score":
        from .api import FastCommParams, run_fastcomm

        result = run_fastcomm(
            FastCommParams(
                expression=args.expression,
                metadata=args.metadata,
                h5ad=args.h5ad,
                output=args.output,
                lr_table=args.lr_table,
                response_matrix=None if args.no_response else args.response_matrix,
                lr_sources=args.lr_source,
                state_pair_output=args.state_pair_output,
                state_expression_output=args.state_expression_output,
                state_key=args.state_key,
                layer=args.layer,
                gene_symbol_col=args.gene_symbol_col,
                species=args.species,
                min_cells=args.min_cells,
                min_ligand_expr=args.min_ligand_expr,
                min_receptor_expr=args.min_receptor_expr,
                min_lr_expression_score=args.min_lr_expression_score,
                max_lr_candidates_per_state_pair=args.max_lr_candidates_per_state_pair,
                include_self_edges=not args.exclude_self_edges,
                baseline_state=args.baseline_state,
                top_n=args.top_n,
            )
        )
        if args.summary_json:
            args.summary_json.parent.mkdir(parents=True, exist_ok=True)
            args.summary_json.write_text(json.dumps(result.summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result.summary, indent=2, sort_keys=True))
        return 0

    if args.command == "benchmark-h5ad":
        from .benchmark import FastCommBenchmarkParams, run_benchmark

        summary = run_benchmark(
            FastCommBenchmarkParams(
                h5ad=args.h5ad,
                output_dir=args.output_dir,
                state_key=args.state_key,
                split_key=args.split_key,
                lr_table=args.lr_table,
                response_matrix=None if args.no_response else args.response_matrix,
                lr_sources=args.lr_source,
                layer=args.layer,
                gene_symbol_col=args.gene_symbol_col,
                species=args.species,
                min_cells=args.min_cells,
                min_ligand_expr=args.min_ligand_expr,
                min_receptor_expr=args.min_receptor_expr,
                min_lr_expression_score=args.min_lr_expression_score,
                max_lr_candidates_per_state_pair=args.max_lr_candidates_per_state_pair,
                include_self_edges=args.include_self_edges,
                top_n_stability=args.top_n_stability,
            )
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if args.command == "build-response-matrix":
        from .training import ResponseTrainingParams, build_response_matrix

        manifest = build_response_matrix(
            ResponseTrainingParams(
                input=args.input,
                output=args.output,
                manifest=args.manifest,
                id_col=args.id_col,
                gene_col=args.gene_col,
                score_col=args.score_col,
                input_format=args.input_format,
                top_genes=args.top_genes,
                min_abs_score=args.min_abs_score,
                l2_normalize=not args.no_l2_normalize,
            )
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0

    if args.command == "build-upstream-resources":
        from .upstream_resources import UpstreamResourceBuildParams, build_upstream_resources

        manifest = build_upstream_resources(
            UpstreamResourceBuildParams(
                species=args.species,
                output_dir=args.output_dir,
                raw_dir=args.raw_dir,
                top_targets_per_signature=args.top_targets_per_signature,
                min_target_weight=args.min_target_weight,
                include_nichenet_lr_network=not args.skip_nichenet_lr_network,
                overwrite=args.overwrite,
            )
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0

    if args.command == "exemplar-report":
        from .reporting import ExemplarReportParams, write_exemplar_report

        manifest = write_exemplar_report(
            ExemplarReportParams(
                scores=args.scores,
                output_tsv=args.output_tsv,
                output_md=args.output_md,
                split_stability=args.split_stability,
                top_n=args.top_n,
                min_score=args.min_score,
                one_per_state_pair=not args.allow_multiple_per_state_pair,
                title=args.title,
            )
        )
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
