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
    benchmark.add_argument("--candidate-edges", type=Path, default=None,
                           help="TSV of sender_state, receiver_state, ligand, receptor: score exactly these edges in every "
                                "split with the expression thresholds and LR cap off (see candidate-edges).")

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

    pb = subparsers.add_parser(
        "sample-pseudobulk",
        help="Per-sample scores (benchmark-h5ad split_scores_long.tsv) as a sample x receiver-state h5ad for cellHarmony_differential.",
    )
    pb.add_argument("--split-scores", type=Path, required=True, help="split_scores_long.tsv from benchmark-h5ad.")
    pb.add_argument("--output", type=Path, required=True, help="Output h5ad.")
    pb.add_argument("--score", default="fastcomm_score", help="Score column used as the value.")
    pb.add_argument("--sample-col", default="Sample", help="obs column for the sample (split) id.")
    pb.add_argument("--population-col", default="cell_type", help="obs column for the receiver state.")
    pb.add_argument("--covariates", type=Path, help="Optional TSV with --sample-col and --covariate-col columns.")
    pb.add_argument("--covariate-col", default="Condition", help="Covariate column in --covariates.")
    pb.add_argument("--covariates-output", type=Path, help="Covariates file for cellHarmony_differential (default beside --output).")
    pb.add_argument("--summary-json", type=Path, help="Summary JSON (default beside --output).")

    ce = subparsers.add_parser(
        "candidate-edges",
        help="Union of per-split edges from a thresholded benchmark-h5ad run, senders present in every split, for rescoring.",
    )
    ce.add_argument("--split-scores", type=Path, required=True, help="split_scores_long.tsv from benchmark-h5ad.")
    ce.add_argument("--h5ad", type=Path, required=True, help="The h5ad the split scores came from.")
    ce.add_argument("--state-key", required=True, help="obs state column used for scoring.")
    ce.add_argument("--split-key", required=True, help="obs split (sample) column used for scoring.")
    ce.add_argument("--min-cells", type=int, required=True, help="min cells per state used for scoring.")
    ce.add_argument("--output", type=Path, required=True, help="Output candidate edge TSV.")

    ov = subparsers.add_parser(
        "isoform-overlay",
        help="Long-read splicing and isoform changes of the ligand, receptor and downstream genes of differential edges.",
    )
    ov.add_argument("--deg-detailed", type=Path, required=True, help="cellHarmony_differential DEG_detailed TSV of the fastComm differential.")
    ov.add_argument("--candidate-edges", type=Path, required=True, help="candidate_edges.tsv of the rescored run (background).")
    ov.add_argument("--event-dir", type=Path, required=True, help="long-read dPSI-covariate-events folder.")
    ov.add_argument("--isoform-dir", type=Path, required=True, help="long-read diff-covariate-isoform folder.")
    ov.add_argument("--case", required=True); ov.add_argument("--control", required=True)
    ov.add_argument("--lr-table", type=Path, default=Path(__file__).parent / "resources" / "cellchat_nichenet_human" / "ligand_receptor.tsv")
    ov.add_argument("--response-matrix", type=Path,
                    default=Path(__file__).parent / "training_data" / "processed" / "upstream" / "cellchat_nichenet_human" / "response_signatures.tsv",
                    help="NicheNet response signatures; downstream genes are the top --top-downstream targets.")
    ov.add_argument("--switches", type=Path, default=None, help="isoform_state_switches.tsv (optional).")
    ov.add_argument("--alpha", type=float, default=0.1, help="edge significance: FDR (or raw P with --use-rawp) < alpha.")
    ov.add_argument("--use-rawp", action="store_true")
    ov.add_argument("--rawp", type=float, default=0.05, help="long-read evidence: ebayesPval < this.")
    ov.add_argument("--dpsi", type=float, default=0.1, help="splicing events: |dPSI| > this.")
    ov.add_argument("--fold", type=float, default=1.5, help="isoforms: |log2 fold| > log2(fold).")
    ov.add_argument("--min-n", type=int, default=3, help="donors with a value in each group.")
    ov.add_argument("--top-downstream", type=int, default=50, help="response genes per ligand or pathway signature.")
    ov.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "isoform-overlay":
        from .isoform_overlay import overlay

        overlay(str(args.deg_detailed), str(args.candidate_edges), str(args.lr_table), str(args.response_matrix),
                str(args.event_dir), str(args.isoform_dir), args.case, args.control, str(args.output_dir),
                alpha=args.alpha, use_rawp=args.use_rawp, rawp=args.rawp, dpsi=args.dpsi, fold=args.fold,
                min_n=args.min_n, top_downstream=args.top_downstream,
                switches=str(args.switches) if args.switches else None)
        return 0

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
                candidate_edges=args.candidate_edges,
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

    if args.command == "sample-pseudobulk":
        from .sample_pseudobulk import main_from_args

        return main_from_args(args)

    if args.command == "candidate-edges":
        from .sample_pseudobulk import build_candidate_edges

        _, summary = build_candidate_edges(args.split_scores, args.h5ad, args.output, args.state_key, args.split_key,
                                           args.min_cells)
        Path(args.output).with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
