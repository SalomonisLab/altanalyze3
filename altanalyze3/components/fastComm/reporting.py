from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


@dataclass(frozen=True)
class ExemplarReportParams:
    scores: Path
    output_tsv: Path
    output_md: Optional[Path] = None
    split_stability: Optional[Path] = None
    top_n: int = 25
    min_score: float = 0.0
    one_per_state_pair: bool = True
    title: str = "fastComm Exemplar Interactions"


def confidence_label(score: float) -> str:
    if score >= 0.65:
        return "high"
    if score >= 0.40:
        return "moderate"
    if score >= 0.25:
        return "suggestive"
    return "low"


def _load_scores(path: Path) -> pd.DataFrame:
    scores = pd.read_csv(path, sep="\t")
    required = {
        "sender_state",
        "receiver_state",
        "ligand",
        "receptor",
        "fastcomm_score",
        "lr_expression_score_scaled",
        "receiver_response_score",
        "state_promotion_score",
        "response_key",
        "response_support_genes",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(f"Scores file is missing columns: {sorted(missing)}")
    return scores


def _select_exemplars(params: ExemplarReportParams) -> pd.DataFrame:
    scores = _load_scores(params.scores)
    scores = scores.loc[scores["fastcomm_score"].astype(float) >= params.min_score].copy()
    scores = scores.sort_values("fastcomm_score", ascending=False)
    if params.one_per_state_pair:
        scores = scores.drop_duplicates(["sender_state", "receiver_state"], keep="first")
    return scores.head(params.top_n).copy()


def _format_evidence(row: pd.Series) -> str:
    pieces = [
        f"ligand={row['ligand_expr']:.3g}",
        f"receptor={row['receptor_expr']:.3g}",
        f"complex={row['complex_completeness']:.2f}",
        f"LR_pct={row.get('lr_within_lr_percentile', 0.0):.2f}",
        f"response={row['receiver_response_score']:.2f}",
    ]
    support = str(row.get("response_support_genes", "")).strip()
    if support:
        pieces.append(f"support={support}")
    return "; ".join(pieces)


def build_exemplar_table(params: ExemplarReportParams) -> pd.DataFrame:
    exemplars = _select_exemplars(params)
    rows: List[Dict[str, object]] = []
    for row in exemplars.itertuples(index=False):
        series = pd.Series(row._asdict())
        probability = float(series["fastcomm_score"])
        rows.append(
            {
                "rank": int(series.get("score_rank", len(rows) + 1)),
                "sender_state": series["sender_state"],
                "receiver_state": series["receiver_state"],
                "ligand": series["ligand"],
                "receptor": series["receptor"],
                "pathway_or_factor": series["response_key"] or series.get("pathway", ""),
                "interaction_class": series.get("interaction_class", ""),
                "prototype_probability": probability,
                "confidence": confidence_label(probability),
                "fastcomm_percentile": float(series.get("fastcomm_percentile", 0.0)),
                "lr_expression_score": float(series["lr_expression_score_scaled"]),
                "receiver_response_score": float(series["receiver_response_score"]),
                "state_promotion_score": float(series["state_promotion_score"]),
                "ligand_expr_sender": float(series["ligand_expr"]),
                "receptor_expr_receiver": float(series["receptor_expr"]),
                "complex_completeness": float(series["complex_completeness"]),
                "ligand_detection": float(series["ligand_detection"]),
                "receptor_detection": float(series["receptor_detection"]),
                "supporting_response_genes": series.get("response_support_genes", ""),
                "evidence_summary": _format_evidence(series),
            }
        )
    return pd.DataFrame(rows)


def _markdown_table(df: pd.DataFrame) -> str:
    display_cols = [
        "rank",
        "sender_state",
        "receiver_state",
        "ligand",
        "receptor",
        "pathway_or_factor",
        "prototype_probability",
        "confidence",
        "lr_expression_score",
        "receiver_response_score",
        "supporting_response_genes",
    ]
    out = df.loc[:, display_cols].copy()
    for col in ["prototype_probability", "lr_expression_score", "receiver_response_score"]:
        out[col] = out[col].map(lambda value: f"{float(value):.3f}")
    headers = [str(column) for column in out.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in out.itertuples(index=False):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


def _stability_text(path: Optional[Path]) -> str:
    if path is None or not path.exists():
        return ""
    stability = pd.read_csv(path, sep="\t")
    if stability.empty:
        return ""
    lines = ["## Split Stability", ""]
    headers = [str(column) for column in stability.columns]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in stability.itertuples(index=False):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    lines.append("")
    return "\n".join(lines)


def write_exemplar_report(params: ExemplarReportParams) -> Dict[str, object]:
    table = build_exemplar_table(params)
    params.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(params.output_tsv, sep="\t", index=False)

    if params.output_md is not None:
        params.output_md.parent.mkdir(parents=True, exist_ok=True)
        body = [
            f"# {params.title}",
            "",
            "Scores are prototype probability-like confidence scores on a 0-1 scale. They are not yet benchmark-calibrated probabilities.",
            "",
            "The `pathway_or_factor` and supporting genes come from the active receiver-response signature matrix.",
            "",
            "## Top Exemplar Interactions",
            "",
            _markdown_table(table) if not table.empty else "No exemplars passed the filters.",
            "",
            _stability_text(params.split_stability),
        ]
        params.output_md.write_text("\n".join(part for part in body if part is not None), encoding="utf-8")

    manifest = {
        "scores": str(params.scores),
        "output_tsv": str(params.output_tsv),
        "output_md": str(params.output_md) if params.output_md else None,
        "split_stability": str(params.split_stability) if params.split_stability else None,
        "top_n": int(params.top_n),
        "min_score": float(params.min_score),
        "one_per_state_pair": bool(params.one_per_state_pair),
        "n_exemplars": int(table.shape[0]),
    }
    manifest_path = params.output_tsv.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
