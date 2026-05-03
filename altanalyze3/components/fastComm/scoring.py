from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse


_EPS = 1e-9


def clean_labels(values: Iterable[object]) -> List[str]:
    return [str(value).strip() for value in values]


def split_complex(value: object) -> Tuple[str, ...]:
    text = str(value).strip()
    if not text:
        return tuple()
    for delimiter in ("+", "|", ";", ","):
        text = text.replace(delimiter, "+")
    return tuple(part.strip() for part in text.split("+") if part.strip())


def geometric_mean(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return 0.0
    if np.any(array <= 0):
        return 0.0
    return float(np.exp(np.log(array).mean()))


@dataclass(frozen=True)
class PseudobulkState:
    expression: pd.DataFrame
    detection: pd.DataFrame
    sizes: pd.Series


def to_dense_frame(matrix, *, index: Sequence[str], columns: Sequence[str]) -> pd.DataFrame:
    if sparse.issparse(matrix):
        values = matrix.toarray()
    else:
        values = np.asarray(matrix)
    return pd.DataFrame(values, index=clean_labels(index), columns=clean_labels(columns))


def make_state_pseudobulk(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    state_key: str,
    min_cells: int = 1,
) -> PseudobulkState:
    if state_key not in metadata.columns:
        raise KeyError(f"State column {state_key!r} was not found in metadata")

    expression = expression.copy()
    expression.index = clean_labels(expression.index)
    expression.columns = clean_labels(expression.columns)
    metadata = metadata.copy()
    metadata.index = clean_labels(metadata.index)

    common = expression.index.intersection(metadata.index)
    if common.empty:
        raise ValueError("No shared cells between expression index and metadata index")

    expression = expression.loc[common]
    states = metadata.loc[common, state_key].astype(str).str.strip()
    valid = states.ne("") & states.notna()
    expression = expression.loc[valid]
    states = states.loc[valid]

    sizes = states.value_counts().sort_index()
    keep_states = sizes.index[sizes >= min_cells]
    if keep_states.empty:
        raise ValueError(f"No states passed min_cells={min_cells}")

    expression = expression.loc[states.isin(keep_states)]
    states = states.loc[states.isin(keep_states)]
    grouped = expression.groupby(states, sort=True)
    mean_expr = grouped.mean().astype(float)
    detection = expression.gt(0).groupby(states, sort=True).mean().astype(float)

    return PseudobulkState(
        expression=mean_expr,
        detection=detection,
        sizes=sizes.loc[keep_states].sort_index(),
    )


def _gene_value(state_by_gene: pd.DataFrame, state: str, gene: str) -> float:
    if gene not in state_by_gene.columns or state not in state_by_gene.index:
        return 0.0
    return float(state_by_gene.at[state, gene])


def _complex_value(state_by_gene: pd.DataFrame, state: str, members: Sequence[str]) -> float:
    values = [_gene_value(state_by_gene, state, member) for member in members]
    return geometric_mean(values)


def _complex_detection(detection: pd.DataFrame, state: str, members: Sequence[str]) -> float:
    values = [_gene_value(detection, state, member) for member in members]
    if not values:
        return 0.0
    return float(min(values))


def _specificity(state_by_gene: pd.DataFrame, state: str, members: Sequence[str]) -> float:
    if not members:
        return 0.0
    numerator = _complex_value(state_by_gene, state, members)
    totals = []
    for member in members:
        if member in state_by_gene.columns:
            totals.append(float(state_by_gene[member].sum()))
        else:
            totals.append(0.0)
    denominator = geometric_mean(totals)
    return float(numerator / (denominator + _EPS))


def prepare_ligand_receptor_table(lr_table: pd.DataFrame) -> pd.DataFrame:
    required = {"ligand", "receptor"}
    missing = required.difference(lr_table.columns)
    if missing:
        raise ValueError(f"LR table is missing required columns: {sorted(missing)}")

    out = lr_table.copy()
    out["ligand"] = clean_labels(out["ligand"])
    out["receptor"] = clean_labels(out["receptor"])
    if "pathway" not in out.columns:
        out["pathway"] = ""
    if "evidence_weight" not in out.columns:
        out["evidence_weight"] = 1.0
    if "interaction_class" not in out.columns:
        out["interaction_class"] = "unknown"
    out["evidence_weight"] = pd.to_numeric(out["evidence_weight"], errors="coerce").fillna(1.0)
    out = out.loc[out["ligand"].ne("") & out["receptor"].ne("")].copy()
    out["lr_key"] = out["ligand"] + "->" + out["receptor"]
    return out.reset_index(drop=True)


def score_ligand_receptor_expression(
    state: PseudobulkState,
    lr_table: pd.DataFrame,
    *,
    min_ligand_expr: float = 0.0,
    min_receptor_expr: float = 0.0,
    include_self_edges: bool = True,
) -> pd.DataFrame:
    lr = prepare_ligand_receptor_table(lr_table)
    rows: List[Dict[str, object]] = []
    states = list(state.expression.index)

    for edge in lr.itertuples(index=False):
        ligands = split_complex(edge.ligand)
        receptors = split_complex(edge.receptor)
        for sender in states:
            ligand_expr = _complex_value(state.expression, sender, ligands)
            if ligand_expr < min_ligand_expr:
                continue
            ligand_detection = _complex_detection(state.detection, sender, ligands)
            ligand_specificity = _specificity(state.expression, sender, ligands)

            for receiver in states:
                if not include_self_edges and sender == receiver:
                    continue
                receptor_expr = _complex_value(state.expression, receiver, receptors)
                if receptor_expr < min_receptor_expr:
                    continue
                receptor_detection = _complex_detection(state.detection, receiver, receptors)
                receptor_specificity = _specificity(state.expression, receiver, receptors)
                complex_completeness = min(ligand_detection, receptor_detection)
                expression_score = np.log1p(ligand_expr) * np.log1p(receptor_expr)
                specificity_score = np.sqrt(max(ligand_specificity, 0.0) * max(receptor_specificity, 0.0))
                lr_expression_score = (
                    expression_score
                    * (0.5 + 0.5 * complex_completeness)
                    * (0.5 + 0.5 * specificity_score)
                    * float(edge.evidence_weight)
                )
                rows.append(
                    {
                        "sender_state": sender,
                        "receiver_state": receiver,
                        "ligand": edge.ligand,
                        "receptor": edge.receptor,
                        "lr_key": edge.lr_key,
                        "pathway": edge.pathway,
                        "interaction_class": edge.interaction_class,
                        "ligand_expr": ligand_expr,
                        "receptor_expr": receptor_expr,
                        "ligand_detection": ligand_detection,
                        "receptor_detection": receptor_detection,
                        "complex_completeness": complex_completeness,
                        "ligand_specificity": ligand_specificity,
                        "receptor_specificity": receptor_specificity,
                        "evidence_weight": float(edge.evidence_weight),
                        "lr_expression_score": float(lr_expression_score),
                    }
                )

    if not rows:
        return pd.DataFrame(
            columns=[
                "sender_state",
                "receiver_state",
                "ligand",
                "receptor",
                "lr_key",
                "pathway",
                "interaction_class",
                "ligand_expr",
                "receptor_expr",
                "ligand_detection",
                "receptor_detection",
                "complex_completeness",
                "ligand_specificity",
                "receptor_specificity",
                "evidence_weight",
                "lr_expression_score",
            ]
        )
    return pd.DataFrame(rows).sort_values("lr_expression_score", ascending=False).reset_index(drop=True)


def load_response_matrix(path: str) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep=None, engine="python", index_col=0)
    matrix.index = clean_labels(matrix.index)
    matrix.columns = clean_labels(matrix.columns)
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def make_receiver_delta(
    state_expression: pd.DataFrame,
    *,
    baseline_state: Optional[str] = None,
    baseline_map: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    if baseline_map:
        rows = {}
        for state_name, base_name in baseline_map.items():
            if state_name in state_expression.index and base_name in state_expression.index:
                rows[state_name] = state_expression.loc[state_name] - state_expression.loc[base_name]
        if not rows:
            raise ValueError("None of the baseline_map state pairs were found in the expression matrix")
        return pd.DataFrame.from_dict(rows, orient="index").astype(float)

    if baseline_state:
        if baseline_state not in state_expression.index:
            raise KeyError(f"Baseline state {baseline_state!r} was not found")
        baseline = state_expression.loc[baseline_state]
    else:
        baseline = state_expression.median(axis=0)
    return state_expression.subtract(baseline, axis=1).astype(float)


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= _EPS:
        return 0.0
    return float(np.dot(a, b) / denom)


def _empirical_upper_tail(values: pd.Series, *, group: Optional[pd.Series] = None) -> pd.Series:
    values = values.astype(float)
    if values.empty:
        return pd.Series(dtype=float, index=values.index)
    if group is None:
        ranks = values.rank(method="min", ascending=False)
        return ranks / float(values.shape[0])

    out = pd.Series(index=values.index, dtype=float)
    for _, idx in group.groupby(group, sort=False).groups.items():
        subset = values.loc[idx]
        ranks = subset.rank(method="min", ascending=False)
        out.loc[idx] = ranks / float(subset.shape[0])
    return out



def add_receiver_response_scores(
    edges: pd.DataFrame,
    receiver_delta: pd.DataFrame,
    response_matrix: Optional[pd.DataFrame],
    *,
    response_key_order: Sequence[str] = ("lr_key", "ligand", "pathway"),
    max_support_genes: int = 8,
) -> pd.DataFrame:
    out = edges.copy()
    if response_matrix is None or response_matrix.empty or out.empty:
        out["receiver_response_score"] = 0.0
        out["state_promotion_score"] = 0.0
        out["response_key"] = ""
        out["response_genes_used"] = 0
        out["response_support_genes"] = ""
        return out

    response_matrix = response_matrix.copy()
    response_matrix.index = clean_labels(response_matrix.index)
    response_matrix.columns = clean_labels(response_matrix.columns)
    common_genes = receiver_delta.columns.intersection(response_matrix.columns)
    if common_genes.empty:
        out["receiver_response_score"] = 0.0
        out["state_promotion_score"] = 0.0
        out["response_key"] = ""
        out["response_genes_used"] = 0
        out["response_support_genes"] = ""
        return out

    response_sub = response_matrix.loc[:, common_genes]
    delta_sub = receiver_delta.loc[:, common_genes]
    response_scores: List[float] = []
    promotion_scores: List[float] = []
    response_keys: List[str] = []
    genes_used: List[int] = []
    support_genes: List[str] = []

    for edge in out.itertuples(index=False):
        key = ""
        expected = None
        for field in response_key_order:
            candidate = str(getattr(edge, field, "")).strip()
            if candidate and candidate in response_sub.index:
                key = candidate
                expected = response_sub.loc[candidate].to_numpy(dtype=float)
                break
        if expected is None or edge.receiver_state not in delta_sub.index:
            response_scores.append(0.0)
            promotion_scores.append(0.0)
            response_keys.append("")
            genes_used.append(0)
            support_genes.append("")
            continue

        observed = delta_sub.loc[edge.receiver_state].to_numpy(dtype=float)
        cosine = _cosine(observed, expected)
        contribution = observed * expected
        if max_support_genes > 0:
            nonzero = np.where((expected != 0) & np.isfinite(contribution))[0]
            ranked = sorted(nonzero, key=lambda idx: contribution[idx], reverse=True)
            support = [
                f"{common_genes[idx]}:{contribution[idx]:.3g}"
                for idx in ranked[:max_support_genes]
                if contribution[idx] > 0
            ]
        else:
            support = []
        response_scores.append(max(cosine, 0.0))
        promotion_scores.append(cosine)
        response_keys.append(key)
        genes_used.append(int(np.count_nonzero(expected)))
        support_genes.append(";".join(support))

    out["receiver_response_score"] = response_scores
    out["state_promotion_score"] = promotion_scores
    out["response_key"] = response_keys
    out["response_genes_used"] = genes_used
    out["response_support_genes"] = support_genes
    return out


def finalize_scores(
    edges: pd.DataFrame,
    *,
    lr_weight: float = 0.55,
    response_weight: float = 0.30,
    promotion_weight: float = 0.15,
) -> pd.DataFrame:
    if edges.empty:
        out = edges.copy()
        out["fastcomm_score"] = []
        out["score_rank"] = []
        return out

    out = edges.copy()
    lr = out["lr_expression_score"].astype(float)
    lr_scaled = lr / (float(lr.max()) + _EPS)
    response = out.get("receiver_response_score", pd.Series(0.0, index=out.index)).astype(float).clip(lower=0.0)
    promotion = out.get("state_promotion_score", pd.Series(0.0, index=out.index)).astype(float).clip(lower=0.0)
    raw = lr_weight * lr_scaled + response_weight * response + promotion_weight * promotion
    out["lr_expression_score_scaled"] = lr_scaled
    out["fastcomm_score"] = raw.clip(lower=0.0, upper=1.0)
    out["lr_empirical_p"] = _empirical_upper_tail(out["lr_expression_score"], group=out["lr_key"])
    out["fastcomm_empirical_p"] = _empirical_upper_tail(out["fastcomm_score"])
    out["lr_within_lr_percentile"] = (1.0 - out["lr_empirical_p"]).clip(lower=0.0, upper=1.0)
    out["fastcomm_percentile"] = (1.0 - out["fastcomm_empirical_p"]).clip(lower=0.0, upper=1.0)
    out["score_rank"] = out["fastcomm_score"].rank(method="first", ascending=False).astype(int)
    return out.sort_values("fastcomm_score", ascending=False).reset_index(drop=True)
