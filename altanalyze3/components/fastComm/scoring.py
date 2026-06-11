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


def make_state_pseudobulk_from_matrix(
    matrix,
    *,
    index: Sequence[str],
    columns: Sequence[str],
    metadata: pd.DataFrame,
    state_key: str,
    min_cells: int = 1,
) -> PseudobulkState:
    if state_key not in metadata.columns:
        raise KeyError(f"State column {state_key!r} was not found in metadata")

    obs_index = pd.Index(clean_labels(index))
    metadata = metadata.copy()
    metadata.index = clean_labels(metadata.index)
    common = obs_index.intersection(metadata.index)
    if common.empty:
        raise ValueError("No shared cells between expression index and metadata index")

    row_positions = obs_index.get_indexer(common)
    states = metadata.loc[common, state_key].astype(str).str.strip()
    valid = states.ne("") & states.notna()
    row_positions = row_positions[valid.to_numpy()]
    states = states.loc[valid]

    sizes = states.value_counts().sort_index()
    keep_states = sizes.index[sizes >= min_cells]
    if keep_states.empty:
        raise ValueError(f"No states passed min_cells={min_cells}")

    state_names = list(keep_states.sort_values())
    matrix_subset = matrix[row_positions]
    if sparse.issparse(matrix_subset):
        matrix_subset = matrix_subset.tocsr()

    mean_rows: List[np.ndarray] = []
    detection_rows: List[np.ndarray] = []
    state_values = states.to_numpy(dtype=object)
    for state_name in state_names:
        mask = state_values == state_name
        state_matrix = matrix_subset[mask]
        if sparse.issparse(state_matrix):
            mean_rows.append(np.asarray(state_matrix.mean(axis=0)).ravel())
            detection_rows.append(np.asarray((state_matrix > 0).mean(axis=0)).ravel())
        else:
            mean_rows.append(np.asarray(state_matrix, dtype=float).mean(axis=0))
            detection_rows.append((np.asarray(state_matrix) > 0).mean(axis=0))

    column_labels = clean_labels(columns)
    return PseudobulkState(
        expression=pd.DataFrame(mean_rows, index=state_names, columns=column_labels).astype(float),
        detection=pd.DataFrame(detection_rows, index=state_names, columns=column_labels).astype(float),
        sizes=sizes.loc[state_names],
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


def _complex_arrays(
    state: PseudobulkState,
    complexes: Iterable[Tuple[str, ...]],
) -> Dict[Tuple[str, ...], Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    expression = state.expression.to_numpy(dtype=float, copy=False)
    detection = state.detection.to_numpy(dtype=float, copy=False)
    gene_to_idx = {gene: idx for idx, gene in enumerate(clean_labels(state.expression.columns))}
    gene_totals = expression.sum(axis=0)
    cache: Dict[Tuple[str, ...], Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    n_states = expression.shape[0]

    for members in complexes:
        if not members:
            zeros = np.zeros(n_states, dtype=float)
            cache[members] = (zeros, zeros, zeros)
            continue

        idxs = [gene_to_idx.get(member) for member in members]
        if any(idx is None for idx in idxs):
            zeros = np.zeros(n_states, dtype=float)
            cache[members] = (zeros, zeros, zeros)
            continue

        idx_array = np.asarray(idxs, dtype=int)
        if idx_array.size == 1:
            expr_values = expression[:, idx_array[0]].astype(float, copy=True)
            detection_values = detection[:, idx_array[0]].astype(float, copy=True)
            denominator = float(gene_totals[idx_array[0]])
        else:
            member_expr = expression[:, idx_array]
            positive = np.all(member_expr > 0, axis=1)
            expr_values = np.zeros(n_states, dtype=float)
            if np.any(positive):
                expr_values[positive] = np.exp(np.log(member_expr[positive]).mean(axis=1))
            detection_values = detection[:, idx_array].min(axis=1).astype(float, copy=False)
            totals = gene_totals[idx_array]
            denominator = geometric_mean(totals)
        specificity_values = expr_values / (denominator + _EPS)
        cache[members] = (expr_values, detection_values, specificity_values)

    return cache


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
    min_lr_expression_score: float = 0.0,
    include_self_edges: bool = True,
) -> pd.DataFrame:
    lr = prepare_ligand_receptor_table(lr_table)
    output_columns = [
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
    if lr.empty:
        return pd.DataFrame(
            columns=output_columns
        )

    ligand_members = [split_complex(value) for value in lr["ligand"]]
    receptor_members = [split_complex(value) for value in lr["receptor"]]
    complex_cache = _complex_arrays(state, set(ligand_members) | set(receptor_members))
    states = np.asarray(clean_labels(state.expression.index), dtype=object)
    chunks: Dict[str, List[np.ndarray]] = {column: [] for column in output_columns}

    ligand_values = lr["ligand"].to_numpy(dtype=object)
    receptor_values = lr["receptor"].to_numpy(dtype=object)
    lr_key_values = lr["lr_key"].to_numpy(dtype=object)
    pathway_values = lr["pathway"].to_numpy(dtype=object)
    class_values = lr["interaction_class"].to_numpy(dtype=object)
    evidence_values = lr["evidence_weight"].to_numpy(dtype=float)

    for idx in range(lr.shape[0]):
        ligand_expr, ligand_detection, ligand_specificity = complex_cache[ligand_members[idx]]
        receptor_expr, receptor_detection, receptor_specificity = complex_cache[receptor_members[idx]]
        senders = np.flatnonzero(ligand_expr >= min_ligand_expr)
        receivers = np.flatnonzero(receptor_expr >= min_receptor_expr)
        if senders.size == 0 or receivers.size == 0:
            continue

        sender_idx = np.repeat(senders, receivers.size)
        receiver_idx = np.tile(receivers, senders.size)
        if not include_self_edges:
            keep = sender_idx != receiver_idx
            if not np.any(keep):
                continue
            sender_idx = sender_idx[keep]
            receiver_idx = receiver_idx[keep]

        ligand_expr_values = ligand_expr[sender_idx]
        receptor_expr_values = receptor_expr[receiver_idx]
        ligand_detection_values = ligand_detection[sender_idx]
        receptor_detection_values = receptor_detection[receiver_idx]
        ligand_specificity_values = ligand_specificity[sender_idx]
        receptor_specificity_values = receptor_specificity[receiver_idx]
        complex_completeness = np.minimum(ligand_detection_values, receptor_detection_values)
        expression_score = np.log1p(ligand_expr_values) * np.log1p(receptor_expr_values)
        specificity_score = np.sqrt(
            np.maximum(ligand_specificity_values, 0.0)
            * np.maximum(receptor_specificity_values, 0.0)
        )
        lr_expression_score = (
            expression_score
            * (0.5 + 0.5 * complex_completeness)
            * (0.5 + 0.5 * specificity_score)
            * evidence_values[idx]
        )
        if min_lr_expression_score > 0:
            keep = lr_expression_score >= min_lr_expression_score
            if not np.any(keep):
                continue
            sender_idx = sender_idx[keep]
            receiver_idx = receiver_idx[keep]
            ligand_expr_values = ligand_expr_values[keep]
            receptor_expr_values = receptor_expr_values[keep]
            ligand_detection_values = ligand_detection_values[keep]
            receptor_detection_values = receptor_detection_values[keep]
            ligand_specificity_values = ligand_specificity_values[keep]
            receptor_specificity_values = receptor_specificity_values[keep]
            complex_completeness = complex_completeness[keep]
            lr_expression_score = lr_expression_score[keep]

        n_rows = sender_idx.size
        chunks["sender_state"].append(states[sender_idx])
        chunks["receiver_state"].append(states[receiver_idx])
        chunks["ligand"].append(np.repeat(ligand_values[idx], n_rows))
        chunks["receptor"].append(np.repeat(receptor_values[idx], n_rows))
        chunks["lr_key"].append(np.repeat(lr_key_values[idx], n_rows))
        chunks["pathway"].append(np.repeat(pathway_values[idx], n_rows))
        chunks["interaction_class"].append(np.repeat(class_values[idx], n_rows))
        chunks["ligand_expr"].append(ligand_expr_values)
        chunks["receptor_expr"].append(receptor_expr_values)
        chunks["ligand_detection"].append(ligand_detection_values)
        chunks["receptor_detection"].append(receptor_detection_values)
        chunks["complex_completeness"].append(complex_completeness)
        chunks["ligand_specificity"].append(ligand_specificity_values)
        chunks["receptor_specificity"].append(receptor_specificity_values)
        chunks["evidence_weight"].append(np.repeat(evidence_values[idx], n_rows))
        chunks["lr_expression_score"].append(lr_expression_score)

    if not chunks["lr_expression_score"]:
        return pd.DataFrame(columns=output_columns)

    return pd.DataFrame(
        {column: np.concatenate(values) for column, values in chunks.items()},
        columns=output_columns,
    ).sort_values("lr_expression_score", ascending=False).reset_index(drop=True)


def limit_lr_candidates_per_state_pair(edges: pd.DataFrame, max_candidates: Optional[int]) -> pd.DataFrame:
    if edges.empty or max_candidates is None or max_candidates <= 0:
        return edges
    return (
        edges.sort_values("lr_expression_score", ascending=False)
        .groupby(["sender_state", "receiver_state"], sort=False, group_keys=False)
        .head(int(max_candidates))
        .reset_index(drop=True)
    )


def load_response_matrix(path: str) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
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
    response_index = pd.Index(response_sub.index.astype(str))
    response_key_values = np.full(out.shape[0], "", dtype=object)
    unresolved = np.ones(out.shape[0], dtype=bool)
    for field in response_key_order:
        if field not in out.columns or not np.any(unresolved):
            continue
        candidates = out[field].astype(str).str.strip()
        matched = unresolved & candidates.isin(response_index).to_numpy()
        if np.any(matched):
            response_key_values[matched] = candidates.loc[matched].to_numpy(dtype=object)
            unresolved[matched] = False

    receiver_index = pd.Index(delta_sub.index.astype(str))
    receiver_values = out["receiver_state"].astype(str).to_numpy(dtype=object)
    receiver_pos = receiver_index.get_indexer(receiver_values)
    response_pos = response_index.get_indexer(response_key_values)
    valid = (receiver_pos >= 0) & (response_pos >= 0)

    response_scores = np.zeros(out.shape[0], dtype=float)
    promotion_scores = np.zeros(out.shape[0], dtype=float)
    genes_used = np.zeros(out.shape[0], dtype=int)
    support_genes = np.full(out.shape[0], "", dtype=object)

    if np.any(valid):
        delta_values = delta_sub.to_numpy(dtype=float, copy=False)
        response_values = response_sub.to_numpy(dtype=float, copy=False)
        delta_norm = np.linalg.norm(delta_values, axis=1)
        response_norm = np.linalg.norm(response_values, axis=1)
        denom = delta_norm[:, None] * response_norm[None, :]
        cosine = np.divide(
            delta_values @ response_values.T,
            denom,
            out=np.zeros((delta_values.shape[0], response_values.shape[0]), dtype=float),
            where=denom > _EPS,
        )
        valid_idx = np.flatnonzero(valid)
        valid_cosine = cosine[receiver_pos[valid_idx], response_pos[valid_idx]]
        promotion_scores[valid_idx] = valid_cosine
        response_scores[valid_idx] = np.maximum(valid_cosine, 0.0)
        response_gene_counts = np.count_nonzero(response_values, axis=1)
        genes_used[valid_idx] = response_gene_counts[response_pos[valid_idx]]

        unique_pairs = pd.DataFrame(
            {
                "receiver_pos": receiver_pos[valid_idx],
                "response_pos": response_pos[valid_idx],
            }
        ).drop_duplicates()
        # Support strings are diagnostic only. Keep them for moderate workloads,
        # but avoid an O(unique pairs * genes) Python pass on full web-scale runs.
        if max_support_genes > 0 and unique_pairs.shape[0] * len(common_genes) <= 50_000_000:
            support_by_pair: Dict[Tuple[int, int], str] = {}
            common_gene_values = np.asarray(common_genes, dtype=object)
            for pair in unique_pairs.itertuples(index=False):
                contribution = delta_values[pair.receiver_pos] * response_values[pair.response_pos]
                eligible = np.flatnonzero((response_values[pair.response_pos] != 0) & np.isfinite(contribution) & (contribution > 0))
                if eligible.size == 0:
                    support_by_pair[(pair.receiver_pos, pair.response_pos)] = ""
                    continue
                if eligible.size > max_support_genes:
                    top_local = np.argpartition(-contribution[eligible], max_support_genes - 1)[:max_support_genes]
                    ranked = eligible[top_local[np.argsort(-contribution[eligible][top_local])]]
                else:
                    ranked = eligible[np.argsort(-contribution[eligible])]
                support_by_pair[(pair.receiver_pos, pair.response_pos)] = ";".join(
                    f"{common_gene_values[idx]}:{contribution[idx]:.3g}" for idx in ranked
                )
            support_genes[valid_idx] = [
                support_by_pair.get((int(receiver), int(response)), "")
                for receiver, response in zip(receiver_pos[valid_idx], response_pos[valid_idx])
            ]

    out["receiver_response_score"] = response_scores
    out["state_promotion_score"] = promotion_scores
    out["response_key"] = response_key_values
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
