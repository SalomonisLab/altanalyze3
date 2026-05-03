"""Synthetic aberrant-expression test for ADT predictors.

Background: a model that's useful for surfacing disease biology must
faithfully convert cell-level RNA changes into protein predictions, even
when those changes contradict the cell's normal cell-state baseline.

This test explicitly probes that property:

  1. Pick a non-T-cell cluster (e.g. Pro-B-1) where CD3 RNA and protein are
     near-zero in the reference.
  2. Take held-out cells from that cluster, set their RNA partners for a
     chosen target ADT (e.g. CD3D / CD3E / CD3G for Hu.CD3) to a high
     value chosen from the high tail of the cell type that *does* normally
     express that ADT (e.g. T cells).
  3. Run predict() on the boosted cells.
  4. Compute the lift in predicted protein vs the unboosted baseline.

A passing model shows a large lift; a centroid-anchored model shows none
(its prediction is determined by the cell's cluster).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Hand-picked test cases of (target ADT, normal-source-cluster, target-host-cluster).
# These pairs are biologically reliable: the target host cluster does NOT
# normally express the target protein in a healthy bone marrow atlas.
DEFAULT_PROBES: List[Tuple[str, List[str], List[str], List[str]]] = [
    # (adt_name,            host_cluster_substrings,  source_cluster_substrings, partner_genes)
    ("Hu.CD3",              ["Pro-B", "Pre-B", "B-cell", "Mono", "MEP", "HSC"],
                            ["T-", "T cell", "CD4", "CD8", "Treg", "NK"],
                            ["CD3D", "CD3E", "CD3G", "CD247"]),
    ("Hu.CD19",             ["T-", "Mono", "MEP", "HSC", "MPP", "GMP"],
                            ["B-cell", "Pro-B", "Pre-B", "Plasma"],
                            ["CD19"]),
    ("Hu.CD14_M5E2",        ["T-", "B-cell", "Pro-B", "Pre-B", "MEP", "HSC"],
                            ["Mono", "Macro", "Classical-Mono"],
                            ["CD14"]),
    ("Hu.CD34",             ["T-", "B-cell", "Mono", "Plasma"],
                            ["HSC", "MPP", "MultiLin", "GMP", "MEP", "CLP"],
                            ["CD34"]),
    ("Hu.CD8",              ["Pro-B", "Pre-B", "Mono", "MEP", "HSC", "B-cell"],
                            ["CD8", "T-CD8"],
                            ["CD8A", "CD8B"]),
]


@dataclass
class AberrantTestReport:
    adt_name: str
    host_cluster: str
    source_cluster: str
    partner_genes: List[str]
    n_host_cells: int
    baseline_pred_mean: float          # what the model predicts on unmodified host cells
    boosted_pred_mean: float           # what it predicts after boosting partner-gene RNA
    source_truth_mean: float           # what the ADT actually is in source cells
    host_truth_mean: float             # what the ADT actually is in host cells (should be low)
    # A useful summary: ratio = (boosted - host_truth) / (source_truth - host_truth)
    # 0.0 = no lift (centroid-style failure), 1.0 = full lift to source-cluster level
    lift_fraction: float


def _match_clusters(all_clusters: Sequence[str], substrings: Sequence[str]) -> List[str]:
    keep: List[str] = []
    for cluster in pd.unique(all_clusters):
        s = str(cluster)
        if any(sub.lower() in s.lower() for sub in substrings):
            keep.append(s)
    return keep


def run_aberrant_test(
    *,
    model,
    adata: ad.AnnData,
    rna_genes: Sequence[str],          # the model's expected RNA columns
    adt_names: Sequence[str],          # the model's output ADTs
    cluster_col: str,
    probes: Optional[List[Tuple[str, List[str], List[str], List[str]]]] = None,
    max_host_cells: int = 1000,
    boost_quantile: float = 0.95,
    seed: int = 0,
) -> List[AberrantTestReport]:
    """Run the aberrant-expression battery against a single fitted model.

    The model must implement ``predict(X)`` taking an (n_cells, len(rna_genes))
    array. Cluster_col must exist on ``adata.obs`` so we can split host vs
    source populations.
    """
    if probes is None:
        probes = DEFAULT_PROBES
    rng = np.random.default_rng(seed)

    var_names = np.array([str(v) for v in adata.var_names])
    rna_pos = {g: i for i, g in enumerate(rna_genes)}
    full_pos = {g: i for i, g in enumerate(var_names)}
    adt_pos = {n: i for i, n in enumerate(var_names)}
    adt_target_pos = {n: i for i, n in enumerate(adt_names)}

    # Pre-cache cluster column.
    if cluster_col not in adata.obs.columns:
        raise KeyError(f"obs[{cluster_col!r}] missing")
    clusters = adata.obs[cluster_col].astype(str).to_numpy()

    reports: List[AberrantTestReport] = []
    for adt_name, host_subs, source_subs, partners in probes:
        if adt_name not in adt_target_pos:
            continue
        host_clusters = _match_clusters(clusters, host_subs)
        source_clusters = _match_clusters(clusters, source_subs)
        if not host_clusters or not source_clusters:
            continue
        host_mask = np.isin(clusters, host_clusters)
        source_mask = np.isin(clusters, source_clusters)
        host_cells = np.where(host_mask)[0]
        source_cells = np.where(source_mask)[0]
        if len(host_cells) == 0 or len(source_cells) == 0:
            continue
        if len(host_cells) > max_host_cells:
            host_cells = rng.choice(host_cells, size=max_host_cells, replace=False)
            host_cells.sort()

        # Build the host-cell RNA input
        X_host = _slice_rna(adata, host_cells, rna_genes, rna_pos)
        baseline_pred = _safe_predict(model, X_host)
        adt_idx_in_pred = adt_target_pos[adt_name]

        # Pick the boost level: high-quantile of partner-gene expression in source cells.
        boost_levels = {}
        for partner in partners:
            if partner not in full_pos or partner not in rna_pos:
                continue
            v = adata.X[source_cells, full_pos[partner]]
            if sp.issparse(v):
                v = np.asarray(v.todense()).ravel()
            else:
                v = np.asarray(v).ravel()
            if v.size == 0:
                continue
            boost_levels[partner] = float(np.quantile(v, boost_quantile))
        if not boost_levels:
            continue

        # Apply the boost
        X_boosted = X_host.copy()
        for partner, level in boost_levels.items():
            X_boosted[:, rna_pos[partner]] = level
        boosted_pred = _safe_predict(model, X_boosted)

        # Truth values
        adt_full_idx = adt_pos[adt_name]
        host_truth = adata.X[host_cells, adt_full_idx]
        source_truth = adata.X[source_cells, adt_full_idx]
        host_truth = _flatten(host_truth)
        source_truth = _flatten(source_truth)

        host_truth_mean = float(host_truth.mean()) if host_truth.size else 0.0
        source_truth_mean = float(source_truth.mean()) if source_truth.size else 0.0
        baseline_mean = float(baseline_pred[:, adt_idx_in_pred].mean())
        boosted_mean = float(boosted_pred[:, adt_idx_in_pred].mean())

        denom = max(source_truth_mean - host_truth_mean, 1e-6)
        lift = (boosted_mean - host_truth_mean) / denom
        reports.append(AberrantTestReport(
            adt_name=adt_name,
            host_cluster=", ".join(host_clusters[:3]) + ("..." if len(host_clusters) > 3 else ""),
            source_cluster=", ".join(source_clusters[:3]) + ("..." if len(source_clusters) > 3 else ""),
            partner_genes=list(boost_levels.keys()),
            n_host_cells=int(len(host_cells)),
            baseline_pred_mean=baseline_mean,
            boosted_pred_mean=boosted_mean,
            source_truth_mean=source_truth_mean,
            host_truth_mean=host_truth_mean,
            lift_fraction=float(lift),
        ))
    return reports


def _slice_rna(adata: ad.AnnData, cell_idx: np.ndarray,
               rna_genes: Sequence[str],
               rna_pos: Dict[str, int]) -> np.ndarray:
    var_names = np.array([str(v) for v in adata.var_names])
    var_pos = {g: i for i, g in enumerate(var_names)}
    cols = [var_pos[g] for g in rna_genes if g in var_pos]
    if len(cols) != len(rna_genes):
        # Some training-RNA genes may be missing; fill missing as zeros.
        out = np.zeros((len(cell_idx), len(rna_genes)), dtype=np.float32)
        # Vectorize the present columns
        present_mask = [g in var_pos for g in rna_genes]
        present_cols = np.array([var_pos[g] for g, m in zip(rna_genes, present_mask) if m])
        present_targets = np.array([i for i, m in enumerate(present_mask) if m])
        if present_cols.size:
            sub = adata.X[cell_idx][:, present_cols]
            if sp.issparse(sub):
                sub = np.asarray(sub.todense())
            out[:, present_targets] = np.asarray(sub, dtype=np.float32)
        return out
    sub = adata.X[cell_idx][:, cols]
    if sp.issparse(sub):
        sub = np.asarray(sub.todense())
    return np.asarray(sub, dtype=np.float32)


def _safe_predict(model, X: np.ndarray) -> np.ndarray:
    try:
        return np.asarray(model.predict(X, cluster_labels=None), dtype=np.float32)
    except TypeError:
        return np.asarray(model.predict(X), dtype=np.float32)


def _flatten(x) -> np.ndarray:
    if sp.issparse(x):
        return np.asarray(x.todense()).ravel()
    return np.asarray(x).ravel()


def report_to_dataframe(reports: List[AberrantTestReport]) -> pd.DataFrame:
    rows = []
    for r in reports:
        rows.append({
            "adt": r.adt_name,
            "host_clusters": r.host_cluster,
            "source_clusters": r.source_cluster,
            "partners_boosted": ",".join(r.partner_genes),
            "n_host_cells": r.n_host_cells,
            "host_adt_truth_mean": round(r.host_truth_mean, 3),
            "source_adt_truth_mean": round(r.source_truth_mean, 3),
            "baseline_pred_mean": round(r.baseline_pred_mean, 3),
            "boosted_pred_mean": round(r.boosted_pred_mean, 3),
            "lift_fraction": round(r.lift_fraction, 3),
        })
    return pd.DataFrame(rows)
