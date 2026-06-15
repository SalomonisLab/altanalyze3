"""rna2grn.training — fit the per-edge regulon-activity regression and write the
bundle.

Each TF->target edge's activity is regressed (ridge) on three standardized
features of the pseudobulk: the target gene's expression, the TF's expression,
and the TF's regulon mean. The fit is vectorized across all edges. The resulting
bundle is tiny (per-edge coefficients + standardization stats + a sparse
regulon-averaging matrix), well under 1 MB.
"""
from __future__ import annotations

import gzip
import pickle
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp

from .model import RegulonEdgeRegressor

DEFAULT_RIDGE = 1.0


def _build_regulon_matrix(tf_list, tf_to_target_cols, n_feat):
    rows, cols, data = [], [], []
    for ti, tf in enumerate(tf_list):
        cols_t = tf_to_target_cols[tf]
        w = 1.0 / len(cols_t)
        for c in cols_t:
            rows.append(ti); cols.append(c); data.append(w)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(tf_list), n_feat))


def build_bundle(
    dataset_npz: str | Path,
    out_path: str | Path,
    *,
    ridge_lambda: float = DEFAULT_RIDGE,
    exclude_samples: Optional[Sequence[str]] = None,
    include_controls: bool = True,
    created_at: str = "",
) -> dict:
    d = np.load(str(dataset_npz), allow_pickle=True)
    X = d["X"].astype(np.float64)
    Y = d["Y"].astype(np.float64)
    genes = d["genes"].astype(str)
    edge_ids = d["edge_ids"].astype(str)
    edge_tf = d["edge_tf"].astype(str)
    edge_gene = d["edge_gene"].astype(str)
    sample = d["sample"].astype(str)
    cell_state = d["cell_state"].astype(str)
    group = d["group"].astype(str)

    keep = np.ones(len(sample), dtype=bool)
    if exclude_samples:
        keep &= ~np.isin(sample, list(exclude_samples))
    if not include_controls:
        keep &= ~np.isin(group, ["Multiome_control", "TEA_control"])
    X, Y = X[keep], Y[keep]

    # ---- feature gene space: TF u target genes present in RNA ----
    gene_pos = {g: i for i, g in enumerate(genes)}
    feat_genes = sorted(set(edge_tf) | set(edge_gene))
    feat_genes = [g for g in feat_genes if g in gene_pos]
    fpos = {g: i for i, g in enumerate(feat_genes)}
    # X is CP10k over the whole transcriptome + log1p (the normalization cellHarmony
    # produces), so the model consumes cellHarmony-normalized expression directly with no
    # re-normalization at inference. All-gene CP10k is also slightly more accurate than
    # feature-gene CP10k (LOSO R^2 0.773 vs 0.731) at the same speed.
    Xf = X[:, np.array([gene_pos[g] for g in feat_genes])]   # (n, n_feat) feature space
    n_feat = len(feat_genes)
    E = len(edge_tf)

    target_col = np.array([fpos[g] for g in edge_gene])
    tf_col = np.array([fpos[g] for g in edge_tf])
    tf_list = sorted(set(edge_tf))
    tf_pos = {t: i for i, t in enumerate(tf_list)}
    tf_of_edge = np.array([tf_pos[t] for t in edge_tf])
    tf_to_target_cols = {t: [] for t in tf_list}
    for e in range(E):
        tf_to_target_cols[edge_tf[e]].append(int(target_col[e]))
    regmat = _build_regulon_matrix(tf_list, tf_to_target_cols, n_feat)

    # ---- per-edge features ----
    t = Xf[:, target_col]                                    # (n, E)
    f = Xf[:, tf_col]
    r = np.asarray(Xf @ regmat.T)[:, tf_of_edge]
    feats = [t, f, r]
    n = Xf.shape[0]
    feat_mu = np.zeros((E, 3)); feat_sd = np.ones((E, 3))
    cols = [np.ones((n, E))]
    for j, A in enumerate(feats):
        mu = A.mean(0); sd = A.std(0); sd[sd < 1e-8] = 1.0
        feat_mu[:, j] = mu; feat_sd[:, j] = sd
        cols.append((A - mu) / sd)
    Dsg = np.stack(cols, axis=2)                             # (n, E, 4)
    AtA = np.einsum("nef,neg->efg", Dsg, Dsg)
    Aty = np.einsum("nef,ne->ef", Dsg, Y)
    reg = ridge_lambda * np.eye(4)[None].repeat(E, 0); reg[:, 0, 0] = 0.0
    coefs = np.linalg.solve(AtA + reg, Aty)                  # (E, 4)

    model = RegulonEdgeRegressor(
        feature_genes=feat_genes, edge_tf=edge_tf, edge_gene=edge_gene,
        tf_col=tf_col, target_col=target_col, tf_of_edge=tf_of_edge,
        regulon_matrix=regmat, coefs=coefs, feat_mu=feat_mu, feat_sd=feat_sd,
        tf_list=tf_list,
    )

    metadata = {
        "modality": "grn",
        "approach": "per_edge_regulon_regression",
        "features": ["target_expr", "tf_expr", "regulon_mean"],
        "ridge_lambda": ridge_lambda,
        "normalization": "cp10k_log1p",
        "n_train_pseudobulks": int(Xf.shape[0]),
        "n_feature_genes": n_feat,
        "n_edges": int(E),
        "n_tfs": len(tf_list),
        "include_controls": include_controls,
        "excluded_samples": list(exclude_samples or []),
        "created_at": created_at,
        "target_scaling": {"mode": "none"},
        "edge_tf": edge_tf.tolist(),
        "edge_gene": edge_gene.tolist(),
        "tf_list": tf_list,
    }
    bundle = {
        "model": model,
        "scaler_x": None,
        "scaler_y": None,
        "X_columns": list(feat_genes),
        "Y_columns": list(edge_ids),
        "metadata": metadata,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if out_path.suffix == ".gz" else open
    with opener(out_path, "wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return bundle
