#!/usr/bin/env python3.11
"""rna2grn — build the matched (X=RNA, Y=GRN) training dataset.

Outputs (to the Dropbox rna2grn/matched dir):
  dataset.npz   : X (pseudobulk x gene, CP10k+log1p), Y (pseudobulk x edge),
                  plus gene names, edge ids, and per-row metadata arrays
  meta.csv      : per-row metadata (sample, cellstate, group, annotation, n_cells)
  match_table.csv, unmatched.csv : provenance / coverage
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from altanalyze3.components.rna2grn import grn_io  # noqa: E402

GRN = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/tf_to_gene_connection_scores_log10-NOT_ordered_clusters_COMBINED_AML_Multiome_RC2.txt"
RNA = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_counts_hashed.h5ad"
OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/matched"
CS_COL = "Hs-BM-titrated-reference-centroid"

HOLDOUT_SAMPLES = {"AML-7_CCHMC"}  # held out of the materialized training table


def cp10k_log1p(counts):
    """CP10k + log1p on a (rows x genes) matrix (dense ndarray)."""
    tot = counts.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(counts / tot * 1e4).astype(np.float32)


def main():
    os.makedirs(OUT, exist_ok=True)
    print("loading GRN ...")
    grn = grn_io.load_grn_matrix(GRN)
    print(f"  GRN: {grn.n_edges} edges x {len(grn.columns)} columns")

    print("opening RNA (backed) ...")
    A = ad.read_h5ad(RNA, backed="r")
    cs_vocab = list(pd.unique(A.obs[CS_COL].astype(str)))
    rna_samples = list(pd.unique(A.obs["Sample"].astype(str)))
    pb_index = A.obs_names.tolist()
    pb_pos = {k: i for i, k in enumerate(pb_index)}
    genes = A.var_names.tolist()
    print(f"  RNA: {A.n_obs} pseudobulks x {A.n_vars} genes; {len(cs_vocab)} cell states")

    match = grn_io.match_columns(grn.columns, cs_vocab, rna_samples, pb_index)
    tbl = match.table
    tbl.to_csv(f"{OUT}/match_table.csv", index=False)
    pd.Series(match.unmatched, name="unmatched_grn_column").to_csv(f"{OUT}/unmatched.csv", index=False)
    print(f"  matched rows: {len(tbl)}  (incl {len(match.control_columns)} control cols); "
          f"unmatched: {len(match.unmatched)}")

    # ---- direct (non-control) matches -> unique RNA pseudobulk rows ----
    direct = tbl[tbl["rna_pseudobulk"].notna()].copy()
    # average GRN replicate columns that map to the same RNA pseudobulk
    grn_col_pos = {c: i for i, c in enumerate(grn.columns)}
    rows = []
    for pb, grp in direct.groupby("rna_pseudobulk"):
        cols = grp["grn_column"].tolist()
        yvec = grn.values[:, [grn_col_pos[c] for c in cols]].mean(axis=1)
        r0 = grp.iloc[0]
        rows.append(dict(key=pb, rna_sample=r0["rna_sample"], cell_state=r0["cell_state"],
                         group=r0["group"], n_grn_cols=len(cols), Y=yvec))
    direct_rows = pd.DataFrame(rows)
    print(f"  unique direct RNA pseudobulks: {len(direct_rows)}")

    # pull RNA counts for these pseudobulks (backed row slice)
    pos = [pb_pos[k] for k in direct_rows["key"]]
    Xc = A[pos].to_memory().X
    Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)
    Xn = cp10k_log1p(Xc.astype(np.float64))
    Yd = np.vstack(direct_rows["Y"].values).astype(np.float32)

    # ---- control aggregates (Multiome_WT, RC2_TEA) ----
    # control RNA = summed counts across CCHMC Control-annotated samples per cell state
    obs = A.obs
    ctrl_mask = (obs["Annotation"].astype(str) == "Control") & (obs["Dataset"].astype(str) == "CCHMC")
    ctrl_obs = obs[ctrl_mask]
    print(f"  CCHMC control samples: {ctrl_obs['Sample'].nunique()} ; pseudobulks: {len(ctrl_obs)}")
    # map cell_state -> summed counts vector
    ctrl_pos_by_cs = {}
    for cs, grp in ctrl_obs.groupby(CS_COL, observed=True):
        ctrl_pos_by_cs[str(cs)] = [pb_pos[k] for k in grp.index]
    ctrl_rows = []
    control_tbl = tbl[tbl["rna_pseudobulk"].isna()].copy()
    for _, r in control_tbl.iterrows():
        cs = r["cell_state"]
        if cs not in ctrl_pos_by_cs:
            continue  # no control RNA for this cell state
        ctrl_rows.append(dict(grn_column=r["grn_column"], sample_token=r["sample_token"],
                              cell_state=cs, group=r["group"]))
    # build matrices for controls
    cX_rows, cY_rows, cmeta = [], [], []
    for cr in ctrl_rows:
        positions = ctrl_pos_by_cs[cr["cell_state"]]
        sub = A[positions].to_memory().X
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        summed = sub.sum(axis=0, keepdims=True).astype(np.float64)
        cX_rows.append(summed)
        cY_rows.append(grn.values[:, grn_col_pos[cr["grn_column"]]])
        cmeta.append(cr)
    if cX_rows:
        cXc = np.vstack(cX_rows)
        cXn = cp10k_log1p(cXc)
        cY = np.vstack(cY_rows).astype(np.float32)
    else:
        cXn = np.zeros((0, len(genes)), np.float32); cY = np.zeros((0, grn.n_edges), np.float32)

    # ---- assemble metadata ----
    direct_meta = pd.DataFrame(dict(
        key=direct_rows["key"], rna_sample=direct_rows["rna_sample"],
        cell_state=direct_rows["cell_state"], group=direct_rows["group"],
        n_grn_cols=direct_rows["n_grn_cols"], source="direct",
    ))
    control_meta = pd.DataFrame([dict(
        key=f"{cr['sample_token']}|{cr['cell_state']}", rna_sample=cr["sample_token"],
        cell_state=cr["cell_state"], group=cr["group"], n_grn_cols=1, source="control_agg",
    ) for cr in cmeta]) if cmeta else pd.DataFrame(columns=direct_meta.columns)

    X = np.vstack([Xn, cXn]).astype(np.float32)
    Y = np.vstack([Yd, cY]).astype(np.float32)
    meta = pd.concat([direct_meta, control_meta], ignore_index=True)
    meta["heldout_demo"] = meta["rna_sample"].isin(HOLDOUT_SAMPLES)

    assert X.shape[0] == Y.shape[0] == len(meta), (X.shape, Y.shape, len(meta))
    print(f"\nFINAL dataset: X {X.shape}  Y {Y.shape}  rows {len(meta)}")
    print("rows by group:\n", meta.groupby(["group", "source"]).size())
    print("n unique samples:", meta["rna_sample"].nunique(),
          " n unique cellstates:", meta["cell_state"].nunique())

    np.savez_compressed(
        f"{OUT}/dataset.npz",
        X=X, Y=Y,
        genes=np.array(genes, dtype=object),
        edge_ids=np.array(grn.edge_ids, dtype=object),
        edge_tf=grn.edges["TF"].values.astype(object),
        edge_gene=grn.edges["Gene"].values.astype(object),
        keys=meta["key"].values.astype(object),
        sample=meta["rna_sample"].values.astype(object),
        cell_state=meta["cell_state"].values.astype(object),
        group=meta["group"].values.astype(object),
        source=meta["source"].values.astype(object),
        heldout_demo=meta["heldout_demo"].values.astype(bool),
    )
    meta.to_csv(f"{OUT}/meta.csv", index=False)
    print(f"\nwrote {OUT}/dataset.npz and meta.csv")


if __name__ == "__main__":
    main()
