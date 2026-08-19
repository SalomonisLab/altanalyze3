#!/usr/bin/env python3.11
"""cellHarmony-web — validate the imputation modalities wired into scALABLE.

Checks the four modalities the viewer exposes today:
  lung   : grn   (rna2grn lung bundle) on the 4 human lung references
  marrow : metabolite, lipid (AML), grn (leukemia bundle) on hs_bm_reference

Sections
  A  config: every reference's declared modalities resolve to an existing bundle
  B  lung GRN end-to-end: `_build_imputed_grn_adata` on a REAL lung query built from
     /Users/saljh8/Dropbox/LungMAP/code/ILD/output_file_with_umap.h5ad
  C  lung statistic exactness: the builder's pseudobulk equals the training formula
     (mean over cells of log1p(CP10k), pseudobulk_only.py line 37) to the bit
  D  bone-marrow regression: metabolite / lipid / GRN builders still run, and the
     leukemia GRN path is byte-identical to the pre-change call signature
  E  guard rails: a count matrix sent down the mean-of-log path raises, and an
     unknown statistic raises

Nothing here writes into a job directory or modifies an input file.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")))
from altanalyze3.components.cellHarmony.flask import pipeline  # noqa: E402
from altanalyze3.components.rna2grn.api import load_bundle as load_rna2grn_bundle  # noqa: E402

REGISTRY = Path(__file__).resolve().parents[1] / "reference_config.json"
LUNG_H5AD = "/Users/saljh8/Dropbox/LungMAP/code/ILD/output_file_with_umap.h5ad"
BM_H5AD = ("/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/"
           "final/pseudobulk/pseudobulk_counts_hashed.h5ad")
OUT = "/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/validation"
LUNG_REFS = ("hs_lung_natri_reference", "hs_lung_hlca_reference",
             "hs_lung_cellref_reference", "hs_lung_bpd_sun_reference")
LUNG_CLUSTER_COL = "manual_annotation_1"
BM_CLUSTER_COL = "Hs-BM-titrated-reference-centroid"


# --------------------------------------------------------------------------- #
def _build_lung_query(n_per_group: int, seed: int = 0):
    """Real lung cells, X = the atlas's CP10k+log1p values, obs shaped like a
    cellHarmony job (a `Sample` column plus a cell-state column)."""
    f = h5py.File(LUNG_H5AD, "r")
    def cat(col):
        g = f["obs"][col]
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
        return cats[g["codes"][:]]
    ann, dia, smp = cat("manual_annotation_1"), cat("Diagnosis"), cat("Sample_Name")
    genes = [x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]
    keep_group = np.isin(dia, ["Control", "IPF"])
    rng = np.random.default_rng(seed)
    picked = []
    for grp in ("Control", "IPF"):
        idx = np.where(keep_group & (dia == grp))[0]
        picked.append(rng.choice(idx, size=min(n_per_group, len(idx)), replace=False))
    sel = np.sort(np.concatenate(picked))
    X = f["X"]; indptr = X["indptr"][:]; n_var = int(X.attrs["shape"][1])
    rows, cols, data = [], [], []
    for r, i in enumerate(sel):
        a, b = indptr[i], indptr[i + 1]
        c = X["indices"][a:b]; d = X["data"][a:b]
        rows.append(np.full(len(c), r)); cols.append(c); data.append(d)
    M = sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(sel), n_var), dtype=np.float32)
    obs = pd.DataFrame({"Sample": smp[sel], "Diagnosis": dia[sel],
                        LUNG_CLUSTER_COL: ann[sel]},
                       index=[f"cell{i}" for i in sel])
    f.close()
    q = ad.AnnData(X=M, obs=obs, var=pd.DataFrame(index=pd.Index(genes, dtype=str)))
    q.obsm["X_umap"] = np.zeros((q.n_obs, 2), np.float32)
    return q


def _build_bm_query(n_rows: int, seed: int = 0):
    """Bone-marrow query with RAW COUNTS in layers['counts'] and CP10k+log1p in X —
    the shape cellHarmony_lite produces (cellHarmony_lite.py lines 554-556)."""
    A = ad.read_h5ad(BM_H5AD, backed="r")
    rng = np.random.default_rng(seed)
    sel = np.sort(rng.choice(A.n_obs, size=min(n_rows, A.n_obs), replace=False))
    sub = A[sel].to_memory()
    A.file.close()
    counts = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    counts = np.asarray(counts, dtype=np.float32)
    tot = counts.sum(1, keepdims=True); tot[tot == 0] = 1.0
    q = ad.AnnData(X=np.log1p(counts / tot * 1e4).astype(np.float32),
                   obs=sub.obs.copy(), var=pd.DataFrame(index=sub.var_names.astype(str)))
    q.layers["counts"] = counts
    q.obs_names = [str(x) for x in sub.obs_names]
    q.obsm["X_umap"] = np.zeros((q.n_obs, 2), np.float32)
    return q


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--lung-cells", type=int, default=6000, help="cells per group")
    ap.add_argument("--bm-rows", type=int, default=1200)
    ap.add_argument("--skip-bm", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    report = {"started": datetime.now().isoformat(timespec="seconds"),
              "registry": str(REGISTRY)}
    fail: list = []

    def check(name, ok, detail=""):
        print(f"    {'PASS' if ok else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))
        if not ok:
            fail.append(name)
        return ok

    # ------------------------------------------------------------------ A
    print("\n[A] reference config -> bundle resolution")
    registry = json.loads(REGISTRY.read_text())
    resolved_refs = {}
    rows = []
    for sp_entry in registry["species"]:
        for ref in sp_entry["references"]:
            mods = pipeline._reference_impute_modalities(ref)
            if not mods:
                continue
            entry = pipeline._lookup_reference(sp_entry["id"], ref["id"], REGISTRY)
            resolved_refs[ref["id"]] = entry
            for mod in mods:
                cfg = pipeline._reference_impute_config(entry, mod)
                bp = cfg.get("bundle_path")
                exists = bool(bp) and Path(bp).exists()
                rows.append(dict(species=sp_entry["id"], reference=ref["id"], modality=mod,
                                 bundle_path=bp or "(module default)",
                                 exists=exists if bp else True,
                                 pseudobulk_statistic=cfg.get("pseudobulk_statistic", "")))
                check(f"{ref['id']}/{mod}", exists if bp else True,
                      f"{Path(bp).name if bp else 'module default'}")
    pd.DataFrame(rows).to_csv(out / "modality_bundle_matrix.csv", index=False)
    report["A_bundles"] = rows
    for rid in LUNG_REFS:
        cfg = pipeline._reference_impute_config(resolved_refs[rid], "grn")
        check(f"{rid} grn declares lung statistic",
              cfg.get("pseudobulk_statistic") == "mean_over_cells_of_log1p_cp10k")
        check(f"{rid} grn points at the lung bundle",
              str(cfg.get("bundle_path", "")).endswith("rna2grn_lung_bundle.pkl.gz"))
    bm_grn_cfg = pipeline._reference_impute_config(resolved_refs["hs_bm_reference"], "grn")
    check("hs_bm_reference grn declares NO statistic (keeps sum_counts)",
          not bm_grn_cfg.get("pseudobulk_statistic"))

    # ------------------------------------------------------------------ B
    print(f"\n[B] lung GRN end-to-end ({args.lung_cells} cells/group from {LUNG_H5AD})")
    q = _build_lung_query(args.lung_cells)
    print(f"    query: {q.n_obs} cells x {q.n_vars} genes; "
          f"{q.obs[LUNG_CLUSTER_COL].nunique()} cell states; "
          f"{q.obs['Sample'].nunique()} samples; X max {float(q.X.max()):.2f}")
    entry = resolved_refs["hs_lung_natri_reference"]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tf_adata, edges_adata, summary = pipeline._build_imputed_grn_adata(
            q, entry, cluster_obs_col=LUNG_CLUSTER_COL)
        counts_warn = [w for w in caught if issubclass(w.category, RuntimeWarning)
                       and "you passed counts" in str(w.message)]
    print(f"    summary: {json.dumps({k: v for k, v in summary.items() if k != 'neighbors'}, default=str)[:400]}")
    check("statistic used = mean_over_cells_of_log1p_cp10k",
          summary.get("pseudobulk_statistic") == "mean_over_cells_of_log1p_cp10k")
    check("no counts-path warning fired", not counts_warn)
    check("TF activity adata = cells x TFs", tf_adata.shape == (q.n_obs, 221),
          str(tf_adata.shape))
    check("edge adata = pseudobulks x 57307 edges",
          edges_adata.shape[1] == 57307 and edges_adata.shape[0] == summary["n_pseudobulks"],
          str(edges_adata.shape))
    check("edge values finite and non-negative",
          bool(np.isfinite(edges_adata.X).all() and (edges_adata.X >= 0).all()))
    check("TF activity finite", bool(np.isfinite(tf_adata.X).all()))
    check("TF activity varies across cells", float(tf_adata.X.std()) > 0,
          f"sd={float(tf_adata.X.std()):.4f}")
    check("edge adata carries pseudobulk_method for the group test",
          edges_adata.uns.get("pseudobulk_method") == "pseudobulk")
    check("summary names the lung reference", summary.get("grn_reference") == "lung",
          str(summary.get("grn_reference")))
    report["B_lung"] = {k: v for k, v in summary.items() if k != "neighbors"}
    report["B_lung"]["tf_shape"] = list(tf_adata.shape)
    report["B_lung"]["edge_shape"] = list(edges_adata.shape)

    # ------------------------------------------------------------------ C
    print("\n[C] lung statistic exactness vs the training formula")
    bundle = load_rna2grn_bundle(pipeline._reference_impute_config(entry, "grn")["bundle_path"])
    key = pipeline._pseudobulk_group_key(q.obs, LUNG_CLUSTER_COL)
    Xd = q.X.toarray() if sp.issparse(q.X) else np.asarray(q.X)
    groups = pd.unique(key.values)
    manual = np.vstack([Xd[key.values == g].mean(0) for g in groups])   # pseudobulk_only.py:37
    manual_df = pd.DataFrame(manual, index=groups, columns=list(q.var_names))
    manual_pred = bundle.predict_from_dataframe(manual_df, normalized=True).predictions
    builder_pred = pd.DataFrame(edges_adata.X, index=edges_adata.obs_names,
                                columns=edges_adata.var["raw_feature_name"].astype(str))
    common = manual_pred.index.intersection(builder_pred.index)
    check("every pseudobulk matched", len(common) == len(groups), f"{len(common)}/{len(groups)}")
    diff = float(np.abs(manual_pred.loc[common].to_numpy()
                        - builder_pred.loc[common].to_numpy()).max())
    check("builder pseudobulk == mean over cells of log1p(CP10k)", diff < 1e-5,
          f"max|diff| = {diff:.3e} (tol 1e-5, float32 h5ad storage)")
    report["C_statistic_exactness"] = dict(n_pseudobulks=len(common), max_abs_diff=diff)

    # ------------------------------------------------------------------ D
    if not args.skip_bm:
        print(f"\n[D] bone-marrow regression ({args.bm_rows} rows from {BM_H5AD})")
        bm = _build_bm_query(args.bm_rows)
        bm_entry = resolved_refs["hs_bm_reference"]
        cl = BM_CLUSTER_COL if BM_CLUSTER_COL in bm.obs.columns else None
        n_groups = int(pipeline._pseudobulk_group_key(bm.obs, cl).nunique())
        n_samples = int(pipeline._pseudobulk_group_key(bm.obs, None).nunique())
        print(f"    query: {bm.n_obs} rows x {bm.n_vars} genes; counts layer: "
              f"{'counts' in bm.layers}; cell-state col: {cl}")
        print(f"    NOTE this h5ad is itself pseudobulk-level (1 row per sample x cell-state), "
              f"so sample x cell-state grouping yields {n_groups} groups for {bm.n_obs} rows "
              f"(no collapse). Collapse is exercised twice: section B (12,000 lung cells -> "
              f"{summary['n_pseudobulks']} pseudobulks) and the sample-level check below "
              f"({bm.n_obs} rows -> {n_samples} samples).")
        bm_report = {"n_rows": int(bm.n_obs), "n_sample_x_state_groups": n_groups,
                     "n_samples": n_samples}
        for mod, fn in (("metabolite", pipeline._build_imputed_metabolite_adata),
                        ("lipid", pipeline._build_imputed_lipid_aml_adata)):
            v, d, s = fn(bm, bm_entry, cluster_obs_col=cl)
            check(f"{mod} viewer adata = rows x features", v.shape[0] == bm.n_obs, str(v.shape))
            check(f"{mod} differential rows == unique sample x cell-state groups",
                  d.uns.get("pseudobulk_method") == "pseudobulk" and d.shape[0] == n_groups,
                  f"{d.shape[0]} rows vs {n_groups} groups")
            check(f"{mod} values finite", bool(np.isfinite(v.X).all()))
            bm_report[mod] = dict(viewer_shape=list(v.shape), diff_shape=list(d.shape))
        # real collapse through the shared helper: group by SAMPLE only
        vC, dC, _ = pipeline._build_imputed_metabolite_adata(bm, bm_entry, cluster_obs_col=None)
        check("collapse: sample-level differential rows == unique samples",
              dC.shape[0] == n_samples and n_samples < bm.n_obs,
              f"{bm.n_obs} rows -> {dC.shape[0]} rows, {n_samples} samples")
        check("collapse: viewer stays per-row (broadcast back to cells)",
              vC.shape[0] == bm.n_obs, str(vC.shape))
        bm_report["collapse_sample_level"] = dict(rows=int(bm.n_obs),
                                                  groups=int(dC.shape[0]))
        gtf, gedges, gsum = pipeline._build_imputed_grn_adata(bm, bm_entry, cluster_obs_col=cl)
        check("marrow GRN keeps sum_counts",
              gsum.get("pseudobulk_statistic") == "sum_counts", str(gsum.get("pseudobulk_statistic")))
        check("marrow GRN edges = 7486", gedges.shape[1] == 7486, str(gedges.shape))
        check("marrow GRN TFs = 217", gtf.shape[1] == 217, str(gtf.shape))
        bm_report["grn"] = dict(tf_shape=list(gtf.shape), edge_shape=list(gedges.shape),
                                statistic=gsum.get("pseudobulk_statistic"))

        # byte-identity: the new kwarg must not change the pre-change call
        bmb = load_rna2grn_bundle(pipeline._reference_impute_config(bm_entry, "grn")["bundle_path"])
        bm.obs["_pb_group"] = pipeline._pseudobulk_group_key(bm.obs, cl).values
        old_call = bmb.predict_from_adata(bm, groupby="_pb_group", layer="counts")
        new_call = bmb.predict_from_adata(bm, groupby="_pb_group", layer="counts",
                                          pseudobulk_statistic="sum_counts")
        del bm.obs["_pb_group"]
        identical = bool(np.array_equal(old_call.predictions.to_numpy(),
                                        new_call.predictions.to_numpy()))
        check("leukemia GRN path bit-identical to the pre-change signature", identical,
              f"{old_call.predictions.size} values")
        bm_report["bit_identical_values"] = int(old_call.predictions.size)
        report["D_marrow"] = bm_report
    else:
        report["D_marrow"] = "skipped"

    # ------------------------------------------------------------------ E
    print("\n[E] guard rails")
    try:
        bundle.predict_from_adata(_build_bm_query(50), groupby=None, layer="counts",
                                  pseudobulk_statistic="mean_over_cells_of_log1p_cp10k")
        check("count matrix down the mean-of-log path raises", False, "no error")
    except ValueError as exc:
        check("count matrix down the mean-of-log path raises", "CP10k+log1p" in str(exc),
              type(exc).__name__)
    try:
        bundle.predict_from_adata(q[:20], groupby=None, pseudobulk_statistic="nonsense")
        check("unknown statistic raises", False, "no error")
    except ValueError as exc:
        check("unknown statistic raises", "unknown pseudobulk_statistic" in str(exc),
              type(exc).__name__)

    report["finished"] = datetime.now().isoformat(timespec="seconds")
    report["failures"] = fail
    report["status"] = "PASS" if not fail else "FAIL"
    (out / "impute_modalities_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nSTATUS: {report['status']}" + (f"  failures: {fail}" if fail else ""))
    print(f"wrote {out / 'impute_modalities_report.json'}")
    print(f"wrote {out / 'modality_bundle_matrix.csv'}")
    return 0 if not fail else 1


if __name__ == "__main__":
    raise SystemExit(main())
