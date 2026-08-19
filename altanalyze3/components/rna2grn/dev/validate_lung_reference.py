#!/usr/bin/env python3.11
"""rna2grn — validate the LUNG reference bundle.

Runs five checks and writes a report. Nothing here changes the shipped bundle.

  A  structural invariants: bundle edge order == GRN file row order, feature-gene
     space, metadata alignment, finite non-negative predictions
  B  public-API round trip: predictions from the shipped CLI/API path equal the
     predictions from the model called directly on the training matrix
  C  RESUBSTITUTION accuracy on the 39 training pseudobulks. This is the
     configuration that was requested (no internal validation, all pseudobulks
     used for training), so these numbers are optimistic BY CONSTRUCTION.
  D  ADDITION (not requested, does not change the shipped model): leave-one-
     pseudobulk-out held-out accuracy, so the resubstitution number in C has an
     honest counterpart. Uses evaluate._Edges, the module's sanctioned
     evaluation core, first proven to reproduce the shipped fit.
  E  the counts-path RuntimeWarning fires for this already-normalized reference
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from altanalyze3.components.rna2grn import Rna2GrnBundle, evaluate  # noqa: E402

DATASET = "/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/matched/dataset.npz"
GRN = ("/Users/saljh8/Dropbox/LungMAP/GRN/"
       "TF_to_Gene_connection_scores_log10-NOT_ordered_clusters_ALL_GENES-threshold-1.txt")
RNA = ("/Users/saljh8/Dropbox/LungMAP/code/ILD/AltAnalyze-create-cH-reference/"
       "ExpressionInput/exp.pseudobulks-IPF-control.txt")
OUT = "/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/validation"


def _r2(true, pred, baseline):
    ss_res = float(((true - pred) ** 2).sum())
    ss_tot = float(((true - baseline) ** 2).sum())
    return 1.0 - ss_res / ss_tot


def _corr_rows(A, B):
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    num = (A * B).sum(1)
    den = np.sqrt((A ** 2).sum(1) * (B ** 2).sum(1))
    out = np.full(len(num), np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default=DATASET)
    ap.add_argument("--grn", default=GRN)
    ap.add_argument("--rna", default=RNA)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--reference", default="lung")
    ap.add_argument("--skip-loo", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    report: dict = {"started": datetime.now().isoformat(timespec="seconds"),
                    "dataset": args.dataset, "reference": args.reference}
    fail: list = []

    b = Rna2GrnBundle.load(reference=args.reference)
    d = np.load(args.dataset, allow_pickle=True)
    X = d["X"].astype(np.float64); Y = d["Y"].astype(np.float64)
    genes = d["genes"].astype(str); edge_ids = d["edge_ids"].astype(str)
    edge_tf = d["edge_tf"].astype(str); edge_gene = d["edge_gene"].astype(str)
    keys = d["keys"].astype(str); cs = d["cell_state"].astype(str)
    grp = d["group"].astype(str)
    n, E = Y.shape
    print(f"bundle : {b.bundle_path}")
    print(f"dataset: {args.dataset}  X {X.shape}  Y {Y.shape}")

    # ---------------------------------------------------------------- A
    print("\n[A] structural invariants")
    src_head = open(args.grn).readline().rstrip("\n").split("\t")
    src = pd.read_csv(args.grn, sep="\t", usecols=[0, 1])
    src_ids = [f"{t}|{g}" for t, g in zip(src.iloc[:, 0].astype(str), src.iloc[:, 1].astype(str))]
    checks = {
        "bundle_edges_equal_grn_file_rows_in_order": list(b.output_edges) == src_ids,
        "dataset_edges_equal_grn_file_rows_in_order": list(edge_ids) == src_ids,
        "n_edges": E == len(src_ids),
        "metadata_edge_tf_aligned": list(np.asarray(b.metadata["edge_tf"], dtype=str)) == list(edge_tf),
        "metadata_edge_gene_aligned": list(np.asarray(b.metadata["edge_gene"], dtype=str)) == list(edge_gene),
        "feature_genes_are_tf_union_target": (
            sorted(b.input_genes) == sorted(set(edge_tf) | set(edge_gene))),
        "all_feature_genes_in_rna": set(b.input_genes).issubset(set(genes)),
        "n_train_pseudobulks_equals_all_rows": b.metadata["n_train_pseudobulks"] == n,
        "no_holdout_recorded": b.metadata.get("internal_validation", "").startswith("none"),
        "X_finite": bool(np.isfinite(X).all()),
        "Y_finite_nonneg": bool(np.isfinite(Y).all() and (Y >= 0).all()),
    }
    for k, v in checks.items():
        print(f"    {'PASS' if v else 'FAIL'}  {k}")
        if not v:
            fail.append(f"A:{k}")
    report["A_invariants"] = {k: bool(v) for k, v in checks.items()}
    report["A_counts"] = dict(n_pseudobulks=int(n), n_edges=int(E),
                              n_feature_genes=len(b.input_genes),
                              n_tfs=len(set(edge_tf)), n_targets=len(set(edge_gene)),
                              n_grn_columns=len(src_head) - 2,
                              n_cell_states=int(pd.Series(cs).nunique()))

    # ---------------------------------------------------------------- B
    print("\n[B] public-API round trip (real entry point)")
    fpos = {g: i for i, g in enumerate(genes)}
    Xf = X[:, np.array([fpos[g] for g in b.input_genes])]
    direct = b.model.predict(Xf)

    # B1 STRICT: the public API fed exactly the dataset matrix must return exactly
    # what the model returns. Same numbers in, so any difference is an alignment or
    # ordering defect. Tolerance 1e-9 (float64 arithmetic on identical inputs).
    df_npz = pd.DataFrame(X, index=list(keys), columns=list(genes))
    api_npz = b.predict_from_dataframe(df_npz, normalized=True).predictions
    same_idx = list(api_npz.index) == list(keys)
    strict_diff = float(np.abs(api_npz.to_numpy() - direct).max())

    # B2 QUANTIZATION: the same call reading the ORIGINAL text file at full float64
    # precision. dataset.npz stores X as float32, so a nonzero difference is expected
    # here; it must stay at the float32 round-off scale (float32 eps 1.2e-7 on
    # expression values up to ~5.5), not larger. Tolerance 1e-5, value reported.
    rna = pd.read_csv(args.rna, sep="\t", index_col=0)
    rna.index = [str(i).strip() for i in rna.index]
    rna.columns = [str(c).strip() for c in rna.columns]
    df = rna[list(keys)].T                       # pseudobulks x genes, as the CLI builds it
    api_pred = b.predict_from_dataframe(df, normalized=True).predictions
    quant_diff = float(np.abs(api_pred.to_numpy() - direct).max())
    x_quant = float(np.abs(df.to_numpy(dtype=np.float64) - X).max())

    finite = bool(np.isfinite(api_pred.to_numpy()).all())
    nonneg = bool((api_pred.to_numpy() >= 0).all())
    print(f"    index order preserved              : {same_idx}")
    print(f"    B1 max|API(npz X) - model(npz X)|  : {strict_diff:.3e}  (tol 1e-9)")
    print(f"    B2 max|API(text X) - model(npz X)| : {quant_diff:.3e}  (tol 1e-5, "
          f"float32 storage; max|X_text - X_npz| = {x_quant:.3e})")
    print(f"    finite / non-negative              : {finite} / {nonneg}")
    for k, v in [("B:index", same_idx), ("B:strict", strict_diff < 1e-9),
                 ("B:float32_quantization", quant_diff < 1e-5),
                 ("B:finite", finite), ("B:nonneg", nonneg)]:
        if not v:
            fail.append(k)
    report["B_api_roundtrip"] = dict(
        index_order_preserved=same_idx,
        strict_max_abs_diff=strict_diff, strict_tol=1e-9,
        float32_quantization_max_abs_diff=quant_diff, quantization_tol=1e-5,
        x_storage_max_abs_diff=x_quant, finite=finite, nonneg=nonneg)
    P = api_pred.to_numpy()

    # ---------------------------------------------------------------- C
    print("\n[C] RESUBSTITUTION accuracy (optimistic by construction)")
    gmean = np.full_like(Y, Y.mean())
    r2_global = _r2(Y, P, gmean)
    prof_r = _corr_rows(Y, P)
    edge_r = _corr_rows(Y.T, P.T)
    resub = dict(
        n_pseudobulks=int(n), n_edges=int(E),
        r2_vs_global_mean=round(r2_global, 4),
        profile_r_median=round(float(np.nanmedian(prof_r)), 4),
        profile_r_min=round(float(np.nanmin(prof_r)), 4),
        edge_r_median=round(float(np.nanmedian(edge_r)), 4),
        edge_r_frac_nan=round(float(np.isnan(edge_r).mean()), 4),
        mean_abs_error=round(float(np.abs(Y - P).mean()), 4),
        true_mean=round(float(Y.mean()), 4), pred_mean=round(float(P.mean()), 4),
    )
    for k, v in resub.items():
        print(f"    {k:24s} {v}")
    report["C_resubstitution"] = resub
    pd.DataFrame({"pseudobulk": keys, "cell_state": cs, "group": grp,
                  "profile_r_resub": prof_r}).to_csv(
        out / "resubstitution_profile_r.csv", index=False)

    # ---------------------------------------------------------------- D
    if not args.skip_loo:
        print("\n[D] ADDITION: leave-one-pseudobulk-out (held out; shipped model unchanged)")
        ed = evaluate._Edges(genes, edge_tf, edge_gene)
        assert ed.feat == list(b.input_genes), "evaluate._Edges feature space differs from bundle"
        Xf_all = X[:, ed.feat_cols]
        coefs_all, st_all = ed.fit(Xf_all, Y)
        full_pred = ed.predict(coefs_all, st_all, Xf_all)
        # compare against the model on the SAME (npz, float32) matrix, so this is an
        # exact-arithmetic check of the evaluation core, not a storage-precision test
        agree = float(np.abs(full_pred - direct).max())
        print(f"    evaluate._Edges reproduces the shipped fit: max|diff| = {agree:.3e}"
              f"  (tol 1e-9)")
        if agree > 1e-9:
            fail.append("D:core_mismatch")
        loo = np.zeros_like(Y)
        for i in range(n):
            m = np.ones(n, bool); m[i] = False
            c, s = ed.fit(Xf_all[m], Y[m])
            loo[i] = ed.predict(c, s, Xf_all[i:i + 1])[0]
            print(f"      fold {i+1:2d}/{n}  {keys[i]}", flush=True)
        loo_prof = _corr_rows(Y, loo)
        loo_edge = _corr_rows(Y.T, loo.T)
        held = dict(
            r2_vs_global_mean=round(_r2(Y, loo, gmean), 4),
            profile_r_median=round(float(np.nanmedian(loo_prof)), 4),
            profile_r_min=round(float(np.nanmin(loo_prof)), 4),
            edge_r_median=round(float(np.nanmedian(loo_edge)), 4),
            mean_abs_error=round(float(np.abs(Y - loo).mean()), 4),
        )
        for k, v in held.items():
            print(f"    {k:24s} {v}")
        report["D_leave_one_pseudobulk_out"] = held
        pd.DataFrame({"pseudobulk": keys, "cell_state": cs, "group": grp,
                      "profile_r_resub": prof_r, "profile_r_loo": loo_prof}).to_csv(
            out / "loo_profile_r.csv", index=False)
    else:
        report["D_leave_one_pseudobulk_out"] = "skipped"

    # ---------------------------------------------------------------- E
    print("\n[E] counts-path warning")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        b.predict_from_dataframe(df.iloc[:2], normalized=False)
        fired = any(issubclass(x.category, RuntimeWarning) for x in w)
    print(f"    RuntimeWarning raised on counts path: {fired}")
    if not fired:
        fail.append("E:no_warning")
    report["E_counts_path_warning"] = bool(fired)

    report["finished"] = datetime.now().isoformat(timespec="seconds")
    report["failures"] = fail
    report["status"] = "PASS" if not fail else "FAIL"
    (out / "validation_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nSTATUS: {report['status']}"
          + (f"  failures: {fail}" if fail else ""))
    print(f"wrote {out / 'validation_report.json'}")
    print(f"wrote {out / 'resubstitution_profile_r.csv'}")
    if not args.skip_loo:
        print(f"wrote {out / 'loo_profile_r.csv'}")
    return 0 if not fail else 1


if __name__ == "__main__":
    raise SystemExit(main())
