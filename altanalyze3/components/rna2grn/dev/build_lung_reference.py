#!/usr/bin/env python3.11
"""rna2grn — build the LUNG (LungMAP IPF/control) reference: matched dataset + bundle.

Inputs (both are read-only; neither is modified):
  GRN  /Users/saljh8/Dropbox/LungMAP/GRN/TF_to_Gene_connection_scores_log10-NOT_ordered_clusters_ALL_GENES-threshold-1.txt
  RNA  /Users/saljh8/Dropbox/LungMAP/code/ILD/AltAnalyze-create-cH-reference/ExpressionInput/exp.pseudobulks-IPF-control.txt

The RNA matrix is genes x pseudobulks. Its values are the MEAN over cells of
scanpy log1p(CP10k) (built by /Users/saljh8/Dropbox/LungMAP/code/ILD/pseudobulk_only.py
line 37 over an X normalized by /Users/saljh8/Dropbox/LungMAP/code/ILD/normalize.py),
so the matrix is ALREADY in the model's CP10k+log1p input space and is written to
the dataset without re-normalization. The statistic is a mean-of-log, not the
leukemia reference's log-of-summed-counts; the difference is recorded in the
bundle metadata as ``pseudobulk_statistic``.

No internal validation split: every matched pseudobulk is a training row
(``--no-holdout`` is the only supported mode for this reference).

Outputs (all absolute paths printed at the end):
  <out>/matched/dataset.npz        X, Y, genes, edge ids, per-row metadata
  <out>/matched/meta.csv           per-row metadata
  <out>/matched/match_table.csv    GRN column -> pseudobulk provenance
  <out>/matched/unmatched.csv      GRN columns with no pseudobulk partner
  <out>/matched/parameters.json    every input path, parameter and version
  <bundle>                         the trained bundle (default: shipped in-module)
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from altanalyze3.components.rna2grn import grn_io, training  # noqa: E402

GRN = ("/Users/saljh8/Dropbox/LungMAP/GRN/"
       "TF_to_Gene_connection_scores_log10-NOT_ordered_clusters_ALL_GENES-threshold-1.txt")
RNA = ("/Users/saljh8/Dropbox/LungMAP/code/ILD/AltAnalyze-create-cH-reference/"
       "ExpressionInput/exp.pseudobulks-IPF-control.txt")
OUT = "/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn"
BUNDLE = str(Path(__file__).resolve().parents[1] / "rna2grn_lung_bundle.pkl.gz")
PSEUDOBULK_STATISTIC = "mean_over_cells_of_log1p_cp10k"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grn", default=GRN)
    ap.add_argument("--rna", default=RNA)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--bundle", default=BUNDLE)
    ap.add_argument("--ridge-lambda", type=float, default=training.DEFAULT_RIDGE)
    args = ap.parse_args(argv)

    matched_dir = Path(args.out) / "matched"
    matched_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now().isoformat(timespec="seconds")

    print(f"GRN : {args.grn}")
    print(f"      mtime {datetime.fromtimestamp(os.path.getmtime(args.grn)).isoformat(timespec='seconds')}")
    grn = grn_io.load_grn_matrix(args.grn)
    print(f"      {grn.n_edges} edges x {len(grn.columns)} columns")
    n_dup_edges = len(grn.edge_ids) - len(set(grn.edge_ids))
    if n_dup_edges:
        raise ValueError(f"GRN file has {n_dup_edges} duplicate TF|Gene rows")

    print(f"RNA : {args.rna}")
    print(f"      mtime {datetime.fromtimestamp(os.path.getmtime(args.rna)).isoformat(timespec='seconds')}")
    rna = pd.read_csv(args.rna, sep="\t", index_col=0)
    rna.index = [str(i).strip() for i in rna.index]
    rna.columns = [str(c).strip() for c in rna.columns]
    print(f"      {rna.shape[0]} genes x {rna.shape[1]} pseudobulks")
    dup_genes = rna.index[rna.index.duplicated()].tolist()
    if dup_genes:
        raise ValueError(f"RNA matrix has {len(dup_genes)} duplicate gene symbols, "
                         f"e.g. {dup_genes[:5]}")
    if not np.isfinite(rna.to_numpy(dtype=np.float64)).all():
        raise ValueError("RNA matrix contains non-finite values")

    # ---- match GRN columns to pseudobulk columns -------------------------
    match = grn_io.match_lung_columns(grn.columns, rna.columns.tolist())
    tbl = match.table
    tbl.to_csv(matched_dir / "match_table.csv", index=False)
    pd.Series(match.unmatched, name="unmatched_grn_column").to_csv(
        matched_dir / "unmatched.csv", index=False)
    unused = sorted(set(rna.columns) - set(tbl["rna_pseudobulk"]))
    pd.Series(unused, name="pseudobulk_without_grn").to_csv(
        matched_dir / "pseudobulk_without_grn.csv", index=False)
    print(f"\nMATCH: {len(tbl)}/{len(grn.columns)} GRN columns matched "
          f"({100*len(tbl)/len(grn.columns):.1f}%); {len(match.unmatched)} unmatched")
    print(f"       {len(unused)}/{rna.shape[1]} RNA pseudobulks have no GRN partner "
          f"(not used as training rows)")
    if len(tbl) == 0:
        raise ValueError("no GRN column matched a pseudobulk")
    if tbl["rna_pseudobulk"].duplicated().any():
        raise ValueError("two GRN columns matched the same pseudobulk")

    # ---- assemble X (pseudobulks x genes) and Y (pseudobulks x edges) ----
    grn_col_pos = {c: i for i, c in enumerate(grn.columns)}
    X = rna[tbl["rna_pseudobulk"].tolist()].to_numpy(dtype=np.float64).T
    Y = grn.values[:, [grn_col_pos[c] for c in tbl["grn_column"]]].T.astype(np.float64)
    genes = rna.index.tolist()

    exp_rows = len(tbl)
    assert X.shape == (exp_rows, len(genes)), (X.shape, exp_rows, len(genes))
    assert Y.shape == (exp_rows, grn.n_edges), (Y.shape, exp_rows, grn.n_edges)

    # invariant: every row's X column really is the pseudobulk the table names
    probe = tbl.index[:: max(1, len(tbl) // 10)]
    for i in probe:
        pb = tbl.loc[i, "rna_pseudobulk"]
        assert np.allclose(X[i], rna[pb].to_numpy(dtype=np.float64)), pb
        col = tbl.loc[i, "grn_column"]
        assert np.allclose(Y[i], grn.values[:, grn_col_pos[col]]), col

    # ---- gene-space coverage (reported, never silently trimmed) ----------
    gset = set(genes)
    tfs = sorted(set(grn.edges["TF"]))
    tgs = sorted(set(grn.edges["Gene"]))
    tf_in = [t for t in tfs if t in gset]
    tg_in = [t for t in tgs if t in gset]
    edge_ok = grn.edges["TF"].isin(gset) & grn.edges["Gene"].isin(gset)
    print(f"\nCOVERAGE: TFs {len(tf_in)}/{len(tfs)} ({100*len(tf_in)/len(tfs):.1f}%); "
          f"targets {len(tg_in)}/{len(tgs)} ({100*len(tg_in)/len(tgs):.1f}%); "
          f"edges with both present {int(edge_ok.sum())}/{len(edge_ok)} "
          f"({100*edge_ok.mean():.2f}%)")
    if edge_ok.mean() < 0.90:
        print("  WARNING: edge retention below 90% — a finding, not a footnote")

    meta = pd.DataFrame(dict(
        key=tbl["rna_pseudobulk"].values,
        rna_sample=tbl["rna_sample"].values,
        cell_state=tbl["cell_state"].values,
        group=tbl["group"].values,
        n_grn_cols=1,
        source="direct",
        heldout_demo=False,
    ))
    print("\nrows by group:\n", meta.groupby("group").size().to_string())
    print(f"unique cell states: {meta['cell_state'].nunique()}")

    npz_path = matched_dir / "dataset.npz"
    np.savez_compressed(
        npz_path,
        X=X.astype(np.float32), Y=Y.astype(np.float32),
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
    meta.to_csv(matched_dir / "meta.csv", index=False)
    print(f"\nwrote {npz_path}")
    print(f"wrote {matched_dir / 'meta.csv'}")

    # ---- fit the shipped per-edge elastic-net (l1_ratio=0) model ---------
    print(f"\nfitting bundle (ridge_lambda={args.ridge_lambda}, "
          f"all {exp_rows} pseudobulks, no holdout) ...")
    bundle = training.build_bundle(
        npz_path, args.bundle,
        ridge_lambda=args.ridge_lambda,
        exclude_samples=None,
        include_controls=True,
        created_at=started,
    )
    md = bundle["metadata"]
    md["reference"] = "lung"
    md["tissue"] = "human lung (LungMAP IPF vs control)"
    md["pseudobulk_statistic"] = PSEUDOBULK_STATISTIC
    md["input_already_normalized"] = True
    md["grn_source"] = args.grn
    md["rna_source"] = args.rna
    md["internal_validation"] = "none (all pseudobulks used for training, by request)"
    md["reference_sample"] = meta["rna_sample"].tolist()
    md["reference_cell_state"] = meta["cell_state"].tolist()
    md["reference_key"] = meta["key"].tolist()
    import gzip
    import pickle
    with gzip.open(args.bundle, "wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    size_mb = os.path.getsize(args.bundle) / 1e6
    print(f"wrote {args.bundle}  ({size_mb:.2f} MB)")

    params = dict(
        started=started, finished=datetime.now().isoformat(timespec="seconds"),
        script=str(Path(__file__).resolve()),
        python=sys.version.split()[0], platform=platform.platform(),
        numpy=np.__version__, pandas=pd.__version__,
        grn_path=args.grn, rna_path=args.rna,
        grn_mtime=datetime.fromtimestamp(os.path.getmtime(args.grn)).isoformat(timespec="seconds"),
        rna_mtime=datetime.fromtimestamp(os.path.getmtime(args.rna)).isoformat(timespec="seconds"),
        out_dir=str(matched_dir), bundle_path=str(Path(args.bundle).resolve()),
        bundle_size_mb=round(size_mb, 3),
        ridge_lambda=args.ridge_lambda,
        n_grn_columns=len(grn.columns), n_matched=len(tbl),
        n_unmatched=len(match.unmatched),
        n_rna_pseudobulks=int(rna.shape[1]), n_rna_pseudobulks_unused=len(unused),
        n_edges=int(grn.n_edges), n_tfs=len(tfs), n_targets=len(tgs),
        n_feature_genes=int(md["n_feature_genes"]),
        n_train_pseudobulks=int(md["n_train_pseudobulks"]),
        pseudobulk_statistic=PSEUDOBULK_STATISTIC,
        internal_validation="none (all pseudobulks used for training, by request)",
    )
    (matched_dir / "parameters.json").write_text(json.dumps(params, indent=2))
    print(f"wrote {matched_dir / 'parameters.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
