#!/usr/bin/env python3
"""Compare SATAY-UDON enrichments between two UDON runs (e.g. naive pooled vs study-aware/SVA).

The two runs have different clusters, so enrichments are compared at the (covariate, cell type)
level: a covariate "enriches" a cell type if it is significant in ANY cluster of that cell type.
Reports SHARED / NEW (in B only) / LOST (in A only), at both the (covariate x celltype) and the
covariate level. For LOST enrichments, optionally flags whether they sat in a batch-driven cluster
in run A (single-study-dominated -> likely a batch artifact the other run removed).

Usage:
  python satay_compare_runs.py --a NAIVE.tsv --b STUDYAWARE.tsv --labels naive,study_aware \
      --fdr 0.05 --out DIR [--a-clusters udon_clusters.txt --sample-metadata sample_metadata.tsv]
"""
import os, sys, argparse, pandas as pd, numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="run A satay_cluster_enrichment.tsv")
    ap.add_argument("--b", required=True, help="run B satay_cluster_enrichment.tsv")
    ap.add_argument("--labels", default="A,B")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    ap.add_argument("--a-clusters", default=None, help="run A udon_clusters.txt (to flag batch-driven LOST)")
    ap.add_argument("--sample-metadata", default=None, help="Sample,Dataset table (for batch dominance)")
    ap.add_argument("--batch-dominance", type=float, default=0.7)
    args = ap.parse_args()
    la, lb = args.labels.split(",")
    os.makedirs(args.out, exist_ok=True)

    da = pd.read_csv(args.a, sep="\t"); db = pd.read_csv(args.b, sep="\t")
    sa = da[da["fdr"] < args.fdr].copy(); sb = db[db["fdr"] < args.fdr].copy()
    A = set(zip(sa["covariate"], sa["celltype"])); B = set(zip(sb["covariate"], sb["celltype"]))
    new, lost, shared = B - A, A - B, A & B

    print(f"=== SATAY enrichment comparison: {la} (A) vs {lb} (B), FDR<{args.fdr} ===")
    print(f"(covariate x celltype): {la}={len(A)}  {lb}={len(B)}  SHARED={len(shared)}  "
          f"NEW(in {lb})={len(new)}  LOST(from {la})={len(lost)}")
    Ca, Cb = set(sa["covariate"]), set(sb["covariate"])
    print(f"covariates with >=1 sig: {la}={len(Ca)}  {lb}={len(Cb)}")
    print(f"  NEW covariates ({lb} only): {sorted(Cb - Ca)}")
    print(f"  LOST covariates ({la} only): {sorted(Ca - Cb)}")

    # batch-dominance per run-A cluster (to label LOST as batch-artifact vs genuine)
    dom = {}
    if args.a_clusters and args.sample_metadata and os.path.exists(args.a_clusters):
        cl = pd.read_csv(args.a_clusters, sep="\t", index_col=0); cl.columns = ["cluster"]
        cl["Sample"] = [str(i).split("__", 1)[1] if "__" in str(i) else str(i) for i in cl.index]
        sm = pd.read_csv(args.sample_metadata, sep="\t", index_col=0)
        cl["Dataset"] = cl["Sample"].map(sm["Dataset"])
        for c, g in cl.groupby("cluster"):
            vc = g["Dataset"].value_counts(normalize=True)
            dom[str(c)] = (vc.index[0], round(float(vc.iloc[0]), 2))

    def rows(pairs, src):
        out = []
        for cov, ct in sorted(pairs):
            sub = src[(src["covariate"] == cov) & (src["celltype"] == ct)].sort_values("fdr")
            r = sub.iloc[0]
            cl = str(r["cluster"]); td = dom.get(cl, ("", np.nan))
            out.append({"covariate": cov, "celltype": ct, "best_cluster": cl, "fdr": r["fdr"],
                        "clusterA_top_study": td[0], "clusterA_dominance": td[1],
                        "batch_driven": (isinstance(td[1], float) and td[1] >= args.batch_dominance)})
        return pd.DataFrame(out)

    new_df = rows(new, sb); lost_df = rows(lost, sa)
    new_df.to_csv(os.path.join(args.out, f"NEW_in_{lb}.tsv"), sep="\t", index=False)
    lost_df.to_csv(os.path.join(args.out, f"LOST_from_{la}.tsv"), sep="\t", index=False)
    pd.DataFrame(sorted(shared), columns=["covariate", "celltype"]).to_csv(
        os.path.join(args.out, "SHARED.tsv"), sep="\t", index=False)

    if dom and len(lost_df):
        nb = int(lost_df["batch_driven"].sum())
        print(f"\nOf {len(lost_df)} LOST enrichments, {nb} sat in batch-driven {la} clusters "
              f"(>= {args.batch_dominance} one-study) -> likely batch artifacts removed by {lb}; "
              f"{len(lost_df)-nb} were in well-mixed clusters (genuine losses).")
    print(f"\nwrote NEW_in_{lb}.tsv, LOST_from_{la}.tsv, SHARED.tsv -> {args.out}")
    if len(new_df):
        print("\ntop NEW (study-aware-only) enrichments:")
        print(new_df.sort_values("fdr").head(10).to_string(index=False))


if __name__ == "__main__":
    main()
