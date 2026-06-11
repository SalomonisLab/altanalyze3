#!/usr/bin/env python3
"""Compare/contrast SATAY-UDON results: which (cell-type x covariate) associations are
SHARED across analyses vs UNIQUE to one. Uses the cell-type-resolved SATAY test
(build_matrix) for: per_study (union of 8 studies), matched_sex_platform, udon_core,
study_aware_final. Flags whether unique merged findings sit in batch-driven clusters.
Pure Python. No R."""
import os, glob, numpy as np, pandas as pd
from satay_heatmap import (build_features, add_age, load_run, load_final, build_matrix,
                           auto_covars, PB, OUT)
from satay_metadata_enrichment import bh_fdr

FDR = 0.05


def to_long(name, rows, M, covars):
    recs = []
    for i, (c, ct) in enumerate(rows):
        for j, (key, disp, grp) in enumerate(covars):
            recs.append([name, str(c), ct, disp, grp, M[i, j]])
    df = pd.DataFrame(recs, columns=["result", "cluster", "celltype", "covariate", "group", "p"])
    df["fdr"] = bh_fdr(df["p"].values)
    return df


def main():
    feat, defined = build_features(); add_age(feat, defined)
    COV = auto_covars(feat, defined)
    longs = {}
    for name, rdir in [("matched_sex_platform", os.path.join(PB, "UDON/matched_sex_platform/udon_core")),
                       ("udon_core", os.path.join(PB, "UDON/udon_core"))]:
        cl = load_run(os.path.join(rdir, "udon_clusters.txt"))
        rows, M = build_matrix(cl, feat, defined, COV, 15); longs[name] = to_long(name, rows, M, COV)
    cl = load_final(os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv"))
    rows, M = build_matrix(cl, feat, defined, COV, 15); longs["study_aware_final"] = to_long("study_aware_final", rows, M, COV)
    ps = []
    for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt"))):
        st = os.path.basename(os.path.dirname(d)); cl = load_run(d)
        for mpb in [6, 4, 3, 2]:
            rows, M = build_matrix(cl, feat, defined, COV, mpb)
            if rows: break
        if rows: ps.append(to_long(f"per_study:{st}", rows, M, COV))
    longs["per_study"] = pd.concat(ps, ignore_index=True)

    # batch dominance for merged clusters (to flag artifacts among unique findings)
    meta = pd.read_csv(os.path.join(PB, "evaluation/outputs/sample_metadata.tsv"), sep="\t", index_col=0)
    def dom(path):
        cl = pd.read_csv(path, sep="\t", index_col=0); cl.columns = ["cluster"]
        cl["D"] = meta.reindex([i.split("__", 1)[1] for i in cl.index])["Dataset"].values
        return {str(c): g["D"].value_counts(normalize=True).iloc[0] for c, g in cl.groupby("cluster")}
    domc = {"matched_sex_platform": dom(os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt")),
            "udon_core": dom(os.path.join(PB, "UDON/udon_core/udon_clusters.txt"))}

    RES = ["per_study", "matched_sex_platform", "udon_core", "study_aware_final"]
    # significant (celltype, covariate) sets per result
    pair = {}; cov = {}
    for r in RES:
        s = longs[r][longs[r]["fdr"] < FDR]
        pair[r] = set(zip(s["celltype"], s["covariate"]))
        cov[r] = set(s["covariate"])

    rep = open(os.path.join(OUT, "SATAY_COMPARE.txt"), "w")
    def out(m): print(m); rep.write(m + "\n")

    out(f"=== detectable COVARIATES per result (FDR<{FDR}, any cell type/cluster) ===")
    for r in RES:
        out(f"  {r:22s}: {len(cov[r]):2d} covariates | {len(pair[r]):3d} (celltype x covariate) findings")
    shared_cov = set.intersection(*cov.values())
    out(f"\nshared by ALL results (covariate): {sorted(shared_cov)}")
    for r in RES:
        uniq = cov[r] - set.union(*[cov[o] for o in RES if o != r])
        out(f"covariate UNIQUE to {r}: {sorted(uniq) if uniq else '(none)'}")

    out(f"\n=== UNIQUE (cell-type x covariate) findings per result ===")
    for r in RES:
        others = set.union(*[pair[o] for o in RES if o != r])
        uniq = pair[r] - others
        out(f"\n--- {r}: {len(uniq)} unique (celltype x covariate) associations ---")
        if not uniq:
            out("   (none)"); continue
        s = longs[r][longs[r]["fdr"] < FDR].copy()
        s["pk"] = list(zip(s["celltype"], s["covariate"]))
        s = s[s["pk"].isin(uniq)].sort_values("fdr")
        seen = set()
        for _, x in s.iterrows():
            k = (x["celltype"], x["covariate"])
            if k in seen: continue
            seen.add(k)
            flag = ""
            if r in domc:
                dd = domc[r].get(x["cluster"], np.nan)
                flag = f"  [cluster {x['cluster']} {100*dd:.0f}% one-study]" if dd > 0.9 else \
                       (f"  [cluster {x['cluster']} {100*dd:.0f}% one-study]" if dd > 0.5 else "")
            out(f"   {x['covariate']:10s} in {x['celltype']:22s} fdr={x['fdr']:.1e}{flag}")
            if len(seen) >= 18: break

    # pairwise overlap matrix
    out(f"\n=== pairwise SHARED (celltype x covariate) findings ===")
    out("           " + " ".join(f"{r[:10]:>11s}" for r in RES))
    for r in RES:
        out(f"  {r[:10]:10s} " + " ".join(f"{len(pair[r]&pair[o]):11d}" for o in RES))
    rep.close()
    pd.concat(longs.values(), ignore_index=True).to_csv(os.path.join(OUT, "satay_celltype_resolved_all.tsv"), sep="\t", index=False)
    print("\nwrote SATAY_COMPARE.txt + satay_celltype_resolved_all.tsv")


if __name__ == "__main__":
    main()
