#!/usr/bin/env python3
"""Unified SATAY-UDON: ONE canonical donor-level cell-type-resolved enrichment table per
UDON run, and EVERY view (cluster-grouped heatmap, metadata-field heatmap, comparison,
summary) is derived from it -- so they are consistent by construction.

Results are written next to each UDON result, in <run_dir>/satay/ with consistent names:
  satay_celltype_resolved.tsv   canonical (cluster x celltype x covariate, donor-level p+FDR)
  satay_cluster_enrichment.tsv  program-level (cluster x covariate, donor-level)
  satay_heatmap.png/.pdf        cluster-grouped (top cell types per cluster)
  satay_byfield.png/.pdf        metadata-field-grouped (all cell types)
Cross-run comparison + harmonization doc go to UDON/satay_metadata/.

Donor-level (>=3 distinct donors), all features automatically (no curation). Pure Python."""
import os, sys, glob, numpy as np, pandas as pd
UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR); os.chdir(UDON_DIR)
import satay_metadata_enrichment as SME
from satay_metadata_enrichment import build_features, load_run, enrich, _units, MIN_UNITS, bh_fdr
from satay_heatmap import auto_covars, load_final, compute_pmap, add_age, plot, MIN_PB, TOP_CT
from satay_heatmap_byfield import plot_byfield, P_THR
from scipy.stats import fisher_exact

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
TOPOUT = os.path.join(PB, "UDON")          # cross-run files at UDON root (satay_metadata is ignored)
SATAY_DIRNAME = "SATAY-UDON"               # per-run folder name, written next to each marker_heatmap
MIN_PB_CT = 8
import pseudobulk_protocol as _PP
_CT, _GF, _SP, _SFX = _PP.udon_restriction()   # honours --cell-type via UDON_CELL_TYPE


def celltype_resolved_long(cl, feat, defined, covars, min_pb=MIN_PB_CT):
    """THE canonical table: donor-level Fisher for every (cluster, celltype>=min_pb, covariate).
    Donors are collapsed once per (celltype, covariate) into a donor x cluster membership matrix,
    then the 2x2 is vectorized over clusters -> fast and identical to unit_fisher per cluster."""
    clusters = sorted(cl["cluster"].unique(), key=lambda x: (len(str(x)), str(x)))
    cl_idx = {c: j for j, c in enumerate(clusters)}
    cts = [ct for ct, n in cl["celltype"].value_counts().items() if n >= min_pb]
    recs = []
    for ct in cts:
        T = cl[cl["celltype"] == ct]
        for (key, disp, grp) in covars:
            if key not in feat: continue
            Tg = T[T["Sample"].isin(defined[key])]
            if len(Tg) < min_pb: continue
            fpos = Tg["Sample"].map(feat[key]).fillna(0).astype(int).values
            units = _units(Tg["Sample"].values)
            uniq, uinv = np.unique(units, return_inverse=True)
            member = np.zeros((len(uniq), len(clusters)), bool)              # donor x cluster
            member[uinv, [cl_idx[c] for c in Tg["cluster"].values]] = True
            donor_f = np.zeros(len(uniq), int); np.maximum.at(donor_f, uinv, fpos)  # donor feature+
            Fp = int(donor_f.sum()); ntot = len(uniq)
            if Fp < MIN_UNITS: continue                                       # <3 positive donors -> untestable
            fpos_d = donor_f == 1
            for j, c in enumerate(clusters):
                inc = member[:, j]
                a = int((inc & fpos_d).sum())
                if a < MIN_UNITS:
                    recs.append([str(c), ct, disp, grp, 1.0, a]); continue
                b = int((inc & ~fpos_d).sum()); cc = Fp - a; d = (ntot - (a + b)) - cc
                _, p = fisher_exact([[a, b], [cc, d]], alternative="greater")
                recs.append([str(c), ct, disp, grp, p, a])
    df = pd.DataFrame(recs, columns=["cluster", "celltype", "covariate", "group", "p", "n_donors"])
    df["fdr"] = np.nan
    for g, idx in df.groupby("group").groups.items():
        df.loc[idx, "fdr"] = bh_fdr(df.loc[idx, "p"].values)
    return df


def long_to_format1(long, cl, covars, min_pb=MIN_PB, top_ct=TOP_CT):
    clusters = sorted(cl["cluster"].unique(), key=lambda x: (len(str(x)), str(x)))
    lut = {(r.cluster, r.celltype, r.covariate): r.p for r in long.itertuples()}
    rows = []
    for c in clusters:
        vc = cl[cl["cluster"] == c]["celltype"].value_counts()
        for ct in [ct for ct in vc.index if vc[ct] >= min_pb][:top_ct]:
            rows.append((c, ct))
    M = np.ones((len(rows), len(covars)))
    for i, (c, ct) in enumerate(rows):
        for j, (key, disp, grp) in enumerate(covars):
            M[i, j] = lut.get((str(c), ct, disp), 1.0)
    return rows, M


def long_to_byfield_mats(long, covars):
    clusters = sorted(long["cluster"].unique(), key=lambda x: (len(str(x)), str(x)))
    out = {}
    for (key, disp, grp) in covars:
        sub = long[long["covariate"] == disp]
        if sub.empty: continue
        piv = sub.pivot_table(index="celltype", columns="cluster", values="p", aggfunc="min").reindex(columns=clusters)
        M = piv.fillna(1.0).values
        keep = [i for i in range(len(piv.index)) if (M[i] < P_THR).any()]
        if keep:
            out[disp] = (list(np.array(piv.index)[keep]), M[keep], clusters)
    return out


def run_one(name, run_dir, marker_dir, cl, feat, defined, covars, pmap):
    satay = os.path.join(run_dir, SATAY_DIRNAME); os.makedirs(satay, exist_ok=True)
    long = celltype_resolved_long(cl, feat, defined, covars)
    long.to_csv(os.path.join(satay, "satay_celltype_resolved.tsv"), sep="\t", index=False)
    enr = enrich(cl, feat, defined, name)
    enr.to_csv(os.path.join(satay, "satay_cluster_enrichment.tsv"), sep="\t", index=False)
    rows, M = long_to_format1(long, cl, covars)
    if rows:
        plot(rows, M, covars, marker_dir, f"SATAY-UDON: {name}", os.path.join(satay, "satay_heatmap"),
             pmap=pmap if name == "study_aware_final" else None)
    plot_byfield(long_to_byfield_mats(long, covars), f"SATAY-UDON by metadata field: {name}",
                 os.path.join(satay, "satay_byfield"))
    print(f"  {name}: {len(long)} celltype-resolved tests -> {satay}")
    return long, enr


def runs_list():
    # (name, run_dir = dir of the marker_heatmap, marker_dir = dir of marker_genes, cl)
    sa = os.path.join(PB, "UDON", "study_aware" + _SFX)
    if _CT:                                            # cell-type run: naive cell-type UDON + study-aware
        naive = os.path.join(PB, "UDON", "celltype_" + _CT.replace("/", "-"))
        R = [(f"celltype_{_CT}_naive", naive, os.path.join(naive, "udon_core"),
              load_run(os.path.join(naive, "udon_core/udon_clusters.txt")))]
    elif _SFX:                                         # tagged variant (e.g. gene-filtered): study-aware only
        R = []
    else:
        R = [("matched_sex_platform", os.path.join(PB, "UDON/matched_sex_platform"),
              os.path.join(PB, "UDON/matched_sex_platform/udon_core"),
              load_run(os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt"))),
             ("udon_core", os.path.join(PB, "UDON"), os.path.join(PB, "UDON/udon_core"),   # heatmap is at UDON/
              load_run(os.path.join(PB, "UDON/udon_core/udon_clusters.txt")))]
    R.append(("study_aware_final", sa, sa,
              load_final(os.path.join(sa, "final_program_assignments.tsv"))))
    for d in sorted(glob.glob(os.path.join(sa, "per_study/*/udon_clusters.txt"))):
        st = os.path.basename(os.path.dirname(d)); rd = os.path.dirname(d)
        R.append((f"per_study_{st}", rd, rd, load_run(d)))
    return R


def main():
    feat, defined = build_features(); add_age(feat, defined)
    covars = auto_covars(feat, defined)
    pmap = compute_pmap()
    print(f"SATAY-UDON unified | {len(covars)} covariates | unit={SME.UNIT_MODE} | >= {SME.MIN_UNITS} donors")
    longs = {}
    for name, run_dir, marker_dir, cl in runs_list():
        try:
            longs[name], _ = run_one(name, run_dir, marker_dir, cl, feat, defined, covars, pmap)
        except Exception as e:
            print(f"  {name}: FAILED {e}")

    # harmonization reference doc (next to the cross-run comparison, at UDON root)
    pd.DataFrame([{"field": cat, "original": o, "harmonized": (n if n else "(dropped)")}
                  for cat, m in SME.NORMALIZE.items() for o, n in m.items() if o != n]
                 ).to_csv(os.path.join(TOPOUT, "SATAY-UDON_metadata_harmonization.tsv"), sep="\t", index=False)

    # cross-run comparison from the SAME canonical tables
    rep = open(os.path.join(TOPOUT, "SATAY-UDON_crossrun_comparison.txt"), "w")
    def out(m): print(m); rep.write(m + "\n")
    out(f"=== SATAY-UDON cross-run comparison (FDR<0.05 celltype-resolved, from canonical tables) ===")
    sig = {}
    for name, lg in longs.items():
        s = lg[lg["fdr"] < 0.05]
        sig[name] = set(zip(s["celltype"], s["covariate"]))
        out(f"  {name:24s}: {len(s)} sig (celltype x covariate) | {s['covariate'].nunique()} covariates")
    keys = [k for k in ["study_aware_final", "matched_sex_platform", "udon_core"] if k in sig]
    for r in keys:
        uniq = sig[r] - set().union(*[sig[o] for o in keys if o != r])
        out(f"  UNIQUE to {r}: {sorted(set(c for _, c in uniq))[:12]}")
    rep.close()
    print(f"\nwrote per-run {SATAY_DIRNAME}/ folders (next to each marker_heatmap) + cross-run files at {TOPOUT}")


if __name__ == "__main__":
    main()
