#!/usr/bin/env python3
"""Did the metadata update (N/A->0, +12 gene columns) lose any prior associations?
Run the IDENTICAL enrichment on the OLD vs NEW metadata and diff the significant sets
(FDR<0.05), restricted to features present in BOTH. Reports lost / gained / shared."""
import os, glob, pandas as pd
import satay_metadata_enrichment as SME

PB = SME.PB
OLD = os.path.join(PB, "AML_harmonized_metadata.xlsx")
NEW = os.path.join(PB, "UDON", "AML_harmonized_metadata.xlsx")
FDR = 0.05


def feats(path):
    SME.XLSX = path
    return SME.build_features()


def sigset(feat, defined):
    runs = [("matched_sex_platform", os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt")),
            ("udon_core", os.path.join(PB, "UDON/udon_core/udon_clusters.txt"))]
    runs += [(f"per_study:{os.path.basename(os.path.dirname(d))}", d)
             for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt")))]
    S = set()
    for name, p in runs:
        cl = SME.load_run(p)
        e = SME.enrich(cl, feat, defined, name)
        s = e[e["fdr"] < FDR]
        rt = "per_study" if name.startswith("per_study") else name
        for _, r in s.iterrows():
            S.add((rt, r["category"], r["feature"]))     # feature-level (cluster ids identical across old/new)
    return S


fo, do = feats(OLD); fn, dn = feats(NEW)
shared = set(fo) & set(fn)                                # features in BOTH metadata versions
print(f"shared features: {len(shared)} | old-only: {len(set(fo)-set(fn))} | new-only: {len(set(fn)-set(fo))}")

So = sigset(fo, do); Sn = sigset(fn, dn)
# restrict to shared features
shared_names = {k[1] for k in shared}   # mutation gene / clinical level names
def restrict(S):
    return {(rt, cat, f) for (rt, cat, f) in S if cat == "celltype" or f in shared_names}
So_r, Sn_r = restrict(So), restrict(Sn)

lost = So_r - Sn_r
gained = Sn_r - So_r
print(f"\nsignificant (run, category, feature) associations on SHARED features:")
print(f"  OLD: {len(So_r)} | NEW: {len(Sn_r)} | shared: {len(So_r & Sn_r)}")
print(f"  LOST (in OLD, gone in NEW): {len(lost)}")
print(f"  GAINED (new in NEW): {len(gained)}")

def show(tag, S):
    print(f"\n{tag}:")
    for rt, cat, f in sorted(S):
        if cat == "mutation" or cat == "clinical":
            print(f"   {rt:22s} {cat:9s} {f}")
    cts = [(rt, f) for (rt, cat, f) in S if cat == "celltype"]
    if cts: print(f"   (+ {len(cts)} celltype associations)")

if lost: show("LOST associations", lost)
else: print("\n** No mutation/clinical associations lost. **")
if gained: show("GAINED associations (from re-curation / new genes)", gained)
