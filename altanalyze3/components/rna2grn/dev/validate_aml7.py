#!/usr/bin/env python3.11
"""Rebuild bundles (corrected weighting) and validate on held-out AML-7.

AML-7 was excluded from the demo reference; we predict its 11 cell-state GRNs
from RNA and compare to the true AML-7 GRN, against the cell-state-mean baseline.
"""
import sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn import training, Rna2GrnBundle

DATA = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
REPO = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2grn"
npz = f"{DATA}/matched/dataset.npz"

# rebuild with corrected cosine weighting
training.build_bundle(npz, f"{REPO}/rna2grn_bundle.pkl.gz", k=15, metric="cosine",
                      include_controls=True, created_at="2026-06-13")
training.build_bundle(npz, f"{DATA}/artifacts/rna2grn_bundle_noAML7.pkl.gz", k=15, metric="cosine",
                      exclude_samples=["AML-7_CCHMC"], include_controls=True, created_at="2026-06-13")
# also a no-controls, no-AML7 bundle to mirror the best benchmark config
training.build_bundle(npz, f"{DATA}/artifacts/rna2grn_bundle_noAML7_noctrl.pkl.gz", k=15, metric="cosine",
                      exclude_samples=["AML-7_CCHMC"], include_controls=False, created_at="2026-06-13")

d = np.load(npz, allow_pickle=True)
genes = d["genes"].astype(str); X = d["X"]; Y = d["Y"]
sample = d["sample"].astype(str); cs = d["cell_state"].astype(str)
group = d["group"].astype(str); heldout = d["heldout_demo"]

a7 = heldout  # AML-7 rows
Xa, Ya, csa = X[a7], Y[a7], cs[a7]
print(f"AML-7 held-out pseudobulks: {a7.sum()}  cell states: {list(csa)}")


def col_pearson(A, B):
    Az = A - A.mean(0, keepdims=True); Bz = B - B.mean(0, keepdims=True)
    num = (Az * Bz).sum(0); den = np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    out = np.full(A.shape[1], np.nan); ok = den > 1e-12; out[ok] = num[ok] / den[ok]; return out


def row_pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / (np.sqrt((a**2).sum() * (b**2).sum()) + 1e-12))


# baselines from the FULL reference excluding AML-7
ref = ~a7
csmean = {c: Y[ref & (cs == c)].mean(0) for c in pd.unique(cs[ref])}
gm = Y[ref].mean(0)
CSb = np.vstack([csmean.get(c, gm) for c in csa])
GMb = np.repeat(gm[None, :], a7.sum(), 0)

for label, bundle_path in [("knn(+ctrl)", f"{DATA}/artifacts/rna2grn_bundle_noAML7.pkl.gz"),
                           ("knn(-ctrl)", f"{DATA}/artifacts/rna2grn_bundle_noAML7_noctrl.pkl.gz")]:
    b = Rna2GrnBundle.load(bundle_path)
    df = pd.DataFrame(Xa, index=[f"AML7|{c}" for c in csa], columns=genes)
    t0 = time.time()
    res = b.predict_from_dataframe(df, normalized=True)
    dt = time.time() - t0
    P = res.predictions.values
    # metrics
    def ss(a, b): return ((a - b) ** 2).sum()
    skill_cs = 1 - ss(Ya, P) / ss(Ya, CSb)
    r2 = 1 - ss(Ya, P) / ss(Ya, GMb)
    prof = np.array([row_pearson(Ya[i], P[i]) for i in range(len(Ya))])
    prof_cs = np.array([row_pearson(Ya[i], CSb[i]) for i in range(len(Ya))])
    edge = np.nanmedian(col_pearson(Ya, P))
    # within-cellstate residual
    Yr = Ya - CSb; Pr = P - CSb
    wedge = np.nanmedian(col_pearson(Yr, Pr))
    print(f"\n=== {label} held-out AML-7 ({dt*1000:.0f} ms for {a7.sum()} pseudobulks) ===")
    print(f"  profile r (model)        : median={np.median(prof):.3f}")
    print(f"  profile r (cellstate mean): median={np.median(prof_cs):.3f}  [baseline]")
    print(f"  R2 vs global mean        : {r2:+.3f}")
    print(f"  skill vs cellstate-mean  : {skill_cs:+.3f}   (>0 beats baseline)")
    print(f"  per-edge across-cs r     : {edge:.3f}")
    print(f"  within-cellstate edge r  : {wedge:.3f}")

import os
print("\nbundle sizes:")
for f in [f"{REPO}/rna2grn_bundle.pkl.gz"]:
    print("  %.2f MB  %s" % (os.path.getsize(f) / 1e6, os.path.basename(f)))
