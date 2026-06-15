#!/usr/bin/env python3.11
"""How similar are Multiome vs TEA-seq vs AML(3') GRNs?

Directly compares GRN connection-score profiles across protocols using:
  (a) Multiome_WT vs RC2_TEA on shared cell states (control-vs-control, same biology)
  (b) within-protocol cell-state-to-cell-state correlation baseline
This answers whether a model trained on one protocol can transfer to another.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

GRN = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/tf_to_gene_connection_scores_log10-NOT_ordered_clusters_COMBINED_AML_Multiome_RC2.txt"
OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn/reports"

df = pd.read_csv(GRN, sep="\t")
cols = df.columns[2:].tolist()
V = df[cols].to_numpy(float)
col_idx = {c: i for i, c in enumerate(cols)}

# parse columns into (group, cellstate) using the same cell-state vocab approach
import h5py
f = h5py.File("/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/pseudobulk_counts_hashed.h5ad", "r")
cs_vocab = [x.decode() if isinstance(x, bytes) else x for x in f['obs']['Hs-BM-titrated-reference-centroid']['categories'][:]]
f.close()
cs_us = sorted(((c.replace(" ", "_"), c) for c in cs_vocab), key=lambda x: -len(x[0]))


def parse(col):
    for us, orig in cs_us:
        if col == us or col.endswith("_" + us):
            return (col[:len(col) - len(us) - 1] if col != us else "", orig)
    return (None, None)


parsed = {c: parse(c) for c in cols}

def group_of(stok):
    if stok == "Multiome_WT": return "Multiome"
    if stok == "RC2_TEA": return "TEA"
    if stok.startswith("AML_"): return "AML3p"
    return "patient3p"

prof = {}  # (group, cs) -> column vector
for c, (stok, cs) in parsed.items():
    prof.setdefault((group_of(stok), cs), []).append(V[:, col_idx[c]])
prof = {k: np.mean(np.vstack(v), 0) for k, v in prof.items()}

def pcorr(a, b):
    m = (a > 0) | (b > 0)
    if m.sum() < 10: return np.nan
    return np.corrcoef(a, b)[0, 1]

# (a) Multiome vs TEA on shared cell states
mt_cs = sorted({cs for (g, cs) in prof if g == "Multiome"} & {cs for (g, cs) in prof if g == "TEA"})
print("=== Multiome_WT vs RC2_TEA on shared cell states ===")
rows = []
for cs in mt_cs:
    a, b = prof[("Multiome", cs)], prof[("TEA", cs)]
    r = pcorr(a, b); rho = spearmanr(a, b).correlation
    rows.append((cs, r, rho))
    print(f"  {cs:24s} pearson={r:.3f} spearman={rho:.3f}")
mt = pd.DataFrame(rows, columns=["cell_state", "pearson", "spearman"])
print(f"  MEAN cross-protocol (same cell state): pearson={mt.pearson.mean():.3f} spearman={mt.spearman.mean():.3f}")

# (b) within-Multiome: same cell state has how much higher corr than different cell state?
mult_cs = [cs for (g, cs) in prof if g == "Multiome"]
import itertools
same = mt.pearson.mean()
diffs = []
mult = {cs: prof[("Multiome", cs)] for (g, cs) in prof if g == "Multiome"}
keys = list(mult)
for a, b in itertools.combinations(keys, 2):
    diffs.append(pcorr(mult[a], mult[b]))
diffs = np.array(diffs)
print(f"\n=== Reference: within-Multiome BETWEEN different cell states ===")
print(f"  mean pearson across different cell states: {np.nanmean(diffs):.3f} (median {np.nanmedian(diffs):.3f})")
print(f"\nINTERPRETATION:")
print(f"  cross-protocol same-cellstate r = {same:.3f}")
print(f"  cross-cellstate within-Multiome r = {np.nanmean(diffs):.3f}")
print(f"  -> if cross-protocol r >> cross-cellstate r, protocols are reconcilable;")
print(f"     cell-state identity dominates the GRN profile over protocol.")
mt.to_csv(f"{OUT}/protocol_similarity_multiome_vs_tea.csv", index=False)
print(f"\nwrote {OUT}/protocol_similarity_multiome_vs_tea.csv")
