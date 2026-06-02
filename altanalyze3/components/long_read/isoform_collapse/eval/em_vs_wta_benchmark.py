#!/usr/bin/env python3
"""Benchmark: EM vs winner-takes-all (WTA) isoform collapse, using splicing (junctions) as a
positive control, within a single sample (WM34_CD34).

Hypothesis: a correct isoform-abundance estimate should track the abundance of the splice junctions
that DEFINE that isoform, across cells. So, per cell, an isoform's expression should correlate with
its constituent junctions' expression. Whichever collapse method yields isoform-by-cell matrices
that correlate BETTER with their junctions assigns ambiguous (substring) reads more correctly.

KEY POINT (why the per-cell matrices differ): both methods assign child structures -> parent
isoform(s); each child structure maps uniquely to MOLECULE IDs in the source sample. We PROPAGATE
those molecule associations down to the per-cell matrix:
  - WTA: each child structure has ONE parent -> its molecules' reads go entirely to that parent.
  - EM : each child structure has FRACTIONAL weights over parents (em_weights) -> its molecules'
         reads split fractionally across parents in each cell.
(The first version wrongly re-keyed EM with the hard assignment, so the matrices were identical.)
"""
import os, sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.long_read.isoform_collapse import pipeline as P

BASE = "/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM"
SAMPLE = ("WM34_CD34", f"{BASE}/WM34/gff-output/transcript_associations.txt")
MOL_H5AD = f"{BASE}/WM34/WM34_CD34.h5ad"
JUNC_H5AD = f"{BASE}/WM34/WM34_CD34-junction.h5ad"
ENST = f"{BASE}/gff-output/ENST_reference_structures.tsv"
OUTDIR = f"{BASE}/analysis/em_vs_wta"; os.makedirs(OUTDIR, exist_ok=True)
MIN_READS = 200


def collapse(method):
    gr, _, _, _ = P.stage1_collapse([SAMPLE], nproc=4, enst_cache=ENST,
                                    tier1_dir=f"{OUTDIR}/tier1", collapse_method=method)
    catalog, kept = P.stage2_outliers(gr, min_total=3)
    return gr, catalog, kept


def fractional_var2final(ta_path, kept, soft):
    """molecule var ('gene:mol') -> {final_iso_id: weight}. Propagates the structure->parent
    association (and EM weights) DOWN to each molecule of that structure in the sample.
    WTA: one parent, weight 1. EM: em_weights of the molecule's structure across parents.

    CRITICAL: both methods are restricted to the SAME post-stage2 kept exemplar set (so the EM
    fractional weights to any parent dropped by the <min_total outlier filter are removed and the
    remaining weights renormalized). Otherwise EM would silently bypass stage-2 filtering and emit
    every pre-filter structure as its own column (the bug that gave 1.18M EM 'isoforms', 100% reads).
    """
    # only-kept exemplars (post-stage2), as full "gene:exemplar" ids
    kept_final = {f"{g}:{ex}" for g, m in kept.items() for s, ex in m.items()}
    s2e = {(g, s): f"{g}:{ex}" for g, m in kept.items() for s, ex in m.items()}
    s2soft = {}
    if soft:
        for g, sm in soft.items():
            for s, wmap in sm.items():
                # keep only weights to surviving (kept) exemplars, then renormalize
                kw = {f"{g}:{ex}": w for ex, w in wmap.items() if f"{g}:{ex}" in kept_final}
                tot = sum(kw.values())
                if tot > 0:
                    s2soft[(g, s)] = {k: v / tot for k, v in kw.items()}
    var2final = {}
    with open(ta_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            g, struct, mol = p[0], p[2], p[3]
            if soft and (g, struct) in s2soft:
                var2final[f"{g}:{mol}"] = s2soft[(g, struct)]
            elif (g, struct) in s2e:
                var2final[f"{g}:{mol}"] = {s2e[(g, struct)]: 1.0}
    return var2final


def build_matrix(var2final):
    a = ad.read_h5ad(MOL_H5AD)
    final_ids = sorted({fid for wm in var2final.values() for fid in wm})
    fidx = {f: j for j, f in enumerate(final_ids)}
    rows, cols, data = [], [], []
    for i, v in enumerate(a.var_names.astype(str).values):
        wm = var2final.get(v)
        if not wm:
            continue
        for fid, w in wm.items():
            rows.append(i); cols.append(fidx[fid]); data.append(w)
    G = sp.coo_matrix((data, (rows, cols)), shape=(a.shape[1], len(final_ids))).tocsr()
    X = (a.X.astype(np.float64)) @ G
    return ad.AnnData(X=X.tocsr(), obs=a.obs.copy(), var=pd.DataFrame(index=final_ids))


def exemplar_structure(kept):
    ex_struct = {}
    for g, sm in kept.items():
        for struct, ex in sm.items():
            fid = f"{g}:{ex}"
            cur = ex_struct.get(fid)
            if cur is None or len(struct) > len(cur):
                ex_struct[fid] = struct
    return ex_struct


def run(method):
    """Collapse + per-cell re-key via the PACKAGE stage3 (structure-based, kept_soft for EM).
    Returns (isoform_adata, ex_struct map). Writes the alternative h5ad."""
    gr, _, _, _ = P.stage1_collapse([SAMPLE], nproc=4, enst_cache=ENST,
                                    tier1_dir=f"{OUTDIR}/tier1_{method}", collapse_method=method)
    cat, kept, ksoft = P.stage2_outliers(gr, min_total=3, return_soft=True)
    odir = f"{OUTDIR}/{method}"; os.makedirs(odir, exist_ok=True)
    out_path, n_final, raw, final, _, _ = P.stage3_rekey_h5ad(
        SAMPLE[0], MOL_H5AD, SAMPLE[1], kept, outdir=odir,
        kept_soft=(ksoft if method == "em" else None))
    # copy to the comparison filename next to the sample
    dest = f"{BASE}/WM34/WM34_CD34-isoform-{method}.h5ad"
    import shutil; shutil.copy(out_path, dest)
    print(f"[{method}] catalog={len(cat):,} isoforms={n_final:,} reads {raw:,}->{int(final):,} "
          f"({100*final/raw:.1f}%) -> {dest}")
    return ad.read_h5ad(dest), exemplar_structure(kept)


print("=== WTA ===")
I_w, exs_w = run("wta")
print("=== EM ===")
I_e, exs_e = run("em")

# CONSTRAINT: WTA and EM must produce the SAME isoform set (counts differ, not the feature set).
set_w, set_e = set(I_w.var_names), set(I_e.var_names)
assert set_w == set_e, (f"WTA/EM isoform sets DIFFER: |WTA|={len(set_w)} |EM|={len(set_e)} "
                        f"only_wta={len(set_w-set_e)} only_em={len(set_e-set_w)}")
print(f"[validate] WTA & EM share the SAME isoform set: {len(set_w):,} isoforms")

J = ad.read_h5ad(JUNC_H5AD)
Jx = J.X.tocsc() if sp.issparse(J.X) else sp.csc_matrix(J.X)
jbc = J.obs_names.astype(str).values
jtot = np.asarray(Jx.sum(0)).ravel()
jname2col = {}
for c, v in enumerate(J.var_names):
    if jtot[c] < MIN_READS:
        continue
    g = v.split(":")[0]; ee = v.split(":", 1)[1].split("=")[0]
    jname2col.setdefault(g, {})[ee] = c
print(f"junctions >= {MIN_READS} reads: {sum(len(d) for d in jname2col.values()):,}")


def junctions_of(struct):
    """The splice junctions an isoform's structure ACTUALLY contains: each ADJACENT pair of exon
    tokens in the structure (Ea immediately followed by Eb) is the exon-exon boundary that isoform
    splices. A junction is only valid for an isoform if both its exon tokens are consecutive here --
    i.e. the isoform structurally uses that junction. This enforces a biologically valid match (the
    junction is part of the isoform) rather than an ambiguous gene-level co-occurrence."""
    t = struct.split("|")
    return {f"{t[k]}-{t[k+1]}" for k in range(len(t) - 1)}


def best_corr(I, ex_struct, label):
    Ix = I.X.tocsc() if sp.issparse(I.X) else sp.csc_matrix(I.X)
    ibc = I.obs_names.astype(str).values
    jset = set(jbc); common = [b for b in ibc if b in jset]
    ib = {b: k for k, b in enumerate(ibc)}; jb = {b: k for k, b in enumerate(jbc)}
    ci = [ib[b] for b in common]; cj = [jb[b] for b in common]
    itot = np.asarray(Ix.sum(0)).ravel()
    rows = []
    for col, v in enumerate(I.var_names):
        if itot[col] < MIN_READS:
            continue
        g = v.split(":")[0]; struct = ex_struct.get(v)
        if not struct:
            continue
        # REQUIRED: only junctions that are within this isoform's structure (adjacent exon tokens)
        # AND are real quantified junctions (>=MIN_READS). No gene-level / ambiguous junction matches.
        struct_juncs = junctions_of(struct)
        juncs = [jj for jj in struct_juncs if jj in jname2col.get(g, {})]
        if not juncs:
            continue
        ie = np.asarray(Ix[ci, col].todense()).ravel()
        if ie.std() == 0:
            continue
        best = -2.0; bj = None
        for jj in juncs:
            je = np.asarray(Jx[cj, jname2col[g][jj]].todense()).ravel()
            if je.std() == 0:
                continue
            rho = spearmanr(ie, je).statistic
            if rho is not None and rho > best:
                best = rho; bj = jj
        if bj is not None:
            rows.append((v, g, best, bj, float(itot[col])))
    df = pd.DataFrame(rows, columns=["isoform", "gene", "best_junction_corr", "best_junction", "iso_reads"])
    df.to_csv(f"{OUTDIR}/best_corr_{label}.tsv", sep="\t", index=False)
    print(f"[{label}] isoforms scored: {len(df)}  mean={df['best_junction_corr'].mean():.4f}  "
          f"median={df['best_junction_corr'].median():.4f}")
    return df


wta_df = best_corr(I_w, exs_w, "wta")
em_df = best_corr(I_e, exs_e, "em")
m = wta_df.merge(em_df, on="isoform", suffixes=("_wta", "_em"))
print(f"\n=== PAIRED (n={len(m)} shared isoforms) ===")
print(f"  mean   WTA={m['best_junction_corr_wta'].mean():.4f}  EM={m['best_junction_corr_em'].mean():.4f}")
print(f"  median WTA={m['best_junction_corr_wta'].median():.4f}  EM={m['best_junction_corr_em'].median():.4f}")
em_better = int((m['best_junction_corr_em'] > m['best_junction_corr_wta']).sum())
wta_better = int((m['best_junction_corr_wta'] > m['best_junction_corr_em']).sum())
print(f"  EM better: {em_better}   WTA better: {wta_better}   tie: {len(m)-em_better-wta_better}")
try:
    diff = (m['best_junction_corr_em'] - m['best_junction_corr_wta'])
    if (diff != 0).any():
        print(f"  Wilcoxon paired p={wilcoxon(m['best_junction_corr_em'], m['best_junction_corr_wta']).pvalue:.3e}")
except Exception as e:
    print("  wilcoxon:", e)
m.to_csv(f"{OUTDIR}/paired_em_vs_wta.tsv", sep="\t", index=False)
print("BENCHMARK_DONE")
