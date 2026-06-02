#!/usr/bin/env python3
"""Paired statistical test of junction-correlation: WTA vs EM, using BOTH Pearson and Spearman.
Reuses the already-written per-cell isoform h5ads (WM34_CD34-isoform-{wta,em}.h5ad) and the junction
h5ad. Only junctions WITHIN an isoform's structure (adjacent exon tokens) and with >=200 reads are
used; isoforms need >=200 reads. For each isoform take the best within-structure junction correlation;
compare WTA vs EM paired over shared isoforms with Wilcoxon signed-rank + paired t-test, for both
correlation types."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr, pearsonr, wilcoxon, ttest_rel
B="/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM"
JUNC=f"{B}/WM34/WM34_CD34-junction.h5ad"; MIN=200
# ex_struct maps (from the kept structure->exemplar of each method); recompute via FINAL maps written by stage3? 
# Simpler: re-derive each isoform's structure from the catalog structure_to_exemplar is not saved per-method here,
# so use the junction-of-structure via the var name set per gene is not enough -> need ex_struct.
# We rebuild ex_struct by re-running stage1/stage2 quickly (cheap relative to full) OR load from saved.
import sys; sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.long_read.isoform_collapse import pipeline as P
SAMPLE=("WM34_CD34", f"{B}/WM34/gff-output/transcript_associations.txt"); ENST=f"{B}/gff-output/ENST_reference_structures.tsv"

def ex_struct_for(method):
    gr,_,_,_=P.stage1_collapse([SAMPLE],nproc=4,enst_cache=ENST,tier1_dir=f"{B}/analysis/em_vs_wta/tier1_{method}",collapse_method=method)
    _,kept,_=P.stage2_outliers(gr,min_total=3,return_soft=True)
    d={}
    for g,sm in kept.items():
        for struct,ex in sm.items():
            fid=f"{g}:{ex}"; 
            if fid not in d or len(struct)>len(d[fid]): d[fid]=struct
    return d

J=ad.read_h5ad(JUNC); Jx=J.X.tocsc() if sp.issparse(J.X) else sp.csc_matrix(J.X); jbc=J.obs_names.astype(str).values
jtot=np.asarray(Jx.sum(0)).ravel(); jname2col={}
for c,v in enumerate(J.var_names):
    if jtot[c]<MIN: continue
    g=v.split(":")[0]; ee=v.split(":",1)[1].split("=")[0]; jname2col.setdefault(g,{})[ee]=c

def juncs_of(s):
    t=s.split("|"); return {f"{t[k]}-{t[k+1]}" for k in range(len(t)-1)}

def best(method, ex_struct, corr):
    I=ad.read_h5ad(f"{B}/WM34/WM34_CD34-isoform-{method}.h5ad")
    Ix=I.X.tocsc() if sp.issparse(I.X) else sp.csc_matrix(I.X); ibc=I.obs_names.astype(str).values
    jset=set(jbc); common=[b for b in ibc if b in jset]
    ib={b:k for k,b in enumerate(ibc)}; jb={b:k for k,b in enumerate(jbc)}
    ci=[ib[b] for b in common]; cj=[jb[b] for b in common]
    itot=np.asarray(Ix.sum(0)).ravel(); out={}
    for col,v in enumerate(I.var_names):
        if itot[col]<MIN: continue
        g=v.split(":")[0]; st=ex_struct.get(v)
        if not st: continue
        js=[j for j in juncs_of(st) if j in jname2col.get(g,{})]
        if not js: continue
        ie=np.asarray(Ix[ci,col].todense()).ravel()
        if ie.std()==0: continue
        bestc=-2
        for j in js:
            je=np.asarray(Jx[cj,jname2col[g][j]].todense()).ravel()
            if je.std()==0: continue
            r=(pearsonr(ie,je)[0] if corr=="pearson" else spearmanr(ie,je).statistic)
            if r is not None and r>bestc: bestc=r
        if bestc>-2: out[v]=bestc
    return out

esw=ex_struct_for("wta"); ese=ex_struct_for("em")
for corr in ("pearson","spearman"):
    w=best("wta",esw,corr); e=best("em",ese,corr)
    shared=[k for k in w if k in e]
    wv=np.array([w[k] for k in shared]); ev=np.array([e[k] for k in shared])
    print(f"\n===== {corr.upper()} (paired n={len(shared)}) =====")
    print(f"  WTA mean={wv.mean():.4f} median={np.median(wv):.4f} | EM mean={ev.mean():.4f} median={np.median(ev):.4f}")
    print(f"  WTA better: {(wv>ev).sum()}  EM better: {(ev>wv).sum()}  tie: {(wv==ev).sum()}")
    try:
        print(f"  Wilcoxon signed-rank (EM-WTA) p={wilcoxon(ev,wv).pvalue:.3e}")
    except Exception as ex_: print("  wilcoxon:",ex_)
    print(f"  paired t-test (EM-WTA) p={ttest_rel(ev,wv).pvalue:.3e}  mean diff(EM-WTA)={ (ev-wv).mean():+.4f}")
    pd.DataFrame({"isoform":shared,"wta":wv,"em":ev}).to_csv(f"{B}/analysis/em_vs_wta/paired_{corr}.tsv",sep="\t",index=False)
print("\nPEARSON_TEST_DONE")
