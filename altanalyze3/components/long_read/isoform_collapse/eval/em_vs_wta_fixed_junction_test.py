#!/usr/bin/env python3
"""Paired Pearson test, WTA vs EM, with the SAME junction fixed per isoform.

For each isoform present in both methods (>=200 reads), select ONE constituent junction independently
of method -- the highest-read junction that is within the isoform's structure (adjacent exon tokens)
and has >=200 reads. Correlate (Pearson, per cell) that one fixed junction against BOTH the WTA and
the EM isoform values. Paired t-test of the two correlation vectors. This controls the junction so the
only thing differing between paired values is the collapse method's per-cell isoform estimate."""
import sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import pearsonr, ttest_rel, wilcoxon
sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.long_read.isoform_collapse import pipeline as P
B="/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM"; MIN=200
SAMPLE=("WM34_CD34", f"{B}/WM34/gff-output/transcript_associations.txt"); ENST=f"{B}/gff-output/ENST_reference_structures.tsv"

def ex_struct_for(method):
    gr,_,_,_=P.stage1_collapse([SAMPLE],nproc=4,enst_cache=ENST,tier1_dir=f"{B}/analysis/em_vs_wta/tier1_{method}",collapse_method=method)
    _,kept,_=P.stage2_outliers(gr,min_total=3,return_soft=True)
    d={}
    for g,sm in kept.items():
        for st,ex in sm.items():
            fid=f"{g}:{ex}"
            if fid not in d or len(st)>len(d[fid]): d[fid]=st
    return d

esw=ex_struct_for("wta"); ese=ex_struct_for("em")  # identical isoform set; use WTA's (same keys)
J=ad.read_h5ad(f"{B}/WM34/WM34_CD34-junction.h5ad"); Jx=J.X.tocsc() if sp.issparse(J.X) else sp.csc_matrix(J.X)
jbc=J.obs_names.astype(str).values; jtot=np.asarray(Jx.sum(0)).ravel()
jname2col={}
for c,v in enumerate(J.var_names):
    if jtot[c]<MIN: continue
    g=v.split(":")[0]; ee=v.split(":",1)[1].split("=")[0]; jname2col.setdefault(g,{})[ee]=(c,jtot[c])
def juncs_of(s):
    t=s.split("|"); return {f"{t[k]}-{t[k+1]}" for k in range(len(t)-1)}

Iw=ad.read_h5ad(f"{B}/WM34/WM34_CD34-isoform-wta.h5ad"); Ie=ad.read_h5ad(f"{B}/WM34/WM34_CD34-isoform-em.h5ad")
Iwx=Iw.X.tocsc() if sp.issparse(Iw.X) else sp.csc_matrix(Iw.X)
Iex=Ie.X.tocsc() if sp.issparse(Ie.X) else sp.csc_matrix(Ie.X)
wcol={v:k for k,v in enumerate(Iw.var_names)}; ecol={v:k for k,v in enumerate(Ie.var_names)}
ibc=Iw.obs_names.astype(str).values; jset=set(jbc)
common=[b for b in ibc if b in jset]
ib={b:k for k,b in enumerate(ibc)}; jb={b:k for k,b in enumerate(jbc)}
ci=[ib[b] for b in common]; cj=[jb[b] for b in common]
itw=np.asarray(Iwx.sum(0)).ravel()

rows=[]
for v in Iw.var_names:
    if v not in ecol: continue
    if itw[wcol[v]]<MIN: continue
    g=v.split(":")[0]; st=esw.get(v)
    if not st: continue
    cand=[(jname2col[g][jj][1], jj) for jj in juncs_of(st) if jj in jname2col.get(g,{})]
    if not cand: continue
    _, fixed_j = max(cand)              # SAME junction for both methods = highest-read constituent
    je=np.asarray(Jx[cj, jname2col[g][fixed_j][0]].todense()).ravel()
    if je.std()==0: continue
    iw=np.asarray(Iwx[ci,wcol[v]].todense()).ravel(); ie=np.asarray(Iex[ci,ecol[v]].todense()).ravel()
    if iw.std()==0 or ie.std()==0: continue
    rw=pearsonr(iw,je)[0]; re=pearsonr(ie,je)[0]
    rows.append((v,fixed_j,rw,re))
df=pd.DataFrame(rows,columns=["isoform","fixed_junction","pearson_wta","pearson_em"])
df.to_csv(f"{B}/analysis/em_vs_wta/fixed_junction_pearson.tsv",sep="\t",index=False)
w=df["pearson_wta"].values; e=df["pearson_em"].values
print(f"n isoforms (same junction both methods) = {len(df):,}")
print(f"  Pearson mean: WTA={w.mean():.4f}  EM={e.mean():.4f}   diff(EM-WTA)={ (e-w).mean():+.4f}")
print(f"  paired t-test p = {ttest_rel(e,w).pvalue:.3e}")
print("FIXED_DONE")
