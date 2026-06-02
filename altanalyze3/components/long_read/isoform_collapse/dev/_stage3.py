import sys, os, time, resource
import numpy as np, scipy.sparse as sp, anndata as ad, pandas as pd
def peak_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)

MAP='/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/FINAL_structure_to_exemplar.tsv'
CAT='/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/FINAL_isoform_catalog.tsv'
ITGA2B='ENSG00000005961'
SAMPLES=[('Ctrl','/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/EP_SRSF2_Ctrl.h5ad',
          '/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt'),
         ('Aza','/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/EP_SFSF2_AZA.h5ad',
          '/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/gff-output/transcript_associations.txt')]

# (gene,structure) -> final exemplar
s2e={}
with open(MAP) as f:
    next(f)
    for line in f:
        g,struct,ex=line.rstrip("\n").split("\t"); s2e[(g,struct)]=ex
itg_finals=set()
with open(CAT) as f:
    next(f)
    for line in f:
        g,ex,blk,tot,b=line.rstrip("\n").split("\t")
        if g==ITGA2B: itg_finals.add(ex)
print(f"loaded {len(s2e):,} surviving structures", flush=True)

for SAMPLE,H5,TA in SAMPLES:
    t=time.time()
    # OPTIMIZED: build var-id -> final exemplar only for surviving molecules (no 8.9M tuple dict)
    var2final={}
    with open(TA) as f:
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)<5: continue
            ex=s2e.get((p[0],p[2]))
            if ex is not None: var2final[f"{p[0]}:{p[3]}"]=ex
    a=ad.read_h5ad(H5)
    final_ids=sorted(set(var2final.values())); fidx={f:j for j,f in enumerate(final_ids)}
    rows=[]; cols=[]
    for i,v in enumerate(a.var_names):
        ex=var2final.get(v)
        if ex is not None: rows.append(i); cols.append(fidx[ex])
    G=sp.coo_matrix((np.ones(len(rows),np.int64),(rows,cols)),shape=(a.shape[1],len(final_ids))).tocsr()
    Xf=(a.X.astype(np.int64))@G
    out=ad.AnnData(X=Xf.tocsr(), obs=a.obs.copy(), var=pd.DataFrame(index=final_ids))
    OUT=f"{os.path.dirname(H5)}/{os.path.basename(H5).split('.h5ad')[0]}-final_isoform.h5ad"
    out.write_h5ad(OUT, compression='gzip')
    raw=int(a.X.sum()); fin=int(Xf.sum())
    print(f"\n[{SAMPLE}] {a.shape[0]} cells x {a.shape[1]:,} mol -> {len(final_ids):,} final | {time.time()-t:.1f}s peakRSS {peak_mb():.0f}MB")
    print(f"   reads raw={raw:,} -> final={fin:,} ({100*fin/raw:.0f}% preserved)  wrote {os.path.basename(OUT)}")
    present=[f for f in final_ids if f in itg_finals]
    if present:
        cs=np.asarray(Xf.sum(axis=0)).ravel()
        print(f"   ITGA2B top finals: "+", ".join(f"{f}={int(cs[fidx[f]])}" for f in sorted(present,key=lambda x:-cs[fidx[x]])[:4]))
