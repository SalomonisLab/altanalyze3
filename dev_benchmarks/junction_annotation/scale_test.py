import sys, time, random
import pandas as pd, numpy as np, scipy.sparse as sp
sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
import altanalyze3.components.long_read.gff_process as gp
EXON="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
ec,gd,sd = gp.importEnsemblGenes(EXON, include_introns=True)

# isolated cost of one find_completely_novel_annot call (the whole-gene scan)
def fcna(pos,strand):
    for gene in gd:
        if gd[gene][0][0] <= pos <= gd[gene][-1][1] and sd[gene]==strand:
            return gene, gp.findNovelSpliceSite(gene,pos,strand)
    return None,None
rnd=random.Random(3)
for strand,label in (('+','hit-possible (+)'),('-','always-miss (-)')):
    ps=[rnd.randint(1_000_000,200_000_000) for _ in range(60)]
    t0=time.time(); [fcna(p,strand) for p in ps]; dt=(time.time()-t0)/60
    print(f"find_completely_novel_annot  {label:<18} {1e3*dt:7.2f} ms/call  ({len(gd):,} genes scanned)")

# iterrows vs itertuples on a var-shaped frame
n=200_000
var=pd.DataFrame({"chr":["chr1"]*n,"start":np.arange(n),"end":np.arange(n)+100,"strand":["+"]*n})
t0=time.time()
for i,(idx,row) in enumerate(var.iterrows()):
    _=row["chr"],row["start"],row["end"],row["strand"]
t_ir=time.time()-t0
t0=time.time()
for r in zip(var["chr"].to_numpy(),var["start"].to_numpy(),var["end"].to_numpy(),var["strand"].to_numpy()):
    _=r
t_np=time.time()-t0
print(f"\nrow iteration over {n:,} rows: iterrows={t_ir:6.2f} s   numpy zip={t_np:6.3f} s   ratio={t_ir/max(t_np,1e-9):.0f}x")

# export_dense_matrix cost model: toarray + to_csv
for nj,ns in ((200_000,20),):
    X=sp.random(ns,nj,density=0.15,format="csr",dtype=np.int32)
    t0=time.time(); D=X.toarray().T; t_d=time.time()-t0
    df=pd.DataFrame(D,columns=[f"s{i}" for i in range(ns)])
    df.insert(0,"UID",[f"g:E1.1-E2.1=chr1:{i}-{i+100}" for i in range(nj)])
    t0=time.time(); df.to_csv("/dev/null",sep="\t",index=False); t_c=time.time()-t0
    print(f"export_dense_matrix  {nj:,} junctions x {ns} samples: toarray={t_d:5.2f} s  to_csv={t_c:6.2f} s  dense RAM={D.nbytes/1e9:.2f} GB")
