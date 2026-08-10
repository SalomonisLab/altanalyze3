import sys, time, random, importlib.util
import pandas as pd
REPO="/Users/saljh8/Documents/GitHub/altanalyze3"; BENCH=f"{REPO}/dev_benchmarks/junction_annotation"
EXON="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
sys.path.insert(0,REPO)
import altanalyze3.components.long_read.gff_process as gp
from altanalyze3.components.aggregate.annotate import annotate_junctions as new_a
spec=importlib.util.spec_from_file_location("fz",f"{BENCH}/annotate_reference_frozen.py")
fz=importlib.util.module_from_spec(spec); spec.loader.exec_module(fz)
class A:
    def __init__(s,v): s.var=v
ec,gd,sd=gp.importEnsemblGenes(EXON,include_introns=True)
from collections import defaultdict
don,acc=defaultdict(list),defaultdict(list)
for (c,p,st,site),v in ec.items(): (don if site==2 else acc)[(c,st)].append(p)
keys=[k for k in don if k in acc and len(don[k])>200]
rnd=random.Random(5)
def build(n,frac_novel):
    rows=[]
    for i in range(n):
        c,st=rnd.choice(keys)
        d=rnd.choice(don[(c,st)]); a=rnd.choice(acc[(c,st)])
        if i < n*frac_novel: d+=rnd.randint(1,40); a+=rnd.randint(1,40)
        rows.append((c,d,a,st))
    v=pd.DataFrame(rows,columns=["chr","start","end","strand"]); v.index=[str(i) for i in range(n)]
    return v
for n,frac in ((20_000,0.15),):
    var=build(n,frac)
    t0=time.time(); o=A(var.copy()); fz.annotate_junctions(o,EXON); t_old=time.time()-t0
    t0=time.time(); m=A(var.copy()); new_a(m,EXON,novel_gene_mode="legacy"); t_new=time.time()-t0
    ident=sum(1 for x,y in zip(o.var["annotation"],m.var["annotation"]) if x==y)
    LOAD=2.19
    print(f"n={n:,}  fully-novel fraction={frac:.0%}")
    print(f"  OLD total {t_old:8.2f} s   NEW total {t_new:8.2f} s   total speedup {t_old/t_new:6.1f}x")
    print(f"  loop only (minus {LOAD:.2f}s reference load): OLD {t_old-LOAD:8.2f} s  NEW {t_new-LOAD:6.2f} s"
          f"   loop speedup {(t_old-LOAD)/max(t_new-LOAD,1e-6):,.0f}x")
    print(f"  identity: {ident:,} of {n:,} annotations identical  (mismatches {n-ident})")
    print(f"\n  projection to 500,000 junctions at {frac:.0%} novel:")
    print(f"    OLD {(t_old-LOAD)*500_000/n/60:8.1f} min loop + {LOAD:.1f}s load")
    print(f"    NEW {(t_new-LOAD)*500_000/n:8.1f} s   loop + {LOAD:.1f}s load")
