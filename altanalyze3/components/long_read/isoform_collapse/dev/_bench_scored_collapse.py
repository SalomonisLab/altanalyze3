import sys, os, time, collections, resource, subprocess, importlib.util
spec=importlib.util.spec_from_file_location("scored_collapse",
  os.path.join(os.path.dirname(__file__),"scored_collapse.py"))
sc=importlib.util.module_from_spec(spec); spec.loader.exec_module(sc)
from altanalyze3.components.long_read.isoform_collapse_utils import structure_tokens_for_containment

def peak_mb():
    kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb/(1024*1024) if sys.platform=='darwin' else kb/1024
def cur_mb():
    o=subprocess.check_output(["ps","-o","rss=","-p",str(os.getpid())],text=True).strip()
    return float(o)/1024 if o else 0

SAMPLE=sys.argv[1]
PATH=f"/Users/saljh8/Dropbox/Revio/BAMs/iPSC/{SAMPLE}/gff-output/transcript_associations.txt"
OUTDIR=os.path.dirname(PATH)

# ---- Step 1: LOAD (parse transcript_associations -> gene->{structure:reads}, first molecule) ----
t=time.time()
gene_struct=collections.defaultdict(collections.Counter); gene_mol=collections.defaultdict(dict)
nreads=0
with open(PATH) as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        if len(p)<5: continue
        g,strand,struct,mol,src=p[:5]
        gene_struct[g][struct]+=1
        d=gene_mol[g]
        if struct not in d: d[struct]=mol
        nreads+=1
t_load=time.time()-t; rss_load=cur_mb()
n_genes=len(gene_struct)
n_unique=sum(len(v) for v in gene_struct.values())

# ---- Step 2: TOKENIZE (structure -> token tuple, per gene) ----
t=time.time()
gene_tok={}
for g,sr in gene_struct.items():
    gene_tok[g]={s:tuple(structure_tokens_for_containment(s,g)) for s in sr}
t_tok=time.time()-t; rss_tok=cur_mb()

# ---- Step 3: COLLAPSE (per gene) ----
t=time.time()
n_final=0; n_struct_ge2=0; rows=[]
for g,sr in gene_struct.items():
    res=sc.collapse_gene(sr, gene_tok[g], min_reads=1)
    reps=set(res['long_reps'])|set(res['other_reps'])
    n_final+=len(reps)
    n_struct_ge2+=len(sr)
    longset=set(res['long_reps'])
    incoming=collections.defaultdict(list)
    for c,par in res['assignment'].items(): incoming[par].append(c)
    for r in reps:
        subs=incoming.get(r,[]); total=sr[r]+sum(sr[s] for s in subs)
        rows.append((g, gene_mol[g][r], res['blocks'][r], f"{res['score'][r]:.2f}", sr[r], len(subs), total, "long" if r in longset else "other"))
t_collapse=time.time()-t; rss_collapse=cur_mb()

OUT=f"{OUTDIR}/{SAMPLE}_scored_collapse_summary.tsv"
with open(OUT,"w") as o:
    o.write("gene\texemplar_id\texon_blocks\tscore\toriginal_reads\tn_collapsed_in\ttotal_reads_after\tbin\n")
    for row in rows: o.write("\t".join(map(str,row))+"\n")

print(f"========== {SAMPLE} ==========")
print(f"  total reads:             {nreads:,}")
print(f"  genes:                   {n_genes:,}")
print(f"  unique isoforms (>=1rd): {n_unique:,}")
print(f"  structures (kept, all): {n_struct_ge2:,}")
print(f"  FINAL collapsed isoforms:{n_final:,}")
print(f"  --- time per step ---")
print(f"  load:     {t_load:7.1f}s   (RSS {rss_load:6.0f}MB)")
print(f"  tokenize: {t_tok:7.1f}s   (RSS {rss_tok:6.0f}MB)")
print(f"  collapse: {t_collapse:7.1f}s   (RSS {rss_collapse:6.0f}MB)")
print(f"  TOTAL:    {t_load+t_tok+t_collapse:7.1f}s   PEAK RSS {peak_mb():.0f}MB")
print(f"  wrote {OUT} ({len(rows):,} rows)")
