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

SAMPLES={'Ctrl':'/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt',
         'Aza':'/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/gff-output/transcript_associations.txt'}
OUTDIR='/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output'

# ---- Step 1: LOAD ALL SAMPLES COMBINED -> gene->{structure:pooled_reads}, first molecule (sample-tagged) ----
t=time.time()
gene_struct=collections.defaultdict(collections.Counter); gene_mol=collections.defaultdict(dict)
nreads=0
for SAMPLE,PATH in SAMPLES.items():
    with open(PATH) as f:
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)<5: continue
            g,strand,struct,mol,src=p[:5]
            gene_struct[g][struct]+=1            # POOLED across samples
            d=gene_mol[g]
            if struct not in d: d[struct]=f"{SAMPLE}:{mol}"
            nreads+=1
t_load=time.time()-t; rss_load=cur_mb()
n_genes=len(gene_struct); n_unique=sum(len(v) for v in gene_struct.values())

# ---- Step 2: TOKENIZE ----
t=time.time()
gene_tok={g:{s:tuple(structure_tokens_for_containment(s,g)) for s in sr} for g,sr in gene_struct.items()}
t_tok=time.time()-t; rss_tok=cur_mb()

# ---- Step 3: COLLAPSE (integrated, all structures, min_reads=1) ----
t=time.time()
n_final=0; rows=[]
for g,sr in gene_struct.items():
    res=sc.collapse_gene(sr, gene_tok[g], min_reads=1)
    reps=set(res['long_reps'])|set(res['other_reps']); longset=set(res['long_reps'])
    incoming=collections.defaultdict(list)
    for c,par in res['assignment'].items(): incoming[par].append(c)
    for r in reps:
        subs=incoming.get(r,[]); total=sr[r]+sum(sr[s] for s in subs)
        rows.append((g, gene_mol[g][r], res['blocks'][r], sr[r], len(subs), total, "long" if r in longset else "other"))
    n_final+=len(reps)
t_collapse=time.time()-t; rss_collapse=cur_mb()

# ---- Step 4: REMOVE OUTLIER ISOFORMS (cross-sample low-read filter on FINAL total reads) ----
t=time.time()
def kept(min_total): return sum(1 for r in rows if r[5]>=min_total)
filt={m:kept(m) for m in [1,2,3,5]}
t_filt=time.time()-t

OUT=f"{OUTDIR}/INTEGRATED_scored_collapse_summary.tsv"
with open(OUT,"w") as o:
    o.write("gene\texemplar_id\texon_blocks\toriginal_reads\tn_collapsed_in\ttotal_reads_after\tbin\n")
    for row in rows: o.write("\t".join(map(str,row))+"\n")

print("========== INTEGRATED (Ctrl + Aza combined) ==========")
print(f"  total reads (pooled):    {nreads:,}")
print(f"  genes:                   {n_genes:,}")
print(f"  unique isoforms (pooled):{n_unique:,}")
print(f"  FINAL collapsed isoforms:{n_final:,}")
print(f"  --- after OUTLIER removal (filter on final total reads across all samples) ---")
for m in [1,2,3,5]:
    print(f"    keep total_reads >= {m}: {filt[m]:,} isoforms")
print(f"  --- time per step ---")
print(f"  load:     {t_load:7.1f}s   (RSS {rss_load:6.0f}MB)")
print(f"  tokenize: {t_tok:7.1f}s   (RSS {rss_tok:6.0f}MB)")
print(f"  collapse: {t_collapse:7.1f}s   (RSS {rss_collapse:6.0f}MB)")
print(f"  filter:   {t_filt:7.1f}s")
print(f"  TOTAL:    {t_load+t_tok+t_collapse+t_filt:7.1f}s   PEAK RSS {peak_mb():.0f}MB")
print(f"  wrote {OUT} ({len(rows):,} rows)")
