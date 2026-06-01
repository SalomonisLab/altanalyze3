import sys, os, time, collections, resource, importlib.util
import multiprocessing as mp
spec=importlib.util.spec_from_file_location("scored_collapse",
  os.path.join(os.path.dirname(__file__),"scored_collapse.py"))
sc=importlib.util.module_from_spec(spec); spec.loader.exec_module(sc)
from altanalyze3.components.long_read.isoform_collapse_utils import structure_tokens_for_containment

REF='/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'
SAMPLES={'Ctrl':'/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt',
         'Aza':'/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/gff-output/transcript_associations.txt'}

def self_peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)  # darwin: bytes

def load_gene_chrom():
    gc={}
    with open(REF) as f:
        next(f)
        for line in f:
            t=line.rstrip("\n").split("\t")
            if len(t)>=3 and t[0] not in gc: gc[t[0]]=t[2]
    return gc

def worker(args):
    chrom, genes_for_chrom = args
    # each worker re-reads only the rows for its genes (low memory: holds one chromosome)
    gene_struct=collections.defaultdict(collections.Counter); gene_mol=collections.defaultdict(dict)
    gset=set(genes_for_chrom)
    for SAMPLE,PATH in SAMPLES.items():
        with open(PATH) as f:
            for line in f:
                i=line.find("\t")
                if i<0: continue
                g=line[:i]
                if g not in gset: continue
                p=line.rstrip("\n").split("\t")
                if len(p)<5: continue
                struct,mol=p[2],p[3]
                gene_struct[g][struct]+=1
                d=gene_mol[g]
                if struct not in d: d[struct]=f"{SAMPLE}:{mol}"
    n_final=0; n_unique=0; rows=[]
    for g,sr in gene_struct.items():
        n_unique+=len(sr)
        tok={s:tuple(structure_tokens_for_containment(s,g)) for s in sr}
        res=sc.collapse_gene(sr,tok,min_reads=1)
        reps=set(res['long_reps'])|set(res['other_reps'])
        n_final+=len(reps)
    return chrom, len(gset), n_unique, n_final, self_peak_mb(), len(gene_struct)

if __name__=='__main__':
    nproc=int(sys.argv[1]) if len(sys.argv)>1 else 8
    t0=time.time()
    gc=load_gene_chrom()
    # which genes appear in samples
    sgenes=set()
    for PATH in SAMPLES.values():
        with open(PATH) as f:
            for line in f:
                i=line.find("\t")
                if i>0: sgenes.add(line[:i])
    # partition genes by chromosome; unmapped -> 'unplaced'
    chrom_genes=collections.defaultdict(list)
    for g in sgenes:
        chrom_genes[gc.get(g,'unplaced')].append(g)
    parts=sorted(chrom_genes.items(), key=lambda kv:-len(kv[1]))  # biggest first for load balance
    print(f"partitions: {len(parts)} (by chromosome); top sizes: {[(c,len(gs)) for c,gs in parts[:6]]}", flush=True)
    t_setup=time.time()-t0
    t1=time.time()
    with mp.Pool(nproc) as pool:
        results=pool.map(worker, parts)
    t_run=time.time()-t1
    results.sort(key=lambda r:-r[2])
    tot_unique=sum(r[2] for r in results); tot_final=sum(r[3] for r in results)
    max_worker_rss=max(r[4] for r in results)
    print(f"\n========== PARALLEL by chromosome (nproc={nproc}) ==========")
    print(f"  partitions(chromosomes): {len(parts)}")
    print(f"  total unique isoforms:   {tot_unique:,}")
    print(f"  total FINAL isoforms:    {tot_final:,}")
    print(f"  setup(partition) time:   {t_setup:.1f}s")
    print(f"  parallel run wall-clock: {t_run:.1f}s")
    print(f"  TOTAL wall-clock:        {t_setup+t_run:.1f}s")
    print(f"  MAX per-worker peak RSS: {max_worker_rss:.0f}MB   <-- memory gain vs 1424MB single-process")
    print(f"\n  per-chromosome (chrom, genes, unique, final, peakRSS_MB):")
    for chrom,ng,nu,nf,rss,ngs in results[:8]:
        print(f"    {chrom:18s} genes={ng:5d} unique={nu:7,} final={nf:7,} rss={rss:5.0f}MB")
