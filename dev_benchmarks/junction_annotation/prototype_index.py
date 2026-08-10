"""Prototype: replace the 68,008-gene linear scan with a precomputed sweep index.

Contract preserved exactly: return the FIRST gene in geneData insertion order for which
    geneData[g][0][0] <= pos <= geneData[g][-1][1] and strandData[g] == strand
Chromosome is NOT tested, matching the current code.
"""
import sys, time, random, bisect, heapq
sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
import altanalyze3.components.long_read.gff_process as gp
EXON="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
ec,gd,sd = gp.importEnsemblGenes(EXON, include_introns=True)

def brute(pos,strand):
    for gene in gd:
        if gd[gene][0][0] <= pos <= gd[gene][-1][1] and sd[gene]==strand:
            return gene
    return None

def build_index(geneData, strandData):
    per_strand={}
    for rank,g in enumerate(geneData):
        lo=geneData[g][0][0]; hi=geneData[g][-1][1]
        if lo>hi:                      # empty interval, can never contain a position
            continue
        per_strand.setdefault(strandData[g],[]).append((lo,hi,rank,g))
    index={}
    for strand,ivs in per_strand.items():
        events=[]
        for lo,hi,rank,g in ivs:
            events.append((lo,0,rank,g)); events.append((hi+1,1,rank,g))
        events.sort(key=lambda e:(e[0],e[1]))
        starts=[]; winners=[]
        heap=[]; alive={}
        i=0
        while i<len(events):
            x=events[i][0]
            while i<len(events) and events[i][0]==x:
                _,kind,rank,g=events[i]
                if kind==0: heapq.heappush(heap,(rank,g)); alive[rank]=alive.get(rank,0)+1
                else: alive[rank]=alive.get(rank,0)-1
                i+=1
            while heap and alive.get(heap[0][0],0)<=0: heapq.heappop(heap)
            starts.append(x); winners.append(heap[0][1] if heap else None)
        index[strand]=(starts,winners)
    return index

t0=time.time(); idx=build_index(gd,sd); t_build=time.time()-t0
print(f"index build: {t_build:.2f} s   segments: " + ", ".join(f"{s}={len(idx[s][0]):,}" for s in sorted(idx)))

def fast(pos,strand):
    e=idx.get(strand)
    if e is None: return None
    starts,winners=e
    j=bisect.bisect_right(starts,pos)-1
    return winners[j] if j>=0 else None

# --- identity check: brute vs fast on a broad position set ---
rnd=random.Random(11)
probes=[]
for g in rnd.sample(list(gd),400):                      # real gene boundaries and just outside them
    lo=gd[g][0][0]; hi=gd[g][-1][1]
    for p in (lo-1,lo,lo+1,hi-1,hi,hi+1,(lo+hi)//2):
        for s in ('+','-'): probes.append((p,s))
for _ in range(600):
    probes.append((rnd.randint(1,250_000_000), rnd.choice(['+','-'])))
mismatch=[(p,s,brute(p,s),fast(p,s)) for p,s in probes if brute(p,s)!=fast(p,s)]
print(f"identity check: {len(probes):,} probes (gene boundaries +/-1, midpoints, random), mismatches = {len(mismatch)}")
if mismatch: print("  first 5:",mismatch[:5])

# --- speed ---
qs=[(rnd.randint(1,250_000_000), rnd.choice(['+','-'])) for _ in range(200)]
t0=time.time(); [brute(p,s) for p,s in qs]; t_b=(time.time()-t0)/200
qs2=[(rnd.randint(1,250_000_000), rnd.choice(['+','-'])) for _ in range(200_000)]
t0=time.time(); [fast(p,s) for p,s in qs2]; t_f=(time.time()-t0)/200_000
print(f"\nper call:  brute={1e3*t_b:8.3f} ms   indexed={1e6*t_f:8.2f} us   speedup={t_b/t_f:,.0f}x")
