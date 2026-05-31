import collections, statistics
from altanalyze3.components.long_read.isoform_collapse_utils import (
    is_contiguous_subsequence, structure_tokens_for_containment)
def nb(tk): return len({t.split(':')[-1].split('.')[0].split('_')[0] for t in tk if t.split(':')[-1].startswith('E')})
def contig(a,b): return len(a)<len(b) and is_contiguous_subsequence(a,b)
def prefix_or_suffix(a,b):  # a is a 5' prefix or 3' suffix of b (a shorter)
    n,m=len(a),len(b); return n<m and (b[:n]==a or b[m-n:]==a)

def bin_gene(counts, known=None, gene=""):
    known=known or {}; total=sum(counts.values()); Sset=list(counts)
    tok={s:structure_tokens_for_containment(s,gene) for s in Sset}
    tokset={s:set(tok[s]) for s in Sset}
    by_len=sorted(Sset,key=lambda s:len(tok[s]))

    # 1. CLUSTER by single-linkage containment
    parent={s:s for s in Sset}
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for i,a in enumerate(by_len):
        for b in by_len[i+1:]:
            if len(tok[b])==len(tok[a]): continue
            if tokset[a]<=tokset[b] and is_contiguous_subsequence(tok[a],tok[b]):
                ra,rb=find(a),find(b)
                if ra!=rb: parent[ra]=rb
    clusters=collections.defaultdict(list)
    for s in Sset: clusters[find(s)].append(s)

    long_iso=[];unique=[];ambiguous=[];orphan=[]
    for root,members in clusters.items():
        # 2. Identify LONG ISOFORMS, expression-aware so a RARE long extension never demotes an
        #    abundant shorter form. Build terminal-depth equivalence classes: s ~ t iff one is a
        #    prefix/suffix (5'/3' depth variant) of the other AND that pairing is UNAMBIGUOUS --
        #    the shorter is a prefix/suffix of no OTHER longer member (else it's a shared end of
        #    two distinct isoforms -> not a depth variant, must stay mappable as ambiguous).
        msort=sorted(members,key=lambda s:len(tok[s]))
        p2={s:s for s in members}
        def f2(x):
            while p2[x]!=x: p2[x]=p2[p2[x]]; x=p2[x]
            return x
        for i,a in enumerate(msort):
            longer_pref_suf=[b for b in members if len(tok[b])>len(tok[a]) and prefix_or_suffix(tok[a],tok[b])]
            if len(longer_pref_suf)==1:          # unambiguous extension -> same isoform
                ra,rb=f2(a),f2(longer_pref_suf[0])
                if ra!=rb: p2[ra]=rb
        classes=collections.defaultdict(list)
        for s in members: classes[f2(s)].append(s)
        cobjs=[]
        for r,mem in classes.items():
            rep=max(mem,key=lambda m:counts[m]); longest=max(mem,key=lambda m:len(tok[m]))
            cobjs.append({'rep':rep,'reads':sum(counts[m] for m in mem),'members':set(mem),
                'tk':tok[longest],'is_known':next((known[m] for m in mem if m in known),None),
                'blocks':nb(tok[longest])})
        # ABSORB rare end-extensions: a class C that strictly CONTAINS a more-abundant class D
        # (D is a contiguous subseq of C) and is itself much rarer than D is just D sequenced with
        # untrusted extra ends -> fold C into D (D keeps representative; C's reads + tk extent move in).
        # "much rarer" = C.reads < D.reads (any longer-but-less-abundant container folds to the
        # most-abundant class it extends). Iterate to a fixed point on the small class list.
        changed=True
        while changed:
            changed=False
            for c in list(cobjs):
                if c not in cobjs: continue
                # candidate hosts D: abundant classes that c strictly contains
                hosts=[d for d in cobjs if d is not c and contig(d['tk'],c['tk']) and d['reads']>c['reads']]
                if hosts:
                    d=max(hosts,key=lambda x:(x['reads'],len(x['tk'])))
                    d['members']|=c['members']; d['reads']+=c['reads']
                    if len(c['tk'])>len(d['tk']): d['tk']=c['tk']; d['blocks']=nb(c['tk'])
                    cobjs.remove(c); changed=True
        # LONG BIN = classes not a contiguous subseq of a >=-expressed class (after absorption).
        cluster_long=[c for c in cobjs
                      if not any(d is not c and contig(c['tk'],d['tk']) and d['reads']>=c['reads']
                                 for d in cobjs)]
        long_iso.extend(cluster_long)
        long_mem=set().union(*[c['members'] for c in cluster_long]) if cluster_long else set()
        # 3. MAP members of non-long classes against the long isoforms.
        for c in cobjs:
            if c in cluster_long: continue
            for s in c['members']:
                hits=[d for d in cluster_long if contig(tok[s],d['tk'])]
                if not hits: orphan.append((s,counts[s]))
                elif len(hits)==1: unique.append((s,counts[s]))
                else:
                    kh=[h for h in hits if h['is_known']]; pool=kh if kh else hits
                    pool.sort(key=lambda h:(-h['reads'],-h['blocks']))
                    ambiguous.append((s,counts[s],len(hits)))
    long_keys=set().union(*[m['members'] for m in long_iso]) if long_iso else set()
    res=dict(total=total,long_iso=long_iso,unique=unique,ambiguous=ambiguous,orphan=orphan,
             counts=counts,tok=tok,n_clusters=len(clusters))
    _check(res,gene); return res

def _check(r,gene):
    counts,tok,li=r['counts'],r['tok'],r['long_iso']
    lk=set().union(*[m['members'] for m in li]) if li else set()
    assert sum(m['reads'] for m in li)+sum(n for _,n in r['unique'])+sum(n for _,n,_ in r['ambiguous'])+sum(n for _,n in r['orphan'])==r['total'],f"{gene}:reads"
    part=lk|{s for s,_ in r['unique']}|{s for s,_,_ in r['ambiguous']}|{s for s,_ in r['orphan']}
    assert part==set(counts),f"{gene}:partition"
    if counts:
        top=max(counts,key=counts.get); assert top in lk,f"{gene}:top {counts[top]}rd not long"
    for s,_,_ in r['ambiguous']: assert sum(1 for m in li if contig(tok[s],m['tk']))>=2,f"{gene}:amb"

if __name__=="__main__":
    PATH="/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt"
    def load(g):
        c=collections.Counter()
        with open(PATH) as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                if len(p)>=5 and p[0]==g: c[p[2]]+=1
        return c
    genes={"ENSG00000130429":"ARF","ENSG00000169756":"LIMS1","ENSG00000133112":"TPT1",
           "ENSG00000196565":"HBG2","ENSG00000092841":"MYL6","ENSG00000166710":"B2M",
           "ENSG00000005961":"big1","ENSG00000187109":"big2","ENSG00000035403":"VCL",
           "ENSG00000012048":"BRCA1"}
    fails=0
    for g,n in genes.items():
        c=load(g)
        if not c: continue
        try:
            r=bin_gene(c,known={},gene=g)
            b=[mi['blocks'] for mi in r['long_iso']]
            ar=sum(x[1] for x in r['ambiguous']); ur=sum(x[1] for x in r['unique']); orr=sum(x[1] for x in r['orphan'])
            print(f"OK {n:6s}: {r['total']:7,}rd {len(c):4d}st {r['n_clusters']:3d}cl -> LONG {len(r['long_iso']):3d}iso(blk{min(b)}-{max(b)} mn{statistics.mean(b):.1f}) | uniq {len(r['unique'])}/{ur}rd  AMB {len(r['ambiguous'])}/{ar}rd  orph {len(r['orphan'])}/{orr}rd")
        except AssertionError as e:
            fails+=1; print(f"FAIL {n}: {e}")
    print(f"\n{'ALL INVARIANTS PASS' if fails==0 else str(fails)+' FAILED'}")
