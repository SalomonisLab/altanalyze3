import sys; sys.path.insert(0,'/tmp')
from cluster_binner2 import bin_gene
def S(*ex): return "|".join(ex)
DOM      = S("E2.1","E3.1","E4.1","E5.1","E6.1")               # DOMINANT, abundant (the real isoform)
RARE_EXT = S("E1.1","E2.1","E3.1","E4.1","E5.1","E6.1","E7.1") # rare long 5'/3' EXTENSION of DOM (artifact-ish, 1 read)
DOM_5t   = S("E4.1","E5.1","E6.1")                              # 3' suffix of DOM only
SKIP     = S("E2.1","E4.1","E5.1","E6.1")                       # exon-skip, distinct long iso
B        = S("E1.1","E8.1","E9.1")
counts={DOM:5000, RARE_EXT:1, DOM_5t:50, SKIP:300, B:200}
# CORRECT answer:
#  - DOM and RARE_EXT are the SAME isoform (RARE_EXT is a 5'/3' extension of DOM) -> ONE long iso,
#    representative = DOM (5000 reads), fullest extent = RARE_EXT. The dominant must be 'long'.
#  - SKIP = separate long iso. B = separate long iso (own cluster).
#  - DOM_5t (E4.1|E5.1|E6.1) is a 3' suffix of DOM and of SKIP and of RARE_EXT -> AMBIGUOUS.
r=bin_gene(counts,known={},gene="SYN")
lk=set().union(*[m['members'] for m in r['long_iso']])
reps={m['rep'] for m in r['long_iso']}
ok=True
def chk(n,c):
    global ok; print(f"  [{'PASS' if c else 'FAIL'}] {n}"); ok=ok and c
print("LONG isoforms:")
for m in r['long_iso']:
    print(f"   rep_reads={counts[m['rep']]:5d} fullest_blocks={m['blocks']} members={len(m['members'])} rep={m['rep']}")
print("AMB:",{s for s,_,_ in r['ambiguous']}, "UNIQUE:",{s for s,_ in r['unique']},"ORPH:",{s for s,_ in r['orphan']})
print("\n--- checks (models real data: rare long extension must NOT demote dominant) ---")
chk("DOM (5000rd) is a long isoform", DOM in lk)
chk("DOM and RARE_EXT in SAME long iso", any(DOM in m['members'] and RARE_EXT in m['members'] for m in r['long_iso']))
chk("DOM is the representative (not the rare ext)", any(m['rep']==DOM for m in r['long_iso'] if RARE_EXT in m['members']))
chk("SKIP separate long iso", SKIP in reps)
chk("B separate long iso", B in reps)
chk("top structure (DOM) is long", max(counts,key=counts.get) in lk)
print(f"\n{'SYNTHETIC PASSES' if ok else 'SYNTHETIC FAILED'}")
