import sys
sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
import altanalyze3.components.long_read.gff_process as gp
EXON="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
ec,gd,sd = gp.importEnsemblGenes(EXON, include_introns=True)

# claim A: some genes have no surviving key in exonCoordinates (key collisions overwrite them)
in_ec=set(v[0] for v in ec.values())
missing=[g for g in gd if g not in in_ec]
print(f"A) genes in geneData = {len(gd):,}; genes reachable from exonCoordinates = {len(in_ec):,}; "
      f"unreachable = {len(missing)}")
print("   examples:", missing[:5])

# claim B: for + strand genes, geneData[g][0][0] == min coord and geneData[g][-1][1] == max coord
badlo=badhi=0
for g,ex in gd.items():
    if sd[g]!='+': continue
    lo=min(min(a,b) for a,b,_ in ex); hi=max(max(a,b) for a,b,_ in ex)
    if ex[0][0]!=lo: badlo+=1
    if ex[-1][1]!=hi: badhi+=1
npl=sum(1 for g in gd if sd[g]=='+')
print(f"B) + strand genes = {npl:,}; where [0][0] != min coord: {badlo}; where [-1][1] != max coord: {badhi}")

# claim C: minus-strand true span vs the legacy expression
bad=0
for g,ex in gd.items():
    if sd[g]!='-': continue
    if ex[0][0] <= ex[-1][1]: bad+=1
print(f"C) - strand genes = {sum(1 for g in gd if sd[g]=='-'):,}; where legacy lo<=hi holds: {bad}")

# claim D: do genes span more than one chromosome in the reference?
from collections import defaultdict
gc=defaultdict(set)
with open(EXON) as f:
    next(f)
    for line in f:
        t=line.rstrip('\n').split('\t')
        gc[t[0]].add(t[2])
multi=[g for g,v in gc.items() if len(v)>1]
print(f"D) genes with >1 chromosome in the reference file: {len(multi)}; total genes in file: {len(gc):,}")
