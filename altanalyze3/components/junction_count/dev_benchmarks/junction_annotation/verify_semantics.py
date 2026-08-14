import sys
sys.path.insert(0,"/Users/saljh8/Documents/GitHub/altanalyze3")
import altanalyze3.components.long_read.gff_process as gp
EXON="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
ec,gd,sd = gp.importEnsemblGenes(EXON, include_introns=True)

# 1. chromosome of ENSG00000081913 as recorded in the reference file
import subprocess
out = subprocess.run(["grep","-m","2","-P","^ENSG00000081913\t",EXON],capture_output=True,text=True).stdout
print("ENSG00000081913 reference rows:\n", out.strip())

# 2. does find_completely_novel_annot ignore chr?  reproduce its scan for the sample case
pos, strand = 62747405, '+'
hits=[g for g in gd if gd[g][0][0] <= pos <= gd[g][-1][1] and sd[g]==strand]
print(f"\nscan for pos={pos} strand={strand}: {len(hits):,} genes satisfy the range+strand test (chr never tested)")
print("first hit returned by dict iteration order:", hits[0] if hits else None)

# 3. minus-strand reachability
minus=[g for g in gd if sd[g]=='-']
ok=[g for g in minus if gd[g][0][0] <= gd[g][-1][1]]
print(f"\nminus-strand genes: {len(minus):,}; of these, lo<=hi holds for {len(ok):,}")

# 4. per-gene exon counts -> cost of findNovelSpliceSite
import statistics
ln=[len(v) for v in gd.values()]
print(f"\nexons+introns per gene: median={statistics.median(ln):.0f} mean={statistics.mean(ln):.1f} max={max(ln):,}")

# 5. how many distinct chromosomes / genes per chromosome
from collections import Counter
chrs=Counter()
seen=set()
for (c,p,s,site) in ec:
    g=ec[(c,p,s,site)][0]
    if g not in seen: seen.add(g); chrs[c]+=1
main={f"chr{i}" for i in range(1,23)}|{"chrX","chrY"}
print(f"\ndistinct chromosomes in reference: {len(chrs):,}; genes on main chromosomes: {sum(v for k,v in chrs.items() if k in main):,} of {sum(chrs.values()):,}")
print("largest per-chromosome gene count:", chrs.most_common(3))
