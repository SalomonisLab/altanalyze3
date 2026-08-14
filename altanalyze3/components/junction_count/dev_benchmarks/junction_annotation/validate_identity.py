"""Prove the new annotate_junctions reproduces the frozen original, branch by branch.

Compares, elementwise, the annotation string from:
  OLD  dev_benchmarks/junction_annotation/annotate_reference_frozen.py  (copy taken before editing)
  NEW  altanalyze3/components/aggregate/annotate.py  with novel_gene_mode='legacy'

Then reports how novel_gene_mode='corrected' differs from 'legacy' on the same junctions.
Writes a per-branch table and the full mismatch list. Reads the reference file only.
"""
import sys, time, random, importlib.util
from collections import Counter, defaultdict
import pandas as pd

REPO = "/Users/saljh8/Documents/GitHub/altanalyze3"
BENCH = f"{REPO}/dev_benchmarks/junction_annotation"
EXON = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
sys.path.insert(0, REPO)

import altanalyze3.components.long_read.gff_process as gff_process
from altanalyze3.components.aggregate.annotate import annotate_junctions as new_annotate

spec = importlib.util.spec_from_file_location("frozen_ref", f"{BENCH}/annotate_reference_frozen.py")
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
old_annotate = frozen.annotate_junctions


class FakeAdata:
    def __init__(self, var):
        self.var = var


# ---------------------------------------------------------------- build a branch-covering set
ec, gd, sd = gff_process.importEnsemblGenes(EXON, include_introns=True)
gene_chr = gff_process.gene_chr

donors, acceptors = defaultdict(list), defaultdict(list)
for (chrom, pos, strand, site), v in ec.items():
    (donors if site == 2 else acceptors)[(chrom, strand)].append(pos)

# genes with an intron long enough to place a mid-intron and two near-edge probes
introns = []
for g, exons in gd.items():
    for lo, hi, name in exons:
        a, b = min(lo, hi), max(lo, hi)
        if name.startswith("I") and name.endswith(".1") and b - a > 400:
            introns.append((g, sd[g], gene_chr[g], a, b))

rnd = random.Random(7)
rows, branch = [], []


def add(chrom, s, e, strand, label):
    rows.append((chrom, s, e, strand))
    branch.append(label)


N = 60  # >= 50 rows per branch, per CLAUDE.md RULE 5

keys_pos = [k for k in donors if k in acceptors and k[1] == '+' and len(donors[k]) > 200]
keys_neg = [k for k in donors if k in acceptors and k[1] == '-' and len(donors[k]) > 200]

for keys, tag in ((keys_pos, '+'), (keys_neg, '-')):
    for _ in range(N):
        chrom, strand = rnd.choice(keys)
        d, a = rnd.choice(donors[(chrom, strand)]), rnd.choice(acceptors[(chrom, strand)])
        add(chrom, d, a, strand, f"both sites known ({tag})")
    for _ in range(N):
        chrom, strand = rnd.choice(keys)
        d, a = rnd.choice(donors[(chrom, strand)]), rnd.choice(acceptors[(chrom, strand)])
        add(chrom, d, a + rnd.randint(1, 40), strand, f"acceptor novel ({tag})")
    for _ in range(N):
        chrom, strand = rnd.choice(keys)
        d, a = rnd.choice(donors[(chrom, strand)]), rnd.choice(acceptors[(chrom, strand)])
        add(chrom, d + rnd.randint(1, 40), a, strand, f"donor novel ({tag})")
    for _ in range(N):
        chrom, strand = rnd.choice(keys)
        d, a = rnd.choice(donors[(chrom, strand)]), rnd.choice(acceptors[(chrom, strand)])
        add(chrom, d + rnd.randint(1, 40), a + rnd.randint(1, 40), strand, f"both sites novel ({tag})")

# positions placed deliberately inside introns, near each edge and in the middle
for tag, want in (('+', '+'), ('-', '-')):
    pool = [x for x in introns if x[1] == want]
    for offset_label, mk in (
        ("intron start <50nt", lambda a, b: a + rnd.randint(1, 45)),
        ("intron end <50nt", lambda a, b: b - rnd.randint(1, 45)),
        ("intron middle", lambda a, b: (a + b) // 2),
    ):
        for _ in range(N):
            g, strand, chrom, a, b = rnd.choice(pool)
            p = mk(a, b)
            add(chrom, p, p + rnd.randint(50, 5000), strand, f"{offset_label} ({tag})")

# before the first exon and after the last exon of a gene
for tag in ('+', '-'):
    pool = [g for g in gd if sd[g] == tag]
    for _ in range(N):
        g = rnd.choice(pool)
        exons = gd[g]
        lo = min(min(x, y) for x, y, _ in exons)
        add(gene_chr[g], lo - rnd.randint(1, 5000), lo - rnd.randint(1, 500), tag, f"upstream of gene ({tag})")
    for _ in range(N):
        g = rnd.choice(pool)
        exons = gd[g]
        hi = max(max(x, y) for x, y, _ in exons)
        add(gene_chr[g], hi + rnd.randint(1, 500), hi + rnd.randint(501, 5000), tag, f"downstream of gene ({tag})")

# scaffold contigs
scaff = [k for k in donors if k[0].startswith("chrCHR") and k in acceptors and len(donors[k]) > 20]
for _ in range(N):
    chrom, strand = rnd.choice(scaff)
    add(chrom, rnd.choice(donors[(chrom, strand)]) + rnd.randint(0, 30),
        rnd.choice(acceptors[(chrom, strand)]) + rnd.randint(0, 30), strand, "scaffold contig")

var = pd.DataFrame(rows, columns=["chr", "start", "end", "strand"])
var.index = [f"{c}:{s}-{e}:{st}#{i}" for i, (c, s, e, st) in enumerate(rows)]
print(f"junction set: {len(var):,} rows across {len(set(branch))} branches\n")

# ---------------------------------------------------------------- run all three
def run(fn, **kw):
    a = FakeAdata(var.copy())
    t0 = time.time()
    fn(a, EXON, **kw)
    return a.var["annotation"].tolist(), time.time() - t0

old, t_old = run(old_annotate)
new_legacy, t_new = run(new_annotate, novel_gene_mode="legacy")
new_corrected, t_corr = run(new_annotate, novel_gene_mode="corrected")

print(f"OLD (frozen original)          {t_old:8.2f} s")
print(f"NEW legacy                     {t_new:8.2f} s   speedup {t_old/t_new:6.1f}x")
print(f"NEW corrected                  {t_corr:8.2f} s\n")

# ---------------------------------------------------------------- identity, per branch
print(f"{'branch':<28} {'n':>5} {'identical':>10} {'corrected differs':>18}")
print("-" * 64)
total_mismatch = 0
per_branch = defaultdict(lambda: [0, 0, 0])
for b, o, nl, nc in zip(branch, old, new_legacy, new_corrected):
    per_branch[b][0] += 1
    if o == nl:
        per_branch[b][1] += 1
    else:
        total_mismatch += 1
    if nl != nc:
        per_branch[b][2] += 1
for b in sorted(per_branch):
    n, ident, diff = per_branch[b]
    flag = "" if ident == n else "   <-- MISMATCH"
    print(f"{b:<28} {n:>5} {ident:>10} {diff:>18}{flag}")
print("-" * 64)
print(f"{'TOTAL':<28} {len(old):>5} {len(old)-total_mismatch:>10} "
      f"{sum(v[2] for v in per_branch.values()):>18}")

print(f"\nOLD vs NEW-legacy mismatches: {total_mismatch} of {len(old):,}")
if total_mismatch:
    print("first 10 mismatches (branch | old | new):")
    shown = 0
    for b, o, nl in zip(branch, old, new_legacy):
        if o != nl and shown < 10:
            print(f"  {b} | {o} | {nl}")
            shown += 1

# ---------------------------------------------------------------- what corrected changes
changed = [(b, o, c) for b, o, c in zip(branch, new_legacy, new_corrected) if o != c]
print(f"\ncorrected vs legacy: {len(changed)} of {len(old):,} junctions change")
kinds = Counter()
for b, l, c in changed:
    lg = l.split(":")[0]
    cg = c.split(":")[0]
    if lg == "None" and cg != "None":
        kinds["legacy None -> corrected assigns a gene"] += 1
    elif lg != "None" and cg == "None":
        kinds["legacy assigned a gene -> corrected None"] += 1
    else:
        kinds["both assigned, different gene or exon"] += 1
for k, v in kinds.most_common():
    print(f"  {v:>5}  {k}")
print("\nexamples (legacy -> corrected):")
for b, l, c in changed[:6]:
    print(f"  [{b}]\n     legacy    {l}\n     corrected {c}")

sys.exit(1 if total_mismatch else 0)
