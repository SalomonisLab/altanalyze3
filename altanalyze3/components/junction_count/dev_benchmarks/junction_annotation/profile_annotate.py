"""Profile altanalyze3 junction annotation. Read-only. Writes nothing except stdout.

Entry point measured: altanalyze3.components.aggregate.annotate.annotate_junctions
Reference file:       /Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt

The harness builds a synthetic var table because no real juncounts BED exists on this Mac.
Junctions are drawn from the real reference coordinates, so the dictionary lookups and the
gene scans run against real data volumes.
"""
import sys, time, random, cProfile, pstats, io
import pandas as pd

REPO = "/Users/saljh8/Documents/GitHub/altanalyze3"
EXON = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt"
sys.path.insert(0, REPO)

import altanalyze3.components.long_read.gff_process as gff_process
from altanalyze3.components.aggregate import annotate as ann


class FakeAdata:
    def __init__(self, var):
        self.var = var


def load_ref():
    t0 = time.time()
    exonCoordinates, geneData, strandData = gff_process.importEnsemblGenes(EXON, include_introns=True)
    t1 = time.time()
    print(f"importEnsemblGenes: {t1-t0:.2f} s   "
          f"exonCoordinates={len(exonCoordinates):,}  geneData={len(geneData):,} genes")
    return exonCoordinates, geneData, strandData


def build_var(exonCoordinates, geneData, strandData, n_known, n_halfnovel, n_fullnovel, seed=0):
    """Build a var table with three junction classes, in the columns annotate_junctions reads."""
    rnd = random.Random(seed)
    # index donor (site 2) and acceptor (site 1) coordinates per (chr, strand)
    donors, acceptors = {}, {}
    for (chrom, pos, strand, site), v in exonCoordinates.items():
        (donors if site == 2 else acceptors).setdefault((chrom, strand), []).append(pos)
    keys = [k for k in donors if k in acceptors and len(donors[k]) > 100 and len(acceptors[k]) > 100]

    rows = []
    for _ in range(n_known):
        chrom, strand = rnd.choice(keys)
        rows.append((chrom, rnd.choice(donors[(chrom, strand)]), rnd.choice(acceptors[(chrom, strand)]), strand))
    for _ in range(n_halfnovel):
        chrom, strand = rnd.choice(keys)
        rows.append((chrom, rnd.choice(donors[(chrom, strand)]),
                     rnd.choice(acceptors[(chrom, strand)]) + rnd.randint(1, 40), strand))
    for _ in range(n_fullnovel):
        chrom, strand = rnd.choice(keys)
        rows.append((chrom, rnd.choice(donors[(chrom, strand)]) + rnd.randint(1, 40),
                     rnd.choice(acceptors[(chrom, strand)]) + rnd.randint(1, 40), strand))
    var = pd.DataFrame(rows, columns=["chr", "start", "end", "strand"])
    var.index = [f"{c}:{s}-{e}" for c, s, e, _ in rows]
    return var


def timed(label, var):
    a = FakeAdata(var.copy())
    t0 = time.time()
    ann.annotate_junctions(a, EXON)
    dt = time.time() - t0
    n = len(var)
    print(f"{label:<28} n={n:>7,}  total={dt:8.2f} s   per-junction={1e3*dt/n:8.3f} ms")
    return a.var["annotation"].tolist(), dt


if __name__ == "__main__":
    exonCoordinates, geneData, strandData = load_ref()

    # How often can find_completely_novel_annot succeed at all?
    plus = sum(1 for g in geneData if strandData[g] == '+')
    minus = sum(1 for g in geneData if strandData[g] == '-')
    satisfiable = sum(1 for g in geneData if geneData[g][0][0] <= geneData[g][-1][1])
    print(f"genes: +strand={plus:,}  -strand={minus:,}  "
          f"range-check lo<=hi satisfiable for {satisfiable:,} of {len(geneData):,}")

    print()
    # class-by-class cost; the reference load is repeated inside annotate_junctions each call
    v_known = build_var(exonCoordinates, geneData, strandData, 2000, 0, 0)
    v_half = build_var(exonCoordinates, geneData, strandData, 0, 2000, 0, seed=1)
    v_full = build_var(exonCoordinates, geneData, strandData, 0, 0, 200, seed=2)
    ann_known, _ = timed("known donor+acceptor", v_known)
    ann_half, _ = timed("one novel splice site", v_half)
    ann_full, _ = timed("both sites novel", v_full)

    print()
    print("sample annotations:")
    for lbl, arr in (("known", ann_known), ("half-novel", ann_half), ("full-novel", ann_full)):
        print(f"  {lbl:<12} {arr[:3]}")

    print()
    print("cProfile on 2,000 known + 200 full-novel:")
    v_mix = pd.concat([v_known, v_full])
    a = FakeAdata(v_mix.copy())
    pr = cProfile.Profile()
    pr.enable()
    ann.annotate_junctions(a, EXON)
    pr.disable()
    s = io.StringIO()
    pstats.Stats(pr, stream=s).sort_stats("cumulative").print_stats(22)
    print(s.getvalue())
