"""End-to-end isoform-collapse pipeline (Stages 1-3), chromosome-parallel integrated.

STAGE 1  (integrated collapse, chromosome-parallel)
  - Pool reads across ALL samples per gene (identical structures unify automatically).
  - Per gene: scored collapse (collapse_gene) -> final isoforms (reps) + molecule->exemplar map.
  - Partition genes by chromosome; run in parallel. Bit-identical to single-process (genes are
    chromosome-local). Each worker holds one chromosome -> low memory.
  - Keeps ALL structures (no per-sample read filter).

STAGE 2  (cross-sample outlier removal + final naming)
  - Drop final isoforms whose TOTAL reads across all samples < min_total (default 2; "observed only
    once across all samples is ignored").
  - Final isoform id = exemplar molecule id (sample-tagged first molecule of the rep structure).
    (Ensembl prioritization, if desired, plugs in here; not required per the scored model.)

STAGE 3  (per-sample h5ad re-keying)
  - For each sample, map its per-read molecule ids -> final isoform id via the structure->exemplar
    map, and aggregate counts into a final-isoform x cell sparse matrix (one h5ad per sample).

Assessment gene: ITGA2B = ENSG00000005961.
"""

from __future__ import annotations

import collections
import os
import sys
import time
import resource
import importlib.util
import multiprocessing as mp

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location("scored_collapse", os.path.join(_HERE, "scored_collapse.py"))
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)
from altanalyze3.components.long_read.isoform_collapse_utils import structure_tokens_for_containment

REF = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'
ITGA2B = 'ENSG00000005961'


def _self_peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def load_gene_chrom(ref=REF):
    gc = {}
    with open(ref) as f:
        next(f)
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 3 and t[0] not in gc:
                gc[t[0]] = t[2]
    return gc


# ---------------------------------------------------------------- STAGE 1 ----
def _collapse_chrom(args):
    """Worker: pooled integrated collapse for all genes on one chromosome across all samples.
    Returns per-gene results: gene -> dict(final_exemplar, struct2exemplar, rep_reads, blocks)."""
    chrom, genes, sample_paths = args
    gset = set(genes)
    gene_struct = collections.defaultdict(collections.Counter)
    gene_mol = collections.defaultdict(dict)  # structure -> first molecule (sample-tagged)
    for sample, path in sample_paths:
        with open(path) as f:
            for line in f:
                i = line.find("\t")
                if i < 0:
                    continue
                g = line[:i]
                if g not in gset:
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                struct, mol = p[2], p[3]
                gene_struct[g][struct] += 1
                d = gene_mol[g]
                if struct not in d:
                    d[struct] = f"{sample}:{mol}"
    out = {}
    for g, sr in gene_struct.items():
        tok = {s: tuple(structure_tokens_for_containment(s, g)) for s in sr}
        res = sc.collapse_gene(sr, tok, min_reads=1)
        reps = set(res['long_reps']) | set(res['other_reps'])
        final = {s: s for s in reps}
        for c, par in res['assignment'].items():
            final[c] = par
        # structure -> final exemplar id ; rep -> total reads (incl. collapsed-in)
        incoming = collections.defaultdict(list)
        for c, par in res['assignment'].items():
            incoming[par].append(c)
        rep_total = {r: sr[r] + sum(sr[s] for s in incoming.get(r, [])) for r in reps}
        struct2exemplar = {s: gene_mol[g][final[s]] for s in final}
        out[g] = dict(
            struct2exemplar=struct2exemplar,
            exemplar_blocks={gene_mol[g][r]: res['blocks'][r] for r in reps},
            exemplar_total={gene_mol[g][r]: rep_total[r] for r in reps},
            exemplar_bin={gene_mol[g][r]: ('long' if r in set(res['long_reps']) else 'other') for r in reps},
            exemplar_orig={gene_mol[g][r]: sr[r] for r in reps},
        )
    return chrom, out, _self_peak_mb()


def stage1(sample_paths, nproc=8):
    gc = load_gene_chrom()
    sgenes = set()
    for _, path in sample_paths:
        with open(path) as f:
            for line in f:
                i = line.find("\t")
                if i > 0:
                    sgenes.add(line[:i])
    chrom_genes = collections.defaultdict(list)
    for g in sgenes:
        chrom_genes[gc.get(g, 'unplaced')].append(g)
    parts = sorted(chrom_genes.items(), key=lambda kv: -len(kv[1]))
    args = [(chrom, genes, sample_paths) for chrom, genes in parts]
    t = time.time()
    with mp.Pool(nproc) as pool:
        results = pool.map(_collapse_chrom, args)
    dt = time.time() - t
    gene_results = {}
    max_rss = 0
    for chrom, out, rss in results:
        gene_results.update(out)
        max_rss = max(max_rss, rss)
    return gene_results, dt, max_rss, len(parts)


# ---------------------------------------------------------------- STAGE 2 ----
def stage2(gene_results, min_total=2):
    """Cross-sample outlier removal: keep final isoforms with total reads >= min_total.
    Returns final_catalog: list of (gene, exemplar_id, blocks, total_reads, bin) and the set of
    kept exemplar ids per gene, plus the surviving structure->exemplar map."""
    catalog = []
    kept_struct2exemplar = {}  # gene -> {structure: exemplar}  (only surviving exemplars)
    for g, r in gene_results.items():
        kept_ex = {ex for ex, tot in r['exemplar_total'].items() if tot >= min_total}
        s2e = {s: ex for s, ex in r['struct2exemplar'].items() if ex in kept_ex}
        kept_struct2exemplar[g] = s2e
        for ex in kept_ex:
            catalog.append((g, ex, r['exemplar_blocks'][ex], r['exemplar_total'][ex], r['exemplar_bin'][ex]))
    return catalog, kept_struct2exemplar


# ---------------------------------------------------------------- driver -----
if __name__ == '__main__':
    SAMPLES = [
        ('Ctrl', '/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output/transcript_associations.txt'),
        ('Aza', '/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Aza/gff-output/transcript_associations.txt'),
    ]
    OUTDIR = '/Users/saljh8/Dropbox/Revio/BAMs/iPSC/Ctrl/gff-output'
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    min_total = int(sys.argv[2]) if len(sys.argv) > 2 else 2

    t0 = time.time()
    gene_results, t1, max_rss, nparts = stage1(SAMPLES, nproc=nproc)
    n_final = sum(len(r['exemplar_total']) for r in gene_results.values())
    print("===== STAGE 1: integrated collapse (chromosome-parallel) =====")
    print(f"  genes={len(gene_results):,}  final isoforms(all)={n_final:,}")
    print(f"  wall-clock={t1:.1f}s  max per-worker RSS={max_rss:.0f}MB  partitions={nparts}")

    catalog, kept = stage2(gene_results, min_total=min_total)
    print(f"\n===== STAGE 2: cross-sample outlier removal (total reads >= {min_total}) =====")
    print(f"  final catalog isoforms={len(catalog):,}  (dropped {n_final-len(catalog):,} outliers)")

    # write final catalog
    catpath = f"{OUTDIR}/FINAL_isoform_catalog.tsv"
    with open(catpath, 'w') as o:
        o.write("gene\tfinal_isoform_id\texon_blocks\ttotal_reads\tbin\n")
        for g, ex, blk, tot, b in sorted(catalog, key=lambda r: (r[0], -r[3])):
            o.write(f"{g}\t{ex}\t{blk}\t{tot}\t{b}\n")
    print(f"  wrote {catpath}")

    # write surviving structure -> final exemplar map (for Stage 3 h5ad re-keying)
    mappath = f"{OUTDIR}/FINAL_structure_to_exemplar.tsv"
    with open(mappath, 'w') as o:
        o.write("gene\tstructure\tfinal_isoform_id\n")
        for g, s2e in kept.items():
            for struct, ex in s2e.items():
                o.write(f"{g}\t{struct}\t{ex}\n")
    print(f"  wrote {mappath}")

    # -------- ITGA2B assessment --------
    print("\n===== ASSESSMENT: ITGA2B (ENSG00000005961) =====")
    r = gene_results.get(ITGA2B)
    if r:
        all_ex = r['exemplar_total']
        kept_ex = {ex: tot for ex, tot in all_ex.items() if tot >= min_total}
        print(f"  final isoforms (pre-filter): {len(all_ex)}   after outlier removal: {len(kept_ex)}")
        print(f"  top final isoforms (exemplar | blocks | bin | total reads):")
        for ex in sorted(kept_ex, key=lambda e: -kept_ex[e])[:10]:
            print(f"    {ex:>18} | {r['exemplar_blocks'][ex]:2d}blk | {r['exemplar_bin'][ex]:5s} | {kept_ex[ex]:6d}")
    print(f"\nTOTAL pipeline wall-clock (S1+S2)={time.time()-t0:.1f}s")
