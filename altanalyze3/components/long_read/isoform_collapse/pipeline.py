"""Production long-read isoform-collapse pipeline (Stages 1-3).

Importable, deployable entry point. See PIPELINE_README.md (in the sibling isoform-collapse/ dir,
also copied here) for the full algorithm description, measured benchmarks, and ITGA2B assessment.

Public API:
    run_pipeline(samples, outdir, nproc=8, min_total=3, write_h5ad=True) -> result dict
        samples: list of (sample_name, h5ad_path, transcript_associations_path)
    stage1_collapse(sample_ta_paths, nproc) -> gene_results
    stage2_outliers(gene_results, min_total) -> (catalog, kept_struct2exemplar)
    stage3_rekey_h5ad(sample, h5ad_path, ta_path, struct2exemplar, outdir) -> out_path

The multiprocessing worker (_collapse_chrom) is a module-level function so it pickles correctly
(this requires the package to be importable -- hence the isoform_collapse/ underscore package).
"""

from __future__ import annotations

import collections
import os
import resource
import time

from ..isoform_collapse_utils import structure_tokens_for_containment
from .scored_collapse import collapse_gene

DEFAULT_REF = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'


def _self_peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def load_gene_chrom(ref):
    gc = {}
    with open(ref) as f:
        next(f)
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 3 and t[0] not in gc:
                gc[t[0]] = t[2]
    return gc


def _run_parallel_or_serial(args, nproc):
    """Cross-platform (Unix fork / Mac+Windows spawn) parallel map of _collapse_chrom over args,
    with an automatic SERIAL fallback whenever multiprocessing is unavailable or fails.

    - nproc <= 1 or a single partition  -> serial.
    - Else try a Pool; prefer the 'fork' start method on Unix (fast, no re-import), fall back to the
      platform default ('spawn' on Mac/Windows). Any failure (PicklingError, OSError, BrokenPool,
      frozen exe, etc.) -> log and run serially. Serial always produces identical results.
    """
    if nproc <= 1 or len(args) <= 1:
        return [_collapse_chrom(a) for a in args]
    try:
        import multiprocessing as mp
        ctx = None
        try:
            if 'fork' in mp.get_all_start_methods():
                ctx = mp.get_context('fork')      # Unix: fast, avoids re-import / spawn pitfalls
        except Exception:
            ctx = None
        pool_factory = ctx.Pool if ctx is not None else mp.Pool
        with pool_factory(min(nproc, len(args))) as pool:
            return pool.map(_collapse_chrom, args)
    except Exception as e:  # any mp failure -> deterministic serial fallback
        try:
            import warnings
            warnings.warn(f"[isoform_collapse] multiprocessing failed ({type(e).__name__}: {e}); "
                          f"falling back to single-processor mode.")
        except Exception:
            pass
        return [_collapse_chrom(a) for a in args]


# ------------------------------------------------------------------ STAGE 1 --
def _load_enst_for_genes(enst_cache, gset):
    """Load {gene: {structure: ENST}} for the given gene set from the cached reference TSV."""
    ref = collections.defaultdict(dict)
    if not enst_cache or not os.path.exists(enst_cache):
        return ref
    with open(enst_cache) as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            g, struct, enst = parts
            if g in gset:
                ref[g][struct] = enst
    return ref


def _collapse_chrom(args):
    """Worker: pooled integrated collapse for all genes on one chromosome across all samples."""
    chrom, genes, sample_paths, enst_cache = args
    gset = set(genes)
    gene_struct = collections.defaultdict(collections.Counter)
    gene_mol = collections.defaultdict(dict)
    gene_mol_sample = collections.defaultdict(dict)   # gene -> {exemplar_mol: defining_sample}
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
                struct, mol, source = p[2], p[3], p[4]
                gene_struct[g][struct] += 1
                d = gene_mol[g]
                if struct not in d:
                    # Exemplar id = "molecule.source_sample" (e.g. 6020262.EP_SRSF2_Ctrl) so the
                    # source sample is encoded in the name (molecule ids can collide across samples).
                    # The final isoform id is then "gene:molecule.sample" -> still exactly one colon,
                    # so the existing annotation's feature.split(':') gives (gene, molecule.sample).
                    # source = col 5 of transcript_associations = the original gff basename.
                    d[struct] = f"{mol}.{source}"
                    gene_mol_sample[g][f"{mol}.{source}"] = source
    gene_enst = _load_enst_for_genes(enst_cache, gset)
    out = {}
    for g, sr in gene_struct.items():
        enst = gene_enst.get(g)
        all_structs = list(sr) + ([s for s in enst if s not in sr] if enst else [])
        tok = {s: tuple(structure_tokens_for_containment(s, g)) for s in all_structs}
        res = collapse_gene(sr, tok, min_reads=1, enst=enst)
        enst_ids = res.get('enst', {})
        reps = set(res['long_reps']) | set(res['other_reps'])
        final = {s: s for s in reps}
        for c, par in res['assignment'].items():
            final[c] = par
        incoming = collections.defaultdict(list)
        for c, par in res['assignment'].items():
            incoming[par].append(c)
        rep_total = {r: sr.get(r, 0) + sum(sr.get(s, 0) for s in incoming.get(r, [])) for r in reps}
        long_rep_set = set(res['long_reps'])

        # exemplar id: ENST if the representative is a reference transcript, else its molecule id.
        def exid(r):
            return enst_ids[r] if r in enst_ids else gene_mol[g][r]

        out[g] = dict(
            struct2exemplar={s: exid(final[s]) for s in final},
            exemplar_blocks={exid(r): res['blocks'][r] for r in reps},
            exemplar_total={exid(r): rep_total[r] for r in reps},
            exemplar_bin={exid(r): ('long' if r in long_rep_set else 'other') for r in reps},
            exemplar_orig={exid(r): sr.get(r, 0) for r in reps},
            exemplar_sample={exid(r): ('ENSEMBL' if r in enst_ids
                                       else gene_mol_sample[g].get(gene_mol[g].get(r))) for r in reps},
            exemplar_known={exid(r): (r in enst_ids) for r in reps},
        )
    return chrom, out, _self_peak_mb()


def stage1_collapse(sample_ta_paths, nproc=8, ref=DEFAULT_REF, enst_cache=None):
    """Chromosome-parallel integrated collapse across all samples.
    sample_ta_paths: list of (sample_name, transcript_associations_path).
    enst_cache: optional path to the cached ENST_reference_structures.tsv (gene, structure, ENST).
      When given, reference transcripts are injected and exactly-matching observed isoforms are
      named/flagged by their ENST.
    Returns (gene_results, wallclock_s, max_worker_rss_mb, n_partitions)."""
    gc = load_gene_chrom(ref)
    sgenes = set()
    for _, path in sample_ta_paths:
        with open(path) as f:
            for line in f:
                i = line.find("\t")
                if i > 0:
                    sgenes.add(line[:i])
    chrom_genes = collections.defaultdict(list)
    for g in sgenes:
        chrom_genes[gc.get(g, 'unplaced')].append(g)
    parts = sorted(chrom_genes.items(), key=lambda kv: -len(kv[1]))
    args = [(chrom, genes, sample_ta_paths, enst_cache) for chrom, genes in parts]
    t = time.time()
    results = _run_parallel_or_serial(args, nproc)
    dt = time.time() - t
    gene_results = {}
    max_rss = 0.0
    for chrom, out, rss in results:
        gene_results.update(out)
        max_rss = max(max_rss, rss)
    return gene_results, dt, max_rss, len(parts)


# ------------------------------------------------------------------ STAGE 2 --
def stage2_outliers(gene_results, min_total=3):
    """Cross-sample outlier removal: keep final isoforms with total reads >= min_total.
    Default min_total=3 removes isoforms with <=2 total reads across all samples (keeps >=3).
    Returns (catalog rows, kept_struct2exemplar)."""
    catalog = []
    kept = {}
    for g, r in gene_results.items():
        # Keep rule:
        #  - known (ENST) isoform: kept at ANY positive read count (a known transcript an observed
        #    isoform matched is real even at low depth) -- but ONLY if reads actually collapsed into
        #    it (tot > 0). A 0-read ENST clustered with reads but nothing matched/collapsed into it,
        #    so it is NOT an observed isoform and must NOT appear anywhere (catalog/h5ad/GFF).
        #  - novel isoform: kept at the cross-sample min_total threshold.
        known_map = r.get('exemplar_known', {})
        kept_ex = {ex for ex, tot in r['exemplar_total'].items()
                   if (known_map.get(ex) and tot > 0) or (not known_map.get(ex) and tot >= min_total)}
        kept[g] = {s: ex for s, ex in r['struct2exemplar'].items() if ex in kept_ex}
        for ex in kept_ex:
            catalog.append((g, ex, r['exemplar_blocks'][ex], r['exemplar_total'][ex],
                            r['exemplar_bin'][ex], known_map.get(ex, False)))
    return catalog, kept


# ---------------------------------------------------- PROTEIN PREDICTION ------
def stage_protein(gene_results, kept_struct2exemplar, sample_gff_paths, outdir,
                  genome_fasta, ref_gff=None, log=print):
    """Build the gff_process-style outputs for the FINAL collapsed isoforms and run the existing
    isoform_translation.gff_translate to produce protein_summary.txt (protein length / NMD per
    isoform) -- the file the downstream differential annotation consumes.

    Mirrors gff_process.consolidateLongReadGFFs' combined.gff: copy each kept exemplar molecule's
    genomic exon/transcript records from the sample GFF that DEFINED it, into gff-output/combined.gff.
    Also writes gff-output/transcript_associations.txt for the final isoforms (gene, strand,
    structure, exemplar_mol, source) so gff_translate's transcript->gene / intron-retention lookups
    work. sample_gff_paths: dict sample_name -> raw .gff(.gz) path.
    """
    import gzip
    gff_dir = os.path.join(outdir, 'gff-output')
    os.makedirs(gff_dir, exist_ok=True)

    # 1. exemplars to retain, grouped by their defining sample: {sample: set(exemplar_mol)}
    keep_by_sample = collections.defaultdict(set)
    exemplar_gene = {}
    for g, r in gene_results.items():
        for ex in r['exemplar_total']:
            if ex in {v for v in kept_struct2exemplar.get(g, {}).values()}:
                keep_by_sample[r['exemplar_sample'][ex]].add(ex)
                exemplar_gene[ex] = g

    # 2. combined.gff: copy genomic records for kept exemplar transcript_ids from each sample GFF.
    combined = os.path.join(gff_dir, 'combined.gff')
    written = 0
    with open(combined, 'w') as out:
        for sample, keep in keep_by_sample.items():
            path = sample_gff_paths.get(sample)
            if not path or not os.path.exists(path):
                log(f"[protein] WARN: no raw GFF for sample {sample}; its exemplars skipped")
                continue
            opener = gzip.open if path.endswith('.gz') else open
            with opener(path, 'rt') as fh:
                for line in fh:
                    if not line or line[0] == '#':
                        continue
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) != 9 or cols[2] == 'CDS' or cols[2] == 'gene':
                        continue
                    info = cols[8]
                    ti = next((i for i, x in enumerate(info.split(';')) if 'transcript_id' in x), -1)
                    if ti < 0:
                        continue
                    td = '=' if '=' in info else '"'
                    try:
                        tid = info.split(';')[ti].split(td)[1]
                    except Exception:
                        continue
                    if tid in keep:
                        out.write(line)
                        written += 1
    log(f"[protein] wrote {combined} ({written:,} genomic records for {sum(len(s) for s in keep_by_sample.values()):,} final exemplars)")

    # 3. transcript_associations for final isoforms (gene, strand, structure, exemplar, source)
    ta = os.path.join(gff_dir, 'transcript_associations.txt')
    with open(ta, 'w') as out:
        for g, r in gene_results.items():
            kept = kept_struct2exemplar.get(g, {})
            # one row per final isoform: use the exemplar's own structure
            seen = set()
            for struct, ex in kept.items():
                if ex in seen:
                    continue
                seen.add(ex)
                # exemplar's structure is the structure that maps to itself
                ex_struct = next((s for s, e in kept.items() if e == ex and r['struct2exemplar'].get(s) == ex), struct)
                out.write(f"{g}\t.\t{ex_struct}\t{ex}\t{r['exemplar_sample'][ex]}\n")
    log(f"[protein] wrote {ta}")

    # 4. run the EXISTING gff_translate -> protein_summary.txt (+ fastas) in outdir
    from .. import isoform_translation as isot
    cwd = os.getcwd()
    try:
        os.chdir(outdir)
        cds, transcripts, proteins = isot.gff_translate(combined, genome_fasta, ref_gff, ta)
        try:
            from Bio import SeqIO
            with open('protein_sequences.fasta', 'w') as f:
                SeqIO.write(proteins, f, 'fasta')
            with open('transcript_sequences.fasta', 'w') as f:
                SeqIO.write(transcripts, f, 'fasta')
            with open('orf_sequences.fasta', 'w') as f:
                SeqIO.write(cds, f, 'fasta')
        except Exception as e:
            log(f"[protein] fasta export note: {type(e).__name__}: {e}")
    finally:
        os.chdir(cwd)
    log(f"[protein] gff_translate done -> {os.path.join(outdir, 'protein_summary.txt')} "
        f"({len(proteins)} proteins)")
    return combined, ta


# ------------------------------------------------------------------ STAGE 3 --
def stage3_rekey_h5ad(sample, h5ad_path, ta_path, kept_struct2exemplar, outdir):
    """Re-key one sample's per-read h5ad to final isoform ids via sparse grouping matrix.
    Memory-optimized: builds var->final only for surviving molecules (no full molecule->struct map).
    Returns (out_path, n_final, raw_reads, final_reads, peak_rss_mb, elapsed_s)."""
    import numpy as np
    import scipy.sparse as sp
    import anndata as ad
    import pandas as pd
    t = time.time()
    # flatten kept map to (gene,structure) -> "gene:exemplar".  Final isoform ids are namespaced
    # by gene so the output is directly consumable by the existing pseudo_cluster_counts /
    # compute_differentials (which do var_name.split(":")[0] to recover the gene for ratios).
    s2e = {}
    for g, m in kept_struct2exemplar.items():
        for struct, ex in m.items():
            s2e[(g, struct)] = f"{g}:{ex}"
    var2final = {}
    with open(ta_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            ex = s2e.get((p[0], p[2]))
            if ex is not None:
                var2final[f"{p[0]}:{p[3]}"] = ex
    del s2e
    a = ad.read_h5ad(h5ad_path)
    final_ids = sorted(set(var2final.values()))
    fidx = {f: j for j, f in enumerate(final_ids)}
    rows = []
    cols = []
    for i, v in enumerate(a.var_names):
        ex = var2final.get(v)
        if ex is not None:
            rows.append(i)
            cols.append(fidx[ex])
    G = sp.coo_matrix((np.ones(len(rows), np.int64), (rows, cols)),
                      shape=(a.shape[1], len(final_ids))).tocsr()
    Xf = (a.X.astype(np.int64)) @ G
    out = ad.AnnData(X=Xf.tocsr(), obs=a.obs.copy(), var=pd.DataFrame(index=final_ids))
    base = os.path.basename(h5ad_path).split('.h5ad')[0]
    # Output name matches the existing convention (exportConsensusIsoformMatrix writes
    # "<base>-isoform.h5ad"), so this is a drop-in replacement for the isoform matrix.
    out_path = os.path.join(outdir or os.path.dirname(h5ad_path), f"{base}-isoform.h5ad")
    out.write_h5ad(out_path, compression='gzip')
    return (out_path, len(final_ids), int(a.X.sum()), int(Xf.sum()), _self_peak_mb(), time.time() - t)


# ------------------------------------------------------------------ DRIVER ---
def run_pipeline(samples, outdir, nproc=8, min_total=3, ref=DEFAULT_REF, write_h5ad=True,
                 genome_fasta=None, ref_gff=None, enst_cache=None, log=print):
    """End-to-end. samples: list of (name, h5ad_path, transcript_associations_path[, raw_gff_path]).
    If genome_fasta is given AND samples include a raw_gff_path, protein prediction is run
    (combined.gff for final isoforms -> existing gff_translate -> protein_summary.txt).
    If enst_cache (ENST_reference_structures.tsv) is given, Ensembl transcripts are injected and
    exactly-matching observed isoforms are named/flagged by their ENST (known)."""
    os.makedirs(outdir, exist_ok=True)
    sample_ta = [(s[0], s[2]) for s in samples]
    sample_gff_paths = {s[0]: s[3] for s in samples if len(s) > 3}
    gene_results, t1, max_rss, nparts = stage1_collapse(sample_ta, nproc=nproc, ref=ref,
                                                        enst_cache=enst_cache)
    n_final = sum(len(r['exemplar_total']) for r in gene_results.values())
    log(f"[stage1] genes={len(gene_results):,} final={n_final:,} {t1:.1f}s maxRSS={max_rss:.0f}MB parts={nparts}")

    catalog, kept = stage2_outliers(gene_results, min_total=min_total)
    n_known = sum(1 for row in catalog if len(row) > 5 and row[5])
    log(f"[stage2] catalog={len(catalog):,} (dropped {n_final-len(catalog):,} outliers, "
        f"min_total={min_total}); known(ENST)={n_known:,} novel={len(catalog)-n_known:,}")

    catpath = os.path.join(outdir, "FINAL_isoform_catalog.tsv")
    with open(catpath, 'w') as o:
        o.write("gene\tfinal_isoform_id\texon_blocks\ttotal_reads\tbin\tknown\n")
        for row in sorted(catalog, key=lambda r: (r[0], -r[3])):
            g, ex, blk, tot, b = row[:5]
            known = 'known' if (len(row) > 5 and row[5]) else 'novel'
            o.write(f"{g}\t{ex}\t{blk}\t{tot}\t{b}\t{known}\n")
    mappath = os.path.join(outdir, "FINAL_structure_to_exemplar.tsv")
    with open(mappath, 'w') as o:
        o.write("gene\tstructure\tfinal_isoform_id\n")
        for g, s2e in kept.items():
            for struct, ex in s2e.items():
                o.write(f"{g}\t{struct}\t{ex}\n")
    log(f"[stage2] wrote {catpath} and {mappath}")

    # PROTEIN PREDICTION: combined.gff for final isoforms -> existing gff_translate -> protein_summary.txt
    protein = None
    if genome_fasta and sample_gff_paths:
        try:
            protein = stage_protein(gene_results, kept, sample_gff_paths, outdir,
                                    genome_fasta, ref_gff=ref_gff, log=log)
        except Exception as e:
            log(f"[protein] FAILED ({type(e).__name__}: {e}); continuing without protein prediction")

    h5ad_results = []
    if write_h5ad:
        for s in samples:
            name, h5, ta = s[0], s[1], s[2]
            r = stage3_rekey_h5ad(name, h5, ta, kept, outdir=None)
            log(f"[stage3:{name}] -> {r[1]:,} final isoforms, reads {r[2]:,}->{r[3]:,} "
                f"({100*r[3]/r[2]:.0f}%), {r[5]:.1f}s peakRSS={r[4]:.0f}MB")
            h5ad_results.append((name,) + r)
    return dict(gene_results=gene_results, catalog=catalog, kept=kept, protein=protein,
                stage1=(t1, max_rss, nparts), h5ad=h5ad_results)
