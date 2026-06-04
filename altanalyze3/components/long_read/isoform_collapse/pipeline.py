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
from .scored_collapse import collapse_gene, collapse_gene_em

# Default gene->chrom reference = the gzipped exon annotation BUNDLED in the program directory
# (components/long_read/resources/), resolved relative to this package so it works on any machine
# without an external --exon_annot. load_gene_chrom() opens it gz-aware. Human default; pass
# --exon_annot for mouse/custom. (Was a hardcoded developer macOS path that existed nowhere else.)
DEFAULT_REF = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),  # components/long_read/
    "resources", "Hs_Ensembl_exon.txt.gz",
)


def _self_peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def load_gene_chrom(ref):
    # ``ref`` is the Ensembl exon annotation (gene->chrom). Callers should pass the resolved
    # --exon_annot path; only fall back to DEFAULT_REF if nothing was supplied. Fail loudly with a
    # clear, actionable message rather than a bare FileNotFoundError deep in the collapse.
    ref = str(ref) if ref else DEFAULT_REF
    if not os.path.exists(ref):
        raise FileNotFoundError(
            f"Exon annotation (gene->chrom reference) not found: {ref}. Pass a valid --exon_annot "
            f"(e.g. Hs_Ensembl_exon.txt), or ensure the bundled "
            f"components/long_read/resources/Hs_Ensembl_exon.txt.gz ships with the package."
        )
    import gzip
    opener = (lambda p: gzip.open(p, 'rt')) if ref.endswith('.gz') else open
    gc = {}
    with opener(ref) as f:
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
    """Tier 2 worker: merge the per-sample Tier-1 exemplar tables for all genes on one chromosome,
    then run the validated containment collapse on the cross-sample SUMMED structure counts.

    Input ``sample_paths`` are the Tier-1 ``<sample>.exemplars.tsv`` files
    (gene, strand, structure, count, rep_molecule_id). Summing the Tier-1 ``count`` per structure
    across samples reproduces exactly the per-read structure counts the previously-validated pooled
    collapse fed to ``collapse_gene`` -- so the result is identical, but memory is bounded by distinct
    structures (Tier 1 already reduced reads->structures per sample)."""
    # collapse_method: 'wta' (default winner-takes-all) or 'em' (soft EM read allocation). The 5-tuple
    # form carries the method; the 4-tuple form (older callers) defaults to 'wta'.
    if len(args) == 5:
        chrom, genes, sample_paths, enst_cache, collapse_method = args
    else:
        chrom, genes, sample_paths, enst_cache = args
        collapse_method = 'wta'
    gset = set(genes)
    gene_struct = collections.defaultdict(collections.Counter)
    gene_mol = collections.defaultdict(dict)
    gene_mol_sample = collections.defaultdict(dict)   # gene -> {exemplar_mol: defining_sample}
    for sample, path in sample_paths:
        with open(path) as f:
            next(f, None)  # header
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
                # Tier-1 exemplars: gene, strand, structure, count, rep_molecule_id
                struct, count, mol = p[2], p[3], p[4]
                gene_struct[g][struct] += int(count)
                d = gene_mol[g]
                if struct not in d:
                    # Exemplar id = "rep_molecule.sample" so the source sample is encoded in the name
                    # (molecule ids can collide across samples). Final isoform id = "gene:molecule.sample"
                    # -> one colon, so the existing annotation's feature.split(':') -> (gene, molecule.sample).
                    d[struct] = f"{mol}.{sample}"
                    gene_mol_sample[g][f"{mol}.{sample}"] = sample
    gene_enst = _load_enst_for_genes(enst_cache, gset)
    out = {}
    for g, sr in gene_struct.items():
        enst = gene_enst.get(g)
        all_structs = list(sr) + ([s for s in enst if s not in sr] if enst else [])
        tok = {s: tuple(structure_tokens_for_containment(s, g)) for s in all_structs}
        if collapse_method == 'em':
            res = collapse_gene_em(sr, tok, min_reads=1, enst=enst)
        else:
            res = collapse_gene(sr, tok, min_reads=1, enst=enst)
        enst_ids = res.get('enst', {})
        reps = set(res['long_reps']) | set(res['other_reps'])
        final = {s: s for s in reps}
        for c, par in res['assignment'].items():
            final[c] = par
        incoming = collections.defaultdict(list)
        for c, par in res['assignment'].items():
            incoming[par].append(c)
        # HARD total = own reads + reads of children hard-assigned to it (the winner-takes-all count).
        # This is the SHARED catalog-membership criterion: the kept isoform SET is decided from the
        # hard totals identically for WTA and EM, so both methods yield the SAME isoform set; EM only
        # redistributes reads WITHIN that fixed set (it must not add/remove isoforms). EM additionally
        # reports its soft abundance (em_reads) for the catalog's reported read count.
        hard_total = {r: sr.get(r, 0) + sum(sr.get(s, 0) for s in incoming.get(r, [])) for r in reps}
        if collapse_method == 'em':
            em_reads = res.get('em_reads', {})
            rep_total = {r: em_reads.get(r, sr.get(r, 0)) for r in reps}
        else:
            rep_total = hard_total
        long_rep_set = set(res['long_reps'])

        # exemplar id: ENST if the representative is a reference transcript, else its molecule id.
        def exid(r):
            return enst_ids[r] if r in enst_ids else gene_mol[g][r]

        entry = dict(
            struct2exemplar={s: exid(final[s]) for s in final},
            exemplar_blocks={exid(r): res['blocks'][r] for r in reps},
            exemplar_total={exid(r): rep_total[r] for r in reps},
            exemplar_hardtotal={exid(r): hard_total[r] for r in reps},
            exemplar_bin={exid(r): ('long' if r in long_rep_set else 'other') for r in reps},
            exemplar_orig={exid(r): sr.get(r, 0) for r in reps},
            exemplar_sample={exid(r): ('ENSEMBL' if r in enst_ids
                                       else gene_mol_sample[g].get(gene_mol[g].get(r))) for r in reps},
            exemplar_known={exid(r): (r in enst_ids) for r in reps},
        )
        # EM: carry the STRUCTURE-keyed soft map forward unchanged (child_structure -> {parent_structure:
        # weight}). It is NOT converted to ids here -- structure is the canonical, cross-sample-consistent
        # key (the same key stage3 uses), so the molecule->final re-link stays correct across all samples.
        # ids are resolved only at the final re-key, exactly like the hard struct2exemplar path.
        if collapse_method == 'em':
            entry['struct2struct_soft'] = res.get('em_weights', {})
        out[g] = entry
    return chrom, out, _self_peak_mb()


def stage1_collapse(sample_ta_paths, nproc=8, ref=DEFAULT_REF, enst_cache=None, tier1_dir=None,
                    collapse_method='wta', log=print):
    """Two-tier cross-sample collapse.

    collapse_method: 'wta' (default winner-takes-all -- each ambiguous substring's reads go to its
      single highest-score long-bin parent) or 'em' (soft EM -- each substring's reads split across
      compatible parents in proportion to their abundance). Same representative set either way; only
      the per-isoform read totals differ.

    TIER 1 (per sample, in series): reduce each sample's molecule-level transcript_associations.txt
      to its distinct structures + read counts + representative molecule id (collapse_sample). Bounded
      by a sample's distinct structures, never its read count.
    TIER 2 (chromosome-parallel): merge the per-sample exemplar tables (sum counts across samples per
      structure) and run the validated containment collapse (_collapse_chrom). Summing Tier-1 counts
      reproduces the pooled per-read counts exactly -> identical result to the single pooled collapse.

    sample_ta_paths: list of (sample_name, transcript_associations_path).
    Returns (gene_results, wallclock_s, max_worker_rss_mb, n_partitions)."""
    from .collapse_sample import collapse_sample
    t0 = time.time()

    # --- TIER 1: per-sample reduction (series; each sample is bounded + independent) -------------
    if tier1_dir is None:
        tier1_dir = os.path.join(os.path.dirname(sample_ta_paths[0][1]) or '.', 'tier1')
    os.makedirs(tier1_dir, exist_ok=True)
    sample_exemplars = []   # (sample_name, exemplars_tsv_path)
    for sample, ta in sample_ta_paths:
        ex_path, _mol = collapse_sample(sample, ta, out_dir=tier1_dir, log=log)
        sample_exemplars.append((sample, ex_path))

    # --- TIER 2: merge per-sample exemplars + collapse, partitioned by chromosome ----------------
    gc = load_gene_chrom(ref)
    sgenes = set()
    for _, path in sample_exemplars:
        with open(path) as f:
            next(f, None)  # header
            for line in f:
                i = line.find("\t")
                if i > 0:
                    sgenes.add(line[:i])
    chrom_genes = collections.defaultdict(list)
    for g in sgenes:
        chrom_genes[gc.get(g, 'unplaced')].append(g)
    parts = sorted(chrom_genes.items(), key=lambda kv: -len(kv[1]))
    args = [(chrom, genes, sample_exemplars, enst_cache, collapse_method) for chrom, genes in parts]
    t = time.time()
    results = _run_parallel_or_serial(args, nproc)
    dt = time.time() - t
    gene_results = {}
    max_rss = 0.0
    for chrom, out, rss in results:
        gene_results.update(out)
        max_rss = max(max_rss, rss)
    log(f"[stage1] tier1+tier2 total {time.time()-t0:.1f}s (tier2 {dt:.1f}s)")
    return gene_results, dt, max_rss, len(parts)


# ------------------------------------------------------------------ STAGE 2 --
def stage2_outliers(gene_results, min_total=3, return_soft=False):
    """Cross-sample outlier removal: keep final isoforms with total reads >= min_total.
    Default min_total=3 removes isoforms with <=2 total reads across all samples (keeps >=3).
    Returns (catalog rows, kept_struct2exemplar) -- or, with return_soft=True, additionally a
    structure-keyed EM soft map kept_soft (gene -> {child_structure: {parent_structure: weight}})
    restricted to kept parents and renormalized (only populated for EM gene_results)."""
    catalog = []
    kept = {}
    kept_soft = {}
    for g, r in gene_results.items():
        # Keep rule:
        #  - known (ENST) isoform: kept at ANY positive read count (a known transcript an observed
        #    isoform matched is real even at low depth) -- but ONLY if reads actually collapsed into
        #    it (tot > 0). A 0-read ENST clustered with reads but nothing matched/collapsed into it,
        #    so it is NOT an observed isoform and must NOT appear anywhere (catalog/h5ad/GFF).
        #  - novel isoform: kept at the cross-sample min_total threshold.
        known_map = r.get('exemplar_known', {})
        # Catalog membership uses the HARD (winner-takes-all) totals so WTA and EM keep the SAME
        # isoform set. EM only redistributes reads within that set; it must not change which isoforms
        # exist. (Falls back to exemplar_total when hardtotal absent, e.g. older results.)
        decide = r.get('exemplar_hardtotal', r['exemplar_total'])
        kept_ex = {ex for ex, tot in decide.items()
                   if (known_map.get(ex) and tot > 0) or (not known_map.get(ex) and tot >= min_total)}
        s2e = r['struct2exemplar']
        kept[g] = {s: ex for s, ex in s2e.items() if ex in kept_ex}
        for ex in kept_ex:
            catalog.append((g, ex, r['exemplar_blocks'][ex], r['exemplar_total'][ex],
                            r['exemplar_bin'][ex], known_map.get(ex, False)))
        # EM structure-keyed soft map, filtered to kept structures & kept parents, renormalized.
        soft = r.get('struct2struct_soft')
        if soft:
            kept_structs = set(kept[g])
            gs = {}
            for child_struct, pw in soft.items():
                if child_struct not in kept_structs:
                    continue
                fw = {ps: w for ps, w in pw.items() if ps in kept_structs}
                tot = sum(fw.values())
                if tot > 0:
                    gs[child_struct] = {ps: w / tot for ps, w in fw.items()}
            if gs:
                kept_soft[g] = gs
    if return_soft:
        return catalog, kept, kept_soft
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
    # ``outdir`` is already the gff-output directory; write directly into it (do not nest another
    # gff-output/ under it).
    gff_dir = outdir
    os.makedirs(gff_dir, exist_ok=True)

    # 1. exemplars to retain. A NOVEL exemplar id is "<molecule>.<sample>" (genomic records live in
    # that sample's read GFF, transcript_id = bare molecule). A KNOWN exemplar is an ENST (version
    # stripped); its genomic records live in the reference GFF (transcript_id = ENST.version). We copy
    # both into combined.gff so gff_translate produces proteins for known AND novel isoforms.
    # For each kept NOVEL exemplar we record (sample, bare_molecule_id, (gene, structure), exemplar_id).
    # The genomic coordinates are retrieved by an indexed point lookup on each sample's
    # structure_coords.sqlite (written by gff_process at the structure-formation stage) keyed by
    # (gene, structure) -- avoiding a full re-scan of the source GFFs. If a sample's store is absent
    # (older run), we fall back to scanning that sample's read GFF by transcript_id.
    keep_mol_to_ex = collections.defaultdict(dict)        # sample -> {bare_molecule_id: exemplar_id}
    keep_struct_by_sample = collections.defaultdict(list)  # sample -> [(gene, structure, exemplar_id)]
    keep_enst = set()                                      # version-stripped ENST exemplar ids (known)
    for g, r in gene_results.items():
        kept_ex = set(kept_struct2exemplar.get(g, {}).values())
        known_map = r.get('exemplar_known', {})
        s2e = r.get('struct2exemplar', {})
        ex2struct = {}
        for s, e in s2e.items():
            ex2struct.setdefault(e, s)                     # exemplar -> its own (representative) structure
        for ex in r['exemplar_total']:
            if ex not in kept_ex:
                continue
            if known_map.get(ex):
                keep_enst.add(ex)                          # ENST id (already version stripped)
                continue
            sample = r['exemplar_sample'].get(ex)
            if not sample or sample == 'ENSEMBL':
                continue
            mol = ex[:-(len(sample) + 1)] if ex.endswith('.' + sample) else ex
            keep_mol_to_ex[sample][mol] = ex
            struct = ex2struct.get(ex)
            if struct is not None:
                # stamp the BARE molecule id as transcript_id -- byte-identical to the source GFF /
                # golden combined.gff (which carry transcript_id "<molecule>", not "<molecule>.<sample>").
                keep_struct_by_sample[sample].append((g, struct, mol))

    combined_gz = os.path.join(gff_dir, 'combined.gff.gz')
    written = 0
    n_novel = sum(len(m) for m in keep_mol_to_ex.values())
    n_from_store = 0
    n_from_scan = 0
    try:
        from ...bam import coord_store
    except Exception:
        try:
            from bam import coord_store
        except Exception:
            coord_store = None

    with gzip.open(combined_gz, 'wt') as out:
        # 2a. novel exemplars: emit genomic records grouped by sample (each sample's structures are
        # already chromosome-local in its store; we group output per chromosome within a sample).
        for sample, keep in keep_mol_to_ex.items():
            path = sample_gff_paths.get(sample)
            store_path = None
            if path:
                # structure_coords.sqlite lives next to the sample's gff-output/transcript_associations.txt
                cand = os.path.join(os.path.dirname(path), 'gff-output', 'structure_coords.sqlite')
                if os.path.exists(cand):
                    store_path = cand
            conn = coord_store.open_reader(store_path) if (coord_store and store_path) else None

            if conn is not None:
                # FAST PATH: indexed (gene, structure) lookups; stamp exemplar id as transcript_id.
                pairs = [(g, s) for (g, s, _mol) in keep_struct_by_sample.get(sample, [])]
                mol_by_pair = {(g, s): mol for (g, s, mol) in keep_struct_by_sample.get(sample, [])}
                # collect blocks then emit grouped by chromosome (genes group naturally per chrom)
                by_chrom = collections.defaultdict(list)
                for g, s, chrom, block in coord_store.fetch_structures(
                        conn, pairs, transcript_id_for=lambda gg, ss: mol_by_pair[(gg, ss)]):
                    by_chrom[chrom].append(block)
                conn.close()
                for chrom in sorted(by_chrom):
                    for block in by_chrom[chrom]:
                        out.write(block)
                        written += block.count('\n')
                        n_from_store += 1
                continue

            # FALLBACK: scan the sample's read GFF by transcript_id (no store available).
            if not path or not os.path.exists(path):
                log(f"[protein] WARN: no raw GFF or coord store for sample {sample}; its exemplars skipped")
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
                        n_from_scan += 1

        # 2b. known (ENST) exemplars: copy genomic records from the reference GFF, matching on the
        # version-stripped transcript_id (reference is "ENST.version", exemplar is "ENST").
        n_enst_records = 0
        if keep_enst and ref_gff and os.path.exists(ref_gff):
            opener = gzip.open if str(ref_gff).endswith('.gz') else open
            with opener(ref_gff, 'rt') as fh:
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
                    tid_bare = tid.split('.')[0]
                    if tid_bare in keep_enst:   # version-stripped match
                        # Rewrite transcript_id to the version-stripped form so combined.gff,
                        # transcript_associations.txt, and the catalog all key ENST identically
                        # (our exemplar/TA ids are version-stripped; the reference is versioned).
                        out.write(line.replace(tid, tid_bare, 1))
                        written += 1
                        n_enst_records += 1
    log(f"[protein] wrote {combined_gz} ({n_novel:,} novel exemplars "
        f"[{n_from_store:,} via coord store, {n_from_scan:,} via GFF scan] + "
        f"{len(keep_enst):,} known ENST [{n_enst_records:,} ref records]; {written:,} lines)")

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

    # 4. run the EXISTING gff_translate -> protein_summary.txt (+ fastas) inside gff-output/.
    # The downstream comparisons.compute_differentials reads these from gff-output/ (PROTEIN_DIR), so
    # producer and consumer agree on the location. combined/ta are absolute paths, so cwd does not
    # matter for reading them.
    from .. import isoform_translation as isot
    cwd = os.getcwd()
    try:
        os.chdir(outdir)
        cds, transcripts, proteins = isot.gff_translate(combined_gz, genome_fasta, ref_gff, ta)
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
    return combined_gz, ta


# ------------------------------------------------------------------ STAGE 3 --
def stage3_rekey_h5ad(sample, h5ad_path, ta_path, kept_struct2exemplar, outdir,
                      barcode_clusters=None, kept_soft=None):
    """Re-key one sample's per-read h5ad to final isoform ids via sparse grouping matrix.
    Memory-optimized: builds var->final only for surviving molecules (no full molecule->struct map).

    The re-link is STRUCTURE-based (the canonical, cross-sample-consistent key): each molecule maps to
    its structure, and the structure maps to the final isoform id. This catches every sample's
    molecules for a shared structure, regardless of which sample 'won' the id at Tier 2.

    kept_soft: optional EM structure-keyed soft map (gene -> {child_structure: {parent_structure:
      weight}}). When given, a molecule of a child structure is split FRACTIONALLY across the parent
      structures' final ids (the EM per-cell analogue of em_reads); without it (WTA) each molecule
      maps to its single final id (weight 1).

    barcode_clusters: optional mapping cell_barcode -> cluster; restricts to annotated cells and
      attaches obs['cluster'] for downstream pseudobulk / PSI.
    Returns (out_path, n_final, raw_reads, final_reads, peak_rss_mb, elapsed_s)."""
    import numpy as np
    import scipy.sparse as sp
    import anndata as ad
    import pandas as pd
    # Coerce path-like args to str: callers may pass pathlib.Path (e.g. argparse-resolved or
    # extraction-returned), but the string ops below (ta_path.endswith, basename().split('.h5ad'))
    # require plain strings.
    h5ad_path = str(h5ad_path); ta_path = str(ta_path)
    if outdir is not None:
        outdir = str(outdir)
    t = time.time()
    # (gene,structure) -> "gene:exemplar". Final isoform ids are gene-namespaced so the output is
    # directly consumable by pseudo_cluster_counts / compute_differentials (var.split(":")[0] = gene).
    s2e = {}
    for g, m in kept_struct2exemplar.items():
        for struct, ex in m.items():
            s2e[(g, struct)] = f"{g}:{ex}"
    # EM: (gene,child_structure) -> {final_id: weight}, resolving parent structures to their exemplar ids.
    s2soft = {}
    if kept_soft:
        for g, gm in kept_soft.items():
            for child_struct, pw in gm.items():
                fw = {}
                for parent_struct, w in pw.items():
                    fid = s2e.get((g, parent_struct))
                    if fid is not None:
                        fw[fid] = fw.get(fid, 0.0) + w
                if fw:
                    s2soft[(g, child_struct)] = fw

    # molecule var ('gene:mol') -> {final_id: weight}
    var2final = {}
    with open(ta_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            key = (p[0], p[2])
            if s2soft and key in s2soft:
                var2final[f"{p[0]}:{p[3]}"] = s2soft[key]
            else:
                ex = s2e.get(key)
                if ex is not None:
                    var2final[f"{p[0]}:{p[3]}"] = {ex: 1.0}
    del s2e, s2soft
    a = ad.read_h5ad(h5ad_path)
    final_ids = sorted({fid for wm in var2final.values() for fid in wm})
    fidx = {f: j for j, f in enumerate(final_ids)}
    rows = []
    cols = []
    data = []
    for i, v in enumerate(a.var_names):
        wm = var2final.get(v)
        if wm:
            for fid, w in wm.items():
                rows.append(i)
                cols.append(fidx[fid])
                data.append(w)
    # integer (weight all 1) for WTA -> int matrix; fractional for EM -> float matrix.
    is_frac = bool(kept_soft)
    G = sp.coo_matrix((np.asarray(data, dtype=(np.float64 if is_frac else np.int64)),
                       (rows, cols)), shape=(a.shape[1], len(final_ids))).tocsr()
    Xf = (a.X.astype(np.float64 if is_frac else np.int64)) @ G
    obs = a.obs.copy()

    # Attach cluster annotations (cellHarmony) and restrict to annotated barcodes, so the final
    # clean-isoform h5ad is directly groupable by cluster downstream (pseudobulk / PSI).
    if barcode_clusters is not None:
        cl = barcode_clusters
        if not isinstance(cl, pd.Series):
            cl = pd.Series(cl)
        cl.index = cl.index.astype(str)
        bc_index = obs.index.astype(str)
        mapped = bc_index.map(cl)
        n_fwd = int(pd.notna(mapped).sum())
        # Long-read cell barcodes are often the REVERSE COMPLEMENT of the (standard-orientation)
        # cellHarmony annotations (metadata reverse=TRUE). The junction path applies this via
        # reverse_complement_seq; the isoform re-key must do the same or EVERY barcode misses ->
        # obs['cluster'] all-NaN -> pseudobulk columns become 'nan.<sample>' and the cluster-group
        # filter drops all data columns. Auto-detect: if the forward orientation matches nothing (or
        # far fewer), retry with reverse-complemented barcodes and keep whichever matches more.
        if n_fwd < len(obs):
            from ..isoform_matrix import reverse_complement_seq
            def _rc(b):
                try:
                    return reverse_complement_seq(b)
                except Exception:
                    return b
            rc_index = pd.Index([_rc(b) for b in bc_index])
            mapped_rc = rc_index.map(cl)
            n_rev = int(pd.notna(mapped_rc).sum())
            if n_rev > n_fwd:
                print(f"[stage3:{sample}] barcodes are reverse-complemented vs annotations: "
                      f"{n_rev}/{len(obs)} match revcomp (vs {n_fwd}/{len(obs)} forward); using revcomp.")
                mapped = mapped_rc
        obs['cluster'] = mapped.to_numpy() if hasattr(mapped, 'to_numpy') else mapped
        # Cells WITHOUT a cluster assignment (no matching barcode in the annotation) must be EXCLUDED
        # downstream -- never carried as 'nan'. Always restrict to annotated cells. If NONE match, fail
        # loudly rather than silently emitting unannotated cells (the old behaviour wrote all cells
        # without restriction, which leaked nan-cluster cells into the pseudobulk).
        keep = pd.notna(obs['cluster']).to_numpy()
        n_keep = int(keep.sum())
        if n_keep == 0:
            raise ValueError(
                f"[stage3:{sample}] 0/{len(obs)} barcodes matched cluster annotations (tried forward "
                f"AND reverse-complement). Refusing to write unannotated cells. Check that the "
                f"--cell_annot sample name matches the library '{sample}' and that its barcodes "
                f"overlap this sample's h5ad.")
        if n_keep < len(obs):
            print(f"[stage3:{sample}] excluding {len(obs) - n_keep}/{len(obs)} cells with no cluster "
                  f"assignment (no matching barcode); keeping {n_keep} annotated cells.")
        Xf = Xf.tocsr()[keep]
        obs = obs.iloc[keep].copy()

    out = ad.AnnData(X=Xf.tocsr(), obs=obs, var=pd.DataFrame(index=final_ids))
    base = os.path.basename(h5ad_path).split('.h5ad')[0]
    # Strip a trailing "-isoform" so re-running on an already-named "<lib>-isoform.h5ad" does not
    # produce "<lib>-isoform-isoform.h5ad". Output is always "<library>-isoform.h5ad" -- the name the
    # downstream get_valid_h5ad expects (drop-in for exportConsensusIsoformMatrix's isoform matrix).
    if base.endswith('-isoform'):
        base = base[:-len('-isoform')]
    out_path = os.path.join(outdir or os.path.dirname(h5ad_path), f"{base}-isoform.h5ad")
    out.write_h5ad(out_path, compression='gzip')
    return (out_path, len(final_ids), int(a.X.sum()), int(Xf.sum()), _self_peak_mb(), time.time() - t)


# ------------------------------------------------------------------ DRIVER ---
def run_pipeline(samples, outdir, nproc=8, min_total=3, ref=DEFAULT_REF, write_h5ad=True,
                 genome_fasta=None, ref_gff=None, enst_cache=None, barcode_clusters=None,
                 collapse_method='wta', log=print):
    """End-to-end. samples: list of (name, h5ad_path, transcript_associations_path[, raw_gff_path]).
    If genome_fasta is given AND samples include a raw_gff_path, protein prediction is run
    (final-isoform records -> existing gff_translate -> protein_summary.txt + *sequences.fasta).
    If enst_cache (ENST_reference_structures.tsv) is given, Ensembl transcripts are injected and
    exactly-matching observed isoforms are named/flagged by their ENST (known).
    barcode_clusters: optional {sample_name: {barcode: cluster}} so each re-keyed -isoform.h5ad
    carries obs['cluster'] (restricted to annotated cells) for downstream pseudobulk / PSI."""
    os.makedirs(outdir, exist_ok=True)
    sample_ta = [(s[0], s[2]) for s in samples]
    sample_gff_paths = {s[0]: s[3] for s in samples if len(s) > 3}
    gene_results, t1, max_rss, nparts = stage1_collapse(sample_ta, nproc=nproc, ref=ref,
                                                        enst_cache=enst_cache,
                                                        collapse_method=collapse_method)
    n_final = sum(len(r['exemplar_total']) for r in gene_results.values())
    log(f"[stage1] genes={len(gene_results):,} final={n_final:,} {t1:.1f}s maxRSS={max_rss:.0f}MB parts={nparts}")

    catalog, kept, kept_soft = stage2_outliers(gene_results, min_total=min_total, return_soft=True)
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

    # Persist the EM soft-weight map so a SEPARATE per-sample re-key job (4-phase cluster mode) can
    # reproduce the fractional EM allocation -- otherwise these weights live only in this process and
    # are lost when the collapse job ends, leaving EM impossible to re-key downstream. WTA needs only
    # FINAL_structure_to_exemplar.tsv (already written); EM additionally needs this. Mirrors the
    # in-memory kept_soft structure exactly: gene -> {child_structure: {parent_structure: weight}}.
    # Empty for WTA (no soft map) -- the file is still written (header only) so its presence/absence is
    # unambiguous and the loader is uniform.
    softpath = os.path.join(outdir, "FINAL_structure_to_exemplar_soft.tsv")
    with open(softpath, 'w') as o:
        o.write("gene\tchild_structure\tparent_structure\tweight\n")
        for g, gm in (kept_soft or {}).items():
            for child_struct, pw in gm.items():
                for parent_struct, w in pw.items():
                    o.write(f"{g}\t{child_struct}\t{parent_struct}\t{w:.10g}\n")
    if kept_soft:
        log(f"[stage2] wrote EM soft map {softpath} ({sum(len(v) for v in kept_soft.values()):,} child structures)")

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
            bc = (barcode_clusters or {}).get(name)
            r = stage3_rekey_h5ad(name, h5, ta, kept, outdir=None, barcode_clusters=bc,
                                  kept_soft=(kept_soft if collapse_method == 'em' else None))
            raw_reads, final_reads = r[2], r[3]
            kept_frac = final_reads / raw_reads if raw_reads else 0.0
            log(f"[stage3:{name}] -> {r[1]:,} final isoforms, reads {raw_reads:,}->{final_reads:,} "
                f"({100*kept_frac:.0f}%), {r[5]:.1f}s peakRSS={r[4]:.0f}MB")
            # LOGICAL VALIDATION (read conservation): collapsing molecule reads onto final isoforms
            # must conserve reads except for the deliberate <min_total outlier removal + cells dropped
            # for lacking a cluster annotation. So final_reads must be (a) never MORE than raw, and
            # (b) not collapse to a tiny fraction (which would signal a barcode/structure mapping bug,
            # like the -1 mismatch that once produced 0% overlap). Warn loudly if violated.
            if final_reads > raw_reads:
                log(f"[stage3:{name}] *** VALIDATION FAIL: final reads ({final_reads:,}) exceed raw "
                    f"({raw_reads:,}) -- impossible; mapping is double-counting.")
            elif raw_reads > 0 and kept_frac < 0.5:
                log(f"[stage3:{name}] *** VALIDATION WARNING: only {100*kept_frac:.1f}% of reads "
                    f"mapped to final isoforms. Expected ~> the non-outlier fraction. Check barcode "
                    f"convention (molecule h5ad vs annotation) and structure->exemplar mapping.")
            h5ad_results.append((name,) + r)

        # Aggregate read-conservation validation across all samples.
        tot_raw = sum(r[3] for r in h5ad_results)
        tot_final = sum(r[4] for r in h5ad_results)
        if tot_raw:
            log(f"[validate] read conservation: {tot_final:,}/{tot_raw:,} "
                f"({100*tot_final/tot_raw:.1f}%) of molecule reads retained in final isoform h5ads "
                f"(remainder = <{min_total}-read outliers + unannotated cells)")
    # ---- workflow summary tables + bar charts (summary/ folder) ----
    # All numbers come from the live collapse results above + per-sample artifacts on disk.
    try:
        from ..summary_report import write_summary
        # per-sample junction h5ad sits next to each sample's molecule h5ad as "<stem>-junction.h5ad";
        # reads come from each sample's transcript_associations.txt (already used for collapse).
        jh5, tap = {}, {}
        for s in samples:
            name, h5, ta = s[0], str(s[1]), s[2]
            stem = h5[:-len('.h5ad')] if h5.endswith('.h5ad') else h5
            cand = stem + '-junction.h5ad'
            if os.path.exists(cand):
                jh5[name] = cand
            if ta and os.path.exists(str(ta)):
                tap[name] = str(ta)
        psum = os.path.join(outdir, 'protein_summary.txt')   # outdir is gff-output/
        # Write the summary/ folder at the RUN-DIR level (sibling of gff-output) for discoverability.
        summary_root = os.path.dirname(os.path.abspath(outdir)) or outdir
        write_summary(summary_root, collapse_method=collapse_method, min_total=min_total,
                      n_structures=n_final, catalog=catalog, h5ad_results=h5ad_results,
                      sample_junction_h5ads=jh5, sample_ta_paths=tap,
                      protein_summary_path=psum, log=log)
    except Exception as e:
        log(f"[summary] table generation skipped ({type(e).__name__}: {e})")

    return dict(gene_results=gene_results, catalog=catalog, kept=kept, protein=protein,
                stage1=(t1, max_rss, nparts), h5ad=h5ad_results)


# ---------------------------------------------------------------------------
# 4-PHASE CLUSTER SUPPORT: build the catalog once (P3), then re-key each sample
# in a SEPARATE job (P4) by loading the persisted maps from disk. This is the
# memory-efficient split -- nothing from the collapse process needs to stay
# resident; the structure->final-id map (WTA) and the EM soft-weight map are
# read back from the two FINAL_structure_to_exemplar*.tsv files.
# ---------------------------------------------------------------------------

def load_struct2exemplar(outdir):
    """Load FINAL_structure_to_exemplar.tsv -> kept {gene: {structure: final_id}} (the WTA map)."""
    path = os.path.join(outdir, "FINAL_structure_to_exemplar.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Collapse catalog map missing: {path}. Run the collapse (P3) first.")
    kept = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            kept.setdefault(p[0], {})[p[1]] = p[2]
    return kept


def load_struct2exemplar_soft(outdir):
    """Load FINAL_structure_to_exemplar_soft.tsv -> kept_soft
    {gene: {child_structure: {parent_structure: weight}}}. Returns None if the file is absent or has
    no rows (i.e. a WTA collapse), so callers fall back to the hard WTA map automatically."""
    path = os.path.join(outdir, "FINAL_structure_to_exemplar_soft.tsv")
    if not os.path.exists(path):
        return None
    kept_soft = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            g, child, parent, w = p[0], p[1], p[2], p[3]
            kept_soft.setdefault(g, {}).setdefault(child, {})[parent] = float(w)
    return kept_soft or None


def rekey_one_sample(name, h5ad_path, ta_path, outdir, barcode_clusters=None,
                     collapse_method='wta', write_dir=None, log=print):
    """P4 per-sample entry: re-key ONE sample's molecule h5ad onto the final catalog, loading the
    structure->final-id map (and EM soft weights) from disk -- no dependency on the collapse process's
    memory. Writes <h5ad_path stem>-isoform.h5ad next to the sample (stage3's default).

    collapse_method: 'wta' uses only the hard map; 'em' additionally loads the persisted soft weights
    and splits each molecule fractionally (identical to the in-process EM re-key). If 'em' is asked
    but no soft map exists on disk, it falls back to WTA (and warns) rather than failing."""
    outdir = str(outdir); h5ad_path = str(h5ad_path); ta_path = str(ta_path)
    if write_dir is not None:
        write_dir = str(write_dir)
    kept = load_struct2exemplar(outdir)
    kept_soft = None
    if collapse_method == 'em':
        kept_soft = load_struct2exemplar_soft(outdir)
        if kept_soft is None:
            log(f"[rekey:{name}] WARNING: --collapse_method em but no EM soft map on disk "
                f"({outdir}/FINAL_structure_to_exemplar_soft.tsv); falling back to WTA re-key.")
    r = stage3_rekey_h5ad(name, h5ad_path, ta_path, kept, outdir=write_dir,
                          barcode_clusters=barcode_clusters, kept_soft=kept_soft)
    raw_reads, final_reads = r[2], r[3]
    kept_frac = final_reads / raw_reads if raw_reads else 0.0
    log(f"[rekey:{name}] -> {r[1]:,} final isoforms, reads {raw_reads:,}->{final_reads:,} "
        f"({100*kept_frac:.0f}%), {r[5]:.1f}s peakRSS={r[4]:.0f}MB")
    if final_reads > raw_reads:
        log(f"[rekey:{name}] *** VALIDATION FAIL: final reads ({final_reads:,}) exceed raw ({raw_reads:,}).")
    elif raw_reads > 0 and kept_frac < 0.5:
        log(f"[rekey:{name}] *** VALIDATION WARNING: only {100*kept_frac:.1f}% of reads mapped.")
    return (name,) + r
