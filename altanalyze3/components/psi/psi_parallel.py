"""Chromosome-parallel PSI with output identical to the serial psi_single.main.

psi_single.main (1) drops junctions whose max count over all columns is below min_read, (2) stable-sorts
junctions by gene (text before the first ':'), (3) computes PSI one gene at a time; the overlap clique, PSI
and event filters never look outside the gene. So the work splits exactly into chunks of WHOLE genes.

This module:
  * reads the junction pseudobulk h5ad once and slices it into chunks by CHROMOSOME (every junction of a
    gene goes to the chromosome of that gene's first junction); a chromosome holding more than 1/n_jobs of
    the junctions is sub-split into whole-gene pieces so one large chromosome does not set the run time;
  * writes each chunk as an h5ad whose obs (the sample columns) is untouched -- never rebuilt with concat;
  * runs the UNMODIFIED psi_single.main on every chunk in parallel;
  * verifies every chunk header equals the input sample columns (names and order), every row has
    1 + n_samples fields, and no gene appears in two chunks; any failure raises;
  * merges chunk rows back in the serial gene-sorted order, so the file matches the serial file.

run_psi_parallel(junction_h5ad, outdir, min_read=None, n_jobs=None, **psi_kwargs) -> outdir
"""
import asyncio
import os
import sys
import time

import numpy as np

from . import psi_single as psi


def _log(m):
    print(f"[psi-parallel {time.strftime('%H:%M:%S')}] {m}", flush=True)


def _run_chunk(args):
    """Worker: the serial PSI on one chunk h5ad."""
    idx, chunk_h5ad, chunk_txt, min_read, psi_kwargs = args
    t = time.time()
    asyncio.run(psi.main(junction_path=chunk_h5ad, query_gene=None, outdir=chunk_txt, min_read=min_read,
                         **psi_kwargs))
    return idx, chunk_txt, time.time() - t


def _chunks_by_chromosome(features, n_jobs, keep_mask):
    """Assign whole genes to chunks by chromosome; sub-split large chromosomes by gene. Returns
    (chunk_of_feature array, chunk labels)."""
    genes = np.array([f.split(':', 1)[0] for f in features])
    chroms = np.array([f.split('=', 1)[1].split(':', 1)[0] if '=' in f else 'NA' for f in features])
    # a gene's chromosome = chromosome of its first junction (genes do not span chromosomes; if one did,
    # keeping it whole preserves the serial result, which groups by gene only)
    first_chrom = {}
    for g, c in zip(genes, chroms):
        first_chrom.setdefault(g, c)
    gene_chrom = np.array([first_chrom[g] for g in genes])
    n_split_genes = len({g for g, c in zip(genes, chroms) if c != first_chrom[g]})
    if n_split_genes:
        _log(f"{n_split_genes} gene(s) have junctions on >1 chromosome; kept whole in their first chromosome's chunk")
    work = keep_mask.sum()
    cap = max(1, int(np.ceil(work / max(n_jobs, 1))))
    chunk_of = np.empty(len(features), dtype=np.int64)
    labels = []
    for c in sorted(set(gene_chrom)):
        idx = np.where(gene_chrom == c)[0]
        ugenes = sorted(set(genes[idx]))
        per_gene = {}
        for i in idx:
            per_gene[genes[i]] = per_gene.get(genes[i], 0) + int(keep_mask[i])
        part, load, gene_part = 0, 0, {}
        for g in ugenes:                                   # whole genes, sorted, greedy to the cap
            if load >= cap and load > 0:
                part, load = part + 1, 0
            gene_part[g] = part
            load += per_gene[g]
        for p in range(part + 1):
            labels.append(f"{c}" if part == 0 else f"{c}.part{p + 1}")
        base = len(labels) - (part + 1)
        for i in idx:
            chunk_of[i] = base + gene_part[genes[i]]
    return chunk_of, labels, genes


def _default_jobs():
    """Workers: ALTANALYZE3_PSI_JOBS, else the LSF slot count (LSB_DJOB_NUMPROC), else the CPUs this process may
    use (sched_getaffinity), else cpu_count; capped at 16 so a small allocation is never oversubscribed."""
    for key in ("ALTANALYZE3_PSI_JOBS", "LSB_DJOB_NUMPROC"):
        try:
            v = int(os.environ.get(key, "0"))
        except ValueError:
            v = 0
        if v > 0:
            return min(v, 16)
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    return max(1, min(n, 16))


def run_psi_parallel(junction_h5ad, outdir, min_read=None, n_jobs=None, **psi_kwargs):
    """Chromosome-parallel equivalent of asyncio.run(psi_single.main(junction_h5ad, None, outdir, min_read)).
    outdir = the PSI TSV path (same meaning as psi_single.main). The chunk h5ads and per-chunk PSI files are
    KEPT in a new folder '<outdir>.chunks_<timestamp>' (never deleted by this code) so every merged row can be
    traced back to its chunk."""
    import anndata as ad
    from multiprocessing import get_context
    T0 = time.time()
    n_jobs = n_jobs or _default_jobs()
    work = f"{os.path.abspath(outdir)}.chunks_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    os.makedirs(work)                                      # new, unique folder; nothing pre-existing is touched

    adata = ad.read_h5ad(junction_h5ad)
    samples = [str(s) for s in adata.obs_names]
    features = np.array([str(f) for f in adata.var_names])
    if len(set(samples)) != len(samples):
        raise ValueError("duplicate sample column names in the junction h5ad; PSI columns would be ambiguous")
    X = adata.X
    colmax = np.asarray(X.max(axis=0).todense()).ravel() if hasattr(X, "todense") else np.nanmax(np.asarray(X), axis=0)
    keep = colmax > (min_read - 1) if min_read is not None else np.ones(len(features), bool)
    kept_genes = {f.split(':', 1)[0] for f in features[keep]}
    if len(kept_genes) <= 1:
        # psi_single writes a second header line when the whole input is one gene; run it serially so the
        # output stays identical to the serial path in that edge case.
        _log(f"{len(kept_genes)} gene(s) pass min_read -> serial psi_single (identical by definition)")
        asyncio.run(psi.main(junction_path=junction_h5ad, query_gene=None, outdir=outdir, min_read=min_read, **psi_kwargs))
        return outdir
    chunk_of, labels, genes = _chunks_by_chromosome(features, n_jobs, keep)
    _log(f"{len(features):,} junctions ({int(keep.sum()):,} pass min_read={min_read}), {len(samples):,} sample columns, "
         f"{len(labels)} chunks, {n_jobs} workers")
    jobs = []
    for k, lab in enumerate(labels):
        idx = np.where(chunk_of == k)[0]                   # ascending = original var order
        if len(idx) == 0:
            continue
        sub = adata[:, idx]
        sub = ad.AnnData(X=sub.X.copy(), obs=adata.obs.copy(), var=sub.var.copy())  # obs = the input columns, untouched
        if list(map(str, sub.obs_names)) != samples:
            raise RuntimeError(f"chunk {lab}: sample columns changed while slicing")
        p = f"{work}/chunk_{k:03d}_{lab}.h5ad"
        sub.write_h5ad(p)
        jobs.append((k, p, f"{work}/psi_{k:03d}_{lab}.txt", min_read, psi_kwargs))
    del adata, X
    _log(f"split into {len(jobs)} chunk files in {time.time() - T0:.1f}s")

    T1 = time.time()
    ctx = get_context("fork" if sys.platform.startswith("linux") else "spawn")
    with ctx.Pool(min(n_jobs, len(jobs))) as pool:
        results = sorted(pool.map(_run_chunk, jobs))
    _log(f"PSI on {len(results)} chunks in {time.time() - T1:.1f}s (slowest {max(r[2] for r in results):.1f}s)")

    # ---- integrity checks + merge in the serial gene order ----
    expected_header = "\t".join(samples) + "\n"
    # psi_single writes a SECOND header ('uid-bg_uid' + samples) inside a chunk that holds a single gene;
    # verify it carries the input columns in order, then drop it (a multi-gene serial run never writes it).
    secondary_header = "uid-bg_uid\t" + "\t".join(samples) + "\n"
    n_secondary = 0
    ncol = 1 + len(samples)
    gene_blocks, gene_owner, n_rows = {}, {}, 0
    for k, txt, _ in results:
        with open(txt) as fh:
            h = fh.readline()
            if h != expected_header:
                raise RuntimeError(f"{txt}: header differs from the input sample columns (names/order)")
            for line in fh:
                if line == expected_header:
                    raise RuntimeError(f"{txt}: header repeated inside the chunk")
                if line.startswith("uid-bg_uid\t"):
                    if line != secondary_header:
                        raise RuntimeError(f"{txt}: psi_single secondary header differs from the input sample columns")
                    n_secondary += 1
                    continue
                if line.count("\t") != ncol - 1:
                    raise RuntimeError(f"{txt}: row with {line.count(chr(9)) + 1} fields, expected {ncol}")
                g = line.split(":", 1)[0]
                if gene_owner.setdefault(g, k) != k:
                    raise RuntimeError(f"gene {g} produced rows in two chunks ({gene_owner[g]} and {k})")
                gene_blocks.setdefault(g, []).append(line)
                n_rows += 1
    with open(outdir, "w") as fo:
        fo.write(expected_header)
        for g in sorted(gene_blocks):                      # same order as psi_single's stable gene sort
            fo.writelines(gene_blocks[g])
    if n_secondary:
        _log(f"dropped {n_secondary} psi_single secondary header line(s) from single-gene chunks (columns verified)")
    _log(f"merged {n_rows:,} PSI rows from {len(gene_blocks):,} genes -> {outdir}; chunks kept in {work}; "
         f"total {time.time() - T0:.1f}s")
    return outdir
