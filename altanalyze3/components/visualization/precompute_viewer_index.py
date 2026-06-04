#!/usr/bin/env python3
"""Pre-build the isoform_structure_viewer's per-cell-type counts caches so ANY query composes
instantly with NO rebuild.

WHY: the viewer keys its isoform-counts cache on the *exact* cell_states selection
(``_counts_cache_key``). Every new cell-type / sample-group combination is therefore a cache miss
-> a full per-sample h5ad scan (the lag). The viewer ALSO has a composer
(``_try_compose_isoform_counts_cache`` + ``_find_disjoint_cover``) that builds any requested
selection by SUMMING existing per-subset caches -- but nothing pre-builds the atomic pieces.

THIS SCRIPT builds, for every sample h5ad, one ATOMIC counts cache per INDIVIDUAL cell type. Once
present, every ``plot_isoform_structures_by_conditions`` call -- for any ``cell_states`` selection
AND any ``conditions`` (sample group) -- is answered by composing those atomic caches, with no h5ad
re-read. (Validated: composed counts == direct mask+sum, 0 mismatches.)

Atomic granularity is per (cell_type x sample h5ad): the viewer iterates samples for ``conditions``
and composes cell types for ``cell_states``, so per-(cell_type x sample) pieces cover any
combination of both.

Usage (mirrors your existing viewer setup -- pass the SAME sample_dict + barcode_sample_dict +
groupby_rev you use for plotting):

    from altanalyze3.components.visualization import precompute_viewer_index as pvi
    pvi.precompute(sample_dict, barcode_sample_dict, groupby_rev=True)

Idempotent: an atomic cache that already exists (and is not stale) is skipped, so re-running only
fills gaps. Safe to run repeatedly / after adding samples.
"""
import os
import time

import anndata as ad
import pandas as pd

from . import isoform_structure_view as isv


def _sample_name(sample, h5ad_path):
    # Prefer an explicit library/gff_name field (this is the sample key used by
    # barcode_sample_dict). Fall back to the h5ad basename, stripping the '-isoform' suffix that
    # export_sample_isoform appends to the molecule h5ad ('<stem>-isoform.h5ad') so the inferred
    # name still matches the barcode annotation key '<stem>'.
    for key in ("library", "gff_name"):
        val = sample.get(key) if isinstance(sample, dict) else None
        if val:
            return str(val)
    base = os.path.basename(str(h5ad_path))
    name = base.split(".h5ad")[0]
    if name.endswith("-isoform"):
        name = name[: -len("-isoform")]
    return name


def _detect_orientation(adata, barcode_series, cell_types, probe=8):
    """Return the reverse_complement flag (True/False) whose masks actually match cells for this
    sample, or None if NEITHER orientation matches anything. Probes up to ``probe`` cell types and
    sums matched cells under each orientation; picks the larger. Long-read h5ads are typically
    reverse-complemented vs the annotation, short-read are forward -- this avoids hardcoding it."""
    fwd = rev = 0
    for ct in list(cell_types)[:probe]:
        mf = isv._load_barcode_mask(adata, barcode_series, [ct], reverse_complement=False)
        mr = isv._load_barcode_mask(adata, barcode_series, [ct], reverse_complement=True)
        fwd += int(mf.sum()) if mf is not None else 0
        rev += int(mr.sum()) if mr is not None else 0
    if fwd == 0 and rev == 0:
        return None
    return rev > fwd


def _build_all_celltype_caches_one_pass(adata, barcode_series, cell_types, h5ad_path,
                                        groupby_rev, mask_rev, groupby_sample=None, log=print):
    """Build ALL per-cell-type atomic caches in ONE pass over the matrix (instead of N mask+sum
    passes). Sums isoform counts per cell type via a single sparse cells x cell_types group matrix
    multiply, then writes each cell type's column as its atomic cache (keyed for groupby_rev so the
    viewer composes from them). Returns the number of caches written.

    Equivalent to calling _build_isoform_counts_cache once per cell type, but ~N times faster: the
    big X matrix is read once and reduced by a single G^T @ X.
    """
    import numpy as np
    import scipy.sparse as sp

    obs_names = np.asarray(adata.obs_names, dtype=str)
    # map each obs barcode -> its cell type, applying the detected orientation to the annotation side
    bser = barcode_series
    if not isinstance(bser, pd.Series):
        bser = pd.Series(bser)
    bser = bser.dropna()
    bser.index = bser.index.astype(str)
    if mask_rev:
        # annotation barcodes are reverse-complemented relative to obs -> rc the annotation keys
        rc = {isv._reverse_complement_barcode(b): v for b, v in bser.items()}
        ann = pd.Series(rc)
    else:
        ann = bser
    ct_index = {ct: j for j, ct in enumerate(cell_types)}
    rows, cols = [], []
    for i, bc in enumerate(obs_names):
        ct = ann.get(bc)
        if ct is not None and ct in ct_index:
            rows.append(i)
            cols.append(ct_index[ct])
    if not rows:
        return 0
    n_cells = adata.n_obs
    G = sp.csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)),
                      shape=(n_cells, len(cell_types)))
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    # cell_types x isoforms = G^T @ X. Keep CT SPARSE (CSR) -- densifying it would be
    # n_cell_types x n_isoforms (58 x 35M = ~8GB) which is the bloat we are eliminating. We densify
    # only ONE cell type's row at a time inside the loop, transiently, then free it.
    CT = (G.T @ X).tocsr()
    var_names = np.asarray(adata.var_names, dtype=str)
    written = 0
    for ct, j in ct_index.items():
        row = CT.getrow(j)
        if row.nnz == 0:
            continue  # no cells of this type -> nothing to cache (composition omits it)
        # IMPORTANT: key by the SAME groupby_sample the viewer passes (sample_name) or the compose
        # path filters every candidate out and the viewer rebuilds instead of composing.
        cache_path = isv._h5ad_counts_cache_path(str(h5ad_path), "cell_states", [ct], groupby_rev, groupby_sample)
        # idempotent: skip if a (sparse or dense) cache already exists and is current
        if isv._counts_cache_ready(cache_path, [str(h5ad_path)]):
            continue
        col = np.asarray(row.todense()).ravel().astype(np.float32)  # transient dense, freed each iter
        isv._save_isoform_counts_cache_arrays(cache_path, col, var_names)  # writes SPARSE on disk
        isv._write_isoform_counts_cache_metadata(cache_path, str(h5ad_path), "cell_states", [ct],
                                                 groupby_rev, groupby_sample)
        del col
        written += 1
    return written


def _prebuild_structure_indexes(h5ad_path, sname, log=print):
    """Pre-build the TWO non-counts indexes the viewer needs so the FIRST render is not a slow 'cold'
    build: (a) the h5ad gene-index (.gene_index.npz), and (b) the per-sample transcripts.db (the
    transcript-structure SQLite, built from the multi-GB transcript_associations.txt -- the dominant
    cold cost). Both are idempotent (skipped if current). With these prebuilt, cold render == warm."""
    try:
        tg = time.time()
        gi_path = isv._h5ad_gene_index_path(h5ad_path)
        if isv._is_index_stale(gi_path, [h5ad_path]):
            isv._build_h5ad_gene_index(h5ad_path, gi_path)
        # var_names sidecar: reading the var index out of a large h5ad (e.g. 35M isoforms) is a ~10s
        # HDF5 read paid on every render that loads a counts cache. Persist it once here so later
        # renders (any process) read a small .npy instead -- this was the dominant remaining cold cost.
        isv._build_var_names_sidecar(h5ad_path)
        assoc = os.path.join(os.path.dirname(h5ad_path), 'gff-output', 'transcript_associations.txt')
        assoc = isv._resolve_transcript_associations_path(assoc) if os.path.exists(assoc) else None
        if assoc:
            isv._ensure_transcripts_db(assoc)   # builds <assoc_dir>/gene_indexes_v2/transcripts.db if stale
        log(f"[precompute] {sname}: gene-index + var_names + transcripts.db ready ({time.time()-tg:.1f}s)")
    except Exception as _e:
        log(f"[precompute] {sname}: structure-index prebuild skipped ({type(_e).__name__}: {_e}); "
            f"first render will build it")


def _cell_types_for(barcode_series):
    """Distinct cell-type labels present for a sample (the atomic units to pre-build)."""
    if barcode_series is None:
        return []
    s = barcode_series
    if not isinstance(s, pd.Series):
        s = pd.Series(s)
    vals = pd.unique(s.dropna().astype(str))
    return [v for v in vals if v and v.lower() != "nan"]


def precompute(sample_dict, barcode_sample_dict, groupby_rev=True, index_dir=None,
               only_samples=None, log=print):
    """Build atomic per-cell-type isoform-counts caches for every sample h5ad.

    sample_dict          : as used by the viewer (uid -> [sample dicts with 'matrix','groups',...]).
    barcode_sample_dict  : sample_name -> Series(barcode -> cell_type) (iso.import_barcode_clusters).
    groupby_rev          : MUST match what you pass to plot_isoform_structures_by_conditions
                           (viewer default True -- long-read barcodes are reverse-complemented).
    only_samples         : optional iterable of sample_names to restrict to.
    Returns a summary dict {sample_name: n_cell_types_built}.
    """
    built = {}
    t0 = time.time()
    for uid, samples in sample_dict.items():
        for sample in samples:
            h5ad_path = sample.get("matrix")
            if not h5ad_path or not os.path.exists(str(h5ad_path)):
                log(f"[precompute] skip (no matrix h5ad): uid={uid}")
                continue
            sname = _sample_name(sample, h5ad_path)
            if only_samples and sname not in set(only_samples):
                continue
            bser = barcode_sample_dict.get(sname)
            if bser is None:
                log(f"[precompute] {sname}: no barcode->cell_type annotation; skipping")
                continue
            cell_types = _cell_types_for(bser)
            if not cell_types:
                log(f"[precompute] {sname}: 0 cell types in annotation; skipping")
                continue

            # Pre-build the structure indexes (gene-index + transcripts.db) FIRST, before the counts
            # early-skip -- so a sample whose counts caches exist but whose transcripts.db does not still
            # gets it (and so it is built regardless of the counts path). Both are idempotent.
            _prebuild_structure_indexes(str(h5ad_path), sname, log)

            # EARLY SKIP (index-once): if EVERY cell type's atomic cache already exists and is
            # current, this sample is fully indexed -- do NOT load the h5ad or recompute the matmul.
            # Without this, re-running precompute (or the workflow) re-reads the big h5ad and redoes
            # the whole G^T@X every time, only to discard it at the per-cell-type skip -- i.e. it
            # "re-indexes every time". The index must be built ONCE and then be a near-instant no-op.
            def _cache_current(ct):
                cp = isv._h5ad_counts_cache_path(str(h5ad_path), "cell_states", [ct], groupby_rev, sname)
                return isv._counts_cache_ready(cp, [str(h5ad_path)])  # sparse OR dense
            if all(_cache_current(ct) for ct in cell_types):
                log(f"[precompute] {sname}: index already complete ({len(cell_types)} cell types); "
                    f"skipping (no h5ad load).")
                built[sname] = 0
                continue

            # Open once (backed) and build one atomic cache per cell type. Skip cell types whose
            # cache already exists & is current (idempotent).
            adata = None
            n_built = 0
            ts = time.time()
            try:
                adata = ad.read_h5ad(str(h5ad_path))  # full load (masks need it)
                # AUTO-DETECT barcode orientation for THIS sample: the h5ad barcodes may match the
                # annotation forward OR reverse-complemented (long-read vs short-read convention). We
                # build masks with whichever orientation actually matches cells -- otherwise every
                # cell type masks 0 cells and nothing is cached. The cache is KEYED by the viewer's
                # requested ``groupby_rev`` (so the viewer's compose path finds it), but the mask uses
                # the DETECTED orientation so the counts are real.
                mask_rev = _detect_orientation(adata, bser, cell_types)
                if mask_rev is None:
                    log(f"[precompute] {sname}: 0 barcodes match in EITHER orientation; skipping "
                        f"(check that the annotation sample name matches the library).")
                    built[sname] = 0
                    continue
                if mask_rev != groupby_rev:
                    log(f"[precompute] {sname}: data barcode orientation (rev={mask_rev}) differs from "
                        f"requested groupby_rev={groupby_rev}; building masks with the matching "
                        f"orientation, keyed for groupby_rev={groupby_rev}. Pass groupby_rev={mask_rev} "
                        f"to the viewer for this sample, OR ensure consistency.")
                # ONE pass over the matrix builds ALL per-cell-type atomic caches (G^T @ X), instead
                # of N separate mask+sum passes. Same cache files, ~N times faster.
                # groupby_sample = sname: the viewer keys its cache requests by sample_name, so the
                # atomic caches MUST carry the same value or the compose path won't match them.
                n_built = _build_all_celltype_caches_one_pass(
                    adata, bser, cell_types, str(h5ad_path), groupby_rev, mask_rev,
                    groupby_sample=sname, log=log)
            finally:
                del adata
            built[sname] = n_built
            log(f"[precompute] {sname}: built {n_built}/{len(cell_types)} atomic cell-type caches "
                f"({time.time() - ts:.1f}s)")
    log(f"[precompute] done: {sum(built.values())} atomic caches across {len(built)} samples "
        f"({time.time() - t0:.1f}s). Viewer queries now compose without rebuilds.")
    return built


if __name__ == "__main__":
    raise SystemExit(
        "Import and call precompute(sample_dict, barcode_sample_dict, groupby_rev=...). "
        "Build sample_dict via isoform_automate.import_metadata and barcode_sample_dict via "
        "isoform_matrix.import_barcode_clusters, exactly as for the viewer.")
