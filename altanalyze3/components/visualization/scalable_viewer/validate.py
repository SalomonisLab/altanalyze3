"""Prove a bundle reproduces its source h5ad.

Run this after every precompute. It re-reads the source h5ad and checks the bundle
against it. Any failure raises; nothing is silently repaired.

  PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
  /opt/homebrew/opt/python@3.11/bin/python3.11 \
    -m altanalyze3.components.visualization.scalable_viewer.validate \
    --bundle /path/to/bundles/MyDataset --prefix My-Dataset [--n-cells 200] [--n-genes 40]

Checks, each reported with its denominator:
  A  nnz conservation between the source layer and the gene-major store
  B  every non-zero of N random cells is present in the store with the same value
  C  cell indices ascend strictly inside each gene
  D  stats_mean / stats_frac recomputed from the store match the stored matrices
  E  stats_n matches the cell-state counts in obs
  F  the embedding is finite and has one row per cell
  G  the legacy TSVs have exactly one row per cell or per gene
  H  the centroid matrix columns are the canonical cell-state order
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import h5py
import numpy as np

from . import bundle as B


def _fail(msg: str) -> None:
    raise AssertionError(msg)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Validate a scalable_viewer bundle against its h5ad.")
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--h5ad", default=None, help="default: source_h5ad from the metadata")
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--n-genes", type=int, default=40)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)

    p = B.BundlePaths(a.bundle, a.prefix)
    missing = p.missing()
    if missing:
        _fail(f"bundle is incomplete: {missing}")
    meta = B.read_metadata(p.metadata)
    sv = B.viewer_block(meta) or _fail("metadata has no scalable_viewer block")
    src = a.h5ad or meta["source_h5ad"]
    if not os.path.isfile(src):
        _fail(f"source h5ad not found: {src}")

    print(f"bundle : {p.bundle_dir}")
    print(f"prefix : {a.prefix}")
    print(f"source : {src}")
    print(f"built  : {sv.get('built_utc')}  layer={sv.get('layer')}  dtype={sv.get('expr_dtype')}")
    print("-" * 78)

    f = h5py.File(src, "r")
    layer = sv.get("layer", "lognorm")
    grp = f["X"] if layer == "X" else f["layers"][layer]
    sip = grp["indptr"][:].astype(np.int64)
    indptr = np.load(p.expr_indptr)
    ind = np.load(p.expr_indices, mmap_mode="r")
    val = np.load(p.expr_data, mmap_mode="r")
    cells = np.load(p.cells)
    sc = cells["state_code"].astype(np.int64)
    emb = cells["embedding"]
    sn = np.load(p.stats_n)
    mean = np.load(p.stats_mean)
    frac = np.load(p.stats_frac)
    states = sv["states"]
    n_cells, n_genes = len(sip) - 1, len(indptr) - 1
    S = len(states)
    rng = np.random.RandomState(a.seed)
    ok = True

    # A ---------------------------------------------------------------------
    src_nnz, store_nnz = int(sip[-1]), int(indptr[-1])
    good = src_nnz == store_nnz
    ok &= good
    print(f"A nnz conservation            : source {src_nnz:,} vs store {store_nnz:,}  "
          f"{'PASS' if good else 'FAIL'}")

    # B ---------------------------------------------------------------------
    t0 = time.time()
    sel = np.sort(rng.choice(n_cells, min(a.n_cells, n_cells), replace=False))
    nchk = nbad = 0
    maxerr = 0.0
    for c in sel:
        s, e = int(sip[c]), int(sip[c + 1])
        gs = grp["indices"][s:e].astype(np.int64)
        vs = np.asarray(grp["data"][s:e], dtype=np.float32)
        for gi, vv in zip(gs, vs):
            lo, hi = int(indptr[gi]), int(indptr[gi + 1])
            seg = np.asarray(ind[lo:hi])
            k = int(np.searchsorted(seg, c))
            nchk += 1
            if k >= seg.size or seg[k] != c:
                nbad += 1
                continue
            err = abs(float(val[lo + k]) - float(vv))
            maxerr = max(maxerr, err)
            if err > 0:
                nbad += 1
    good = nbad == 0
    ok &= good
    print(f"B store vs source, {len(sel)} cells : {nchk:,} non-zeros checked, {nbad} mismatched, "
          f"max abs error {maxerr:.3e}  [{time.time() - t0:.1f}s]  {'PASS' if good else 'FAIL'}")

    # C ---------------------------------------------------------------------
    gsel = rng.choice(n_genes, min(200, n_genes), replace=False)
    unsorted = 0
    for gi in gsel:
        lo, hi = int(indptr[gi]), int(indptr[gi + 1])
        seg = np.asarray(ind[lo:hi])
        if seg.size and not np.all(np.diff(seg) > 0):
            unsorted += 1
    good = unsorted == 0
    ok &= good
    print(f"C CSC sorted, {len(gsel)} genes      : {unsorted} genes not strictly ascending  "
          f"{'PASS' if good else 'FAIL'}")

    # D ---------------------------------------------------------------------
    me = fe = 0.0
    for gi in rng.choice(n_genes, min(a.n_genes, n_genes), replace=False):
        lo, hi = int(indptr[gi]), int(indptr[gi + 1])
        col = np.zeros(n_cells, dtype=np.float32)
        col[np.asarray(ind[lo:hi], dtype=np.int64)] = np.asarray(val[lo:hi], dtype=np.float32)
        s1 = np.bincount(sc, weights=col.astype(np.float64), minlength=S)
        c1 = np.bincount(sc, weights=(col > 0).astype(np.float64), minlength=S)
        me = max(me, float(np.abs(s1 / sn - mean[gi]).max()))
        fe = max(fe, float(np.abs(c1 / sn - frac[gi]).max()))
    good = me < 1e-4 and fe < 1e-4
    ok &= good
    print(f"D statistics recomputed        : max |mean diff| {me:.3e}, max |frac diff| {fe:.3e}  "
          f"{'PASS' if good else 'FAIL'}")

    # E ---------------------------------------------------------------------
    node = f["obs"][meta["cluster_key"]]
    cats = [x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]]
    codes = node["codes"][:]
    src_n = np.array([int((codes == cats.index(s)).sum()) for s in states])
    good = np.array_equal(src_n, sn) and int(sn.sum()) == n_cells
    ok &= good
    print(f"E state counts vs obs          : {S} states, total {int(sn.sum()):,} of {n_cells:,}  "
          f"{'PASS' if good else 'FAIL'}")

    # F ---------------------------------------------------------------------
    nbadf = int((~np.isfinite(emb)).sum())
    good = emb.shape == (n_cells, 2) and nbadf == 0
    ok &= good
    print(f"F embedding                    : shape {emb.shape}, {nbadf} non-finite  "
          f"{'PASS' if good else 'FAIL'}")

    # G ---------------------------------------------------------------------
    def nlines(path):
        n = 0
        with open(path, "r") as fh:
            for _ in fh:
                n += 1
        return n
    lu, lc, lg = nlines(p.umap), nlines(p.clusters), nlines(p.genes)
    good = lu == n_cells + 1 and lc == n_cells + 1 and lg == n_genes + 1
    ok &= good
    print(f"G legacy TSV rows              : umap {lu:,}, clusters {lc:,}, genes {lg:,} "
          f"(expect {n_cells + 1:,}/{n_cells + 1:,}/{n_genes + 1:,})  {'PASS' if good else 'FAIL'}")

    # H ---------------------------------------------------------------------
    with open(p.centroids, "r") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    good = hdr[1:] == states
    ok &= good
    print(f"H centroid columns             : {len(hdr) - 1} columns match the canonical order  "
          f"{'PASS' if good else 'FAIL'}")

    print("-" * 78)
    warns = sv.get("warnings", [])
    print(f"precompute warnings recorded   : {len(warns)}")
    for w in warns:
        print(f"  ! {w}")
    print(f"RESULT: {'ALL CHECKS PASS' if ok else 'FAILURES PRESENT'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
