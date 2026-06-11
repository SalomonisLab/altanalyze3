#!/usr/bin/env python3
"""Build the ISV web read-level RETRIEVAL index: a per-sample SQLite mapping
``gene -> [(molecule_id, cell_type, count)]`` so the molecule (read-level) view can fetch one gene's
molecules in ~ms instead of loading each sample's whole-matrix counts cache (~12 s for the largest
sample) on every render. Clustering is NOT precomputed -- it stays dynamic at visualization time; this
only makes molecule RETRIEVAL fast.

It rides on the same matrix pass precompute_viewer_index already does: build G (cells x cell_types),
``CT = G^T @ X`` (cell_types x molecules), then PARTITION CT by gene (using the existing
``*.gene_index.npz``) and write the per-gene blocks. The molecule structures themselves come from the
per-sample ``transcripts.db`` at retrieval time (already ~10 ms/gene), so they are not duplicated here.

Idempotent: a sample's db that already exists and is newer than its h5ad is skipped.

Usage (mirrors precompute_viewer_index):
    from altanalyze3.components.visualization.isv_web import precompute_reads_index as pri
    pri.build_reads_index(viewer_sample_dict, barcode_sample_dict, out_dir)
"""
import os
import sqlite3
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .. import isoform_structure_view as isv
from .. import precompute_viewer_index as pvi


def reads_db_path(out_dir, library):
    return os.path.join(out_dir, f"{library}.reads.db")


def _current(db_path, h5ad_path):
    return (os.path.exists(db_path) and os.path.exists(h5ad_path)
            and os.path.getmtime(db_path) >= os.path.getmtime(h5ad_path))


def _cell_types(bser):
    s = bser if isinstance(bser, pd.Series) else pd.Series(bser)
    return [v for v in pd.unique(s.dropna().astype(str)) if v and v.lower() != "nan"]


def _build_group_matrix(adata, bser, cell_types, mask_rev):
    """G (n_obs x n_cell_types) one-hot of each annotated barcode -> its cell type (orientation applied)."""
    obs = np.asarray(adata.obs_names, dtype=str)
    s = bser if isinstance(bser, pd.Series) else pd.Series(bser)
    s = s.dropna(); s.index = s.index.astype(str)
    ann = pd.Series({isv._reverse_complement_barcode(b): v for b, v in s.items()}) if mask_rev else s
    ct_index = {ct: j for j, ct in enumerate(cell_types)}
    rows, cols = [], []
    for i, bc in enumerate(obs):
        ct = ann.get(bc)
        if ct is not None and ct in ct_index:
            rows.append(i); cols.append(ct_index[ct])
    G = sp.csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)),
                      shape=(adata.n_obs, len(cell_types)))
    return G, ct_index


def build_one(h5ad_path, bser, db_path, log=print):
    import anndata as ad
    cell_types = _cell_types(bser)
    if not cell_types:
        log(f"[reads-index] {os.path.basename(h5ad_path)}: 0 cell types; skip"); return 0
    t0 = time.time()
    adata = ad.read_h5ad(h5ad_path)
    mask_rev = pvi._detect_orientation(adata, bser, cell_types)
    if mask_rev is None:
        log(f"[reads-index] {os.path.basename(h5ad_path)}: barcodes match neither orientation; skip")
        del adata; return 0
    G, ct_index = _build_group_matrix(adata, bser, cell_types, mask_rev)
    ct_labels = [None] * len(ct_index)
    for ct, j in ct_index.items():
        ct_labels[j] = ct
    X = adata.X
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    CT = (G.T @ X).tocsc()                                  # cell_types x molecules
    var_names = np.asarray(adata.var_names, dtype=str)
    del adata, X, G

    g2p, indptr, indices = isv._load_h5ad_gene_index(h5ad_path)

    tmp = db_path + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    cur = con.cursor()
    cur.execute("PRAGMA journal_mode=OFF"); cur.execute("PRAGMA synchronous=OFF")
    cur.execute("CREATE TABLE mol (gene TEXT, mol TEXT, ct TEXT, cnt REAL)")
    n_rows = 0
    for gene, pos in g2p.items():
        cols = np.asarray(indices[int(indptr[pos]):int(indptr[pos + 1])], dtype=np.int64)
        if cols.size == 0:
            continue
        sub = CT[:, cols].tocoo()
        if sub.nnz == 0:
            continue
        # molecule id = var token after the gene prefix ('ENSG..:12963' -> '12963')
        molids = [var_names[c].split(":", 1)[-1] for c in cols]
        batch = [(gene, molids[lc], ct_labels[r], float(v))
                 for r, lc, v in zip(sub.row, sub.col, sub.data) if v > 0]
        if batch:
            cur.executemany("INSERT INTO mol VALUES (?,?,?,?)", batch)
            n_rows += len(batch)
    cur.execute("CREATE INDEX ix_gene ON mol(gene)")
    con.commit(); con.close()
    os.replace(tmp, db_path)
    log(f"[reads-index] {os.path.basename(h5ad_path)}: {n_rows:,} rows, {len(ct_index)} cell types "
        f"({time.time() - t0:.1f}s)")
    return n_rows


def build_reads_index(viewer_sample_dict, barcode_sample_dict, out_dir, log=print):
    """Build/refresh the per-sample reads index for every sample. Returns {library: n_rows}."""
    os.makedirs(out_dir, exist_ok=True)
    built = {}
    t0 = time.time()
    for lib, entries in viewer_sample_dict.items():
        h5 = entries[0].get("matrix")
        bser = barcode_sample_dict.get(lib)
        if not h5 or not os.path.exists(str(h5)) or bser is None:
            log(f"[reads-index] {lib}: missing h5ad or barcodes; skip"); continue
        db = reads_db_path(out_dir, lib)
        if _current(db, h5):
            log(f"[reads-index] {lib}: up to date; skip"); built[lib] = -1; continue
        built[lib] = build_one(h5, bser, db, log=log)
    log(f"[reads-index] done ({time.time() - t0:.1f}s) for {len(built)} samples -> {out_dir}")
    return built


def _build_sample_dict(metadata, cell_annot, root, only_uid=None):
    """Construct sample_dict (library -> [{matrix, groups, library}]) + barcode_sample_dict EXACTLY as
    the workflow's own isoform_automate.precompute_isoform_viewer_index does, so matrix resolution is
    authoritative: the MOLECULE matrix is '<stem>-isoform.h5ad', with a '<root>/<uid>-isoform.h5ad'
    fallback for multi-BAM / hashed uids (the convention these datasets use)."""
    from ...long_read import isoform_automate as isoa
    from ...long_read import cli as _cli
    from ...long_read.isoform_matrix import import_barcode_clusters
    sd_meta = isoa.import_metadata(str(metadata), include_hashed_samples=True)
    barcode_sample_dict = import_barcode_clusters(_cli._discover_barcode_cluster_dirs(metadata, cell_annot))
    sample_dict = {}
    items = sd_meta.items() if only_uid is None else [(only_uid, sd_meta.get(only_uid, []))]
    for uid, libs in items:
        for s in libs:
            lib = s.get("library") or uid
            mol = isoa._isoform_molecule_h5ad(s.get("matrix"))
            if not os.path.exists(mol):
                alt = isoa._isoform_molecule_h5ad(os.path.join(root, str(uid)))
                if os.path.exists(alt):
                    mol = alt
            sample_dict[lib] = [{"matrix": str(mol), "groups": s.get("groups"), "library": lib}]
    return sample_dict, barcode_sample_dict


def main(argv=None):
    """Cluster entry point. Build the isv_web retrieval index (and optionally the engine viewer
    indexes) for the dataset's samples; --uid restricts to one metadata uid for parallel bsub jobs.

      python3 -m altanalyze3.components.visualization.isv_web.precompute_reads_index \
          --metadata M --gene_model EXON --gene_symbol SYM [--cell_annot CA] --root R \
          [--uid <uid>] [--viewer_index]
    """
    import argparse
    ap = argparse.ArgumentParser(description="Build isv_web per-gene molecule retrieval index.")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--gene_model", required=False)     # accepted for parity with the workflow scripts
    ap.add_argument("--gene_symbol", required=False)
    ap.add_argument("--cell_annot", default=None)
    ap.add_argument("--root", default=None)
    ap.add_argument("--uid", default=None, help="restrict to one metadata uid (for per-sample parallel jobs)")
    ap.add_argument("--viewer_index", action="store_true",
                    help="also build the engine viewer indexes (gene_index + transcripts.db + counts "
                         "caches) for these samples first (precompute_viewer_index)")
    a = ap.parse_args(argv)
    root = a.root or os.path.dirname(os.path.abspath(a.metadata))
    sample_dict, barcode_sample_dict = _build_sample_dict(a.metadata, a.cell_annot, root, only_uid=a.uid)
    if not sample_dict:
        print(f"[reads-index] no samples for uid={a.uid}; nothing to do"); return
    print(f"[reads-index] samples: {list(sample_dict.keys())}", flush=True)

    if a.viewer_index:
        from .. import precompute_viewer_index as pvi
        print("[reads-index] building engine viewer indexes first (precompute_viewer_index) ...", flush=True)
        pvi.precompute(sample_dict, barcode_sample_dict, groupby_rev=True,
                       only_samples=set(sample_dict.keys()))

    out_dir = os.path.join(root, "_isv_web_cache", "mol_index")
    build_reads_index(sample_dict, barcode_sample_dict, out_dir)

    # also build the molecule -> FINAL collapsed-isoform index (for the "group reads by final isoform"
    # molecule view); no-op if FINAL_structure_to_exemplar.tsv / tier1 mol2struct files are absent.
    try:
        from . import precompute_final_isoform_index as pfi
        final_tsv, lib2m2s = pfi.discover(root)
        pfi.build_mol2final_index(lib2m2s, final_tsv, out_dir)
    except Exception as e:
        print(f"[reads-index] mol2final index skipped: {e}")


if __name__ == "__main__":
    main()
