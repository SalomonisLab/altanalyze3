"""Build a scalable_viewer bundle from an h5ad.

The server never opens an h5ad. This script does all the scanning once and writes the
bundle described in bundle.py.

  PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
  /opt/homebrew/opt/python@3.11/bin/python3.11 \
    -m altanalyze3.components.visualization.scalable_viewer.precompute \
    --h5ad   /path/to/data.h5ad \
    --out    /path/to/bundles/MyDataset \
    --prefix My-Dataset \
    --label  "My dataset" \
    --cluster-key cell_state \
    --layer  lognorm \
    --markers /path/to/markers.tsv \
    --order   /path/to/canonical_order.tsv \
    --deg     /path/to/differential

Method notes that the README repeats, so a reader never has to guess:

* Expression store. layers[--layer] is a cell-major CSR in the h5ad. This script
  transposes it to a gene-major (CSC) triple of .npy files with a two-pass counting
  sort. Values are copied at their source dtype (float32), so the store is bit-exact.
  Pass --expr-dtype float16 to halve the size and accept the rounding.
* Per-state statistics. mean is the UNWEIGHTED mean over the cells (or metacells) of a
  state, not a cell-count-weighted mean. frac is the fraction of those cells with a
  value above zero. The denominator per state is written to <prefix>_stats_n.npy.
* Highly variable genes. dispersion = variance / mean computed on the layer values.
  Genes are put in 20 equal-count bins of mean, the dispersion is z-scored inside each
  bin, and the top --n-hvg genes by that z-score are kept. This is a described method;
  it is NOT a call into scanpy.pp.highly_variable_genes.
* Embedding. PCA (randomised SVD) on the z-scored HVG matrix, then UMAP. The source
  h5ad may carry no obsm embedding; --embedding-from obsm:<key> reuses one when it does.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

from . import bundle as B

_LOG_T0 = time.time()


def log(msg: str) -> None:
    print(f"[precompute {time.time() - _LOG_T0:8.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- h5ad readers


def _decode(arr) -> List[str]:
    return [x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]


def read_obs_column(obs: h5py.Group, name: str):
    """Return ('categorical', codes int32, categories list) or ('numeric', values float64, None)
    or ('string', codes int32, categories list). Raises on an unreadable column."""
    node = obs[name]
    if isinstance(node, h5py.Group):                       # anndata categorical
        cats = _decode(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int32)
        return "categorical", codes, cats
    vals = node[:]
    if vals.dtype.kind in "OSU":                           # plain string column
        strs = _decode(vals)
        cats = sorted(set(strs))
        lut = {c: i for i, c in enumerate(cats)}
        return "string", np.asarray([lut[s] for s in strs], dtype=np.int32), cats
    if vals.dtype.kind == "b":
        return "categorical", vals.astype(np.int32), ["False", "True"]
    return "numeric", np.asarray(vals, dtype=np.float64), None


def read_uns_json(uns: h5py.Group, key: str) -> Optional[dict]:
    if key not in uns:
        return None
    raw = uns[key][()]
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    try:
        return json.loads(raw)
    except (ValueError, TypeError):
        return None


def read_uns_strlist(uns: h5py.Group, key: str) -> Optional[List[str]]:
    if key not in uns:
        return None
    node = uns[key]
    if isinstance(node, h5py.Group):
        return None
    return _decode(node[:])


# ---------------------------------------------------- CSR -> CSC + statistics


def build_expression_store(
    grp: h5py.Group,
    n_cells: int,
    n_genes: int,
    state_code: np.ndarray,
    n_states: int,
    paths: B.BundlePaths,
    expr_dtype: str,
    row_block: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Transpose the cell-major CSR layer into a gene-major store and accumulate stats.

    Returns (sum_gs, cnt_gs, gene_sum, gene_sumsq, nnz).
    sum_gs / cnt_gs are (n_genes, n_states) float64 / int64.
    """
    indptr_ds, indices_ds, data_ds = grp["indptr"], grp["indices"], grp["data"]
    nnz = int(indptr_ds[-1])
    log(f"layer nnz = {nnz:,}  density = {nnz / (n_cells * n_genes):.4f}")

    # ---- pass 1: nnz per gene -------------------------------------------------
    counts = np.zeros(n_genes, dtype=np.int64)
    step = 64_000_000
    for s in range(0, nnz, step):
        e = min(s + step, nnz)
        counts += np.bincount(indices_ds[s:e].astype(np.int64), minlength=n_genes)
    if counts.sum() != nnz:
        raise RuntimeError(f"pass 1 count mismatch: {counts.sum()} != {nnz}")
    log(f"pass 1 done. genes with zero nnz: {int((counts == 0).sum()):,} of {n_genes:,}")

    out_indptr = np.zeros(n_genes + 1, dtype=np.int64)
    np.cumsum(counts, out=out_indptr[1:])
    if out_indptr[-1] != nnz:
        raise RuntimeError("indptr tail does not equal nnz")

    # ---- allocate the gene-major store ---------------------------------------
    np.save(paths.expr_indptr, out_indptr)
    dt = np.float16 if expr_dtype == "float16" else np.float32
    out_idx = np.lib.format.open_memmap(paths.expr_indices, mode="w+", dtype=np.uint32, shape=(nnz,))
    out_val = np.lib.format.open_memmap(paths.expr_data, mode="w+", dtype=dt, shape=(nnz,))
    log(f"allocated {paths.expr_indices} ({nnz * 4 / 2**30:.2f} GiB) and "
        f"{paths.expr_data} ({nnz * np.dtype(dt).itemsize / 2**30:.2f} GiB)")

    # ---- pass 2: scatter + statistics ----------------------------------------
    write_ptr = out_indptr[:-1].copy()
    gs_len = n_genes * n_states
    sum_gs = np.zeros(gs_len, dtype=np.float64)
    cnt_gs = np.zeros(gs_len, dtype=np.int64)
    gene_sumsq = np.zeros(n_genes, dtype=np.float64)
    state_code64 = state_code.astype(np.int64)
    full_indptr = indptr_ds[:].astype(np.int64)

    n_blocks = (n_cells + row_block - 1) // row_block
    for bi, r0 in enumerate(range(0, n_cells, row_block)):
        r1 = min(r0 + row_block, n_cells)
        s, e = int(full_indptr[r0]), int(full_indptr[r1])
        if e == s:
            continue
        gi = indices_ds[s:e].astype(np.int64)
        dv = np.asarray(data_ds[s:e], dtype=np.float32)
        rowlens = np.diff(full_indptr[r0:r1 + 1])
        ci = np.repeat(np.arange(r0, r1, dtype=np.int64), rowlens)

        # statistics on the unsorted chunk
        key = gi * n_states + state_code64[ci]
        sum_gs += np.bincount(key, weights=dv.astype(np.float64), minlength=gs_len)
        cnt_gs += np.bincount(key, minlength=gs_len)
        gene_sumsq += np.bincount(gi, weights=(dv.astype(np.float64) ** 2), minlength=n_genes)
        del key

        # stable sort by gene keeps cell indices ascending inside every gene, and
        # blocks are visited in ascending cell order, so the CSC comes out sorted.
        order = np.argsort(gi, kind="stable")
        gs_sorted = gi[order]
        uniq, first_idx, cnts = np.unique(gs_sorted, return_index=True, return_counts=True)
        dest = (np.repeat(write_ptr[uniq], cnts)
                + (np.arange(gs_sorted.size, dtype=np.int64) - np.repeat(first_idx, cnts)))
        out_idx[dest] = ci[order].astype(np.uint32)
        out_val[dest] = dv[order].astype(dt)
        write_ptr[uniq] += cnts
        del gi, dv, ci, order, gs_sorted, uniq, first_idx, cnts, dest
        log(f"pass 2 block {bi + 1}/{n_blocks} cells {r0:,}-{r1:,}")

    if not np.array_equal(write_ptr, out_indptr[1:]):
        bad = int((write_ptr != out_indptr[1:]).sum())
        raise RuntimeError(f"pass 2 did not fill the store: {bad} genes have a wrong write pointer")
    if cnt_gs.sum() != nnz:
        raise RuntimeError(f"statistics count {cnt_gs.sum()} != nnz {nnz}")
    out_idx.flush(); out_val.flush()
    del out_idx, out_val
    log("pass 2 done, store flushed, invariants held")

    gene_sum = sum_gs.reshape(n_genes, n_states).sum(axis=1)
    return (sum_gs.reshape(n_genes, n_states), cnt_gs.reshape(n_genes, n_states),
            gene_sum, gene_sumsq, nnz)


# ------------------------------------------------------------------ HVG + PCA


def select_hvg(gene_mean: np.ndarray, gene_var: np.ndarray, n_top: int, n_bins: int = 20) -> np.ndarray:
    """Binned-dispersion HVG on the layer values. See the module docstring."""
    ok = (gene_mean > 0) & (gene_var > 0)
    disp = np.full(gene_mean.shape, np.nan, dtype=np.float64)
    disp[ok] = np.log(gene_var[ok] / gene_mean[ok])
    logmean = np.full(gene_mean.shape, np.nan, dtype=np.float64)
    logmean[ok] = np.log(gene_mean[ok])

    score = np.full(gene_mean.shape, -np.inf, dtype=np.float64)
    idx_ok = np.flatnonzero(ok)
    if idx_ok.size == 0:
        raise RuntimeError("no gene has a positive mean and variance; cannot pick HVGs")
    qs = np.quantile(logmean[idx_ok], np.linspace(0, 1, n_bins + 1))
    qs[0] -= 1e-9
    binid = np.digitize(logmean[idx_ok], qs[1:-1], right=True)
    for b in range(n_bins):
        sel = idx_ok[binid == b]
        if sel.size < 2:
            continue
        d = disp[sel]
        sd = d.std()
        score[sel] = (d - d.mean()) / sd if sd > 0 else 0.0
    n_top = min(n_top, int(np.isfinite(score).sum()))
    hvg = np.argsort(-score, kind="stable")[:n_top]
    return np.sort(hvg)


def dense_from_csc(paths: B.BundlePaths, gene_idx: np.ndarray, n_cells: int) -> np.ndarray:
    """Materialise a (n_cells, len(gene_idx)) dense float32 block from the gene-major store."""
    indptr = np.load(paths.expr_indptr, mmap_mode="r")
    ind = np.load(paths.expr_indices, mmap_mode="r")
    val = np.load(paths.expr_data, mmap_mode="r")
    out = np.zeros((n_cells, gene_idx.size), dtype=np.float32)
    for j, g in enumerate(gene_idx):
        s, e = int(indptr[g]), int(indptr[g + 1])
        if e > s:
            out[np.asarray(ind[s:e], dtype=np.int64), j] = np.asarray(val[s:e], dtype=np.float32)
    return out


def compute_embedding(X: np.ndarray, n_pcs: int, n_neighbors: int, min_dist: float,
                      seed: int) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """z-score -> randomised-SVD PCA -> UMAP. Returns (umap (N,2), pcs (N,n_pcs), warnings)."""
    from sklearn.utils.extmath import randomized_svd

    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    X -= mu
    X /= sd
    np.clip(X, -10, 10, out=X)
    log(f"z-scored HVG matrix {X.shape}, running randomized_svd for {n_pcs} components")
    U, S, _Vt = randomized_svd(X, n_components=n_pcs, random_state=seed)
    pcs = (U * S).astype(np.float32)
    log(f"PCA done. top-5 singular values: {np.round(S[:5], 2).tolist()}")

    # umap-learn 0.5.7 imports the optional parametric_umap, which imports tensorflow.
    # tensorflow 2.x on numpy 1.23 raises AttributeError (np.dtypes), which umap's own
    # `except ImportError` does not catch. Forcing an ImportError makes umap take its
    # documented "Tensorflow not installed" path. No installed file is modified.
    sys.modules.setdefault("tensorflow", None)
    import umap  # noqa: E402

    caught: List[str] = []
    log(f"running UMAP on {pcs.shape} (n_neighbors={n_neighbors}, min_dist={min_dist}, seed={seed})")
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        emb = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=2,
                        random_state=seed, verbose=True).fit_transform(pcs)
        for w in wlist:
            caught.append(f"{w.category.__name__}: {w.message}")
    emb = np.asarray(emb, dtype=np.float32)
    if emb.shape != (X.shape[0], 2):
        raise RuntimeError(f"UMAP returned {emb.shape}, expected {(X.shape[0], 2)}")
    n_bad = int((~np.isfinite(emb)).sum())
    if n_bad:
        raise RuntimeError(f"UMAP produced {n_bad} non-finite coordinates")
    return emb, pcs, caught


# ------------------------------------------------------------- side-car tables


def ingest_markers(src: Optional[str], dst: str) -> Dict:
    """Normalise a marker table to gene / cluster / fold / p columns. Returns a summary."""
    if not src:
        return {"available": False, "n_rows": 0, "path": None, "source": None}
    with open(src, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    low = [h.strip().lower() for h in header]

    def find(*names):
        for n in names:
            if n in low:
                return low.index(n)
        return None

    i_gene = find("gene", "uid", "symbol")
    i_clu = find("cluster", "population", "cell_state", "celltype")
    i_fold = find("fold", "log2fc", "logfc", "fold change")
    i_p = find("fdr p-value", "fdr", "adjp", "padj", "fdr_p")
    if i_gene is None or i_clu is None:
        raise RuntimeError(f"marker table {src} has no gene/cluster column. header={header}")
    out = ["gene\tcluster\tfold\tp\t" + "\t".join(header)]
    for r in rows:
        g = r[i_gene] if i_gene < len(r) else ""
        c = r[i_clu] if i_clu < len(r) else ""
        f = r[i_fold] if (i_fold is not None and i_fold < len(r)) else ""
        p = r[i_p] if (i_p is not None and i_p < len(r)) else ""
        out.append("\t".join([g, c, f, p] + r))
    with open(dst, "w") as fh:
        fh.write("\n".join(out) + "\n")
    clusters = sorted({ln.split("\t")[1] for ln in out[1:]})
    return {"available": True, "n_rows": len(rows), "path": dst, "source": os.path.abspath(src),
            "n_clusters": len(clusters), "clusters": clusters,
            "columns": ["gene", "cluster", "fold", "p"] + header}


def ingest_deg(deg_root: Optional[str], deg_dir: str, manifest_path: str,
               modality_roots: Optional[Dict[str, str]] = None) -> Dict:
    """Copy every DEG_detailed_*.tsv / DEG_pooled_overall_*.tsv found under each root.

    Nothing is recomputed. Each file is copied verbatim and indexed, so the browser
    shows the numbers the differential workflow produced.

    `deg_root` holds the RNA differentials. `modality_roots` maps a modality id to the
    root of that modality's own differential runs, so one bundle can serve an ADT, lipid
    or GRN contrast beside the RNA one. A non-RNA file is copied under a modality
    subdirectory, because two modalities of the same contrast share a file name.
    """
    manifest = {"comparisons": []}
    roots: List[Tuple[str, Optional[str]]] = [("rna", deg_root)]
    for modality_id, root in sorted((modality_roots or {}).items()):
        roots.append((modality_id, root))

    for modality_id, root in roots:
        if not root or not os.path.isdir(root):
            continue
        target_dir = deg_dir if modality_id == "rna" else os.path.join(deg_dir, modality_id)
        os.makedirs(target_dir, exist_ok=True)
        for dirpath, _d, filenames in os.walk(root):
            for fn in sorted(filenames):
                if not fn.endswith(".tsv"):
                    continue
                kind = None
                if fn.startswith("DEG_detailed_"):
                    kind = "per_cell_state"
                    comp = fn[len("DEG_detailed_"):-4]
                elif fn.startswith("DEG_pooled_overall_"):
                    kind = "pooled_overall"
                    comp = fn[len("DEG_pooled_overall_"):-4]
                else:
                    continue
                src = os.path.join(dirpath, fn)

                # THE CONTRAST COMES FROM THE PATH, NOT THE FILE NAME.
                #
                # cellHarmony names every comparison directory after the two arm labels, and
                # a covariate file that labels its arms CASE and CONTROL therefore yields
                # CASE_vs_CONTROL for every contrast in the project. Deriving `comp` from the
                # file name alone gave every contrast the same id AND the same destination
                # path, so each copy overwrote the previous one. A CellRef2 bundle built on
                # 2026-09-03 recorded 505 tables, kept 7, and lost 100 of 101 contrasts.
                #
                # The layout is <root>/<contrast>/<comparison>/DEGs/<file>, so the first
                # component of the path relative to the root names the contrast. A flat root
                # with no contrast directory keeps the previous behaviour unchanged.
                rel = os.path.relpath(dirpath, root)
                parts = [p for p in rel.split(os.sep) if p not in ("", ".")]
                contrast = parts[0] if parts else ""

                if contrast:
                    out_dir = os.path.join(target_dir, contrast)
                    os.makedirs(out_dir, exist_ok=True)
                else:
                    out_dir = target_dir
                dst = os.path.join(out_dir, fn)

                with open(src, "r") as fi:
                    text = fi.read()
                with open(dst, "w") as fo:
                    fo.write(text)
                lines = [ln for ln in text.split("\n") if ln.strip()]
                header = lines[0].split("\t") if lines else []
                label = f"{contrast}::{comp}" if contrast else comp
                entry_id = (f"{label}::{kind}" if modality_id == "rna"
                            else f"{modality_id}::{label}::{kind}")
                manifest["comparisons"].append({
                    "id": entry_id, "comparison": comp, "kind": kind,
                    "contrast": contrast,
                    "modality": modality_id,
                    "file": os.path.relpath(dst, deg_dir), "path": os.path.abspath(dst),
                    "source": os.path.abspath(src), "columns": header,
                    "n_rows": max(len(lines) - 1, 0),
                })
    manifest["comparisons"].sort(key=lambda c: (c.get("modality", "rna") != "rna",
                                                c.get("modality", "rna"),
                                                c.get("contrast", ""),
                                                c["comparison"], c["kind"]))
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def ingest_ccc(src: Optional[str], dst: str) -> Dict:
    if not src or not os.path.isfile(src):
        return {"available": False, "n_rows": 0, "path": None, "source": None}
    with open(src, "r") as fi:
        text = fi.read()
    with open(dst, "w") as fo:
        fo.write(text)
    lines = [ln for ln in text.split("\n") if ln.strip()]
    return {"available": True, "n_rows": max(len(lines) - 1, 0), "path": os.path.abspath(dst),
            "source": os.path.abspath(src),
            "columns": lines[0].split("\t") if lines else []}


# ------------------------------------------------------------------- modalities


def _read_feature_matrix(path: str):
    """Read a modality prediction table as (row labels, feature names, float32 matrix).

    rna2adt, rna2lipid and rna2grn all write a delimited table through
    `predictions.to_csv` (rna2adt/cli.py:25), so a `.h5ad` extension on one of those
    outputs names a CSV, not HDF5. The format is decided by the file's own first bytes,
    never by its extension.
    """
    import pandas as pd

    with open(path, "rb") as fh:
        head = fh.read(8)
    if head[:8] == b"\x89HDF\r\n\x1a\n":
        import anndata as ad

        adata = ad.read_h5ad(path)
        matrix = adata.X
        matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        return (np.asarray(adata.obs_names, dtype=str),
                np.asarray(adata.var_names, dtype=str),
                np.asarray(matrix, dtype=np.float32))

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    frame = pd.read_csv(path, sep=sep, index_col=0)
    return (frame.index.astype(str).to_numpy(),
            frame.columns.astype(str).to_numpy(),
            np.ascontiguousarray(frame.to_numpy(dtype=np.float32)))


def _read_display_names(path: str) -> Dict[str, str]:
    """`feature<TAB>display` from a TSV, ignoring any further columns.

    A header line naming the first column `feature` is skipped; a file without one is
    read from its first row, so a hand-made two-column list works.
    """
    out: Dict[str, str] = {}
    with open(path, "r") as fh:
        first = True
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if first:
                first = False
                if parts[0].strip().lower() in ("feature", "id", "gene", "symbol"):
                    continue
            key, shown = parts[0].strip(), parts[1].strip()
            if key and shown:
                out[key] = shown
    if not out:
        raise SystemExit(f"--modality-display-names read no pair from {path}")
    return out


def ingest_modality(
    modality_id: str,
    source: str,
    *,
    paths: B.BundlePaths,
    barcodes: List[str],
    states: List[str],
    state_code: np.ndarray,
    state_n: np.ndarray,
    label: str = "",
    feature_label: str = "feature",
    expr_dtype: str = "float32",
    display_names: Optional[Dict[str, str]] = None,
) -> Dict:
    """Write one modality's sidecars beside the RNA store and report what aligned.

    The table's row labels decide the store kind. Rows that are cell barcodes give a
    `per_cell` store, the same feature-major layout the RNA store uses. Rows that are
    cell-state names give a `per_state` store: the modality was predicted at cell-state
    granularity, so every cell of a state carries that state's value and only the
    feature x cell-state matrix is stored. Nothing is interpolated in either case.
    """
    mpaths = paths.modality(modality_id)
    rows, features, matrix = _read_feature_matrix(source)
    n_features = len(features)
    n_states = len(states)
    n_cells = len(barcodes)
    log(f"[{modality_id}] {source}: {matrix.shape[0]:,} rows x {n_features:,} features")

    barcode_pos = {str(b): i for i, b in enumerate(barcodes)}
    state_pos = {str(s): i for i, s in enumerate(states)}
    hit_cells = sum(1 for r in rows if str(r) in barcode_pos)
    hit_states = sum(1 for r in rows if str(r) in state_pos)

    if hit_cells >= hit_states and hit_cells > 0:
        kind = "per_cell"
    elif hit_states > 0:
        kind = "per_state"
    else:
        raise RuntimeError(
            f"modality '{modality_id}': no row label of {source} matches a bundle barcode "
            f"or a cell state. First rows: {list(rows[:3])}. "
            f"First bundle barcodes: {barcodes[:2]}. First states: {states[:2]}")

    if kind == "per_cell":
        matched = np.full(n_cells, -1, dtype=np.int64)
        for source_row, label_value in enumerate(rows):
            target = barcode_pos.get(str(label_value))
            if target is not None and matched[target] < 0:
                matched[target] = source_row
        n_matched = int((matched >= 0).sum())
        retention = n_matched / n_cells if n_cells else 0.0
        log(f"[{modality_id}] per-cell store: {n_matched:,} of {n_cells:,} bundle cells "
            f"matched a row ({retention:.4f})")
        if retention < 0.90:
            log(f"[WARN] [{modality_id}] retention below 90%: {n_matched:,} of {n_cells:,}")
        dense = np.zeros((n_cells, n_features), dtype=np.float32)
        present = matched >= 0
        dense[present] = matrix[matched[present]]
        # A cell with no row in the table keeps zeros; the count above is the denominator.
    else:
        matched = np.full(n_states, -1, dtype=np.int64)
        for source_row, label_value in enumerate(rows):
            target = state_pos.get(str(label_value))
            if target is not None and matched[target] < 0:
                matched[target] = source_row
        n_matched = int((matched >= 0).sum())
        retention = n_matched / n_states if n_states else 0.0
        log(f"[{modality_id}] per-state store: {n_matched:,} of {n_states:,} cell states "
            f"matched a row ({retention:.4f})")
        unmatched = [states[i] for i in range(n_states) if matched[i] < 0]
        if unmatched:
            log(f"[WARN] [{modality_id}] {len(unmatched)} cell states carry no prediction "
                f"and are served as zero: {unmatched}")
        state_matrix = np.zeros((n_states, n_features), dtype=np.float32)
        present = matched >= 0
        state_matrix[present] = matrix[matched[present]]

    # ---- statistics, in the same definition the RNA store uses ----------------
    counts = np.asarray(state_n, dtype=np.int64)
    denominator = np.maximum(counts, 1)
    if kind == "per_cell":
        codes = np.asarray(state_code, dtype=np.int64)
        sums = np.zeros((n_states, n_features), dtype=np.float64)
        nonzero = np.zeros((n_states, n_features), dtype=np.int64)
        np.add.at(sums, codes, dense.astype(np.float64))
        np.add.at(nonzero, codes, (dense != 0.0))
        mean_fs = (sums / denominator[:, None]).T.astype(np.float32)
        frac_fs = (nonzero / denominator[:, None]).T.astype(np.float32)
    else:
        mean_fs = state_matrix.T.astype(np.float32)
        frac_fs = (state_matrix != 0.0).T.astype(np.float32)

    np.save(mpaths.stats_mean, mean_fs)
    np.save(mpaths.stats_frac, frac_fs)
    # A FOURTH COLUMN CARRIES THE NAME A READER SEES, AND THE KEY NEVER MOVES.
    #
    # Nathan, 2026-09-08: scALABLE should show "CE(18:2)" as "18:2 Cholesterol ester",
    # and the abbreviation must stay searchable. So `symbol` keeps the key every
    # differential table, saved URL and Chat question already uses, and `display` holds
    # the reader's name. data_api indexes BOTH, so either resolves.
    #
    # The column is written only when a mapping is given, so a bundle built without one
    # is byte-identical to before.
    n_display = 0
    with open(mpaths.genes, "w") as fh:
        header = "index\tgene_id\tsymbol"
        fh.write(header + ("\tdisplay\n" if display_names else "\n"))
        for i, name in enumerate(features):
            if display_names:
                shown = str(display_names.get(name) or name)
                if shown != name:
                    n_display += 1
                fh.write(f"{i}\t{name}\t{name}\t{shown}\n")
            else:
                fh.write(f"{i}\t{name}\t{name}\n")
    if display_names:
        missing = [f for f in features if f not in display_names]
        log(f"[{modality_id}] display names: {n_display:,} of {n_features:,} features "
            f"renamed ({n_display / max(1, n_features):.4f})")
        if missing:
            log(f"[WARN] [{modality_id}] {len(missing)} of {n_features:,} features carry "
                f"no display name and keep their key, first: {missing[:5]}")

    nnz = 0
    if kind == "per_cell":
        dtype = np.float16 if expr_dtype == "float16" else np.float32
        columns = []
        indptr = np.zeros(n_features + 1, dtype=np.int64)
        for feature_row in range(n_features):
            column = dense[:, feature_row]
            keep = np.nonzero(column)[0]
            columns.append((keep.astype(np.uint32), column[keep].astype(dtype)))
            indptr[feature_row + 1] = indptr[feature_row] + keep.size
        nnz = int(indptr[-1])
        out_idx = np.lib.format.open_memmap(mpaths.expr_indices, mode="w+",
                                            dtype=np.uint32, shape=(max(nnz, 1),))
        out_val = np.lib.format.open_memmap(mpaths.expr_data, mode="w+",
                                            dtype=dtype, shape=(max(nnz, 1),))
        for feature_row, (keep, values) in enumerate(columns):
            start = int(indptr[feature_row])
            out_idx[start:start + keep.size] = keep
            out_val[start:start + keep.size] = values
        out_idx.flush(); out_val.flush()
        del out_idx, out_val
        np.save(mpaths.expr_indptr, indptr)
        written = int(np.load(mpaths.expr_indptr)[-1])
        if written != nnz:
            raise RuntimeError(f"modality '{modality_id}' store tail {written} != nnz {nnz}")
        log(f"[{modality_id}] feature-major store: {nnz:,} non-zero values "
            f"of {n_cells * n_features:,} ({nnz / max(n_cells * n_features, 1):.4f} dense)")

    info = {
        "id": modality_id,
        "label": label or modality_id.upper(),
        "feature_label": feature_label,
        "kind": kind,
        "source": os.path.abspath(source),
        "source_mtime": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                      time.gmtime(os.path.getmtime(source))),
        "n_features": int(n_features),
        "n_source_rows": int(matrix.shape[0]),
        "n_matched": int(n_matched),
        "n_expected": int(n_cells if kind == "per_cell" else n_states),
        "retention": float(retention),
        "nnz": nnz,
        "expr_dtype": expr_dtype if kind == "per_cell" else None,
        "stats_method": ("mean = unweighted mean over the cells of a state; frac = fraction "
                         "of those cells above zero"
                         if kind == "per_cell" else
                         "the predicted cell-state value; every cell of the state carries it"),
    }
    return info


def _parse_id_value(pairs: Optional[List[str]], flag: str) -> Dict[str, str]:
    """`--flag id=value` pairs into a dict, with the whole pair named on an error."""
    out: Dict[str, str] = {}
    for raw in pairs or []:
        if "=" not in raw:
            raise SystemExit(f"{flag} expects id=value, got {raw!r}")
        key, _, value = raw.partition("=")
        key = key.strip()
        if not key:
            raise SystemExit(f"{flag} expects id=value, got {raw!r}")
        out[key] = value.strip()
    return out


def read_canonical_order(path: Optional[str]) -> Tuple[Optional[List[str]], Dict[str, str]]:
    """Read an order TSV with columns order / cell_state / color."""
    if not path or not os.path.isfile(path):
        return None, {}
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        low = [h.strip().lower() for h in header]
        i_state = low.index("cell_state") if "cell_state" in low else 1
        i_color = low.index("color") if "color" in low else None
        i_ord = low.index("order") if "order" in low else 0
        rows = []
        for ln in fh:
            if not ln.strip():
                continue
            r = ln.rstrip("\n").split("\t")
            try:
                o = int(r[i_ord])
            except (ValueError, IndexError):
                o = len(rows)
            rows.append((o, r[i_state], r[i_color] if (i_color is not None and i_color < len(r)) else None))
    rows.sort(key=lambda x: x[0])
    order = [r[1] for r in rows]
    colors = {r[1]: r[2] for r in rows if r[2]}
    return order, colors


# ------------------------------------------------------------------------ main


def _add_modalities_to_bundle(a, sources: Dict[str, str], labels: Dict[str, str],
                              feature_labels: Dict[str, str],
                              deg_roots: Optional[Dict[str, str]] = None,
                              display_names: Optional[Dict[str, Dict[str, str]]] = None) -> int:
    """Write modality sidecars into a bundle that already exists.

    Reads the bundle's own metadata for the cell barcodes, the cell states and the
    per-state cell counts, so the modality store is aligned to the cells the viewer
    serves and nothing about the RNA store is recomputed or rewritten.
    """
    paths = B.BundlePaths(a.out, a.prefix)
    if not os.path.isfile(paths.metadata):
        raise SystemExit(f"--modalities-only needs a built bundle; {paths.metadata} is absent")
    meta = B.read_metadata(paths.metadata)
    block = B.viewer_block(meta)
    if block is None:
        raise SystemExit(f"{paths.metadata} has no 'scalable_viewer' block")

    states: List[str] = list(block["states"])
    state_n = np.asarray(block["state_n"], dtype=np.int64)
    barcodes: List[str] = []
    codes: List[int] = []
    position = {s: i for i, s in enumerate(states)}
    with open(paths.clusters, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        state_column = 1 if len(header) > 1 else 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            barcodes.append(parts[0])
            codes.append(position.get(parts[state_column], -1))
    state_code = np.asarray(codes, dtype=np.int64)
    unassigned = int((state_code < 0).sum())
    if unassigned:
        raise SystemExit(f"{paths.clusters} holds {unassigned} rows whose cell state is not in "
                         f"the bundle's state list")
    if len(barcodes) != int(meta["n_cells"]):
        raise SystemExit(f"{paths.clusters} holds {len(barcodes)} barcodes but the bundle "
                         f"declares {meta['n_cells']} cells")
    counted = np.bincount(state_code, minlength=len(states))
    if not np.array_equal(counted, state_n):
        raise SystemExit("cells per state read from the clusters table disagree with the "
                         "bundle metadata; the bundle is inconsistent")
    log(f"bundle {paths.bundle_dir}: {len(barcodes):,} cells, {len(states)} cell states")

    existing = dict(block.get("modalities") or {})
    for modality_id, source in sources.items():
        info = ingest_modality(
            modality_id, source, paths=paths, barcodes=barcodes, states=states,
            state_code=state_code, state_n=state_n,
            label=labels.get(modality_id, ""),
            feature_label=feature_labels.get(modality_id, "feature"),
            expr_dtype=a.expr_dtype,
            display_names=(display_names or {}).get(modality_id),
        )
        existing[modality_id] = info
    block["modalities"] = existing

    # Modality differentials are MERGED into the DEG manifest. Only the named modality
    # roots are read here, so the RNA comparisons the bundle already indexes are kept.
    if deg_roots:
        added = ingest_deg(None, paths.deg_dir, paths.deg_manifest, modality_roots=deg_roots)
        manifest = dict(block.get("deg") or {"comparisons": []})
        kept = [c for c in manifest.get("comparisons", [])
                if str(c.get("modality") or "rna") not in deg_roots]
        manifest["comparisons"] = kept + added["comparisons"]
        manifest["comparisons"].sort(key=lambda c: (str(c.get("modality") or "rna") != "rna",
                                                    str(c.get("modality") or "rna"),
                                                    c["comparison"], c["kind"]))
        block["deg"] = manifest
        with open(paths.deg_manifest, "w") as fh:
            json.dump(manifest, fh, indent=2)
        for modality_id in sorted(deg_roots):
            n = sum(1 for c in added["comparisons"] if c.get("modality") == modality_id)
            log(f"[{modality_id}] differential tables indexed: {n}")
        log(f"DEG manifest: {len(manifest['comparisons'])} tables "
            f"({len(kept)} kept, {len(added['comparisons'])} added)")

    block["modalities_updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(paths.metadata, "w") as fh:
        json.dump(meta, fh, indent=2)
    log(f"updated {paths.metadata}: modalities = {sorted(existing)}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build a scalable_viewer bundle from an h5ad.")
    ap.add_argument("--h5ad", default=None,
                    help="source h5ad. Required unless --modalities-only, which reads only "
                         "the built bundle")
    ap.add_argument("--out", required=True, help="bundle directory (created if absent)")
    ap.add_argument("--prefix", required=True, help="file-name prefix inside the bundle")
    ap.add_argument("--label", default=None, help="human label for the catalog")
    ap.add_argument("--dataset-id", default=None, help="catalog id (default: --prefix)")
    ap.add_argument("--study-id", default=None,
                    help="LungMAP study id this bundle belongs to, for example "
                         "lmdata:LMEX0000009416. The viewer's Study tab reads it. Without "
                         "it the Study tab shows no record rather than another study's")
    ap.add_argument("--cluster-key", default="cell_state", help="obs column holding the cell state")
    ap.add_argument("--layer", default="lognorm", help="layer to serve ('X' for adata.X)")
    ap.add_argument("--markers", default=None)
    ap.add_argument("--order", default=None, help="canonical order TSV (order/cell_state/color)")
    ap.add_argument("--deg", default=None, help="directory tree holding the RNA DEG_*.tsv tables")
    ap.add_argument("--deg-modality", action="append", default=None, metavar="ID=DIR",
                    help="directory tree holding one modality's own DEG_*.tsv tables, for "
                         "example `--deg-modality adt=/path/to/adt_differential`. Repeatable.")
    ap.add_argument("--ccc", default=None, help="cell-cell communication TSV")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--min-dist", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--row-block", type=int, default=8192)
    ap.add_argument("--expr-dtype", choices=["float32", "float16"], default="float32")
    ap.add_argument("--embedding-from", default=None,
                    help="obsm:<key> to reuse an existing embedding instead of computing one")
    ap.add_argument("--max-centroid-genes", type=int, default=4000,
                    help="rows written to <prefix>.txt (marker genes first)")
    ap.add_argument("--skip-expr", action="store_true",
                    help="reuse an existing expression store and statistics (rebuild only the rest)")
    ap.add_argument("--modality", action="append", default=None, metavar="ID=PATH",
                    help="add a modality store: `--modality adt=/path/to/predictions.csv`. "
                         "Repeatable. The table's row labels decide the store: cell barcodes "
                         "give a per-cell store, cell-state names give a per-state store "
                         "broadcast to the cells of the state. rna2adt, rna2lipid and rna2grn "
                         "write CSV even when the file is named .h5ad; the format is read from "
                         "the file's own bytes.")
    ap.add_argument("--modality-label", action="append", default=None, metavar="ID=LABEL",
                    help="human label for a modality, for example `--modality-label adt='ADT (imputed)'`")
    ap.add_argument("--modality-feature-label", action="append", default=None, metavar="ID=NOUN",
                    help="what one feature of a modality is called, for example "
                         "`--modality-feature-label lipids=lipid`")
    ap.add_argument("--modality-display-names", action="append", default=None,
                    metavar="ID=TSV",
                    help="a two-column TSV, `feature<TAB>display`, naming what a reader "
                         "sees for each feature of a modality, for example "
                         "`--modality-display-names lipid=/path/lipid_display_names.tsv`. "
                         "The feature name stays the key and stays searchable; only the "
                         "shown label changes. A feature absent from the file keeps its "
                         "key, and the count is logged.")
    ap.add_argument("--modalities-only", action="store_true",
                    help="write ONLY the --modality sidecars into an existing bundle and update "
                         "its metadata. Nothing else is read or rewritten, so adding a modality "
                         "to a built bundle costs minutes instead of the hours a full rebuild takes.")
    a = ap.parse_args(argv)

    modality_sources = _parse_id_value(a.modality, "--modality")
    modality_labels = _parse_id_value(a.modality_label, "--modality-label")
    modality_display = {
        mid: _read_display_names(path)
        for mid, path in _parse_id_value(
            a.modality_display_names, "--modality-display-names").items()}
    modality_feature_labels = _parse_id_value(a.modality_feature_label, "--modality-feature-label")
    for key in list(modality_labels) + list(modality_feature_labels):
        if key not in modality_sources:
            raise SystemExit(f"--modality-label/--modality-feature-label names {key!r}, "
                             f"which no --modality declares")
    modality_deg_roots = _parse_id_value(a.deg_modality, "--deg-modality")
    if a.modalities_only and not (modality_sources or modality_deg_roots):
        raise SystemExit("--modalities-only needs at least one --modality id=path or "
                         "--deg-modality id=dir")

    if a.modalities_only:
        return _add_modalities_to_bundle(a, modality_sources, modality_labels,
                                         modality_feature_labels, modality_deg_roots,
                                         modality_display)
    if not a.h5ad:
        raise SystemExit("--h5ad is required unless --modalities-only is given")

    os.makedirs(a.out, exist_ok=True)
    paths = B.BundlePaths(a.out, a.prefix)
    warn_log: List[str] = []

    log(f"opening {os.path.abspath(a.h5ad)}")
    f = h5py.File(a.h5ad, "r")
    obs, var, uns = f["obs"], f["var"], f.get("uns")

    # ---- cells, genes ---------------------------------------------------------
    obs_index_key = obs.attrs.get("_index", "_index")
    if isinstance(obs_index_key, bytes):
        obs_index_key = obs_index_key.decode()
    barcodes = _decode(obs[obs_index_key][:])
    n_cells = len(barcodes)

    var_index_key = var.attrs.get("_index", "_index")
    if isinstance(var_index_key, bytes):
        var_index_key = var_index_key.decode()
    gene_ids = _decode(var[var_index_key][:])
    gene_syms = _decode(var["gene_symbols"][:]) if "gene_symbols" in var else list(gene_ids)
    n_genes = len(gene_ids)
    log(f"{n_cells:,} cells x {n_genes:,} genes")

    # ---- cell states ----------------------------------------------------------
    if a.cluster_key not in obs:
        raise RuntimeError(f"--cluster-key '{a.cluster_key}' not in obs. columns={list(obs.keys())}")
    kind, codes, cats = read_obs_column(obs, a.cluster_key)
    if kind == "numeric":
        raise RuntimeError(f"--cluster-key '{a.cluster_key}' is numeric; it must be categorical")
    if int((codes < 0).sum()):
        raise RuntimeError(f"{int((codes < 0).sum())} cells have no {a.cluster_key} value")
    present_states = cats
    log(f"{len(present_states)} cell states present in obs['{a.cluster_key}']")

    # ---- canonical order + colors --------------------------------------------
    order_file, colors_file = read_canonical_order(a.order)
    order_uns = read_uns_strlist(uns, "lineage_order") if uns is not None else None
    colors_uns = read_uns_json(uns, "cluster_colors_json") if uns is not None else None
    canonical = order_file or order_uns
    if not canonical:
        raise RuntimeError("no canonical order found: pass --order or store uns['lineage_order']")
    colors = dict(colors_uns or {})
    colors.update(colors_file or {})
    order_source = "--order file" if order_file else "uns['lineage_order']"

    ordered_states = [s for s in canonical if s in set(present_states)]
    extra = sorted(set(present_states) - set(canonical))
    if extra:
        warn_log.append(f"{len(extra)} state(s) present in obs but absent from the canonical order, "
                        f"appended alphabetically: {extra}")
        ordered_states += extra
    if len(ordered_states) != len(present_states):
        raise RuntimeError(f"ordered {len(ordered_states)} states but obs holds {len(present_states)}")
    log(f"canonical order from {order_source}: {len(canonical)} states, "
        f"{len(ordered_states)} of them present in this object")

    # remap obs codes to the canonical position
    pos = {s: i for i, s in enumerate(ordered_states)}
    remap = np.asarray([pos[c] for c in present_states], dtype=np.int16)
    state_code = remap[codes].astype(np.int16)
    n_states = len(ordered_states)
    state_n = np.bincount(state_code.astype(np.int64), minlength=n_states).astype(np.int64)
    if int(state_n.sum()) != n_cells:
        raise RuntimeError(f"state counts sum to {state_n.sum()}, expected {n_cells}")

    missing_color = [s for s in ordered_states if s not in colors]
    if missing_color:
        warn_log.append(f"{len(missing_color)} of {n_states} states have no colour; "
                        f"the client falls back to grey: {missing_color}")

    # ---- covariates -----------------------------------------------------------
    covariates: Dict[str, Dict] = {}
    cov_arrays: Dict[str, np.ndarray] = {}
    for name in obs.keys():
        if name == obs_index_key:
            continue
        try:
            k, vals, cs = read_obs_column(obs, name)
        except Exception as exc:                    # noqa: BLE001 - recorded, never silent
            warn_log.append(f"obs column '{name}' skipped: {type(exc).__name__}: {exc}")
            continue
        if k == "numeric":
            finite = np.isfinite(vals)
            covariates[name] = {"kind": "numeric", "n_missing": int((~finite).sum()),
                                "min": float(vals[finite].min()) if finite.any() else None,
                                "max": float(vals[finite].max()) if finite.any() else None}
            cov_arrays["cov_num_" + name] = vals.astype(np.float32)
        else:
            if len(cs) > 200:
                warn_log.append(f"obs column '{name}' has {len(cs)} categories (>200); not exposed")
                continue
            covariates[name] = {"kind": "categorical", "categories": cs,
                                "n_missing": int((vals < 0).sum())}
            cov_arrays["cov_cat_" + name] = vals.astype(np.int32)
    log(f"{len(covariates)} obs columns exposed as covariates")

    # ---- expression store + statistics ---------------------------------------
    if a.layer == "X":
        grp = f["X"]
    else:
        if "layers" not in f or a.layer not in f["layers"]:
            raise RuntimeError(f"layer '{a.layer}' not found. layers={list(f.get('layers', {}).keys())}")
        grp = f["layers"][a.layer]
    if grp.attrs.get("encoding-type") != "csr_matrix":
        raise RuntimeError(f"layer '{a.layer}' is {grp.attrs.get('encoding-type')}, expected csr_matrix")

    if a.skip_expr and not paths.missing():
        log("--skip-expr: reusing the existing expression store and statistics")
        mean_gs = np.load(paths.stats_mean)
        frac_gs = np.load(paths.stats_frac)
        nnz = int(np.load(paths.expr_indptr)[-1])
        gene_mean = mean_gs.astype(np.float64) @ (state_n / n_cells)
        gene_var = None
    else:
        t0 = time.time()
        sum_gs, cnt_gs, gene_sum, gene_sumsq, nnz = build_expression_store(
            grp, n_cells, n_genes, state_code, n_states, paths, a.expr_dtype, a.row_block)
        log(f"expression store built in {time.time() - t0:.1f}s")
        denom = state_n.astype(np.float64)[None, :]
        mean_gs = (sum_gs / denom).astype(np.float32)
        frac_gs = (cnt_gs / denom).astype(np.float32)
        np.save(paths.stats_mean, mean_gs)
        np.save(paths.stats_frac, frac_gs)
        np.save(paths.stats_n, state_n)
        gene_mean = gene_sum / n_cells
        gene_var = (gene_sumsq - n_cells * gene_mean ** 2) / max(n_cells - 1, 1)
        gene_var[gene_var < 0] = 0.0

    # ---- embedding ------------------------------------------------------------
    umap_warnings: List[str] = []
    if a.embedding_from and a.embedding_from.startswith("obsm:"):
        key = a.embedding_from.split(":", 1)[1]
        if "obsm" not in f or key not in f["obsm"]:
            raise RuntimeError(f"obsm['{key}'] not found. obsm={list(f.get('obsm', {}).keys())}")
        emb = np.asarray(f["obsm"][key][:, :2], dtype=np.float32)
        embedding_method = f"reused obsm['{key}'] from the source h5ad"
        log(embedding_method)
    else:
        if gene_var is None:
            raise RuntimeError("--skip-expr cannot compute an embedding: per-gene variance is gone. "
                               "Rerun without --skip-expr or pass --embedding-from obsm:<key>.")
        hvg = select_hvg(gene_mean, gene_var, a.n_hvg)
        log(f"selected {hvg.size} HVGs of {n_genes} by binned dispersion")
        t0 = time.time()
        Xh = dense_from_csc(paths, hvg, n_cells)
        log(f"HVG matrix materialised {Xh.shape} in {time.time() - t0:.1f}s "
            f"({Xh.nbytes / 2**30:.2f} GiB)")
        emb, _pcs, umap_warnings = compute_embedding(Xh, a.n_pcs, a.n_neighbors, a.min_dist, a.seed)
        del Xh
        embedding_method = (f"binned-dispersion HVG (n={hvg.size}) -> z-score, clip +/-10 -> "
                            f"randomized_svd PCA (n={a.n_pcs}) -> UMAP "
                            f"(n_neighbors={a.n_neighbors}, min_dist={a.min_dist}, seed={a.seed})")
    for w in umap_warnings:
        warn_log.append("UMAP " + w)

    # ---- write the cell sidecar ----------------------------------------------
    np.savez(paths.cells, embedding=emb, state_code=state_code, **cov_arrays)
    log(f"wrote {paths.cells}")

    # ---- legacy bundle files --------------------------------------------------
    with open(paths.umap, "w") as fh:
        fh.write("barcode\tUMAP1\tUMAP2\n")
        for i, bc in enumerate(barcodes):
            fh.write(f"{bc}\t{emb[i, 0]}\t{emb[i, 1]}\n")
    with open(paths.clusters, "w") as fh:
        fh.write(f"barcode\t{a.cluster_key}\tPopulation\n")
        for i, bc in enumerate(barcodes):
            s = ordered_states[int(state_code[i])]
            fh.write(f"{bc}\t{s}\t{s}\n")
    with open(paths.genes, "w") as fh:
        fh.write("index\tgene_id\tsymbol\n")
        for i, (gid, sym) in enumerate(zip(gene_ids, gene_syms)):
            fh.write(f"{i}\t{gid}\t{sym}\n")
    log(f"wrote {paths.umap}, {paths.clusters}, {paths.genes}")

    # ---- markers, DEG, CCC ----------------------------------------------------
    marker_info = ingest_markers(a.markers, paths.markers)
    deg_manifest = ingest_deg(a.deg, paths.deg_dir, paths.deg_manifest,
                              modality_roots=_parse_id_value(a.deg_modality, "--deg-modality"))
    ccc_info = ingest_ccc(a.ccc, paths.ccc)
    log(f"markers: {marker_info['n_rows']} rows; DEG: {len(deg_manifest['comparisons'])} tables; "
        f"CCC available: {ccc_info['available']}")
    if marker_info["available"]:
        unknown = [c for c in marker_info["clusters"] if c not in pos]
        if unknown:
            warn_log.append(f"{len(unknown)} of {len(marker_info['clusters'])} marker clusters do not "
                            f"match a cell state in obs: {unknown}")

    # ---- centroid matrix (<prefix>.txt) --------------------------------------
    sym_to_row: Dict[str, int] = {}
    for i, sym in enumerate(gene_syms):
        sym_to_row.setdefault(sym, i)
    chosen: List[int] = []
    seen = set()
    n_marker_hit = 0
    if marker_info["available"]:
        with open(paths.markers, "r") as fh:
            fh.readline()
            for ln in fh:
                g = ln.split("\t", 1)[0]
                r = sym_to_row.get(g)
                if r is not None:
                    n_marker_hit += 1
                    if r not in seen:
                        seen.add(r); chosen.append(r)
    n_marker_rows = len(chosen)
    if len(chosen) < a.max_centroid_genes:
        extra_order = np.argsort(-mean_gs.max(axis=1))
        for r in extra_order:
            if len(chosen) >= a.max_centroid_genes:
                break
            r = int(r)
            if r not in seen:
                seen.add(r); chosen.append(r)
    with open(paths.centroids, "w") as fh:
        fh.write("UID\t" + "\t".join(ordered_states) + "\n")
        for r in chosen:
            fh.write(gene_syms[r] + "\t" + "\t".join(f"{v:.6g}" for v in mean_gs[r]) + "\n")
    log(f"wrote {paths.centroids}: {len(chosen)} genes x {n_states} states "
        f"({n_marker_rows} unique marker genes matched)")

    # ---- modality stores ------------------------------------------------------
    modality_manifest: Dict[str, Dict] = {}
    for modality_id, source in modality_sources.items():
        modality_manifest[modality_id] = ingest_modality(
            modality_id, source, paths=paths, barcodes=barcodes, states=ordered_states,
            state_code=state_code, state_n=state_n,
            label=modality_labels.get(modality_id, ""),
            feature_label=modality_feature_labels.get(modality_id, "feature"),
            expr_dtype=a.expr_dtype,
            display_names=modality_display.get(modality_id),
        )
    if modality_manifest:
        log(f"modalities: {sorted(modality_manifest)}")

    # ---- metadata + config snippet -------------------------------------------
    ds_id = a.dataset_id or a.prefix
    meta = {
        "cluster_colors": {s: colors.get(s) for s in ordered_states if colors.get(s)},
        "cluster_key": a.cluster_key,
        "lineage_order": canonical,
        "n_cells": n_cells,
        "n_features": n_genes,
        "umap_key": "X_umap",
        "source_h5ad": os.path.abspath(a.h5ad),
        "reference_clusters_tsv": paths.clusters,
        "reference_coords_tsv": paths.umap,
        "scalable_viewer": {
            "bundle_version": B.BUNDLE_VERSION,
            "id": ds_id,
            "label": a.label or ds_id,
            # The LungMAP study this bundle belongs to. `study_ids_for_dataset`
            # (scalable_app.py) reads it, and shows no study record when it is empty.
            "study_id": (a.study_id or "").strip(),
            "prefix": a.prefix,
            "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "layer": a.layer,
            "expr_dtype": a.expr_dtype,
            "nnz": int(nnz),
            "n_states": n_states,
            "states": ordered_states,
            "state_n": state_n.tolist(),
            "canonical_order_source": order_source,
            "canonical_order_file": os.path.abspath(a.order) if a.order else None,
            "states_in_canonical_order_not_present": [s for s in canonical if s not in pos],
            "covariates": covariates,
            "embedding_method": embedding_method,
            "stats_method": ("mean = unweighted mean of the layer over the cells of a state; "
                             "frac = fraction of those cells with a value above zero"),
            "hvg_method": ("dispersion = variance / mean of the layer, 20 equal-count bins of mean, "
                           "z-scored inside each bin, top n kept"),
            "markers": marker_info,
            "deg": deg_manifest,
            "ccc": ccc_info,
            "modalities": modality_manifest,
            "centroid_genes": len(chosen),
            "centroid_marker_genes": n_marker_rows,
            "warnings": warn_log,
        },
    }
    with open(paths.metadata, "w") as fh:
        json.dump(meta, fh, indent=2)
    with open(paths.config_snippet, "w") as fh:
        json.dump({"cluster_key": a.cluster_key, "id": ds_id, "label": a.label or ds_id,
                   "reference_clusters_tsv": paths.clusters,
                   "reference_coords_tsv": paths.umap,
                   "ambient_options": ["default"]}, fh, indent=2)
    log(f"wrote {paths.metadata} and {paths.config_snippet}")

    missing = paths.missing()
    if missing:
        raise RuntimeError(f"bundle incomplete, missing: {missing}")

    log("=" * 70)
    log(f"BUNDLE OK: {paths.bundle_dir} (prefix {a.prefix})")
    log(f"cells {n_cells:,} | genes {n_genes:,} | states {n_states} | nnz {nnz:,}")
    if warn_log:
        log(f"WARNINGS: {len(warn_log)}")
        for w in warn_log:
            log("  ! " + w)
    else:
        log("WARNINGS: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
