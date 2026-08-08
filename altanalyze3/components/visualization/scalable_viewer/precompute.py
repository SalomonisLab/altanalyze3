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


def ingest_deg(deg_root: Optional[str], deg_dir: str, manifest_path: str) -> Dict:
    """Copy every DEG_detailed_*.tsv / DEG_pooled_overall_*.tsv found under deg_root.

    Nothing is recomputed. Each file is copied verbatim and indexed, so the browser
    shows the numbers the differential workflow produced.
    """
    manifest = {"comparisons": []}
    if not deg_root or not os.path.isdir(deg_root):
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        return manifest
    os.makedirs(deg_dir, exist_ok=True)
    for dirpath, _d, filenames in os.walk(deg_root):
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
            dst = os.path.join(deg_dir, fn)
            with open(src, "r") as fi:
                text = fi.read()
            with open(dst, "w") as fo:
                fo.write(text)
            lines = [ln for ln in text.split("\n") if ln.strip()]
            header = lines[0].split("\t") if lines else []
            manifest["comparisons"].append({
                "id": f"{comp}::{kind}", "comparison": comp, "kind": kind,
                "file": os.path.basename(dst), "path": os.path.abspath(dst),
                "source": os.path.abspath(src), "columns": header,
                "n_rows": max(len(lines) - 1, 0),
            })
    manifest["comparisons"].sort(key=lambda c: (c["comparison"], c["kind"]))
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build a scalable_viewer bundle from an h5ad.")
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", required=True, help="bundle directory (created if absent)")
    ap.add_argument("--prefix", required=True, help="file-name prefix inside the bundle")
    ap.add_argument("--label", default=None, help="human label for the catalog")
    ap.add_argument("--dataset-id", default=None, help="catalog id (default: --prefix)")
    ap.add_argument("--cluster-key", default="cell_state", help="obs column holding the cell state")
    ap.add_argument("--layer", default="lognorm", help="layer to serve ('X' for adata.X)")
    ap.add_argument("--markers", default=None)
    ap.add_argument("--order", default=None, help="canonical order TSV (order/cell_state/color)")
    ap.add_argument("--deg", default=None, help="directory tree holding DEG_*.tsv tables")
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
    a = ap.parse_args(argv)

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
    deg_manifest = ingest_deg(a.deg, paths.deg_dir, paths.deg_manifest)
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
