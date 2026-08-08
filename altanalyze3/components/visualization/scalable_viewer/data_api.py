"""Read side of the scalable_viewer.

Every method here reads a precomputed bundle. Nothing opens an h5ad, and nothing
recomputes a statistic that precompute.py already wrote.

Loading is lazy in two steps:
  1. Catalog build reads only <prefix>_metadata.json for each bundle. This stays fast
     with many datasets.
  2. The first request that needs an array memory-maps it. Memory-mapped arrays cost
     no resident memory until a page is touched, so a 3.9 GiB expression store opens
     in microseconds and a single gene read touches a few pages.
"""
from __future__ import annotations

import json
import math
import os
import struct
import threading
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import bundle as B


class Dataset:
    """One precomputed bundle."""

    def __init__(self, paths: B.BundlePaths):
        self.paths = paths
        missing = paths.missing()
        if missing:
            raise FileNotFoundError(f"bundle {paths.bundle_dir} is incomplete: {missing}")
        self.meta = B.read_metadata(paths.metadata)
        blk = B.viewer_block(self.meta)
        if blk is None:
            raise ValueError(f"{paths.metadata} has no 'scalable_viewer' block")
        self.sv = blk
        self.id: str = blk["id"]
        self.label: str = blk.get("label", self.id)
        self.states: List[str] = blk["states"]
        self.state_n: List[int] = blk["state_n"]
        self.cluster_key: str = self.meta.get("cluster_key", "cell_state")
        self.colors: Dict[str, str] = self.meta.get("cluster_colors", {})
        self.n_cells: int = int(self.meta["n_cells"])
        self.n_genes: int = int(self.meta["n_features"])
        self._lock = threading.Lock()
        self._cells = None
        self._indptr = None
        self._indices = None
        self._data = None
        self._mean = None
        self._frac = None
        self._symbols: Optional[List[str]] = None
        self._gene_ids: Optional[List[str]] = None
        self._sym_index: Optional[Dict[str, int]] = None
        self._markers: Optional[List[Dict]] = None

    # ------------------------------------------------------------ lazy loaders

    @property
    def cells(self):
        if self._cells is None:
            with self._lock:
                if self._cells is None:
                    self._cells = np.load(self.paths.cells)
        return self._cells

    @property
    def state_code(self) -> np.ndarray:
        return self.cells["state_code"]

    @property
    def embedding(self) -> np.ndarray:
        return self.cells["embedding"]

    @property
    def indptr(self) -> np.ndarray:
        if self._indptr is None:
            with self._lock:
                if self._indptr is None:
                    self._indptr = np.load(self.paths.expr_indptr, mmap_mode="r")
        return self._indptr

    @property
    def indices(self) -> np.ndarray:
        if self._indices is None:
            with self._lock:
                if self._indices is None:
                    self._indices = np.load(self.paths.expr_indices, mmap_mode="r")
        return self._indices

    @property
    def data(self) -> np.ndarray:
        if self._data is None:
            with self._lock:
                if self._data is None:
                    self._data = np.load(self.paths.expr_data, mmap_mode="r")
        return self._data

    @property
    def stats_mean(self) -> np.ndarray:
        if self._mean is None:
            with self._lock:
                if self._mean is None:
                    self._mean = np.load(self.paths.stats_mean, mmap_mode="r")
        return self._mean

    @property
    def stats_frac(self) -> np.ndarray:
        if self._frac is None:
            with self._lock:
                if self._frac is None:
                    self._frac = np.load(self.paths.stats_frac, mmap_mode="r")
        return self._frac

    def _load_genes(self) -> None:
        if self._symbols is not None:
            return
        with self._lock:
            if self._symbols is not None:
                return
            syms: List[str] = []
            gids: List[str] = []
            with open(self.paths.genes, "r") as fh:
                fh.readline()
                for ln in fh:
                    parts = ln.rstrip("\n").split("\t")
                    gids.append(parts[1])
                    syms.append(parts[2] if len(parts) > 2 else parts[1])
            idx: Dict[str, int] = {}
            for i, s in enumerate(syms):
                idx.setdefault(s.upper(), i)
            for i, g in enumerate(gids):
                idx.setdefault(g.upper(), i)
            self._symbols, self._gene_ids, self._sym_index = syms, gids, idx

    @property
    def symbols(self) -> List[str]:
        self._load_genes()
        return self._symbols                                    # type: ignore[return-value]

    # ------------------------------------------------------------------ genes

    def resolve_gene(self, name: str) -> Optional[int]:
        self._load_genes()
        return self._sym_index.get((name or "").strip().upper())  # type: ignore[union-attr]

    def search_genes(self, q: str, limit: int = 25) -> List[str]:
        self._load_genes()
        q = (q or "").strip().upper()
        syms = self._symbols                                    # type: ignore[assignment]
        if not q:
            return sorted(syms)[:limit]
        pre = [s for s in syms if s.upper().startswith(q)]
        if len(pre) >= limit:
            return sorted(set(pre))[:limit]
        sub = [s for s in syms if q in s.upper() and not s.upper().startswith(q)]
        return sorted(set(pre))[:limit] + sorted(set(sub))[: max(0, limit - len(set(pre)))]

    def gene_column(self, row: int) -> Tuple[np.ndarray, np.ndarray]:
        """Return (cell indices, values) for one gene. One memmap slice, no scan."""
        s, e = int(self.indptr[row]), int(self.indptr[row + 1])
        if e <= s:
            return np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.float32)
        return (np.asarray(self.indices[s:e], dtype=np.uint32),
                np.asarray(self.data[s:e], dtype=np.float32))

    # ------------------------------------------------------------- covariates

    def covariate_names(self) -> Dict[str, Dict]:
        return self.sv.get("covariates", {})

    def covariate_values(self, name: str) -> Tuple[str, np.ndarray, Optional[List[str]]]:
        info = self.covariate_names().get(name)
        if info is None:
            raise KeyError(name)
        if info["kind"] == "numeric":
            return "numeric", self.cells["cov_num_" + name], None
        return "categorical", self.cells["cov_cat_" + name], info["categories"]

    # ------------------------------------------------------------------ plots

    def violin(self, row: int, n_bins: int = 40) -> Dict:
        """Per-state distribution of one gene, computed from the sparse column plus the
        known per-state cell counts. Zeros are never dropped: they are the difference
        between the state size and the number of stored non-zeros."""
        idx, val = self.gene_column(row)
        codes = self.state_code
        n_states = len(self.states)
        vmax = float(val.max()) if val.size else 0.0
        edges = np.linspace(0.0, max(vmax, 1e-6), n_bins + 1)

        st = codes[idx.astype(np.int64)].astype(np.int64) if idx.size else np.empty(0, np.int64)
        out = []
        for s in range(n_states):
            n_tot = int(self.state_n[s])
            sel = val[st == s] if idx.size else np.empty(0, np.float32)
            n_nz = int(sel.size)
            n_zero = n_tot - n_nz
            if n_zero < 0:
                raise RuntimeError(f"state {self.states[s]}: {n_nz} non-zeros exceed {n_tot} cells")
            full = np.concatenate([sel, np.zeros(n_zero, dtype=np.float32)]) if n_zero else sel
            hist = np.histogram(full, bins=edges)[0].astype(int).tolist() if n_tot else []
            q = (np.percentile(full, [5, 25, 50, 75, 95]).tolist() if n_tot else [0, 0, 0, 0, 0])
            out.append({
                "state": self.states[s], "n": n_tot, "n_detected": n_nz,
                "frac": (n_nz / n_tot) if n_tot else 0.0,
                "mean": float(full.mean()) if n_tot else 0.0,
                "p05": q[0], "q25": q[1], "median": q[2], "q75": q[3], "p95": q[4],
                "max": float(full.max()) if n_tot else 0.0,
                "hist": hist,
            })
        return {"gene": self.symbols[row], "bin_edges": edges.tolist(), "states": out,
                "value_max": vmax, "layer": self.sv.get("layer")}

    def dotplot(self, rows: List[int]) -> Dict:
        """Mean and detected fraction per (gene, state), read straight from the stats matrices."""
        m = np.asarray(self.stats_mean[rows, :], dtype=np.float32)
        fr = np.asarray(self.stats_frac[rows, :], dtype=np.float32)
        return {"genes": [self.symbols[r] for r in rows], "states": self.states,
                "state_n": self.state_n, "mean": m.tolist(), "frac": fr.tolist(),
                "layer": self.sv.get("layer")}

    # ---------------------------------------------------------------- markers

    def markers(self) -> List[Dict]:
        if self._markers is not None:
            return self._markers
        with self._lock:
            if self._markers is not None:
                return self._markers
            rows: List[Dict] = []
            if os.path.isfile(self.paths.markers):
                with open(self.paths.markers, "r") as fh:
                    fh.readline()
                    for ln in fh:
                        p = ln.rstrip("\n").split("\t")
                        if len(p) < 4:
                            continue
                        rows.append({"gene": p[0], "cluster": p[1],
                                     "fold": _f(p[2]), "p": _f(p[3])})
            self._markers = rows
            return rows

    def marker_genes_for(self, states: Optional[List[str]], top: int) -> List[str]:
        """Top markers per state in CANONICAL state order, de-duplicated."""
        want = self.states if not states else [s for s in self.states if s in set(states)]
        by_state: Dict[str, List[Dict]] = {s: [] for s in want}
        for m in self.markers():
            if m["cluster"] in by_state:
                by_state[m["cluster"]].append(m)
        out: List[str] = []
        seen = set()
        for s in want:                                          # canonical order, never alphabetical
            ms = sorted(by_state[s], key=lambda r: (-(r["fold"] if r["fold"] is not None else -1e9)))
            for m in ms[:top]:
                if m["gene"] not in seen:
                    seen.add(m["gene"]); out.append(m["gene"])
        return out

    def heatmap_tsv(self, states: Optional[List[str]], top: int, row_scale: bool) -> Tuple[str, Dict]:
        """A Morpheus-ready TSV: genes x cell states, columns in canonical order."""
        cols = self.states if not states else [s for s in self.states if s in set(states)]
        col_idx = [self.states.index(c) for c in cols]
        genes = self.marker_genes_for(states, top)
        rows, missing = [], []
        for g in genes:
            r = self.resolve_gene(g)
            if r is None:
                missing.append(g); continue
            rows.append((g, r))
        lines = ["Name\t" + "\t".join(cols)]
        for g, r in rows:
            v = np.asarray(self.stats_mean[r, :], dtype=np.float64)[col_idx]
            if row_scale:
                sd = v.std()
                v = (v - v.mean()) / sd if sd > 0 else v * 0.0
            lines.append(g + "\t" + "\t".join(f"{x:.5g}" for x in v))
        info = {"n_genes_requested": len(genes), "n_genes_written": len(rows),
                "n_genes_missing": len(missing), "missing": missing[:20],
                "n_states": len(cols), "row_scale": row_scale}
        return "\n".join(lines) + "\n", info

    # -------------------------------------------------------------------- DEG

    def deg_manifest(self) -> Dict:
        return self.sv.get("deg", {"comparisons": []})

    def deg_table(self, comp_id: str, max_rows: int, fdr_max: Optional[float],
                  state: Optional[str]) -> Dict:
        entries = {c["id"]: c for c in self.deg_manifest().get("comparisons", [])}
        c = entries.get(comp_id)
        if c is None:
            raise KeyError(comp_id)
        path = os.path.join(self.paths.deg_dir, c["file"])
        with open(path, "r") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            raw = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        low = [h.strip().lower() for h in header]

        def col(*names):
            for n in names:
                if n in low:
                    return low.index(n)
            return None

        i_g = col("gene", "uid", "symbol")
        i_fc = col("log2fc", "fold", "logfc")
        i_fdr = col("fdr", "adjp", "padj", "fdr p-value")
        i_p = col("pval", "p", "pvalue")
        i_pop = col("population", "cluster", "cell_state")

        n_total = len(raw)
        recs = []
        n_no_fdr = 0
        for r in raw:
            g = r[i_g] if (i_g is not None and i_g < len(r)) else ""
            fc = _f(r[i_fc]) if (i_fc is not None and i_fc < len(r)) else None
            fdr = _f(r[i_fdr]) if (i_fdr is not None and i_fdr < len(r)) else None
            pv = _f(r[i_p]) if (i_p is not None and i_p < len(r)) else None
            pop = r[i_pop] if (i_pop is not None and i_pop < len(r)) else None
            if fdr is None:
                n_no_fdr += 1
            recs.append({"gene": g, "log2fc": fc, "fdr": fdr, "pval": pv, "population": pop})

        n_after_state = len(recs)
        if state and i_pop is not None:
            recs = [r for r in recs if r["population"] == state]
            n_after_state = len(recs)
        n_after_fdr = n_after_state
        if fdr_max is not None:
            recs = [r for r in recs if r["fdr"] is not None and r["fdr"] <= fdr_max]
            n_after_fdr = len(recs)

        recs.sort(key=lambda r: (r["fdr"] if r["fdr"] is not None else 2.0,
                                 -abs(r["log2fc"] or 0.0)))
        shown = recs[:max_rows]
        volcano = [{"gene": r["gene"], "x": r["log2fc"],
                    "y": (-math.log10(r["fdr"]) if (r["fdr"] and r["fdr"] > 0) else None),
                    "population": r["population"]}
                   for r in recs if r["log2fc"] is not None]
        return {"id": comp_id, "comparison": c["comparison"], "kind": c["kind"],
                "source": c.get("source"), "path": path, "columns": header,
                "n_rows_file": n_total, "n_rows_after_filters": len(recs),
                "n_rows_returned": len(shown), "n_rows_no_fdr": n_no_fdr,
                "populations": sorted({r["population"] for r in raw_pops(raw, i_pop)}) if i_pop is not None else [],
                "rows": shown, "volcano": volcano[:20000]}

    # -------------------------------------------------------------------- CCC

    def ccc(self, limit: int = 5000) -> Dict:
        info = self.sv.get("ccc", {"available": False})
        if not info.get("available") or not os.path.isfile(self.paths.ccc):
            return {"available": False, "reason": "no cell-cell communication table in this bundle",
                    "rows": [], "columns": [], "n_rows": 0}
        with open(self.paths.ccc, "r") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = [ln.rstrip("\n").split("\t") for _, ln in zip(range(limit), fh) if ln.strip()]
        return {"available": True, "columns": header,
                "rows": [dict(zip(header, r)) for r in rows],
                "n_rows": info.get("n_rows", len(rows)), "source": info.get("source")}

    # ------------------------------------------------------------- catalog row

    def summary(self) -> Dict:
        return {"id": self.id, "label": self.label, "n_cells": self.n_cells,
                "n_genes": self.n_genes, "n_states": len(self.states),
                "cluster_key": self.cluster_key, "layer": self.sv.get("layer"),
                "built_utc": self.sv.get("built_utc"),
                "bundle_dir": self.paths.bundle_dir, "prefix": self.paths.prefix,
                "has_markers": bool(self.sv.get("markers", {}).get("available")),
                "n_deg_tables": len(self.deg_manifest().get("comparisons", [])),
                "has_ccc": bool(self.sv.get("ccc", {}).get("available")),
                "n_warnings": len(self.sv.get("warnings", []))}


def raw_pops(raw, i_pop):
    for r in raw:
        if i_pop < len(r):
            yield {"population": r[i_pop]}


def _f(s) -> Optional[float]:
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def pack_sparse(idx: np.ndarray, val: np.ndarray) -> bytes:
    """[n: int32 LE][idx: uint32 * n][val: float32 * n] - parsed by a DataView in the client."""
    return (struct.pack("<i", int(idx.size))
            + np.ascontiguousarray(idx, dtype="<u4").tobytes()
            + np.ascontiguousarray(val, dtype="<f4").tobytes())


class Catalog:
    """Datasets discovered once at startup; array pages are mapped on first use."""

    def __init__(self, bundles: List[B.BundlePaths]):
        self._bundles = bundles
        self._loaded: Dict[str, Dataset] = {}
        self._entries: List[Dict] = []
        self._paths_by_id: Dict[str, B.BundlePaths] = {}
        self._lock = threading.Lock()
        self.load_errors: List[Dict] = []
        for bp in bundles:
            try:
                meta = B.read_metadata(bp.metadata)
                blk = B.viewer_block(meta) or {}
                missing = bp.missing()
                if missing:
                    self.load_errors.append({"bundle_dir": bp.bundle_dir, "prefix": bp.prefix,
                                             "error": f"incomplete bundle, missing {missing}"})
                    continue
                ds_id = blk.get("id", bp.prefix)
                if ds_id in self._paths_by_id:
                    self.load_errors.append({
                        "bundle_dir": bp.bundle_dir, "prefix": bp.prefix,
                        "error": f"duplicate dataset id '{ds_id}'; already claimed by "
                                 f"{self._paths_by_id[ds_id].bundle_dir}"})
                    continue
                self._paths_by_id[ds_id] = bp
                self._entries.append({
                    "id": ds_id, "label": blk.get("label", bp.prefix),
                    "prefix": bp.prefix, "bundle_dir": bp.bundle_dir,
                    "n_cells": meta.get("n_cells"), "n_genes": meta.get("n_features"),
                    "n_states": blk.get("n_states"), "cluster_key": meta.get("cluster_key"),
                    "layer": blk.get("layer"), "built_utc": blk.get("built_utc"),
                    "bundle_version": blk.get("bundle_version"),
                    "has_markers": bool(blk.get("markers", {}).get("available")),
                    "n_deg_tables": len(blk.get("deg", {}).get("comparisons", [])),
                    "has_ccc": bool(blk.get("ccc", {}).get("available")),
                    "loaded": False,
                })
            except Exception as exc:                            # noqa: BLE001 - recorded, never silent
                self.load_errors.append({"bundle_dir": bp.bundle_dir, "prefix": bp.prefix,
                                         "error": f"{type(exc).__name__}: {exc}"})

    @property
    def entries(self) -> List[Dict]:
        for e in self._entries:
            e["loaded"] = e["id"] in self._loaded
        return self._entries

    def ids(self) -> List[str]:
        return [e["id"] for e in self._entries]

    def get(self, ds_id: str) -> Dataset:
        """Load a dataset on first request and keep it. This is the lazy step."""
        if ds_id in self._loaded:
            return self._loaded[ds_id]
        with self._lock:
            if ds_id in self._loaded:
                return self._loaded[ds_id]
            bp = self._paths_by_id.get(ds_id)
            if bp is None:
                raise KeyError(ds_id)
            ds = Dataset(bp)
            self._loaded[ds_id] = ds
            return ds


def build_catalog(root: Optional[str] = None, catalog_file: Optional[str] = None) -> Catalog:
    bundles: List[B.BundlePaths] = []
    if catalog_file:
        bundles += B.load_catalog_file(catalog_file)
    if root:
        bundles += B.discover(root)
    seen, uniq = set(), []
    for bp in bundles:
        key = (bp.bundle_dir, bp.prefix)
        if key not in seen:
            seen.add(key); uniq.append(bp)
    return Catalog(uniq)
