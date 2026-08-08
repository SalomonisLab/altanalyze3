"""Bundle layout for the scalable_viewer.

The bundle EXTENDS the existing cellHarmony reference bundle used at
/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/flask/references/Human/Lung/HLCA/
(files Hs-Lung-HLCA.txt, Hs-Lung-HLCA_clusters.tsv, Hs-Lung-HLCA_umap.tsv,
Hs-Lung-HLCA_metadata.json, Hs-Lung-HLCA_config_snippet.json).

Every legacy file keeps its name, its columns and its meaning. The viewer adds
binary sidecars so the server never scans an h5ad at request time.

Legacy files (unchanged contract)
  <prefix>.txt                     centroid matrix. Row 1 = "UID" + cell-state names in
                                   CANONICAL order. Column 1 = gene symbol.
  <prefix>_clusters.tsv            barcode <cluster_key> Population
  <prefix>_umap.tsv                barcode UMAP1 UMAP2
  <prefix>_metadata.json           cluster_colors, cluster_key, lineage_order, n_cells,
                                   n_features, umap_key, source_h5ad ... plus a new
                                   "scalable_viewer" block (legacy readers ignore it).
  <prefix>_config_snippet.json     cellHarmony reference registration snippet

Viewer sidecars (new)
  <prefix>_cells.npz               embedding (N,2) f32, state code (N,) i16, covariate codes
  <prefix>_genes.tsv               row index -> gene id + symbol
  <prefix>_stats_mean.npy          (G,S) f32  mean of layers['lognorm'] per gene per state
  <prefix>_stats_frac.npy          (G,S) f32  fraction of cells with value > 0
  <prefix>_stats_n.npy             (S,)  i64  cells per state (the denominator)
  <prefix>_expr_indptr.npy         (G+1,) i64  gene-major (CSC) index
  <prefix>_expr_indices.npy        (nnz,) u32  cell index
  <prefix>_expr_data.npy           (nnz,) f32  lognorm value
  <prefix>_markers.tsv             marker table (gene, cluster, fold, p ...)
  <prefix>_deg/<comparison>.tsv    one precomputed differential table per comparison
  <prefix>_deg_manifest.json       comparison id -> file, columns, row count
  <prefix>_ccc.tsv                 cell-cell communication edges (optional)

Nothing here reads an h5ad. precompute.py writes the bundle; data_api.py reads it.
"""
from __future__ import annotations

import json
import os
from typing import Dict, List, Optional

# Marker of a bundle the viewer can serve. Bumped only on an incompatible layout change.
BUNDLE_VERSION = 1

# Legacy suffixes, kept byte-for-byte compatible with the HLCA bundle.
SUF_CENTROIDS = ".txt"
SUF_CLUSTERS = "_clusters.tsv"
SUF_UMAP = "_umap.tsv"
SUF_METADATA = "_metadata.json"
SUF_CONFIG = "_config_snippet.json"

# Viewer sidecars.
SUF_CELLS = "_cells.npz"
SUF_GENES = "_genes.tsv"
SUF_STATS_MEAN = "_stats_mean.npy"
SUF_STATS_FRAC = "_stats_frac.npy"
SUF_STATS_N = "_stats_n.npy"
SUF_EXPR_INDPTR = "_expr_indptr.npy"
SUF_EXPR_INDICES = "_expr_indices.npy"
SUF_EXPR_DATA = "_expr_data.npy"
SUF_MARKERS = "_markers.tsv"
SUF_DEG_DIR = "_deg"
SUF_DEG_MANIFEST = "_deg_manifest.json"
SUF_CCC = "_ccc.tsv"


class BundlePaths:
    """Every path in one bundle. All paths are absolute."""

    def __init__(self, bundle_dir: str, prefix: str):
        self.bundle_dir = os.path.abspath(bundle_dir)
        self.prefix = prefix

    def _p(self, suffix: str) -> str:
        return os.path.join(self.bundle_dir, self.prefix + suffix)

    @property
    def centroids(self) -> str: return self._p(SUF_CENTROIDS)

    @property
    def clusters(self) -> str: return self._p(SUF_CLUSTERS)

    @property
    def umap(self) -> str: return self._p(SUF_UMAP)

    @property
    def metadata(self) -> str: return self._p(SUF_METADATA)

    @property
    def config_snippet(self) -> str: return self._p(SUF_CONFIG)

    @property
    def cells(self) -> str: return self._p(SUF_CELLS)

    @property
    def genes(self) -> str: return self._p(SUF_GENES)

    @property
    def stats_mean(self) -> str: return self._p(SUF_STATS_MEAN)

    @property
    def stats_frac(self) -> str: return self._p(SUF_STATS_FRAC)

    @property
    def stats_n(self) -> str: return self._p(SUF_STATS_N)

    @property
    def expr_indptr(self) -> str: return self._p(SUF_EXPR_INDPTR)

    @property
    def expr_indices(self) -> str: return self._p(SUF_EXPR_INDICES)

    @property
    def expr_data(self) -> str: return self._p(SUF_EXPR_DATA)

    @property
    def markers(self) -> str: return self._p(SUF_MARKERS)

    @property
    def deg_dir(self) -> str: return self._p(SUF_DEG_DIR)

    @property
    def deg_manifest(self) -> str: return self._p(SUF_DEG_MANIFEST)

    @property
    def ccc(self) -> str: return self._p(SUF_CCC)

    def required(self) -> List[str]:
        """Files that must exist before the server will serve this bundle."""
        return [self.metadata, self.cells, self.genes, self.stats_mean, self.stats_frac,
                self.stats_n, self.expr_indptr, self.expr_indices, self.expr_data]

    def missing(self) -> List[str]:
        return [p for p in self.required() if not os.path.isfile(p)]


def read_metadata(path: str) -> Dict:
    with open(path, "r") as fh:
        return json.load(fh)


def discover(root: str) -> List[BundlePaths]:
    """Find every viewer bundle under root.

    A directory is a bundle when it holds a *_metadata.json whose JSON carries a
    "scalable_viewer" block. Directories without that block (for example a plain
    cellHarmony reference) are skipped, so pointing the server at a shared reference
    tree is safe.
    """
    found: List[BundlePaths] = []
    root = os.path.abspath(root)
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in sorted(filenames):
            if not fn.endswith(SUF_METADATA):
                continue
            full = os.path.join(dirpath, fn)
            try:
                meta = read_metadata(full)
            except (OSError, ValueError):
                continue
            if not isinstance(meta, dict) or "scalable_viewer" not in meta:
                continue
            found.append(BundlePaths(dirpath, fn[: -len(SUF_METADATA)]))
    return found


def load_catalog_file(path: str) -> List[BundlePaths]:
    """Read an explicit catalog JSON: {"datasets": [{"bundle_dir": ..., "prefix": ...}, ...]}."""
    with open(path, "r") as fh:
        cat = json.load(fh)
    out: List[BundlePaths] = []
    for entry in cat.get("datasets", []):
        bd = entry["bundle_dir"]
        if not os.path.isabs(bd):
            bd = os.path.join(os.path.dirname(os.path.abspath(path)), bd)
        out.append(BundlePaths(bd, entry["prefix"]))
    return out


def viewer_block(meta: Dict) -> Optional[Dict]:
    blk = meta.get("scalable_viewer")
    return blk if isinstance(blk, dict) else None
