#!/usr/bin/env python3
"""Unified multi-modality UMAP visualization (gene / isoform / junction / genomic-variant).

Two capabilities:

1. **UMAP propagation** -- assign UMAP coordinates once and share them across every modality h5ad
   for the SAME cell barcodes, so gene/isoform/junction/variant views are all in the same embedding.
   Coordinates come either from an h5ad that already has ``obsm['X_umap']`` or from cellHarmony-style
   reference projection via :func:`approximate_umap.approximate_umap`.

2. **Feature plotting** -- for any modality, plot a feature (or the top-N features of a gene) on the
   UMAP, either as an **expression scatter** (per-cell colour, q99 vmax, high-on-top) or as a
   **density contour** (expression-weighted KDE; filled + lines). Each modality has its own default
   colour so panels are visually distinguishable:

       gene     -> Greens     isoform  -> Reds
       junction -> Blues      variant  -> Purples

Feature ids are matched by a per-modality key:
  - gene:     var_names start with ``<ENSG>``               (whole-gene = sum of its columns)
  - isoform:  var_names ``<ENSG>:<isoform_id>``             (one column per isoform)
  - junction: var_names ``<ENSG>:Ea-Eb=chr:pos-pos``        (one column per junction)
  - variant:  var_names contain the variant id (chr:pos / rsid); matched by substring/prefix

Designed to subsume the per-modality plotting scripts (grey->red isoform tiles, grey->blue junction
tiles, red expression-density contours).
"""

from __future__ import annotations

import os
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------- palettes ----
# grey -> colour ramps (grey low end keeps non-expressing cells visible / neutral).
_RAMPS = {
    "gene":     ["#d9d9d9", "#e5f5e0", "#a1d99b", "#41ab5d", "#006d2c"],   # greens
    "isoform":  ["#d9d9d9", "#fee0d2", "#fc9272", "#ef3b2c", "#a50f15"],   # reds
    "junction": ["#d9d9d9", "#deebf7", "#9ecae1", "#4292c6", "#08519c"],   # blues
    "variant":  ["#d9d9d9", "#efedf5", "#bcbddc", "#807dba", "#54278f"],   # purples
}
MODALITY_CMAP = {m: LinearSegmentedColormap.from_list(f"grey_{m}", r) for m, r in _RAMPS.items()}
# matplotlib named cmap (for contourf, which wants a colormap with a clean low end)
MODALITY_CONTOUR_CMAP = {"gene": "Greens", "isoform": "Reds", "junction": "Blues", "variant": "Purples"}

DOT = 4.2


# ------------------------------------------------------- UMAP propagation ----
def assign_umap(
    primary_h5ad: str,
    reference: Optional[str] = None,
    reference_clusters: Optional[str] = None,
    query_cluster_key: str = "cluster",
    reference_cluster_key: Optional[str] = None,
    umap_key: str = "X_umap",
    propagate_to: Optional[Sequence[str]] = None,
    write: bool = True,
    log=print,
) -> Dict[str, np.ndarray]:
    """Ensure ``primary_h5ad`` has ``obsm[umap_key]`` and propagate those coords to other modality
    h5ads by shared cell barcode.

    If the primary already carries ``obsm[umap_key]`` it is used as-is. Otherwise, when a reference
    UMAP/clusters TSV pair is given, coordinates are projected with ``approximate_umap`` (the cell
    must have ``obs[query_cluster_key]``). The resulting barcode->coord map is applied to every h5ad
    in ``propagate_to`` (restricted to matched barcodes), so all modalities share one embedding.

    Returns the barcode -> (x, y) coordinate map.
    """
    a = ad.read_h5ad(primary_h5ad)
    if umap_key in a.obsm:
        log(f"[umap] using existing {umap_key} in {os.path.basename(primary_h5ad)} ({a.n_obs} cells)")
    else:
        if not (reference and reference_clusters):
            raise ValueError(
                f"{os.path.basename(primary_h5ad)} has no obsm['{umap_key}'] and no reference "
                f"(reference + reference_clusters) was provided to project one."
            )
        from .approximate_umap import approximate_umap, _load_reference_from_tsv
        ref_key = reference_cluster_key or "Population"
        ref = _load_reference_from_tsv(reference, reference_clusters,
                                       umap_key=umap_key, cluster_key=ref_key)
        if query_cluster_key not in a.obs.columns:
            raise ValueError(f"primary h5ad lacks obs['{query_cluster_key}'] needed for projection.")
        res = approximate_umap(query=a, reference=ref, query_cluster_key=query_cluster_key,
                               reference_cluster_key=ref_key, umap_key=umap_key)
        a = res.query_adata
        log(f"[umap] projected {umap_key} for {a.n_obs} cells via approximate_umap")
        if write:
            a.write_h5ad(primary_h5ad, compression="gzip")

    coords = np.asarray(a.obsm[umap_key])
    bc2xy = dict(zip(a.obs_names.astype(str), coords))

    for other in (propagate_to or []):
        o = ad.read_h5ad(other)
        keep = np.array([str(b) in bc2xy for b in o.obs_names])
        if keep.sum() == 0:
            log(f"[umap] WARN: 0 shared barcodes with {os.path.basename(other)}; skipped")
            continue
        o = o[keep].copy()
        o.obsm[umap_key] = np.array([bc2xy[str(b)] for b in o.obs_names])
        if write:
            o.write_h5ad(other, compression="gzip")
        log(f"[umap] propagated {umap_key} to {os.path.basename(other)} "
            f"({int(keep.sum())}/{len(keep)} barcodes)")
    return bc2xy


# --------------------------------------------------------- feature lookup ----
def _feature_columns(var_names: np.ndarray, modality: str, feature: str) -> list:
    """Column indices for a feature within a modality.

    gene: all columns whose var starts with '<feature>:' or equal '<feature>' (whole-gene sum).
    isoform/junction: columns whose var starts with '<feature>:' (a gene id) OR equal the exact id.
    variant: columns containing the feature id as a substring.
    """
    vn = var_names
    if modality == "variant":
        return [i for i, v in enumerate(vn) if feature in v]
    pref = feature + ":"
    cols = [i for i, v in enumerate(vn) if v == feature or v.startswith(pref)]
    return cols


# ------------------------------------------------------------- plotting ------
def _expr_of(X, cols) -> np.ndarray:
    sub = X[:, cols]
    e = np.asarray(sub.sum(axis=1)).ravel() if len(cols) > 1 else np.asarray(sub.todense()).ravel()
    return e


def plot_feature(
    ax, U: np.ndarray, expr: np.ndarray, modality: str, mode: str = "expression",
    title: str = "", thr: float = 0.08, bw: float = 0.25, nlev: int = 25,
    fill_alpha: float = 0.5, line_alpha: float = 0.7, xlim=None, ylim=None,
):
    """Draw one feature on a UMAP axis. mode='expression' (scatter) or 'contour' (density)."""
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    if mode == "expression":
        nz = expr[expr > 0]
        vmax = float(np.quantile(nz, 0.99)) if nz.size else 1.0
        o = np.argsort(expr)  # high on top
        sc = ax.scatter(U[o, 0], U[o, 1], c=expr[o], cmap=MODALITY_CMAP[modality],
                        vmin=0, vmax=vmax, s=DOT, linewidths=0)
        ax.set_title(title, fontsize=8)
        return sc
    elif mode == "contour":
        from scipy.stats import gaussian_kde
        ax.scatter(U[:, 0], U[:, 1], c="#dddddd", s=2.0, linewidths=0, zorder=1)
        m = expr > 0
        sc = None
        if m.sum() >= 10:
            gx = np.linspace(U[:, 0].min(), U[:, 0].max(), 200)
            gy = np.linspace(U[:, 1].min(), U[:, 1].max(), 200)
            GX, GY = np.meshgrid(gx, gy)
            try:
                kde = gaussian_kde(U[m].T, weights=expr[m], bw_method=bw)
                Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
                Z = Z / Z.max()
                Zm = np.ma.masked_less(Z, thr)
                levels = np.linspace(thr, 1.0, nlev)
                cm = MODALITY_CONTOUR_CMAP[modality]
                ax.contourf(GX, GY, Zm, levels=levels, cmap=cm, alpha=fill_alpha, zorder=2)
                sc = ax.contour(GX, GY, Z, levels=levels, cmap=cm, linewidths=0.6,
                                alpha=line_alpha, zorder=3)
            except Exception:
                pass
        ax.set_title(title, fontsize=8)
        return sc
    raise ValueError(f"mode must be 'expression' or 'contour', got {mode!r}")


def plot_gene_panels(
    h5ad: str, gene: str, modality: str, outdir: str, *,
    label: str = "", mode: str = "expression", top_n: int = 10, umap_key: str = "X_umap",
    gene_symbol: Optional[str] = None, log=print,
) -> Optional[str]:
    """Tile the top-N features of one gene for a modality. Returns the PDF path.

    For ``modality='gene'`` the single whole-gene track is plotted (top_n ignored). For isoform /
    junction / variant, the gene's top-N features (by total reads) are tiled.
    """
    a = ad.read_h5ad(h5ad)
    if umap_key not in a.obsm:
        raise ValueError(f"{os.path.basename(h5ad)} lacks obsm['{umap_key}']; run assign_umap first.")
    U = np.asarray(a.obsm[umap_key])
    X = a.X.tocsc() if sp.issparse(a.X) else sp.csc_matrix(a.X)
    var = np.array(a.var_names)
    xpad = (U[:, 0].max() - U[:, 0].min()) * 0.03
    ypad = (U[:, 1].max() - U[:, 1].min()) * 0.03
    xlim = (U[:, 0].min() - xpad, U[:, 0].max() + xpad)
    ylim = (U[:, 1].min() - ypad, U[:, 1].max() + ypad)
    sym = gene_symbol or gene
    os.makedirs(outdir, exist_ok=True)

    if modality == "gene":
        cols = _feature_columns(var, "gene", gene)
        if not cols:
            log(f"[{modality}] {sym}: not found; skipped"); return None
        expr = _expr_of(X, cols)
        fig, ax = plt.subplots(figsize=(7, 6))
        sc = plot_feature(ax, U, expr, modality, mode=mode,
                          title=f"{sym} ({gene})  {int(expr.sum()):,} reads, {int((expr>0).sum())} cells",
                          xlim=xlim, ylim=ylim)
        if sc is not None:
            cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label(f"{modality} ({'q99' if mode=='expression' else 'density'})", fontsize=8)
        out = os.path.join(outdir, f"{modality.upper()}_{label}_{sym}_{gene}.pdf")
    else:
        cols = _feature_columns(var, modality, gene)
        if not cols:
            log(f"[{modality}] {sym}: no features; skipped"); return None
        tot = np.array([X[:, c].sum() for c in cols])
        order = [cols[k] for k in np.argsort(tot)[::-1][:top_n]]
        n = len(order)
        ncol = 5; nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow), squeeze=False)
        for ax, c in zip(axes.ravel(), order):
            e = np.asarray(X[:, c].todense()).ravel()
            fid = var[c].split(":", 1)[1] if ":" in var[c] else var[c]
            sc = plot_feature(ax, U, e, modality, mode=mode,
                              title=f"{fid}\n{int(e.sum()):,} reads, {int((e>0).sum())} cells",
                              xlim=xlim, ylim=ylim)
            if mode == "expression" and sc is not None:
                cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cb.set_label("q99", fontsize=7)
        for ax in axes.ravel()[n:]:
            ax.axis("off")
        fig.suptitle(f"{label}  {sym} ({gene}) top-{n} {modality}s ({mode})", fontsize=13, y=1.0)
        out = os.path.join(outdir, f"{modality.upper()}_{label}_{sym}_{gene}_top{n}.pdf")

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"[{modality}] {sym} -> {os.path.basename(out)}")
    return out


# --------------------------------------------------------------- module API ----
MODALITY_H5AD_SUFFIX = {
    "gene": "-gene.h5ad", "isoform": "-isoform.h5ad",
    "junction": "-junction.h5ad", "variant": "-variant.h5ad",
}
