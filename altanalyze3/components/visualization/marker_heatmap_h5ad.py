#!/usr/bin/env python3
"""Marker heatmap and MarkerFinder marker table.

This IS the module the cellHarmony-web harness uses for MarkerFinder. It calls
`generate_marker_heatmap_from_adata(..., marker_method="markerfinder")` at
`altanalyze3/components/cellHarmony/flask/pipeline.py:1580`.
`altanalyze3.components.cellHarmony.markerFinder` is imported nowhere in flask/ or webapp/.

Two outputs, two different scopes:

  `*_markers.tsv`   EVERY gene assigned to a population at rho > 0 and FDR <= 0.05.
                    Not bounded by --top-n.
  the heatmap       --top-n markers per cluster, because a figure with one row per
                    assigned gene is unreadable.

Before 2026-09-25 these were one object, so --top-n silently truncated the marker table too.

Harness parameters, for matching a web run from the CLI:
  RNA                 --top-n 50 --cells-per-cluster 100 --export-networks
                      --network-top-n 1000 --skip-expression-tsv
  imputed modalities  --top-n 5 --no-scaling-check --centroid-method mean
  (pipeline.py:1580 and pipeline.py:598)
"""
import argparse
import gzip
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys
import time
import warnings

import matplotlib
import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.errors import PerformanceWarning
import scanpy as sc
from scipy import sparse
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, NoNorm, TwoSlopeNorm, to_rgb, to_rgba
from matplotlib.patches import Rectangle
from altanalyze3.components.visualization import NetPerspective

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
# Compression 0 wrote every rasterized image as raw RGB: no /Filter, no predictor, so a
# 2925 x 3619 heatmap cost 31.8 MB of the 33 MB file. Any non-zero level makes matplotlib
# emit FlateDecode with PNG Predictor 10, the same encoding that keeps the sibling SVG near
# 1 MB. Text stays text because pdf.fonttype 42 is set above and is independent of this.
plt.rcParams['pdf.compression'] = int(os.environ.get("MARKER_HEATMAP_PDF_COMPRESSION", "9"))
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

warnings.filterwarnings(
    'ignore',
    category=PerformanceWarning,
    message='DataFrame is highly fragmented.*',
)

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


def _progress(iterable, **kwargs):
    if tqdm is None:
        return iterable
    return tqdm(iterable, **kwargs)

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def _log_step_timing(label, started_at, timings=None, key=None):
    elapsed = max(0.0, float(time.perf_counter() - started_at))
    print(f"[timing] {label}={elapsed:.2f}s")
    if timings is not None:
        timings[key or label] = round(elapsed, 2)
    return elapsed


def _write_heatmap_cache(cache_path, heatmap_df, row_clusters, column_clusters, ordered_cells,
                         source_PSI=None, imputed_mask=None):
    """Save the plotted fold matrix as an h5ad.

    Earlier versions wrote `np.savez_compressed`. An h5ad holds the same numbers plus named
    obs and var frames, so scanpy, anndata and every other standard reader open the cache
    directly. `_read_heatmap_cache` still reads a `.npz` written before this change.
    """
    matrix = np.asarray(heatmap_df.to_numpy(), dtype=np.float32)
    genes = [str(gene) for gene in heatmap_df.index.tolist()]
    barcodes = [str(barcode) for barcode in ordered_cells]
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(
            f"Heatmap cache matrix is {matrix.shape} but names {len(genes)} genes and "
            f"{len(barcodes)} cells."
        )
    obs = pd.DataFrame(
        {
            "cluster": [str(cluster) for cluster in column_clusters],
            "col_id": [f"{c}:{b}" for c, b in zip(column_clusters, barcodes)],
        },
        index=pd.Index(barcodes, name="barcode"),
    )
    var = pd.DataFrame(
        {
            "cluster": [str(cluster) for cluster in row_clusters],
            "gene": genes,
        },
        index=pd.Index([f"{c}:{g}" for c, g in zip(row_clusters, genes)], name="row_id"),
    )
    cache = AnnData(X=np.ascontiguousarray(matrix.T), obs=obs, var=var)
    cache.uns["marker_heatmap_cache"] = {
        "format": "marker_heatmap_fold_matrix",
        "version": 1,
        "orientation": "cells x markers; transpose for the plotted matrix",
    }
    if source_PSI is not None:
        if not source_PSI.index.equals(heatmap_df.index) or not source_PSI.columns.equals(heatmap_df.columns):
            raise ValueError('Original PSI cache must match the complete plotted identity order')
        cache.layers['source_PSI'] = np.asarray(source_PSI, dtype=np.float32).T.copy()
    if imputed_mask is not None:
        if source_PSI is None or not imputed_mask.index.equals(heatmap_df.index) or not imputed_mask.columns.equals(heatmap_df.columns):
            raise ValueError('Imputation mask requires the aligned original PSI matrix')
        cache.layers['imputed_PSI'] = np.asarray(imputed_mask, dtype=bool).T.copy()
    try:
        cache.write_h5ad(cache_path, compression="lzf")
    except Exception:
        cache.write_h5ad(cache_path, compression="gzip")


def _read_heatmap_cache(cache_path):
    """Rebuild the plotting inputs that `_write_heatmap_cache` saved.

    The cache stores the finished fold matrix, so a completed MarkerFinder result re-renders in
    seconds. This function recomputes nothing. The matrix, the gene order, the cell order and both
    cluster vectors come off disk exactly as the original run wrote them, which is what makes a
    re-render a faithful copy of that run rather than a second analysis of the same data.
    """
    if str(cache_path).lower().endswith(".npz"):
        # Caches written before the h5ad format landed. Kept so a finished run still re-renders.
        with np.load(cache_path, allow_pickle=False) as payload:
            for key in ("matrix", "row_ids", "col_ids", "col_barcodes"):
                if key not in payload:
                    raise ValueError(f"Heatmap cache {cache_path} is missing '{key}'.")
            matrix = np.asarray(payload["matrix"], dtype=np.float32)
            row_ids = [str(x) for x in payload["row_ids"]]
            col_ids = [str(x) for x in payload["col_ids"]]
            col_barcodes = [str(x) for x in payload["col_barcodes"]]
        row_clusters = [entry.partition(":")[0] for entry in row_ids]
        genes = [entry.partition(":")[2] for entry in row_ids]
        column_clusters = [entry.partition(":")[0] for entry in col_ids]
    else:
        cache = sc.read_h5ad(cache_path)
        for frame, column in (("var", "gene"), ("var", "cluster"), ("obs", "cluster")):
            if column not in getattr(cache, frame).columns:
                raise ValueError(f"Heatmap cache {cache_path} is missing {frame}['{column}'].")
        stored = cache.X
        stored = stored.toarray() if sparse.issparse(stored) else np.asarray(stored)
        matrix = np.asarray(stored, dtype=np.float32).T
        genes = [str(gene) for gene in cache.var["gene"].tolist()]
        row_clusters = [str(cluster) for cluster in cache.var["cluster"].tolist()]
        col_barcodes = [str(barcode) for barcode in cache.obs_names.tolist()]
        column_clusters = [str(cluster) for cluster in cache.obs["cluster"].tolist()]

    if matrix.shape != (len(genes), len(col_barcodes)):
        raise ValueError(
            f"Heatmap cache {cache_path} is inconsistent: matrix is {matrix.shape} but it names "
            f"{len(genes)} rows and {len(col_barcodes)} columns."
        )
    if len(column_clusters) != len(col_barcodes):
        raise ValueError(
            f"Heatmap cache {cache_path} names {len(column_clusters)} column clusters for "
            f"{len(col_barcodes)} columns."
        )

    heatmap_df = pd.DataFrame(matrix, index=genes, columns=col_barcodes)

    # The cache is written in plotted order, so each cluster occupies one contiguous block.
    # Counting consecutive runs reproduces the exact block boundaries the original figure drew.
    cluster_counts = []
    for cluster in column_clusters:
        if cluster_counts and cluster_counts[-1][0] == cluster:
            cluster_counts[-1][1] += 1
        else:
            cluster_counts.append([cluster, 1])
    cluster_counts = [(cluster, count) for cluster, count in cluster_counts]
    cluster_order = list(dict.fromkeys(column_clusters))

    return heatmap_df, row_clusters, column_clusters, cluster_order, cluster_counts


def _goelite_terms_for_cache(cache_path, run_dir=None, max_terms=30):
    """Load the GO-Elite term labels that belong to a saved heatmap cache.

    ICGS3 writes the cache to `<run>/MarkerFinder/` and the enrichment table to `<run>/GO-Elite/`,
    so the run directory is the cache's grandparent. This reuses ICGS3's own loader, which means
    the terms, their ordering and the per-cluster cap match what a full run draws.
    """
    if run_dir is None:
        run_dir = os.path.dirname(os.path.dirname(os.path.abspath(cache_path)))
    enrichment = os.path.join(run_dir, "GO-Elite", "icgs3_biomarker_enrichment.tsv")
    if not os.path.exists(enrichment):
        print(f"[INFO] No GO-Elite enrichment found at {enrichment}; drawing without term labels.")
        return {}
    try:
        from altanalyze3.components.clustering.ICGS import ICGS3Config, _load_goelite_heatmap_terms
    except Exception as exc:
        print(f"[WARN] Could not import the ICGS3 GO-Elite loader: {exc}")
        return {}
    config = ICGS3Config(input_paths=[], output_dir=run_dir)
    config.heatmap_goelite_terms = True
    config.heatmap_goelite_max_terms = max_terms
    return _load_goelite_heatmap_terms(run_dir, config)


def _covariates_for_cache(h5ad_path, ordered_cells, covariate_columns):
    """Pull covariate bars for a cache re-render out of the source h5ad.

    The cache stores the plotted cell order but not the cell metadata, so covariate bars need the
    h5ad the run started from. Cells are taken in cache order, so the bars line up with the columns.
    """
    if not h5ad_path or not covariate_columns:
        return None
    if not os.path.exists(h5ad_path):
        print(f"[WARN] Covariate source h5ad not found: {h5ad_path}")
        return None
    import anndata as ad_mod

    adata = ad_mod.read_h5ad(h5ad_path, backed="r")
    try:
        obs = adata.obs.copy()
    finally:
        adata.file.close()
    wanted = [str(c).strip() for c in covariate_columns if str(c).strip()]
    wanted = [c for c in dict.fromkeys(wanted) if c in obs]
    if not wanted:
        print("[INFO] None of the requested covariates are present in the source h5ad obs.")
        return None
    missing = [c for c in ordered_cells if c not in obs.index]
    if missing:
        print(f"[WARN] {len(missing)} of {len(ordered_cells)} plotted cells are absent from "
              f"{h5ad_path}; skipping covariate bars.")
        return None
    retained, skipped = _filter_covariate_columns_for_heatmap(obs, ordered_cells, wanted,
                                                             max_categories=12)
    if skipped:
        skipped_text = ", ".join(f"{name} ({count})" for name, count in skipped)
        print(f"[INFO] Skipping heatmap covariates with >12 displayed categories: {skipped_text}")
    if not retained:
        return None
    print(f"[INFO] Rendering heatmap covariate bars: {', '.join(retained)}")
    return obs.loc[ordered_cells, retained].copy()


def render_heatmap_from_cache(cache_path, output_path, covariate_df=None, go_terms=None,
                              go_terms_max=30, run_dir=None, covariate_h5ad=None,
                              covariate_columns=None, load_go_terms=True,
                              dpi=None, write_svg=True, label_events=None, species=None,
                              label_shared_exons=True, colorbar_label="Norm Exp (z-score)",
                              label_replacements=None, label_shared_introns=False,
                              callout_color='#000000', sample_labels=None,
                              show_pattern_labels=False, caption=None):
    """Redraw a saved heatmap cache with the current plotting code.

    Use this to change how a finished result is drawn -- resolution, colour range, labels --
    without rerunning clustering or MarkerFinder. GO-Elite term labels load automatically from the
    run directory beside the cache. Returns the path written.
    """
    if not os.path.exists(cache_path):
        raise FileNotFoundError(f"Heatmap cache not found: {cache_path}")
    (heatmap_df, row_clusters, column_clusters,
     cluster_order, cluster_counts) = _read_heatmap_cache(cache_path)
    cluster_order = list(dict.fromkeys([*cluster_order, *row_clusters]))
    print(f"Loaded heatmap cache: {heatmap_df.shape[0]} markers x {heatmap_df.shape[1]} cells, "
          f"{len(cluster_order)} clusters")

    if go_terms is None and load_go_terms:
        go_terms = _goelite_terms_for_cache(cache_path, run_dir=run_dir, max_terms=go_terms_max)
    if go_terms:
        print(f"[INFO] GO-Elite term labels: {sum(len(v) for v in go_terms.values())} terms across "
              f"{len(go_terms)} clusters")
    if covariate_df is None:
        covariate_df = _covariates_for_cache(covariate_h5ad, list(heatmap_df.columns),
                                             covariate_columns)

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    _plot_heatmap(
        heatmap_df,
        output_path,
        cluster_counts,
        cluster_order,
        column_clusters,
        row_clusters,
        covariate_df=covariate_df,
        go_terms=go_terms,
        go_terms_max=go_terms_max,
        dpi=dpi,
        write_svg=write_svg,
        label_events=label_events,
        species=species,
        label_shared_exons=label_shared_exons,
        colorbar_label=colorbar_label,
        label_replacements=label_replacements,
        label_shared_introns=label_shared_introns,
        callout_color=callout_color, sample_labels=sample_labels,
        show_pattern_labels=show_pattern_labels, caption=caption,
    )
    return output_path



_ENSEMBL_GENE = re.compile(r"ENS[A-Z]*G\d+")


def _event_uid(label):
    """The splicing-event part of a row label: drops a leading ``cluster:`` prefix before the Ensembl ID."""
    label = str(label)
    match = re.search(r"ENS[A-Z]*G\d+:", label)
    return label[match.start():] if match and match.start() > 0 else label


def _is_event_label(label):
    uid = _event_uid(label)
    return "|" in uid and bool(_ENSEMBL_GENE.match(uid))


def _event_key(label):
    """Orientation-free identity of an event: the set of its junctions (A|B and B|A are one event)."""
    return frozenset(_event_uid(label).split("|"))


def _load_gene_symbols(species=None, annotation_file=None):
    """Ensembl gene ID to symbol, from the bundled ``<Hs|Mm>_Ensembl-annotations.txt.gz``."""
    path = annotation_file
    if path is None:
        prefix = {"human": "Hs", "hs": "Hs", "homo sapiens": "Hs", "mouse": "Mm", "mm": "Mm",
                  "mus musculus": "Mm"}.get(str(species or "human").strip().lower())
        if prefix is None:
            print(f"[WARN] No gene annotation for species {species!r}; event labels keep Ensembl IDs.")
            return {}
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "long_read", "resources", f"{prefix}_Ensembl-annotations.txt.gz")
    if not os.path.exists(path):
        print(f"[WARN] Gene annotation file not found: {path}; event labels keep Ensembl IDs.")
        return {}
    symbols = {}
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 1 and fields[0].startswith("ENS") and fields[1]:
                symbols[fields[0]] = fields[1]
    return symbols


def _format_event_label(label, symbols):
    """``ENSG00000154723:E1.7-E3.1|ENSG00000154723:E1.11-E3.1`` -> ``ATP5J:E1.7-E3.1|E1.11-E3.1``.

    The gene ID becomes its symbol; a later junction drops the repeated gene. A junction into another
    gene keeps that gene, as its symbol. An ID without a symbol stays as the ID.
    """
    parts = _event_uid(label).split("|")
    first_gene = parts[0].split(":", 1)[0]
    out = []
    for index, part in enumerate(parts):
        if index > 0 and part.startswith(first_gene + ":"):
            part = part[len(first_gene) + 1:]
        out.append(_ENSEMBL_GENE.sub(lambda m: symbols.get(m.group(0), m.group(0)), part))
    return "|".join(out)


def _read_label_events(path):
    """Event IDs from a text or TSV file: the ``UID`` column when present, else the first column."""
    table = pd.read_csv(path, sep="\t", dtype=str, header=None)
    first = [str(v).strip() for v in table.iloc[0].tolist()]
    if "UID" in first:
        values = table.iloc[1:, first.index("UID")]
    else:
        values = table.iloc[:, 0]
    events = [str(v).strip() for v in values if isinstance(v, str) and "|" in v]
    print(f"[INFO] Row-label whitelist: {len(events)} events read from {path}")
    return events


def _draw_row_callouts(ax, labels, rows, n_rows, fontsize=6.0, color='#000000'):
    """Label only ``rows``, to the right of the heatmap, with 2-point leader lines spread to avoid overlap."""
    if not rows:
        return
    fig = ax.figure
    fig.canvas.draw()
    height_pt = ax.get_window_extent().height * 72.0 / fig.dpi
    fontsize = min(fontsize, height_pt / (1.2 * max(1, len(rows))))
    original_ylim = ax.get_ylim()
    sep = 1.2 * fontsize * n_rows / max(height_pt, 1.0)
    order = sorted(range(len(rows)), key=lambda i: rows[i])
    target = [float(rows[i]) for i in order]
    placed = []
    for value in target:
        placed.append(max(value, placed[-1] + sep) if placed else max(value, -0.5 + sep / 2.0))
    limit = n_rows - 0.5 - sep / 2.0
    for i in range(len(placed) - 1, -1, -1):
        upper = limit if i == len(placed) - 1 else placed[i + 1] - sep
        placed[i] = min(placed[i], upper)
    transform = ax.get_yaxis_transform()
    artists = []
    for i, y in zip(order, placed):
        artists += ax.plot([1.0, 1.03], [rows[i], y], color=color, linewidth=0.4, transform=transform,
                           clip_on=False)
        artists.append(ax.text(1.035, y, labels[i], transform=transform, ha="left", va="center",
                               fontsize=fontsize, color=color, clip_on=False))
    ax.set_ylim(original_ylim)
    # The page and the heatmap keep their size; a label that would run off the page shrinks the font instead.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    extents = [a.get_window_extent(renderer) for a in artists if hasattr(a, "get_text")]
    left = min(e.x0 for e in extents)
    # Agg hinting can underestimate long labels in the exported PDF. Measure with
    # the standard PDF backend as well, so vector output keeps the complete UID.
    from matplotlib.backends.backend_pdf import RendererPdf
    pdf_renderer = RendererPdf(None, fig.dpi, fig.get_figheight(), fig.get_figwidth())
    pdf_right = [a.get_transform().transform(a.get_position())[0] +
        pdf_renderer.get_text_width_height_descent(a.get_text(), a.get_fontproperties(), False)[0]
        * fig.dpi / 72.0 for a in artists if hasattr(a, 'get_text')]
    right = max([e.x1 for e in extents] + pdf_right)
    overflow = right - (fig.bbox.x1 - 4.0 * fig.dpi / 72.0)
    if overflow > 0 and fontsize > 3.5:
        for artist in artists:
            artist.remove()
        width = right - left
        _draw_row_callouts(ax, labels, rows, n_rows,
                           fontsize=max(3.5, .99 * fontsize * (width - overflow) / width), color=color)


# Splicing PSI defaults (Nathan, 2026-10-11): the settings of the KINNEX-5 healthy PSI run he approved, "appropriate for
# PSI analysis in general". The CLI applies them with --psi; any flag given explicitly overrides its preset value.
# API: generate_marker_heatmap_from_adata(adata, cluster_key=..., out=..., **PSI_DEFAULTS).
PSI_DEFAULTS = {
    "marker_method": "markerfinder",
    "validate_scaling": False,
    "centroid_method": "mean",
    "max_missing_fraction": 0.55,
    "min_rho": 0.25,
    "pval_threshold": 0.01,
    "collapse_shared_exons": True,
    "top_n": 100,
    "max_sample_missing_above_group_median": 0.10,
}


def _drop_high_missing_samples(adata, cluster_key, use_raw, layer, threshold, removed_tsv):
    """Remove samples whose missing-value fraction exceeds their group's median by more than ``threshold``.

    For PSI pseudobulks the missing fraction tracks junction-read depth, so a sample far above its group's median
    is a low-coverage outlier. Writes every sample with its fraction, group median and decision to ``removed_tsv``.
    """
    matrix, _ = _get_expression_matrix(adata, use_raw, layer)
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.issubdtype(values.dtype, np.floating):
        print("[INFO] Sample missing-value filter: the matrix holds no NaN; no sample removed.")
        return adata
    fraction = np.isnan(values).mean(axis=1)
    groups = adata.obs[cluster_key].astype(str).to_numpy()
    group_median = pd.Series(fraction).groupby(groups).transform("median").to_numpy()
    drop = fraction > group_median + threshold
    pd.DataFrame({
        "sample": adata.obs_names.astype(str), cluster_key: groups, "missing_fraction": fraction,
        "group_median": group_median, "removed": drop,
    }).to_csv(removed_tsv, sep="\t", index=False)
    kept = pd.Series(groups[~drop]).value_counts()
    print(f"[INFO] Sample missing-value filter: removed {int(drop.sum())} of {len(drop)} samples missing more than "
          f"their group median + {threshold:g} ({removed_tsv}); smallest group now {int(kept.min()) if len(kept) else 0} "
          f"samples; groups emptied: {sorted(set(groups) - set(kept.index))}.")
    return adata[~drop].copy()


def YellowBlackSky():
    cdict = {
        "red": [(0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)],
        "green": [(0.0, 0.0, 0.8), (0.5, 0.1, 0.0), (1.0, 1.0, 1.0)],
        "blue": [(0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0)],
    }
    return LinearSegmentedColormap("YellowBlackSky", cdict)


def _coerce_lineage_order(raw_order):
    if raw_order is None:
        return None
    try:
        if isinstance(raw_order, dict):
            raw_order = list(raw_order.values())
        if hasattr(raw_order, "tolist"):
            raw_order = raw_order.tolist()
        if isinstance(raw_order, str):
            return None
        return [str(x) for x in raw_order]
    except Exception:
        return None


def _resolve_cluster_order(cluster_series, lineage_order):
    cluster_series = cluster_series.dropna()
    unique_clusters = [str(c) for c in pd.unique(cluster_series)]
    unique_set = set(unique_clusters)

    if lineage_order:
        ordered = [c for c in lineage_order if c in unique_set]
        missing = [c for c in unique_clusters if c not in ordered]
        if ordered:
            return ordered + missing

    if pd.api.types.is_categorical_dtype(cluster_series):
        cats = [str(c) for c in cluster_series.cat.categories if str(c) in unique_set]
        if cats:
            return cats

    return sorted(unique_clusters)


def downsample_cells_per_group(adata, groupby, cells_per_cluster=50, seed=0, group_order=None):
    if cells_per_cluster is None or cells_per_cluster <= 0:
        return adata.copy()

    rng = np.random.default_rng(seed)
    cluster_series = adata.obs[groupby].astype(str)
    if group_order is None:
        group_order = cluster_series.value_counts().index.tolist()

    idx = []
    for group in _progress(group_order, desc="Downsampling clusters"):
        cells = adata.obs_names[cluster_series == str(group)]
        if len(cells) == 0:
            continue
        selected = rng.choice(cells, min(len(cells), cells_per_cluster), replace=False)
        idx.extend(selected.tolist())

    return adata[idx].copy()


def _select_unique_markers(pvals_df, cluster_order, top_n, effect_df=None, pval_threshold=0.05, min_rho=None,
                           ranking="fdr"):
    if ranking not in {"fdr", "rho"}:
        raise ValueError("Marker ranking must be fdr or rho")
    if ranking == "rho":
        if effect_df is None or pval_threshold is not None:
            raise ValueError("Rho-only ranking requires effects and no P/FDR selection threshold")
        effects = effect_df.reindex(columns=cluster_order).replace([np.inf, -np.inf], np.nan)
        effects = effects.where(effects > 0)
        if min_rho is not None:
            effects = effects.where(effects >= min_rho)
        effects = effects.loc[effects.notna().any(axis=1)]
        if effects.empty:
            return pd.DataFrame(columns=["gene", "cluster", "pval", "effect"])
        winner = effects.idxmax(axis=1)
        markers = pd.DataFrame({"gene": effects.index.astype(str), "cluster": winner.to_numpy()})
        markers["effect"] = [effects.loc[g, c] for g, c in zip(effects.index, winner)]
        markers["pval"] = [pvals_df.loc[g, c] for g, c in zip(effects.index, winner)]
        markers["cluster"] = pd.Categorical(markers["cluster"], categories=cluster_order, ordered=True)
        markers = markers.sort_values(["cluster", "effect", "gene"],
                                      ascending=[True, False, True], kind="mergesort")
        if top_n is not None:
            markers = markers.groupby("cluster", sort=False, observed=True).head(top_n)
        return markers.reset_index(drop=True)
    if pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval"])

    pvals_df = pvals_df.reindex(columns=cluster_order)
    pvals_df = pvals_df.replace([np.inf, -np.inf], np.nan)
    mask = ~pvals_df.index.to_series().str.contains('rik', case=False, na=False)
    pvals_df = pvals_df.loc[mask]

    if effect_df is not None:
        effect_df = effect_df.reindex(index=pvals_df.index, columns=cluster_order)
        effect_df = effect_df.replace([np.inf, -np.inf], np.nan)
        pvals_df = pvals_df.where(effect_df > 0)
        if min_rho is not None:
            # Minimum MarkerFinder score: a feature counts for a population only at rho >= min_rho.
            pvals_df = pvals_df.where(effect_df >= min_rho)

    if pval_threshold is not None:
        pvals_df = pvals_df.where(pvals_df <= pval_threshold)

    pvals_df = pvals_df.loc[pvals_df.notna().any(axis=1)]
    if pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval"])

    pvals_df = pvals_df.fillna(1.0)
    min_cluster = pvals_df.idxmin(axis=1)
    min_pval = pvals_df.min(axis=1)

    markers = pd.DataFrame(
        {
            "gene": pvals_df.index.astype(str),
            "cluster": min_cluster.values,
            "pval": min_pval.values,
        }
    )
    markers = markers.dropna(subset=["cluster"])
    markers = markers[markers["cluster"].isin(cluster_order)]
    markers["cluster"] = pd.Categorical(markers["cluster"], categories=cluster_order, ordered=True)

    if effect_df is not None:
        effect_df = effect_df.reindex(index=markers["gene"], columns=cluster_order)
        markers["effect"] = [
            effect_df.loc[g, c] if g in effect_df.index and c in effect_df.columns else np.nan
            for g, c in zip(markers["gene"], markers["cluster"])
        ]
        markers["effect"] = markers["effect"].fillna(-np.inf)
        markers = markers.sort_values(
            ["cluster", "pval", "effect", "gene"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
    else:
        markers = markers.sort_values(
            ["cluster", "pval", "gene"],
            ascending=[True, True, True],
            kind="mergesort",
        )

    if top_n is None:
        # Unique markers are NOT capped. Every gene assignable to a population at rho > 0 and
        # FDR <= pval_threshold is assigned (Nathan's directive, 2026-09-25). Capping here
        # silently discarded genes that met the threshold and wrote no record of them.
        # The heatmap takes its own top_n subset at the call site.
        return markers.reset_index(drop=True)
    selected = markers.groupby("cluster", sort=False).head(top_n)
    return selected.reset_index(drop=True)


def _event_exons(feature, include_introns=False):
    """(gene, exon block) pairs of an AltAnalyze splicing-event ID such as
    ENSG...:E4.5-E5.2|ENSG...:E4.5-E5.1. E4.5 and E4.2 both belong to exon E4. A suffix such as
    _81510496 marks a novel splice site inside the region. Intron regions (I) are
    included only when explicitly requested for exon/intron signature callouts."""
    exons = set()
    for junction in str(feature).split("|"):
        if ":" not in junction:
            continue
        gene, regions = junction.split(":", 1)
        for part in regions.split("-"):
            part_gene = gene
            if ":" in part:
                part_gene, part = part.split(":", 1)
            region = part.split("_", 1)[0]
            if region.startswith("E") or (include_introns and region.startswith("I")):
                exons.add((part_gene, region.split(".", 1)[0]))
    return exons


def _collapse_shared_exon_events(markers):
    """Keep only the most significant of the splicing events that share any exon.

    Events are ranked over all populations by FDR, then by rho (the hook's own order). An event is
    excluded when one of its exons belongs to a more significant event already kept. Returns the
    kept markers in their original order and a table of every excluded event with the kept event
    that shares its exon.
    """
    ranked = markers.assign(rank_order=np.arange(len(markers)))
    sort_cols, ascending = ["pval"], [True]
    if "effect" in ranked.columns:
        sort_cols.append("effect"); ascending.append(False)
    ranked = ranked.sort_values(sort_cols + ["gene"], ascending=ascending + [True], kind="mergesort")
    owner, keep, excluded = {}, [], []
    for row in ranked.itertuples(index=False):
        exons = _event_exons(row.gene)
        shared = sorted(e for e in exons if e in owner)
        if shared:
            excluded.append({"gene": row.gene, "cluster": str(row.cluster), "pval": row.pval,
                             "rho": getattr(row, "effect", np.nan), "kept_event": owner[shared[0]],
                             "shared_exons": ";".join(f"{g}:{e}" for g, e in shared)})
            continue
        keep.append(row.rank_order)
        for exon in exons:
            owner[exon] = row.gene
    kept = markers.iloc[sorted(keep)].reset_index(drop=True)
    return kept, pd.DataFrame(excluded, columns=["gene", "cluster", "pval", "rho", "kept_event", "shared_exons"])


def _select_top_markers_per_cluster(pvals_df, cluster_order, top_n, effect_df=None):
    if pvals_df is None or pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval", "effect"])

    pvals_df = pvals_df.reindex(columns=cluster_order)
    rows = []
    for cluster in cluster_order:
        if cluster not in pvals_df.columns:
            continue
        cluster_pvals = pd.to_numeric(pvals_df[cluster], errors="coerce").fillna(1.0)
        cluster_df = pd.DataFrame(
            {
                "gene": pvals_df.index.astype(str),
                "cluster": str(cluster),
                "pval": cluster_pvals.values,
            }
        )
        if effect_df is not None and cluster in effect_df.columns:
            cluster_effect = pd.to_numeric(effect_df[cluster], errors="coerce")
            cluster_df["effect"] = cluster_effect.reindex(pvals_df.index).values
        else:
            cluster_df["effect"] = np.nan
        cluster_df = cluster_df.sort_values(
            ["pval", "effect", "gene"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        rows.append(cluster_df.head(int(top_n)))
    if not rows:
        return pd.DataFrame(columns=["gene", "cluster", "pval", "effect"])
    return pd.concat(rows, ignore_index=True)



def _zscore_rows(df):
    means = df.mean(axis=1)
    stds = df.std(axis=1, ddof=0).replace(0, np.nan)
    scaled = df.sub(means, axis=0).div(stds, axis=0)
    return scaled.fillna(0.0)


def _zscore_columns(df):
    means = df.mean(axis=0)
    stds = df.std(axis=0, ddof=0).replace(0, np.nan)
    scaled = df.sub(means, axis=1).div(stds, axis=1)
    return scaled.fillna(0.0)


def _bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    finite = np.isfinite(p)
    out = np.full_like(p, np.nan)
    valid = p[finite]
    n = valid.size
    if n == 0:
        return out
    order = np.argsort(valid)
    ranks = np.arange(1, n + 1)
    q = valid[order] * n / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[order] = np.clip(q, 0, 1)
    out[finite] = adjusted
    return out


def _get_expression_matrix(adata, use_raw, layer):
    if use_raw:
        if adata.raw is None:
            raise ValueError("Requested use_raw=True but adata.raw is not set.")
        return adata.raw.X, adata.raw.var_names
    if layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        return adata.layers[layer], adata.var_names
    return adata.X, adata.var_names


def _prepare_marker_stats_aggregates(adata, cluster_key, genes, use_raw, layer):
    ordered_genes = pd.Index(pd.Series(genes, dtype=str)).drop_duplicates().tolist()
    if not ordered_genes:
        return None

    matrix, var_names = _get_expression_matrix(adata, use_raw, layer)
    var_names = pd.Index(var_names).astype(str)
    gene_indexer = var_names.get_indexer(ordered_genes)
    valid_gene_mask = gene_indexer >= 0
    if not np.any(valid_gene_mask):
        return None

    selected_genes = pd.Index(ordered_genes)[valid_gene_mask].astype(str).tolist()
    import h5py
    workspace = getattr(adata, '_matrix_workspace', None)
    disk_matrix = isinstance(matrix, h5py.Dataset) or workspace is not None
    selected_columns = gene_indexer[valid_gene_mask]
    if not disk_matrix:
        matrix = matrix[:, selected_columns]

    cluster_obs = adata.obs[cluster_key]
    cluster_series = cluster_obs.astype(str)
    unique_clusters = set(pd.unique(cluster_series))
    if pd.api.types.is_categorical_dtype(cluster_obs):
        cluster_order = [str(c) for c in cluster_obs.cat.categories if str(c) in unique_clusters]
    else:
        cluster_order = [str(c) for c in pd.unique(cluster_series)]
    if not cluster_order:
        return None

    cluster_codes = pd.Categorical(cluster_series, categories=cluster_order, ordered=True).codes
    valid_cells = cluster_codes >= 0
    if not np.any(valid_cells):
        return None

    row_idx = np.nonzero(valid_cells)[0]
    col_idx = cluster_codes[valid_cells]
    indicator = sparse.csr_matrix(
        (np.ones(len(row_idx), dtype=np.float32), (row_idx, col_idx)),
        shape=(len(cluster_series), len(cluster_order)),
    )

    if disk_matrix:
        unique, inverse = np.unique(selected_columns, return_inverse=True)
        sums_values = np.zeros((len(cluster_order), len(selected_columns)), dtype=matrix.dtype)
        total = np.zeros(len(selected_columns), dtype=matrix.dtype)
        for start in range(0, matrix.shape[0], 512):
            block = matrix[start:start + 512, unique][:, inverse]
            if sparse.issparse(block):
                block = block.tocsr()
                row_groups = cluster_codes[start:start + block.shape[0]]
                keep = valid_cells[start:start + block.shape[0]]
                values = block[keep]
                np.add.at(total, values.indices, values.data)
                groups = np.repeat(row_groups[keep].astype(np.intp), np.diff(values.indptr))
                np.add.at(sums_values.ravel(), groups * sums_values.shape[1] + values.indices, values.data)
                if workspace is not None and start % (512 * 32) == 0:
                    workspace.release_pages()
                continue
            keep = valid_cells[start:start + len(block)]
            values = block[keep]
            np.add.at(sums_values, cluster_codes[start:start + len(block)][keep], values)
            if len(values):
                total = np.cumsum(np.vstack([total, values]), axis=0, dtype=matrix.dtype)[-1]
            if workspace is not None and start % (512 * 32) == 0:
                workspace.release_pages()
        total_sum_values = total.astype(float)
    else:
        measured_values = None
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
        if np.issubdtype(values.dtype, np.floating) and np.isnan(values).any():
            # NaN is a missing value (e.g. PSI). Sum the measured values and count them, so a
            # mean divides by the samples that measured the feature, as MarkerFinder scores it.
            dense = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix, dtype=float)
            measured = ~np.isnan(dense)
            filled = np.where(measured, dense, 0.0)
            sums_values = np.asarray(indicator.T.dot(filled), dtype=float)
            total_sum_values = filled[valid_cells, :].sum(axis=0).astype(float)
            measured_values = np.asarray(indicator.T.dot(measured.astype(np.float64)), dtype=float)
            total_measured_values = measured[valid_cells, :].sum(axis=0).astype(float)
        else:
            sums_matrix = indicator.T.dot(matrix)
            sums_values = sums_matrix.toarray() if sparse.issparse(sums_matrix) else np.asarray(sums_matrix)
            total_sum_values = np.asarray(matrix[valid_cells, :].sum(axis=0)).ravel().astype(float)
    counts_values = np.bincount(col_idx, minlength=len(cluster_order)).astype(float)
    total_count = float(np.sum(valid_cells))

    aggregates = {
        "sums": pd.DataFrame(sums_values, index=cluster_order, columns=selected_genes),
        "counts": pd.Series(counts_values, index=cluster_order),
        "total_sum": pd.Series(total_sum_values, index=selected_genes),
        "total_count": total_count,
    }
    if not disk_matrix and measured_values is not None:
        aggregates["measured_counts"] = pd.DataFrame(measured_values, index=cluster_order, columns=selected_genes)
        aggregates["total_measured"] = pd.Series(total_measured_values, index=selected_genes)
    return aggregates


def _compute_marker_stats(adata, cluster_key, markers_df, fdr_df, use_raw, layer, aggregates=None,
                          effect_column=None):
    """Marker statistics for the output TSV.

    ``effect_column`` names the column that carries ``markers_df["effect"]`` into the output.
    The caller passes "rho" only for ``--marker-method markerfinder``, where the effect is the
    point-biserial Pearson r that MarkerFinder ranks on (markerFinder.py:393). The scanpy path
    puts logfoldchanges in the same slot, so it passes None and the column stays absent.
    Adding the column changes no gene, no row and no existing column.
    """
    if markers_df.empty:
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    if aggregates is None:
        aggregates = _prepare_marker_stats_aggregates(
            adata,
            cluster_key,
            markers_df["gene"].astype(str).tolist(),
            use_raw,
            layer,
        )
    if not aggregates:
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    sums = aggregates["sums"]
    counts = aggregates["counts"]
    total_sum = aggregates["total_sum"]
    total_count = float(aggregates["total_count"])

    sums_values = sums.to_numpy(dtype=float)
    row_clusters = pd.Categorical(markers_df["cluster"].astype(str), categories=sums.index).codes
    col_genes = pd.Categorical(markers_df["gene"].astype(str), categories=sums.columns).codes
    valid = (row_clusters >= 0) & (col_genes >= 0)
    if not np.any(valid):
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    query_sum = np.full(len(markers_df), np.nan, dtype=float)
    query_sum[valid] = sums_values[row_clusters[valid], col_genes[valid]]

    marker_clusters = markers_df["cluster"].astype(str)
    marker_genes = markers_df["gene"].astype(str)
    query_count = counts.reindex(marker_clusters).to_numpy(dtype=float)
    total_gene_sum = total_sum.reindex(marker_genes).to_numpy(dtype=float)

    query_count = np.where(np.isfinite(query_count), query_count, 0.0)
    total_gene_sum = np.where(np.isfinite(total_gene_sum), total_gene_sum, np.nan)
    ref_total_count = np.full(len(markers_df), total_count, dtype=float)
    measured_counts = aggregates.get("measured_counts")
    if measured_counts is not None:
        # Missing values (NaN, e.g. PSI): each mean divides by the samples that measured the feature.
        measured_values = measured_counts.to_numpy(dtype=float)
        query_count = np.zeros(len(markers_df), dtype=float)
        query_count[valid] = measured_values[row_clusters[valid], col_genes[valid]]
        ref_total_count = aggregates["total_measured"].reindex(marker_genes).to_numpy(dtype=float)

    ref_sum = total_gene_sum - query_sum
    ref_count = np.maximum(ref_total_count - query_count, 1.0)
    query_mean = query_sum / np.maximum(query_count, 1.0)
    ref_mean = ref_sum / ref_count

    eps = 1e-9
    fold = np.log2((query_mean + eps) / (ref_mean + eps))

    fdr = np.full(len(markers_df), np.nan, dtype=float)
    if fdr_df is not None and not fdr_df.empty:
        fdr_dedup = fdr_df.loc[~fdr_df.index.duplicated(keep="first")]
        fdr_dedup = fdr_dedup.loc[:, ~fdr_dedup.columns.duplicated(keep="first")]
        fdr_row = pd.Categorical(marker_genes, categories=fdr_dedup.index).codes
        fdr_col = pd.Categorical(marker_clusters, categories=fdr_dedup.columns).codes
        fdr_valid = (fdr_row >= 0) & (fdr_col >= 0)
        if np.any(fdr_valid):
            fdr_values = fdr_dedup.to_numpy(dtype=float)
            fdr[fdr_valid] = fdr_values[fdr_row[fdr_valid], fdr_col[fdr_valid]]

    columns = {
        "Gene": marker_genes,
        "Fold": fold,
        "Query Exp": query_mean,
        "Ref Exp": ref_mean,
        "FDR p-value": fdr,
        "cluster": marker_clusters,
    }
    if effect_column and "effect" in markers_df.columns:
        effect = pd.to_numeric(markers_df["effect"], errors="coerce").to_numpy(dtype=float)
        # _select_unique_markers fills a missing effect with -inf; report that as empty, not -inf.
        columns[effect_column] = np.where(np.isfinite(effect), effect, np.nan)
    out = pd.DataFrame(columns)
    out = out.loc[np.isfinite(out["Query Exp"].to_numpy(dtype=float))].reset_index(drop=True)
    return out


def _order_cells_by_cluster(adata, cluster_key, cluster_order):
    ordered_cells = []
    cluster_counts = []
    cluster_series = adata.obs[cluster_key].astype(str)

    for cluster in cluster_order:
        cells = adata.obs_names[cluster_series == str(cluster)]
        if len(cells) == 0:
            continue
        ordered_cells.extend(cells.tolist())
        cluster_counts.append((str(cluster), len(cells)))

    return ordered_cells, cluster_counts


def _build_heatmap(adata, cluster_key, genes, cluster_order, use_raw, layer, include_column_zscore=True):
    if not genes:
        return pd.DataFrame(), pd.DataFrame(), [], []

    matrix, gene_reference = _get_expression_matrix(adata, use_raw, layer)
    gene_reference = pd.Index(gene_reference).astype(str)
    genes_present = [str(g) for g in genes if str(g) in gene_reference]
    if not genes_present:
        return pd.DataFrame(), pd.DataFrame(), [], []

    ordered_cells, cluster_counts = _order_cells_by_cluster(adata, cluster_key, cluster_order)
    if not ordered_cells:
        return pd.DataFrame(), pd.DataFrame(), [], []

    gene_reference = pd.Index(gene_reference).astype(str)
    cell_indexer = adata.obs_names.get_indexer(ordered_cells)
    gene_indexer = gene_reference.get_indexer(genes_present)
    if np.any(cell_indexer < 0) or np.any(gene_indexer < 0):
        raise ValueError("Failed to align selected cells or genes for heatmap construction.")

    submatrix = matrix[cell_indexer, :][:, gene_indexer]
    if sparse.issparse(submatrix):
        values = submatrix.toarray()
    else:
        values = np.asarray(submatrix)

    heatmap_raw_df = pd.DataFrame(values.T, index=genes_present, columns=ordered_cells)
    heatmap_col_df = _zscore_columns(heatmap_raw_df) if include_column_zscore else pd.DataFrame()
    heatmap_row_df = _zscore_rows(heatmap_raw_df)
    return heatmap_row_df, heatmap_col_df, cluster_counts, ordered_cells


def _build_marker_centroids(adata, cluster_key, genes, cluster_order, use_raw, layer):
    if not genes:
        return pd.DataFrame()

    expr_df = sc.get.obs_df(adata, keys=genes, use_raw=use_raw, layer=layer)
    cluster_series = adata.obs[cluster_key].astype(str)
    expr_df[cluster_key] = cluster_series.values

    ordered_clusters = [cluster for cluster in cluster_order if cluster in set(cluster_series)]
    if not ordered_clusters:
        ordered_clusters = [str(cluster) for cluster in pd.unique(cluster_series)]

    centroid_columns = {}
    for cluster in ordered_clusters:
        subset = expr_df.loc[expr_df[cluster_key] == str(cluster), genes]
        if subset.empty:
            continue
        summed = subset.sum(axis=0).astype(float)
        total_counts = float(summed.sum())
        if total_counts <= 0:
            normalized = pd.Series(0.0, index=genes)
        else:
            normalized = np.log2((summed / total_counts) * 10000.0 + 1.0)
        centroid_columns[str(cluster)] = normalized.reindex(genes)

    return pd.DataFrame(centroid_columns, index=genes)


def _build_marker_centroids_from_aggregates(genes, cluster_order, aggregates, method="log2_cp10k"):
    if method not in {"log2_cp10k", "mean"}:
        raise ValueError(f"Unknown centroid method: {method}")
    if not genes or not aggregates:
        return pd.DataFrame()
    sums = aggregates.get("sums")
    if sums is None or sums.empty:
        return pd.DataFrame()

    genes = [str(g) for g in genes if str(g) in sums.columns]
    if not genes:
        return pd.DataFrame()

    ordered_clusters = [str(cluster) for cluster in cluster_order if str(cluster) in sums.index]
    if not ordered_clusters:
        ordered_clusters = [str(cluster) for cluster in sums.index]
    if not ordered_clusters:
        return pd.DataFrame()

    cluster_gene_sums = sums.reindex(index=ordered_clusters, columns=genes, fill_value=0.0).to_numpy(dtype=float)
    if method == "mean":
        # Imputed lipid, ADT and TF values have a model-defined scale. Retain
        # their signed means; normalizing them as library counts changes it.
        counts = aggregates["counts"].reindex(ordered_clusters).to_numpy(dtype=float)
        means = cluster_gene_sums / np.maximum(counts[:, None], 1.0)
        if "measured_counts" in aggregates:
            # Missing values (NaN, e.g. PSI): mean of the measured samples; NaN where none measured.
            measured = aggregates["measured_counts"].reindex(
                index=ordered_clusters, columns=genes, fill_value=0.0).to_numpy(dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                means = np.where(measured > 0, cluster_gene_sums / measured, np.nan)
        return pd.DataFrame(means.T, index=genes, columns=ordered_clusters)
    # The centroid is a pseudobulk CPM: log2(cluster_sum / cluster_total * 10000 + 1). The
    # formula is defined only for non-negative counts. Summing mean-centred or z-scored values
    # gives negative cluster totals, and log2 of the resulting negative ratio is NaN, which
    # would write a silently unusable cellHarmony reference.
    if np.any(cluster_gene_sums < 0):
        # A CITE-seq matrix carries signed AB_ columns beside non-negative RNA, so a few
        # cluster x feature sums go negative while the RNA block is fine. Clamp those to 0,
        # which is the floor of a background-subtracted protein measurement, and say how many.
        # Nathan's directive, 2026-10-01: "make values compatible (if negative values,
        # truncate to 0)". The check below stays in force for anything the clamp cannot fix.
        n_neg = int((cluster_gene_sums < 0).sum())
        frac = n_neg / float(cluster_gene_sums.size)
        print("[centroid] clamped %d of %d cluster x feature sums (%.4f%%) from negative to 0 "
              "before the log2 CPM; a background-subtracted value has 0 as its floor"
              % (n_neg, cluster_gene_sums.size, 100 * frac))
        cluster_gene_sums = np.clip(cluster_gene_sums, 0.0, None)
    if np.any(cluster_gene_sums < 0):
        raise ValueError(
            "Centroid matrix still holds negative cluster x gene sums after clamping; "
            "the log2 CPM is undefined. Point --layer at a counts matrix."
        )
    totals = cluster_gene_sums.sum(axis=1, keepdims=True)
    normalized = np.zeros_like(cluster_gene_sums, dtype=float)
    valid = totals[:, 0] > 0
    if np.any(valid):
        normalized[valid] = np.log2((cluster_gene_sums[valid] / totals[valid]) * 10000.0 + 1.0)
    return pd.DataFrame(normalized.T, index=genes, columns=ordered_clusters)


def _append_suffix_before_extension(path, suffix):
    root, ext = os.path.splitext(path)
    if ext:
        return f"{root}{suffix}{ext}"
    return f"{path}{suffix}"


def _strip_marker_cluster_prefix(value):
    text = str(value or "")
    if ":" not in text:
        return text
    _, suffix = text.split(":", 1)
    return suffix or text


def _cluster_color_map(cluster_order):
    # `palette(i % palette.N)` wrapped every 20 clusters, so cluster 1 and cluster 21 drew in
    # one colour on the heatmap cluster bar. palettes.categorical_colors keeps every colour
    # distinct into the hundreds; see visualization/palettes.py.
    from altanalyze3.components.visualization.palettes import categorical_colors

    hex_colors = categorical_colors(len(cluster_order))
    colors = [to_rgba(c) for c in hex_colors]
    return {cluster: colors[i] for i, cluster in enumerate(cluster_order)}, ListedColormap(colors)


def _axis_label_fontsize(label_count, axis):
    label_count = max(1, int(label_count))
    if axis == "x":
        if label_count <= 6:
            return 10
        if label_count <= 12:
            return 9
        if label_count <= 20:
            return 8
        if label_count <= 32:
            return 7
        if label_count <= 48:
            return 6
        if label_count <= 72:
            return 5
        if label_count <= 96:
            return 4
        return 3

    if label_count <= 20:
        return 9
    if label_count <= 40:
        return 8
    if label_count <= 80:
        return 7
    if label_count <= 140:
        return 6
    if label_count <= 220:
        return 5
    if label_count <= 320:
        return 4
    return 3


def _format_heatmap_pvalue(pvalue):
    try:
        pval = float(pvalue)
    except Exception:
        return "p=NA"
    if pval < 1e-3:
        return "p={:.1e}".format(pval)
    return "p={:.3f}".format(pval)


def _category_color_map(categories, cmap_name="tab20"):
    # cmap_name stays in the signature for callers that pass it, but a named matplotlib map
    # cannot supply more than its own colour count. Build the covariate colours the same way
    # the cluster bar does, so a 25-library covariate bar no longer repeats 5 colours.
    from altanalyze3.components.visualization.palettes import categorical_colors

    cats = [str(cat) for cat in categories]
    hex_colors = categorical_colors(len(cats))
    return {cat: to_rgb(hex_colors[i]) for i, cat in enumerate(cats)}


def _covariate_palette_name(index):
    palettes = [
        "tab20",
        "Set3",
        "Paired",
        "Dark2",
        "Accent",
        "Pastel1",
        "Pastel2",
        "Set2",
        "tab20b",
        "tab20c",
    ]
    return palettes[int(index) % len(palettes)]


def _covariate_categories(values):
    series = pd.Series(values).astype(str)
    if pd.api.types.is_categorical_dtype(values):
        present = set(series)
        return [str(c) for c in values.cat.categories if str(c) in present]
    return [str(c) for c in pd.unique(series)]


def _build_covariate_track(covariate_df):
    """Return integer category codes, one colormap covering them, and the legends.

    This function used to return an (n_covariates, n_cells, 3) RGB float image. Matplotlib's
    Agg resampler mishandles that array. With a single covariate it must upsample 1 input row
    to the full bar height while downsampling the cells to the raster width, and it writes the
    data into output row 0 and floods every other row with one colour. The BPD run of
    2026-08-23 drew author_disease as 100% 'control' across 236,620 cells for that reason:
    scanline 0 of the rendered bar held the correct 43.8/23.3/22.0/10.9 split and the other
    552 of 553 scanlines held flat pink.

    Integer codes plus a ListedColormap take the same path as the ax_top cluster bar, which
    renders 20 distinct colours correctly at the same width in the same figure.
    """
    if covariate_df is None or covariate_df.empty:
        return None, None, {}
    n_cells = int(covariate_df.shape[0])
    n_covariates = int(covariate_df.shape[1])
    codes = np.zeros((n_covariates, n_cells), dtype=int)
    colors = []
    legends = {}
    for idx, col in enumerate(covariate_df.columns):
        values = covariate_df[col].astype(str)
        categories = _covariate_categories(covariate_df[col])
        color_map = _category_color_map(categories, cmap_name=_covariate_palette_name(idx))
        lookup = {}
        for category in categories:
            lookup[str(category)] = len(colors)
            colors.append(color_map[str(category)])
        unknown_code = len(colors)
        colors.append((0.8, 0.8, 0.8))
        codes[idx] = [lookup.get(str(value), unknown_code) for value in values]
        legends[str(col)] = {
            "categories": categories,
            "colors": color_map,
        }
    if codes.shape != (n_covariates, n_cells):
        raise ValueError(
            f"Covariate track is {codes.shape} but the frame holds {n_covariates} covariates "
            f"and {n_cells} cells."
        )
    return codes, ListedColormap(colors), legends


def _draw_covariate_bars(ax, codes, colors, max_spans=1500):
    """Draw each covariate bar as filled rectangles, one per run of one category.

    imshow cannot draw this bar. Matplotlib rasterizes a thin wide image into the PDF and writes
    the data into output row 0 only, then floods every other row with a single value. Measured on
    the BPD run of 2026-08-23: 1 of 553 rendered scanlines held the correct author_disease split
    of 43.8/23.3/22.0/10.9 and the other 552 held flat pink, so the bar read as 100% control while
    the true composition was 43.5% control of 236,620 cells. Row count (1 to 700), interpolation
    ('none' and 'nearest') and column count (4,096 to 236,620) all reproduce it.

    Rectangles never enter matplotlib's image resampler, so the defect cannot occur. They also
    stay vector in the PDF, which keeps the bar editable.

    Columns bin to at most `max_spans` spans and each span takes the category of the cell at its
    centre. The span count is capped at 1,500 so one span stays about 2 printed pixels wide on a
    600 dpi page. At 4,000 spans a span fell to 0.73 px, isolated spans antialiased away, and only
    long runs kept an exact colour: control then measured 55.2% of the exactly-coloured pixels
    against a true 43.5%. Centre sampling is what a correct nearest-neighbour image downsample does, and it keeps
    the drawn shares unbiased: decimating this run's 236,620 cells to 1,400 and 2,000 spans gave
    control at 43.4% and 42.5% against a true 43.5%. Two rejected alternatives, both measured on
    this run: the most frequent category per bin inflated control to 64.9%, and dividing each bin
    among its categories in proportion produced sub-pixel slivers that antialiased into blends, so
    only 719 pixels of the bar held any exact category colour.
    """
    n_rows, n_cells = codes.shape
    n_spans = int(min(int(max_spans), n_cells))
    if n_spans <= 0:
        raise ValueError("Covariate bar needs at least one cell to draw.")
    edges = np.linspace(0, n_cells, n_spans + 1).astype(int)
    centres = np.clip((edges[:-1] + edges[1:]) // 2, 0, n_cells - 1)
    ax.set_xlim(-0.5, n_cells - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    drawn = 0
    for row in range(n_rows):
        spans = np.asarray(codes[row], dtype=int)[centres]
        start = 0
        for index in range(1, n_spans + 1):
            if index == n_spans or spans[index] != spans[start]:
                left = float(edges[start]) - 0.5
                right = float(edges[index]) - 0.5
                ax.add_patch(Rectangle(
                    (left, row - 0.5), max(right - left, 1.0), 1.0,
                    facecolor=colors[int(spans[start])], edgecolor="none", linewidth=0,
                ))
                drawn += 1
                start = index
    return drawn


def _filter_covariate_columns_for_heatmap(obs, ordered_cells, covariate_columns, max_categories=12):
    retained = []
    skipped = []
    for covariate in covariate_columns:
        categories = _covariate_categories(obs.loc[ordered_cells, covariate])
        if len(categories) > int(max_categories):
            skipped.append((covariate, len(categories)))
        else:
            retained.append(covariate)
    return retained, skipped


def _draw_go_term_labels(fig, ax_terms, ax, row_clusters, go_terms, row_color_map, max_terms):
    if ax_terms is None or not go_terms:
        return
    ax_terms.set_xlim(0, 1)
    ax_terms.set_ylim(ax.get_ylim())
    ax_terms.set_xticks([])
    ax_terms.set_yticks([])
    ax_terms.axis("off")
    ax_terms.patch.set_alpha(0.0)

    cluster_ranges = []
    current = None
    start = 0
    for idx, cluster in enumerate(row_clusters):
        cluster = str(cluster)
        if cluster != current:
            if current is not None:
                cluster_ranges.append((current, start, idx))
            current = cluster
            start = idx
    if current is not None:
        cluster_ranges.append((current, start, len(row_clusters)))

    fig.canvas.draw()
    axis_height_px = ax.get_window_extent().height
    row_height_px = axis_height_px / max(1, len(row_clusters))
    term_font_size = 6
    term_spacing_px = max(10.0, term_font_size * fig.dpi / 72.0 * 3.0)
    gene_scale = max(1.0, len(row_clusters) / 600.0)
    min_rows_per_term = max(1.0, term_spacing_px / max(1.0, row_height_px)) * gene_scale

    entries = []
    for cluster, start, end in cluster_ranges:
        terms = go_terms.get(str(cluster), []) or []
        block_height = max(1, end - start)
        slot_count = max(1, int(block_height / max(1.0, min_rows_per_term)))
        limit = min(len(terms), slot_count)
        if max_terms:
            limit = min(limit, int(max_terms))
        entries.append(
            {
                "cluster": str(cluster),
                "start": start,
                "end": end,
                "block_height": block_height,
                "terms": terms[:limit],
                "positions": [],
            }
        )

    def compute_positions(entry):
        terms = entry["terms"]
        if not terms:
            entry["positions"] = []
            return
        block_height = max(1, entry["end"] - entry["start"])
        available_rows = max(1.0, block_height - 1)
        total_height = (len(terms) - 1) * min_rows_per_term if len(terms) > 1 else 0.0
        total_height = min(total_height, available_rows)
        top_offset = max(0.0, (available_rows - total_height) / 2.0)
        start_y = entry["start"] + 0.5 + top_offset
        positions = []
        for i in range(len(terms)):
            y = start_y + i * min_rows_per_term
            if y > entry["end"] - 0.5:
                break
            positions.append(y)
        entry["positions"] = positions

    for entry in entries:
        compute_positions(entry)

    # Resolve collisions between adjacent cluster blocks, but NEVER prune a block to
    # zero terms. The previous loop dropped the last term from whichever block was
    # taller, with `prev` winning every tie. When all blocks are the same height --
    # the normal case, since each cluster contributes the same top_n markers -- the
    # tie always fell on `prev`, so each block in turn was stripped to empty and only
    # the final cluster kept a label. On the 42-cluster Adams run that rendered 1 of
    # 42 GO-Elite labels. Keeping a floor of one term per block preserves the intended
    # one-label-per-cluster minimum; a lone label is centred in its own block by
    # compute_positions, so it stays inside that cluster's rows.
    def _n_terms(entry):
        return len(entry["terms"])

    for idx in range(1, len(entries)):
        prev = entries[idx - 1]
        curr = entries[idx]
        while prev["positions"] and curr["positions"]:
            if (curr["positions"][0] - prev["positions"][-1]) >= min_rows_per_term:
                break
            # Prune from the block that still has the most terms; stop once both are
            # down to their last label.
            if _n_terms(prev) <= 1 and _n_terms(curr) <= 1:
                break
            if _n_terms(prev) >= _n_terms(curr):
                prev["terms"] = prev["terms"][:-1]
                compute_positions(prev)
            else:
                curr["terms"] = curr["terms"][:-1]
                compute_positions(curr)

    for entry in entries:
        color = row_color_map.get(str(entry["cluster"]), "blue")
        for y, item in zip(entry["positions"], entry["terms"]):
            if isinstance(item, (tuple, list)) and len(item) >= 2:
                label = "{} {}".format(str(item[0]), _format_heatmap_pvalue(item[1]))
            else:
                label = str(item)
            ax_terms.text(
                0.995,
                y,
                label,
                transform=ax_terms.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=term_font_size,
                color=color,
                clip_on=False,
            )


def _draw_covariate_legends(fig, covariate_legends, rect=(0.255, 0.078, 0.690, 0.115)):
    if not covariate_legends:
        return
    import matplotlib.patches as mpatches
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    ax_leg = fig.add_axes(rect)
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")
    fig.canvas.draw()
    axis_width_pt = max(1.0, ax_leg.get_window_extent().width * 72.0 / fig.dpi)
    legend_font = FontProperties(family=plt.rcParams.get("font.family", "sans-serif"), size=8)

    def text_width_axes(label):
        if not str(label):
            return 0.0
        return TextPath((0, 0), str(label), prop=legend_font).get_extents().width / axis_width_pt

    max_x = 1.0
    title_x = 0.245
    item_start_x = 0.275
    sw = 0.022
    sh = 0.115
    row_h = 0.190
    row_idx = 0
    for covariate, info in covariate_legends.items():
        y_center = 0.82 - row_idx * row_h
        ax_leg.text(
            title_x,
            y_center,
            str(covariate),
            fontsize=8,
            ha="right",
            va="center",
            weight="bold",
            transform=ax_leg.transAxes,
        )
        x = item_start_x
        categories = list(info["categories"])
        remaining_width = max(0.05, max_x - item_start_x)
        raw_widths = []
        for cat in categories:
            label = str(cat)
            raw_widths.append(sw + 0.010 + text_width_axes(label) + 0.030)
        min_gap = 0.010
        for cat, raw_item_w in zip(categories, raw_widths):
            label = str(cat)
            # Shrinking spacing without shrinking the text made long legends overlap.
            # Wrap complete entries instead; no category is silently omitted.
            item_w = raw_item_w
            if x > item_start_x and x + item_w > max_x:
                row_idx += 1
                y_center = 0.82 - row_idx * row_h
                x = item_start_x
            patch = mpatches.Rectangle(
                (x, y_center - sh / 2.0),
                sw,
                sh,
                transform=ax_leg.transAxes,
                color=info["colors"][label],
                linewidth=0.25,
                clip_on=False,
            )
            ax_leg.add_patch(patch)
            ax_leg.text(x + sw + min_gap, y_center, label, fontsize=8, ha="left", va="center", transform=ax_leg.transAxes)
            x += item_w
        row_idx += 1


def _plot_heatmap(
    heatmap_df,
    output_path,
    cluster_counts,
    cluster_order,
    column_clusters,
    row_clusters,
    covariate_df=None,
    go_terms=None,
    go_terms_max=30,
    dpi=None,
    write_svg=True,
    label_events=None,
    species=None,
    label_shared_exons=True,
    colorbar_label="Norm Exp (z-score)",
    label_replacements=None,
    figure_title=None,
    caption=None,
    sample_labels=None,
    show_pattern_labels=False,
    label_shared_introns=False,
    callout_color='#000000',
):
    if heatmap_df.empty:
        raise ValueError("Heatmap dataframe is empty.")
    if label_replacements:
        def display_label(value):
            value=str(value)
            for before,after in label_replacements.items(): value=value.replace(str(before),str(after))
            return value
        cluster_counts=[(display_label(c),n) for c,n in cluster_counts]
        cluster_order=[display_label(c) for c in cluster_order]
        if len(set(cluster_order)) != len(cluster_order):
            raise ValueError('Display label replacements may not merge distinct source groups')
        column_clusters=[display_label(c) for c in column_clusters]
        row_clusters=[display_label(c) for c in row_clusters]
        if covariate_df is not None:
            covariate_df=covariate_df.copy()
            for col in covariate_df:
                if isinstance(covariate_df[col].dtype,pd.CategoricalDtype):
                    covariate_df[col]=covariate_df[col].cat.rename_categories(display_label)
                else: covariate_df[col]=covariate_df[col].map(display_label)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    show_go_terms = bool(go_terms)
    covariate_track, covariate_cmap, covariate_legends = _build_covariate_track(covariate_df)

    fixed_width = 10.8 if show_go_terms else 8.3

    # Covariate bars and their legends need room under the heatmap. Taking that room out of the
    # heatmap shortened it by 19% (6.12 in -> 4.95 in) and squeezed the GO-Elite term labels until
    # they overlapped. Growing the page instead keeps the heatmap exactly as tall as it is with no
    # covariates, so the term labels read the same either way. These are inches, and the fractions
    # they produce reproduce the previous layout exactly when no covariate bars are drawn.
    axes_top = 0.82
    heatmap_height_in = 6.12 if show_go_terms else 6.29
    below_axes_in = 2.43 if covariate_track is not None else (0.85 if show_go_terms else 0.68)
    # Optional iterative-figure context is opt-in; ordinary MarkerFinder geometry is unchanged.
    import textwrap
    caption_lines = textwrap.wrap(str(caption), width=135) if caption else []
    caption_height_in = (len(caption_lines) * 0.115 + 0.16) if caption_lines else 0.
    sample_height_in = (max(.65, max((len(str(v)) for v in sample_labels.values()), default=0)
                           * 5.5 * .55 / 72 + .12) if sample_labels is not None else 0.)
    footer_height_in = caption_height_in + sample_height_in
    below_axes_in += footer_height_in
    fixed_height = (heatmap_height_in + below_axes_in) / axes_top
    bottom_frac = below_axes_in / fixed_height

    if show_go_terms:
        fig = plt.figure(figsize=(fixed_width, fixed_height))
        gs = fig.add_gridspec(
            1,
            3,
            width_ratios=[0.35, 0.50, 0.15],
            wspace=0.0,
            left=0.06,
            right=0.98,
            top=axes_top,
            bottom=bottom_frac,
        )
        ax_terms = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_terms)
        ax_labels = fig.add_subplot(gs[0, 2], sharey=ax)
    else:
        fig, ax = plt.subplots(figsize=(fixed_width, fixed_height))
        fig.subplots_adjust(left=0.26 if show_pattern_labels else 0.08,
                            right=0.80, top=axes_top, bottom=bottom_frac)
        ax_terms = None
        ax_labels = None

    contrast = float(os.environ.get("MARKER_HEATMAP_CONTRAST", "1.0"))
    if contrast <= 0:
        contrast = 1.0
    # Row-median-centred log2 values almost never reach +/-3, so a +/-3 scale wastes most of
    # the colour range and washes real differences out. +/-1.5 is the default; override with
    # MARKER_HEATMAP_RANGE (a single number, the half-range) or MARKER_HEATMAP_CONTRAST.
    half_range = float(os.environ.get("MARKER_HEATMAP_RANGE", "1.5"))
    if half_range <= 0:
        half_range = 1.5
    vmin, vmax = -half_range / contrast, half_range / contrast
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    heatmap_cmap = YellowBlackSky()
    # Missing PSI is visually black by explicit figure instruction; retain NaN in all data/caches.
    heatmap_cmap.set_bad('black')
    im = ax.imshow(
        heatmap_df.values,
        aspect="auto",
        cmap=heatmap_cmap,
        norm=norm,
        interpolation="none",
    )
    # Same protocol as fastCNV.py:1140. Without it matplotlib embeds one image pixel per
    # matrix cell: a 289,066-cell heatmap wrote a 289066 x 1160 image, 335 million pixels,
    # into a 1,007 MB PDF. set_rasterized re-renders this artist through the raster backend
    # at the save dpi, so the stamp is sized to the page. Text, axes, labels, the colour bar
    # and the legends stay vector.
    im.set_rasterized(_heatmap_rasterize())

    cluster_colors, cluster_cmap = _cluster_color_map(cluster_order)
    cluster_to_id = {cluster: idx for idx, cluster in enumerate(cluster_order)}

    centers = []
    labels = []
    boundaries = []
    pos = 0
    total = sum(count for _, count in cluster_counts)
    for cluster, count in cluster_counts:
        if count <= 0:
            continue
        start = pos
        end = pos + count - 1
        labels.append(cluster)
        centers.append((start + end) / 2.0)
        pos = end + 1
        if pos < total:
            boundaries.append(end + 0.5)

    x_fontsize = _axis_label_fontsize(len(labels), axis="x")
    if centers:
        ax.set_xticks(centers)
        ax.set_xticklabels(labels, rotation=45, ha="left", fontsize=x_fontsize)
    else:
        ax.set_xticks([])

    ax.xaxis.tick_top()
    ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, length=0, pad=7.5)

    gene_fontsize = _axis_label_fontsize(heatmap_df.shape[0], axis="y")
    y_positions = np.arange(heatmap_df.shape[0])
    raw_labels = [str(g) for g in heatmap_df.index]
    # Splicing events read as SYMBOL:junction|junction (Nathan, 2026-10-11); other labels are unchanged.
    symbols = _load_gene_symbols(species) if any(_is_event_label(g) for g in raw_labels) else {}
    y_labels = [_format_event_label(g, symbols) if symbols and _is_event_label(g) else g for g in raw_labels]
    callout_rows = None
    if label_events is not None:
        # A row is labelled when it is a listed event (either junction order) or, unless
        # label_shared_exons is False, when it shares an exon block with one: the redundancy rule of
        # the shared-exon collapse (Nathan, 2026-10-11).
        listed = {_event_key(e): _event_uid(e) for e in label_events}
        listed_exons = {uid: _event_exons(uid, include_introns=label_shared_introns)
                        for uid in listed.values()}
        records = []
        for i, g in enumerate(raw_labels):
            if _event_key(g) in listed:
                records.append((i, "listed event", listed[_event_key(g)], ""))
                continue
            if not label_shared_exons:
                continue
            row_exons = _event_exons(_event_uid(g), include_introns=label_shared_introns)
            for uid, exons in listed_exons.items():
                shared = sorted(row_exons & exons)
                if shared:
                    records.append((i, ("shares an exon/intron with a listed event" if label_shared_introns
                                        else "shares an exon with a listed event"), uid,
                                    ";".join(f"{gene}:{block}" for gene, block in shared)))
        callout_rows = sorted({r[0] for r in records})
        n_exact = len({r[0] for r in records if r[1] == "listed event"})
        print(f"[INFO] Row-label whitelist: {len(callout_rows)} of {len(raw_labels)} heatmap rows labelled for "
              f"{len(listed)} listed events: {n_exact} listed events (either junction order), "
              f"{len(callout_rows) - n_exact} sharing an exon with one")
        label_table = os.path.splitext(output_path)[0] + "_label_events.tsv"
        pd.DataFrame(
            [{"heatmap_row": i + 1, "event": _event_uid(raw_labels[i]), "label": y_labels[i],
              "cluster": row_clusters[i] if i < len(row_clusters) else "", "match": kind,
              "listed_event": uid, "shared_exons": exons} for i, kind, uid, exons in records],
            columns=["heatmap_row", "event", "label", "cluster", "match", "listed_event", "shared_exons"],
        ).to_csv(label_table, sep="\t", index=False)
        print(f"[INFO] Labelled rows and their listed events: {label_table}")
    if ax_labels is None and callout_rows is not None:
        ax.set_yticks([])
        ax.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)
    elif ax_labels is None:
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=gene_fontsize)
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", left=False, right=True, labelleft=False, labelright=True, pad=1)
    else:
        ax.set_yticks([])
        ax.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)

    for boundary in boundaries:
        ax.axvline(boundary, color="white", linewidth=0.5)

    if show_pattern_labels:
        start = 0
        while start < len(row_clusters):
            end = start + 1
            while end < len(row_clusters) and row_clusters[end] == row_clusters[start]:
                end += 1
            label = str(row_clusters[start]).replace('up-', 'up ', 1).replace('down-', 'down ', 1)
            ax.text(-0.045, (start + end - 1) / 2., label,
                    transform=ax.get_yaxis_transform(), ha='right', va='center',
                    fontsize=6.5, clip_on=False)
            if end < len(row_clusters):
                ax.axhline(end - .5, color='white', linewidth=.6)
            start = end

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="1.125%", pad=0.02)
    ax_left = divider.append_axes("left", size="1.5%", pad=0.02)
    col_ids = np.array([cluster_to_id[c] for c in column_clusters], dtype=int)[None, :]
    row_ids = np.array([cluster_to_id[c] for c in row_clusters], dtype=int)[:, None]

    # Categorical integer IDs index the full shared palette directly. Automatic
    # normalization stretched a cohort-only subset across row-pattern colors,
    # contradicting the same cohort's covariate legend in consolidated figures.
    ax_top.imshow(col_ids, aspect="auto", cmap=cluster_cmap, norm=NoNorm(),
                  interpolation="none").set_rasterized(_heatmap_rasterize())
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.set_ylabel("")

    ax_left.imshow(row_ids, aspect="auto", cmap=cluster_cmap, norm=NoNorm(),
                   interpolation="none").set_rasterized(_heatmap_rasterize())
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    ax_left.set_xlabel("")

    if covariate_track is not None:
        cov_rows = int(covariate_track.shape[0])
        ax_cov = divider.append_axes("bottom", size=f"{max(4.0, cov_rows * 2.0)}%", pad=0.040)
        _draw_covariate_bars(ax_cov, covariate_track, list(covariate_cmap.colors))
        ax_cov.set_aspect("auto")
        ax_cov.set_xticks([])
        ax_cov.set_yticks([])
        for idx, name in enumerate(covariate_df.columns):
            ax_cov.text(
                -0.006,
                idx,
                str(name),
                transform=ax_cov.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=7,
                weight="bold",
                clip_on=False,
            )
        ax_cov.tick_params(axis="y", length=0)
        for spine in ax_cov.spines.values():
            spine.set_visible(False)
        ax_cov.set_facecolor("none")
        for boundary in boundaries:
            ax_cov.axvline(boundary, color="white", linewidth=0.35)
        for row_boundary in np.arange(0.5, cov_rows - 0.5 + 1e-9, 1.0):
            ax_cov.axhline(row_boundary, color="white", linewidth=0.35)
        # The default rect is a fraction of a 9.0 in page. Convert to the same absolute inches
        # so a taller page does not stretch the legend block away from the covariate bars.
        _draw_covariate_legends(
            fig,
            covariate_legends,
            rect=(0.255, (0.702 + caption_height_in) / fixed_height, 0.690, 1.035 / fixed_height),
        )
        cbar_anchor_ax = ax_cov
    else:
        cbar_anchor_ax = ax

    if ax_terms is not None and ax_labels is not None:
        fig.canvas.draw()
        heatmap_bbox = ax.get_position()
        ax_terms.set_position([ax_terms.get_position().x0, heatmap_bbox.y0, ax_terms.get_position().width, heatmap_bbox.height])
        ax_labels.set_position([ax_labels.get_position().x0, heatmap_bbox.y0, ax_labels.get_position().width, heatmap_bbox.height])
        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(ax.get_ylim())
        ax_labels.set_xticks([])
        ax_labels.set_yticks([])
        ax_labels.axis("off")
        ax_labels.patch.set_alpha(0.0)

        fig.canvas.draw()
        row_height_px = ax_labels.get_window_extent().height / max(1, heatmap_df.shape[0])
        row_height_pt = row_height_px * 72.0 / fig.dpi
        fitted_gene_fontsize = max(2.5, min(gene_fontsize, row_height_pt * 0.80))
        if callout_rows is not None:
            # The whitelist uses the same spaced leader-line callouts with or without
            # GO terms. Fitting them to every matrix row made sparse labels overlap
            # and shrink on large heatmaps despite the empty rows between callouts.
            _draw_row_callouts(ax, [y_labels[i] for i in callout_rows],
                               callout_rows, heatmap_df.shape[0], color=callout_color)
        else:
            for pos, gene in enumerate(y_labels):
                ax_labels.text(
                    0.02,
                    pos,
                    gene,
                    transform=ax_labels.get_yaxis_transform(),
                    ha="left",
                    va="center",
                    fontsize=fitted_gene_fontsize,
                    fontstyle="italic",
                    clip_on=False,
                )
        _draw_go_term_labels(fig, ax_terms, ax, row_clusters, go_terms, cluster_colors, go_terms_max)

    fig.canvas.draw()
    heatmap_bbox = ax.get_position()
    anchor_bbox = cbar_anchor_ax.get_position()
    cbar_width = heatmap_bbox.width * 0.35
    cbar_height = 0.012
    cbar_x = heatmap_bbox.x0 + (heatmap_bbox.width - cbar_width) / 2.0
    cbar_y = max(0.02, anchor_bbox.y0 - (0.026 if covariate_track is not None else 0.055))
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.set_xlim(vmin, vmax)
    cbar.set_label(colorbar_label, fontsize=8, labelpad=2)
    cbar.set_ticks([])
    # :g not :.0f -- the default half-range is 1.5, which :.0f would print as "-2"/"2"
    cbar.ax.text(-0.08, 0.5, f"{vmin:g}", ha="right", va="center", transform=cbar.ax.transAxes)
    cbar.ax.text(1.08, 0.5, f"{vmax:g}", ha="left", va="center", transform=cbar.ax.transAxes)

    if ax_labels is None and callout_rows is not None:
        _draw_row_callouts(ax, [y_labels[i] for i in callout_rows], callout_rows,
                           heatmap_df.shape[0], color=callout_color)

    if sample_labels is not None:
        labels = [str(sample_labels[c]) for c in heatmap_df.columns]
        # Put donor identities below the covariates; the heatmap retains every ordered column.
        identity_axis = cbar_anchor_ax.secondary_xaxis('bottom')
        identity_axis.set_xticks(np.arange(len(labels)))
        identity_axis.set_xticklabels(labels, rotation=90, ha='center', fontsize=5.5)
        identity_axis.tick_params(length=0, pad=2)
        identity_axis.spines['bottom'].set_visible(False)
        cax.set_position([cbar_x, cbar_y - sample_height_in / fixed_height,
                          cbar_width, cbar_height])
    if figure_title:
        fig.text(.5, .985, str(figure_title), ha='center', va='top', fontsize=12, weight='bold')
    if caption_lines:
        fig.text(.045, .035 / fixed_height, '\n'.join(caption_lines),
                 ha='left', va='bottom', fontsize=7, linespacing=1.15)

    # Iterative review labels can contain long literal multi-gene junction IDs.
    # Preserve their full text by expanding the exported page, without shrinking
    # fonts, clipping identifiers, or changing the heatmap data/axis geometry.
    _save_figure(fig, output_path, dpi=dpi, write_svg=write_svg,
                 bbox_inches="tight" if show_pattern_labels else None)
    plt.close(fig)


def _heatmap_rasterize():
    """Whether the heatmap image and cluster strips go through the raster backend.

    True, the default, re-renders them through the raster backend at the save dpi. At the default
    2400 dpi a 27,997-cell heatmap draws 11,700 image columns and 13,835 rows, so every marker row
    gets about 6 pixels of height instead of the single pixel it gets at full matrix resolution.
    False keeps matplotlib's own behaviour, one image pixel per matrix cell: exact per cell, but
    one pixel per marker row. Override with MARKER_HEATMAP_RASTERIZE=0.
    """
    return os.environ.get("MARKER_HEATMAP_RASTERIZE", "1").strip().lower() not in ("0", "false", "no")


def _heatmap_dpi():
    """Raster resolution for the heatmap stamp, in dots per inch.

    Only the heatmap image and the cluster strips rasterize; every label, axis and legend stays
    vector at any setting. Raising this sharpens the image and grows the file roughly with the
    square of the value. The default is 2400, which held a 27,997-cell UPenn heatmap to 4.1 MB
    with PDF compression on. Override with MARKER_HEATMAP_DPI or the --dpi flag.
    """
    try:
        dpi = float(os.environ.get("MARKER_HEATMAP_DPI", "2400"))
    except (TypeError, ValueError):
        return 2400.0
    return dpi if dpi > 0 else 2400.0


def _save_figure(fig, output_path, dpi=None, write_svg=True, bbox_inches=None):
    """Write the heatmap and, by default, a sibling `.svg`.

    Illustrator opens the SVG with the text still text, because `svg.fonttype='none'` keeps glyphs
    as characters instead of converting them to outlines, which matches what `pdf.fonttype=42`
    does for the PDF. Returns the paths written.

    An SVG failure never costs the caller the PDF: the primary write happens first, and a failed
    SVG is reported rather than raised.
    """
    written = [output_path]
    dpi = float(dpi) if dpi else _heatmap_dpi()
    fig.savefig(output_path, facecolor="white", transparent=False, dpi=dpi, bbox_inches=bbox_inches)
    base, ext = os.path.splitext(output_path)
    svg_path = base + ".svg"
    if write_svg and ext.lower() != ".svg":
        prior = plt.rcParams.get("svg.fonttype")
        try:
            plt.rcParams["svg.fonttype"] = "none"      # keep text editable
            fig.savefig(svg_path, facecolor="white", transparent=False, dpi=dpi, bbox_inches=bbox_inches)
            written.append(svg_path)
            print(f"Saved marker heatmap SVG to: {svg_path}")
        except Exception as exc:                        # never lose the primary figure over the SVG
            print(f"[WARN] SVG export failed for {svg_path}: {exc}")
        finally:
            plt.rcParams["svg.fonttype"] = prior
    return written


def _load_marker_finder():
    try:
        from cellHarmony import markerFinder as ch_marker_finder

        return ch_marker_finder
    except ImportError:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        components_dir = os.path.dirname(script_dir)
        if components_dir not in sys.path:
            sys.path.insert(0, components_dir)
        try:
            from cellHarmony import markerFinder as ch_marker_finder

            return ch_marker_finder
        except ImportError as exc:
            raise ImportError(
                "Unable to import cellHarmony.markerFinder. Ensure the altanalyze3 components directory is on PYTHONPATH."
            ) from exc


def _markerfinder_stats(adata, cluster_key, use_raw, layer, scale_data=False,
                        scale_factor=None, validate_scaling=True, max_missing_fraction=0.15):
    """MarkerFinder correlations for a cells x genes matrix.

    ``validate_scaling`` defaults to True. MarkerFinder correlates each gene against a 0/1
    cluster indicator, which is sensitive to per-cell sequencing depth, so raw counts and
    log1p-without-normalization layers give markers ranked partly by cell depth.
    """
    marker_finder = _load_marker_finder()
    matrix, gene_names = _get_expression_matrix(adata, use_raw, layer)
    clusters = adata.obs[cluster_key].astype(str).tolist()
    r_df, p_df = marker_finder.marker_finder(
        matrix,
        clusters,
        gene_names=gene_names,
        scale_data=scale_data,
        scale_factor=scale_factor,
        validate_scaling=validate_scaling,
        max_missing_fraction=max_missing_fraction,
    )
    return p_df, r_df


def _generate_single_marker_network(cluster, selected, interactions_df, output_prefix):
    try:
        pdf_path, png_path, tsv_path = NetPerspective.generate_network_for_genes(
            selected,
            interactions_df,
            output_prefix,
            gene_column="gene",
            fold_change_column="log2fc",
            pval_column="fdr" if "fdr" in selected.columns else None,
            max_genes=None,
        )
    except NetPerspective.NetworkGenerationError as exc:
        print(f"[WARN] No marker interaction edges found for {cluster}; skipping network export. {exc}")
        return None
    except ImportError as exc:
        print(f"[WARN] Marker network export disabled: {exc}")
        return None
    return {
        "id": os.path.basename(output_prefix),
        "population": str(cluster),
        "pdf": pdf_path,
        "png": png_path,
        "tsv": tsv_path,
    }


def _export_marker_networks(marker_stats, out_dir, network_top_n=1000, network_jobs=1, species=None):
    if marker_stats is None or marker_stats.empty:
        return []
    try:
        interactions_df = NetPerspective.load_interaction_data(species=species)
    except Exception as exc:
        print(f"[WARN] Marker network export skipped: {exc}")
        return []

    out_dir = str(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    networks = []
    used_ids = set()
    working = marker_stats.copy()
    working["cluster"] = working["cluster"].astype(str)
    working["Gene"] = working["Gene"].astype(str)
    working["Fold"] = pd.to_numeric(working.get("Fold"), errors="coerce")
    working["FDR p-value"] = pd.to_numeric(working.get("FDR p-value"), errors="coerce")

    jobs = []
    for index, (cluster, cluster_df) in enumerate(working.groupby("cluster"), start=1):
        selected = (
            cluster_df.loc[:, [col for col in ("Gene", "Fold", "FDR p-value") if col in cluster_df.columns]]
            .dropna(subset=["Gene", "Fold"])
            .sort_values(["FDR p-value", "Fold", "Gene"], ascending=[True, False, True])
            .drop_duplicates(subset=["Gene"])
            .head(int(network_top_n))
            .rename(columns={"Gene": "gene", "Fold": "log2fc", "FDR p-value": "fdr"})
        )
        if selected.empty or selected["gene"].nunique() < 2:
            continue

        network_id = NetPerspective.safe_component(str(cluster), fallback=f"network_{index}")
        if network_id in used_ids:
            network_id = f"{network_id}_{index}"
        used_ids.add(network_id)
        output_prefix = os.path.join(out_dir, network_id)
        jobs.append((index, cluster, selected, output_prefix))

    if not jobs:
        return networks

    max_workers = max(1, int(network_jobs or 1))
    if max_workers <= 1:
        for _, cluster, selected, output_prefix in jobs:
            network = _generate_single_marker_network(cluster, selected, interactions_df, output_prefix)
            if network is not None:
                networks.append(network)
        return networks

    indexed_outputs = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(_generate_single_marker_network, cluster, selected, interactions_df, output_prefix): idx
            for idx, cluster, selected, output_prefix in jobs
        }
        for future in as_completed(future_map):
            idx = future_map[future]
            try:
                network = future.result()
            except Exception as exc:
                print(f"[WARN] Marker network export task failed: {exc}")
                continue
            if network is not None:
                indexed_outputs.append((idx, network))
    indexed_outputs.sort(key=lambda x: x[0])
    networks.extend([entry for _, entry in indexed_outputs])
    return networks


class NoMarkersSelectedError(ValueError):
    """Valid input produced no features passing marker selection."""


def generate_marker_heatmap_from_adata(
    adata,
    *,
    cluster_key,
    out,
    top_n=5,
    markers_tsv=None,
    heatmap_tsv=None,
    heatmap_column_tsv=None,
    heatmap_cache=None,
    marker_method="scanpy",
    method="wilcoxon",
    use_raw=False,
    layer=None,
    scale_data=False,
    scale_factor=None,
    validate_scaling=True,
    max_missing_fraction=0.15,
    min_rho=None,
    collapse_shared_exons=False,
    max_sample_missing_above_group_median=None,
    label_events=None,
    label_shared_exons=True,
    centroid_method="log2_cp10k",
    cells_per_cluster=50,
    seed=0,
    export_networks=False,
    network_top_n=1000,
    network_jobs=1,
    species=None,
    write_heatmap_tsv=False,
    write_expression_tsv=True,
    write_heatmap_cache=True,
    render_heatmap=True,
    write_svg=True,
    heatmap_dpi=None,
    pval_threshold=0.001,
    covariate_columns=None,
    go_terms=None,
    go_terms_max=30,
    return_statistics=False,
    marker_ranking="fdr",
    missing_psi="preserve",
):
    """Find markers using all input cells; output controls do not change scoring.

    ``render_heatmap=False`` skips static figures while retaining the optional
    interactive cache. Disable all three matrix exports as well to avoid sampling
    and constructing a plotting matrix entirely. Marker tables and compact
    marker-by-state centroids are retained in either mode.
    """
    if heatmap_dpi is not None and (not np.isfinite(heatmap_dpi) or heatmap_dpi <= 0):
        raise ValueError("heatmap_dpi must be a finite positive number.")
    total_started = time.perf_counter()
    timings = {}
    if marker_ranking not in {"fdr", "rho"}:
        raise ValueError("marker_ranking must be fdr or rho")
    if marker_ranking == "rho" and (marker_method != "markerfinder" or pval_threshold is not None):
        raise ValueError("Rho-only ranking requires MarkerFinder and pval_threshold=None")
    if missing_psi not in {"preserve", "row_median"}:
        raise ValueError("missing_psi must be preserve or row_median")
    if missing_psi == "row_median" and (marker_method != "markerfinder" or use_raw or scale_data):
        raise ValueError("Row-median PSI filling requires MarkerFinder, use_raw=False and scale_data=False")
    out = str(out)
    out_dir = os.path.dirname(out) or "."
    base_name = os.path.splitext(os.path.basename(out))[0]
    markers_tsv = markers_tsv or os.path.join(out_dir, f"{base_name}_markers.tsv")
    heatmap_tsv = (
        heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
        if write_heatmap_tsv
        else None
    )
    heatmap_column_tsv = (
        heatmap_column_tsv or os.path.join(out_dir, f"{base_name}_exp_matrix.tsv")
        if write_expression_tsv
        else None
    )
    heatmap_cache = (
        heatmap_cache or os.path.join(out_dir, f"{base_name}_fold_matrix.h5ad")
        if write_heatmap_cache
        else None
    )
    centroids_base = heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
    centroids_tsv = _append_suffix_before_extension(centroids_base, ".centroids")
    redundant_markers_tsv = os.path.join(out_dir, f"{base_name}_redundant_markers.tsv")

    if cluster_key not in adata.obs:
        available_obs = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise KeyError(
            f"Cluster key '{cluster_key}' not found in adata.obs. "
            f"Available obs columns: {available_obs}"
        )

    if use_raw and adata.raw is None:
        raise ValueError("Requested use_raw=True but adata.raw is not set.")

    if layer and layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in adata.layers.")

    if max_sample_missing_above_group_median is not None:
        adata = _drop_high_missing_samples(
            adata, cluster_key, use_raw, layer, max_sample_missing_above_group_median,
            os.path.join(out_dir, f"{base_name}_sample_missing_filter.tsv"),
        )

    source_adata = adata
    imputation = {"method": missing_psi, "filled_cells": 0}
    if missing_psi == "row_median":
        source_matrix, source_features = _get_expression_matrix(adata, False, layer)
        values = source_matrix.toarray() if sparse.issparse(source_matrix) else np.asarray(source_matrix)
        frame = pd.DataFrame(values.T, index=source_features, columns=adata.obs_names)
        filled, filled_mask = fill_PSI_row_median(frame)
        adata = adata.copy()
        if layer:
            adata.layers[layer] = filled.to_numpy().T.copy()
        else:
            adata.X = filled.to_numpy().T.copy()
        imputation.update(original_missing_cells=int(frame.isna().to_numpy().sum()),
                          filled_cells=int(filled_mask.to_numpy().sum()),
                          all_missing_events=int(frame.isna().all(axis=1).sum()))

    resolve_cluster_started = time.perf_counter()
    raw_clusters = adata.obs[cluster_key]
    lineage_order = _coerce_lineage_order(adata.uns.get("lineage_order", None))
    cluster_order = _resolve_cluster_order(raw_clusters, lineage_order)
    clusters = raw_clusters.astype(str)
    print(f"[INFO] Resolved {len(cluster_order)} clusters.")
    _log_step_timing(
        "marker_heatmap.resolve_cluster_order",
        resolve_cluster_started,
        timings=timings,
        key="resolve_cluster_order",
    )

    adata.obs[cluster_key] = pd.Categorical(clusters, categories=cluster_order, ordered=True)

    effect_df = None
    fdr_df = None
    marker_core_started = time.perf_counter()
    if marker_method == "markerfinder":
        print(f"[INFO] Running markerFinder using {'raw' if use_raw else (layer or 'X')} matrix.")
        pvals, effect_df = _markerfinder_stats(
            adata, cluster_key, use_raw, layer,
            scale_data=scale_data, scale_factor=scale_factor,
            validate_scaling=validate_scaling,
            max_missing_fraction=max_missing_fraction,
        )
        fdr_df = pvals.copy() if pvals.empty else pvals.apply(_bh_fdr, axis=0)
    else:
        n_genes = adata.raw.n_vars if use_raw and adata.raw is not None else adata.n_vars
        print(f"[INFO] Running rank_genes_groups for {len(cluster_order)} clusters with {n_genes} genes.")
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=method,
            use_raw=use_raw,
            layer=layer,
            n_genes=n_genes,
        )
        rg_df = sc.get.rank_genes_groups_df(adata, group=None)
        if rg_df.empty:
            raise ValueError("rank_genes_groups produced no results.")
        pval_col = "pvals_adj"
        if pval_col not in rg_df.columns or not rg_df[pval_col].notna().any():
            pval_col = "pvals"
        pvals = rg_df.pivot_table(index="names", columns="group", values=pval_col, aggfunc="first")
        if pval_col == "pvals_adj":
            fdr_df = pvals.copy()
        else:
            fdr_df = pvals.copy() if pvals.empty else pvals.apply(_bh_fdr, axis=0)
        if "logfoldchanges" not in rg_df.columns:
            raise ValueError("logfoldchanges missing; cannot enforce upregulated marker selection.")
        effect_df = rg_df.pivot_table(index="names", columns="group", values="logfoldchanges", aggfunc="first")
    _log_step_timing(
        "marker_heatmap.compute_markers",
        marker_core_started,
        timings=timings,
        key="compute_markers",
    )

    unique_select_started = time.perf_counter()
    # The marker TABLE keeps every assigned gene; the HEATMAP takes top_n per cluster, because
    # a figure with one row per assigned gene is unreadable. The two were the same object, so
    # the figure's size silently defined the table.
    selected_all = _select_unique_markers(
        fdr_df,
        cluster_order,
        None,
        effect_df=effect_df,
        pval_threshold=pval_threshold,
        min_rho=min_rho,
        ranking=marker_ranking,
    )
    if collapse_shared_exons and not selected_all.empty:
        n_before = len(selected_all)
        selected_all, excluded_events = _collapse_shared_exon_events(selected_all)
        excluded_tsv = os.path.join(out_dir, f"{base_name}_shared_exon_excluded.tsv")
        excluded_events.to_csv(excluded_tsv, sep="\t", index=False)
        print(f"[INFO] Shared-exon collapse: kept {len(selected_all)} of {n_before} assigned events; "
              f"{len(excluded_events)} share an exon with a more significant event ({excluded_tsv}).")
    selected = (
        selected_all.groupby("cluster", sort=False).head(top_n).reset_index(drop=True)
        if top_n is not None else selected_all
    )
    selection_description = "rho-only selection" if marker_ranking == "rho" else "FDR/effect filtering"
    print(f"[INFO] Assigned {selected_all.shape[0]} unique markers after {selection_description} "
          f"(written to the marker TSV); {selected.shape[0]} go to the heatmap (top_n={top_n}).")
    _log_step_timing(
        "marker_heatmap.select_unique_markers",
        unique_select_started,
        timings=timings,
        key="select_unique_markers",
    )
    if selected_all.empty:
        if return_statistics and marker_method == "markerfinder":
            return {"pdf": None, "rho": effect_df.reindex(adata.var_names),
                    "pval": pvals.reindex(adata.var_names), "fdr": fdr_df.reindex(adata.var_names),
                    "sample_ids": adata.obs_names.tolist(), "status": "NO_MARKERS_SELECTED"}
        raise NoMarkersSelectedError("No markers were selected. Check inputs and parameters.")

    redundant_select_started = time.perf_counter()
    redundant_selected = _select_top_markers_per_cluster(
        fdr_df,
        cluster_order,
        250,
        effect_df=effect_df,
    )
    print(f"[INFO] Selected {redundant_selected.shape[0]} redundant markers for network export.")
    _log_step_timing(
        "marker_heatmap.select_redundant_markers",
        redundant_select_started,
        timings=timings,
        key="select_redundant_markers",
    )

    selected_genes = selected["gene"].tolist()
    needs_plot_matrix = render_heatmap or write_heatmap_cache or write_heatmap_tsv or write_expression_tsv
    if needs_plot_matrix:
        if cells_per_cluster and cells_per_cluster > 0:
            print(f"[INFO] Downsampling to {cells_per_cluster} cells per cluster (seed={seed}).")
        else:
            print("[INFO] Using all cells (no downsampling).")
        downsample_started = time.perf_counter()
        adata_plot = downsample_cells_per_group(
            adata,
            cluster_key,
            cells_per_cluster=cells_per_cluster,
            seed=seed,
            group_order=cluster_order,
        )
        _log_step_timing(
            "marker_heatmap.downsample_cells",
            downsample_started,
            timings=timings,
            key="downsample_cells",
        )
        print(f"[INFO] Heatmap will use {adata_plot.n_obs} cells.")
        adata_plot.obs[cluster_key] = pd.Categorical(
            adata_plot.obs[cluster_key].astype(str),
            categories=cluster_order,
            ordered=True,
        )

        print(f"[INFO] Building heatmap matrix for {len(selected_genes)} markers.")
        build_heatmap_started = time.perf_counter()
        heatmap_df, heatmap_col_df, cluster_counts, ordered_cells = _build_heatmap(
            adata_plot,
            cluster_key,
            selected_genes,
            cluster_order,
            use_raw,
            layer,
            include_column_zscore=write_expression_tsv,
        )
        _log_step_timing(
            "marker_heatmap.build_heatmap_matrix",
            build_heatmap_started,
            timings=timings,
            key="build_heatmap_matrix",
        )
        if heatmap_df.empty:
            raise ValueError("Heatmap data is empty after filtering genes.")

        column_clusters = adata_plot.obs[cluster_key].astype(str).loc[ordered_cells].tolist()
        row_clusters = selected.set_index("gene").loc[heatmap_df.index, "cluster"].astype(str).tolist()
        covariate_df = None
        if covariate_columns:
            covariate_columns = [str(c).strip() for c in covariate_columns if str(c).strip()]
            covariate_columns = [c for c in dict.fromkeys(covariate_columns) if c in adata_plot.obs]
            if covariate_columns:
                retained_covariates, skipped_covariates = _filter_covariate_columns_for_heatmap(
                    adata_plot.obs,
                    ordered_cells,
                    covariate_columns,
                    max_categories=12,
                )
                if skipped_covariates:
                    skipped_text = ", ".join(f"{name} ({count})" for name, count in skipped_covariates)
                    print(f"[INFO] Skipping heatmap covariates with >12 displayed categories: {skipped_text}")
                if retained_covariates:
                    covariate_df = adata_plot.obs.loc[ordered_cells, retained_covariates].copy()
                    print(f"[INFO] Rendering heatmap covariate bars: {', '.join(retained_covariates)}")
                else:
                    print("[INFO] No requested heatmap covariates passed the <=12 category display limit.")
            else:
                print("[INFO] No requested heatmap covariates were present in adata.obs.")
        if render_heatmap:
            print("[INFO] Rendering heatmap.")
            render_heatmap_started = time.perf_counter()
            _plot_heatmap(
                heatmap_df,
                out,
                cluster_counts,
                cluster_order,
                column_clusters,
                row_clusters,
                covariate_df=covariate_df,
                go_terms=go_terms,
                go_terms_max=go_terms_max,
                dpi=heatmap_dpi,
                write_svg=write_svg,
                label_events=label_events,
                species=species,
                label_shared_exons=label_shared_exons,
            )
            _log_step_timing(
                "marker_heatmap.render_heatmap_pdf",
                render_heatmap_started,
                timings=timings,
                key="render_heatmap_pdf",
            )
            print(f"Saved marker heatmap to: {out}")

    else:
        print("[INFO] Skipping plot matrix construction and static heatmap rendering.")

    print("[INFO] Computing marker statistics for TSV.")
    aggregate_stats_started = time.perf_counter()
    stat_genes = pd.Index(
        pd.concat(
            [
                selected_all["gene"].astype(str),
                redundant_selected["gene"].astype(str),
            ],
            ignore_index=True,
        )
    ).drop_duplicates().tolist()
    stats_aggregates = _prepare_marker_stats_aggregates(
        adata,
        cluster_key,
        stat_genes,
        use_raw,
        layer,
    )
    _log_step_timing(
        "marker_heatmap.prepare_stats_aggregates",
        aggregate_stats_started,
        timings=timings,
        key="prepare_stats_aggregates",
    )

    marker_stats_started = time.perf_counter()
    marker_stats = _compute_marker_stats(
        adata,
        cluster_key,
        selected_all,
        fdr_df,
        use_raw,
        layer,
        aggregates=stats_aggregates,
        effect_column="rho" if marker_method == "markerfinder" else None,
    )
    _log_step_timing(
        "marker_heatmap.compute_marker_stats_unique",
        marker_stats_started,
        timings=timings,
        key="compute_marker_stats_unique",
    )
    write_marker_stats_started = time.perf_counter()
    marker_stats.to_csv(markers_tsv, sep="\t", index=False, float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.write_marker_stats_tsv",
        write_marker_stats_started,
        timings=timings,
        key="write_marker_stats_tsv",
    )
    print(f"Saved marker stats TSV to: {markers_tsv}")

    redundant_stats_started = time.perf_counter()
    redundant_marker_stats = _compute_marker_stats(
        adata,
        cluster_key,
        redundant_selected,
        fdr_df,
        use_raw,
        layer,
        aggregates=stats_aggregates,
        effect_column="rho" if marker_method == "markerfinder" else None,
    )
    _log_step_timing(
        "marker_heatmap.compute_marker_stats_redundant",
        redundant_stats_started,
        timings=timings,
        key="compute_marker_stats_redundant",
    )
    write_redundant_stats_started = time.perf_counter()
    redundant_marker_stats.to_csv(redundant_markers_tsv, sep="\t", index=False, float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.write_redundant_marker_stats_tsv",
        write_redundant_stats_started,
        timings=timings,
        key="write_redundant_marker_stats_tsv",
    )
    print(f"Saved redundant marker stats TSV to: {redundant_markers_tsv}")

    write_heatmap_matrices_started = time.perf_counter()
    if write_heatmap_cache and heatmap_cache:
        write_cache_started = time.perf_counter()
        source_PSI = None
        imputed_mask = None
        if missing_psi == "row_median":
            source_matrix, source_features = _get_expression_matrix(source_adata, False, layer)
            values = source_matrix[source_adata.obs_names.get_indexer(ordered_cells), :]
            values = values[:, pd.Index(source_features).get_indexer(heatmap_df.index)]
            values = values.toarray() if sparse.issparse(values) else np.asarray(values)
            source_PSI = pd.DataFrame(values.T, index=heatmap_df.index, columns=ordered_cells)
            imputed_mask = source_PSI.isna() & source_PSI.median(axis=1).notna().to_numpy()[:, None]
        _write_heatmap_cache(
            heatmap_cache,
            heatmap_df,
            row_clusters,
            column_clusters,
            ordered_cells,
            source_PSI=source_PSI,
            imputed_mask=imputed_mask,
        )
        _log_step_timing(
            "marker_heatmap.write_heatmap_cache",
            write_cache_started,
            timings=timings,
            key="write_heatmap_cache",
        )
        print(f"Saved heatmap cache to: {heatmap_cache}")

    if write_heatmap_tsv and heatmap_tsv:
        heatmap_tsv_df = heatmap_df.copy()
        heatmap_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_tsv_df.index)]
        heatmap_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_tsv_df.to_csv(heatmap_tsv, sep="\t", float_format="%.4g")
        print(f"Saved heatmap matrix TSV to: {heatmap_tsv}")
    else:
        print("[INFO] Skipping fold heatmap TSV export.")

    if write_expression_tsv and heatmap_column_tsv:
        heatmap_col_tsv_df = heatmap_col_df.copy()
        heatmap_col_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_col_tsv_df.index)]
        heatmap_col_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_col_tsv_df.to_csv(heatmap_column_tsv, sep="\t", float_format="%.4g")
    else:
        print("[INFO] Skipping expression matrix TSV export.")
    _log_step_timing(
        "marker_heatmap.write_heatmap_tsvs",
        write_heatmap_matrices_started,
        timings=timings,
        key="write_heatmap_tsvs",
    )
    if write_expression_tsv and heatmap_column_tsv:
        print(f"Saved expression matrix TSV to: {heatmap_column_tsv}")

    centroid_started = time.perf_counter()
    centroid_df = _build_marker_centroids_from_aggregates(
        selected_genes,
        cluster_order,
        stats_aggregates,
        method=centroid_method,
    )
    if centroid_method != "mean":
        centroid_df.index = [_strip_marker_cluster_prefix(gene) for gene in centroid_df.index]
    centroid_df.to_csv(centroids_tsv, sep="\t", float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.build_and_write_centroids_tsv",
        centroid_started,
        timings=timings,
        key="build_and_write_centroids_tsv",
    )
    print(f"Saved centroid matrix TSV to: {centroids_tsv}")

    networks = []
    if export_networks:
        network_dir = os.path.join(out_dir, f"{base_name}_networks")
        print(f"[INFO] Exporting NetPerspective networks for up to {int(network_top_n)} markers per cluster.")
        networks_started = time.perf_counter()
        networks = _export_marker_networks(
            redundant_marker_stats,
            network_dir,
            network_top_n=network_top_n,
            network_jobs=network_jobs,
            species=species,
        )
        _log_step_timing(
            "marker_heatmap.export_networks",
            networks_started,
            timings=timings,
            key="export_networks",
        )

    _log_step_timing(
        "marker_heatmap.total",
        total_started,
        timings=timings,
        key="total",
    )

    result = {
        "pdf": out if render_heatmap else None,
        "markers_tsv": markers_tsv,
        "redundant_markers_tsv": redundant_markers_tsv,
        "heatmap_tsv": heatmap_tsv,
        "heatmap_column_tsv": heatmap_column_tsv,
        "heatmap_cache": heatmap_cache,
        "centroids_tsv": centroids_tsv,
        "networks": networks,
        "timings": timings,
        "marker_ranking": marker_ranking,
        "missing_PSI_handling": imputation,
    }
    if return_statistics:
        result.update(rho=effect_df.reindex(adata.var_names), pval=pvals.reindex(adata.var_names),
                      fdr=fdr_df.reindex(adata.var_names), sample_ids=adata.obs_names.tolist(),
                      assignments=selected_all.copy())
    return result


def render_signed_marker_heatmap_from_saved_statistics(
        adata, rho, pval, fdr, *, cluster_key, out, top_n=100, min_rho=.4,
        pval_threshold=1., label_events=None, dpi=300, caption=None):
    """Extend a native single-scheme view with inverse signatures, without rescoring.

    Reuses the native unique-marker selection, row z-score builder and renderer.
    Each event has one winning signed group, selected by the original FDR/rho rule.
    Full original statistics and source observations must be supplied in exact order.
    """
    from pathlib import Path
    for frame in (rho, pval, fdr):
        if not frame.index.equals(adata.var_names) or not frame.columns.equals(rho.columns):
            raise ValueError('Saved statistics must retain every ordered source event and group')
    order = _resolve_cluster_order(adata.obs[cluster_key],
                                  _coerce_lineage_order(adata.uns.get('lineage_order')))
    if set(order) != set(rho.columns):
        raise ValueError('Saved statistic groups differ from the complete source roster')
    roles = [(direction, group) for direction in ('up', 'down') for group in order]
    effects = pd.DataFrame({direction+' '+group: rho[group] * (1 if direction=='up' else -1)
                            for direction, group in roles})
    adjusted = pd.DataFrame({direction+' '+group: fdr[group] for direction, group in roles})
    selected_all = _select_unique_markers(adjusted, list(effects.columns), None,
        effect_df=effects, pval_threshold=pval_threshold, min_rho=min_rho)
    selected = selected_all.groupby('cluster', sort=False, observed=True).head(top_n)
    if selected.empty:
        raise NoMarkersSelectedError('No signed markers pass the original selection rule')
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    selected_all.to_csv(out.with_suffix('.all_assignments.tsv'), sep='\t', index=False)
    selected.to_csv(out.with_suffix('.displayed_rows.tsv'), sep='\t', index=False)
    matrix, _, counts, cells = _build_heatmap(adata, cluster_key,
        selected.gene.tolist(), order, False, None, include_column_zscore=False)
    groups = adata.obs.loc[cells, cluster_key].astype(str).tolist()
    row_groups = selected.cluster.astype(str).tolist()
    color_order = list(dict.fromkeys([*order, *row_groups]))
    _plot_heatmap(matrix, str(out), counts, color_order, groups, row_groups, dpi=dpi,
        species='human', label_events=label_events or [], show_pattern_labels=True,
        colorbar_label='PSI row z-score', caption=caption)
    matrix.to_csv(out.with_suffix('.zscore_matrix.tsv'), sep='\t')
    cache = str(out.with_suffix('.zscore_matrix.h5ad'))
    _write_heatmap_cache(cache, matrix, row_groups, groups, cells)
    return dict(pdf=str(out), heatmap_cache=cache, assigned=len(selected_all),
                plotted=len(selected), pattern_counts=selected.cluster.value_counts().to_dict(),
                normalization='Native row z-score; missing standardized values displayed as zero',
                selection='Native unique-marker FDR then rho ordering across up/down signed groups')


def fill_PSI_row_median(matrix):
    """Fill only missing PSI with its event median over all supplied samples.

    Rows are events, columns are samples. Observed zeros remain measurements;
    all-missing events remain missing because they have no observed median.
    The supplied matrix is never modified.
    """
    if np.isinf(matrix.to_numpy()).any():
        raise ValueError('Infinite PSI is invalid')
    mask=matrix.isna()
    result=matrix.T.fillna(matrix.median(axis=1)).T
    imputed=mask & result.notna()
    return result, imputed


def render_consolidated_marker_heatmap(adata, assignments, *, cluster_key, out,
                                       top_n=100, species=None, half_range=0.2,
                                       covariate_columns=None, dpi=None, label_events=None,
                                       label_replacements=None, figure_title=None, caption=None,
                                       sample_label_key=None, show_pattern_labels=False,
                                       retain_callout_events=True, retain_unassigned_callouts=False,
                                       label_shared_exons=True, normalization='row_median',
                                       missing_psi='preserve'):
    """Render native iterative assignments with an explicit display normalization.

    ``row_median`` preserves the historical median-centred PSI display. ``row_zscore``
    uses the standard native row z-score (population SD), retaining missing cache cells.

    Assignment order is supplied by iterative MarkerFinder, not re-ranked by a renderer.
    Missing measurements remain missing in the matrix/cache and are never converted to PSI zero.
    ``top_n=None`` displays all assignments; the full assignment table is always retained.
    Only supplied discovery assignments are displayed by default, including qualifying
    callouts beyond the cap. Unassigned reference rows require explicit historical opt-in.
    """
    from pathlib import Path
    out = str(out); root = Path(out).parent; root.mkdir(parents=True, exist_ok=True)
    if not adata.var_names.is_unique or not adata.obs_names.is_unique:
        raise ValueError("Consolidation requires unique feature and pseudobulk identities.")
    if not assignments.index.is_unique or not assignments.index.isin(adata.var_names).all():
        raise ValueError("Every assigned event must have exactly one source row.")
    if top_n is not None and (not isinstance(top_n, (int, np.integer)) or top_n <= 0):
        raise ValueError('top_n must be a positive integer or None')
    if sample_label_key is not None and sample_label_key not in adata.obs:
        raise ValueError('The sample-label column must exist in the verified input metadata')
    if normalization not in {'row_median', 'row_zscore'}:
        raise ValueError('normalization must be row_median or row_zscore')
    if missing_psi not in {'preserve', 'row_median'}:
        raise ValueError('missing_psi must be preserve or row_median')
    assignments.to_csv(Path(out).with_suffix('.assignments.tsv'), sep='\t')
    selected = assignments if top_n is None else assignments.groupby('pattern', sort=False).head(top_n)
    callout_features = []
    if label_events and (retain_callout_events or retain_unassigned_callouts):
        listed = {_event_key(e) for e in label_events}
        listed_exons = set().union(*(_event_exons(_event_uid(e)) for e in label_events))
        callout_features = [feature for feature in adata.var_names
            if _event_key(feature) in listed or (label_shared_exons and
                _event_exons(_event_uid(feature)) & listed_exons)]
    if retain_callout_events and callout_features:
        keep = set(selected.index)
        keep.update(set(callout_features) & set(assignments.index))
        selected = assignments.loc[assignments.index.isin(keep)]
    references = [feature for feature in callout_features if feature not in assignments.index]
    selected = selected.copy()
    selected['display_role'] = 'native_discovery_assignment'
    if retain_unassigned_callouts and references:
        reference_rows = pd.DataFrame(index=pd.Index(references, name=assignments.index.name),
            columns=selected.columns)
        reference_rows['pattern'] = 'Callout reference (not discovery-qualified)'
        reference_rows['display_role'] = 'unassigned_callout_reference'
        selected = pd.concat([selected, reference_rows])
    selected.to_csv(Path(out).with_suffix('.displayed_rows.tsv'), sep='\t')
    if selected.empty:
        return {'pdf': None, 'status': 'NO_MARKERS_SELECTED', 'assigned': 0}
    order = _resolve_cluster_order(adata.obs[cluster_key], _coerce_lineage_order(adata.uns.get('lineage_order')))
    cells, counts = _order_cells_by_cluster(adata, cluster_key, order)
    block = adata[cells, selected.index].X
    block = block.toarray() if sparse.issparse(block) else np.asarray(block)
    matrix = pd.DataFrame(block.T, index=selected.index, columns=cells)
    source_PSI=matrix.copy()
    imputed_mask=None
    if missing_psi=='row_median':
        matrix,imputed_mask=fill_PSI_row_median(matrix)
        caption=(caption or '')+(' Missing PSI is filled with the event median over ALL plotted '
            'MDS/control samples before scaling; observed PSI is unchanged. '
            f'{int(imputed_mask.to_numpy().sum())} cells filled; '
            f'{int(matrix.isna().to_numpy().sum())} cells lack an observed median and remain missing. '
            'Original PSI and imputation mask are retained in cache layers.')
    if normalization == 'row_zscore':
        matrix = _zscore_rows(matrix).where(matrix.notna())
        colorbar_label = 'PSI row z-score'
        normalization_description = 'row z-score (population SD); original missing values retained'
    else:
        matrix = matrix.sub(matrix.median(axis=1), axis=0)
        colorbar_label = 'ΔPSI from row median'
        normalization_description = 'row median subtraction; original missing values retained'
    row_groups = selected['pattern'].astype(str).tolist()
    column_groups = adata.obs.loc[cells, cluster_key].astype(str).tolist()
    color_order = list(dict.fromkeys([*order, *row_groups]))
    cov = adata.obs.loc[cells, covariate_columns] if covariate_columns else None
    old_range = os.environ.get('MARKER_HEATMAP_RANGE')
    old_contrast = os.environ.get('MARKER_HEATMAP_CONTRAST')
    try:
        os.environ['MARKER_HEATMAP_RANGE'] = str(half_range)
        os.environ['MARKER_HEATMAP_CONTRAST'] = '1'
        _plot_heatmap(matrix, out, counts, color_order, column_groups, row_groups,
                      covariate_df=cov, species=species, dpi=dpi,
                      colorbar_label=colorbar_label,
                      label_replacements=label_replacements,
                      figure_title=figure_title, caption=caption,
                      sample_labels=adata.obs[sample_label_key].to_dict() if sample_label_key else None,
                      show_pattern_labels=show_pattern_labels,
                      label_shared_exons=label_shared_exons,
                      label_events=[] if label_events is None else label_events)
    finally:
        for key, value in [('MARKER_HEATMAP_RANGE', old_range), ('MARKER_HEATMAP_CONTRAST', old_contrast)]:
            if value is None: os.environ.pop(key, None)
            else: os.environ[key] = value
    cache = str(Path(out).with_suffix('.fold_matrix.h5ad'))
    _write_heatmap_cache(cache, matrix, row_groups, column_groups, cells,
        source_PSI=source_PSI if missing_psi=='row_median' else None,imputed_mask=imputed_mask)
    matrix.to_csv(Path(out).with_suffix('.fold_matrix.tsv'), sep='\t')
    return dict(pdf=out, heatmap_cache=cache, assigned=len(assignments), plotted=len(selected),
                retain_callout_events=retain_callout_events,
                label_shared_exons=label_shared_exons,
                retain_unassigned_callouts=retain_unassigned_callouts,
                unassigned_callout_references=len(references) if retain_unassigned_callouts else 0,
                normalization=normalization_description,missing_PSI_handling=missing_psi,
                original_missing_cells=int(source_PSI.isna().to_numpy().sum()),
                imputed_cells=int(imputed_mask.to_numpy().sum()) if imputed_mask is not None else 0)


def _load_text_expression_inputs(matrix_tsv, groups_tsv, cluster_key):
    print(f"[INFO] Reading expression matrix TSV from: {matrix_tsv}")
    expr_df = pd.read_csv(matrix_tsv, sep="\t", index_col=0)
    if expr_df.empty:
        raise ValueError(f"Expression matrix '{matrix_tsv}' is empty.")

    expr_df.index = expr_df.index.astype(str)
    expr_df.columns = expr_df.columns.astype(str)
    expr_df = expr_df.loc[~expr_df.index.duplicated(keep="first")]
    expr_df = expr_df.loc[:, ~expr_df.columns.duplicated(keep="first")]
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(
        f"[INFO] Loaded expression matrix with {expr_df.shape[0]} genes and "
        f"{expr_df.shape[1]} cells."
    )

    print(f"[INFO] Reading groups file from: {groups_tsv}")
    groups_df = pd.read_csv(groups_tsv, sep="\t", header=None, dtype=str)
    if groups_df.empty:
        raise ValueError(f"Groups file '{groups_tsv}' is empty.")
    if groups_df.shape[1] < 2:
        raise ValueError(
            f"Groups file '{groups_tsv}' must contain at least 2 columns: cell barcode and group label."
        )

    label_col = 2 if groups_df.shape[1] >= 3 else 1
    groups_df = groups_df.iloc[:, [0, label_col]].copy()
    groups_df.columns = ["cell_id", cluster_key]
    groups_df["cell_id"] = groups_df["cell_id"].astype(str)
    groups_df[cluster_key] = groups_df[cluster_key].astype(str)
    groups_df = groups_df.loc[~groups_df["cell_id"].duplicated(keep="first")]

    group_cells = set(groups_df["cell_id"])
    common_cells = [cell for cell in expr_df.columns if cell in group_cells]
    if not common_cells:
        raise ValueError(
            "No overlapping cell identifiers were found between the expression matrix columns "
            "and the groups file."
        )

    missing_group_count = expr_df.shape[1] - len(common_cells)
    if missing_group_count:
        print(
            f"[INFO] Excluding {missing_group_count} matrix cells not present in the groups file."
        )
    extra_group_count = groups_df.shape[0] - len(common_cells)
    if extra_group_count:
        print(
            f"[INFO] Ignoring {extra_group_count} group rows not present in the expression matrix."
        )

    expr_df = expr_df.loc[:, common_cells]
    groups_df = groups_df.set_index("cell_id").loc[common_cells]

    obs = pd.DataFrame(index=common_cells)
    obs[cluster_key] = groups_df[cluster_key].values
    var = pd.DataFrame(index=expr_df.index.astype(str))
    adata = AnnData(X=expr_df.T.to_numpy(dtype=np.float32, copy=True), obs=obs, var=var)
    adata.obs_names = pd.Index(list(map(str, common_cells)))
    adata.var_names = pd.Index(expr_df.index.astype(str))
    print(
        f"[INFO] Built AnnData from text inputs with {adata.n_obs} cells and {adata.n_vars} genes."
    )
    return adata



def main():
    parser = argparse.ArgumentParser(
        description=(
            "Marker heatmap and MarkerFinder marker table. This is the module the "
            "cellHarmony-web harness uses for MarkerFinder (flask/pipeline.py:1580). "
            "The *_markers.tsv carries every gene assigned to a population at rho > 0 and "
            "FDR <= --fdr-threshold (default 0.001); --top-n bounds the HEATMAP only. Harness values: RNA --top-n 50 "
            "--cells-per-cluster 100 --export-networks --network-top-n 1000; imputed "
            "modalities --top-n 5 --no-scaling-check --centroid-method mean; splicing PSI --psi --top-n 100. "
            "The PSI options change nothing for a matrix without NaN or event IDs, so gene analysis keeps its defaults."
        )
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--h5ad", help="Input h5ad file.")
    input_group.add_argument(
        "--render-from-cache",
        help=(
            "Path to a saved *_fold_matrix.h5ad (or legacy .npz) heatmap cache. Redraws that finished result with "
            "the current plotting code in seconds. Nothing is reclustered and MarkerFinder does "
            "not run again."
        ),
    )
    input_group.add_argument(
        "--matrix-tsv",
        help="Log2-normalized expression matrix TSV with genes as rows and cells as columns.",
    )
    parser.add_argument(
        "--groups-tsv",
        default=None,
        help="Groups file for --matrix-tsv input. Column 1 is cell barcode, column 3 (or 2 if absent) is the cluster label.",
    )
    parser.add_argument(
        "--groups-txt",
        dest="groups_tsv",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--cluster-key",
        default=None,
        help="obs column with cluster labels for h5ad input. Ignored for --matrix-tsv, which uses the groups file.",
    )
    parser.add_argument("--top-n", type=int, default=5, help="Markers per cluster (default 5).")
    parser.add_argument("--out", default="marker_heatmap.pdf", help="Output figure path.")
    parser.add_argument(
        "--covariate-h5ad",
        default=None,
        help=(
            "h5ad supplying obs for the covariate bars during --render-from-cache. The cache "
            "stores the plotted cell order but not the cell metadata. Defaults to "
            "icgs3_result.h5ad in the run directory when that file exists."
        ),
    )
    parser.add_argument(
        "--run-dir",
        default=None,
        help=(
            "ICGS3 run directory holding GO-Elite/icgs3_biomarker_enrichment.tsv, used with "
            "--render-from-cache. Defaults to the cache file's parent run directory."
        ),
    )
    parser.add_argument(
        "--no-go-terms",
        action="store_true",
        help="Draw a cache re-render without the GO-Elite enrichment term labels.",
    )
    parser.add_argument(
        "--go-terms-max",
        type=int,
        default=30,
        help="Maximum GO-Elite terms labelled per cluster (default 30).",
    )
    parser.add_argument(
        "--no-rasterize",
        action="store_true",
        help=(
            "Embed the heatmap at full matrix resolution, one image pixel per cell, instead of "
            "re-rendering it at --dpi."
        ),
    )
    parser.add_argument(
        "--rasterize",
        action="store_true",
        help=(
            "Re-render the heatmap image through the raster backend at --dpi. This is the default."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=float,
        default=None,
        help=(
            "Raster resolution of the heatmap image in dots per inch (default 2400, or "
            "MARKER_HEATMAP_DPI). Labels, axes and legends stay vector at every setting."
        ),
    )
    parser.add_argument("--markers-tsv", default=None, help="Output TSV for marker statistics.")
    parser.add_argument("--heatmap-tsv", default=None, help="Output TSV for the row-scaled heatmap matrix.")
    parser.add_argument(
        "--heatmap-column-tsv",
        default=None,
        help="Output TSV for the column-scaled expression matrix.",
    )
    parser.add_argument(
        "--marker-method",
        choices=["scanpy", "markerfinder"],
        default="scanpy",
        help="Marker selection via scanpy rank_genes_groups or cellHarmony markerFinder.",
    )
    parser.add_argument("--method", default="wilcoxon", help="DE method for rank_genes_groups.")
    parser.add_argument("--use-raw", action="store_true", help="Use adata.raw for DE and heatmap.")
    parser.add_argument("--layer", default=None, help="Layer to use for DE and heatmap.")
    parser.add_argument(
        "--scale-data",
        action="store_true",
        help=(
            "Depth-normalize raw counts before MarkerFinder: "
            "log2(counts / total_counts_per_cell * scale_factor + 1). Required when .X or the "
            "chosen --layer holds raw counts."
        ),
    )
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=None,
        help=(
            "Scale factor for --scale-data. Default: the median cell total, matching scanpy "
            "normalize_total(). Keep the default for targeted panels (Xenium, CosMx, MERFISH)."
        ),
    )
    parser.add_argument(
        "--no-scaling-check",
        action="store_true",
        help="Skip the MarkerFinder input-scaling check (not recommended).",
    )
    parser.add_argument(
        "--cells-per-cluster",
        type=int,
        default=50,
        help="Number of cells to display per cluster (default 50). Use 0 for all cells.",
    )
    parser.add_argument(
        "--export-networks",
        action="store_true",
        help="Export NetPerspective networks for each cluster using the top marker genes.",
    )
    parser.add_argument(
        "--network-top-n",
        type=int,
        default=1000,
        help="Top marker genes per cluster to use for network export (default 1000).",
    )
    parser.add_argument(
        "--network-jobs",
        type=int,
        default=1,
        help="Parallel workers for per-cluster network generation (default 1).",
    )
    parser.add_argument(
        "--species",
        default=None,
        help="Species for GRN interaction loading (human or mouse). Defaults to combined resources when omitted.",
    )
    parser.add_argument(
        "--skip-expression-tsv",
        action="store_true",
        help="Skip writing both expression-matrix and fold-heatmap TSV outputs.",
    )
    parser.add_argument("--skip-heatmap-render", action="store_true",
                        help="Skip static PDF/SVG heatmaps; retain requested matrix/cache outputs.")
    parser.add_argument("--skip-svg", action="store_true",
                        help="Skip the companion SVG when generating a static heatmap.")
    parser.add_argument("--markers-only", action="store_true",
                        help="Write marker statistics and state centroids, without heatmap matrices, caches or figures.")
    parser.add_argument(
        "--heatmap-cache",
        default=None,
        help="Output compressed H5AD cache for selected markers and sampled cells.",
    )
    parser.add_argument(
        "--skip-heatmap-cache",
        action="store_true",
        help="Skip writing the compressed marker heatmap cache.",
    )
    parser.add_argument(
        "--heatmap-covariates",
        default=None,
        help="Comma-delimited h5ad obs columns to render as compact bottom covariate bars.",
    )
    parser.add_argument(
        "--centroid-method",
        choices=["log2_cp10k", "mean"],
        default="log2_cp10k",
        help=(
            "Centroid statistic. 'log2_cp10k' is a pseudobulk CPM and is defined only for "
            "non-negative counts. 'mean' keeps signed per-cluster means and is what the "
            "cellHarmony-web harness passes for imputed modalities such as dsb ADT, lipid and "
            "GRN (flask/pipeline.py:601), whose values are legitimately negative."
        ),
    )
    parser.add_argument(
        "--max-missing-fraction",
        type=float,
        default=0.15,
        help=(
            "MarkerFinder treats NaN as a missing value (e.g. PSI) and scores each feature over its "
            "measured samples, as AltAnalyze2 did. A feature missing in this fraction of samples or "
            "more is not considered (default 0.15). No effect on a matrix without NaN."
        ),
    )
    parser.add_argument(
        "--min-rho",
        type=float,
        default=None,
        help=(
            "Minimum MarkerFinder score (Pearson rho) for a feature to be assigned to a population, "
            "on top of the FDR gate. Default: none (rho > 0). AltAnalyze2 used 0.4 for PSI."
        ),
    )
    parser.add_argument(
        "--collapse-shared-exons",
        action="store_true",
        help=(
            "Splicing events: keep only the most significant (FDR, then rho) of the assigned events "
            "that share any exon; the others go to *_shared_exon_excluded.tsv."
        ),
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.001,
        help=(
            "Maximum Benjamini-Hochberg FDR (per population) for a feature to be assigned "
            "(default 0.001, the API default pval_threshold)."
        ),
    )
    parser.add_argument("--rho-only", action="store_true",
                        help="MarkerFinder: assign and rank by positive rho alone; disable P/FDR selection.")
    parser.add_argument("--missing-psi", choices=["preserve", "row_median"], default="preserve",
                        help="Fill missing event PSI on a copy with its median over all input samples before scoring and scaling.")
    parser.add_argument(
        "--max-sample-missing-above-group-median",
        type=float,
        default=None,
        help=(
            "Remove a sample whose missing-value (NaN) fraction exceeds its cluster's median by more than this "
            "fraction, before scoring; the decision per sample goes to *_sample_missing_filter.tsv. Default: off."
        ),
    )
    parser.add_argument(
        "--label-events",
        default=None,
        help=(
            "Text or TSV file of event IDs (column UID, else the first column). Only heatmap rows matching an "
            "event, in either junction order, or sharing an exon block with one, are labelled, as callouts to "
            "the right; *_label_events.tsv records why. Default: label every row."
        ),
    )
    parser.add_argument(
        "--label-events-exact",
        action="store_true",
        help="With --label-events, label only the listed events themselves, not rows that share an exon with them.",
    )
    parser.add_argument(
        "--psi",
        action="store_true",
        help=(
            "Splicing PSI defaults (PSI_DEFAULTS): --marker-method markerfinder --no-scaling-check "
            "--centroid-method mean --max-missing-fraction 0.55 --min-rho 0.25 --fdr-threshold 0.01 "
            "--collapse-shared-exons --top-n 100 --max-sample-missing-above-group-median 0.10. "
            "A flag given explicitly overrides its preset value."
        ),
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed for cell sampling.")
    preset_args, _ = parser.parse_known_args()
    if preset_args.psi:
        parser.set_defaults(
            marker_method=PSI_DEFAULTS["marker_method"],
            no_scaling_check=not PSI_DEFAULTS["validate_scaling"],
            centroid_method=PSI_DEFAULTS["centroid_method"],
            max_missing_fraction=PSI_DEFAULTS["max_missing_fraction"],
            min_rho=PSI_DEFAULTS["min_rho"],
            fdr_threshold=PSI_DEFAULTS["pval_threshold"],
            collapse_shared_exons=PSI_DEFAULTS["collapse_shared_exons"],
            top_n=PSI_DEFAULTS["top_n"],
            max_sample_missing_above_group_median=PSI_DEFAULTS["max_sample_missing_above_group_median"],
        )
    args = parser.parse_args()
    if args.render_from_cache and (args.markers_only or args.skip_heatmap_render):
        parser.error("--render-from-cache cannot be combined with an option that skips rendering.")
    if args.dpi is not None and (not np.isfinite(args.dpi) or args.dpi <= 0):
        parser.error("--dpi must be a finite positive number.")

    if args.render_from_cache:
        if args.dpi:
            os.environ["MARKER_HEATMAP_DPI"] = str(args.dpi)
        if args.no_rasterize:
            os.environ["MARKER_HEATMAP_RASTERIZE"] = "0"
        elif args.rasterize:
            os.environ["MARKER_HEATMAP_RASTERIZE"] = "1"
        covariate_columns = None
        if getattr(args, "heatmap_covariates", None):
            covariate_columns = [c.strip() for c in str(args.heatmap_covariates).split(",") if c.strip()]
        if covariate_columns and not args.covariate_h5ad:
            guess_dir = args.run_dir or os.path.dirname(
                os.path.dirname(os.path.abspath(args.render_from_cache)))
            guess = os.path.join(guess_dir, "icgs3_result.h5ad")
            if os.path.exists(guess):
                args.covariate_h5ad = guess
                print(f"[INFO] Using {guess} for covariate bars.")
        written = render_heatmap_from_cache(
            args.render_from_cache,
            args.out,
            go_terms_max=getattr(args, "go_terms_max", 30) or 30,
            run_dir=args.run_dir,
            covariate_h5ad=args.covariate_h5ad,
            covariate_columns=covariate_columns,
            load_go_terms=not args.no_go_terms,
            dpi=args.dpi,
            write_svg=not args.skip_svg,
            label_events=_read_label_events(args.label_events) if args.label_events else None,
            species=args.species,
            label_shared_exons=not args.label_events_exact,
        )
        print(f"Saved marker heatmap to: {written}")
        return

    out_dir = os.path.dirname(args.out) or "."
    logs_dir = os.path.join(out_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(args.out))[0]
    log_path = os.path.join(logs_dir, f"{base_name}.log")
    markers_tsv = args.markers_tsv or os.path.join(out_dir, f"{base_name}_markers.tsv")
    heatmap_tsv = (
        args.heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
        if not (args.skip_expression_tsv or args.markers_only)
        else None
    )
    heatmap_column_tsv = (
        args.heatmap_column_tsv or os.path.join(out_dir, f"{base_name}_exp_matrix.tsv")
        if not (args.skip_expression_tsv or args.markers_only)
        else None
    )
    heatmap_cache = (
        args.heatmap_cache or os.path.join(out_dir, f"{base_name}_fold_matrix.h5ad")
        if not (args.skip_heatmap_cache or args.markers_only)
        else None
    )
    centroids_base = heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
    centroids_tsv = _append_suffix_before_extension(centroids_base, ".centroids")
    redundant_markers_tsv = os.path.join(out_dir, f"{base_name}_redundant_markers.tsv")

    log_file = open(log_path, "w")
    tee_out = Tee(sys.stdout, log_file)
    tee_err = Tee(sys.stderr, log_file)
    original_stdout, original_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = tee_out, tee_err

    try:
        print("[INFO] Parameters:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print(f"[INFO] Log file: {log_path}")
        print(f"[INFO] Marker TSV: {markers_tsv}")
        print(f"[INFO] Redundant marker TSV: {redundant_markers_tsv}")
        if heatmap_tsv:
            print(f"[INFO] Heatmap TSV (fold matrix): {heatmap_tsv}")
        if heatmap_column_tsv:
            print(f"[INFO] Heatmap TSV (expression matrix): {heatmap_column_tsv}")
        if heatmap_cache:
            print(f"[INFO] Heatmap cache (h5ad): {heatmap_cache}")
        print(f"[INFO] Centroid TSV: {centroids_tsv}")

        layer = args.layer
        if args.use_raw and layer:
            print("Warning: --layer ignored when --use-raw is set.")
            layer = None

        cluster_key = args.cluster_key

        if args.h5ad:
            if not cluster_key:
                raise ValueError("--cluster-key is required when --h5ad is used.")
            print(f"[INFO] Reading AnnData from: {args.h5ad}")
            adata = sc.read_h5ad(args.h5ad)
            print(f"[INFO] Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes.")
        else:
            if not args.groups_tsv:
                raise ValueError("--groups-tsv is required when --matrix-tsv is used.")
            if args.use_raw:
                print("[WARN] --use-raw ignored for text matrix input.")
            if layer:
                print("[WARN] --layer ignored for text matrix input.")
            cluster_key = "cluster"
            adata = _load_text_expression_inputs(args.matrix_tsv, args.groups_tsv, cluster_key)
            layer = None

        if cluster_key not in adata.obs:
            available_obs = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
            raise KeyError(
                f"Cluster key '{cluster_key}' not found in adata.obs. "
                f"Available obs columns: {available_obs}"
            )

        if args.use_raw and adata.raw is None:
            raise ValueError("Requested --use-raw but adata.raw is not set.")

        if layer and layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")

        outputs = generate_marker_heatmap_from_adata(
            adata,
            cluster_key=cluster_key,
            out=args.out,
            top_n=args.top_n,
            markers_tsv=markers_tsv,
            heatmap_tsv=heatmap_tsv,
            heatmap_column_tsv=heatmap_column_tsv,
            heatmap_cache=heatmap_cache,
            marker_method=args.marker_method,
            method=args.method,
            use_raw=args.use_raw,
            layer=layer,
            scale_data=args.scale_data,
            scale_factor=args.scale_factor,
            centroid_method=args.centroid_method,
            validate_scaling=not args.no_scaling_check,
            max_missing_fraction=args.max_missing_fraction,
            min_rho=args.min_rho,
            collapse_shared_exons=args.collapse_shared_exons,
            pval_threshold=None if args.rho_only else args.fdr_threshold,
            marker_ranking="rho" if args.rho_only else "fdr",
            missing_psi=args.missing_psi,
            max_sample_missing_above_group_median=args.max_sample_missing_above_group_median,
            label_events=_read_label_events(args.label_events) if args.label_events else None,
            label_shared_exons=not args.label_events_exact,
            cells_per_cluster=args.cells_per_cluster,
            seed=args.seed,
            export_networks=args.export_networks,
            network_top_n=args.network_top_n,
            network_jobs=args.network_jobs,
            species=args.species,
            write_heatmap_tsv=not (args.skip_expression_tsv or args.markers_only),
            write_expression_tsv=not (args.skip_expression_tsv or args.markers_only),
            write_heatmap_cache=not (args.skip_heatmap_cache or args.markers_only),
            render_heatmap=not (args.skip_heatmap_render or args.markers_only),
            write_svg=not args.skip_svg,
            heatmap_dpi=args.dpi,
            covariate_columns=[c.strip() for c in str(args.heatmap_covariates or "").split(",") if c.strip()],
        )
        timings = outputs.get("timings", {}) if isinstance(outputs, dict) else {}
        if timings:
            ordered_keys = [
                "resolve_cluster_order",
                "compute_markers",
                "select_unique_markers",
                "select_redundant_markers",
                "downsample_cells",
                "build_heatmap_matrix",
                "render_heatmap_pdf",
                "prepare_stats_aggregates",
                "compute_marker_stats_unique",
                "write_marker_stats_tsv",
                "compute_marker_stats_redundant",
                "write_redundant_marker_stats_tsv",
                "write_heatmap_cache",
                "write_heatmap_tsvs",
                "build_and_write_centroids_tsv",
                "export_networks",
                "total",
            ]
            summary = " | ".join(
                f"{k}={float(timings[k]):.2f}s" for k in ordered_keys if k in timings
            )
            if summary:
                print(f"[timing] marker_heatmap.summary {summary}")

        print(f"Saved log to: {log_path}")
    finally:
        sys.stdout, sys.stderr = original_stdout, original_stderr
        log_file.close()


if __name__ == "__main__":

    main()
