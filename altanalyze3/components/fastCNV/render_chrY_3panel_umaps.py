#!/usr/bin/env python3
"""Render 3-panel chrY UMAPs (LOY vs WT scatter + gene expression in LOY +
gene expression in WT) for an AnnData with precomputed UMAP coordinates.

Panel layout (single row):
  [1] UMAP scatter: WT cells grey, LOY cells red
  [2] gene expression in LOY cells only (grey -> red colormap, q99 vmax)
  [3] gene expression in WT cells only  (same colormap and vmax)

This is the canonical visualization used for fastCNV LOY detection QC.
LOY membership is derived from a fastCNV.cnv_intervals.tsv (chr == 'chrY' &
call == 'loss'); one or more interval files can be passed to support
pooled/per-group cohorts.

CLI:
    python -m altanalyze3.components.fastCNV.render_chrY_3panel_umaps \
        --h5ad PATH \
        --cnv-intervals PATH [PATH ...] \
        --genes RPS4Y1 KDM5D UTY ... \
        --outdir DIR \
        [--label cohort_name]
        [--library-key Library]
        [--per-sample-intervals SAMPLE=PATH ...]
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Sequence
import argparse

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

DEFAULT_CHRY_GENES = [
    'RPS4Y1', 'KDM5D', 'UTY', 'DDX3Y', 'EIF1AY', 'USP9Y', 'ZFY',
    'TMSB4Y', 'NLGN4Y', 'RPS4Y2', 'TBL1Y', 'AMELY', 'TSPY1', 'PRY', 'PRKY',
]
GREY_RED_CMAP = LinearSegmentedColormap.from_list(
    'grey_red',
    ['#d9d9d9', '#fee0d2', '#fc9272', '#ef3b2c', '#a50f15'],
)
LOY_RED = '#d62728'
WT_GREY = '#e0e0e0'
DOT = 4.2
DOT_BG = 3.6


def _load_loy_barcodes(cnv_intervals_paths: Iterable[str | Path]) -> set[str]:
    """Union of chrY-loss barcodes across one or more fastCNV.cnv_intervals.tsv files."""
    loy: set[str] = set()
    for p in cnv_intervals_paths:
        ci = pd.read_csv(p, sep='\t')
        chry = ci[(ci['chr'] == 'chrY') & (ci['call'] == 'loss')]
        loy |= set(chry['CellBarcode'].astype(str).unique())
    return loy


def _load_per_sample_loy(
    per_sample_intervals: dict[str, str | Path],
) -> dict[str, set[str]]:
    """Per-sample LOY barcode sets, keyed by library/sample ID."""
    return {
        sample: _load_loy_barcodes([path])
        for sample, path in per_sample_intervals.items()
    }


def render(
    h5ad: str | Path,
    outdir: str | Path,
    cnv_intervals: Sequence[str | Path] | None = None,
    per_sample_intervals: dict[str, str | Path] | None = None,
    library_key: str = 'Library',
    genes: Sequence[str] = DEFAULT_CHRY_GENES,
    label: str = 'cohort',
    umap_key: str = 'X_umap',
) -> list[Path]:
    """Render chrY 3-panel UMAPs.

    Either pass `cnv_intervals` (a list of TSV paths whose chrY-loss
    barcodes are unioned) OR `per_sample_intervals` (a mapping of
    sample-id -> TSV path, where each cell is matched against the intervals
    for its own sample via obs[library_key]). Use per-sample mode whenever
    the same barcode can repeat across samples.

    Returns the list of written PDF paths.
    """
    h5ad = Path(h5ad)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(h5ad)
    if umap_key not in a.obsm:
        raise ValueError(f"obsm['{umap_key}'] not found in {h5ad}")
    umap = np.asarray(a.obsm[umap_key])
    bcs = a.obs_names.astype(str).values
    var_names = a.var_names.astype(str).values

    if per_sample_intervals:
        if library_key not in a.obs.columns:
            raise ValueError(
                f"obs['{library_key}'] not found; needed for per-sample LOY matching"
            )
        lib = a.obs[library_key].astype(str).values
        sample_loy = _load_per_sample_loy(per_sample_intervals)
        missing = [s for s in np.unique(lib) if s not in sample_loy]
        if missing:
            print(f'  [warn] no LOY intervals supplied for samples: {missing}')
            for s in missing:
                sample_loy[s] = set()
        is_loy = np.array(
            [bc in sample_loy.get(ll, set()) for bc, ll in zip(bcs, lib)]
        )
    elif cnv_intervals:
        loy_set = _load_loy_barcodes(cnv_intervals)
        is_loy = np.isin(bcs, list(loy_set))
    else:
        raise ValueError(
            'Must provide either cnv_intervals (list of TSVs) or '
            'per_sample_intervals (mapping)'
        )
    is_wt = ~is_loy
    n_loy = int(is_loy.sum()); n_wt = int(is_wt.sum()); n_total = umap.shape[0]
    print(f'[{label}] cells={n_total}  LOY={n_loy}  WT={n_wt}')

    xpad = (umap[:, 0].max() - umap[:, 0].min()) * 0.03
    ypad = (umap[:, 1].max() - umap[:, 1].min()) * 0.03
    xlim = (umap[:, 0].min() - xpad, umap[:, 0].max() + xpad)
    ylim = (umap[:, 1].min() - ypad, umap[:, 1].max() + ypad)

    written: list[Path] = []
    for gene in genes:
        if gene not in var_names:
            print(f'  skip {gene} (absent in var)')
            continue
        gi = int(np.where(var_names == gene)[0][0])
        col = a.X[:, gi]
        if sp.issparse(col): col = col.toarray()
        expr = np.asarray(col).ravel()
        nz = expr[expr > 0]
        vmax = float(np.quantile(nz, 0.99)) if nz.size else 1.0

        fig = plt.figure(figsize=(15, 4.5))
        width_ratios = [4.0, 0.6, 4.0, 0.22, 0.6, 4.0, 0.22]
        gs = GridSpec(1, len(width_ratios), figure=fig,
                      width_ratios=width_ratios, wspace=0.05)
        fig.suptitle(f'{label}  -  {gene} expression', fontsize=11, y=1.02)

        ax1 = fig.add_subplot(gs[0, 0])
        ax1.scatter(umap[is_wt, 0], umap[is_wt, 1], c=WT_GREY,
                    s=DOT_BG, linewidths=0)
        ax1.scatter(umap[is_loy, 0], umap[is_loy, 1], c=LOY_RED,
                    s=DOT, linewidths=0)
        ax1.set_title(f'{label}  LOY ({n_loy}/{n_total})', fontsize=9)
        ax1.set_xlim(xlim); ax1.set_ylim(ylim)
        ax1.set_xticks([]); ax1.set_yticks([])
        for sp_ in ax1.spines.values(): sp_.set_visible(False)

        ax2 = fig.add_subplot(gs[0, 2])
        sub = umap[is_loy]; sub_e = expr[is_loy]
        order = np.argsort(sub_e)
        sc2 = ax2.scatter(sub[order, 0], sub[order, 1], c=sub_e[order],
                          cmap=GREY_RED_CMAP, vmin=0, vmax=vmax,
                          s=DOT, linewidths=0)
        ax2.set_title(f'{label}  {gene} in LOY (n={n_loy})', fontsize=9)
        ax2.set_xlim(xlim); ax2.set_ylim(ylim)
        ax2.set_xticks([]); ax2.set_yticks([])
        for sp_ in ax2.spines.values(): sp_.set_visible(False)
        cax2 = fig.add_subplot(gs[0, 3])
        fig.colorbar(sc2, cax=cax2).set_label(f'{gene} (log1p CP10K)', fontsize=7)

        ax3 = fig.add_subplot(gs[0, 5])
        sub = umap[is_wt]; sub_e = expr[is_wt]
        order = np.argsort(sub_e)
        sc3 = ax3.scatter(sub[order, 0], sub[order, 1], c=sub_e[order],
                          cmap=GREY_RED_CMAP, vmin=0, vmax=vmax,
                          s=DOT, linewidths=0)
        ax3.set_title(f'{label}  {gene} in WT (n={n_wt})', fontsize=9)
        ax3.set_xlim(xlim); ax3.set_ylim(ylim)
        ax3.set_xticks([]); ax3.set_yticks([])
        for sp_ in ax3.spines.values(): sp_.set_visible(False)
        cax3 = fig.add_subplot(gs[0, 6])
        fig.colorbar(sc3, cax=cax3).set_label(f'{gene} (log1p CP10K)', fontsize=7)

        out = outdir / f'UMAP_{label}_{gene}_3panel.pdf'
        fig.savefig(out, dpi=200, bbox_inches='tight')
        plt.close(fig)
        pct_loy = (expr[is_loy] > 0).mean() * 100 if n_loy else 0
        pct_wt = (expr[is_wt] > 0).mean() * 100 if n_wt else 0
        print(f'  Wrote {out.name}  LOY %+ ={pct_loy:.2f}  WT %+ ={pct_wt:.2f}')
        written.append(out)
    return written


def _parse_per_sample(values: Sequence[str] | None) -> dict[str, str]:
    if not values:
        return {}
    out: dict[str, str] = {}
    for v in values:
        if '=' not in v:
            raise argparse.ArgumentTypeError(
                f"--per-sample-intervals expects SAMPLE=PATH, got: {v!r}"
            )
        sample, path = v.split('=', 1)
        out[sample] = path
    return out


def main(argv: Sequence[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description='Render 3-panel chrY UMAPs from an AnnData with X_umap '
                    'and fastCNV.cnv_intervals.tsv LOY calls.',
    )
    p.add_argument('--h5ad', required=True, help='Input h5ad with obsm[X_umap]')
    p.add_argument('--outdir', required=True, help='Where to write PDFs')
    p.add_argument(
        '--cnv-intervals', nargs='+',
        help='One or more fastCNV.cnv_intervals.tsv files; chrY-loss '
             'barcodes are unioned across files. Mutually exclusive with '
             '--per-sample-intervals.',
    )
    p.add_argument(
        '--per-sample-intervals', nargs='+', metavar='SAMPLE=PATH',
        help='Per-sample intervals as SAMPLE=PATH pairs (e.g. 657=/path/cnv_intervals.tsv). '
             'Each cell is matched against the intervals for its sample, '
             'using obs[--library-key] to identify the sample.',
    )
    p.add_argument('--library-key', default='Library',
                   help="obs column with sample/library ID (default 'Library')")
    p.add_argument('--genes', nargs='+', default=DEFAULT_CHRY_GENES,
                   help='Genes to plot (defaults to the standard chrY panel)')
    p.add_argument('--label', default='cohort',
                   help="Cohort label used in filenames and titles")
    p.add_argument('--umap-key', default='X_umap',
                   help="obsm key with UMAP coords (default X_umap)")
    args = p.parse_args(argv)

    if not args.cnv_intervals and not args.per_sample_intervals:
        p.error('Must pass --cnv-intervals or --per-sample-intervals')

    render(
        h5ad=args.h5ad,
        outdir=args.outdir,
        cnv_intervals=args.cnv_intervals,
        per_sample_intervals=_parse_per_sample(args.per_sample_intervals) or None,
        library_key=args.library_key,
        genes=args.genes,
        label=args.label,
        umap_key=args.umap_key,
    )


if __name__ == '__main__':
    main()
