#!/usr/bin/env python3
"""Stacked isoform-usage bar/line plot with ISV-web-faithful structure tracks -- a reusable iso2function
callable (ported from the KINNEX-5 ``stack_lineage_struct_cohort_figset.py`` analysis driver).

One panel of per-cell-state (or per-lineage) isoform composition (stacked bar renormalised to 100% among
the shown isoforms, or a mean+/-SE usage line), stacked ABOVE one ISV-web structure track per isoform
(:mod:`.isoform_struct_view`). Layout: the isoform structure tracks sit BELOW the (enlarged, Arial)
cell-state x-axis labels via a GridSpec spacer row so the two never overlap.

Public API:
    plot_stacked_isoform_bar(gene, ensg, isoforms, pseudobulk, out_path, gene_model=…, gff_dir=…, …)
    load_pseudobulk_slice(path)   -> a lightweight AnnData-like slice from a .npz cache or .h5ad
"""
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from scipy.stats import mannwhitneyu      # noqa: E402
import scipy.sparse as sp                 # noqa: E402

from .isoform_struct_view import StructRenderer  # noqa: E402

# distinct, professional (RGB) palette -- one colour per isoform
PALETTE = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#1B9E77",
           "#666666", "#D95F02"]


class _Slice:
    """Minimal AnnData-like holder: ``obs_names`` (list[str]), ``X`` (2-D array/sparse), ``var_names`` (list[str])."""
    __slots__ = ("obs_names", "X", "var_names")

    def __init__(self, obs_names, X, var_names):
        self.obs_names = list(obs_names)
        self.X = X
        self.var_names = list(var_names)


def load_pseudobulk_slice(path, ensg=None):
    """Load a pseudobulk slice from a ``.npz`` cache (keys ``X``/``obs`` plus ``var`` or ``iso``) or an
    ``.h5ad`` (read once). ``ensg`` is only needed to reconstruct var names from a single-gene ``iso`` npz.
    Returns a :class:`_Slice`."""
    if str(path).endswith(".npz"):
        z = np.load(path, allow_pickle=True)
        obs = [str(o) for o in z["obs"]]
        if "var" in z.files:
            var = [str(v) for v in z["var"]]
        else:
            var = [f"{ensg}:{i}" for i in z["iso"]]
        return _Slice(obs, z["X"], var)
    import anndata as ad
    a = ad.read_h5ad(path)
    return _Slice(list(a.obs_names), a.X, list(np.asarray(a.var_names, dtype=str)))


def _sample_of(obs_name):
    """'<cellstate>.<sample>-isoform' -> sample token (matches the pseudobulk obs naming)."""
    return obs_name.split(".", 1)[1].rsplit("-isoform", 1)[0] if "." in obs_name else ""


def _cellstate_of(obs_name):
    return obs_name.split(".", 1)[0]


def _in_cohort(sample, cohort_samples):
    if not cohort_samples:
        return True
    return any(sample == c or sample.startswith(c + "_") or sample.startswith(c + "-") for c in cohort_samples)


def plot_stacked_isoform_bar(gene, ensg, isoforms, pseudobulk, out_path, *,
                             gene_model, gff_dir, level="cellstate", plottype="bar",
                             cohort_samples=None, lineage_map=None, cluster_order=None,
                             clone_map=None, colors=None, intron_scale=0.2,
                             min_group_reads=20, min_donors=2, xlabel_fontsize=None,
                             struct_renderer=None, title=None):
    """Render a stacked isoform-composition plot with structure tracks and write ``out_path`` (.pdf + .png).

    Parameters
    ----------
    gene, ensg : gene symbol and ENSG id.
    isoforms   : ordered list of isoform ids to show (the token after ':' in var_names, e.g. an ENST or a
                 novel '<num>.<sample>' id). Order = left-to-right stack / structure-track order.
    pseudobulk : an AnnData-like object (``.obs_names``/``.X``/``.var_names``) or a path (see
                 :func:`load_pseudobulk_slice`). obs names are '<cellstate>.<sample>-isoform'.
    out_path   : output base path (".pdf"/".png" appended).
    gene_model : Ensembl exon coords TSV (for structure tracks). gff_dir : the dataset gff-output dir.
    level      : 'cellstate' (per cell state) or 'lineage' (aggregate via ``lineage_map``).
    plottype   : 'bar' (100% stacked composition) or 'line' (mean+/-SE usage).
    cohort_samples : keep only these sample tokens (prefix match); None = all samples.
    lineage_map    : {cell_state: lineage}; required for level='lineage'.
    cluster_order  : ordered x-axis groups (cell states or lineages). None = observed order.
    clone_map      : {isoform_id: label} -> legend shows '<label>: <isoform>'.
    colors         : {isoform_id: hex}; None -> distinct PALETTE.
    xlabel_fontsize: x-tick label size; None -> 2x default (12 for cellstate, 23 for lineage).
    """
    if isinstance(pseudobulk, str):
        pseudobulk = load_pseudobulk_slice(pseudobulk, ensg=ensg)
    a = pseudobulk
    obs = np.asarray(a.obs_names, dtype=str)
    vn = np.asarray(a.var_names, dtype=str)
    gof = np.array([v.split(":", 1)[0] for v in vn])
    sm = np.array([_sample_of(o) for o in obs])
    cs = np.array([_cellstate_of(o) for o in obs])
    keep = np.array([_in_cohort(s, cohort_samples) for s in sm])
    if level == "lineage":
        if not lineage_map:
            raise ValueError("level='lineage' requires lineage_map {cell_state: lineage}")
        grp = np.array([lineage_map.get(c, "") for c in cs])
    else:
        grp = cs
    grpk = grp[keep]; smk = sm[keep]
    GRPORDER = list(cluster_order) if cluster_order else sorted(set(g for g in grpk if g))

    cols = np.where(gof == ensg)[0]
    if len(cols) < 2:
        print(f"  ! {gene}: <2 isoforms for {ensg}; skip"); return None
    X = a.X[:, cols]
    X = (np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X))[keep]
    short = np.array([vn[c].split(":", 1)[1] for c in cols])
    order = []
    for f in isoforms:
        w = np.where(short == f)[0]
        if not len(w):
            print(f"  ! {gene}: isoform {f} absent in data; skip gene"); return None
        order.append(int(w[0]))
    isos = [short[i] for i in order]
    colmap = {short[i]: (colors.get(short[i]) if colors else PALETTE[k % len(PALETTE)])
              for k, i in enumerate(order)}

    grps = [L for L in GRPORDER if (grpk == L).sum() > 0]
    data = {i: {L: [] for L in grps} for i in order}
    for L in grps:
        for s in sorted(set(smk[grpk == L])):
            msk = (grpk == L) & (smk == s); v = X[msk].sum(0); t = v.sum()
            if t < min_group_reads:
                continue
            for i in order:
                data[i][L].append(100 * v[i] / t)
    keepL = [L for L in grps if len(data[order[0]][L]) >= min_donors]
    if not keepL:
        print(f"  ! {gene}: no group with >={min_donors} donors; skip"); return None

    R = struct_renderer or StructRenderer(gene_model=gene_model, gff=gff_dir, intron_scale=intron_scale)
    R.load_genes([ensg]); xlim = R.shared_xlim([(ensg, i) for i in isos])

    _fw = (0.95 if level == "lineage" else 0.42)
    nstruct = len(isos)
    _spacer = (1.6 if level == "cellstate" else 1.0)   # empty row so the enlarged rotated x-labels sit ABOVE tracks
    xfs = xlabel_fontsize if xlabel_fontsize is not None else (23 if level == "lineage" else 12)
    fig = plt.figure(figsize=(max(8, _fw * len(keepL) + 3), 3.6 + 0.46 * nstruct + _spacer))
    gs = GridSpec(2 + nstruct, 1, height_ratios=[3.4, _spacer] + [0.5] * nstruct, hspace=0.45)
    ax = fig.add_subplot(gs[0]); xpos = np.arange(len(keepL))

    if plottype == "bar":
        muM = np.array([[np.mean(data[i][L]) if data[i][L] else 0.0 for L in keepL] for i in order])
        colsum = muM.sum(axis=0); colsum[colsum == 0] = 1.0
        frac = 100 * muM / colsum
        bottom = np.zeros(len(keepL)); bw = 0.8 if level == "lineage" else 0.9
        for k, i in enumerate(order):
            lbl = (clone_map or {}).get(short[i], "")
            leg = (lbl + ": " if lbl else "") + short[i]
            ax.bar(xpos, frac[k], bottom=bottom, width=bw, color=colmap[short[i]], edgecolor="#ffffff",
                   linewidth=0.4, label=leg, zorder=3)
            bottom = bottom + frac[k]
        ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(keepL) - 0.4); ax.patch.set_visible(False)
        ax.set_ylabel("% composition among shown isoforms", fontsize=12)
    else:
        for i in order:
            mu = [np.mean(data[i][L]) if data[i][L] else np.nan for L in keepL]
            se = [np.std(data[i][L], ddof=1) / np.sqrt(len(data[i][L])) if len(data[i][L]) > 1 else 0 for L in keepL]
            lbl = (clone_map or {}).get(short[i], "")
            leg = (lbl + ": " if lbl else "") + short[i]
            ax.errorbar(xpos, mu, yerr=se, marker="o", ms=7.0, color=colmap[short[i]], lw=2.4, capsize=4.5,
                        capthick=1.2, elinewidth=1.0, label=leg, zorder=3, clip_on=False)
        ax.set_ylim(0, 100); ax.set_xlim(-0.4, len(keepL) - 0.6); ax.patch.set_visible(False)
        ax.set_ylabel("% isoform usage (mean ± SE)", fontsize=12)

    ax.set_xticks(xpos)
    ax.set_xticklabels(keepL, rotation=40, ha="right", fontsize=xfs, fontname="Arial")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(fontsize=8.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.34), ncol=min(3, nstruct))

    for k, i in enumerate(isos):                            # structure tracks BELOW the spacer (gs index 2+k)
        axs = fig.add_subplot(gs[2 + k]); R.draw(axs, ensg, i, colmap[i], xlim=xlim)

    fig.suptitle(title or f"{gene} ({ensg}) isoform usage by {level} + structures", fontsize=13, y=0.997)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path + ".pdf", bbox_inches="tight")
    fig.savefig(out_path + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  DONE {gene}: {os.path.basename(out_path)}.pdf (+.png)  groups={len(keepL)} isoforms={len(isos)}")
    return out_path + ".pdf"
