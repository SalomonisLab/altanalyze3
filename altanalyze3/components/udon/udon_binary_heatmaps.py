#!/usr/bin/env python3
"""Standard UDON binary composition heatmaps (produced for every UDON result):

  donor_cluster_heatmap    rows = donors,     ordered by DOMINANT cluster (left bar = study)
  celltype_cluster_heatmap rows = cell types, clustered within lineage (left bar = lineage)

ONE column per PSEUDOBULK, ordered by cluster exactly like the MarkerFinder heatmap (top bar =
the same rainbow cluster palette). Cell = 1 (black) / 0 (white): does that pseudobulk belong to
the row (donor / cell type). So each column has a single black cell, and each cluster's column
block shows its donor / cell-type composition -- i.e. whether a cluster's pseudobulks come from
many donors (shared) or few (donor-specific / batch-like).

donor_specificity_per_cluster.tsv quantifies this: per cluster, n_donors, top-donor fraction,
and a specificity call (donor-specific if top donor >=50%, few-donor if <=3 donors, else shared).
Donor rows are NOT hierarchically clustered (one-hot rows are orthogonal -> meaningless); they are
ordered by their dominant cluster. Figure height grows with row count so labels stay legible.

Pure Python. Editable TrueType fonts; Adobe-safe PDFs (compression off + rasterized imshow)."""
import os, re, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

try:
    from visualizations import assign_rainbow_colors_to_groups
except Exception:
    assign_rainbow_colors_to_groups = None

_QUAL = ["#3B6FB0", "#E07B39", "#2E8B57", "#B0413E", "#7E3CB0", "#C9A227",
         "#1B9E9E", "#8B5A2B", "#D4549A", "#5B5B5B", "#6AA84F", "#A64D79"]


def _natnum(x):
    m = re.search(r"\d+", str(x)); return int(m.group()) if m else 0


def _lineage(ct):
    c = str(ct).lower()
    if any(k in c for k in ["hsc", "mpp", "lmpp", "multilin", "cmp", "gmp", "mdp", "bmcp", "cmop", "mkp", "hspc", "prog"]):
        return "Progenitor"
    if any(k in c for k in ["eryth", "erp", "mep", "megak", "platelet", "rbc"]):
        return "Erythroid/Mega"
    if any(k in c for k in ["mono", "dc", "neu", "mast", "macro", "myeloid", "baso", "eos", "asdc"]):
        return "Myeloid"
    if any(k in c for k in ["cd4", "cd8", "nk", "mait", "ilc", "plasma", "tcm", "tem", "treg", "naive", "pro-b", "t ", "b "]):
        return "Lymphoid"
    return "Other"


def _cluster_colors(clusters):
    if assign_rainbow_colors_to_groups is not None:
        try:
            return assign_rainbow_colors_to_groups(groups=np.array(clusters))
        except Exception:
            pass
    import matplotlib.cm as cm
    return {c: cm.rainbow(i / max(len(clusters) - 1, 1)) for i, c in enumerate(clusters)}


def _order_within_groups(binary, groups):
    order = []
    for g in sorted(set(groups.values())):
        rows = [r for r in binary.index if groups[r] == g]
        sub = binary.loc[rows]
        if len(rows) > 2 and sub.values.sum() > 0:
            try:
                d = np.nan_to_num(pdist(sub.values.astype(float), metric="jaccard"), nan=1.0)
                rows = [sub.index[i] for i in leaves_list(linkage(d, method="average"))]
            except Exception:
                pass
        order += list(rows)
    return order


def _plot(M, col_clusters, groups, agg, outbase, title, row_fontsize=6, row_order=None):
    """M: rows (donor/celltype) x PSEUDOBULK columns, one-hot (1 = that pseudobulk belongs to the row).
    col_clusters: cluster per pseudobulk column (aligned to M.columns). row_order: explicit row order
    (e.g. donors by dominant cluster); else hierarchical clustering of `agg` within each left-bar group."""
    seq = np.asarray(col_clusters)
    uniq = sorted(set(seq), key=_natnum)
    ccolors = _cluster_colors(uniq)
    gcolors = {g: _QUAL[i % len(_QUAL)] for i, g in enumerate(sorted(set(groups.values())))}
    M = M.loc[row_order if row_order is not None else _order_within_groups(agg, groups)]
    nrow, npb = M.shape
    fig = plt.figure(figsize=(13.0, max(4.0, 0.12 * nrow + 2.0)))
    gs = GridSpec(2, 2, width_ratios=(0.5, 40), height_ratios=(1, 22), wspace=0.01, hspace=0.02)
    axt = fig.add_subplot(gs[0, 1]); axl = fig.add_subplot(gs[1, 0]); ax = fig.add_subplot(gs[1, 1])

    im = ax.imshow(M.values, cmap=ListedColormap(["#FFFFFF", "#000000"]), vmin=0, vmax=1,
                   aspect="auto", interpolation="nearest"); im.set_rasterized(True)
    ax.set_xticks([]); ax.set_yticks(range(nrow)); ax.yaxis.tick_right()
    ax.set_yticklabels(M.index, fontsize=row_fontsize); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)

    # top cluster color bar (one cell per pseudobulk) + P# labels at block centers
    axt.imshow(np.array([[to_rgb(ccolors[c]) for c in seq]]), aspect="auto")
    axt.set_yticks([]); axt.set_xticks([]); axt.set_xlim(-0.5, npb - 0.5)
    for sp in axt.spines.values(): sp.set_visible(False)
    s = 0
    for i in range(1, npb + 1):
        if i == npb or seq[i] != seq[s]:
            axt.text((s + i - 1) / 2.0, -0.7, str(seq[s]), ha="center", va="bottom", fontsize=7, rotation=90)
            s = i

    axl.imshow(np.array([[to_rgb(gcolors[groups[r]])] for r in M.index]), aspect="auto")
    axl.set_xticks([]); axl.set_yticks([])
    for sp in axl.spines.values(): sp.set_visible(False)
    gseq = [groups[r] for r in M.index]; s = 0
    for i in range(1, nrow + 1):
        if i == nrow or gseq[i] != gseq[s]:
            axl.text(-0.8, (s + i - 1) / 2.0, gseq[s], transform=axl.get_yaxis_transform(),
                     ha="right", va="center", fontsize=8, color=gcolors[gseq[s]], fontweight="bold")
            s = i

    ax.set_title(title, fontsize=13, pad=24)
    fig.legend(handles=[Patch(fc="#000000", label="1 (pseudobulk of this row)"),
                        Patch(fc="#FFFFFF", ec="#888888", label="0")],
               loc="lower left", bbox_to_anchor=(0.0, -0.01), ncol=2, fontsize=9, frameon=False)
    _prev = matplotlib.rcParams['pdf.compression']; matplotlib.rcParams['pdf.compression'] = 0
    fig.savefig(outbase + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(outbase + ".pdf", bbox_inches="tight")
    matplotlib.rcParams['pdf.compression'] = _prev; plt.close(fig)
    print(f"  wrote {os.path.basename(outbase)}  ({nrow} rows x {npb} pseudobulks)")


def make_udon_binary_heatmaps(clusters_df, outdir, donor_of=None, study_of=None, cluster_col="cluster"):
    """clusters_df: index = 'celltype__Sample', column `cluster_col`. ONE column per pseudobulk
    (ordered by cluster, like the marker heatmap). Writes celltype_cluster_heatmap (always) and
    donor_cluster_heatmap (if donor_of+study_of given)."""
    os.makedirs(outdir, exist_ok=True)
    df = pd.DataFrame({"cluster": clusters_df[cluster_col].astype(str).values}, index=clusters_df.index)
    df["Sample"] = [str(i).split("__", 1)[1] if "__" in str(i) else str(i) for i in df.index]
    df["celltype"] = [str(i).split("__", 1)[0] for i in df.index]
    # pseudobulk order = by cluster (numeric, stable within cluster), as in the marker heatmap
    df = df.iloc[np.argsort([_natnum(c) for c in df["cluster"]], kind="stable")]
    col_clusters = df["cluster"].values

    ct = pd.get_dummies(df["celltype"]).T                       # celltype x pseudobulk (one-hot)
    ct_agg = (pd.crosstab(df["celltype"], df["cluster"]) > 0).astype(int)
    _plot(ct.loc[:, df.index], col_clusters, {c: _lineage(c) for c in ct.index}, ct_agg,
          os.path.join(outdir, "celltype_cluster_heatmap"), "UDON cell-type composition (one column per pseudobulk)")

    if donor_of is not None and study_of is not None:
        df["Donor"] = df["Sample"].map(lambda s: donor_of.get(str(s), str(s)))
        cnt = pd.crosstab(df["Donor"], df["cluster"])           # donor x cluster pseudobulk COUNTS
        clusters = sorted(cnt.columns, key=_natnum); cnt = cnt[clusters]
        dn = pd.get_dummies(df["Donor"]).T                      # donor x pseudobulk (one-hot)
        d2s = df.drop_duplicates("Donor").set_index("Donor")["Sample"].map(lambda s: study_of.get(str(s), "NA"))
        groups = {d: d2s.get(d, "NA") for d in cnt.index}

        # SPECIFICITY of donor pseudobulks per cluster: is a cluster driven by few donors?
        tot = cnt.sum(0); top = cnt.max(0)
        spec = pd.DataFrame({"n_pseudobulks": tot.astype(int),
                             "n_donors": (cnt > 0).sum(0).astype(int),
                             "n_studies": df.groupby("cluster")["Sample"].apply(
                                 lambda s: pd.Series([study_of.get(str(x), "NA") for x in s]).nunique()).reindex(clusters),
                             "top_donor": cnt.idxmax(0),
                             "top_donor_frac": (top / tot).round(3),
                             "specificity": np.where(top / tot >= 0.5, "donor-specific",
                                            np.where((cnt > 0).sum(0) <= 3, "few-donor", "shared"))})
        spec.index.name = "cluster"
        spec.to_csv(os.path.join(outdir, "donor_specificity_per_cluster.tsv"), sep="\t")

        # order donors by their dominant cluster (reveals which clusters are donor-specific), then study
        dom = cnt.idxmax(axis=1)
        row_order = sorted(cnt.index, key=lambda d: (_natnum(dom[d]), groups[d], str(d)))
        _plot(dn.loc[:, df.index], col_clusters, groups, None,
              os.path.join(outdir, "donor_cluster_heatmap"),
              "UDON donor pseudobulks per cluster (donors ordered by dominant cluster)", row_order=row_order)


def _mappings(PB):
    donor_of = study_of = None
    try:
        clin = pd.read_excel(os.path.join(PB, "UDON", "AML_harmonized_metadata.xlsx"), sheet_name="Clinical_Metadata")
        donor_of = dict(zip(clin["Sample"].astype(str), clin["Donor_ID"].astype(str)))
    except Exception:
        pass
    try:
        sm = pd.read_csv(os.path.join(PB, "evaluation/outputs/sample_metadata.tsv"), sep="\t", index_col=0)
        study_of = dict(zip(sm.index.astype(str), sm["Dataset"].astype(str)))
    except Exception:
        pass
    return donor_of, study_of


def generate_all(PB="/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"):
    """Standard binary heatmaps for every UDON result, written next to each marker_heatmap."""
    import glob
    donor_of, study_of = _mappings(PB)

    def load_clusters(path):
        c = pd.read_csv(path, sep="\t", index_col=0); c.columns = ["cluster"]; return c

    jobs = []
    fa = os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv")
    if os.path.exists(fa):
        f = pd.read_csv(fa, sep="\t")
        jobs.append((pd.DataFrame({"cluster": f["final_program"].values}, index=f["pseudobulk"]),
                     os.path.join(PB, "UDON/study_aware")))
    for cl, rd in [("UDON/matched_sex_platform/udon_core/udon_clusters.txt", "UDON/matched_sex_platform"),
                   ("UDON/udon_core/udon_clusters.txt", "UDON")]:
        p = os.path.join(PB, cl)
        if os.path.exists(p): jobs.append((load_clusters(p), os.path.join(PB, rd)))
    for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt"))):
        jobs.append((load_clusters(d), os.path.dirname(d)))

    for cdf, outdir in jobs:
        try:
            make_udon_binary_heatmaps(cdf, outdir, donor_of=donor_of, study_of=study_of)
        except Exception as e:
            print(f"  {outdir}: FAILED {e}")


if __name__ == "__main__":
    generate_all()
