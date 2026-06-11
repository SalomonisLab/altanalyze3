#!/usr/bin/env python3
"""SATAY-UDON heatmaps: rows = UDON cluster x its cell types, columns = AML covariates
(demographic / clinical / mutation), color = enrichment p-value (purple-on-black,
matching the SATAY-UDON format). SATAY test: WITHIN each cell type, is the cluster's
slice enriched for the covariate (Fisher's exact, one-sided; background = all pseudobulks
of that cell type)? -- so cell-of-origin is controlled for.

Pure Python (pandas/scipy/matplotlib). No R. Editable PDF fonts (Arial).
"""
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

UDON_DIR = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, UDON_DIR)
import satay_metadata_enrichment as SME
from satay_metadata_enrichment import build_features, load_run, unit_fisher, _units, MIN_UNITS

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42; plt.rcParams['ps.fonttype'] = 42

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
OUT = os.path.join(PB, "UDON", "satay_metadata"); os.makedirs(OUT, exist_ok=True)
MIN_PB = 15; TOP_CT = 5

GROUP_COLORS = {"Demographic": "#3B6FB0", "Clinical": "#2E8B57", "Mutation": "#8B3A3A"}
_GORDER = {"Demographic": 0, "Clinical": 1, "Mutation": 2}


def auto_covars(feat, defined, min_pos=1):
    """Assemble the covariate list AUTOMATICALLY from every metadata feature
    (all mutations + all clinical levels + demographic) with >= min_pos positive
    samples. No hand selection / no AI curation -- the columns shown are then
    filtered data-drivenly to those actually enriched (see plot)."""
    rows = []
    for key, smap in feat.items():
        pos = [s for s, v in smap.items() if v]
        if len(set(_units(pos))) < min_pos:          # >= min_pos distinct positive donors/samples
            continue
        if key[0] == "mutation":
            grp, disp = "Mutation", key[1]
        elif key[0] == "clinical":
            grp, disp = "Clinical", key[1].replace("disease_category=", "").replace("timepoint=", "").replace("=", ":")
        else:
            grp, disp = "Demographic", ("Age>=60" if str(key[0]).startswith("age") else str(key[0]))
        rows.append((key, disp, grp))
    rows.sort(key=lambda r: (_GORDER.get(r[2], 3), r[1]))
    return rows


def add_age(feat, defined):
    c = pd.read_excel(os.path.join(PB, "UDON", "AML_harmonized_metadata.xlsx"), sheet_name="Clinical_Metadata")
    c = c.set_index(c["Sample"].astype(str))
    age = pd.to_numeric(c["age"], errors="coerce").groupby(level=0).first()   # dedupe dup samples
    dfn = set(age.dropna().index)
    feat[("age>=60",)] = {s: int(age[s] >= 60) for s in dfn}
    defined[("age>=60",)] = dfn


def top_marker(run_dir, cluster):
    for fn in ("marker_genes_all.txt", "marker_genes.txt", "final_program_markers.txt"):
        p = os.path.join(run_dir, fn)
        if not os.path.exists(p): continue
        try:
            m = pd.read_csv(p, sep="\t")
            g = m[m["top_cluster"].astype(str) == str(cluster)].sort_values("pearson_r", ascending=False)
            if len(g): return ",".join(g["marker"].astype(str).head(2).tolist())
        except Exception:
            pass
    return ""


def build_matrix(cl, feat, defined, covars, min_pb=MIN_PB):
    clusters = sorted(cl["cluster"].unique(), key=lambda x: (len(x), x))
    rows = []
    for c in clusters:
        sub = cl[cl["cluster"] == c]
        vc = sub["celltype"].value_counts()
        cts = [ct for ct in vc.index if vc[ct] >= min_pb][:TOP_CT]
        for ct in cts:
            rows.append((c, ct))
    M = np.ones((len(rows), len(covars)))
    for i, (c, ct) in enumerate(rows):
        T = cl[cl["celltype"] == ct]
        for j, (key, disp, grp) in enumerate(covars):
            if key not in feat: continue
            smap = feat[key]; dfn = defined[key]
            Tg = T[T["Sample"].isin(dfn)]
            if len(Tg) < 8: continue
            fpos = Tg["Sample"].map(smap).fillna(0).astype(int).values
            p, _orr, a, _cn, _fn = unit_fisher(Tg["Sample"].values, (Tg["cluster"] == c).values, fpos)
            M[i, j] = p            # p=1 automatically when < MIN_UNITS distinct donors back it
    return rows, M


def compute_pmap():
    """P# (study-aware final program) -> dominant merged matched_sex_platform U#."""
    try:
        fa = pd.read_csv(os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv"),
                         sep="\t").set_index("pseudobulk")
        mu = pd.read_csv(os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt"),
                         sep="\t", index_col=0); mu.columns = ["U"]
        ct = pd.crosstab(fa.join(mu, how="inner")["final_program"], mu["U"])
        return {p: ct.loc[p].idxmax() for p in ct.index}
    except Exception:
        return {}


def plot(rows, M, covars, run_dir, title, outbase, pmap=None):
    # data-driven: keep only covariate columns that are enriched (p<0.1) in >=1 row
    keepc = [j for j in range(len(covars)) if (M[:, j] < 0.1).any()]
    if not rows or not keepc:
        print(f"  {os.path.basename(outbase)}: no enriched covariates, skipped"); return
    M = M[:, keepc]; covars = [covars[j] for j in keepc]
    nlp = -np.log10(np.clip(M, 1e-12, 1))
    cmap = LinearSegmentedColormap.from_list("satay_purple",
            ["#E4D6EE", "#C9A8DE", "#A86BC9", "#7E3CB0", "#4B1A86"])      # light -> dark violet
    cmap.set_under("#000000")                                            # p>=0.1 -> black (n.s.)
    nrow, ncol = len(rows), len(covars)
    fig, ax = plt.subplots(figsize=(0.46 * ncol + 8.0, 0.34 * nrow + 3.0))
    im = ax.imshow(nlp, aspect="auto", cmap=cmap, vmin=1.0, vmax=4.0, interpolation="nearest")
    im.set_rasterized(True)                                             # embed raster in PDF (fixes vector render)

    clusters = [r[0] for r in rows]
    ax.set_yticks(range(nrow)); ax.set_yticklabels([r[1] for r in rows], fontsize=10)
    ax.set_xticks(range(ncol)); ax.set_xticklabels([c[1] for c in covars], rotation=90, fontsize=10)
    bounds = []; start = 0
    for k in range(1, nrow + 1):
        if k == nrow or clusters[k] != clusters[start]:
            bounds.append((start, k - 1))
            if k < nrow: ax.axhline(k - 0.5, color="white", lw=2)
            start = k
    for (s, e) in bounds:
        c = clusters[s]; mk = top_marker(run_dir, c)
        lab = f"{c}:{pmap[c]}" if (pmap and c in pmap) else f"{c}"      # P#:U# mapping when available
        ax.text(-0.06 * ncol - 3.6, (s + e) / 2, lab, ha="right", va="center",
                fontsize=12, fontweight="bold", color="#1A3A6B")
        if mk:   # enriched gene set listed on the RIGHT of the heatmap
            ax.text(ncol - 0.3, (s + e) / 2, mk, ha="left", va="center",
                    fontsize=9, color="#3B6FB0", style="italic")
    for grp, col in GROUP_COLORS.items():
        idx = [j for j, c in enumerate(covars) if c[2] == grp]
        if idx:
            ax.plot([min(idx) - 0.4, max(idx) + 0.4], [-1.1, -1.1], color=col, lw=4, clip_on=False)
            ax.text((min(idx) + max(idx)) / 2, -1.7, grp, ha="center", va="bottom",
                    fontsize=11, color=col, fontweight="bold")
    ax.set_xlim(-0.5, ncol - 0.5); ax.set_ylim(nrow - 0.5, -2.2)
    ax.set_title(title, fontsize=15, pad=28)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)

    leg = [Patch(fc="#4B1A86", label="<0.0001"), Patch(fc="#7E3CB0", label="0.001"),
           Patch(fc="#A86BC9", label="0.01"), Patch(fc="#E4D6EE", label="0.1"),
           Patch(fc="#000000", label="n.s.")]
    fig.legend(handles=leg, title="P-value", loc="lower center", bbox_to_anchor=(0.5, -0.015),
               ncol=5, fontsize=11, title_fontsize=11, frameon=False)
    fig.savefig(outbase + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(outbase + ".pdf", dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {os.path.basename(outbase)}.png/.pdf  ({nrow} rows x {ncol} covariates)")


def load_final(path):
    """final_program_assignments.tsv -> cl frame (cluster=final_program, Sample, celltype)."""
    fa = pd.read_csv(path, sep="\t")
    cl = pd.DataFrame(index=fa["pseudobulk"])
    cl["cluster"] = fa["final_program"].astype(str).values
    cl["Sample"] = fa["Sample"].astype(str).values
    cl["celltype"] = [str(p).split("__", 1)[0] for p in fa["pseudobulk"]]   # celltype before first "__"
    return cl


def generate_all(feat=None, defined=None):
    """Default SATAY-UDON heatmap export: one heatmap per UDON result."""
    if feat is None:
        feat, defined = build_features()
    if ("age>=60",) not in feat:
        add_age(feat, defined)
    jobs = []   # (name, cl, marker_run_dir)
    # merged runs
    for name, rdir in [("matched_sex_platform", os.path.join(PB, "UDON/matched_sex_platform/udon_core")),
                       ("udon_core", os.path.join(PB, "UDON/udon_core"))]:
        f = os.path.join(rdir, "udon_clusters.txt")
        if os.path.exists(f): jobs.append((name, load_run(f), rdir))
    # per-study runs
    import glob
    for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt"))):
        st = os.path.basename(os.path.dirname(d))
        jobs.append((f"per_study_{st}", load_run(d), os.path.dirname(d)))
    # final study-aware programs
    fa = os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv")
    if os.path.exists(fa):
        jobs.append(("study_aware_final", load_final(fa), os.path.join(PB, "UDON/study_aware")))

    PMAP = compute_pmap()                                            # P# -> merged U# mapping
    COV = auto_covars(feat, defined)                                 # automatic: every feature, no curation
    for name, cl, rdir in jobs:
        try:
            # per-study clusters are small; step the cell-type floor down until rows exist
            thresholds = [6, 4, 3, 2] if name.startswith("per_study") else [MIN_PB]
            rows, M, used = [], None, None
            for mpb in thresholds:
                rows, M = build_matrix(cl, feat, defined, COV, min_pb=mpb); used = mpb
                if rows: break
            if not rows:
                print(f"  {name}: no cell-type rows even at floor=2, skipped"); continue
            note = f" (cell-type floor={used} pb)" if used != MIN_PB else ""
            plot(rows, M, COV, rdir, f"SATAY-UDON: {name}{note}", os.path.join(OUT, f"satay_heatmap_{name}"),
                 pmap=PMAP if name == "study_aware_final" else None)
        except Exception as e:
            print(f"  {name}: FAILED {e}")


if __name__ == "__main__":
    generate_all()
