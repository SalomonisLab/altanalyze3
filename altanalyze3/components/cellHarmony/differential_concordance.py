"""Concordance of differential effects between populations, from a cellHarmony fixed-order heatmap.

Inputs are the two tables build_fixed_order_heatmap and the differential write: the heatmap TSV (rows
'<group>:<feature>', columns '<population>:<population>', values = log2 fold or dPSI, blank where not tested) and
DEG_detailed (gene, population, log2fc, sig_metric). For every pair of populations with >= --min_hits hits:
  shared        hits in both
  jaccard       shared / hits in either
  same_sign     fraction of shared hits that change in the same direction
  pearson_r     Pearson r of the values over the hits of either population that are tested in both
                (>= --min_pairs features, else missing)
Populations are clustered by average linkage on 1 - r (missing r counts as 0) and cut into k groups (--k, several
allowed). For each k, each group's hits are the features that are hits in >= 1 of its populations; a hit is unique
to a group when it is a hit in no population of another group, and 'unique, tested elsewhere' when it is also tested
in >= 1 population of every other group (so its absence there is observed, not missing). Gene lists of unique hits
(symbol = first token of the feature label) are written for enrichment, and GO-Elite runs on them with --goelite.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd


def load_heatmap(path):
    h = pd.read_csv(path, sep="\t", index_col=0)
    h.index = [str(i).split(":", 1)[1] if ":" in str(i) else str(i) for i in h.index]
    h.columns = [str(c).split(":", 1)[0] for c in h.columns]
    if not h.index.is_unique:
        raise ValueError(f"{path}: duplicate feature labels after removing the group prefix")
    return h.astype(float)


def load_hits(deg_detailed, features, populations):
    d = pd.read_csv(deg_detailed, sep="\t")
    hits = pd.DataFrame(False, index=features, columns=populations)
    known = d[d.gene.isin(hits.index) & d.population.astype(str).isin(populations)]
    missing = int((~d.gene.isin(hits.index)).sum())
    for pop, g in known.groupby(known.population.astype(str)):
        hits.loc[g.gene, pop] = True
    return hits, missing


def pairwise(values, hits, min_hits=100, min_pairs=30):
    pops = [p for p in hits.columns if hits[p].sum() >= min_hits]
    rows = []
    R = pd.DataFrame(np.nan, index=pops, columns=pops)
    for i, a in enumerate(pops):
        R.loc[a, a] = 1.0
        for b in pops[i + 1:]:
            ha, hb = hits[a], hits[b]
            both = ha & hb
            either = ha | hb
            va, vb = values[a], values[b]
            tested = either & va.notna() & vb.notna()
            r = float(np.corrcoef(va[tested], vb[tested])[0, 1]) if tested.sum() >= min_pairs else np.nan
            same = (np.sign(va[both]) == np.sign(vb[both])).mean() if both.sum() else np.nan
            R.loc[a, b] = R.loc[b, a] = r
            rows.append(dict(population_a=a, population_b=b, hits_a=int(ha.sum()), hits_b=int(hb.sum()),
                             shared=int(both.sum()), jaccard=float(both.sum() / either.sum()) if either.sum() else np.nan,
                             same_sign=float(same), pearson_r=r, features_for_r=int(tested.sum())))
    return R, pd.DataFrame(rows), pops


def cluster(R, ks, method="average"):
    """average: average linkage on 1 - r. profile: average linkage on the correlation distance between populations'
    rows of r (populations whose concordance with every other population is alike group together)."""
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
    from scipy.spatial.distance import pdist, squareform
    D = 1 - R.fillna(0.0).to_numpy()
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    if method == "profile":
        P = R.fillna(float(np.nanmedian(R.to_numpy()))).to_numpy()
        Z = linkage(pdist(P, "correlation"), "average")
    else:
        Z = linkage(squareform(np.clip(D, 0, None), checks=False), "average")
    order = [R.index[i] for i in leaves_list(Z)]
    groups = {k: pd.Series(fcluster(Z, k, "maxclust"), index=R.index) for k in ks}
    # silhouette on the same distances
    sil = {}
    for k, g in groups.items():
        if g.nunique() < 2:
            sil[k] = np.nan
            continue
        s = []
        for i, p in enumerate(R.index):
            same = [j for j, q in enumerate(R.index) if g[q] == g[p] and j != i]
            if not same:
                s.append(0.0)
                continue
            a = D[i, same].mean()
            b = min(D[i, [j for j, q in enumerate(R.index) if g[q] == c]].mean() for c in set(g) if c != g[p])
            s.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
        sil[k] = float(np.mean(s))
    return order, groups, sil


def group_unique(values, hits, groups):
    """Per group: hits, unique hits, unique hits tested in every other group; and per-feature membership."""
    names = sorted(set(groups))
    in_group = pd.DataFrame({g: hits[groups.index[groups == g]].any(axis=1) for g in names})
    tested_group = pd.DataFrame({g: values[groups.index[groups == g]].notna().any(axis=1) for g in names})
    in_group = in_group[in_group.any(axis=1)]
    tested_group = tested_group.loc[in_group.index]
    summary, membership = [], in_group.copy()
    for g in names:
        others = [o for o in names if o != g]
        uniq = in_group[g] & ~in_group[others].any(axis=1)
        observed = uniq & tested_group[others].all(axis=1)
        summary.append(dict(group=g, populations=";".join(groups.index[groups == g]), hits=int(in_group[g].sum()),
                            unique=int(uniq.sum()), unique_tested_in_every_other_group=int(observed.sum())))
        membership[f"unique_{g}"] = uniq
        membership[f"unique_tested_elsewhere_{g}"] = observed
    regions = in_group.apply(lambda r: "&".join(str(g) for g in names if r[g]), axis=1).value_counts()
    return pd.DataFrame(summary), membership, regions


def symbol(label):
    return str(label).split(" ", 1)[0]


def plot_correlation(R, order, groups, outfile, value_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    M = R.loc[order, order].to_numpy()
    n = len(order)
    fig, ax = plt.subplots(figsize=(0.12 * n + 2.2, 0.12 * n + 1.6))
    im = ax.pcolormesh(np.arange(n + 1), np.arange(n + 1), M, cmap=plt.get_cmap("Spectral_r"), vmin=-1, vmax=1,
                       edgecolors="none")
    ax.set_xlim(0, n); ax.set_ylim(n, 0)
    ax.set_xticks(np.arange(n) + 0.5); ax.set_xticklabels(order, rotation=90, fontsize=5)
    ax.set_yticks(np.arange(n) + 0.5); ax.set_yticklabels(order, fontsize=5)
    g = groups.loc[order].to_numpy()
    for i in range(1, n):
        if g[i] != g[i - 1]:
            ax.plot([i, i], [0, n], color="#000000", lw=0.8)
            ax.plot([0, n], [i, i], color="#000000", lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(f"Pearson r of {value_label}", fontsize=6); cb.ax.tick_params(labelsize=5)
    fig.savefig(outfile, bbox_inches="tight"); plt.close(fig)


def plot_venn(summary, regions, outfile, labels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    names = [str(g) for g in summary.group]
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    if len(names) in (2, 3):
        from matplotlib_venn import venn2, venn3
        keys = ["10", "01", "11"] if len(names) == 2 else ["100", "010", "110", "001", "101", "011", "111"]
        sub = {}
        for key in keys:
            sel = "&".join(n for n, bit in zip(names, key) if bit == "1")
            sub[key] = int(regions.get(sel, 0))
        (venn2 if len(names) == 2 else venn3)(subsets=sub, set_labels=[labels[n] for n in names], ax=ax,
                                              set_colors=["#4884BC", "#E9B940", "#59A14F"][:len(names)])
        for t in ax.texts:
            t.set_fontsize(6)
    else:
        r = regions.sort_values(ascending=True)
        ax.barh(np.arange(len(r)), r.to_numpy(), color="#7F7F7F", height=0.7)
        ax.set_yticks(np.arange(len(r))); ax.set_yticklabels([" & ".join(labels[x] for x in k.split("&")) for k in r.index], fontsize=5)
        ax.set_xlabel("Hits", fontsize=6); ax.tick_params(axis="x", labelsize=5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.savefig(outfile, bbox_inches="tight"); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--heatmap_tsv", required=True)
    ap.add_argument("--deg_detailed", required=True)
    ap.add_argument("--k", type=int, nargs="+", default=[2, 3, 4], help="numbers of groups to cut the tree into")
    ap.add_argument("--min_hits", type=int, default=100, help="populations with fewer hits are left out")
    ap.add_argument("--min_pairs", type=int, default=30, help="features needed for a Pearson r")
    ap.add_argument("--method", choices=["average", "profile"], default="average", help="clustering (see cluster)")
    ap.add_argument("--groups_file", default=None,
                    help="TSV of population and group (named groups) used instead of clustering; written as k = number "
                         "of groups")
    ap.add_argument("--value_label", default="dPSI (aged - young)")
    ap.add_argument("--goelite_species", default=None, help="run GO-Elite on each group's unique genes (e.g. human)")
    ap.add_argument("--goelite_background", default=None, help="newline gene list for the GO-Elite background")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    values = load_heatmap(a.heatmap_tsv)
    hits, missing = load_hits(a.deg_detailed, values.index, list(values.columns))
    R, pairs, pops = pairwise(values, hits, a.min_hits, a.min_pairs)
    R.to_csv(os.path.join(a.outdir, "concordance_pearson_r.tsv"), sep="\t")
    pairs.to_csv(os.path.join(a.outdir, "concordance_pairs.tsv"), sep="\t", index=False)
    order, groups, sil = cluster(R, a.k, method=a.method)
    if a.groups_file:
        gf = pd.read_csv(a.groups_file, sep="\t", dtype=str)
        gmap = dict(zip(gf.iloc[:, 0], gf.iloc[:, 1]))
        keep = [p for p in pops if p in gmap]
        groups = {len(set(gmap[p] for p in keep)): pd.Series([gmap[p] for p in keep], index=keep)}
        report_left = [p for p in pops if p not in gmap]
        R, pops = R.loc[keep, keep], keep
        order = [p for p in order if p in keep]
        sil = {}
    else:
        report_left = []
    report = dict(heatmap_tsv=os.path.abspath(a.heatmap_tsv), deg_detailed=os.path.abspath(a.deg_detailed),
                  features=int(len(values)), hits_not_in_heatmap=missing,
                  populations_used=pops, populations_left_out={p: int(hits[p].sum()) for p in hits.columns if p not in pops},
                  silhouette=sil, method=a.method, groups_file=a.groups_file, populations_without_group=report_left,
                  groups={})
    for k, g in groups.items():
        d = os.path.join(a.outdir, f"k{k}"); os.makedirs(d, exist_ok=True)
        g.rename("group").rename_axis("population").to_csv(os.path.join(d, "population_groups.tsv"), sep="\t")
        summ, memb, regions = group_unique(values[pops], hits[pops], g)
        labels = {str(x): f"G{x} ({', '.join(g.index[g == x][:3])}{', ...' if (g == x).sum() > 3 else ''})" for x in set(g)}
        summ.to_csv(os.path.join(d, "group_unique_summary.tsv"), sep="\t", index=False)
        memb.to_csv(os.path.join(d, "group_membership_by_feature.tsv"), sep="\t")
        regions.rename("hits").rename_axis("groups").to_csv(os.path.join(d, "venn_regions.tsv"), sep="\t")
        plot_correlation(R, order, g, os.path.join(d, "concordance_pearson_r.pdf"), a.value_label)
        regions.index = regions.index.map(lambda s: "&".join(str(x) for x in s.split("&")))
        plot_venn(summ.assign(group=summ.group.astype(str)), regions, os.path.join(d, "venn_unique_hits.pdf"), labels)
        for x in sorted(set(g)):
            for col, tag in ((f"unique_{x}", "unique"), (f"unique_tested_elsewhere_{x}", "unique_tested_elsewhere")):
                genes = sorted({symbol(f) for f in memb.index[memb[col]]} - {"", "nan"})
                path = os.path.join(d, f"genes_{tag}_G{x}.txt")
                with open(path, "w") as fh:
                    fh.write("\n".join(genes) + "\n")
                if a.goelite_species and genes:
                    cmd = [sys.executable, "-m", "altanalyze3.components.goelite.main", "--species", a.goelite_species,
                           "--query", path, "--outdir", os.path.join(d, f"goelite_{tag}_G{x}"), "--name", f"{tag}_G{x}"]
                    if a.goelite_background:
                        cmd += ["--background", a.goelite_background]
                    subprocess.run(cmd, check=True)
        report["groups"][str(k)] = dict(populations={str(x): list(g.index[g == x]) for x in sorted(set(g))},
                                        summary=summ.to_dict("records"), regions=regions.astype(int).to_dict())
    with open(os.path.join(a.outdir, "concordance_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(json.dumps({k: v for k, v in report.items() if k in ("populations_used", "silhouette")}, indent=1, default=str))
    return report


if __name__ == "__main__":
    main()
