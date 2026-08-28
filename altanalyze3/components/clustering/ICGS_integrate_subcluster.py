"""Sub-cluster replacement for ICGS_integrate: let a finer dataset replace a coarse state.

ICGS_integrate has 2 outcomes for a candidate cluster: admit it as a new state, or judge it
redundant against an existing state and drop its cells (ICGS_integrate.py:1487-1489 keeps only
cells whose cluster is a state member). It has no third outcome, so the seed dataset's
granularity becomes a ceiling: on 7 human lung datasets, CellRef's single `AM` state absorbed all
3 of HLCA's alveolar macrophage states.

This module supplies the missing outcome. When dataset D publishes k separate cell states and all
k are judged redundant against ONE reference state S, the audit already records that S is
under-resolved: D distinguishes k populations there and the reference distinguishes 1. The
evidence is combinatorial, so no cell of D is ever compared with a cell of S. That matters,
because ICGS_INTEGRATE_METHODS.md section 1.2 rules out cross-study cell comparison, and 2
earlier prototypes that ignored that rule failed: their negative controls scored more extreme
than their true positives.

Usage:
    python -m altanalyze3.components.clustering.ICGS_integrate_subcluster \
        <integration_dir> <runs_root> [--apply]

`runs_root` holds one directory per dataset, each with ICGS3_run/MarkerFinder/.
`--apply` writes harmonized_states.subcluster_replaced.tsv beside the input.
"""
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

RHO = 0.3
MIN_MARKERS = 1
FDR = 0.05
TOP_N = 60
MIN_OVERLAP = 10


def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def markers_for(root, ds):
    f = os.path.join(root, ds, "ICGS3_run", "MarkerFinder", "icgs3_markers_all_correlations.tsv")
    d = pd.read_csv(f, sep="\t")
    d["top_cluster"] = d["top_cluster"].astype(str)
    d["marker"] = d["marker"].astype(str)
    keep = d[d["pearson_r"] >= RHO]
    top = {c: list(g.sort_values("pearson_r", ascending=False)["marker"].head(TOP_N))
           for c, g in keep.groupby("top_cluster")}
    return top, set(d["marker"])


def detect(intdir, root):
    a = pd.read_csv(os.path.join(intdir, "hierarchical_audit.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(intdir, "harmonized_states.tsv"), sep="\t")
    name = {r["state"]: (r["name"] if isinstance(r.get("name"), str) else r["state"])
            for _, r in s.iterrows()}
    members = {r["state"]: str(r["members"]) for _, r in s.iterrows()}
    exc = a[a["action"] == "excluded_redundant"].copy()
    n_all = len(exc)
    exc = exc[(exc["fdr_r025"] <= FDR) & (exc["overlap_r025"] >= MIN_OVERLAP)]
    print("exclusions carrying significant gate-1 evidence: %d of %d (dropped %d with "
          "FDR > %.2f or overlap < %d)" % (len(exc), n_all, n_all - len(exc), FDR, MIN_OVERLAP))
    cache, rows = {}, []
    for (ds, st), grp in exc.groupby(["step", "nearest_state_r025"]):
        clusters = sorted(grp["cluster"].astype(str))
        if len(clusters) < 2:
            continue
        if ds not in cache:
            cache[ds] = markers_for(root, ds)
        top, universe = cache[ds]
        N = len(universe)
        ok = [c for c in clusters if len(top.get(c, [])) >= MIN_MARKERS]
        pairs, pv = list(combinations(ok, 2)), []
        for x, y in pairs:
            A, B = set(top[x]), set(top[y])
            k = len(A & B)
            pv.append(hypergeom.sf(k - 1, N, len(A), len(B)) if k else 1.0)
        q = bh(np.array(pv)) if pv else np.array([])
        distinct = {c: True for c in ok}
        for (x, y), qq in zip(pairs, q):
            if qq <= FDR:
                distinct[x] = distinct[y] = False
        survivors = [c for c in ok if distinct[c]]
        rows.append({"reference_state": st, "reference_name": name.get(st, st),
                     "reference_members": members.get(st, ""),
                     "finer_dataset": ds, "n_clusters_collapsed": len(clusters),
                     "clusters": ",".join(clusters),
                     "with_enough_markers": len(ok),
                     "mutually_distinct": len(survivors),
                     "survivors": ",".join(survivors),
                     "verdict": "REPLACE" if len(survivors) >= 2 else "keep"})
    return pd.DataFrame(rows).sort_values(["verdict", "n_clusters_collapsed"],
                                          ascending=[True, False])


def replacement_flags(table):
    """The two --force-admit / --suppress-cluster values a second ICGS_integrate pass needs.

    Writing a revised state table would not be a result. A state carries markers, a centroid, an
    annotation, a row in the cell table and a column in the heatmap, and none of those can be
    produced by editing a TSV. So this function emits the flags instead, and the replacement runs
    through the validated entry point, which writes the full standard output set.
    """
    rep = table[table["verdict"] == "REPLACE"]
    admit, suppress = [], []
    for _, r in rep.iterrows():
        for c in str(r["survivors"]).split(","):
            admit.append("%s|%s" % (r["finer_dataset"], c))
        for m in str(r["reference_members"]).split(","):
            if m:
                suppress.append(m)
    admit = sorted(set(admit))
    suppress = sorted(set(suppress))
    return ",".join(admit), ",".join(suppress)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    intdir, root = args[0], args[1]
    t = detect(intdir, root)
    if t.empty:
        print("no reference state is collapsed by 2 or more clusters of one dataset")
        sys.exit(0)
    rep = t[t["verdict"] == "REPLACE"]
    print("candidate splits: %d reference states collapse >=2 clusters of one dataset; "
          "%d pass the test" % (len(t), len(rep)))
    print("net states if enforced: +%d" % int((rep["mutually_distinct"] - 1).sum()))
    print()
    for _, r in t.iterrows():
        print("  %-8s %-34s <- %-8s %d clusters, %d distinct  %s" % (
            r["reference_state"], str(r["reference_name"])[:34], r["finer_dataset"],
            r["n_clusters_collapsed"], r["mutually_distinct"], r["verdict"]))
        print("        %s" % r["clusters"][:150])
    out = os.path.join(intdir, "subcluster_replacement_candidates.tsv")
    t.to_csv(out, sep="\t", index=False)
    print("\nwrote %s" % out)
    admit, suppress = replacement_flags(t)
    fl = os.path.join(intdir, "subcluster_replacement_flags.txt")
    with open(fl, "w") as h:
        h.write("--force-admit\n%s\n\n--suppress-cluster\n%s\n" % (admit, suppress))
    print("\nwrote %s" % fl)
    print("\nSecond pass: re-run ICGS_integrate with")
    print("  --force-admit '%s'" % admit)
    print("  --suppress-cluster '%s'" % suppress)
