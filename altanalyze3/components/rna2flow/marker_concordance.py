"""The strict success measure: does a transfer reproduce each population's antibody markers?

ARI against FlowSOM asks whether the partitions agree. That is weak: FlowSOM is itself one
clustering of 20 channels, and a transfer can be biologically right while partitioning
differently. The stronger question is whether the SAME populations are defined by the SAME
antibodies, in a similar order, on both platforms.

Procedure, per label set:
  1. MarkerFinder on the CITE-seq panel with the original labels   -> reference ranking
  2. MarkerFinder on the flow panel with the TRANSFERRED labels    -> query ranking
  3. For each population present in both, compare the shared-marker rankings.

A population PASSES when its top flow marker appears in the CITE top-k for that population and
the two rankings correlate positively. The headline is how many populations pass; the target is
at least 7.

MarkerFinder is altanalyze3's own point-biserial statistic
(components/cellHarmony/markerFinder.marker_finder), not a reimplementation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from ..cellHarmony.markerFinder import marker_finder

__all__ = ["rank_markers", "concordance", "MIN_POPULATIONS"]

MIN_POPULATIONS = 7


def rank_markers(X, labels, feature_names, min_cells: int = 20):
    """Per-population ranking of features by MarkerFinder rho. Returns {population: Series}."""
    labels = np.asarray(labels).astype(str)
    keep_pop = {p for p, n in zip(*np.unique(labels, return_counts=True)) if n >= min_cells}
    mask = np.isin(labels, list(keep_pop))
    if mask.sum() == 0 or len(keep_pop) < 2:
        return {}
    r_df, _ = marker_finder(np.asarray(X, dtype=np.float64)[mask], list(labels[mask]),
                            gene_names=list(feature_names), validate_scaling=False)
    return {p: r_df[p].sort_values(ascending=False) for p in r_df.columns}


def concordance(cite_X, cite_labels, flow_X, flow_labels, pairs, top_k: int = 3,
                min_cells: int = 20):
    """Compare per-population marker rankings across platforms.

    `pairs` maps flow marker -> cite feature, so the two rankings share one index.
    """
    cite_rank = rank_markers(cite_X, cite_labels, list(pairs["cite"]), min_cells)
    flow_rank = rank_markers(flow_X, flow_labels, list(pairs["flow"]), min_cells)
    f2c = dict(zip(pairs["flow"], pairs["cite"]))

    rows = []
    for pop in sorted(set(cite_rank) & set(flow_rank)):
        c = cite_rank[pop]
        f = flow_rank[pop].rename(index=f2c)          # put flow on CITE feature names
        shared = [i for i in f.index if i in c.index]
        if len(shared) < 3:
            continue
        rho = spearmanr(c.loc[shared].values, f.loc[shared].values).statistic
        top_flow = f.loc[shared].idxmax()
        cite_top_k = list(c.loc[shared].sort_values(ascending=False).index[:top_k])
        rows.append({
            "population": pop,
            "n_shared_markers": len(shared),
            "top_flow_marker": top_flow,
            "cite_top_%d" % top_k: ", ".join(cite_top_k),
            "top_marker_in_cite_top_k": top_flow in cite_top_k,
            "rank_spearman": None if rho is None or np.isnan(rho) else round(float(rho), 4),
            "passes": bool(top_flow in cite_top_k and rho is not None and not np.isnan(rho) and rho > 0),
        })
    df = pd.DataFrame(rows)
    summary = {
        "populations_compared": int(len(df)),
        "populations_passing": int(df["passes"].sum()) if len(df) else 0,
        "median_rank_spearman": float(df["rank_spearman"].median()) if len(df) else float("nan"),
        "meets_target": bool(len(df) and df["passes"].sum() >= MIN_POPULATIONS),
    }
    return df, summary
