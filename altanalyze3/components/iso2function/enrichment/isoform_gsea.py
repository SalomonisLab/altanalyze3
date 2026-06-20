"""Layer 4 (enrichment) — isoform-focused gene-set enrichment.

Tests whether the interaction partners / DNA-target genes that are GAINED or LOST across an isoform
switch (from ``network.differential_interactome``) are enriched for a pathway / GO term / interactor
module — i.e. whether the switch concentrates its interaction change on a coherent functional unit.

Key honesty choice: the background is the **assayed space** (the co-tested target universe), not the
whole genome. Enrichment of a gained/lost partner set is only meaningful relative to what *could* have
been gained/lost given what was tested. Callers pass that background explicitly.

Statistics use scipy only (no R): hypergeometric tail probability + Benjamini-Hochberg FDR.
"""

import logging

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)


def benjamini_hochberg(pvalues):
    """Return BH-FDR adjusted q-values for a 1-D array-like of p-values (NaNs preserved)."""
    p = np.asarray(pvalues, dtype=float)
    mask = ~np.isnan(p)
    q = np.full_like(p, np.nan, dtype=float)
    pv = p[mask]
    if pv.size == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    m = pv.size
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]   # enforce monotonicity
    out = np.empty(m, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    q[mask] = out
    return q


def hypergeometric_test(query, gene_set, background):
    """One gene set vs a query, restricted to a background universe.

    Returns dict(k, n, K, N, fold_enrichment, p_value):
      N = |background|, K = |gene_set ∩ background|, n = |query ∩ background|,
      k = |query ∩ gene_set ∩ background|. p = P(X >= k) under the hypergeometric null.
    """
    bg = set(map(str, background))
    q = set(map(str, query)) & bg
    gs = set(map(str, gene_set)) & bg
    N, K, n = len(bg), len(gs), len(q)
    k = len(q & gs)
    if N == 0 or K == 0 or n == 0:
        return {"k": k, "n": n, "K": K, "N": N, "fold_enrichment": np.nan, "p_value": np.nan}
    p = float(hypergeom.sf(k - 1, N, K, n))  # P(X >= k)
    expected = n * K / N
    fold = (k / expected) if expected > 0 else np.nan
    return {"k": k, "n": n, "K": K, "N": N, "fold_enrichment": fold, "p_value": p}


def enrich(query_genes, gene_sets, background, min_overlap=1):
    """Test a query gene list against many gene sets within a background universe.

    query_genes : iterable of partner/target gene symbols (e.g. the 'gained' set of a switch).
    gene_sets   : dict {set_name: iterable_of_genes} (pathways/GO/interactor modules).
    background  : iterable defining the universe (e.g. co-tested targets) — REQUIRED, no genome default.
    Returns a DataFrame sorted by p_value with a BH q_value column; rows with overlap < min_overlap
    are dropped before FDR so the correction reflects only testable sets.
    """
    rows = []
    for name, gs in gene_sets.items():
        res = hypergeometric_test(query_genes, gs, background)
        res["gene_set"] = name
        rows.append(res)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df[df["k"] >= min_overlap].copy()
    if df.empty:
        return df
    df = df.sort_values("p_value", kind="mergesort").reset_index(drop=True)
    df["q_value"] = benjamini_hochberg(df["p_value"].values)
    cols = ["gene_set", "k", "n", "K", "N", "fold_enrichment", "p_value", "q_value"]
    return df[cols]


def gene_sets_from_goelite():
    """Hook for pathway/GO gene sets via the bundled GO-Elite engine (``components.goelite``).

    Deferred wiring: GO-Elite owns the curated pathway/ontology resources used elsewhere in
    altanalyze3. The plan is to load its gene-set membership and feed it to :func:`enrich`, rather than
    duplicating any annotation here. Until wired, callers supply ``gene_sets`` explicitly (e.g.
    interactor modules from the atlas itself).
    """
    raise NotImplementedError(
        "GO-Elite gene-set loading not yet wired; pass `gene_sets` explicitly (see GOALS.md P4). "
        "Interactor-module gene sets can be built directly from the atlas partner lists in the meantime."
    )
