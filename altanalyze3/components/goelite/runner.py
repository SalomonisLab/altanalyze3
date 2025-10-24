"""
High-level orchestration for GO-Elite enrichment.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence

import logging
import numpy as np

try:
    from scipy.stats import hypergeom
except ImportError:  # pragma: no cover
    hypergeom = None

from .parser import ParsedGO
from .prio import PrioritizationSettings, prioritize_terms
from .structures import EnrichmentResult, GOTree, compute_z_score


def _bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    v = np.asarray(p_values, dtype=float)
    n = v.size
    if n == 0:
        return np.array([])
    order = np.argsort(v)
    ranked = np.empty_like(order, dtype=float)
    ranked[order] = np.arange(1, n + 1)
    fdr = v * n / ranked
    fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
    adjusted = np.empty_like(fdr_sorted)
    adjusted[order] = np.minimum(fdr_sorted, 1.0)
    return adjusted


@dataclass
class EnrichmentSettings:
    min_term_size: int = 5
    max_term_size: int = 2000
    prioritization: PrioritizationSettings = field(default_factory=PrioritizationSettings)


class GOEliteRunner:
    def __init__(
        self,
        go_data: ParsedGO,
        logger=None,
        settings: Optional[EnrichmentSettings] = None,
    ) -> None:
        self.go_tree: GOTree = go_data.tree
        self.term_to_genes = go_data.term_to_genes
        self.settings = settings or EnrichmentSettings()
        if logger is None:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger("go_elite")
        else:
            self.logger = logger

    def run(
        self,
        query_genes: Iterable[str],
        background_genes: Iterable[str],
        *,
        apply_prioritization: bool = True,
    ) -> List[EnrichmentResult]:
        query = {g.upper() for g in query_genes}
        background = {g.upper() for g in background_genes}

        if not background:
            raise ValueError("Background gene set is empty.")

        q_total = len(query)
        bg_total = len(background)
        results: List[EnrichmentResult] = []

        for term_id, gene_set in self.term_to_genes.items():
            genes = {g.upper() for g in gene_set if g.upper() in background}
            term_size = len(genes)
            if term_size < self.settings.min_term_size or term_size > self.settings.max_term_size:
                continue
            hits = len(query & genes)
            if hits == 0:
                continue

            z_score = compute_z_score(hits, q_total, term_size, bg_total)
            p_val = self._hypergeom_p(hits, q_total, term_size, bg_total)
            results.append(
                EnrichmentResult(
                    term_id=term_id,
                    z_score=z_score,
                    p_value=p_val,
                    fdr=1.0,
                    overlap=hits,
                    total_genes=term_size,
                    background=bg_total,
                )
            )

        if not results:
            return []

        fdr = _bh_fdr([res.p_value for res in results])
        for res, adj in zip(results, fdr):
            res.fdr = float(adj)

        if apply_prioritization:
            results = prioritize_terms(self.go_tree, results, settings=self.settings.prioritization)

        return results

    @staticmethod
    def _hypergeom_p(k: int, query_total: int, term_total: int, background_total: int) -> float:
        if hypergeom is None:
            from math import comb

            max_x = min(query_total, term_total)
            numerator = 0.0
            denominator = comb(background_total, query_total)
            for x in range(k, max_x + 1):
                numerator += comb(term_total, x) * comb(background_total - term_total, query_total - x)
            if denominator == 0:
                return 1.0
            return float(min(1.0, numerator / denominator))

        rv = hypergeom(M=background_total, n=term_total, N=query_total)
        return float(rv.sf(k - 1))
