"""
GO DAG-based term prioritisation for GO-Elite.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from .structures import EnrichmentResult, GOTree


@dataclass
class PrioritizationSettings:
    min_z: float = 1.96
    max_fdr: float = 0.1
    min_overlap: int = 2
    delta_z: float = 0.5


def prioritize_terms(
    tree: GOTree,
    results: Sequence[EnrichmentResult],
    settings: Optional[PrioritizationSettings] = None,
) -> List[EnrichmentResult]:
    cfg = settings or PrioritizationSettings()
    selected: Dict[str, EnrichmentResult] = {}
    result_map = {res.term_id: res for res in results}
    ordered = sorted(results, key=lambda res: res.significance_tuple())

    for res in ordered:
        if abs(res.z_score) < cfg.min_z or res.fdr > cfg.max_fdr or res.overlap < cfg.min_overlap:
            continue

        ancestors = tree.ancestors(res.term_id)
        blocker = _blocking_ancestor(res, ancestors, selected, cfg.delta_z)
        if blocker:
            res.blocked_by = blocker
            continue

        res.selected = True
        selected[res.term_id] = res

        for child_id in tree.descendants(res.term_id):
            child = result_map.get(child_id)
            if not child or child.selected:
                continue
            if abs(child.z_score) + cfg.delta_z <= abs(res.z_score):
                child.blocked_by = res.term_id

    return list(result_map.values())


def _blocking_ancestor(
    res: EnrichmentResult,
    ancestors: Iterable[str],
    selected: Dict[str, EnrichmentResult],
    delta_z: float,
) -> Optional[str]:
    for ancestor_id in ancestors:
        ancestor = selected.get(ancestor_id)
        if not ancestor:
            continue
        if abs(ancestor.z_score) >= abs(res.z_score) - delta_z:
            return ancestor.term_id
    return None
