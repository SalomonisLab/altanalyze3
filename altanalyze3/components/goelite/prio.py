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
    min_overlap: int = 3
    delta_z: float = 0.5
    parent_min_overlap_ratio: float = 1.25


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

        descendants = tree.descendants(res.term_id)
        blocker_desc = _blocking_descendant(
            res,
            descendants,
            selected,
            cfg.parent_min_overlap_ratio,
            cfg.delta_z,
        )
        if blocker_desc:
            res.blocked_by = blocker_desc
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
            if child.p_value >= res.p_value * (1 + cfg.delta_z):
                child.blocked_by = res.term_id

    _prune_selected_parents(tree, selected, cfg.parent_min_overlap_ratio, cfg.delta_z)
    return list(result_map.values())


def _blocking_descendant(
    res: EnrichmentResult,
    descendants: Iterable[str],
    selected: Dict[str, EnrichmentResult],
    min_overlap_ratio: float,
    delta_ratio: float,
) -> Optional[str]:
    res_overlap = res.overlap or 0
    for desc_id in descendants:
        desc = selected.get(desc_id)
        if not desc:
            continue
        if desc.p_value <= res.p_value * (1 + delta_ratio):
            min_required = (desc.overlap or 0) * min_overlap_ratio
            if res_overlap <= min_required:
                return desc.term_id
    return None


def _prune_selected_parents(
    tree: GOTree,
    selected: Dict[str, EnrichmentResult],
    min_overlap_ratio: float,
    delta_ratio: float,
) -> None:
    for term_id in list(selected.keys()):
        res = selected.get(term_id)
        if res is None or not res.selected:
            continue
        res_overlap = res.overlap or 0
        for desc_id in tree.descendants(term_id):
            desc = selected.get(desc_id)
            if not desc:
                continue
            if desc.p_value <= res.p_value * (1 + delta_ratio):
                min_required = (desc.overlap or 0) * min_overlap_ratio
                if res_overlap <= min_required:
                    res.selected = False
                    res.blocked_by = desc.term_id
                    selected.pop(term_id, None)
                    break


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
        if ancestor.p_value <= res.p_value * (1 + delta_z):
            return ancestor.term_id
    return None
