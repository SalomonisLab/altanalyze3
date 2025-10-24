"""
Core data structures for the GO-Elite component.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Iterable, Mapping, Optional, Sequence, Set, Tuple

import numpy as np


@dataclass(frozen=True)
class GOTermNode:
    term_id: str
    name: str
    namespace: str
    parents: Tuple[str, ...] = field(default_factory=tuple)
    children: Tuple[str, ...] = field(default_factory=tuple)
    depth: int = 0
    genes: frozenset[str] = field(default_factory=frozenset)

    def with_children(self, children: Iterable[str]) -> "GOTermNode":
        return GOTermNode(
            term_id=self.term_id,
            name=self.name,
            namespace=self.namespace,
            parents=self.parents,
            children=tuple(children),
            depth=self.depth,
            genes=self.genes,
        )

    def with_parents(self, parents: Iterable[str]) -> "GOTermNode":
        return GOTermNode(
            term_id=self.term_id,
            name=self.name,
            namespace=self.namespace,
            parents=tuple(parents),
            children=self.children,
            depth=self.depth,
            genes=self.genes,
        )

    def with_depth(self, depth: int) -> "GOTermNode":
        return GOTermNode(
            term_id=self.term_id,
            name=self.name,
            namespace=self.namespace,
            parents=self.parents,
            children=self.children,
            depth=depth,
            genes=self.genes,
        )

    def with_genes(self, genes: Iterable[str]) -> "GOTermNode":
        return GOTermNode(
            term_id=self.term_id,
            name=self.name,
            namespace=self.namespace,
            parents=self.parents,
            children=self.children,
            depth=self.depth,
            genes=frozenset(genes),
        )


@dataclass
class GOTree:
    nodes: Dict[str, GOTermNode]
    roots: Tuple[str, ...] = field(default_factory=tuple)
    namespace_index: Mapping[str, Tuple[str, ...]] = field(default_factory=dict)

    def get(self, term_id: str) -> Optional[GOTermNode]:
        return self.nodes.get(term_id)

    def parents(self, term_id: str) -> Tuple[str, ...]:
        node = self.get(term_id)
        return node.parents if node else tuple()

    def children(self, term_id: str) -> Tuple[str, ...]:
        node = self.get(term_id)
        return node.children if node else tuple()

    def descendants(self, term_id: str) -> Set[str]:
        desc: Set[str] = set()
        stack = list(self.children(term_id))
        while stack:
            current = stack.pop()
            if current in desc:
                continue
            desc.add(current)
            stack.extend(self.children(current))
        return desc

    def ancestors(self, term_id: str) -> Set[str]:
        anc: Set[str] = set()
        stack = list(self.parents(term_id))
        while stack:
            current = stack.pop()
            if current in anc:
                continue
            anc.add(current)
            stack.extend(self.parents(current))
        return anc

    def ensure_depths(self) -> None:
        queue = deque((root, 0) for root in self.roots)
        visited: Set[str] = set()
        while queue:
            term_id, depth = queue.popleft()
            if term_id in visited:
                continue
            visited.add(term_id)
            node = self.nodes.get(term_id)
            if not node:
                continue
            self.nodes[term_id] = node.with_depth(depth)
            for child in node.children:
                queue.append((child, depth + 1))

    def subgraph(self, terms: Iterable[str]) -> "GOTree":
        term_set = set(terms)
        return GOTree(
            nodes={tid: node for tid, node in self.nodes.items() if tid in term_set},
            roots=tuple(t for t in self.roots if t in term_set),
            namespace_index={
                ns: tuple(t for t in tids if t in term_set)
                for ns, tids in self.namespace_index.items()
            },
        )


@dataclass
class EnrichmentResult:
    term_id: str
    z_score: float
    p_value: float
    fdr: float
    overlap: int
    total_genes: int
    background: int
    selected: bool = False
    blocked_by: Optional[str] = None

    def significance_tuple(self) -> Tuple[float, float]:
        return (-abs(self.z_score), self.p_value)


def compute_z_score(hit_count: int, query_total: int, term_total: int, background_total: int) -> float:
    if query_total == 0 or background_total == 0:
        return 0.0
    term_fraction = term_total / background_total
    expected = term_fraction * query_total
    variance = query_total * term_fraction * (1 - term_fraction)
    if variance <= 0:
        return 0.0
    return (hit_count - expected) / np.sqrt(variance)
