"""
Ontology and annotation parsing utilities for GO-Elite.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from .structures import GOTermNode, GOTree


@dataclass
class ParsedGO:
    tree: GOTree
    term_to_genes: Dict[str, Set[str]]


def _open_text(path: str) -> io.TextIOWrapper:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")


def parse_obo(path: str) -> GOTree:
    nodes: Dict[str, GOTermNode] = {}
    current: Dict[str, List[str]] = {}
    roots: Set[str] = set()
    namespace_index: Dict[str, List[str]] = {}

    with _open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line == "[Term]":
                if "id" in current:
                    _commit_term(current, nodes, namespace_index, roots)
                current = {}
                continue
            if line.startswith("id:"):
                current.setdefault("id", []).append(line.split("id:")[1].strip())
            elif line.startswith("name:"):
                current.setdefault("name", []).append(line.split("name:")[1].strip())
            elif line.startswith("namespace:"):
                current.setdefault("namespace", []).append(line.split("namespace:")[1].strip())
            elif line.startswith("is_a:"):
                parent = line.split("is_a:")[1].split("!")[0].strip()
                current.setdefault("is_a", []).append(parent)
            elif line.startswith("relationship:"):
                relation = line.split("relationship:")[1].strip()
                if relation.startswith("part_of"):
                    parent = relation.split("part_of")[1].split("!")[0].strip()
                    current.setdefault("part_of", []).append(parent)
    if "id" in current:
        _commit_term(current, nodes, namespace_index, roots)

    for term_id, node in list(nodes.items()):
        children = [
            child_id for child_id, child in nodes.items()
            if term_id in child.parents
        ]
        nodes[term_id] = node.with_children(children)

    tree = GOTree(
        nodes=nodes,
        roots=tuple(roots),
        namespace_index={ns: tuple(ids) for ns, ids in namespace_index.items()},
    )
    tree.ensure_depths()
    return tree


def _commit_term(
    current: Dict[str, List[str]],
    nodes: Dict[str, GOTermNode],
    namespace_index: Dict[str, List[str]],
    roots: Set[str],
) -> None:
    term_id = current["id"][0]
    node = GOTermNode(
        term_id=term_id,
        name=current.get("name", [""])[0],
        namespace=current.get("namespace", ["unknown"])[0],
        parents=tuple(current.get("is_a", []) + current.get("part_of", [])),
    )
    nodes[term_id] = node
    namespace_index.setdefault(node.namespace, []).append(term_id)
    if not node.parents:
        roots.add(term_id)


def parse_gaf(path: str, evidence: Optional[Set[str]] = None) -> pd.DataFrame:
    columns = [
        "db", "db_object_id", "db_object_symbol", "qualifier", "go_id",
        "db_reference", "evidence_code", "with_from", "aspect",
        "db_object_name", "db_object_synonym", "db_object_type",
        "taxon", "date", "assigned_by", "annotation_extension", "gene_product_form_id",
    ]
    records: List[Tuple[str, str, str]] = []

    with _open_text(path) as handle:
        for raw in handle:
            if raw.startswith("!"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < len(columns):
                continue
            go_id = parts[4]
            if not go_id.startswith("GO:"):
                continue
            evidence_code = parts[6]
            if evidence and evidence_code not in evidence:
                continue
            gene_id = parts[1] or parts[2]
            gene_symbol = parts[2]
            records.append((gene_id, gene_symbol, go_id))

    df = pd.DataFrame(records, columns=["gene_id", "gene_symbol", "go_id"])
    return df.drop_duplicates()


def build_tree_with_annotations(
    obo_path: str,
    gaf_path: str,
    acceptable_evidence: Optional[Set[str]] = None,
) -> ParsedGO:
    tree = parse_obo(obo_path)
    annotations = parse_gaf(gaf_path, evidence=acceptable_evidence)

    term_to_genes: Dict[str, Set[str]] = {}
    for term_id, group in annotations.groupby("go_id"):
        term_to_genes[term_id] = set(group["gene_symbol"].astype(str))

    for term_id, node in list(tree.nodes.items()):
        if term_id in term_to_genes:
            tree.nodes[term_id] = node.with_genes(term_to_genes[term_id])

    return ParsedGO(tree=tree, term_to_genes=term_to_genes)
