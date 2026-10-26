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


REGULATES_RELATIONS = ("regulates", "positively_regulates", "negatively_regulates")


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
                elif relation.split()[0] in REGULATES_RELATIONS:
                    current.setdefault("regulates", []).append(relation.split("!")[0].split()[1].strip())
    if "id" in current:
        _commit_term(current, nodes, namespace_index, roots)

    # One pass over the nodes (file order kept), not one scan of every node per term: the former
    # quadratic loop cost 48,340^2 membership tests, minutes per build (2026-10-24 fix).
    children_of: Dict[str, List[str]] = {}
    for child_id, child in nodes.items():
        for parent in child.parents:
            children_of.setdefault(parent, []).append(child_id)
    for term_id, node in list(nodes.items()):
        nodes[term_id] = node.with_children(children_of.get(term_id, []))

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
        regulates=tuple(current.get("regulates", [])),
    )
    nodes[term_id] = node
    namespace_index.setdefault(node.namespace, []).append(term_id)
    if not node.parents:
        roots.add(term_id)


def parse_gaf(path: str, evidence: Optional[Set[str]] = None, exclude_not: bool = True) -> pd.DataFrame:
    """Gene-to-term rows of a GAF file. Rows whose qualifier carries NOT are negative statements
    ("this gene does NOT have this function") and are excluded by default (2026-10-24 fix)."""
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
            if exclude_not and "NOT" in parts[3].upper().split("|"):
                continue
            gene_id = parts[1] or parts[2]
            gene_symbol = parts[2]
            records.append((gene_id, gene_symbol, go_id))

    df = pd.DataFrame(records, columns=["gene_id", "gene_symbol", "go_id"])
    return df.drop_duplicates()


def parse_gaf_synonyms(path: str) -> Dict[str, str]:
    """Map every synonym (GAF column 11, pipe-separated) and alias to the current db_object_symbol,
    upper-cased, so a query written with an older gene symbol still reaches its annotations.
    A synonym that is itself a current symbol is never remapped."""
    current: Set[str] = set()
    pairs: Dict[str, Set[str]] = {}
    with _open_text(path) as handle:
        for raw in handle:
            if raw.startswith("!"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            symbol = parts[2].strip()
            if not symbol:
                continue
            current.add(symbol.upper())
            for syn in parts[10].split("|"):
                syn = syn.strip().upper()
                if syn and syn != symbol.upper():
                    pairs.setdefault(syn, set()).add(symbol)
    return {syn: sorted(symbols)[0] for syn, symbols in pairs.items()
            if syn not in current and len(symbols) == 1}


def propagate_annotations(tree: GOTree, term_to_genes: Dict[str, Set[str]], use_regulates: bool = True) -> Dict[str, Set[str]]:
    """GO true-path rule: a gene annotated to a term is annotated to every ancestor (is_a, part_of) and,
    with ``use_regulates`` (the default since 2026-10-24), to the terms its term regulates, positively
    regulates or negatively regulates, transitively. ToppFun propagates over regulates ("leukocyte
    activation" holds 1,224 genes there and 512 with is_a/part_of alone); Nathan's benchmark is ToppFun.

    Without propagation a broad term such as "regulation of leukocyte activation" holds only its direct
    annotations and never ranks, while leaf terms are over-weighted (2026-10-24 fix). Idempotent.
    """
    propagated: Dict[str, Set[str]] = {tid: set(genes) for tid, genes in term_to_genes.items()}
    ancestor_cache: Dict[str, Set[str]] = {}

    def _ancestors(term_id: str) -> Set[str]:
        out: Set[str] = set()
        node = tree.nodes[term_id]
        stack = list(node.parents) + (list(node.regulates) if use_regulates else [])
        while stack:
            current = stack.pop()
            if current in out:
                continue
            out.add(current)
            cur = tree.nodes.get(current)
            if cur is None:
                continue
            stack.extend(cur.parents)
            if use_regulates:
                stack.extend(cur.regulates)
        return out

    for term_id, genes in term_to_genes.items():
        if term_id not in tree.nodes:
            continue
        if term_id not in ancestor_cache:
            ancestor_cache[term_id] = _ancestors(term_id)
        for anc in ancestor_cache[term_id]:
            propagated.setdefault(anc, set()).update(genes)
    return propagated


def parse_gene_orthologs(path: str, taxon_id: int, other_taxon_id: int, gene_id_to_symbol: Dict[str, str],
                         other_gene_id_to_symbol: Dict[str, str]) -> Dict[str, Set[str]]:
    """NCBI gene_orthologs (tax_id, GeneID, relationship, Other_tax_id, Other_GeneID), either row direction,
    relationship Ortholog only. Returns other-species symbol -> set of this-species symbols. GeneIDs absent
    from a gene_info table are counted in ``attrs``-like keys of the returned dict under ``__unmapped__``."""
    this_tax, other_tax = str(taxon_id), str(other_taxon_id)
    out: Dict[str, Set[str]] = {}
    unmapped = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[2] != "Ortholog":
                continue
            if f[0] == this_tax and f[3] == other_tax:
                this_id, other_id = f[1], f[4]
            elif f[0] == other_tax and f[3] == this_tax:
                this_id, other_id = f[4], f[1]
            else:
                continue
            this_sym, other_sym = gene_id_to_symbol.get(this_id), other_gene_id_to_symbol.get(other_id)
            if this_sym is None or other_sym is None:
                unmapped += 1
                continue
            out.setdefault(other_sym, set()).add(this_sym)
    out["__unmapped__"] = {str(unmapped)}
    return out


def transfer_ortholog_annotations(frame: pd.DataFrame, orthologs: Dict[str, Set[str]]) -> pd.DataFrame:
    """Annotations of the other species (gene_symbol, go_id) carried to this species' ortholog symbols.
    One other-species gene with several orthologs annotates all of them. Returns gene_symbol/go_id rows."""
    rows = []
    for other_sym, go_id in frame[["gene_symbol", "go_id"]].drop_duplicates().itertuples(index=False):
        for this_sym in orthologs.get(other_sym, ()):
            rows.append((this_sym, go_id))
    df = pd.DataFrame(rows, columns=["gene_symbol", "go_id"]).drop_duplicates().reset_index(drop=True)
    df.attrs["source_rows"] = int(len(frame))
    df.attrs["transferred_genes"] = int(df["gene_symbol"].nunique()) if len(df) else 0
    return df


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
    term_to_genes = propagate_annotations(tree, term_to_genes)

    for term_id, node in list(tree.nodes.items()):
        if term_id in term_to_genes:
            tree.nodes[term_id] = node.with_genes(term_to_genes[term_id])

    return ParsedGO(tree=tree, term_to_genes=term_to_genes)


def build_flat_gene_sets(path: str, label: str = "WikiPathways") -> ParsedGO:
    """A flat gene-set collection (no hierarchy) as a ParsedGO, for pathway enrichment.

    Accepted layouts, tab-delimited with a header:
      * AltAnalyze GO-Elite ``gene-mapp`` files (``UID``, ``SystemCode``, ``MAPP``), where MAPP is
        ``<pathway name>:<WP id>`` (EnsMart100 ``Ensembl-MAPP.txt`` is the WikiPathways collection);
      * two columns ``gene`` and ``set`` (optional third column ``name``).
    Term ids are the WP id when present, else the set string; every node is a root, so the
    GO-Elite hierarchy pruning has nothing to prune and each pathway is reported on its own.
    """
    path = str(path)
    term_to_genes: Dict[str, Set[str]] = {}
    names: Dict[str, str] = {}
    with _open_text(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        cols = {c.strip().lower(): i for i, c in enumerate(header)}
        if "uid" in cols and "mapp" in cols:
            gi, si, ni = cols["uid"], cols["mapp"], None
        elif "gene" in cols and "set" in cols:
            gi, si, ni = cols["gene"], cols["set"], cols.get("name")
        else:
            raise ValueError("Gene-set file needs columns UID/SystemCode/MAPP or gene/set[/name]: " + path)
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= max(gi, si) or not parts[gi].strip() or not parts[si].strip():
                continue
            gene, raw = parts[gi].strip(), parts[si].strip()
            if ni is not None:
                term_id, name = raw, (parts[ni].strip() if len(parts) > ni and parts[ni].strip() else raw)
            elif ":WP" in raw:
                name, _, wp = raw.rpartition(":")
                term_id = wp
            else:
                term_id, name = raw, raw
            term_to_genes.setdefault(term_id, set()).add(gene)
            names.setdefault(term_id, name)
    if not term_to_genes:
        raise ValueError("Gene-set file holds no gene/set rows: " + path)
    nodes = {tid: GOTermNode(term_id=tid, name=names[tid], namespace=label, genes=frozenset(genes))
             for tid, genes in term_to_genes.items()}
    tree = GOTree(nodes=nodes, roots=tuple(sorted(nodes)), namespace_index={label: tuple(sorted(nodes))})
    return ParsedGO(tree=tree, term_to_genes=term_to_genes)


def parse_gene_info(path: str, taxon_id: Optional[int] = None) -> Tuple[Dict[str, str], Dict[str, str]]:
    """NCBI gene_info: GeneID -> Symbol, and synonym (upper) -> Symbol (Synonyms column, pipe-separated)."""
    gene_id_to_symbol: Dict[str, str] = {}
    synonyms: Dict[str, Set[str]] = {}
    with _open_text(str(path)) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            if taxon_id is not None and parts[0] != str(taxon_id):
                continue
            gene_id, symbol = parts[1], parts[2].strip()
            if not symbol:
                continue
            gene_id_to_symbol[gene_id] = symbol
            for syn in parts[4].split("|"):
                syn = syn.strip().upper()
                if syn and syn != "-" and syn != symbol.upper():
                    synonyms.setdefault(syn, set()).add(symbol)
    current = {v.upper() for v in gene_id_to_symbol.values()}
    return gene_id_to_symbol, {k: sorted(v)[0] for k, v in synonyms.items() if k not in current and len(v) == 1}


def parse_gene2go(path: str, gene_id_to_symbol: Dict[str, str], taxon_id: int, exclude_not: bool = True) -> pd.DataFrame:
    """NCBI gene2go rows of one taxon as gene_symbol / go_id (GeneIDs mapped through gene_info).
    Rows whose Qualifier carries NOT are excluded by default. Unmapped GeneIDs are counted in
    ``attrs['unmapped_gene_ids']`` of the returned frame."""
    records: List[Tuple[str, str]] = []
    unmapped: Set[str] = set()
    taxon = str(taxon_id)
    prefix = taxon + "\t"
    with _open_text(str(path)) as handle:
        for raw in handle:
            # gene2go holds every taxon (tens of millions of lines); test the taxon prefix before splitting
            if not raw.startswith(prefix):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            if exclude_not and "NOT" in parts[4].upper().replace(" ", "|").split("|"):
                continue
            symbol = gene_id_to_symbol.get(parts[1])
            if symbol is None:
                unmapped.add(parts[1])
                continue
            records.append((symbol, parts[2]))
    df = pd.DataFrame(records, columns=["gene_symbol", "go_id"]).drop_duplicates()
    df.attrs["unmapped_gene_ids"] = len(unmapped)
    return df


def build_tree_with_annotation_frames(obo_path: str, frames: Dict[str, pd.DataFrame],
                                      use_regulates: bool = True) -> Tuple[ParsedGO, Dict[str, int]]:
    """Union of several gene_symbol/go_id annotation frames (GOA GAF, NCBI gene2go, ortholog-transferred
    annotations), then the true-path propagation (over regulates too unless ``use_regulates`` is False).
    Returns the ParsedGO and the per-source and union row counts."""
    tree = parse_obo(obo_path)
    counts: Dict[str, int] = {}
    union = []
    for name, df in frames.items():
        d = df[["gene_symbol", "go_id"]].drop_duplicates()
        counts[name] = int(len(d))
        union.append(d)
    if not union:
        raise ValueError("No annotation source available to build GO-Elite resources")
    all_rows = pd.concat(union, ignore_index=True).drop_duplicates()
    counts["union"] = int(len(all_rows))
    term_to_genes: Dict[str, Set[str]] = {}
    for term_id, group in all_rows.groupby("go_id"):
        if term_id in tree.nodes:
            term_to_genes[term_id] = set(group["gene_symbol"].astype(str))
    term_to_genes = propagate_annotations(tree, term_to_genes, use_regulates=use_regulates)
    for term_id, node in list(tree.nodes.items()):
        if term_id in term_to_genes:
            tree.nodes[term_id] = node.with_genes(term_to_genes[term_id])
    return ParsedGO(tree=tree, term_to_genes=term_to_genes), counts
