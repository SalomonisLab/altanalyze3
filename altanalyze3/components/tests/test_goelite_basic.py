"""
Basic GO-Elite tests to verify prioritisation and enrichment logic.
"""

from __future__ import annotations

from pathlib import Path

from altanalyze3.components.goelite.prio import PrioritizationSettings, prioritize_terms
from altanalyze3.components.goelite.runner import EnrichmentSettings, GOEliteRunner
from altanalyze3.components.goelite.structures import EnrichmentResult, GOTermNode, GOTree
from altanalyze3.components.goelite.parser import ParsedGO
from altanalyze3.components.goelite.resources import (
    load_cached_resources,
    prepare_species_resources,
)


def _mock_parsed_go() -> ParsedGO:
    nodes = {
        "GO:0001": GOTermNode("GO:0001", "root", "bp"),
        "GO:0002": GOTermNode("GO:0002", "childA", "bp", parents=("GO:0001",)),
        "GO:0003": GOTermNode("GO:0003", "childB", "bp", parents=("GO:0001",)),
    }
    nodes["GO:0001"] = nodes["GO:0001"].with_children(("GO:0002", "GO:0003")).with_depth(0)
    nodes["GO:0002"] = nodes["GO:0002"].with_depth(1)
    nodes["GO:0003"] = nodes["GO:0003"].with_depth(1)
    tree = GOTree(nodes=nodes, roots=("GO:0001",))
    term_to_genes = {
        "GO:0001": {"A", "B", "C", "D"},
        "GO:0002": {"A", "B"},
        "GO:0003": {"C", "D"},
    }
    for term_id, genes in term_to_genes.items():
        tree.nodes[term_id] = tree.nodes[term_id].with_genes(genes)
    return ParsedGO(tree=tree, term_to_genes=term_to_genes)


def test_prioritize_terms_blocks_child_terms():
    parsed = _mock_parsed_go()
    results = [
        EnrichmentResult(term_id="GO:0001", z_score=5.0, p_value=1e-4, fdr=0.001, overlap=4, total_genes=4, background=100),
        EnrichmentResult(term_id="GO:0002", z_score=3.0, p_value=1e-3, fdr=0.01, overlap=2, total_genes=2, background=100),
        EnrichmentResult(term_id="GO:0003", z_score=2.5, p_value=5e-3, fdr=0.02, overlap=2, total_genes=2, background=100),
    ]
    prioritised = prioritize_terms(
        parsed.tree,
        results,
        PrioritizationSettings(min_z=2.0, max_fdr=0.05, min_overlap=1, delta_z=0.5),
    )

    selected = {res.term_id for res in prioritised if res.selected}
    blocked = {res.term_id: res.blocked_by for res in prioritised if res.blocked_by}

    assert "GO:0001" in selected
    assert blocked["GO:0002"] == "GO:0001"
    assert blocked["GO:0003"] == "GO:0001"


def test_runner_enrichment_scores_terms():
    parsed = _mock_parsed_go()
    runner = GOEliteRunner(parsed, settings=EnrichmentSettings(min_term_size=1, max_term_size=10))
    results = runner.run({"A", "B"}, {"A", "B", "C", "D"}, apply_prioritization=False)
    term_ids = {res.term_id for res in results}
    assert {"GO:0001", "GO:0002"}.issubset(term_ids)
    z_map = {res.term_id: res.z_score for res in results}
    assert z_map["GO:0002"] > z_map["GO:0001"]


def _write_minimal_go_files(base: Path) -> tuple[Path, Path]:
    obo_text = """format-version: 1.2
data-version: releases/2024-01-01

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: child process
namespace: biological_process
is_a: GO:0000001 ! root process
"""
    obo_path = base / "mini.obo"
    obo_path.write_text(obo_text, encoding="utf-8")

    # Construct minimal GAF rows (17 columns)
    row_template = [
        "UniProtKB",
        "P12345",
        "GENEA",
        "",
        "GO:0000001",
        "PMID:1",
        "IDA",
        "",
        "P",
        "Gene A",
        "",
        "protein",
        "taxon:9606",
        "20240101",
        "GO_Curator",
        "",
        "",
    ]
    row_child = row_template.copy()
    row_child[4] = "GO:0000002"
    row_child[2] = "GENEA"

    gaf_lines = [
        "!gaf-version: 2.2",
        "\t".join(row_template),
        "\t".join(row_child),
    ]
    gaf_path = base / "mini.gaf"
    gaf_path.write_text("\n".join(gaf_lines) + "\n", encoding="utf-8")
    return obo_path, gaf_path


def test_resource_cache_roundtrip(tmp_path):
    obo_path, gaf_path = _write_minimal_go_files(tmp_path)
    cache_dir = tmp_path / "cache"
    parsed = prepare_species_resources(
        "human",
        cache_dir=str(cache_dir),
        version="testrun",
        obo_path=str(obo_path),
        gaf_path=str(gaf_path),
        force=True,
    )
    assert "GO:0000001" in parsed.tree.nodes
    assert parsed.term_to_genes["GO:0000002"] == {"GENEA"}

    # Remove source files to ensure cache is used on the next call
    obo_path.unlink()
    gaf_path.unlink()

    cached = prepare_species_resources("human", cache_dir=str(cache_dir), version="testrun")
    assert "GO:0000002" in cached.tree.nodes
    assert cached.term_to_genes["GO:0000002"] == {"GENEA"}

    loaded = load_cached_resources("human", cache_dir=str(cache_dir), version="testrun")
    assert loaded.term_to_genes["GO:0000002"] == {"GENEA"}
