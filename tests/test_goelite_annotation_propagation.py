"""GO-Elite annotation handling fixed 2026-10-24: true-path propagation to parent terms, exclusion of
NOT-qualified GAF rows, and mapping of outdated query symbols through GAF synonyms."""
import logging

from altanalyze3.components.goelite.parser import (build_tree_with_annotations, parse_gaf, parse_gaf_synonyms,
                                                   propagate_annotations)
from altanalyze3.components.goelite.runner import EnrichmentSettings, GOEliteRunner

OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: regulation of activation
namespace: biological_process
is_a: GO:0000001 ! root process

[Term]
id: GO:0000003
name: regulation of leukocyte activation
namespace: biological_process
is_a: GO:0000002 ! regulation of activation

[Term]
id: GO:0000004
name: part term
namespace: biological_process
relationship: part_of GO:0000003 ! regulation of leukocyte activation
"""


def _gaf_line(symbol, go, qualifier="involved_in", synonyms=""):
    cols = ["UniProtKB", "P" + symbol, symbol, qualifier, go, "PMID:1", "IDA", "", "P", symbol + " name",
            synonyms, "protein", "taxon:9606", "20240101", "UniProt", "", ""]
    return "\t".join(cols)


def _write(tmp_path, gaf_lines):
    obo = tmp_path / "go.obo"; obo.write_text(OBO)
    gaf = tmp_path / "goa.gaf"; gaf.write_text("!gaf-version: 2.2\n" + "\n".join(gaf_lines) + "\n")
    return obo, gaf


def test_direct_annotation_reaches_every_ancestor_through_is_a_and_part_of(tmp_path):
    obo, gaf = _write(tmp_path, [_gaf_line("A", "GO:0000004"), _gaf_line("B", "GO:0000003"), _gaf_line("C", "GO:0000002")])
    parsed = build_tree_with_annotations(str(obo), str(gaf))
    assert parsed.term_to_genes["GO:0000004"] == {"A"}
    assert parsed.term_to_genes["GO:0000003"] == {"A", "B"}
    assert parsed.term_to_genes["GO:0000002"] == {"A", "B", "C"}
    assert parsed.term_to_genes["GO:0000001"] == {"A", "B", "C"}
    # idempotent
    again = propagate_annotations(parsed.tree, parsed.term_to_genes)
    assert again == parsed.term_to_genes


def test_not_qualified_rows_are_excluded(tmp_path):
    obo, gaf = _write(tmp_path, [_gaf_line("A", "GO:0000003"), _gaf_line("B", "GO:0000003", qualifier="NOT|involved_in")])
    df = parse_gaf(str(gaf))
    assert set(df["gene_symbol"]) == {"A"}
    assert set(parse_gaf(str(gaf), exclude_not=False)["gene_symbol"]) == {"A", "B"}


def test_synonyms_map_old_symbols_to_current_ones(tmp_path):
    obo, gaf = _write(tmp_path, [_gaf_line("RO60", "GO:0000003", synonyms="TROVE2|SSA2"),
                                 _gaf_line("MTREX", "GO:0000003", synonyms="SKIV2L2"),
                                 _gaf_line("SKIV2L2", "GO:0000002")])   # a current symbol is never remapped
    syn = parse_gaf_synonyms(str(gaf))
    assert syn["TROVE2"] == "RO60" and syn["SSA2"] == "RO60"
    assert "SKIV2L2" not in syn


def test_propagated_terms_rank_broad_terms_the_hypergeometric_way(tmp_path):
    lines = [_gaf_line(f"L{i}", "GO:0000003") for i in range(10)]
    lines += [_gaf_line(f"P{i}", "GO:0000004") for i in range(5)]
    lines += [_gaf_line(f"O{i}", "GO:0000002") for i in range(30)]
    obo, gaf = _write(tmp_path, lines)
    parsed = build_tree_with_annotations(str(obo), str(gaf))
    runner = GOEliteRunner(parsed, logger=logging.getLogger("t"), settings=EnrichmentSettings(min_term_size=2, max_term_size=2000))
    background = [f"L{i}" for i in range(10)] + [f"P{i}" for i in range(5)] + [f"O{i}" for i in range(30)] + [f"X{i}" for i in range(55)]
    results = runner.run(["L0", "L1", "P0", "P1"], background, apply_prioritization=False)
    by_id = {r.term_id: r for r in results}
    assert by_id["GO:0000003"].overlap == 4 and by_id["GO:0000003"].total_genes == 15   # L + P through part_of
    assert by_id["GO:0000004"].overlap == 2 and by_id["GO:0000004"].total_genes == 5
