"""GO-Elite default build protocol (2026-10-24): union of the GOA GAF and NCBI gene2go, plus the ortholog species'
union carried over through NCBI gene_orthologs, propagated over is_a/part_of/regulates; one source alone when the
other is missing; never a directory removal; the richest cache preferred."""
import gzip
import json

import pytest

from altanalyze3.components.goelite import resources as res
from altanalyze3.components.goelite.parser import (build_tree_with_annotation_frames, parse_gene2go, parse_gene_info, parse_gaf,
                                                   parse_gene_orthologs, transfer_ortholog_annotations, propagate_annotations, parse_obo)

OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: parent
namespace: biological_process
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: child
namespace: biological_process
is_a: GO:0000002 ! parent
"""


def _gaf_line(symbol, go, qualifier="involved_in"):
    return "\t".join(["UniProtKB", "P" + symbol, symbol, qualifier, go, "PMID:1", "IDA", "", "P", symbol, "", "protein",
                      "taxon:9606", "20240101", "UniProt", "", ""])


def _fixture(tmp_path):
    (tmp_path / "go.obo").write_text(OBO)
    (tmp_path / "goa.gaf").write_text("!gaf-version: 2.2\n!date-generated: 2026-10-11\n" + _gaf_line("A", "GO:0000003") + "\n" + _gaf_line("B", "GO:0000002") + "\n")
    with gzip.open(tmp_path / "gene_info.gz", "wt") as fh:
        fh.write("#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n")
        fh.write("9606\t1\tA\t-\tOLDA|AA\t-\n9606\t2\tC\t-\t-\t-\n9606\t3\tD\t-\tOLDD\t-\n10090\t4\tMouse\t-\t-\t-\n")
    with gzip.open(tmp_path / "gene2go.gz", "wt") as fh:
        fh.write("#tax_id\tGeneID\tGO_ID\tEvidence\tQualifier\tGO_term\tPubMed\tCategory\n")
        fh.write("9606\t2\tGO:0000003\tIEA\tinvolved_in\tchild\t-\tProcess\n")       # C: gene2go only
        fh.write("9606\t3\tGO:0000002\tIEA\tNOT involved_in\tparent\t-\tProcess\n")   # D: NOT, excluded
        fh.write("9606\t1\tGO:0000002\tIDA\tinvolved_in\tparent\t-\tProcess\n")       # A again, parent
        fh.write("9606\t99\tGO:0000002\tIEA\tinvolved_in\tparent\t-\tProcess\n")      # unmapped GeneID
        fh.write("10090\t4\tGO:0000003\tIEA\tinvolved_in\tchild\t-\tProcess\n")       # other taxon
    return tmp_path


def test_gene2go_rows_map_through_gene_info_and_exclude_not_and_other_taxa(tmp_path):
    f = _fixture(tmp_path)
    ids, syn = parse_gene_info(str(f / "gene_info.gz"), 9606)
    assert ids == {"1": "A", "2": "C", "3": "D"} and syn == {"OLDA": "A", "AA": "A", "OLDD": "D"}
    df = parse_gene2go(str(f / "gene2go.gz"), ids, 9606)
    assert set(map(tuple, df[["gene_symbol", "go_id"]].values)) == {("C", "GO:0000003"), ("A", "GO:0000002")}
    assert df.attrs["unmapped_gene_ids"] == 1


def test_union_then_propagation(tmp_path):
    f = _fixture(tmp_path)
    ids, _ = parse_gene_info(str(f / "gene_info.gz"), 9606)
    frames = {"goa_gaf": parse_gaf(str(f / "goa.gaf")), "ncbi_gene2go": parse_gene2go(str(f / "gene2go.gz"), ids, 9606)}
    parsed, counts = build_tree_with_annotation_frames(str(f / "go.obo"), frames)
    assert counts == {"goa_gaf": 2, "ncbi_gene2go": 2, "union": 4}      # GAF: A-child, B-parent; gene2go adds C-child, A-parent
    assert parsed.term_to_genes["GO:0000003"] == {"A", "C"}
    assert parsed.term_to_genes["GO:0000002"] == {"A", "B", "C"}
    assert parsed.term_to_genes["GO:0000001"] == {"A", "B", "C"}


def test_build_uses_union_by_default_and_falls_back_to_one_source(tmp_path, monkeypatch):
    f = _fixture(tmp_path)
    root = tmp_path / "cache"
    parsed = res._build_species_resources("human", cache_root=root, user_version=None, obo_path=str(f / "go.obo"),
                                          gaf_path=str(f / "goa.gaf"), acceptable_evidence=None, use_orthologs=False, use_regulates=False,
                                          gene2go_path=str(f / "gene2go.gz"), gene_info_path=str(f / "gene_info.gz"))
    vdir = root / "Hs" / "20261011_union"
    meta = json.load(open(vdir / "metadata.json"))
    assert meta["annotation_protocol"] == "union(goa_gaf,ncbi_gene2go)" and meta["annotation_rows"]["union"] == 4
    assert parsed.term_to_genes["GO:0000002"] == {"A", "B", "C"}
    syn = dict(line.split("\t") for line in gzip.open(vdir / "data" / "synonyms.tsv.gz", "rt").read().strip().split("\n")[1:])
    assert syn["OLDA"] == "A"
    # gene2go unobtainable -> GAF only, recorded, no exception
    parsed2 = res._build_species_resources("human", cache_root=root, user_version=None, obo_path=str(f / "go.obo"),
                                           gaf_path=str(f / "goa.gaf"), acceptable_evidence=None, use_orthologs=False, use_regulates=False,
                                           gene2go_path=str(f / "missing.gz"), gene_info_path=str(f / "gene_info.gz"))
    meta2 = json.load(open(root / "Hs" / "20261011_goa_gaf" / "metadata.json"))
    assert meta2["annotation_protocol"] == "single(goa_gaf)" and meta2["build_problems"]
    assert parsed2.term_to_genes["GO:0000002"] == {"A", "B"}
    # GAF unobtainable -> gene2go only
    parsed3 = res._build_species_resources("human", cache_root=root, user_version=None, obo_path=str(f / "go.obo"),
                                           gaf_path=str(f / "missing.gaf"), acceptable_evidence=None, use_orthologs=False, use_regulates=False,
                                           gene2go_path=str(f / "gene2go.gz"), gene_info_path=str(f / "gene_info.gz"))
    assert parsed3.term_to_genes["GO:0000003"] == {"C"} and parsed3.term_to_genes["GO:0000002"] == {"A", "C"}
    # nothing -> error
    with pytest.raises(RuntimeError):
        res._build_species_resources("human", cache_root=root, user_version=None, obo_path=str(f / "go.obo"),
                                     gaf_path=str(f / "missing.gaf"), acceptable_evidence=None, use_orthologs=False, use_regulates=False,
                                     gene2go_path=str(f / "missing.gz"), gene_info_path=str(f / "gene_info.gz"))
    # rebuilding the same version moves the old one under delete/, never removes it
    res._build_species_resources("human", cache_root=root, user_version=None, obo_path=str(f / "go.obo"),
                                 gaf_path=str(f / "goa.gaf"), acceptable_evidence=None, use_orthologs=False, use_regulates=False,
                                 gene2go_path=str(f / "gene2go.gz"), gene_info_path=str(f / "gene_info.gz"))
    assert any(p.name.startswith("20261011_union_") for p in (root / "Hs" / "delete").iterdir())
    # a union cache is preferred over a GAF-only cache in another root
    monkeypatch.setattr(res, "candidate_cache_dirs", lambda cache_dir=None: (root / "other", root))
    (root / "other" / "Hs").mkdir(parents=True)
    json.dump({"version": "20261011_goa_gaf"}, open(root / "other" / "Hs" / "latest.json", "w"))
    import shutil as _sh
    _sh.copytree(root / "Hs" / "20261011_goa_gaf", root / "other" / "Hs" / "20261011_goa_gaf")
    loaded = res.load_cached_resources("human")
    assert loaded.version_dir.endswith("20261011_union")


OBO_REG = OBO + """
[Term]
id: GO:0000010
name: regulation of child
namespace: biological_process
is_a: GO:0000001 ! root
relationship: negatively_regulates GO:0000003 ! child
"""


def _ortholog_fixture(tmp_path):
    f = _fixture(tmp_path)
    (f / "go.obo").write_text(OBO_REG)
    (f / "goa_mouse.gaf").write_text("!gaf-version: 2.2\n!date-generated: 2026-10-11\n"
                                     + _gaf_line("Mouse", "GO:0000003").replace("taxon:9606", "taxon:10090") + "\n"
                                     + _gaf_line("Mreg", "GO:0000010").replace("taxon:9606", "taxon:10090") + "\n")
    with gzip.open(f / "Mus_musculus.gene_info.gz", "wt") as fh:
        fh.write("#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n10090\t4\tMouse\t-\t-\t-\n10090\t5\tMreg\t-\t-\t-\n10090\t6\tMnone\t-\t-\t-\n")
    with gzip.open(f / "gene_orthologs.gz", "wt") as fh:
        fh.write("#tax_id\tGeneID\trelationship\tOther_tax_id\tOther_GeneID\n")
        fh.write("9606\t3\tOrtholog\t10090\t4\n")      # D <- Mouse (D only had a NOT row: gains child through the ortholog)
        fh.write("10090\t5\tOrtholog\t9606\t2\n")      # C <- Mreg, reversed row direction
        fh.write("9606\t1\tOrtholog\t10090\t999\n")    # unmapped mouse GeneID
        fh.write("9606\t1\tOrtholog\t7227\t50\n")      # other taxon, ignored
    return f


def test_ortholog_parsing_and_transfer(tmp_path):
    f = _ortholog_fixture(tmp_path)
    hid, _ = parse_gene_info(str(f / "gene_info.gz"), 9606)
    mid, _ = parse_gene_info(str(f / "Mus_musculus.gene_info.gz"), 10090)
    orth = parse_gene_orthologs(str(f / "gene_orthologs.gz"), 9606, 10090, hid, mid)
    assert orth.pop("__unmapped__") == {"1"} and orth == {"Mouse": {"D"}, "Mreg": {"C"}}
    mouse = parse_gaf(str(f / "goa_mouse.gaf"))
    moved = transfer_ortholog_annotations(mouse, orth)
    assert set(map(tuple, moved[["gene_symbol", "go_id"]].values)) == {("D", "GO:0000003"), ("C", "GO:0000010")}
    assert moved.attrs["transferred_genes"] == 2


def test_regulates_propagation_is_optional_and_transitive(tmp_path):
    f = _ortholog_fixture(tmp_path)
    tree = parse_obo(str(f / "go.obo"))
    assert tree.nodes["GO:0000010"].regulates == ("GO:0000003",) and tree.parents("GO:0000010") == ("GO:0000001",)
    direct = {"GO:0000010": {"R"}}
    assert propagate_annotations(tree, direct, use_regulates=False) == {"GO:0000010": {"R"}, "GO:0000001": {"R"}}
    with_reg = propagate_annotations(tree, direct, use_regulates=True)
    assert with_reg["GO:0000003"] == {"R"} and with_reg["GO:0000002"] == {"R"}   # regulated term and its is_a ancestors


def test_default_build_transfers_orthologs_and_propagates_regulates(tmp_path, monkeypatch):
    f = _ortholog_fixture(tmp_path)
    root = tmp_path / "cache"
    kw = dict(cache_root=root, user_version=None, obo_path=str(f / "go.obo"), gaf_path=str(f / "goa.gaf"),
              acceptable_evidence=None, gene2go_path=str(f / "gene2go.gz"), gene_info_path=str(f / "gene_info.gz"),
              gene_orthologs_path=str(f / "gene_orthologs.gz"), ortholog_gaf_path=str(f / "goa_mouse.gaf"),
              ortholog_gene_info_path=str(f / "Mus_musculus.gene_info.gz"))
    parsed = res._build_species_resources("human", **kw)
    vdir = root / "Hs" / "20261011_union_orth_reg"
    meta = json.load(open(vdir / "metadata.json"))
    assert meta["annotation_protocol"] == "union(goa_gaf,ncbi_gene2go)+ortholog(mouse:union(goa_gaf,ncbi_gene2go))"
    assert meta["propagation"] == "is_a,part_of,regulates" and meta["build_problems"] == []
    assert meta["annotation_rows"]["ortholog_mouse"] == 2 and meta["ortholog_transfer"]["transferred_genes"] == 2
    # Mouse gene2go row (Mouse -> GO:0000003) and GAF rows both arrive: D gains child, C gains regulation-of-child
    assert parsed.term_to_genes["GO:0000003"] == {"A", "C", "D"}      # C through regulates from GO:0000010
    assert parsed.term_to_genes["GO:0000010"] == {"C"}
    assert parsed.term_to_genes["GO:0000002"] == {"A", "B", "C", "D"}
    # the tree JSON keeps the regulates edges, and the loader re-propagates with the recorded relations
    loaded = res.load_cached_resources("human", cache_dir=str(root))
    assert loaded.tree.nodes["GO:0000010"].regulates == ("GO:0000003",)
    assert loaded.term_to_genes["GO:0000003"] == {"A", "C", "D"}
    # ortholog files unobtainable -> build continues without transfer, records it, no regulates in the name either
    parsed2 = res._build_species_resources("human", **{**kw, "gene_orthologs_path": str(f / "missing_orth.gz"), "use_regulates": False})
    meta2 = json.load(open(root / "Hs" / "20261011_union" / "metadata.json"))
    assert meta2["annotation_protocol"] == "union(goa_gaf,ncbi_gene2go)" and any("ortholog transfer" in p for p in meta2["build_problems"])
    assert parsed2.term_to_genes["GO:0000003"] == {"A", "C"}
    # the ortholog cache is preferred over the plain union cache whatever latest.json says
    json.dump({"version": "20261011_union"}, open(root / "Hs" / "latest.json", "w"))
    assert res.load_cached_resources("human", cache_dir=str(root)).version_dir.endswith("20261011_union_orth_reg")
