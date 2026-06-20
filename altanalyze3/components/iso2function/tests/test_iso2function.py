"""Unit tests for iso2function. The pure-logic tests need no external data; the integration test is
skipped automatically when the source supplementary zips are not present."""

import os

import pandas as pd
import pytest

from altanalyze3.components.iso2function import config
from altanalyze3.components.iso2function.crosswalk import structure_key, clone_map
from altanalyze3.components.iso2function.network import build
from altanalyze3.components.iso2function.enrichment import isoform_gsea
from altanalyze3.components.iso2function.ingest import paper2
from altanalyze3.components.iso2function.crosswalk import sequence_map


# --------------------------------------------------------------------------- structure_key
def test_normalize_structure():
    assert structure_key.normalize_structure("E1.1|E2.1|I2.1") == ("E1.1", "E2.1", "I2.1")
    assert structure_key.normalize_structure("  E1.1 | E2.1 |") == ("E1.1", "E2.1")
    assert structure_key.normalize_structure("") == ()
    assert structure_key.normalize_structure(None) == ()


def test_match_structures_exact_and_containment():
    a = "E2.1|E3.1|E4.1"
    assert structure_key.match_structures(a, a) == "exact"
    # a is a contiguous block of the longer b -> contained_in
    b = "E1.1|E2.1|E3.1|E4.1|E5.1"
    assert structure_key.match_structures(a, b) == "contained_in"
    # reversed direction -> superset_of
    assert structure_key.match_structures(b, a) == "superset_of"
    # non-contiguous overlap -> none
    assert structure_key.match_structures("E1.1|E5.1", b) == "none"


def test_is_assigned():
    assert structure_key.is_assigned("exact")
    assert structure_key.is_assigned("contained_in")
    assert not structure_key.is_assigned("superset_of")
    assert not structure_key.is_assigned("novel")


# --------------------------------------------------------------------------- clone id parsing
def test_parse_gff_clone_id():
    p = clone_map.parse_gff_clone_id("AEBP2|1/3|11E07")
    assert p == {"gene_symbol": "AEBP2", "iso_index": 1, "iso_total": 3, "well": "11E07"}
    assert clone_map.parse_gff_clone_id("garbage")["gene_symbol"] is None


def test_parse_paper2_clone_id():
    assert clone_map.parse_paper2_clone_id("AEBP2-2") == {"gene_symbol": "AEBP2", "iso_num": 2}
    # trailing -N absent -> whole string is the gene, iso_num None
    assert clone_map.parse_paper2_clone_id("FOO")["iso_num"] is None


# --------------------------------------------------------------------------- differential interactome
def _toy_edges():
    # up isoform U detects {P1,P2}; down isoform D detects {P2,P3}; all co-tested
    rows = []
    for src, partners in (("U", {"P1": True, "P2": True, "P3": False}),
                          ("D", {"P1": False, "P2": True, "P3": True})):
        for p, det in partners.items():
            rows.append({"source": src, "target": p, "edge_type": "ppi", "assay": "y2h",
                         "detected": det, "tested": True, "score": ""})
    return pd.DataFrame(rows)


def test_differential_interactome_respects_tested_space():
    d = build.differential_interactome("U", "D", _toy_edges(), edge_type="ppi")
    assert d["gained"] == ["P1"]   # positive in U, tested-negative in D
    assert d["lost"] == ["P3"]     # positive in D, tested-negative in U
    assert d["shared"] == ["P2"]
    assert set(d["co_tested_targets"]) == {"P1", "P2", "P3"}


# --------------------------------------------------------------------------- enrichment
def test_hypergeometric_and_bh():
    background = [f"G{i}" for i in range(100)]
    gene_set = {"G0", "G1", "G2", "G3", "G4"}      # K=5
    query = ["G0", "G1", "G2", "G90", "G91"]        # n=5, k=3 overlap
    res = isoform_gsea.hypergeometric_test(query, gene_set, background)
    assert res["N"] == 100 and res["K"] == 5 and res["n"] == 5 and res["k"] == 3
    assert 0.0 <= res["p_value"] <= 1.0
    assert res["fold_enrichment"] > 1  # 3 observed vs 0.25 expected

    df = isoform_gsea.enrich(query, {"setA": gene_set, "setB": {"G50", "G51"}}, background)
    assert "q_value" in df.columns
    assert (df["q_value"].dropna().between(0, 1)).all()


def test_bh_monotonic():
    q = isoform_gsea.benjamini_hochberg([0.01, 0.02, 0.03, 0.5])
    assert all(q[i] <= q[i + 1] + 1e-9 for i in range(len(q) - 1))


# --------------------------------------------------------------------------- thresholds (binary-only)
def test_m1h_call_uses_paper_cutoff():
    # verbatim paper rule: M1H signal >= 1 -> activator, <= -1 -> repressor, else neither
    assert paper2._m1h_call(1.0) == "activator"
    assert paper2._m1h_call(2.5) == "activator"
    assert paper2._m1h_call(-1.0) == "repressor"
    assert paper2._m1h_call(0.3) == "neither"
    assert paper2._m1h_call(float("nan")) == ""
    assert paper2._m1h_call(None) == ""
    # constants are the documented values, not invented
    assert config.M1H_ACTIVATOR_MIN == 1.0 and config.M1H_REPRESSOR_MAX == -1.0


def test_n2h_excluded_from_called_edges():
    # binary-only policy: N2H has no numeric cutoff in metadata -> never a called edge
    nodes, edges = build.build_graph()
    if len(edges):
        assert "n2h" not in set(edges["assay"].astype(str))
        assert set(edges["assay"].astype(str)) <= {"y2h", "ey1h"}


# --------------------------------------------------------------------------- sequence crosswalk (paper1/UniProt)
def test_translate_orf():
    # ATG GCA TGA -> M A (stop trimmed)
    assert sequence_map.translate_orf("ATGGCATGA") == "MA"
    assert sequence_map.translate_orf("") == ""
    # lowercase + trailing partial codon handled
    assert sequence_map.translate_orf("atggca tg".replace(" ", "")) == "MA"


def test_pct_identity_is_gap_aware():
    # an N-terminal offset (indel) must NOT collapse identity to ~0 (the old ungapped bug).
    # "MAAAAK" vs "XMAAAAK" share 6/7 -> difflib ratio ~0.92; ungapped-from-N-term would be ~0.
    assert sequence_map.pct_identity("MAAAAK", "XMAAAAK") > 85
    assert sequence_map.pct_identity("MAAAAK", "MAAAAK") == 100.0
    assert sequence_map.pct_identity("MAAAAK", "WWWWWW") < 30


def test_resolve_structure_tries_all_ensts():
    from altanalyze3.components.iso2function.crosswalk import reference_structures
    e2s = {"ENST_OLD": {"structure": "Sx", "gene": "ENSG9"}}
    # first ENST (newest build) is absent; the second one carries the structure -> must be found
    enst, struct, ensg = reference_structures.resolve_structure(["ENST_NEW", "ENST_OLD"], e2s)
    assert (enst, struct, ensg) == ("ENST_OLD", "Sx", "ENSG9")
    # none present -> first enst, empty structure
    assert reference_structures.resolve_structure(["ENST_A"], e2s) == ("ENST_A", "", "")


def test_sequence_crosswalk_exact_and_fuzzy():
    target = sequence_map.build_target_index_from_records([
        ("MAAAK", {"transcript_id": "ENST_A", "structure": "Sa", "enst": "ENST_A"}),
        ("MBBBK", {"transcript_id": "ENST_B", "structure": "Sb", "enst": "ENST_B"}),
    ])
    # exact match
    meta, pid = sequence_map.match_one("MAAAK", target)
    assert meta["structure"] == "Sa" and pid == 100.0
    # no exact, exact-only mode -> unresolved
    meta, pid = sequence_map.match_one("MAAAR", target)
    assert meta is None
    # fuzzy within allowed same-gene targets
    meta, pid = sequence_map.match_one("MAAAR", target, min_identity=80.0, allowed_targets=["MAAAK"])
    assert meta["structure"] == "Sa" and pid >= 80.0

    res = sequence_map.crosswalk_sequences({"q1": "MAAAK", "q2": "ZZZZ"}, target)
    by = {r["query_id"]: r for r in res}
    assert by["q1"]["matched"] and by["q1"]["structure"] == "Sa"
    assert not by["q2"]["matched"]


# --------------------------------------------------------------------------- Leg A resolver (gene+index)
def test_resolve_clone_to_structure_joins_by_gene_and_index():
    struct = pd.DataFrame([
        {"gene_symbol": "AEBP2", "ensg": "ENSG1", "iso_index": 1, "structure": "Sa",
         "match_type": "exact", "assigned": True, "final_isoform_id": "ENST1", "final_known": "known",
         "enst": "ENST1", "protein_length": "503", "nmd_status": "Not-NMD", "gff_clone_id": "AEBP2|1/3|x"},
        {"gene_symbol": "AEBP2", "ensg": "ENSG1", "iso_index": 2, "structure": "Sb",
         "match_type": "exact", "assigned": True, "final_isoform_id": "ENST2", "final_known": "known",
         "enst": "ENST2", "protein_length": "268", "nmd_status": "Not-NMD", "gff_clone_id": "AEBP2|2/3|y"},
    ])
    clones = pd.DataFrame([
        {"gene_symbol": "AEBP2", "clone_id": "AEBP2-2", "iso_num": 2, "sources": "ppi_y2h"},  # -> Sb
        {"gene_symbol": "AEBP2", "clone_id": "AEBP2-9", "iso_num": 9, "sources": "ppi_y2h"},  # idx absent
        {"gene_symbol": "ZZZ",   "clone_id": "ZZZ-1",   "iso_num": 1, "sources": "ppi_y2h"},  # gene absent
    ])
    resolved, unresolved = clone_map.resolve_clone_to_structure(struct, clones)
    assert len(resolved) == 1
    row = resolved.iloc[0]
    assert row["clone_id"] == "AEBP2-2" and row["structure"] == "Sb" and row["enst"] == "ENST2"
    reasons = set(unresolved["reason"])
    assert any("not present" in r for r in reasons)      # AEBP2-9 index out of range
    assert any("not in TFIso bridge" in r for r in reasons)  # ZZZ gene absent


# --------------------------------------------------------------------------- integration (guarded)
def _zips_present():
    try:
        return all(os.path.exists(config.resolve_shared(os.path.join(config.PAPER2_DIR, z)))
                   for z, _, _ in config.PAPER2_TABLES.values())
    except Exception:
        return False


@pytest.mark.skipif(not _zips_present(), reason="paper2 supplementary zips not available")
def test_ingest_qc_passes(tmp_path):
    from altanalyze3.components.iso2function.ingest import paper2
    manifest = paper2.parse_all(out_dir=str(tmp_path))
    assert manifest["qc_ok"].all(), manifest[~manifest["qc_ok"]].to_string(index=False)


def _pipeline_outputs_present():
    d = config.DATA_DIR
    return all(os.path.exists(os.path.join(d, f)) for f in
               ("switch_function.tsv", "clone_to_structure.tsv", "ppi_y2h.tsv"))


@pytest.mark.skipif(not _pipeline_outputs_present(),
                    reason="run `iso2function all` first (needs data/ outputs)")
def test_my_differential_agrees_with_authors_ppi_category():
    # Validation: an independent binary gained/lost computation should match the authors' Data_S8
    # PPI_category direction -- "PPI loss" switches lose more than they gain; "no PPI change" ~0.
    from altanalyze3.components.iso2function.enrichment import switch_enrichment
    df = switch_enrichment.differential_interactions_by_switch()
    df = df[df["n_ppi_cotested"] > 0]
    loss = df[df["PPI_category"].astype(str).str.startswith("PPI loss")]
    nochange = df[df["PPI_category"].astype(str).str.startswith("no PPI change")]
    if len(loss):
        assert loss["n_ppi_lost"].mean() > loss["n_ppi_gained"].mean()
    if len(nochange):
        assert nochange["n_ppi_lost"].mean() < 1.0  # near-zero turnover
