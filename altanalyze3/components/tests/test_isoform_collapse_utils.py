"""Unit tests for the cross-sample isoform collapse token/containment primitives.

Covers the five biological separation cases the unified design mandates
(``isoform-collapse/UNIFIED_isoform_collapse_design.md`` §11):

  1. terminal truncation  -> FOLD
  2. exon skip            -> stay separate
  3. alternative donor    -> stay separate
  4. alternative acceptor -> stay separate
  5. intron retention     -> stay separate

plus tokenizer / containment edge cases and the gene-level fold driver.

The structure strings here use the real AltAnalyze token grammar observed in the
iPSC data, e.g. ``"E3.1|E3.2|E4.2|E4.3|E7.1|E9.1"`` and intron-retention tokens
like ``"I7.1"``.
"""

from __future__ import annotations

from altanalyze3.components.long_read.isoform_collapse_utils import (
    tokenize,
    exon_block_count,
    is_contiguous_subsequence,
    structure_tokens_for_containment,
    contiguous_containment_fold,
    split_token,
    strip_terminal_coords,
    filter_exon_intron_tokens,
)


# ---------------------------------------------------------------------------
# tokenize
# ---------------------------------------------------------------------------
def test_tokenize_basic():
    assert tokenize("E1.1|E2.2|E2.4|") == ["E1.1", "E2.2", "E2.4"]


def test_tokenize_no_trailing_pipe():
    assert tokenize("E1.1|E2.2") == ["E1.1", "E2.2"]


def test_tokenize_empty():
    assert tokenize("") == []
    assert tokenize("|") == []


def test_tokenize_single():
    assert tokenize("E7.1") == ["E7.1"]


# ---------------------------------------------------------------------------
# exon_block_count — length metric for binning (sub-exon versions collapsed)
# ---------------------------------------------------------------------------
def test_exon_block_count_collapses_subexons():
    # E1.1|E1.2|E1.3 are sub-versions of one block E1
    assert exon_block_count("E1.1|E1.2|E1.3|E2.1") == 2


def test_exon_block_count_excludes_introns():
    # I* intron-retention tokens are not exon blocks
    assert exon_block_count("E1.1|E2.1|I2.1|E3.1") == 3


def test_exon_block_count_trans_spliced_prefix():
    # gene-prefixed (trans-spliced) tokens still reduce to their block
    assert exon_block_count("ENSG9:E1.1|ENSG9:E1.2|E2.1") == 2


def test_exon_block_count_empty():
    assert exon_block_count("") == 0


def test_exon_block_count_distinct_blocks():
    assert exon_block_count("E1.1|E2.1|E3.1|E4.1|E5.1") == 5


def test_exon_block_count_is_length_metric_not_identity():
    # Two DISTINCT sub-exon structures can share the same block count; the metric
    # is intentionally lossy for *length* only — they remain distinct strings.
    a = "E1.1|E2.1"
    b = "E1.2|E2.2"
    assert a != b
    assert exon_block_count(a) == exon_block_count(b) == 2


# ---------------------------------------------------------------------------
# is_contiguous_subsequence — the core test
# ---------------------------------------------------------------------------
PARENT = ["E1.1", "E2.1", "E3.1", "E4.1", "E5.1"]


def test_contiguous_prefix_folds():
    # 3' truncation keeps a contiguous prefix
    assert is_contiguous_subsequence(["E1.1", "E2.1", "E3.1"], PARENT)


def test_contiguous_suffix_folds():
    # 5' truncation keeps a contiguous suffix
    assert is_contiguous_subsequence(["E3.1", "E4.1", "E5.1"], PARENT)


def test_contiguous_interior_folds():
    # both ends truncated -> interior contiguous run
    assert is_contiguous_subsequence(["E2.1", "E3.1", "E4.1"], PARENT)


def test_exon_skip_does_not_fold():
    # E2.1 missing in the middle -> ordered but NON-contiguous -> must NOT fold
    assert not is_contiguous_subsequence(["E1.1", "E3.1", "E4.1", "E5.1"], PARENT)


def test_full_equality_is_contiguous_but_strictness_handled_by_caller():
    # The window test itself returns True for an equal sequence; strict
    # containment (pn > cn) is enforced by the fold driver, not here.
    assert is_contiguous_subsequence(PARENT, PARENT)


def test_longer_than_parent_never_contained():
    assert not is_contiguous_subsequence(PARENT + ["E6.1"], PARENT)


def test_empty_fragment_not_contained():
    assert not is_contiguous_subsequence([], PARENT)


def test_token_boundary_no_false_substring_match():
    # The legacy string bug: "E1|E12" is a *string* substring of "E1|E12|E23"
    # but here we operate on tokens so E1.1 != E1.10 etc. Confirm no spurious
    # match when a token is a string-prefix of another token.
    parent = ["E1.1", "E12.1", "E23.1"]
    assert not is_contiguous_subsequence(["E1.1", "E1"], parent)
    assert is_contiguous_subsequence(["E12.1", "E23.1"], parent)


# ---------------------------------------------------------------------------
# structure_tokens_for_containment — interior novel sites & IR handling
# ---------------------------------------------------------------------------
def test_interior_novel_exon_site_dropped():
    # A novel E*_coord interior site is normalised away so a read differing only
    # by a novel annotation still compares equal on the shared junctions.
    toks = structure_tokens_for_containment("E4.2|E4.3|E7.1_135203|E9.1")
    assert toks == ["E4.2", "E4.3", "E9.1"]


def test_intron_retention_token_retained():
    # I* tokens are real distinguishing features and MUST be kept.
    toks = structure_tokens_for_containment("E3.1|E3.2|E4.2|E4.3|E7.1|I7.1")
    assert "I7.1" in toks


# ---------------------------------------------------------------------------
# contiguous_containment_fold — the five separation cases at gene level
# ---------------------------------------------------------------------------
def _s(structure, count=1, is_ref=False, **extra):
    d = {"structure": structure, "count": count, "is_ref": is_ref}
    d.update(extra)
    return d


def test_fold_terminal_truncation():
    parent = "E1.1|E2.1|E3.1|E4.1|E5.1"
    frag5 = "E3.1|E4.1|E5.1"   # 5' truncation
    frag3 = "E1.1|E2.1|E3.1"   # 3' truncation
    winners, fold_map = contiguous_containment_fold([
        _s(parent, count=10),
        _s(frag5, count=3),
        _s(frag3, count=2),
    ])
    assert fold_map.get(frag5) == parent
    assert fold_map.get(frag3) == parent
    assert len(winners) == 1
    w = winners[0]
    assert w["structure"] == parent
    # absorbed counts accumulate onto the parent
    assert w["count"] == 10 + 3 + 2
    assert set(w["absorbed"]) == {frag5, frag3}


def test_exon_skip_stays_separate():
    full = "E1.1|E2.1|E3.1|E4.1|E5.1"
    skip = "E1.1|E3.1|E4.1|E5.1"   # E2.1 skipped
    winners, fold_map = contiguous_containment_fold([
        _s(full, count=10),
        _s(skip, count=8),
    ])
    assert skip not in fold_map
    assert {w["structure"] for w in winners} == {full, skip}


def test_alt_donor_stays_separate():
    # Alternative donor = different version of the donor exon (E3.1 vs E3.2);
    # the differing token breaks contiguity with the parent's run.
    iso_a = "E1.1|E2.1|E3.1|E5.1"
    iso_b = "E1.1|E2.1|E3.2|E5.1"
    winners, fold_map = contiguous_containment_fold([
        _s(iso_a, count=5),
        _s(iso_b, count=5),
    ])
    assert not fold_map  # neither folds into the other
    assert len(winners) == 2


def test_alt_acceptor_stays_separate():
    iso_a = "E1.1|E2.1|E3.1|E4.1"
    iso_b = "E1.1|E2.2|E3.1|E4.1"   # alternative acceptor on exon 2
    winners, fold_map = contiguous_containment_fold([
        _s(iso_a, count=5),
        _s(iso_b, count=4),
    ])
    assert not fold_map
    assert len(winners) == 2


def test_intron_retention_stays_separate():
    spliced = "E3.1|E3.2|E4.2|E4.3|E7.1|E9.1"
    retained = "E3.1|E3.2|E4.2|E4.3|E7.1|I7.1"   # IR after E7.1 instead of E9.1
    winners, fold_map = contiguous_containment_fold([
        _s(spliced, count=9),
        _s(retained, count=4),
    ])
    assert retained not in fold_map
    assert spliced not in fold_map
    assert len(winners) == 2


# ---------------------------------------------------------------------------
# fold driver — priority and reference protection
# ---------------------------------------------------------------------------
def test_reference_never_folds_as_fragment():
    # A short reference isoform must survive even when a longer observed
    # structure contains it contiguously.
    long_obs = "E1.1|E2.1|E3.1|E4.1|E5.1"
    short_ref = "E2.1|E3.1|E4.1"
    winners, fold_map = contiguous_containment_fold([
        _s(long_obs, count=50, is_ref=False),
        _s(short_ref, count=2, is_ref=True, ref_id="ENST00000001"),
    ])
    assert short_ref not in fold_map      # ref protected
    structures = {w["structure"] for w in winners}
    assert short_ref in structures and long_obs in structures


def test_fragment_prefers_reference_parent():
    # Two parents contain the fragment contiguously; the reference one wins.
    frag = "E2.1|E3.1"
    ref_parent = "E1.1|E2.1|E3.1|E4.1"
    obs_parent = "E2.1|E3.1|E4.1|E5.1"
    winners, fold_map = contiguous_containment_fold([
        _s(ref_parent, count=5, is_ref=True, ref_id="ENST00000002"),
        _s(obs_parent, count=40, is_ref=False),
        _s(frag, count=3),
    ])
    assert fold_map.get(frag) == ref_parent


def test_fragment_prefers_longer_then_higher_count_when_no_ref():
    frag = "E2.1|E3.1"
    longer = "E1.1|E2.1|E3.1|E4.1"
    shorter_parent = "E2.1|E3.1|E4.1"
    winners, fold_map = contiguous_containment_fold([
        _s(longer, count=5),
        _s(shorter_parent, count=50),
        _s(frag, count=3),
    ])
    # length dominates count in parent priority -> folds into the longer one
    assert fold_map.get(frag) == longer


def test_single_structure_no_fold():
    only = "E1.1|E2.1"
    winners, fold_map = contiguous_containment_fold([_s(only, count=7)])
    assert fold_map == {}
    assert len(winners) == 1
    assert winners[0]["structure"] == only


def test_chained_fragments_resolve_to_top_parent():
    # frag1 ⊂ mid ⊂ top ; both fragments should resolve to a surviving winner,
    # never to a folded intermediate.
    top = "E1.1|E2.1|E3.1|E4.1|E5.1"
    mid = "E2.1|E3.1|E4.1"
    frag = "E3.1|E4.1"
    winners, fold_map = contiguous_containment_fold([
        _s(top, count=20),
        _s(mid, count=5),
        _s(frag, count=2),
    ])
    survivors = {w["structure"] for w in winners}
    # the only surviving winner is top; mid and frag both map to a survivor
    assert top in survivors
    for child in (mid, frag):
        assert fold_map[child] in survivors


# ---------------------------------------------------------------------------
# re-exported token helpers behave as in gff_process
# ---------------------------------------------------------------------------
def test_split_token_with_coord():
    gene, base, coord = split_token("E2.1_135203", default_gene="ENSG1")
    assert base == "E2.1" and coord == 135203 and gene == "ENSG1"


def test_split_token_trans_spliced_prefix():
    gene, base, coord = split_token("ENSG9:E1.1", default_gene="ENSG1")
    assert gene == "ENSG9" and base == "E1.1" and coord is None


def test_strip_terminal_coords_drops_terminal_novel_only():
    toks = ["E2.1_100", "E3.1", "E4.1_200"]
    assert strip_terminal_coords(toks, "ENSG1") == ["E3.1"]


def test_strip_terminal_coords_keeps_interior_novel():
    toks = ["E2.1", "E3.1_150", "E4.1"]
    # interior novel coord is NOT a terminal -> retained
    assert strip_terminal_coords(toks, "ENSG1") == ["E2.1", "E3.1_150", "E4.1"]


def test_filter_exon_intron_drops_novel_exon_keeps_intron():
    toks = ["E1.1", "E2.1_100", "I3.1"]
    out = filter_exon_intron_tokens(toks, "ENSG1")
    assert out == ["E1.1", "I3.1"]
