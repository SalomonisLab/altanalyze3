"""Tests for the scored isoform-collapse core (isoform_collapse.scored_collapse).

Covers the hard invariants and the key behaviors validated on real ITGA2B data:
  - one-parent partition (every structure is a representative XOR has exactly one parent)
  - perfect-substring-only collapse edges
  - long bin = blocks > 0.8*max; long isoforms never collapse into each other
  - score = (blocks/max)*reads; a fragment collapses into the highest-scoring long isoform it fits
  - an abundant short isoform below the 0.8 gate is NOT a sink (does not swallow full-length forms)
  - a rare-but-long structure does not demote an abundant one in the long bin
"""

from __future__ import annotations

from altanalyze3.components.long_read.isoform_collapse.scored_collapse import (
    collapse_gene, exon_blocks, LONG_BIN_FRACTION,
)


def _tok(structs):
    return {s: tuple(t for t in s.split('|') if t) for s in structs}


def test_exon_blocks_collapses_subexons():
    assert exon_blocks(("E1.1", "E1.2", "E2.1", "I2.1")) == 2  # E1,E2; intron excluded


def test_one_parent_partition_invariant():
    # collapse_gene asserts invariants internally; a clean call must not raise.
    sr = {"E1.1|E2.1|E3.1|E4.1|E5.1": 50, "E2.1|E3.1|E4.1": 5, "E3.1|E4.1": 3}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    reps = set(res['long_reps']) | set(res['other_reps'])
    # partition: each structure is rep xor child
    for s in sr:
        assert (s in reps) != (s in res['assignment'])


def test_truncation_collapses_into_full_length():
    full = "E1.1|E2.1|E3.1|E4.1|E5.1"
    sr = {full: 50, "E2.1|E3.1|E4.1": 5, "E1.1|E2.1|E3.1": 4}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    # full is the only long isoform; the contiguous fragments collapse into it
    assert full in res['long_reps']
    assert res['assignment'].get("E2.1|E3.1|E4.1") == full
    assert res['assignment'].get("E1.1|E2.1|E3.1") == full


def test_exon_skip_stays_separate():
    full = "E1.1|E2.1|E3.1|E4.1|E5.1"
    skip = "E1.1|E3.1|E4.1|E5.1"  # non-contiguous (E2.1 skipped) -> NOT a substring
    sr = {full: 50, skip: 30}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    # both are long (5 vs 4 blocks, both > 0.8*5=4? 4 is not >4) -> at least neither collapses
    assert skip not in res['assignment'] or res['assignment'][skip] != full
    # the exon-skip is never a child of the full (not a contiguous substring)
    assert res['assignment'].get(skip) != full


def test_abundant_short_is_not_a_sink():
    # a 2-block, very abundant structure must NOT swallow a 5-block full-length form.
    full = "E1.1|E2.1|E3.1|E4.1|E5.1"   # 5 blocks
    short = "E1.1|E2.1"                  # 2 blocks, huge reads
    sr = {full: 10, short: 5000}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    # full (5 blocks > 0.8*5) is a long representative; short is below the gate.
    assert full in res['long_reps']
    # full must never collapse INTO the short abundant one.
    assert res['assignment'].get(full) != short


def test_rare_long_does_not_demote_abundant_in_long_bin():
    # two long isoforms; both kept as distinct reps regardless of read counts (no merging).
    a = "E1.1|E2.1|E3.1|E4.1|E5.1"   # abundant
    b = "E1.1|E2.1|E3.1|E4.1|E5.2"   # alt last exon, rare, same length
    sr = {a: 1000, b: 2}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    reps = set(res['long_reps'])
    assert a in reps and b in reps          # both long, neither collapses into the other
    assert a not in res['assignment'] and b not in res['assignment']


def test_score_picks_most_expressed_long_parent():
    # a fragment that fits two long isoforms goes to the higher-scoring (more reads) one.
    p1 = "E1.1|E2.1|E3.1|E4.1|E5.1"   # 5 blocks, 100 reads -> score 100
    p2 = "E1.1|E2.1|E3.1|E4.1|E6.1"   # 5 blocks, 10 reads  -> score 10
    frag = "E1.1|E2.1|E3.1|E4.1"      # contiguous prefix of BOTH
    sr = {p1: 100, p2: 10, frag: 3}
    res = collapse_gene(sr, _tok(sr), min_reads=1)
    assert res['assignment'].get(frag) == p1   # highest score wins
