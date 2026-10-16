"""Exact boundary/query parity with the independent scan-based predicate."""
import itertools
import random
import pytest
from altanalyze3.components.long_read.isoform_collapse.splice_evidence import (
    merge_exons, splice_spans_compatible,
)
from altanalyze3.components.long_read.isoform_collapse.splice_span_index import make_checker


def test_index_matches_scan_on_closed_bounds_and_random_geometry():
    fast = make_checker(merge_exons)
    spans = [{(a, b)} for a in range(-4, 12) for b in range(-4, 12)]
    for a in range(-2, 7):
        for b in range(a + 2, 10):
            parent = [(b, b + 3), (a - 3, a)]
            for source, slop in itertools.product(spans, (0, 1, 3, 1.5)):
                assert fast(source, parent, slop) == splice_spans_compatible(source, parent, slop)
    randomizer = random.Random(20261016)
    for _ in range(5000):
        source = {(randomizer.randrange(-50, 200), randomizer.randrange(-50, 200))
                  for _ in range(randomizer.randrange(1, 60))}
        parent = []
        for _ in range(randomizer.randrange(1, 12)):
            left = randomizer.randrange(-50, 200)
            parent.append((left, left + randomizer.randrange(0, 40)))
        slop = randomizer.choice((0, 1, 2, 3, 4, 1.5))
        assert fast(source, parent, slop) == splice_spans_compatible(source, parent, slop)


def test_empty_evidence_preserves_parent_validation_and_errors():
    fast = make_checker(merge_exons)
    for source in (set(), [], (), frozenset()):
        for parent in ([], [(3, 1)]):
            assert fast(source, parent, 3) == splice_spans_compatible(source, parent, 3)
    for parent in ([], [(3, 1)]):
        with pytest.raises(ValueError) as original:
            splice_spans_compatible({(1, 9)}, parent, 3)
        with pytest.raises(ValueError, match=str(original.value)):
            fast({(1, 9)}, parent, 3)


def test_cache_reuses_geometry_only_within_one_gene_checker():
    merges = []
    def measured(parent):
        merges.append(parent)
        return merge_exons(parent)
    source = ((10, 30),)
    parent = ((0, 13), (33, 40))
    a = make_checker(measured)
    for slop in (0, 3, 4, 0):
        assert a(source, parent, slop) == splice_spans_compatible(source, parent, slop)
    assert len(merges) == 1
    b = make_checker(measured)
    assert b(source, parent, 3) == splice_spans_compatible(source, parent, 3)
    assert len(merges) == 2
