"""Exact per-gene index for positive splice-span compatibility.

For parent gap (a,b), reject iff a source span (x,y) satisfies both
abs(x-a) <= slop and abs(y-b) <= slop. Binary searches implement exactly
these closed bounds. Inputs remain immutable during one collapse_prepared
call; the checker and its strong identity references are released per gene.
"""
from bisect import bisect_left, bisect_right


def make_checker(merge_exons):
    span_cache = {}
    parent_cache = {}

    def compatible(spans, parent_exons, splice_slop):
        if not spans:
            return True
        parent_id = id(parent_exons)
        entry = parent_cache.get(parent_id)
        if entry is None:
            blocks = merge_exons(parent_exons)
            gaps = tuple((a[1], b[0]) for a, b in zip(blocks, blocks[1:]))
            # Retain the actual object so its identity cannot be recycled.
            parent_cache[parent_id] = (parent_exons, gaps)
        else:
            gaps = entry[1]
        if not gaps:
            return True
        span_id = id(spans)
        entry = span_cache.get(span_id)
        if entry is None:
            grouped = {}
            for left, right in spans:
                grouped.setdefault(left, []).append(right)
            lefts = tuple(sorted(grouped))
            rights = tuple(tuple(sorted(grouped[left])) for left in lefts)
            span_cache[span_id] = (spans, lefts, rights)
        else:
            _, lefts, rights = entry
        for a, b in gaps:
            lo = bisect_left(lefts, a - splice_slop)
            hi = bisect_right(lefts, a + splice_slop)
            lower, upper = b - splice_slop, b + splice_slop
            for i in range(lo, hi):
                values = rights[i]
                j = bisect_left(values, lower)
                if j < len(values) and values[j] <= upper:
                    return False
        return True

    return compatible
