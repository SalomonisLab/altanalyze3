"""Indexed gene lookup for junctions whose two splice sites both miss the reference.

``annotate_junctions`` previously scanned all 68,008 genes for every such junction, at 3.40 ms per
call on the plus strand and 14.78 ms on the minus strand. This module answers the same question by
one ``bisect`` into a precomputed sweep, at 0.24 us per call.

Two modes:

``legacy``
    Reproduces the previous scan exactly, including two behaviours that are almost certainly wrong.
    The scan never tests the chromosome, so a chr2 position can return a chr18 gene. The scan also
    uses ``geneData[gene][0][0] <= pos <= geneData[gene][-1][1]``, and the loader stores minus-strand
    coordinates reversed, so that test reads ``high <= pos <= low`` and holds for 0 of 32,995
    minus-strand genes. Default, so existing output does not move.

``corrected``
    Tests the chromosome, and spans each gene from its lowest to its highest coordinate, so
    minus-strand genes become reachable. Differs from ``legacy`` on minus-strand junctions, on
    cross-chromosome matches, and on 376 plus-strand genes whose exon blocks overlap.

Both modes break ties the same way: the gene that appears first in ``geneData`` insertion order,
which is first-appearance order in the exon reference file. ``legacy`` therefore returns the same
gene the loop returned.
"""

import bisect
import heapq


LEGACY = "legacy"
CORRECTED = "corrected"


def _sweep(intervals):
    """Map every coordinate to the lowest-rank interval covering it.

    ``intervals`` holds ``(lo, hi, rank, gene)`` with ``hi`` inclusive. Returns the segment start
    coordinates and the winning gene per segment, both ordered by coordinate.
    """
    events = []
    for lo, hi, rank, gene in intervals:
        events.append((lo, 0, rank, gene))
        events.append((hi + 1, 1, rank, gene))
    events.sort(key=lambda e: (e[0], e[1]))

    starts, winners = [], []
    heap, alive = [], {}
    i = 0
    while i < len(events):
        coord = events[i][0]
        while i < len(events) and events[i][0] == coord:
            _, kind, rank, gene = events[i]
            if kind == 0:
                heapq.heappush(heap, (rank, gene))
                alive[rank] = alive.get(rank, 0) + 1
            else:
                alive[rank] = alive.get(rank, 0) - 1
            i += 1
        while heap and alive.get(heap[0][0], 0) <= 0:
            heapq.heappop(heap)
        starts.append(coord)
        winners.append(heap[0][1] if heap else None)
    return starts, winners


class NovelGeneIndex(object):
    """Answer "which gene covers this position" without scanning every gene."""

    def __init__(self, geneData, strandData, gene_chr=None, mode=LEGACY):
        if mode not in (LEGACY, CORRECTED):
            raise ValueError(f"mode must be {LEGACY!r} or {CORRECTED!r}, got {mode!r}")
        if mode == CORRECTED and not gene_chr:
            raise ValueError(
                "corrected mode needs gene_chr; call gff_process.importEnsemblGenes first and pass "
                "gff_process.gene_chr"
            )
        self.mode = mode
        self._index = {}

        buckets = {}
        for rank, gene in enumerate(geneData):
            exons = geneData[gene]
            strand = strandData[gene]
            if mode == LEGACY:
                lo, hi = exons[0][0], exons[-1][1]
                if lo > hi:
                    continue                      # empty span, covers no position
                key = strand
            else:
                lo = min(min(a, b) for a, b, _ in exons)
                hi = max(max(a, b) for a, b, _ in exons)
                key = (gene_chr[gene], strand)
            buckets.setdefault(key, []).append((lo, hi, rank, gene))

        for key, intervals in buckets.items():
            self._index[key] = _sweep(intervals)

        self.n_genes_indexed = sum(len(v) for v in buckets.values())
        self.n_genes_total = len(geneData)

    def find(self, chrom, pos, strand):
        """Return the covering gene, or None. ``legacy`` mode ignores ``chrom``."""
        key = strand if self.mode == LEGACY else (chrom, strand)
        entry = self._index.get(key)
        if entry is None:
            return None
        starts, winners = entry
        i = bisect.bisect_right(starts, pos) - 1
        if i < 0:
            return None
        return winners[i]
