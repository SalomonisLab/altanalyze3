"""Post-collapse consolidation of exact physical splice-junction chains.

This operates on an existing, fully reconciled retained catalog. It does not
repeat source preparation, infer missing junctions, change splice coordinates,
or rescue reads excluded by the original catalog. Non-splice terminal extents
are deliberately ignored. Consequently a pooled abundance describes a splice
structure, not necessarily a unique full-length RNA or protein sequence.

The caller supplies the frozen native parent-ranking inputs; catalog totals are
not a substitute for a novel representative's original own-class support.
Every input identity remains in the returned crosswalk. Keep the original
geometry/sequence/protein tables alongside that crosswalk as provenance.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from numbers import Integral


def physical_blocks(exons):
    """Normalize inclusive, ascending genomic exon blocks without inventing gaps."""
    result = []
    for a, b in sorted(exons):
        if not isinstance(a, Integral) or not isinstance(b, Integral) or a > b:
            raise ValueError('Verified ascending integer exon coordinates required')
        if result and a <= result[-1][1] + 1:
            result[-1] = (result[-1][0], max(result[-1][1], b))
        else:
            result.append((int(a), int(b)))
    if not result:
        raise ValueError('Every catalog identity requires physical exon geometry')
    return tuple(result)


def junction_chain(exons):
    """Exact genomic donor/acceptor boundaries; outer endpoints are excluded."""
    blocks = physical_blocks(exons)
    return tuple((a[1], b[0]) for a, b in zip(blocks, blocks[1:]))


@dataclass(frozen=True)
class CatalogModel:
    gene: str
    transcript: str
    chromosome: str
    strand: str
    exons: tuple
    total_reads: int
    known: bool
    # Known: (-native full-token length, native reference cache key).
    # Novel: (-native own-class support, complete evidence-class key).
    native_rank: tuple

    @property
    def identity(self):
        return self.gene, self.transcript


def plan_cleanup(models, *, expected_identities):
    """Return a complete direct crosswalk and summed retained-catalog counts.

    All already-retained known IDs survive. Novel matches prefer the first
    native-ranked known ID; novel-only groups use native own-class ranking.
    Single-exon models always survive independently. There is no approximate
    matching, fragment containment, or transitive closure.
    """
    rows = list(models)
    expected = list(expected_identities)
    identities = [r.identity for r in rows]
    if len(set(expected)) != len(expected) or len(set(identities)) != len(identities):
        raise ValueError('Duplicate required or candidate catalog identities')
    if set(identities) != set(expected):
        raise ValueError('Complete instructed catalog identity coverage required')
    by_id = {r.identity: r for r in rows}
    groups = defaultdict(list)
    mapping = {k: k for k in identities}
    for r in rows:
        if (not r.gene or not r.transcript or not r.chromosome or r.strand not in ('+', '-')
                or type(r.known) is not bool or not isinstance(r.total_reads, Integral)
                or r.total_reads <= 0):
            raise ValueError('Verified positive-support catalog metadata required')
        if (len(r.native_rank) != 2 or not isinstance(r.native_rank[0], Integral)
                or r.native_rank[0] > 0 or not isinstance(r.native_rank[1], str)
                or not r.native_rank[1]):
            raise ValueError('Explicit frozen native ranking inputs required')
        chain = junction_chain(r.exons)
        if chain:
            groups[r.gene, r.chromosome, r.strand, chain].append(r)
    merged_groups = []
    for shape, members in groups.items():
        known = [r for r in members if r.known]
        # Exact ties would make source identity depend on input row order.
        choices = known or members
        ranks = [r.native_rank for r in choices]
        if len(ranks) != len(set(ranks)):
            raise ValueError('Native representative ranking is not uniquely resolved')
        target = min(choices, key=lambda r: r.native_rank)
        changed = []
        for r in members:
            if not r.known and r.identity != target.identity:
                mapping[r.identity] = target.identity
                changed.append(r.identity)
        if changed:
            merged_groups.append(dict(target=target.identity, children=changed,
                                      junctions=shape[-1], known_targets=len(known)))
    totals = Counter()
    for key, target in mapping.items():
        if mapping[target] != target:
            raise ValueError('Cleanup crosswalk must have direct surviving targets')
        if by_id[key].known and target != key:
            raise ValueError('An existing known identity was redirected')
        totals[target] += by_id[key].total_reads
    if sum(totals.values()) != sum(r.total_reads for r in rows):
        raise ValueError('Retained-catalog read counts were not conserved')
    return dict(mapping=mapping, totals=dict(totals), merged_groups=merged_groups,
                all_original_identities=len(rows), retained_identities=len(totals),
                original_known_identities=sum(r.known for r in rows),
                merged_novel_identities=sum(k != v for k, v in mapping.items()),
                original_RNA_support=sum(r.total_reads for r in rows))


def terminal_merge_preserves_evidence(evidence, original_exons, target_exons, *, splice_slop):
    """Check one original source class under the endpoint-insensitive rule.

    The original parent must pass the unchanged native evidence checks. Both
    parents must have the same complete nonempty chain. Recheck every measured
    junction and continuous span and reject any internal covered-gap conflict.
    Only terminal coverage/parent-overlap differences are waived. Do not clip
    measured coverage into artificial short intervals at the terminal splice.
    """
    from ..isoform_collapse_utils import intronic_edge_compatible
    from .splice_evidence import splice_spans_compatible

    original = physical_blocks(original_exons)
    target = physical_blocks(target_exons)
    chain = junction_chain(original)
    if not chain or chain != junction_chain(target):
        raise ValueError('Source evidence waiver requires the same complete nonempty chain')
    junctions = evidence.get('protected_junctions', ())
    coverage = evidence.get('protected_coverage', ())
    spans = evidence.get('protected_splice_spans', ())
    if not (intronic_edge_compatible(junctions, coverage, original, splice_slop=splice_slop)
            and splice_spans_compatible(spans, original, splice_slop)):
        return False
    if not (intronic_edge_compatible(junctions, (), target, splice_slop=splice_slop)
            and splice_spans_compatible(spans, target, splice_slop)):
        return False
    return not any(max(a, left+1+splice_slop) <= min(b, right-1-splice_slop)
                   for a, b in coverage for left, right in chain)


def remap_sparse_counts(matrix, source_ids, crosswalk, *, output_ids=None):
    """Sum raw count columns and verify every ordered observation's exact sum.

    ``source_ids`` is the entire ordered feature roster of this matrix; its
    members must all occur in an already-validated catalog crosswalk. Output
    IDs default to sorted observed targets, as in native WTA quantification.
    CPM, ratios and other transformed values must be recomputed afterwards.
    """
    import numpy as np
    from scipy import sparse

    ids = list(source_ids)
    if not sparse.issparse(matrix) or matrix.shape[1] != len(ids) or len(set(ids)) != len(ids):
        raise ValueError('Full unique ordered sparse-matrix feature roster required')
    if any(k not in crosswalk for k in ids):
        raise ValueError('A required matrix feature has no catalog crosswalk entry')
    targets = [crosswalk[k] for k in ids]
    if any(crosswalk.get(k) != k for k in targets):
        raise ValueError('Every crosswalk target must survive directly')
    output = sorted(set(targets)) if output_ids is None else list(output_ids)
    if len(set(output)) != len(output) or set(output) != set(targets):
        raise ValueError('Complete unique matrix target roster required')
    data = matrix.data
    if (not np.isfinite(data).all() or (data < 0).any()
            or (data >= 2**63).any() or not np.equal(data, np.floor(data)).all()):
        raise ValueError('Finite nonnegative integer raw counts required')
    # Integer arithmetic avoids float32 rounding when merged counts grow.
    counts = matrix.astype(np.int64)
    positions = {k: i for i, k in enumerate(output)}
    projection = sparse.csr_matrix((np.ones(len(ids), dtype=np.int64),
                                    (np.arange(len(ids)), [positions[k] for k in targets])),
                                   shape=(len(ids), len(output)))
    result = (counts @ projection).tocsr()
    result.sum_duplicates()
    if (result.data < 0).any() or not np.array_equal(
            np.asarray(counts.sum(axis=1)), np.asarray(result.sum(axis=1))):
        raise ValueError('Per-observation raw counts were not conserved')
    return result, output
