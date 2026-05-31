"""Scored per-sample isoform collapse (Stage 1 core).

Algorithm (per gene), confirmed with the user on the ITGA2B exemplar:

  1. Take the gene's distinct exon structures (terminal-stripped AltAnalyze strings) with read
     counts. ALL structures are kept (no per-sample read-count filter) -- low-read/singleton
     removal happens later, cross-sample (Stage 2), since a structure rare here may be abundant
     in another sample.
  2. LONG BIN = structures whose exon-BLOCK count > 0.8 * (cluster max blocks). Every long-bin
     structure is a DISTINCT representative isoform -- long isoforms NEVER collapse into each other.
  3. score(s) = (blocks(s) / max_blocks) * reads(s).  A non-long structure collapses into the
     HIGHEST-SCORING long-bin isoform it is a perfect CONTIGUOUS SUBSTRING of (each child -> exactly
     ONE parent).
  4. OTHER BIN = the structures not consumed by the long bin (neither a long rep nor collapsed in).
     Apply the SAME scoring (no 0.8 gate) within the other bin: process by score desc; a structure is
     a representative unless it is a perfect contiguous substring of a higher-scoring other-bin
     structure, in which case it collapses there (again exactly one parent).
  5. A structure that is a substring of nothing remains its own isoform.

Hard invariants (enforced):
  * one-parent: every structure maps to exactly one final exemplar (itself if a representative).
  * perfect-substring-only: a collapse edge exists only for a true contiguous token-subsequence.
  * representatives never collapse into each other within a bin.

Exon BLOCKS (E1.1|E1.2 -> E1) are used ONLY to size bins / score; collapse identity uses FULL
structures. No Ensembl needed at this (per-sample) stage. See
[[isoform_longbin_algorithm]] / the ITGA2B worked example.

Efficiency: candidate parents for a structure are found via a token-posting index (a structure can
only be a substring of structures that contain ALL its tokens), so the contiguous-subsequence test
runs only against a few candidates instead of all n -- avoids the O(n^2) brute force.
"""

from __future__ import annotations

import collections
from typing import Dict, List, Sequence, Tuple

try:  # works as a package import
    from ..isoform_collapse_utils import is_contiguous_subsequence
except ImportError:  # works when loaded by file path (hyphenated dir is not a package)
    from altanalyze3.components.long_read.isoform_collapse_utils import is_contiguous_subsequence

LONG_BIN_FRACTION = 0.8
# Stage 1 (per-sample) keeps ALL structures regardless of read count. Low-read / singleton
# filtering is a CROSS-SAMPLE (Stage 2) operation done AFTER collapsing all samples, because a
# structure rare in one sample may be abundant in another. Do NOT drop low-read structures here.
DEFAULT_MIN_READS = 1


def exon_blocks(tokens: Sequence[str]) -> int:
    """Distinct exon blocks (E1.1|E1.2 -> E1; introns excluded)."""
    seen = set()
    for t in tokens:
        base = t.split(':')[-1]
        if base.startswith('E'):
            seen.add(base.split('.')[0].split('_')[0])
    return len(seen)


def _build_token_postings(structs: List[str], tok: Dict[str, tuple]):
    """token -> set(struct indices containing it). Used to prune parent candidates."""
    postings: Dict[str, set] = {}
    for i, s in enumerate(structs):
        for t in set(tok[s]):
            postings.setdefault(t, set()).add(i)
    return postings


def _candidate_supersequences(child_idx, structs, tok, postings, longer_than):
    """Indices of structures that contain ALL of child's tokens and are strictly longer.
    Intersection of postings = necessary condition for being a supersequence."""
    ctoks = tok[structs[child_idx]]
    if not ctoks:
        return []
    cand = None
    for t in ctoks:
        p = postings.get(t)
        if not p:
            return []
        cand = p if cand is None else (cand & p)
        if not cand:
            return []
    return [i for i in cand if i != child_idx and len(tok[structs[i]]) > longer_than]


def _assign_to_reps(child_indices, rep_indices, structs, tok, postings, score, rep_set):
    """Each child -> the highest-score representative it is a perfect contiguous substring of.
    Returns {child_idx: parent_idx}; children that fit no rep are omitted (stay unassigned)."""
    # reps sorted by score desc once
    reps_sorted = sorted(rep_indices, key=lambda i: (-score[i], -len(tok[structs[i]]), structs[i]))
    assigned: Dict[int, int] = {}
    for ci in child_indices:
        clen = len(tok[structs[ci]])
        cand = set(_candidate_supersequences(ci, structs, tok, postings, clen))
        if not cand:
            continue
        for pi in reps_sorted:                 # highest score first
            if pi not in cand or pi not in rep_set:
                continue
            if is_contiguous_subsequence(tok[structs[ci]], tok[structs[pi]]):
                assigned[ci] = pi
                break
    return assigned


def collapse_gene(struct_reads: Dict[str, int], tok: Dict[str, tuple],
                  min_reads: int = DEFAULT_MIN_READS):
    """Collapse one gene's structures. Returns dict with:
        long_reps, other_reps : list[str] (representative structures)
        assignment            : {child_struct: parent_struct}  (one parent each)
        blocks, score         : per-struct
    All structures (>= min_reads) are partitioned: each is either a representative or has exactly
    one parent."""
    structs = [s for s, n in struct_reads.items() if n >= min_reads]
    if not structs:
        return dict(long_reps=[], other_reps=[], assignment={}, blocks={}, score={})

    idx = {s: i for i, s in enumerate(structs)}
    blocks = {i: exon_blocks(tok[structs[i]]) for i in range(len(structs))}
    maxb = max(blocks.values())
    score = {i: (blocks[i] / maxb) * struct_reads[structs[i]] for i in range(len(structs))}
    postings = _build_token_postings(structs, tok)

    # ---- LONG BIN: blocks > 0.8 * maxb ; every long struct is a distinct representative ----
    long_idx = [i for i in range(len(structs)) if blocks[i] > LONG_BIN_FRACTION * maxb]
    long_set = set(long_idx)
    non_long = [i for i in range(len(structs)) if i not in long_set]
    long_assign = _assign_to_reps(non_long, long_idx, structs, tok, postings, score, long_set)

    # ---- OTHER BIN: structures not consumed by the long bin ----
    consumed = long_set | set(long_assign)
    others = [i for i in range(len(structs)) if i not in consumed]
    other_assign: Dict[int, int] = {}
    other_reps: List[int] = []
    # process by score desc; a struct is a rep unless it's a substring of a higher-score other rep
    other_ranked = sorted(others, key=lambda i: (-score[i], -len(tok[structs[i]]), structs[i]))
    claimed = set()
    other_rep_set: set = set()
    for i in other_ranked:
        if i in claimed:
            continue
        other_reps.append(i)
        other_rep_set.add(i)
    # assign each non-rep other struct to its highest-score other-rep parent
    other_children = [i for i in others if i not in other_rep_set]
    other_assign = _assign_to_reps(other_children, other_reps, structs, tok, postings, score, other_rep_set)
    # any other-child that fit no other-rep becomes its own representative (substring of nothing here)
    for i in other_children:
        if i not in other_assign:
            other_reps.append(i)

    assignment = {}
    for ci, pi in long_assign.items():
        assignment[structs[ci]] = structs[pi]
    for ci, pi in other_assign.items():
        assignment[structs[ci]] = structs[pi]

    res = dict(
        long_reps=[structs[i] for i in long_idx],
        other_reps=[structs[i] for i in other_reps],
        assignment=assignment,
        blocks={structs[i]: blocks[i] for i in range(len(structs))},
        score={structs[i]: score[i] for i in range(len(structs))},
        structs=structs,
    )
    _check_invariants(res, struct_reads, tok)
    return res


def _check_invariants(res, struct_reads, tok):
    structs = res['structs']
    reps = set(res['long_reps']) | set(res['other_reps'])
    assignment = res['assignment']
    # 1. partition: every struct is either a rep or has exactly one parent
    for s in structs:
        is_rep = s in reps
        has_parent = s in assignment
        assert is_rep != has_parent, f"struct neither/both rep and child: {s[:40]}"
    # 2. one-parent: assignment maps each child to a single structure
    for c, p in assignment.items():
        assert c != p, "self-parent"
        # 3. perfect-substring: edge is a true contiguous subsequence
        assert is_contiguous_subsequence(tok[c], tok[p]), f"non-substring collapse {c[:30]}->{p[:30]}"
    # 4. reps don't collapse into each other
    for r in reps:
        assert r not in assignment, "representative has a parent"
