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

from ..isoform_collapse_utils import is_contiguous_subsequence, structure_tokens_for_containment

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


def exon_regions(tokens: Sequence[str]) -> int:
    """Distinct exon REGIONS (E6.2 vs E6.3 are distinct; novel _coord suffix stripped;
    introns excluded). Finer-grained than :func:`exon_blocks`, which would count both as E6.
    Used as the region-resolution factor in the score so that, among structures with the same
    block count, the one covering more regions ranks higher (becomes the collapse parent)."""
    seen = set()
    for t in tokens:
        base = t.split(':')[-1]
        if base.startswith('E'):
            seen.add(base.split('_')[0])   # E6.2_56659293 -> E6.2 ; keep block.region
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


def _assign_to_reps_by_length(child_indices, rep_indices, structs, tok, postings, rep_set, prefer=None):
    """Each child -> a representative it is a perfect contiguous substring of. Used by the OTHER bin,
    which collapses by region-length (NOT score). Parent preference order:
        1. an ENST container (indices in `prefer`) -- known transcript wins as parent, so a non-ENST
           that fits inside an ENST folds into the ENST;
        2. otherwise the LONGEST containing rep (region-length), ties broken by structure string.
    Returns {child_idx: parent_idx}; children fitting no rep are omitted."""
    prefer = prefer or set()
    # ENST containers first (longest among them), then non-ENST by length.
    reps_sorted = sorted(rep_indices,
                         key=lambda i: (0 if i in prefer else 1, -len(tok[structs[i]]), structs[i]))
    assigned: Dict[int, int] = {}
    for ci in child_indices:
        clen = len(tok[structs[ci]])
        cand = set(_candidate_supersequences(ci, structs, tok, postings, clen))
        if not cand:
            continue
        for pi in reps_sorted:                 # ENST containers first, then longest
            if pi not in cand or pi not in rep_set:
                continue
            if is_contiguous_subsequence(tok[structs[ci]], tok[structs[pi]]):
                assigned[ci] = pi
                break
    return assigned


def collapse_gene(struct_reads: Dict[str, int], tok: Dict[str, tuple],
                  min_reads: int = DEFAULT_MIN_READS, enst: Dict[str, str] = None):
    """Collapse one gene's structures. Returns dict with:
        long_reps, other_reps : list[str] (representative structures)
        assignment            : {child_struct: parent_struct}  (one parent each)
        blocks, score         : per-struct
    All structures (>= min_reads) are partitioned: each is either a representative or has exactly
    one parent."""
    obs_structs = [s for s, n in struct_reads.items() if n >= min_reads]
    if not obs_structs:
        return dict(long_reps=[], other_reps=[], assignment={}, blocks={}, score={}, structs=[],
                    enst={})

    # ENST INJECTION: inject reference structures ONLY when they belong to this gene's cluster --
    # i.e. an ENST is kept only if some observed structure is contiguity-related to it (the ENST is
    # a superstring of, equal to, or substring of an observed structure). An ENST related to nothing
    # observed is dropped (we only catalog observed isoforms + their known names). Each injected ENST
    # is later given a score strictly greater than any observed long isoform, so within its cluster
    # it wins as the representative (Ensembl-prioritized naming).
    enst = enst or {}
    obs_tok = {s: tok[s] for s in obs_structs}
    kept_enst = {}
    if enst:
        for es, eid in enst.items():
            etk = tok.get(es)
            if etk is None:
                etk = tuple(structure_tokens_for_containment(es))
                tok[es] = etk
            if es in obs_tok:
                kept_enst[es] = eid  # exact match to an observed structure
                continue
            for os_ in obs_structs:                       # related to some observed structure?
                ot = obs_tok[os_]
                if (len(etk) > len(ot) and is_contiguous_subsequence(ot, etk)) or \
                   (len(etk) < len(ot) and is_contiguous_subsequence(etk, ot)):
                    kept_enst[es] = eid
                    break

    structs = list(obs_structs) + [s for s in kept_enst if s not in struct_reads]
    enst_ids = {s: kept_enst[s] for s in structs if s in kept_enst}

    idx = {s: i for i, s in enumerate(structs)}
    blocks = {i: exon_blocks(tok[structs[i]]) for i in range(len(structs))}
    regions = {i: exon_regions(tok[structs[i]]) for i in range(len(structs))}
    # maxb / maxr are computed from OBSERVED (non-ENST) structures ONLY. An injected reference transcript
    # that is not observed must NOT inflate the gene-max block/region count -- otherwise it raises the
    # long-bin gate (0.8*maxb) and demotes real observed structures out of the long bin (a 4-block
    # observed isoform would land in "other" while a 2-block ENST sits in "long"). Reference structures
    # are still forced into the long bin separately; they just don't move the gate for observed data.
    obs_idx = [i for i in range(len(structs)) if structs[i] not in enst_ids]
    maxb = max((blocks[i] for i in obs_idx), default=0) or max(blocks.values())
    maxr = max((regions[i] for i in obs_idx), default=0) or max(regions.values()) or 1
    # Score:
    #   non-ENST: 2 * (blocks/maxb) * reads * (regions/maxr) -- read-weighted; the region factor breaks
    #             ties between structures with the same block count.
    #   ENST:     1_000_000 * (blocks/maxb) * (regions/maxr) -- purely STRUCTURAL (no read term), so a
    #             reference transcript always outranks observed structures and wins as the representative.
    score = {i: (1_000_000.0 * (blocks[i] / maxb) * (regions[i] / maxr)
                 if structs[i] in enst_ids
                 else 2.0 * (blocks[i] / maxb) * struct_reads.get(structs[i], 0) * (regions[i] / maxr))
             for i in range(len(structs))}
    postings = _build_token_postings(structs, tok)

    # ---- LONG BIN: blocks > 0.8 * maxb ----
    # Rule (per user):
    #   * An ENST long structure is ALWAYS a representative -- never folded away.
    #   * A NON-ENST long structure folds into a LONGER long-bin structure it fits inside (the longer
    #     one may be non-ENST or ENST). Score is NOT a gate here: containment + length drives the fold,
    #     so a high-read structure still folds into a longer container it is a substring of.
    #   * When it fits inside SEVERAL longer containers, it folds into the HIGHEST-SCORING one.
    # The parent must be a SURVIVING long representative (ENST, or a non-ENST that itself did not fold).
    # We process by length descending so every potential (longer) parent is decided before its children.
    enst_struct = set(enst_ids)                         # structures that ARE reference transcripts
    # ENSTs are ALWAYS long-bin (regardless of block count); their structural 1,000,000x score governs
    # collapse there so every reference transcript is a long-bin representative.
    long_bin = [i for i in range(len(structs))
                if blocks[i] > LONG_BIN_FRACTION * maxb or structs[i] in enst_struct]
    # longest first; ties by score then string for determinism
    long_ranked = sorted(long_bin, key=lambda i: (-len(tok[structs[i]]), -score[i], structs[i]))
    long_idx: List[int] = []            # surviving long representatives
    long_long_assign: Dict[int, int] = {}
    for i in long_ranked:
        if structs[i] in enst_struct:
            long_idx.append(i)          # ENST always survives as a representative
            continue
        ti = tok[structs[i]]
        # all SURVIVING long reps that are strictly longer and contain i
        containers = [r for r in long_idx
                      if len(tok[structs[r]]) > len(ti) and is_contiguous_subsequence(ti, tok[structs[r]])]
        if containers:
            parent = max(containers, key=lambda r: (score[r], len(tok[structs[r]]), structs[r]))
            long_long_assign[i] = parent
        else:
            long_idx.append(i)          # distinct long isoform (substring of no longer rep)
    long_set = set(long_idx)
    # non-long structures (shorts) fold into the surviving long reps as before. ENSTs are EXCLUDED here
    # -- an ENST must never become a child (ENST always survives), so any non-long ENST is handled in the
    # other bin where it is forced to be a representative.
    non_long = [i for i in range(len(structs))
                if i not in long_set and i not in long_long_assign and structs[i] not in enst_struct]
    long_assign = _assign_to_reps(non_long, long_idx, structs, tok, postings, score, long_set)

    # ---- OTHER BIN: structures not consumed by the long bin ----
    consumed = long_set | set(long_long_assign) | set(long_assign)
    others = [i for i in range(len(structs)) if i not in consumed]
    other_assign: Dict[int, int] = {}
    other_reps: List[int] = []
    # OTHER BIN collapses PURELY BY REGION-LENGTH (no score): process LONGEST-first; a struct is a rep
    # UNLESS it is a perfect contiguous subsequence (substring) of an already-chosen LONGER other rep,
    # in which case it is a CHILD of that rep. Score is NOT used here -- a short high-read fragment must
    # still fold into its superset (e.g. HOPX 'E6.2|...|E8.2' folds into 'E6.1|E6.2|...|E8.2' even though
    # the shorter one has more reads). Region-length = number of region-resolution tokens = len(tok[...]).
    # ENSTs ALWAYS survive: an ENST is never folded away in the other bin either (same rule as the long
    # bin). Process ENSTs first so they are guaranteed representatives; a non-ENST that is a contiguous
    # substring of an ENST therefore folds INTO the ENST.
    other_enst = {i for i in others if structs[i] in enst_ids}
    other_ranked = sorted(others, key=lambda i: (0 if i in other_enst else 1,
                                                 -len(tok[structs[i]]), structs[i]))
    claimed = set()
    other_rep_set: set = set()
    for i in other_ranked:
        if i in claimed:
            continue
        if i in other_enst:
            other_reps.append(i)          # ENST always a representative -- never demoted to a child
            other_rep_set.add(i)
            continue
        # non-ENST: it is a child UNLESS it is a contiguous substring of an existing (longer) rep.
        ti = tok[structs[i]]
        is_child = False
        for r in other_reps:
            tr = tok[structs[r]]
            if len(ti) < len(tr) and is_contiguous_subsequence(ti, tr):
                is_child = True
                break
        if is_child:
            claimed.add(i)
            continue
        other_reps.append(i)
        other_rep_set.add(i)
    # assign each non-rep other struct to its containing other-rep parent. Prefer an ENST container
    # (known transcript wins as parent); otherwise the LONGEST containing rep (region-length, not score).
    other_children = [i for i in others if i not in other_rep_set]
    other_assign = _assign_to_reps_by_length(other_children, other_reps, structs, tok, postings,
                                             other_rep_set, prefer=other_enst)
    # any other-child that fit no other-rep becomes its own representative (substring of nothing here)
    for i in other_children:
        if i not in other_assign:
            other_reps.append(i)

    assignment = {}
    for ci, pi in long_long_assign.items():
        assignment[structs[ci]] = structs[pi]
    for ci, pi in long_assign.items():
        assignment[structs[ci]] = structs[pi]
    for ci, pi in other_assign.items():
        assignment[structs[ci]] = structs[pi]

    long_reps_s = [structs[i] for i in long_idx]
    other_reps_s = [structs[i] for i in other_reps]

    # Drop injected ENST representatives that have NO children -- an ENST that clusters (so it was
    # injected) but ends up representing only itself, with no observed structure collapsing into it,
    # is not an observed isoform and must not appear in the catalog. Only CHILDLESS, 0-read ENST reps
    # are dropped, so no observed child is ever orphaned (an ENST with observed children is kept and
    # names those isoforms).
    if enst_ids:
        n_children = collections.Counter(assignment.values())
        drop = {s for s in list(long_reps_s) + list(other_reps_s)
                if s in enst_ids and struct_reads.get(s, 0) == 0 and n_children.get(s, 0) == 0}
        if drop:
            long_reps_s = [s for s in long_reps_s if s not in drop]
            other_reps_s = [s for s in other_reps_s if s not in drop]
            structs = [s for s in structs if s not in drop]
            # dropped ENSTs are childless, so no assignment edge references them; just prune entries.
            assignment = {c: p for c, p in assignment.items() if c not in drop}
            enst_ids = {s: e for s, e in enst_ids.items() if s not in drop}

    res = dict(
        long_reps=long_reps_s,
        other_reps=other_reps_s,
        assignment=assignment,
        blocks={s: exon_blocks(tok[s]) for s in structs},
        score={s: score[idx[s]] for s in structs},
        structs=structs,
        enst=enst_ids,          # {structure: ENST id} for structures that ARE reference transcripts
    )
    _check_invariants(res, struct_reads, tok)
    return res


def _candidate_parents_for_child(ci, structs, tok, postings, rep_indices, rep_set):
    """All representatives (in rep_set) that the child ci is a perfect contiguous subsequence of."""
    clen = len(tok[structs[ci]])
    cand = set(_candidate_supersequences(ci, structs, tok, postings, clen))
    parents = []
    for pi in rep_indices:
        if pi in cand and pi in rep_set and is_contiguous_subsequence(tok[structs[ci]], tok[structs[pi]]):
            parents.append(pi)
    return parents


def collapse_gene_em(struct_reads: Dict[str, int], tok: Dict[str, tuple],
                     min_reads: int = DEFAULT_MIN_READS, enst: Dict[str, str] = None,
                     n_iter: int = 50, tol: float = 1e-4):
    """EM variant of collapse_gene (OPTIONAL alternative to winner-takes-all).

    Same representative SET and the same long-bin / containment rules as ``collapse_gene`` -- the
    ONLY difference is how an ambiguous child substring's reads are distributed among the
    representatives that contain it. Winner-takes-all gives all of a child's reads to its single
    highest-score parent; here we soft-assign fractionally, proportional to the parents' current read
    abundance, iterating to convergence (an RSEM-style EM):

        E-step: child c's reads split across its compatible parents p in proportion to abundance(p).
        M-step: abundance(p) = own_reads(p) + sum over children of their fraction assigned to p.

    Per your spec, if A=200, B=100, C=50 reads all contain a child's substring, the child's reads go
    mostly to A, less to B, least to C (200:100:50), and the proportions are refined as abundances
    update. A child compatible with exactly one parent contributes all its reads to that parent
    (identical to WTA). Efficient: parents-per-child are precomputed once; each iteration is O(edges).

    Returns the same dict as collapse_gene, PLUS ``em_reads`` = {rep_structure: EM-estimated reads}.
    The catalog/quantification can use ``em_reads`` in place of the hard ``exemplar_total`` when the
    EM option is selected. The hard ``assignment`` (argmax parent) is still returned so the partition
    invariants and the h5ad re-key remain well-defined.
    """
    # Reuse the WTA collapse to get the rep set, bins, scores, ENST handling, and the hard assignment.
    res = collapse_gene(struct_reads, tok, min_reads=min_reads, enst=enst)
    structs = res['structs']
    reps = list(set(res['long_reps']) | set(res['other_reps']))
    if not structs:
        res['em_reads'] = {}
        return res

    idx = {s: i for i, s in enumerate(structs)}
    rep_set = {idx[s] for s in reps}
    rep_indices = list(rep_set)
    postings = _build_token_postings(structs, tok)

    own = {i: float(struct_reads.get(structs[i], 0)) for i in range(len(structs))}

    # children = non-rep structures; precompute each child's compatible parents (rep indices).
    child_parents = {}
    for s in structs:
        ci = idx[s]
        if ci in rep_set:
            continue
        parents = _candidate_parents_for_child(ci, structs, tok, postings, rep_indices, rep_set)
        if not parents:
            # WTA leaves these as their own reps too; collapse_gene already made them reps if childless.
            continue
        child_parents[ci] = parents

    # initialise rep abundance with own reads (ENST floor reads are 0, fine).
    abund = {pi: own[pi] for pi in rep_set}
    for _ in range(max(1, n_iter)):
        new = {pi: own[pi] for pi in rep_set}
        for ci, parents in child_parents.items():
            w = own[ci]
            if w <= 0:
                continue
            denom = sum(abund[pi] for pi in parents)
            if denom <= 0:
                # all parents currently zero -> split evenly (rare; happens iter 0 if parents have 0 own)
                share = w / len(parents)
                for pi in parents:
                    new[pi] += share
            else:
                for pi in parents:
                    new[pi] += w * (abund[pi] / denom)
        # convergence on total movement
        delta = sum(abs(new[pi] - abund[pi]) for pi in rep_set)
        abund = new
        if delta < tol:
            break

    res['em_reads'] = {structs[pi]: abund[pi] for pi in rep_set}

    # Final soft assignment weights: child_structure -> {parent_structure: fraction}, from the
    # converged parent abundances. Used to propagate a child's molecules FRACTIONALLY to parents in
    # the per-cell re-key (the per-cell analogue of em_reads). A rep maps to itself with weight 1.
    em_weights = {}
    for ci, parents in child_parents.items():
        denom = sum(abund[pi] for pi in parents)
        if denom <= 0:
            frac = {structs[pi]: 1.0 / len(parents) for pi in parents}
        else:
            frac = {structs[pi]: abund[pi] / denom for pi in parents}
        em_weights[structs[ci]] = frac
    res['em_weights'] = em_weights
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
