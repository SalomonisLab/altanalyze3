"""Canonical isoform identity: the Ensembl91 structure string and how to match two of them.

The structure string is the pipe-delimited exon/intron token list, e.g.
``E1.2|E1.3|I1.1|E2.1|...|E11.1`` (Ensembl91 reference). It is the cross-project isoform key (see
GOALS.md sec. 2). This module is PURE (no I/O) so it is trivially testable and reusable.

Match semantics mirror the prior TFIso annotation (``tfiso_to_final_assignment``):
  - exact         : identical token sequence
  - token_exact   : identical token *set/order* after normalization (e.g. whitespace/case), used when
                    raw strings differ only cosmetically
  - contained_in  : A is a contiguous token subsequence of B (A folds into the longer isoform B)
  - superset_of   : B is a contiguous token subsequence of A (only a fragment of A is present in B)
  - none          : no containment relationship

Containment is *contiguous* (a run of consecutive tokens), matching splice-structure semantics: a
shorter isoform observed as a 5'/3' truncation of a longer one is a contiguous block of its tokens.
"""

SEP = "|"


def normalize_structure(structure):
    """Return the canonical token tuple for a structure string. Strips whitespace around the string
    and around each token; drops empty tokens; preserves order (order is meaningful). Returns () for
    None/empty input."""
    if structure is None:
        return ()
    s = str(structure).strip()
    if not s:
        return ()
    return tuple(tok.strip() for tok in s.split(SEP) if tok.strip())


def _is_contiguous_subsequence(short, long_):
    """True iff tuple ``short`` appears as a run of consecutive elements within tuple ``long_``.
    Empty ``short`` is never a match (no information)."""
    n, m = len(short), len(long_)
    if n == 0 or n > m:
        return False
    first = short[0]
    for i in range(m - n + 1):
        if long_[i] == first and long_[i:i + n] == short:
            return True
    return False


def match_structures(a, b):
    """Classify the relationship of structure ``a`` to structure ``b`` (both raw strings). Returns one
    of: 'exact', 'token_exact', 'contained_in', 'superset_of', 'none'. Directionality: 'contained_in'
    means a is contained in b; 'superset_of' means a contains b."""
    ta, tb = normalize_structure(a), normalize_structure(b)
    if not ta or not tb:
        return "none"
    if ta == tb:
        return "exact"
    # token_exact: same multiset of tokens in the same order after normalization but raw strings
    # differed (already covered by exact above on normalized tuples); reserved for callers that
    # normalize differently (e.g. dropping terminal-coordinate suffixes). Kept for parity with the
    # prior pipeline's vocabulary.
    if _is_contiguous_subsequence(ta, tb):
        return "contained_in"
    if _is_contiguous_subsequence(tb, ta):
        return "superset_of"
    return "none"


# match types that constitute a confident isoform-level assignment (vs 'superset_of'/'none'/'novel')
ASSIGNED_MATCH_TYPES = frozenset({"exact", "token_exact", "contained_in"})


def is_assigned(match_type):
    """True if ``match_type`` represents a confident link to an observed isoform."""
    return match_type in ASSIGNED_MATCH_TYPES
