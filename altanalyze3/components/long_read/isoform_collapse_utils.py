"""Shared structure-token utilities for the cross-sample isoform collapse pipeline.

This module is the single, importable home for the AltAnalyze structure-string
primitives used by both the legacy ``gff_process`` collapse path and the new
two-tier collapse (Tier 1 ``collapse_sample`` / Tier 2 ``merge_catalog``).

Design reference: ``isoform-collapse/UNIFIED_isoform_collapse_design.md`` (Areas 3+4,
Milestone 1).

Two classes of helper live here:

1. Token parsing helpers that mirror the canonical AltAnalyze semantics already
   implemented in :mod:`gff_process` (``_split_token``, ``_strip_terminal_coords``,
   ``_filter_exon_intron_tokens``).  To guarantee there is exactly one definition
   of those semantics, this module *re-exports* them from ``gff_process`` rather
   than re-implementing them.  If that import is unavailable (e.g. the module is
   exercised in isolation without the heavy ``gff_process`` dependencies), an
   identical local fallback is used and is covered by the same unit tests.

2. The collapse primitives the unified design mandates:
   ``tokenize`` and ``is_contiguous_subsequence`` plus the gene-level
   ``contiguous_containment_fold`` driver.

CONTIGUOUS, not ordered.  A 5'/3'-truncated read shares a *contiguous* run of
interior tokens with its parent and may fold into it.  An exon-skip, intron
retention, or alternative donor/acceptor isoform is a *non-contiguous*
subsequence and is a biologically distinct isoform that must NOT fold.  See the
worked example in ``is_contiguous_subsequence``.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# 1. Canonical token helpers — re-exported from gff_process (single source of
#    truth) with an identical local fallback for isolated use.
# ---------------------------------------------------------------------------
try:  # pragma: no cover - exercised implicitly; both branches are equivalent
    from .gff_process import (
        _split_token as split_token,
        _strip_terminal_coords as strip_terminal_coords,
        _filter_exon_intron_tokens as filter_exon_intron_tokens,
    )
    _USING_GFF_PROCESS = True
except Exception:  # pragma: no cover - fallback path
    _USING_GFF_PROCESS = False

    def split_token(token, default_gene):
        """Parse ``ENSG..:E1.1_12345`` -> (gene_id, base='E1.1', coord=12345).

        Identical to ``gff_process._split_token``.  ``coord`` is ``None`` when the
        token carries no trailing ``_<int>`` novel coordinate.
        """
        gene_id = default_gene
        core = token
        if ':' in token:
            gene_id, core = token.split(':', 1)
        coord = None
        base = core
        if '_' in core:
            base, coord_str = core.rsplit('_', 1)
            try:
                coord = int(coord_str)
            except ValueError:
                base = core
                coord = None
        return gene_id, base, coord

    def strip_terminal_coords(tokens, default_gene):
        """Drop a leading/trailing exon token that carries a novel coordinate.

        Identical to ``gff_process._strip_terminal_coords``: only the first and
        last positions are considered, and only ``E`` tokens with a coordinate
        are dropped (uncertain 5'/3' ends).
        """
        if not tokens:
            return []
        trimmed = list(tokens)
        drop_indices = set()
        for idx in (0, len(trimmed) - 1):
            _, base, coord = split_token(trimmed[idx], default_gene)
            if coord is not None and base.startswith('E'):
                drop_indices.add(idx)
        if drop_indices:
            return [tok for i, tok in enumerate(trimmed) if i not in drop_indices]
        return trimmed

    def filter_exon_intron_tokens(tokens, default_gene):
        """Keep only annotated E/I tokens, dropping novel ``E*_coord`` sites.

        Identical to ``gff_process._filter_exon_intron_tokens``.
        """
        kept_tokens = []
        for tok in tokens:
            gene_id, base, _ = split_token(tok, default_gene)
            if base.startswith(('E', 'I')):
                if base.startswith('E') and '_' in tok:
                    continue
                if gene_id != default_gene and ':' not in tok:
                    tok = f"{gene_id}:{base}"
                kept_tokens.append(tok)
        return kept_tokens


# ---------------------------------------------------------------------------
# 2. Collapse primitives mandated by the unified design.
# ---------------------------------------------------------------------------
def tokenize(structure: str) -> List[str]:
    """Split an AltAnalyze structure string into its non-empty tokens.

    The structure strings are pipe-delimited, e.g. ``"E1.1|E2.2|E2.4|"``.  Trailing
    empties (from a terminal ``|``) are dropped.  Tokenising before any containment
    test is essential: it prevents the legacy substring bug where ``"E1|E12"`` would
    spuriously match inside the *string* ``"E1|E12|E23"`` at a non-token boundary.

    >>> tokenize("E1.1|E2.2|E2.4|")
    ['E1.1', 'E2.2', 'E2.4']
    >>> tokenize("")
    []
    """
    if not structure:
        return []
    return [t for t in structure.split('|') if t]


def exon_block_count(structure: str) -> int:
    """Number of DISTINCT exon blocks in a structure (length metric for binning).

    Sub-exon versions are collapsed: ``E1.1``, ``E1.2``, ``E1.3`` all count as the
    single exon block ``E1``.  Intron-retention tokens (``I*``) are NOT exon blocks
    and are excluded.  A trans-spliced gene prefix (``ENSG..:E1.1``) is dropped
    before taking the block.

    IMPORTANT — scope: this is used ONLY to compare the *number of exons* between
    candidate isoform models when binning "the longest" (see the read-weighted
    exon-block mode method).  It is NOT used as the collapse key and NOT used for
    contiguity tests — sub-exon structures remain fully distinct isoforms for
    collapsing.  (Confirmed design decision: "unique exons but ONLY for comparing
    number of exons in different models and NOT when comparing isoform structures
    for collapsing.")

    Empirically (real iPSC data), counting exon blocks instead of sub-exon tokens
    sharply unimodalises the read-weighted length distribution — e.g. MYL6/B2M/TPT1
    go from ~22-28% mode mass (sub-feature) to 78-88% (blocks) — because read mass
    that was split across sub-exon splice-boundary variants of the *same* block
    isoform now coheres.  See [[isoform_longest_bin_finding]].

    >>> exon_block_count("E1.1|E1.2|E2.1|I2.1|E3.1")
    3
    """
    seen = set()
    for tok in tokenize(structure):
        base = tok.split(':')[-1]          # drop trans-spliced gene prefix
        if not base.startswith('E'):       # exons only (skip I* introns)
            continue
        block = base.split('.')[0]         # E1.2 -> E1
        seen.add(block)
    return len(seen)


def is_contiguous_subsequence(short: Sequence[str], long: Sequence[str]) -> bool:
    """Return True iff ``short`` appears as a CONTIGUOUS run of tokens in ``long``.

    This is the single collapse test for terminal-truncation folding.  It is
    deliberately *contiguous* (a consecutive window), NOT an ordered
    subsequence, because contiguity is exactly what distinguishes a truncated
    read from a structurally distinct isoform:

        parent:                E1.1|E2.1|E3.1|E4.1|E5.1
        E2.1|E3.1|E4.1     ->  contiguous run            => FOLD (5'/3' truncation)
        E1.1|E3.1|E4.1|E5.1 -> E2.1 skipped, NOT a run   => DO NOT FOLD (exon skip)

    Contiguity is self-sufficient: an exon-skip / intron-retention / alt-donor /
    alt-acceptor variant can never be a contiguous window of the full structure,
    so no separate "junction-absent" guard is required.

    An empty ``short`` is treated as not-contained (there is no meaningful empty
    fragment to fold); this also keeps the function total for degenerate inputs.
    """
    n, m = len(short), len(long)
    if n == 0 or n > m:
        return False
    short = list(short)
    long = list(long)
    first = short[0]
    # Scan only windows whose first token matches, then compare the slice.
    for i in range(m - n + 1):
        if long[i] == first and long[i:i + n] == short:
            return True
    return False


def structure_tokens_for_containment(structure: str, default_gene: str = "") -> List[str]:
    """Return the token list used for containment tests for one structure.

    The collapse key has already had its uncertain terminal exon coordinates
    stripped upstream (``exon_str`` / :func:`strip_terminal_coords`).  For
    containment we additionally normalise away any *interior* novel ``E*_coord``
    sites and keep only annotated exon/intron tokens via
    :func:`filter_exon_intron_tokens`, so that two reads that differ only by a
    novel-site annotation on a shared junction still compare as the same run.
    Intron-retention tokens (``I*``) are RETAINED — they are real, distinguishing
    features and their presence is what keeps an IR isoform from folding into a
    spliced parent.
    """
    return filter_exon_intron_tokens(tokenize(structure), default_gene)


def contiguous_containment_fold(
    structures: Sequence[dict],
    default_gene: str = "",
) -> Tuple[List[dict], Dict[str, str]]:
    """Fold contiguous-subsequence fragments into their best parent within a gene.

    Parameters
    ----------
    structures:
        An iterable of dicts, one per distinct structure for a single gene.  Each
        must provide:
          - ``structure``  : the AltAnalyze structure string (collapse key)
          - ``count``      : aggregate read support (int)
          - ``is_ref``     : True if this structure is a reference (ENST/NM/..) model
        Any other keys (e.g. ``ref_id``, sample provenance) are passed through
        untouched on the winner dicts.
    default_gene:
        Gene id used to resolve bare/trans-spliced tokens.

    Returns
    -------
    (winners, fold_map)
        ``winners`` is the list of surviving (non-folded) structure dicts, each
        annotated with ``absorbed`` (list of folded child structure strings) and a
        ``count`` increased by the children it absorbed.  ``fold_map`` maps every
        folded child structure string -> the winner structure string it folded
        into.  Structures that did not fold appear only as winners (identity is
        implied; they are intentionally absent from ``fold_map``).

    Rules (from the unified design):
      * Parent priority is ``(is_ref desc, token-length desc, count desc)``;
        a fragment folds into the highest-priority parent that contains it as a
        contiguous token run.
      * A reference structure is NEVER folded as a fragment (references are
        canonical even when short).
      * Only strict containment folds: a structure never folds into itself, and
        equal-length structures never fold into each other (equal structures were
        already merged by exact-key aggregation upstream).
    """
    # Precompute tokens once.
    enriched = []
    for s in structures:
        toks = structure_tokens_for_containment(s["structure"], default_gene)
        enriched.append((s, toks, len(toks)))

    # Parent-priority order: references first, then longer, then better supported.
    order = sorted(
        range(len(enriched)),
        key=lambda i: (
            not enriched[i][0].get("is_ref", False),  # is_ref True sorts first
            -enriched[i][2],                           # longer first
            -int(enriched[i][0].get("count", 0)),      # higher count first
        ),
    )

    # Token-posting index: token -> set of structure indices containing it.
    # Used as a pre-filter so we only run the full contiguity test against
    # candidate parents that contain every token of the fragment.
    postings: Dict[str, set] = {}
    for i, (_s, toks, _n) in enumerate(enriched):
        for t in set(toks):
            postings.setdefault(t, set()).add(i)

    folded_into: Dict[int, int] = {}  # child idx -> parent idx (direct)

    # Iterate fragments shortest-first so a fragment prefers the most specific
    # (highest-priority) parent; we test candidates in parent-priority order.
    frag_order = sorted(range(len(enriched)), key=lambda i: enriched[i][2])
    for ci in frag_order:
        child, ctoks, cn = enriched[ci]
        if cn == 0:
            continue
        if child.get("is_ref", False):
            continue  # references are never folded as fragments
        if ci in folded_into:
            continue
        # Candidate parents: structures containing all of the child's tokens.
        candidate = None
        for t in ctoks:
            posting = postings.get(t)
            if not posting:
                candidate = set()
                break
            candidate = posting if candidate is None else (candidate & posting)
            if not candidate:
                break
        if not candidate:
            continue
        # Test candidate parents in parent-priority order.
        for pi in order:
            if pi == ci or pi not in candidate:
                continue
            parent, ptoks, pn = enriched[pi]
            if pn <= cn:
                continue  # strict containment only
            if pi in folded_into:
                continue  # a folded structure cannot be a parent
            if is_contiguous_subsequence(ctoks, ptoks):
                folded_into[ci] = pi
                break

    # Resolve to ultimate winners (no chained parents survive because parents that
    # were themselves folded are skipped above, but resolve defensively).
    def resolve(idx: int) -> int:
        seen = set()
        while idx in folded_into and idx not in seen:
            seen.add(idx)
            idx = folded_into[idx]
        return idx

    winners: List[dict] = []
    winner_by_idx: Dict[int, dict] = {}
    fold_map: Dict[str, str] = {}

    for i, (s, _toks, _n) in enumerate(enriched):
        if i not in folded_into:
            w = dict(s)
            w.setdefault("absorbed", [])
            winners.append(w)
            winner_by_idx[i] = w

    for ci in folded_into:
        pi = resolve(ci)
        child_struct = enriched[ci][0]["structure"]
        parent_struct = enriched[pi][0]["structure"]
        fold_map[child_struct] = parent_struct
        w = winner_by_idx.get(pi)
        if w is not None:
            w["absorbed"].append(child_struct)
            w["count"] = int(w.get("count", 0)) + int(enriched[ci][0].get("count", 0))

    return winners, fold_map


__all__ = [
    "split_token",
    "strip_terminal_coords",
    "filter_exon_intron_tokens",
    "tokenize",
    "exon_block_count",
    "is_contiguous_subsequence",
    "structure_tokens_for_containment",
    "contiguous_containment_fold",
]
