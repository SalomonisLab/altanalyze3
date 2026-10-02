"""Match flow antibody names to CITE-seq panel features.

Flow names a channel "<detector> :: <antibody>" (FlowJo) and a CITE-seq panel names a feature
by its own convention: "mouse_CD4", "AB_CD4", "CD4-TotalSeqC". Neither is canonical, so this
module normalizes both sides to a comparison key and then applies an explicit alias table for
the cases normalization cannot reach (IL-7R and CD127 are one protein; Int_B7 is integrin
beta-7; Sca-1 is Ly6a).

Nothing here guesses silently: `build_crosswalk` returns the matched pairs AND the unmatched
names on both sides, so a caller can print what failed rather than quietly map fewer markers.
"""
from __future__ import annotations

import re
import pandas as pd

__all__ = ["normalize_marker", "ALIASES", "build_crosswalk"]

# One protein, several written forms. Keys and values are normalized forms.
ALIASES = {
    "il7r": "cd127",          # IL-7R alpha chain is CD127; the panel ships both channels
    "intb7": "itgb7",         # integrin beta-7
    "integrinb7": "itgb7",
    "sca1": "ly6a",           # Sca-1 is Ly6a
    "nkp46": "cd335",         # NKp46 is CD335 / Ncr1
    "tcrb": "tcrbeta",
    "tcrbchain": "tcrbeta",
    "iaie": "mhcii",
    "ly6ae": "ly6a",
    "cd90.2": "cd90",
    "thy1.2": "cd90",
    "b220": "cd45r",
    "cd45rb220": "cd45r",
    # Documented decisions, not guesses:
    "tcrchain": "tcrbeta",   # the Chinese panel writes "mouse_TCR_chain"; the Greek beta was
                             # lost in the vendor name. Flow stains TCRb. Stated, not silent.
    "sca1": "ly6a",
}

_STRIP = re.compile(r"(^(mouse|human|rat|anti|ms|hu)[_\-\s]+)|([_\-\s]*totalseq[a-z]?$)", re.I)
_PREFIX = re.compile(r"^(ab|adt|prot)[_\-]", re.I)
_NONALNUM = re.compile(r"[^a-z0-9]")
# A CD number is the most reliable identity token a panel carries. "mouse_CD117_c_kit" and
# flow "CD117" are the same protein; "mouse_CD86" must NOT collapse onto "CD8". The optional
# trailing a/b is the chain (CD8a vs CD8b), kept so the two chains stay distinct.
_CDNUM = re.compile(r"cd(\d+)([ab])?(?![0-9])")


def normalize_marker(name: str) -> str:
    """A comparison key: lowercase, strip species/vendor decoration and punctuation."""
    s = str(name)
    if "::" in s:                      # FlowJo "detector :: antibody"
        s = s.split("::", 1)[1]
    s = s.strip()
    s = _PREFIX.sub("", s)
    prev = None
    while prev != s:                   # species prefixes can stack: "human_mouse_CD44"
        prev = s
        s = _STRIP.sub("", s).strip()
    for sp in ("mouse", "human", "rat"):
        s = re.sub(r"^%s[_\-\s]+" % sp, "", s, flags=re.I)
    s = _NONALNUM.sub("", s.lower())
    s = ALIASES.get(s, s)
    if s in ALIASES:
        return ALIASES[s]
    m = _CDNUM.search(s)
    if m:
        s = "cd" + m.group(1) + (m.group(2) or "")
    elif "sca1" in s or s.startswith("ly6a"):
        s = "ly6a"
    return ALIASES.get(s, s)


def build_crosswalk(flow_names, cite_names, extra_aliases=None):
    """Return (pairs_df, unmatched_flow, unmatched_cite).

    pairs_df columns: flow, cite, key. A flow marker matching several CITE features keeps all
    of them, because a panel legitimately carries two clones of one protein; the caller
    decides how to collapse.
    """
    al = dict(ALIASES)
    al.update({normalize_marker(k): normalize_marker(v) for k, v in (extra_aliases or {}).items()})

    def key(n):
        k = normalize_marker(n)
        return al.get(k, k)

    cite_by_key = {}
    for c in cite_names:
        cite_by_key.setdefault(key(c), []).append(c)

    rows, unmatched_flow = [], []
    for f in flow_names:
        k = key(f)
        if k in cite_by_key:
            for c in cite_by_key[k]:
                rows.append({"flow": f, "cite": c, "key": k})
        else:
            unmatched_flow.append(f)
    matched_keys = {r["key"] for r in rows}
    unmatched_cite = [c for c in cite_names if key(c) not in matched_keys]
    return pd.DataFrame(rows, columns=["flow", "cite", "key"]), unmatched_flow, unmatched_cite
