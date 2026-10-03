#!/usr/bin/env python3
"""Score a derived gating strategy against the gates a human actually drew.

The manual FlowJo tree is the reference. Two questions get separate answers, because they
fail differently:

  RECOVERY   for each manual population, how well does the best derived population
             reproduce it? A manual population with no good match is a miss.
  PRECISION  for each derived population, which manual population does it land on, and
             how cleanly? A derived population that spreads over many manual populations
             is not a gate, it is a smear.

Both are reported per population with its own denominator, never as a single headline.
F1 is used rather than ARI: ARI rewards a collapsed partition, which is the failure mode
this comparison exists to catch.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["pairwise_f1", "recovery_table", "precision_table", "summarize"]


def _masks_from_codes(codes, levels, min_n=1):
    codes = np.asarray(codes)
    return {lev: (codes == i) for i, lev in enumerate(levels)
            if (codes == i).sum() >= min_n}


def pairwise_f1(ref: dict, qry: dict):
    """F1 for every reference x query population pair. Returns a DataFrame (ref rows)."""
    rk, qk = list(ref), list(qry)
    R = np.array([ref[k] for k in rk], bool)
    Q = np.array([qry[k] for k in qk], bool)
    inter = (R.astype(np.int32) @ Q.T.astype(np.int32)).astype(float)
    rn = R.sum(1)[:, None].astype(float)
    qn = Q.sum(1)[None, :].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.where((rn + qn) > 0, 2.0 * inter / (rn + qn), 0.0)
    return pd.DataFrame(f1, index=rk, columns=qk)


def recovery_table(F: pd.DataFrame, ref_n: dict):
    rows = []
    for r in F.index:
        s = F.loc[r]
        j = s.values.argmax()
        rows.append(dict(manual_population=r, manual_n=int(ref_n[r]),
                         best_match=str(s.index[j]), f1=round(float(s.values[j]), 4),
                         runner_up=round(float(np.sort(s.values)[-2]) if len(s) > 1 else 0.0, 4)))
    return pd.DataFrame(rows).sort_values("f1", ascending=False).reset_index(drop=True)


def precision_table(F: pd.DataFrame, qry_n: dict):
    rows = []
    for c in F.columns:
        s = F[c]
        i = s.values.argmax()
        rows.append(dict(derived_population=c, derived_n=int(qry_n[c]),
                         lands_on=str(s.index[i]), f1=round(float(s.values[i]), 4)))
    return pd.DataFrame(rows).sort_values("f1", ascending=False).reset_index(drop=True)


def summarize(rec: pd.DataFrame, thresholds=(0.25, 0.5, 0.75)):
    out = {"manual_populations": int(len(rec)), "median_best_f1": round(float(rec.f1.median()), 4)}
    for t in thresholds:
        out["recovered_f1_ge_%.2f" % t] = int((rec.f1 >= t).sum())
    return out
