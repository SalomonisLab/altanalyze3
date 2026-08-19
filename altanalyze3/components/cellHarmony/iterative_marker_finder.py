#!/usr/bin/env python3
"""
iterative_marker_finder.py

Consolidated MarkerFinder across several grouping schemes.

The AltAnalyze2 module `iterativeMarkerFinder.py` ran MarkerFinder once per grouping file, then
kept the single strongest (rho, state, direction) per gene. That collapse hides every program that
two conditions share. This module keeps the full feature-by-scheme correlation matrix, codes each
feature into a sign pattern, and separates unique, shared and global programs the way
`cellHarmony_differential.py` separates local, co-regulated and global genes.

The statistics come from `altanalyze3.components.cellHarmony.markerFinder.marker_finder`.
The figure comes from `altanalyze3.components.visualization.marker_heatmap_h5ad._plot_heatmap`,
the established MarkerFinder heatmap renderer. This module adds no statistic and no plotting code
of its own.

Stages
  1. score_schemes      - one marker_finder call per scheme, full signed rho and p retained
  2. code_signs         - -1 / 0 / +1 per (feature, scheme) under a rho gate and an FDR gate
  3. assign_blocks      - unique / shared / global / discordant, one block per feature
  4. pooled_retest      - each shared and global block re-tested at the sample level
  5. consolidated_heatmap - the established renderer, rows in block order
  6. goelite_per_block  - the altanalyze3 GO-Elite CLI, one query per block

All code edits must be validated within altanalyze3.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

if __package__ in (None, ""):
    _repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    if _repo_root not in sys.path:
        sys.path.insert(0, _repo_root)

from altanalyze3.components.cellHarmony.markerFinder import marker_finder
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh

CONTRAST_SEP = "::"
BLOCK_KIND_ORDER = ["unique", "shared", "global", "single-coverage", "discordant"]


# --------------------------------------------------------------------------------------- inputs


@dataclass
class Scheme:
    """One grouping scheme: the AltAnalyze2 equivalent of one `exp.`/`groups.` file pair.

    name        : scheme label, used in every output and in the figure
    adata       : samples x features. `adata.X` must hold no NaN, because
                  `marker_finder` drops any feature that returns a NaN correlation.
    group_key   : `adata.obs` column holding the levels of this scheme
    level_order : the levels in report order. With exactly 2 levels the correlation of the
                  first is the negative of the second, so the module keeps one signed column
                  and a positive value means higher in `level_order[-1]`.
    testable    : optional per-feature bool Series. False marks a feature this scheme could
                  not assess, which stays distinct from "assessed and not significant".
    sample_key  : optional `adata.obs` column naming the biological sample. The pooled
                  re-test needs it to keep the sample, not the scheme, as the unit.
    """

    name: str
    adata: AnnData
    group_key: str
    level_order: Sequence[str]
    testable: Optional[pd.Series] = None
    sample_key: Optional[str] = None
    measured: Optional[pd.DataFrame] = None

    def levels(self) -> List[str]:
        return [str(level) for level in self.level_order]

    def groups(self) -> List[str]:
        return self.adata.obs[self.group_key].astype(str).tolist()

    def features(self) -> pd.Index:
        return pd.Index(self.adata.var_names.astype(str))


@dataclass
class IterativeResults:
    rho: pd.DataFrame
    pval: pd.DataFrame
    fdr: pd.DataFrame
    testable: pd.DataFrame
    signs: pd.DataFrame
    assignments: pd.DataFrame
    blocks: pd.DataFrame
    diagnostics: List[dict] = field(default_factory=list)


# ------------------------------------------------------------------------------- stage 1: score


def _bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg over the finite entries. NaN stays NaN and never enters the rank."""
    values = pd.to_numeric(pvals, errors="coerce")
    finite = values.dropna()
    if finite.empty:
        return pd.Series(np.nan, index=values.index)
    adjusted = pd.Series(mh._bh_fdr(finite.to_numpy(dtype=float)), index=finite.index)
    return adjusted.reindex(values.index)


def score_schemes(
    schemes: Sequence[Scheme],
    features: Optional[Sequence[str]] = None,
    complementary_tol: float = 1e-8,
    min_range: float = 1e-9,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, List[dict]]:
    """One `marker_finder` call per scheme. Returns rho, p, FDR, the testable mask and diagnostics.

    Columns are the scheme name for a 2-level scheme, or `scheme::level` otherwise. The rows are
    the union of the schemes' features unless `features` fixes the set. A feature that a scheme
    never scored comes back as NaN, never as 0.
    """
    if not schemes:
        raise ValueError("score_schemes needs at least one scheme.")

    index = pd.Index([str(f) for f in features]) if features is not None else None
    if index is None:
        union: List[str] = []
        seen = set()
        for scheme in schemes:
            for feature in scheme.features():
                if feature not in seen:
                    seen.add(feature)
                    union.append(feature)
        index = pd.Index(union)
    if index.has_duplicates:
        raise ValueError("The feature index holds duplicates; MarkerFinder rows must be unique.")

    rho_cols: Dict[str, pd.Series] = {}
    p_cols: Dict[str, pd.Series] = {}
    testable_cols: Dict[str, pd.Series] = {}
    diagnostics: List[dict] = []

    for scheme in schemes:
        levels = scheme.levels()
        observed = sorted(set(scheme.groups()))
        missing = [level for level in levels if level not in observed]
        if missing:
            raise ValueError(f"{scheme.name}: level(s) {missing} absent from obs['{scheme.group_key}'].")
        if len(levels) < 2:
            raise ValueError(f"{scheme.name}: a scheme needs at least 2 levels.")

        matrix = scheme.adata.X
        # float64 for the correlation. `marker_finder` forms the covariance as
        # `sum_xy - sum_x*sum_y/n` (markerFinder.py:160), a subtraction of 2 nearly equal numbers,
        # then divides by the standard deviation. In float32 that cancellation leaves an error
        # near 1e-5 for a feature with a small range. The cast changes no value, only the
        # precision the sum is carried at.
        dense = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        dense = np.asarray(dense, dtype=np.float64)
        matrix = dense
        if not np.isfinite(dense).all():
            raise ValueError(
                f"{scheme.name}: the matrix holds NaN or inf. marker_finder drops those features "
                "silently, so fill or remove them before scoring and report what you did."
            )

        # A feature that does not vary inside this scheme has an undefined correlation: the engine
        # divides by a zero standard deviation. It returns a meaningless ratio of 2 rounding
        # errors, or +/-inf, and its `dropna` (markerFinder.py:180) keeps both. Such a value can
        # exceed any rho gate and rank first. Identify these features from the data BEFORE the
        # test, mark them unassessable in this scheme, and count them.
        feature_range = dense.max(axis=0) - dense.min(axis=0)
        constant = pd.Series(feature_range <= float(min_range), index=scheme.features())
        n_constant = int(constant.sum())

        r_df, p_df = marker_finder(matrix, scheme.groups(), gene_names=scheme.features())
        n_dropped = int(len(scheme.features()) - r_df.shape[0])

        degenerate = constant.reindex(r_df.index).fillna(False)
        if degenerate.any():
            r_df.loc[degenerate.to_numpy(), :] = np.nan
            p_df.loc[degenerate.to_numpy(), :] = np.nan
        nonfinite_mask = ~np.isfinite(r_df.to_numpy(dtype=float))
        n_nonfinite = int(nonfinite_mask.any(axis=1).sum())
        if n_nonfinite:
            r_df = r_df.where(np.isfinite(r_df))
            p_df = p_df.where(r_df.notna())

        if len(levels) == 2:
            pair = (r_df[levels[0]].to_numpy(dtype=float) + r_df[levels[1]].to_numpy(dtype=float))
            pair = pair[np.isfinite(pair)]
            complementary = float(np.max(np.abs(pair))) if pair.size else 0.0
            if complementary > complementary_tol:
                raise ValueError(
                    f"{scheme.name}: the 2 levels are not complementary "
                    f"(max |r0 + r1| = {complementary:.3g} > {complementary_tol:g})."
                )
            named = {scheme.name: levels[-1]}
            diagnostics.append(
                {
                    "scheme": scheme.name,
                    "check": "complementary_levels",
                    "max_abs_r_sum": complementary,
                    "positive_means": levels[-1],
                }
            )
        else:
            named = {f"{scheme.name}{CONTRAST_SEP}{level}": level for level in levels}

        for column, level in named.items():
            rho_cols[column] = r_df[level].reindex(index)
            p_cols[column] = p_df[level].reindex(index)
            if scheme.testable is None:
                mask = pd.Series(True, index=index)
            else:
                mask = scheme.testable.reindex(index).fillna(False).astype(bool)
            mask = mask & rho_cols[column].notna()
            testable_cols[column] = mask

        diagnostics.append(
            {
                "scheme": scheme.name,
                "check": "marker_finder_rows",
                "n_features_in": int(len(scheme.features())),
                "n_features_scored": int(r_df.shape[0]),
                "n_features_dropped_nan_rho": n_dropped,
                "n_features_constant_in_scheme": n_constant,
                "n_features_nonfinite_rho": n_nonfinite,
                "n_samples": int(scheme.adata.n_obs),
                "levels": ";".join(levels),
            }
        )

    rho = pd.DataFrame(rho_cols, index=index)
    pval = pd.DataFrame(p_cols, index=index)
    testable = pd.DataFrame(testable_cols, index=index)
    fdr = pd.DataFrame({column: _bh_fdr(pval[column].where(testable[column])) for column in pval.columns},
                       index=index)
    return rho, pval, fdr, testable, diagnostics


# ---------------------------------------------------------------- stage 1b: the broad contrast


def donor_level_matrix(
    schemes: Sequence[Scheme],
    features: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Average each sample's MEASURED values over the schemes it appears in.

    The broad contrast asks one question of the whole data set, so its unit must be the sample.
    Pooling the scheme columns instead would count a donor measured in 5 schemes 5 times and
    inflate every p-value. Returns the sample-by-feature matrix, the level of each sample, and the
    number of schemes that measured each cell.
    """
    if features is None:
        union: List[str] = []
        seen = set()
        for scheme in schemes:
            for feature in scheme.features():
                if feature not in seen:
                    seen.add(feature)
                    union.append(feature)
        features = union
    features = pd.Index([str(f) for f in features])

    samples: List[str] = []
    levels: Dict[str, str] = {}
    for scheme in schemes:
        if not scheme.sample_key:
            raise ValueError(f"{scheme.name}: the broad contrast needs Scheme.sample_key.")
        for sample, level in zip(
            scheme.adata.obs[scheme.sample_key].astype(str),
            scheme.adata.obs[scheme.group_key].astype(str),
        ):
            if sample in levels and levels[sample] != level:
                raise ValueError(f"Sample {sample} carries 2 levels: {levels[sample]} and {level}.")
            if sample not in levels:
                samples.append(sample)
            levels[sample] = level

    total = pd.DataFrame(0.0, index=pd.Index(samples), columns=features)
    count = pd.DataFrame(0.0, index=pd.Index(samples), columns=features)
    for scheme in schemes:
        present = [f for f in features if f in set(scheme.features())]
        if not present:
            continue
        sub = scheme.adata[:, present]
        values = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
        frame = pd.DataFrame(np.asarray(values, dtype=np.float64),
                             index=sub.obs[scheme.sample_key].astype(str), columns=present)
        if scheme.measured is not None:
            mask = scheme.measured.reindex(index=frame.index, columns=present).fillna(False)
            mask = mask.astype(bool)
        else:
            mask = pd.DataFrame(True, index=frame.index, columns=present)
        frame = frame.where(mask, 0.0)
        total.loc[frame.index, present] = total.loc[frame.index, present].to_numpy() + frame.to_numpy()
        count.loc[frame.index, present] = count.loc[frame.index, present].to_numpy() + mask.to_numpy().astype(float)

    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = pd.DataFrame(total.to_numpy() / count.to_numpy(), index=total.index, columns=features)
    pooled = pooled.where(count.to_numpy() > 0)
    return pooled, pd.Series({sample: levels[sample] for sample in samples}), count


def broad_contrast(
    schemes: Sequence[Scheme],
    level_order: Sequence[str],
    *,
    features: Optional[Sequence[str]] = None,
    min_samples_per_level: Optional[Dict[str, int]] = None,
    min_range: float = 1e-9,
) -> dict:
    """MarkerFinder on the broad contrast, over every scheme at once, one row per sample.

    Returns the correlation, the p-value, the BH value, the mean difference between the levels,
    and the mask of features enough samples measured. A positive correlation means higher in
    `level_order[-1]`.
    """
    levels = [str(level) for level in level_order]
    if len(levels) != 2:
        raise ValueError("The broad contrast needs exactly 2 levels.")
    pooled, sample_levels, count = donor_level_matrix(schemes, features=features)
    features_index = pooled.columns

    measured = pooled.notna()
    per_level = {
        level: measured.loc[sample_levels[sample_levels == level].index].sum(axis=0)
        for level in levels
    }
    floors = dict(min_samples_per_level or {})
    for level in levels:
        floors.setdefault(level, int((sample_levels == level).sum()))
    testable = pd.Series(True, index=features_index)
    for level in levels:
        testable &= per_level[level] >= int(floors[level])

    event_mean = pooled.mean(axis=0, skipna=True)
    filled = pooled.fillna(event_mean)
    n_filled = int(pooled.isna().to_numpy().sum())
    n_filled_testable = int(pooled.loc[:, testable[testable].index].isna().to_numpy().sum())

    feature_range = filled.max(axis=0) - filled.min(axis=0)
    constant = feature_range <= float(min_range)
    testable &= ~constant

    groups = [sample_levels[sample] for sample in filled.index]
    r_df, p_df = marker_finder(filled.to_numpy(dtype=np.float64), groups, gene_names=features_index)
    r_df = r_df.reindex(features_index)
    p_df = p_df.reindex(features_index)
    reference = levels[-1]
    rho = r_df[reference].where(np.isfinite(r_df[reference]) & testable)
    pval = p_df[reference].where(rho.notna())
    fdr = _bh_fdr(pval)

    positive = filled.loc[[s for s in filled.index if sample_levels[s] == levels[-1]]].mean(axis=0)
    negative = filled.loc[[s for s in filled.index if sample_levels[s] == levels[0]]].mean(axis=0)
    effect = (positive - negative).where(testable)

    return {
        "rho": rho,
        "pval": pval,
        "fdr": fdr,
        "effect": effect,
        "testable": testable,
        "pooled": filled,
        "sample_levels": sample_levels,
        "n_samples": int(filled.shape[0]),
        "n_samples_per_level": {level: int((sample_levels == level).sum()) for level in levels},
        "n_testable": int(testable.sum()),
        "n_constant": int(constant.sum()),
        "n_filled_cells": n_filled,
        "n_filled_cells_in_tested_features": n_filled_testable,
        "positive_means": reference,
    }


def residual_schemes(
    schemes: Sequence[Scheme],
    effect: pd.Series,
    positive_level: str,
) -> List[Scheme]:
    """Remove the broad effect from every scheme, so what remains is scheme-specific.

    Each sample of `positive_level` loses the broad mean difference of its feature. A feature whose
    scheme-specific difference equals the broad one leaves no residual difference. A feature whose
    difference is larger, smaller or opposite keeps one, and MarkerFinder can then find it.
    """
    out: List[Scheme] = []
    for scheme in schemes:
        adata = scheme.adata.copy()
        values = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        values = np.asarray(values, dtype=np.float64)
        shift = effect.reindex(scheme.features()).to_numpy(dtype=float)
        shift = np.nan_to_num(shift, nan=0.0)
        indicator = (adata.obs[scheme.group_key].astype(str) == str(positive_level)).to_numpy(dtype=float)
        adata.X = (values - np.outer(indicator, shift)).astype(np.float64)
        out.append(
            Scheme(
                name=scheme.name,
                adata=adata,
                group_key=scheme.group_key,
                level_order=scheme.level_order,
                testable=scheme.testable,
                sample_key=scheme.sample_key,
                measured=scheme.measured,
            )
        )
    return out


# -------------------------------------------------------------------------------- stage 2: sign


def code_signs(
    rho: pd.DataFrame,
    fdr: pd.DataFrame,
    testable: pd.DataFrame,
    rho_min: float,
    alpha: float,
    pval: Optional[pd.DataFrame] = None,
    use_rawp: bool = False,
) -> pd.DataFrame:
    """-1, 0 or +1 per (feature, scheme). NaN marks a scheme that could not assess the feature.

    A feature is significant when it clears BOTH gates: |rho| >= rho_min and the significance
    metric < alpha. `use_rawp` swaps the BH value for the raw p-value and needs `pval`.
    """
    if use_rawp:
        if pval is None:
            raise ValueError("use_rawp=True needs the raw p-value frame.")
        metric = pval
    else:
        metric = fdr

    significant = (rho.abs() >= float(rho_min)) & (metric < float(alpha)) & testable
    signs = pd.DataFrame(0.0, index=rho.index, columns=rho.columns)
    signs = signs.where(~significant, np.sign(rho))
    return signs.where(testable, np.nan)


# ------------------------------------------------------------------------------ stage 3: blocks


def _direction_name(sign: float, level_names: Tuple[str, str]) -> str:
    return level_names[1] if sign > 0 else level_names[0]


def assign_blocks(
    signs: pd.DataFrame,
    level_names: Tuple[str, str] = ("down", "up"),
    min_features_per_pattern: int = 10,
    scheme_order: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Assign every significant feature to exactly one block.

    unique     - significant in 1 of its testable schemes
    shared     - significant in >= 2 testable schemes, one sign, but not in all of them
    global     - significant in every scheme that could test it, and in >= 2
    discordant - significant in >= 2 schemes with opposite signs; never merged into a program

    A shared pattern carried by fewer than `min_features_per_pattern` features moves to a
    pooled `shared-other` block of the same size class. Nothing is discarded.
    """
    columns = list(scheme_order) if scheme_order else list(signs.columns)
    missing = [column for column in columns if column not in signs.columns]
    if missing:
        raise ValueError(f"scheme_order names columns absent from the sign matrix: {missing}")
    frame = signs[columns]

    rows = []
    for feature, series in frame.iterrows():
        testable_cols = [column for column in columns if not pd.isna(series[column])]
        positive = [column for column in testable_cols if series[column] > 0]
        negative = [column for column in testable_cols if series[column] < 0]
        n_sig = len(positive) + len(negative)
        if n_sig == 0:
            continue

        if positive and negative:
            kind = "discordant"
            members = [column for column in testable_cols if series[column] != 0]
            direction = "mixed"
        else:
            members = positive if positive else negative
            direction = _direction_name(1.0 if positive else -1.0, level_names)
            if n_sig == 1 and len(testable_cols) == 1:
                # No other scheme could test this feature, so nothing here is specific to the one
                # scheme that could. Calling it unique would read a coverage difference as biology.
                kind = "single-coverage"
            elif n_sig == 1:
                kind = "unique"
            elif n_sig == len(columns) and len(testable_cols) == len(columns):
                kind = "global"
            else:
                kind = "shared"

        rows.append(
            {
                "feature": feature,
                "kind": kind,
                "direction": direction,
                "members": "+".join(members),
                "n_schemes_significant": n_sig,
                "n_schemes_testable": len(testable_cols),
            }
        )

    assignments = pd.DataFrame(
        rows,
        columns=["feature", "kind", "direction", "members", "n_schemes_significant", "n_schemes_testable"],
    )
    if assignments.empty:
        return assignments.assign(block=pd.Series(dtype=str)), pd.DataFrame(
            columns=["block", "kind", "direction", "members", "n_members", "n_features"]
        )

    def _label(row):
        if row["kind"] == "discordant":
            return "discordant|mixed"
        if row["kind"] == "single-coverage":
            return f"{row['members']}-1scheme|{row['direction']}"
        if row["kind"] == "unique":
            return f"{row['members']}-only|{row['direction']}"
        if row["kind"] == "global":
            return f"all-{len(columns)}|{row['direction']}"
        return f"{row['members']}|{row['direction']}"

    assignments["block"] = assignments.apply(_label, axis=1)

    counts = assignments["block"].value_counts()
    small = {
        block
        for block in counts.index
        if counts[block] < int(min_features_per_pattern)
        and assignments.loc[assignments["block"] == block, "kind"].iloc[0] == "shared"
    }
    if small:
        mask = assignments["block"].isin(small)
        assignments.loc[mask, "block"] = (
            "shared-other-"
            + assignments.loc[mask, "n_schemes_significant"].astype(str)
            + "|"
            + assignments.loc[mask, "direction"]
        )

    blocks = (
        assignments.groupby("block")
        .agg(
            kind=("kind", "first"),
            direction=("direction", "first"),
            n_members=("n_schemes_significant", "max"),
            n_features=("feature", "size"),
        )
        .reset_index()
    )
    blocks["members"] = [
        ";".join(sorted(set(assignments.loc[assignments["block"] == block, "members"])))
        for block in blocks["block"]
    ]

    order = order_blocks(blocks, columns, level_names)
    blocks["display_order"] = [order.index(block) for block in blocks["block"]]
    blocks = blocks.sort_values("display_order").reset_index(drop=True)
    return assignments, blocks


def order_blocks(
    blocks: pd.DataFrame,
    scheme_order: Sequence[str],
    level_names: Tuple[str, str],
) -> List[str]:
    """Row order for the figure: unique, then shared, then global, then discordant.

    The order follows `build_fixed_order_heatmap` in cellHarmony_differential.py, which lists the
    local genes first, then the co-regulated patterns, then the global ones.
    """
    directions = [level_names[0], level_names[1]]

    def _direction_rank(value):
        return directions.index(value) if value in directions else len(directions)

    ordered: List[str] = []
    frame = blocks.set_index("block")

    for direction in directions:
        for scheme in scheme_order:
            label = f"{scheme}-only|{direction}"
            if label in frame.index:
                ordered.append(label)

    shared = frame[frame["kind"] == "shared"].copy()
    shared = shared[~shared.index.isin(ordered)]
    if not shared.empty:
        shared["direction_rank"] = shared["direction"].map(_direction_rank)
        shared["is_other"] = shared.index.str.startswith("shared-other-")
        shared = shared.sort_values(
            ["is_other", "n_members", "direction_rank", "n_features"],
            ascending=[True, False, True, False],
        )
        ordered.extend(shared.index.tolist())

    for direction in directions:
        for label in frame.index[(frame["kind"] == "global") & (frame["direction"] == direction)]:
            ordered.append(label)

    for direction in directions:
        for scheme in scheme_order:
            label = f"{scheme}-1scheme|{direction}"
            if label in frame.index:
                ordered.append(label)

    for label in frame.index[frame["kind"] == "discordant"]:
        ordered.append(label)

    for label in frame.index:
        if label not in ordered:
            ordered.append(label)
    return ordered


def assign_blocks_broad_first(
    broad_signs: pd.Series,
    residual_signs: pd.DataFrame,
    level_names: Tuple[str, str] = ("down", "up"),
    min_features_per_pattern: int = 10,
    scheme_order: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Broad first, then whatever survives it.

    The broad contrast decides first. A feature it calls significant belongs to the broad program,
    whatever the individual schemes say, because the broad difference already explains it. Only a
    feature the broad contrast does NOT call keeps a scheme-specific block, and only when its
    residual, after the broad effect is removed, is still significant in that scheme.

    A broad feature that also carries a residual scheme effect keeps that fact in the column
    `also_specific_in`, so an amplified program stays visible without leaving the broad block.
    """
    columns = list(scheme_order) if scheme_order else list(residual_signs.columns)
    missing = [column for column in columns if column not in residual_signs.columns]
    if missing:
        raise ValueError(f"scheme_order names columns absent from the residual sign matrix: {missing}")
    frame = residual_signs[columns]
    broad_signs = broad_signs.reindex(frame.index)

    # Only a feature the broad contrast or some scheme calls significant can enter a block, so the
    # loop below skips the rest instead of walking every feature.
    carries_signal = (broad_signs.fillna(0).abs() > 0) | (frame.abs().sum(axis=1) > 0)
    frame = frame.loc[carries_signal.to_numpy()]

    rows = []
    for feature, series in frame.iterrows():
        testable_cols = [column for column in columns if not pd.isna(series[column])]
        positive = [column for column in testable_cols if series[column] > 0]
        negative = [column for column in testable_cols if series[column] < 0]
        residual_members = positive + negative
        broad = broad_signs.get(feature, np.nan)

        if pd.notna(broad) and broad != 0:
            rows.append(
                {
                    "feature": feature,
                    "kind": "broad",
                    "direction": _direction_name(float(broad), level_names),
                    "members": "all-schemes",
                    "n_schemes_significant": len(residual_members),
                    "n_schemes_testable": len(testable_cols),
                    "also_specific_in": "+".join(residual_members),
                }
            )
            continue

        if not residual_members:
            continue

        if positive and negative:
            kind = "discordant"
            direction = "mixed"
            members = residual_members
        else:
            members = positive if positive else negative
            direction = _direction_name(1.0 if positive else -1.0, level_names)
            kind = "scheme-specific" if len(members) == 1 else "shared-specific"

        rows.append(
            {
                "feature": feature,
                "kind": kind,
                "direction": direction,
                "members": "+".join(members),
                "n_schemes_significant": len(members),
                "n_schemes_testable": len(testable_cols),
                "also_specific_in": "",
            }
        )

    assignments = pd.DataFrame(
        rows,
        columns=["feature", "kind", "direction", "members", "n_schemes_significant",
                 "n_schemes_testable", "also_specific_in"],
    )
    if assignments.empty:
        return assignments.assign(block=pd.Series(dtype=str)), pd.DataFrame(
            columns=["block", "kind", "direction", "members", "n_members", "n_features"]
        )

    def _label(row):
        if row["kind"] == "broad":
            return f"broad|{row['direction']}"
        if row["kind"] == "discordant":
            return "discordant|mixed"
        return f"{row['members']}-specific|{row['direction']}"

    assignments["block"] = assignments.apply(_label, axis=1)

    counts = assignments["block"].value_counts()
    small = {
        block
        for block in counts.index
        if counts[block] < int(min_features_per_pattern)
        and assignments.loc[assignments["block"] == block, "kind"].iloc[0] == "shared-specific"
    }
    if small:
        mask = assignments["block"].isin(small)
        assignments.loc[mask, "block"] = (
            "shared-specific-other-"
            + assignments.loc[mask, "n_schemes_significant"].astype(str)
            + "|"
            + assignments.loc[mask, "direction"]
        )

    blocks = (
        assignments.groupby("block")
        .agg(
            kind=("kind", "first"),
            direction=("direction", "first"),
            n_members=("n_schemes_significant", "max"),
            n_features=("feature", "size"),
        )
        .reset_index()
    )
    blocks["members"] = [
        ";".join(sorted(set(assignments.loc[assignments["block"] == block, "members"])))
        for block in blocks["block"]
    ]

    order = order_blocks_broad_first(blocks, columns, level_names)
    blocks["display_order"] = [order.index(block) for block in blocks["block"]]
    blocks = blocks.sort_values("display_order").reset_index(drop=True)
    return assignments, blocks


def order_blocks_broad_first(
    blocks: pd.DataFrame,
    scheme_order: Sequence[str],
    level_names: Tuple[str, str],
) -> List[str]:
    """Broad blocks first, then the shared-specific ones, then one scheme at a time."""
    directions = [level_names[0], level_names[1]]
    ordered: List[str] = []
    frame = blocks.set_index("block")

    for direction in directions:
        label = f"broad|{direction}"
        if label in frame.index:
            ordered.append(label)

    shared = frame[(frame["kind"] == "shared-specific") & (~frame.index.isin(ordered))].copy()
    if not shared.empty:
        shared["direction_rank"] = shared["direction"].map(
            lambda value: directions.index(value) if value in directions else len(directions))
        shared["is_other"] = shared.index.str.startswith("shared-specific-other-")
        shared = shared.sort_values(
            ["is_other", "n_members", "direction_rank", "n_features"],
            ascending=[True, False, True, False],
        )
        ordered.extend(shared.index.tolist())

    for scheme in scheme_order:
        for direction in directions:
            label = f"{scheme}-specific|{direction}"
            if label in frame.index and label not in ordered:
                ordered.append(label)

    for label in frame.index[frame["kind"] == "discordant"]:
        if label not in ordered:
            ordered.append(label)

    for label in frame.index:
        if label not in ordered:
            ordered.append(label)
    return ordered


# ------------------------------------------------------------------------ stage 4: pooled re-test


def pooled_retest(
    schemes: Sequence[Scheme],
    assignments: pd.DataFrame,
    blocks: pd.DataFrame,
    alpha: float = 0.05,
    min_samples_per_level: int = 2,
) -> pd.DataFrame:
    """Re-test every multi-scheme block on its own pooled contrast.

    A shared block claims that several schemes move together. The claim needs its own statistic,
    not only the intersection of separate tests. This pools the member schemes at the SAMPLE
    level: each sample's value is the mean of that feature over the member schemes that measured
    it, so a sample that appears in several schemes still counts once. Pooling the scheme columns
    instead would count the same donor several times and inflate the significance.

    Needs `Scheme.sample_key`. Returns one row per block with the confirmed fraction.
    """
    by_name = {scheme.name: scheme for scheme in schemes}
    rows = []

    for _, block in blocks.iterrows():
        label = block["block"]
        kind = block["kind"]
        features = assignments.loc[assignments["block"] == label, "feature"].tolist()
        member_sets = [
            set(value.split("+"))
            for value in assignments.loc[assignments["block"] == label, "members"].unique()
        ]
        members = sorted(set().union(*member_sets)) if member_sets else []
        members = [name for name in members if name in by_name]

        # A pooled contrast only means something when every feature in the block carries the SAME
        # member set. A `shared-other` block merges several rare patterns, so pooling the union of
        # their members would test a contrast no feature actually showed.
        heterogeneous = len(member_sets) > 1
        if (
            kind in ("unique", "single-coverage", "discordant")
            or len(members) < 2
            or heterogeneous
            or not features
        ):
            if kind in ("unique", "single-coverage", "discordant"):
                note = "not applicable to a single-scheme or discordant block"
            elif heterogeneous:
                note = f"skipped: {len(member_sets)} different member sets in this block"
            else:
                note = "fewer than 2 member schemes"
            rows.append(
                {
                    "block": label,
                    "kind": kind,
                    "n_features": len(features),
                    "n_member_schemes": len(members),
                    "pooled_n_samples": np.nan,
                    "pooled_n_features_tested": 0,
                    "pooled_n_features_confirmed": 0,
                    "pooled_frac_confirmed": np.nan,
                    "pooled_median_abs_rho": np.nan,
                    "note": note,
                }
            )
            continue

        stacks: Dict[str, Dict[str, np.ndarray]] = {}
        sample_level: Dict[str, str] = {}
        for name in members:
            scheme = by_name[name]
            if not scheme.sample_key:
                raise ValueError(f"{name}: pooled_retest needs Scheme.sample_key.")
            present = [feature for feature in features if feature in set(scheme.features())]
            if not present:
                continue
            sub = scheme.adata[:, present]
            values = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
            frame = pd.DataFrame(values, index=sub.obs[scheme.sample_key].astype(str), columns=present)
            for sample, level in zip(sub.obs[scheme.sample_key].astype(str), sub.obs[scheme.group_key].astype(str)):
                if sample in sample_level and sample_level[sample] != level:
                    raise ValueError(f"Sample {sample} carries 2 levels: {sample_level[sample]} and {level}.")
                sample_level[sample] = level
            for sample, row in frame.iterrows():
                stacks.setdefault(sample, {})
                for feature, value in row.items():
                    stacks[sample].setdefault(feature, [])
                    stacks[sample][feature].append(float(value))

        samples = sorted(stacks)
        pooled = pd.DataFrame(index=pd.Index(samples), columns=pd.Index(features), dtype=float)
        for sample in samples:
            for feature, values in stacks[sample].items():
                pooled.at[sample, feature] = float(np.mean(values))
        pooled = pooled.dropna(axis=1, how="any")

        levels = [sample_level[sample] for sample in samples]
        level_counts = pd.Series(levels).value_counts()
        if pooled.shape[1] == 0 or level_counts.min() < int(min_samples_per_level) or len(level_counts) < 2:
            rows.append(
                {
                    "block": label,
                    "kind": kind,
                    "n_features": len(features),
                    "n_member_schemes": len(members),
                    "pooled_n_samples": len(samples),
                    "pooled_n_features_tested": int(pooled.shape[1]),
                    "pooled_n_features_confirmed": 0,
                    "pooled_frac_confirmed": np.nan,
                    "pooled_median_abs_rho": np.nan,
                    "note": "too few samples or features after pooling",
                }
            )
            continue

        r_df, p_df = marker_finder(pooled.to_numpy(dtype=float), levels, gene_names=pooled.columns)
        reference = by_name[members[0]].levels()[-1]
        rho_pooled = r_df[reference]
        p_pooled = p_df[reference]
        confirmed = int((p_pooled < float(alpha)).sum())
        rows.append(
            {
                "block": label,
                "kind": kind,
                "n_features": len(features),
                "n_member_schemes": len(members),
                "pooled_n_samples": len(samples),
                "pooled_n_features_tested": int(r_df.shape[0]),
                "pooled_n_features_confirmed": confirmed,
                "pooled_frac_confirmed": float(confirmed) / float(max(r_df.shape[0], 1)),
                "pooled_median_abs_rho": float(rho_pooled.abs().median()),
                "note": f"sample-level pooling over {len(members)} schemes; positive means {reference}",
            }
        )

    return pd.DataFrame(rows)


# ------------------------------------------------------------------------- stage 5: the heatmap


def build_consolidated_matrix(
    schemes: Sequence[Scheme],
    assignments: pd.DataFrame,
    blocks: pd.DataFrame,
    zscore_within_scheme: bool = True,
    column_grouping: str = "level",
) -> Tuple[pd.DataFrame, List[Tuple[str, int]], List[str], List[str], pd.DataFrame]:
    """Rows in block order. Columns grouped by level (the default) or by scheme.

    `column_grouping="level"` puts every sample of the first level together, then every sample of
    the second, so the broad contrast reads as one split down the middle and a scheme-specific
    block shows its difference inside its own scheme only. Within a level the columns stay in
    scheme order, and the scheme of each column goes on a covariate bar.

    `zscore_within_scheme` centers and scales each feature inside each scheme. Every block comes
    from a within-scheme contrast between the levels, so a z-score taken across all schemes would
    show the baseline difference between schemes instead of the contrast that defined the rows.
    Centering inside a scheme removes that baseline. It cannot create a difference between the
    levels, because it subtracts one number from all samples of the scheme.
    """
    if column_grouping not in ("level", "scheme"):
        raise ValueError("column_grouping must be 'level' or 'scheme'.")
    block_order = blocks.sort_values("display_order")["block"].tolist()
    ordered_features: List[str] = []
    row_blocks: List[str] = []
    for block in block_order:
        features = assignments.loc[assignments["block"] == block, "feature"].tolist()
        ordered_features.extend(features)
        row_blocks.extend([block] * len(features))

    columns: List[str] = []
    column_schemes: List[str] = []
    counts: List[Tuple[str, int]] = []
    covariate_rows = []
    blocks_matrix = []

    for scheme in schemes:
        levels = scheme.levels()
        obs = scheme.adata.obs
        order = sorted(
            range(scheme.adata.n_obs),
            key=lambda i: (levels.index(str(obs[scheme.group_key].iloc[i])), str(scheme.adata.obs_names[i])),
        )
        present = [feature for feature in ordered_features if feature in set(scheme.features())]
        sub = scheme.adata[:, present]
        values = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
        frame = pd.DataFrame(values.T, index=present, columns=list(sub.obs_names.astype(str)))
        frame = frame.reindex(index=ordered_features)
        frame = frame.iloc[:, order]
        names = [f"{scheme.name}|{name}" for name in frame.columns]
        frame.columns = names
        if zscore_within_scheme:
            frame = mh._zscore_rows(frame)
        blocks_matrix.append(frame)
        columns.extend(names)
        column_schemes.extend([scheme.name] * len(names))
        counts.append((scheme.name, len(names)))
        for position in order:
            covariate_rows.append(
                {
                    "column": f"{scheme.name}|{scheme.adata.obs_names[position]}",
                    scheme.group_key: str(obs[scheme.group_key].iloc[position]),
                    "sample": str(obs[scheme.sample_key].iloc[position]) if scheme.sample_key else "",
                }
            )

    matrix = pd.concat(blocks_matrix, axis=1)
    matrix = matrix.loc[ordered_features, columns]
    if matrix.isna().any().any():
        n_missing = int(matrix.isna().sum().sum())
        raise ValueError(
            f"The consolidated matrix holds {n_missing} missing cells. Every feature must carry a "
            "value in every scheme, or the figure would show a gap as a color."
        )

    covariate_df = pd.DataFrame(covariate_rows).set_index("column").loc[columns]
    group_key = schemes[0].group_key
    scheme_of_column = pd.Series(column_schemes, index=columns)

    if column_grouping == "level":
        levels = schemes[0].levels()
        rank = {level: position for position, level in enumerate(levels)}
        order = sorted(
            columns,
            key=lambda name: (rank.get(str(covariate_df.at[name, group_key]), len(levels)),
                              [scheme.name for scheme in schemes].index(scheme_of_column[name]),
                              name),
        )
        matrix = matrix[order]
        covariate_df = covariate_df.loc[order]
        column_groups = [str(covariate_df.at[name, group_key]) for name in order]
        counts = [(level, int(sum(1 for value in column_groups if value == level)))
                  for level in levels if any(value == level for value in column_groups)]
        covariate_df = covariate_df.copy()
        covariate_df.insert(0, "cell_state", [scheme_of_column[name] for name in order])
        covariate_df = covariate_df.drop(columns=[group_key])
    else:
        column_groups = list(column_schemes)

    if "sample" in covariate_df.columns and not covariate_df["sample"].astype(bool).all():
        covariate_df = covariate_df.drop(columns=["sample"])
    return matrix, counts, column_groups, row_blocks, covariate_df


def consolidated_heatmap(
    matrix: pd.DataFrame,
    output_path: str,
    cluster_counts: Sequence[Tuple[str, int]],
    column_group_order: Sequence[str],
    column_groups: Sequence[str],
    row_blocks: Sequence[str],
    block_order: Sequence[str],
    covariate_df: Optional[pd.DataFrame] = None,
    go_terms: Optional[dict] = None,
    go_terms_max: int = 12,
) -> None:
    """Draw with the established MarkerFinder renderer. No plotting code is added here.

    `_plot_heatmap` colors the column bar and the row bar from one label space, so the label list
    it receives is the column groups followed by the blocks. The x tick labels come from
    `cluster_counts`, which holds the column groups only.
    """
    order = list(dict.fromkeys(list(column_group_order) + list(block_order)))
    mh._plot_heatmap(
        matrix,
        output_path,
        list(cluster_counts),
        order,
        list(column_groups),
        list(row_blocks),
        covariate_df=covariate_df,
        go_terms=go_terms,
        go_terms_max=go_terms_max,
    )


# ---------------------------------------------------------------------------- stage 6: GO-Elite


def goelite_for_block(
    block: str,
    query_genes: Sequence[str],
    background_genes: Sequence[str],
    outdir: str,
    species: str = "human",
    repo_root: Optional[str] = None,
    log=print,
) -> pd.DataFrame:
    """Run the altanalyze3 GO-Elite CLI for one block. Returns its table, or an empty frame."""
    os.makedirs(outdir, exist_ok=True)
    safe = block.replace("|", "_").replace("+", "-").replace("/", "-")
    query_path = os.path.join(outdir, f"{safe}_query_genes.txt")
    background_path = os.path.join(outdir, "background_genes.txt")
    with open(query_path, "w") as handle:
        handle.write("\n".join(query_genes) + "\n")
    if not os.path.exists(background_path):
        with open(background_path, "w") as handle:
            handle.write("\n".join(background_genes) + "\n")

    repo_root = repo_root or os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    command = [
        sys.executable, "-m", "altanalyze3.components.goelite.main",
        "--species", species,
        "--query", query_path,
        "--background", background_path,
        "--outdir", outdir,
        "--name", safe,
    ]
    log(f"    [goelite] {' '.join(command)}")
    env = dict(os.environ, PYTHONPATH=repo_root)
    result = subprocess.run(command, cwd=repo_root, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"    [goelite][FAIL] exit {result.returncode}\n{result.stderr[-1500:]}")
        return pd.DataFrame()
    table = os.path.join(outdir, f"{safe}_goelite_results.tsv")
    if not os.path.exists(table):
        log(f"    [goelite][WARN] no results table at {table}")
        return pd.DataFrame()
    frame = pd.read_csv(table, sep="\t")
    log(f"    [goelite] {block}: {len(query_genes)} query genes, {len(background_genes)} background, "
        f"{len(frame)} terms")
    return frame


def terms_for_heatmap(frame: pd.DataFrame, max_terms: int) -> List[Tuple[str, float]]:
    """Rank by FDR then raw p, drop repeats, keep at most `max_terms`. The ICGS rule."""
    if frame is None or frame.empty or "term_name" not in frame.columns:
        return []
    table = frame.copy()
    table["p_value"] = pd.to_numeric(table["p_value"], errors="coerce")
    table["sort_fdr"] = pd.to_numeric(table["fdr"], errors="coerce") if "fdr" in table.columns else np.nan
    if "selected" in table.columns:
        chosen = table[table["selected"].astype(str).str.lower().isin(["true", "1"])]
        if not chosen.empty:
            table = chosen
    table = table.dropna(subset=["p_value"]).sort_values(
        ["sort_fdr", "p_value", "term_name"], na_position="last"
    )
    out: List[Tuple[str, float]] = []
    seen = set()
    for _, row in table.iterrows():
        name = str(row["term_name"]).strip()
        if not name or name in seen:
            continue
        seen.add(name)
        out.append((name, float(row["p_value"])))
        if len(out) >= int(max_terms):
            break
    return out


# ------------------------------------------------------------------------------- orchestration


def run_iterative_marker_finder(
    schemes: Sequence[Scheme],
    outdir: str,
    *,
    rho_min: float = 0.3,
    alpha: float = 0.05,
    use_rawp: bool = False,
    min_features_per_pattern: int = 10,
    level_names: Optional[Tuple[str, str]] = None,
    features: Optional[Sequence[str]] = None,
    prefix: str = "iterative_markerfinder",
    log=print,
) -> IterativeResults:
    """Stages 1 to 4, with every intermediate table written to `outdir`."""
    os.makedirs(outdir, exist_ok=True)
    scheme_order = [scheme.name for scheme in schemes]
    if level_names is None:
        first = schemes[0].levels()
        level_names = (first[0], first[-1])

    log(f"[iterative] {len(schemes)} schemes: {', '.join(scheme_order)}")
    rho, pval, fdr, testable, diagnostics = score_schemes(schemes, features=features)
    log(f"[iterative] score matrix {rho.shape[0]} features x {rho.shape[1]} contrast columns")
    for entry in diagnostics:
        log(f"  [diag] {entry}")

    signs = code_signs(rho, fdr, testable, rho_min=rho_min, alpha=alpha, pval=pval, use_rawp=use_rawp)
    n_significant = int((signs.abs() > 0).sum().sum())
    log(f"[iterative] gates |rho| >= {rho_min} and {'raw p' if use_rawp else 'FDR'} < {alpha}: "
        f"{n_significant} significant (feature, scheme) pairs of "
        f"{int(testable.sum().sum())} testable")

    assignments, blocks = assign_blocks(
        signs,
        level_names=level_names,
        min_features_per_pattern=min_features_per_pattern,
        scheme_order=list(rho.columns),
    )
    log(f"[iterative] {len(assignments)} features assigned to {len(blocks)} blocks")
    for kind in BLOCK_KIND_ORDER:
        subset = assignments[assignments["kind"] == kind]
        log(f"  [blocks] {kind}: {len(subset)} features of {len(assignments)} "
            f"({100.0 * len(subset) / max(len(assignments), 1):.1f}%)")

    # How many schemes could test each feature decides what its block means. A feature called
    # unique after 5 schemes tested it is a different claim from one only 2 schemes could test.
    coverage = (
        pd.crosstab(assignments["kind"], assignments["n_schemes_testable"])
        .reindex([kind for kind in BLOCK_KIND_ORDER if kind in set(assignments["kind"])])
    )
    coverage.to_csv(os.path.join(outdir, f"{prefix}_coverage_qc.tsv"), sep="\t")
    log("  [coverage] features by block kind (rows) and number of schemes able to test them (columns):")
    for line in coverage.fillna(0).astype(int).to_string().splitlines():
        log(f"    {line}")

    rho.to_csv(os.path.join(outdir, f"{prefix}_rho.tsv"), sep="\t")
    pval.to_csv(os.path.join(outdir, f"{prefix}_pvalue.tsv"), sep="\t")
    fdr.to_csv(os.path.join(outdir, f"{prefix}_fdr.tsv"), sep="\t")
    testable.to_csv(os.path.join(outdir, f"{prefix}_testable.tsv"), sep="\t")
    signs.to_csv(os.path.join(outdir, f"{prefix}_signs.tsv"), sep="\t")
    assignments.to_csv(os.path.join(outdir, f"{prefix}_assignments.tsv"), sep="\t", index=False)
    blocks.to_csv(os.path.join(outdir, f"{prefix}_blocks.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, f"{prefix}_diagnostics.json"), "w") as handle:
        json.dump(diagnostics, handle, indent=2)

    return IterativeResults(
        rho=rho,
        pval=pval,
        fdr=fdr,
        testable=testable,
        signs=signs,
        assignments=assignments,
        blocks=blocks,
        diagnostics=diagnostics,
    )


def run_broad_then_specific(
    schemes: Sequence[Scheme],
    outdir: str,
    *,
    rho_min: float = 0.4,
    alpha: float = 0.05,
    use_rawp: bool = False,
    min_features_per_pattern: int = 10,
    level_names: Optional[Tuple[str, str]] = None,
    features: Optional[Sequence[str]] = None,
    min_samples_per_level: Optional[Dict[str, int]] = None,
    prefix: str = "iterative_markerfinder",
    log=print,
) -> dict:
    """The broad contrast first, then the scheme-specific effects that survive it.

    MarkerFinder runs once on the whole data set, one row per sample, to find the broad difference
    between the 2 levels. The broad effect then comes off every scheme, and MarkerFinder runs again
    inside each scheme on what is left. A feature still significant there carries a difference the
    broad program does not explain.
    """
    os.makedirs(outdir, exist_ok=True)
    if level_names is None:
        first = schemes[0].levels()
        level_names = (first[0], first[-1])
    levels = list(level_names)

    log(f"[broad] MarkerFinder on the broad {levels[0]} versus {levels[1]} contrast, "
        f"over all {len(schemes)} schemes at once")
    broad = broad_contrast(schemes, levels, features=features,
                           min_samples_per_level=min_samples_per_level)
    log(f"  [broad] {broad['n_samples']} samples "
        f"({', '.join(f'{level} {count}' for level, count in broad['n_samples_per_level'].items())}); "
        f"{broad['n_testable']} of {len(broad['rho'])} features testable; "
        f"{broad['n_constant']} constant; {broad['n_filled_cells']} filled cells")

    metric = broad["pval"] if use_rawp else broad["fdr"]
    broad_significant = (broad["rho"].abs() >= float(rho_min)) & (metric < float(alpha))
    broad_signs = pd.Series(0.0, index=broad["rho"].index)
    broad_signs = broad_signs.where(~broad_significant, np.sign(broad["rho"]))
    broad_signs = broad_signs.where(broad["testable"], np.nan)
    log(f"  [broad] {int(broad_significant.sum())} features clear |rho| >= {rho_min} and "
        f"{'raw p' if use_rawp else 'FDR'} < {alpha} in the broad contrast")

    log("[residual] removing the broad effect from every scheme, then scoring what is left")
    residual = residual_schemes(schemes, broad["effect"], positive_level=levels[-1])
    r_rho, r_pval, r_fdr, r_testable, diagnostics = score_schemes(
        residual, features=list(broad["rho"].index))
    r_signs = code_signs(r_rho, r_fdr, r_testable, rho_min=rho_min, alpha=alpha,
                         pval=r_pval, use_rawp=use_rawp)
    log(f"  [residual] {int((r_signs.abs() > 0).sum().sum())} significant (feature, scheme) pairs "
        f"of {int(r_testable.sum().sum())} testable")

    assignments, blocks = assign_blocks_broad_first(
        broad_signs, r_signs, level_names=level_names,
        min_features_per_pattern=min_features_per_pattern,
        scheme_order=list(r_rho.columns),
    )
    log(f"[blocks] {len(assignments)} features assigned to {len(blocks)} blocks")
    for kind in ["broad", "shared-specific", "scheme-specific", "discordant"]:
        subset = assignments[assignments["kind"] == kind]
        log(f"  [blocks] {kind}: {len(subset)} features of {len(assignments)} "
            f"({100.0 * len(subset) / max(len(assignments), 1):.1f}%)")
    amplified = assignments[(assignments["kind"] == "broad") & (assignments["also_specific_in"] != "")]
    log(f"  [blocks] broad features that ALSO keep a scheme-specific residual: {len(amplified)} of "
        f"{int((assignments['kind'] == 'broad').sum())}")

    broad_table = pd.DataFrame(
        {
            "rho": broad["rho"],
            "pval": broad["pval"],
            "fdr": broad["fdr"],
            "effect": broad["effect"],
            "testable": broad["testable"],
            "sign": broad_signs,
        }
    )
    broad_table.to_csv(os.path.join(outdir, f"{prefix}_broad_contrast.tsv"), sep="\t")
    broad["pooled"].to_csv(os.path.join(outdir, f"{prefix}_broad_sample_matrix.tsv"), sep="\t")
    r_rho.to_csv(os.path.join(outdir, f"{prefix}_residual_rho.tsv"), sep="\t")
    r_fdr.to_csv(os.path.join(outdir, f"{prefix}_residual_fdr.tsv"), sep="\t")
    r_signs.to_csv(os.path.join(outdir, f"{prefix}_residual_signs.tsv"), sep="\t")
    assignments.to_csv(os.path.join(outdir, f"{prefix}_assignments.tsv"), sep="\t", index=False)
    blocks.to_csv(os.path.join(outdir, f"{prefix}_blocks.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, f"{prefix}_residual_diagnostics.json"), "w") as handle:
        json.dump(diagnostics, handle, indent=2)

    return {
        "broad": broad,
        "broad_signs": broad_signs,
        "residual_schemes": residual,
        "residual_rho": r_rho,
        "residual_pval": r_pval,
        "residual_fdr": r_fdr,
        "residual_testable": r_testable,
        "residual_signs": r_signs,
        "assignments": assignments,
        "blocks": blocks,
        "diagnostics": diagnostics,
    }


def schemes_from_h5ad(
    h5ad_path: str,
    scheme_col: str,
    group_key: str,
    level_order: Sequence[str],
    sample_key: Optional[str] = None,
) -> List[Scheme]:
    """Split one h5ad into a scheme per level of `scheme_col`."""
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path)
    for column in [scheme_col, group_key] + ([sample_key] if sample_key else []):
        if column not in adata.obs.columns:
            raise KeyError(f"obs column '{column}' absent from {h5ad_path}")
    schemes = []
    for name in [str(value) for value in pd.unique(adata.obs[scheme_col].astype(str))]:
        sub = adata[adata.obs[scheme_col].astype(str) == name].copy()
        schemes.append(
            Scheme(name=name, adata=sub, group_key=group_key, level_order=list(level_order),
                   sample_key=sample_key)
        )
    return schemes


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Consolidated MarkerFinder across grouping schemes: unique, shared and global programs."
    )
    parser.add_argument("--h5ad", required=True, help="Input h5ad, samples x features, no NaN in X.")
    parser.add_argument("--scheme-col", required=True, help="obs column that splits the data into schemes.")
    parser.add_argument("--group-key", required=True, help="obs column holding the levels of each scheme.")
    parser.add_argument("--levels", required=True, help="Comma-delimited level order, e.g. young,aged.")
    parser.add_argument("--sample-col", default=None, help="obs column naming the biological sample.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--rho", type=float, default=0.3, help="Absolute MarkerFinder rho gate (default 0.3).")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance gate (default 0.05).")
    parser.add_argument("--use-rawp", action="store_true", help="Gate on the raw p-value, not the BH value.")
    parser.add_argument("--min-features-per-pattern", type=int, default=10)
    parser.add_argument(
        "--mode",
        choices=["broad-first", "per-scheme"],
        default="broad-first",
        help=("'broad-first' (default) scores the broad contrast over every scheme, then the "
              "scheme-specific effects that survive it. 'per-scheme' scores each scheme on its "
              "own and reports unique, shared and global patterns."),
    )
    parser.add_argument("--heatmap", default=None, help="Optional consolidated heatmap PDF path.")
    args = parser.parse_args(argv)

    levels = [value.strip() for value in args.levels.split(",") if value.strip()]
    schemes = schemes_from_h5ad(args.h5ad, args.scheme_col, args.group_key, levels, args.sample_col)

    if args.mode == "broad-first":
        outcome = run_broad_then_specific(
            schemes, args.outdir, rho_min=args.rho, alpha=args.alpha, use_rawp=args.use_rawp,
            min_features_per_pattern=args.min_features_per_pattern,
            level_names=(levels[0], levels[-1]),
        )
        assignments, blocks = outcome["assignments"], outcome["blocks"]
        draw_schemes = schemes
    else:
        results = run_iterative_marker_finder(
            schemes, args.outdir, rho_min=args.rho, alpha=args.alpha, use_rawp=args.use_rawp,
            min_features_per_pattern=args.min_features_per_pattern,
            level_names=(levels[0], levels[-1]),
        )
        assignments, blocks = results.assignments, results.blocks
        draw_schemes = schemes

    if args.heatmap:
        matrix, counts, column_groups, row_blocks, covariate_df = build_consolidated_matrix(
            draw_schemes, assignments, blocks, column_grouping="level"
        )
        consolidated_heatmap(
            matrix,
            args.heatmap,
            counts,
            levels,
            column_groups,
            row_blocks,
            blocks.sort_values("display_order")["block"].tolist(),
            covariate_df=covariate_df,
        )
        print(f"[write] {args.heatmap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
