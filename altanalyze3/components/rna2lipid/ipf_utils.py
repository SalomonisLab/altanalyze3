"""Shared RNA metadata, expression mapping and donor-level statistical helpers."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from scipy.stats import mannwhitneyu, spearmanr, t as student_t
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def metadata(path):
    fields = {}
    with gzip.open(path, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            if row[0] == "!Sample_title":
                titles = row[1:]
            elif row[0] == "!Sample_geo_accession":
                fields["geo_accession"] = row[1:]
            elif row[0] == "!Sample_characteristics_ch1":
                parts = [s.split(": ", 1) for s in row[1:]]
                if len({s[0] for s in parts}) != 1:
                    raise ValueError("Inconsistent GEO metadata field")
                fields[parts[0][0]] = [s[1] for s in parts]
    result = pd.DataFrame(fields, index=titles)
    result.index.name = "sample"
    if not result.index.is_unique:
        raise ValueError("Ambiguous RNA sample IDs")
    return result


def canonical(name):
    """Acyl-chain composition only; preserve ether type and sphingobase.

    Ignores chain position, double-bond stereochemistry and iso index. Ambiguous
    multiple peaks are subsequently excluded, never summed or averaged.
    """
    if not isinstance(name, str):
        return None
    name = re.sub(r"_[PN]$", "", name.strip())
    name = re.sub(r"\[iso\d+\]", "", name)
    name = re.sub(r"\(\d+[EZ](?:,\d+[EZ])*\)", "", name)
    match = re.fullmatch(r"(LysoPE|LysoPC|LPE|LPC|PC|PE|PI|PG|PS|PA|DG|TG|Cer|SM)\s*\(?([^()]+)\)?", name)
    if not match:
        return None
    cls, composition = match.groups()
    cls = {"LysoPE": "LPE", "LysoPC": "LPC"}.get(cls, cls)
    chains = composition.split("/")
    if not all(re.fullmatch(r"(?:[OP]-|d)?\d+:\d+(?:;O2)?", s) for s in chains):
        return None
    # d18:1 and 18:1;O2 denote the same dihydroxy sphingoid base.
    if cls in {"Cer", "SM"}:
        chains = ["d" + s[:-3] if s.endswith(";O2") else s for s in chains]
        if len(chains) != 2 or not any(s.startswith("d") for s in chains):
            return None
    if cls in {"PC", "PE"} and "0:0" in chains:
        cls = "L" + cls
    chains = [s for s in chains if s != "0:0"]
    # A summed species such as TG(56:3) cannot match a resolved three-chain TG.
    return cls + "(" + "/".join(sorted(chains)) + ")"


def map_expression(frame, mapping, logged):
    symbols = pd.Index([mapping.get(str(g).split(".")[0], None) for g in frame.index])
    keep = symbols.notna()
    values = frame.loc[keep].copy()
    values.index = symbols[keep]
    # Multiple Ensembl IDs mapping to one symbol are aggregated in linear space.
    values = np.exp2(values) if logged else values
    values = values.groupby(level=0).sum()
    return (np.log2(values) if logged else values).T


def donor_means(frame, meta, region=None):
    chosen = meta.diseasegroup.isin(["IPF", "NDC"]) & meta.lunglocation.isin(["Apex", "Base"])
    if region:
        chosen &= meta.lunglocation.eq(region)
    rows = meta.loc[chosen]
    # Average technical/anatomical left-right profiles within each region;
    # give apex/base equal weight within a donor, then each donor equal weight.
    values = frame.loc[rows.index].copy()
    values.index = pd.MultiIndex.from_arrays([rows.donorid, rows.lunglocation], names=["donor", "region"])
    return values.groupby(level=[0, 1]).mean().groupby(level=0).mean()


def welch_pvalues(case, control):
    case, control = np.asarray(case), np.asarray(control)
    pooled_variance = case.var(axis=0, ddof=1) + control.var(axis=0, ddof=1)
    magnitude = np.maximum(np.maximum(abs(case).max(axis=0), abs(control).max(axis=0)), 1.)
    testable = pooled_variance > np.finfo(float).eps * magnitude ** 2
    p = np.ones(case.shape[1])
    if testable.any():
        a, b = case[:, testable], control[:, testable]
        va, vb = a.var(axis=0, ddof=1) / len(a), b.var(axis=0, ddof=1) / len(b)
        statistic = (a.mean(axis=0) - b.mean(axis=0)) / np.sqrt(va + vb)
        degrees = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        p[testable] = 2 * student_t.sf(abs(statistic), degrees)
    return np.where(np.isfinite(p), p, 1.), testable


def fold_table(linear, log_values, groups, rng, bootstraps):
    case = linear.loc[groups.eq("IPF")].to_numpy()
    ctrl = linear.loc[groups.eq("NDC")].to_numpy()
    a, b = case.mean(axis=0), ctrl.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        fc = np.log2(a / b)
    # Multinomial subject resampling avoids large B x subjects x lipids arrays.
    wa = rng.multinomial(len(case), np.ones(len(case)) / len(case), size=bootstraps) / len(case)
    wb = rng.multinomial(len(ctrl), np.ones(len(ctrl)) / len(ctrl), size=bootstraps) / len(ctrl)
    with np.errstate(divide="ignore", invalid="ignore"):
        draws = np.log2((wa @ case) / (wb @ ctrl))
    finite_draws = np.isfinite(draws).mean(axis=0)
    lo, hi = np.full(len(a), np.nan), np.full(len(a), np.nan)
    for j in range(len(a)):
        # Require all resampled group means to be positive for a finite CI.
        if finite_draws[j] == 1:
            lo[j], hi[j] = np.quantile(draws[:, j], [.025, .975])
    if log_values is None:
        p = mannwhitneyu(case, ctrl, axis=0, method="asymptotic").pvalue
        geometric = np.full(len(a), np.nan)
    else:
        ca = log_values.loc[groups.eq("IPF")].to_numpy()
        co = log_values.loc[groups.eq("NDC")].to_numpy()
        p, _ = welch_pvalues(ca, co)
        geometric = ca.mean(axis=0) - co.mean(axis=0)
    p = np.where(np.isfinite(p), p, 1.)
    result = pd.DataFrame({"IPF_mean": a, "control_mean": b, "log2FC": fc,
                           "fold_change": np.exp2(fc), "log2FC_CI_low": lo, "log2FC_CI_high": hi,
                           "bootstrap_finite_fraction": finite_draws,
                           "geometric_log2FC": geometric, "pvalue": p,
                           "FDR": multipletests(p, method="fdr_bh")[1],
                           "IPF_detection": (case > 0).mean(axis=0),
                           "control_detection": (ctrl > 0).mean(axis=0)}, index=linear.columns)
    result.index.name = "feature"
    return result


def concordance(table):
    use = table[np.isfinite(table.measured_log2FC) & np.isfinite(table.predicted_log2FC)]
    if len(use) < 2:
        return {"lipids": len(use), "sufficient_overlap": False}
    a, b = use.measured_log2FC.to_numpy(), use.predicted_log2FC.to_numpy()
    substantial = abs(a) >= .5
    slope = float(np.dot(a, b) / np.dot(a, a)) if np.dot(a, a) else None
    result = {"lipids": len(use), "Pearson": float(np.corrcoef(a, b)[0, 1]),
              "Spearman": float(spearmanr(a, b).statistic),
              "RMSE_log2FC": float(np.sqrt(np.mean((a - b) ** 2))),
              "median_absolute_log2FC_error": float(np.median(abs(a - b))),
              "slope_through_origin_predicted_on_measured": slope,
              "direction_agreement": float(np.mean(np.sign(a) == np.sign(b))),
              "substantial_measured_changes": int(substantial.sum()),
              "substantial_direction_agreement": float(np.mean(np.sign(a[substantial]) == np.sign(b[substantial]))) if substantial.any() else None,
              "within_0_5_log2FC": int((abs(a - b) <= .5).sum())}
    # Descriptive acceptance gate, frozen in source before IPF predictions.
    result["magnitude_replication_gate"] = bool(len(use) >= 5 and substantial.sum() >= 5 and
        result["RMSE_log2FC"] <= .5 and slope is not None and .8 <= slope <= 1.25 and
        result["substantial_direction_agreement"] >= .8)
    return result

