"""
variant2clone -- map single-cell variant calls onto inferCNV / fastCNV CNV results
to (1) find the CNV genomic interval(s) most associated with each mutation,
(2) cross-validate those associations across independent genotyping assays
(e.g. GoT vs long-read), and (3) divide cells into discrete CNV clones using the
validated loci as a prior.

Design notes
------------
* Signal = MUTATION PRESENCE. A MUT read means the cell IS mutant; WT reads are
  not used as the negative class, because heterozygous cells carry both alleles
  and either is detected by chance at low sensitivity. The contrast is therefore
  "mutation detected" vs "all other cells".
* Assays are independent measurements of the same mutation: a real CNV locus must
  associate with the mutation in BOTH assays (e.g. GoT and LR). Cross-assay
  concordance is the validation, and it is robust to assay-specific artefacts.
* Multi-scale sliding windows (small and large) so focal associations are not
  merged away by arm-level averaging.
* Test = cell-type-adjusted point-biserial correlation (signed: r<0 -> the
  mutation tracks DOWN-regulation/loss; r>0 -> UP-regulation/gain). Adjusting on
  the reference cell type removes cell-type expression programs; the raw value is
  also returned because a true clonal CNV that coincides with a cell state can be
  partially removed by the adjustment.

Inputs are deliberately generic (a residual/CNV matrix + a gene-order table + a
per-cell variant table), so the same module runs on inferCNV's
``infercnv.observations`` matrix, on a pyInferCNV ``run_infercnv`` residual, or on
a fastCNV scored matrix.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


# --------------------------------------------------------------------------- #
# parameters
# --------------------------------------------------------------------------- #
@dataclass
class Variant2CloneParams:
    # The PRIMARY assay drives locus discovery and clone assignment (long-read, LR).
    # VALIDATION assays (e.g. GoT) only confirm the primary associations and are
    # usually ABSENT in production -- the workflow runs on the primary alone and
    # simply reports that validation was unavailable.
    primary_assay: str = "LR"
    validation_assays: Tuple[str, ...] = ("GoT",)
    window_sizes: Tuple[int, ...] = (3, 8, 20, 50)   # interval sizes in genes (small -> large)
    min_mut_cells: int = 20                          # skip a (gene, assay) with fewer MUT cells
    concordance_p: float = 0.01                      # validation-assay p needed to confirm a locus
    fdr_alpha: float = 0.05
    adjust_celltype: bool = False                    # DO NOT adjust on cell type by default: the
    #   malignant clone is itself a cell state, so adjusting removes the real clonal CNV signal
    #   (e.g. chr6). Cross-assay (GoT vs LR) concordance is the validation instead.
    n_clones: int = 0                                # 0 = auto (n_genes_with_locus + 2, capped 6)
    center_value: float = 1.0                        # neutral value of the CNV matrix (1.0 for FC)
    seed: int = 0


def _match_assay(assays, name):
    """Return the assay in `assays` matching `name` (case-insensitive), else None."""
    nl = str(name).lower()
    for a in assays:
        if str(a).lower() == nl or nl in str(a).lower():
            return a
    return None


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def infer_mutation_specs(columns) -> Dict[str, Dict[str, str]]:
    """Parse variant-table column names into {gene: {assay: column}} for MUT columns.

    Recognises ``GENE-MUT-ASSAY`` (e.g. SRSF2-MUT-GoT), ``GENE_MUT_ASSAY`` and
    bare ``GENE-MUT`` (assay 'NA'). WT columns are ignored (see module docstring).
    """
    specs: Dict[str, Dict[str, str]] = {}
    for col in columns:
        toks = re.split(r"[-_]", str(col))
        up = [t.upper() for t in toks]
        if "MUT" not in up:
            continue
        mi = up.index("MUT")
        gene = "-".join(toks[:mi]) or str(col)
        assay = "-".join(toks[mi + 1:]) or "NA"
        specs.setdefault(gene, {})[assay] = col
    return specs


def _ct_resid(M: np.ndarray, codes: np.ndarray, K: int) -> np.ndarray:
    """Center each column of M (n x p) or vector (n,) by cell-type group mean."""
    cnt = np.bincount(codes, minlength=K).astype(float)
    if M.ndim == 1:
        return M - (np.bincount(codes, weights=M, minlength=K) / cnt)[codes]
    gsum = np.zeros((K, M.shape[1]))
    np.add.at(gsum, codes, M)
    return M - (gsum / cnt[:, None])[codes]


def _point_biserial(S: np.ndarray, y: np.ndarray, codes: Optional[np.ndarray], K: int):
    """Signed point-biserial r + two-sided p of every column of S against binary y.

    If codes is given, both S and y are centered within cell type first
    (cell-type-adjusted partial correlation). Returns (r, p)."""
    if codes is not None:
        Sc = _ct_resid(S, codes, K)
        yc = _ct_resid(y.astype(float), codes, K)
        dof = len(y) - K - 1
    else:
        Sc = S - S.mean(0)
        yc = y.astype(float) - y.mean()
        dof = len(y) - 2
    r = (Sc.T @ yc) / (np.sqrt((Sc * Sc).sum(0) * (yc * yc).sum()) + 1e-12)
    dof = max(dof, 1)
    t = r * np.sqrt(dof / np.maximum(1 - r * r, 1e-12))
    return r, 2 * stats.t.sf(np.abs(t), df=dof)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = (p[order] * n / (np.arange(n) + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# --------------------------------------------------------------------------- #
# 1. break the genome into multi-scale CNV intervals
# --------------------------------------------------------------------------- #
def breakdown_intervals(residual: pd.DataFrame, gene_order: pd.DataFrame,
                        params: Variant2CloneParams = Variant2CloneParams()):
    """Multi-scale sliding windows -> (scores, intervals).

    residual   : genes x cells DataFrame (CNV/FC values), index = gene names.
    gene_order : DataFrame indexed by gene with columns 'chr','start','stop'.
    Returns (scores: cells x intervals float32 ndarray, intervals: DataFrame with
    chr/g0/g1/size_genes/start/stop/genes).
    """
    genes = [g for g in residual.index if g in gene_order.index]
    R = residual.loc[genes].values.astype(np.float32)            # genes x cells
    go = gene_order.loc[genes]
    gchr = go["chr"].values
    gstart = go["start"].values.astype(float)
    gstop = go["stop"].values.astype(float)
    gname = np.asarray(genes)
    nC = R.shape[1]

    cols, meta = [], []
    for ch in pd.unique(gchr):
        gi = np.where(gchr == ch)[0]
        gi = gi[np.argsort(gstart[gi])]
        pref = np.vstack([np.zeros((1, nC)), np.cumsum(R[gi, :], axis=0)])
        for size in params.window_sizes:
            if size > len(gi):
                continue
            step = max(2, size // 2)
            for s in range(0, len(gi) - size + 1, step):
                e = s + size
                cols.append(((pref[e] - pref[s]) / size).astype(np.float32))
                g0, g1 = gi[s], gi[e - 1]
                meta.append((ch, int(g0), int(g1), size,
                             float(gstart[g0]), float(gstop[g1]),
                             ";".join(gname[gi[s:e]][:8])))
    scores = np.array(cols, dtype=np.float32).T                 # cells x intervals
    intervals = pd.DataFrame(meta, columns=["chr", "g0", "g1", "size_genes", "start", "stop", "genes"])
    return scores, intervals


# --------------------------------------------------------------------------- #
# 2. associate every variant (per assay) with every interval
# --------------------------------------------------------------------------- #
def associate_variants(scores: np.ndarray, intervals: pd.DataFrame, variant_df: pd.DataFrame,
                       cell_names, mut_specs: Optional[Dict[str, Dict[str, str]]] = None,
                       cell_types=None, params: Variant2CloneParams = Variant2CloneParams()):
    """Cell-type-adjusted point-biserial of each interval vs MUT-presence, per assay.

    Returns a long DataFrame with one row per (gene, assay, interval):
    r, p (cell-type-adjusted), r_raw, p_raw, fdr, direction.
    """
    if mut_specs is None:
        mut_specs = infer_mutation_specs(variant_df.columns)
    vd = variant_df.reindex(list(cell_names)).fillna(0)
    if cell_types is not None and params.adjust_celltype:
        cats = pd.Categorical(pd.Series(list(cell_types)))
        codes, K = cats.codes.astype(int), len(cats.categories)
    else:
        codes, K = None, 0

    blocks = []
    for gene, assays in mut_specs.items():
        for assay, col in assays.items():
            y = (pd.to_numeric(vd[col], errors="coerce").fillna(0).values > 0)
            if y.sum() < params.min_mut_cells:
                continue
            r, p = _point_biserial(scores, y, codes, K)
            r0, p0 = _point_biserial(scores, y, None, 0)
            b = intervals.copy()
            b["gene"], b["assay"] = gene, assay
            b["r"], b["p"], b["r_raw"], b["p_raw"] = r, p, r0, p0
            b["n_mut"] = int(y.sum())
            blocks.append(b)
    assoc = pd.concat(blocks, ignore_index=True) if blocks else pd.DataFrame()
    if len(assoc):
        assoc["fdr"] = _bh_fdr(assoc["p"].values)
        assoc["direction"] = np.where(assoc["r"] < 0, "down(loss)", "up(gain)")
        assoc["start_Mb"] = (assoc["start"] / 1e6).round(2)
        assoc["stop_Mb"] = (assoc["stop"] / 1e6).round(2)
    return assoc


# --------------------------------------------------------------------------- #
# 2b. MATRIX-LEVEL association: CNV intervals as features predicting each mutation
# --------------------------------------------------------------------------- #
def associate_matrix(scores: np.ndarray, intervals: pd.DataFrame, variant_df: pd.DataFrame,
                     cell_names, mut_specs: Optional[dict] = None, cell_types=None,
                     params: Variant2CloneParams = Variant2CloneParams(), C: float = 0.05):
    """Treat the CNV intervals as a feature matrix (cells x intervals) and fit, per
    mutation, an L1-regularized logistic model predicting MUT-presence. This selects
    the *combination* of CNV loci associated with the mutation (not just the single
    best), and yields a per-cell signature score used downstream for clone calling.

    Validation: the model is trained on the PRIMARY (LR) label with cross-validated
    AUC; if a VALIDATION (GoT) label exists, the same out-of-fold predictions are
    scored against it (independent-assay AUC). Returns (model_summary, signatures)
    where signatures[gene] is the per-cell decision score (clone prior).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    if mut_specs is None:
        mut_specs = infer_mutation_specs(variant_df.columns)
    vd = variant_df.reindex(list(cell_names)).fillna(0)
    X = scores.astype(float)
    if params.adjust_celltype and cell_types is not None:
        cats = pd.Categorical(pd.Series(list(cell_types)))
        X = _ct_resid(X, cats.codes.astype(int), len(cats.categories))
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)

    rows, signatures = [], {}
    for gene, assays in mut_specs.items():
        primary = _match_assay(list(assays), params.primary_assay) or list(assays)[0]
        validation = next((_match_assay(list(assays), v) for v in params.validation_assays
                           if _match_assay(list(assays), v) and _match_assay(list(assays), v) != primary), None)
        y = (pd.to_numeric(vd[assays[primary]], errors="coerce").fillna(0).values > 0)
        if y.sum() < params.min_mut_cells:
            continue
        clf = LogisticRegression(penalty="l1", C=C, solver="liblinear",
                                 class_weight="balanced", max_iter=300)
        try:
            cv = StratifiedKFold(5, shuffle=True, random_state=params.seed)
            proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
            auc_cv = float(roc_auc_score(y, proba))
        except Exception:
            proba, auc_cv = np.full(len(y), y.mean()), np.nan
        clf.fit(X, y)
        signatures[gene] = clf.decision_function(X)
        coef = clf.coef_[0]
        nz = np.where(np.abs(coef) > 1e-8)[0]
        top = sorted(nz, key=lambda j: -abs(coef[j]))[:12]
        feats = [{"chr": intervals.iloc[j]["chr"],
                  "start_Mb": round(intervals.iloc[j]["start"] / 1e6, 2),
                  "stop_Mb": round(intervals.iloc[j]["stop"] / 1e6, 2),
                  "size_genes": int(intervals.iloc[j]["size_genes"]),
                  "coef": round(float(coef[j]), 3), "dir": "gain" if coef[j] > 0 else "loss"} for j in top]
        auc_val = np.nan
        if validation is not None:
            yv = (pd.to_numeric(vd[assays[validation]], errors="coerce").fillna(0).values > 0)
            if yv.sum() >= 5:
                auc_val = float(roc_auc_score(yv, proba))   # LR-trained OOF predictions vs GoT label
        rows.append({"gene": gene, "primary_assay": primary, "n_mut": int(y.sum()),
                     "cv_auc_primary": round(auc_cv, 3) if auc_cv == auc_cv else np.nan,
                     "validation_assay": validation,
                     "validation_auc": round(auc_val, 3) if auc_val == auc_val else np.nan,
                     "n_selected_features": int(len(nz)),
                     "top_features": "; ".join(f"{f['chr']}:{f['start_Mb']}-{f['stop_Mb']}Mb({f['dir']},{f['coef']:+})" for f in feats[:6])})
    return pd.DataFrame(rows), signatures


def assign_clones_matrix(signatures: dict, variant_df: pd.DataFrame, cell_names,
                         mut_specs=None, params: Variant2CloneParams = Variant2CloneParams()):
    """Combinatorial clones from the per-mutation matrix-model signatures (cells x
    mutations). Clusters jointly so a cell high on SRSF2 *and* RUNX1 signatures forms
    a multi-mutation clone. Annotates each clone with per-assay mutation frequencies."""
    from sklearn.cluster import KMeans
    genes = list(signatures)
    if not genes:
        per = pd.DataFrame({"cell": list(cell_names), "clone": "C0"})
        return pd.DataFrame([{"clone": "C0", "n_cells": len(cell_names)}]), per
    Z = np.column_stack([signatures[g] for g in genes])
    Zz = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    k = params.n_clones or min(8, max(3, len(genes) + 2))
    labels = KMeans(n_clusters=k, random_state=params.seed, n_init=10).fit_predict(Zz)
    if mut_specs is None:
        mut_specs = infer_mutation_specs(variant_df.columns)
    vd = variant_df.reindex(list(cell_names)).fillna(0)
    prof = []
    for c in range(k):
        m = labels == c
        rec = {"clone": f"C{c}", "n_cells": int(m.sum())}
        for g in genes:
            rec[f"{g}_signature_mean"] = round(float(Z[m, genes.index(g)].mean()), 3)
        for gene, assays in mut_specs.items():
            for asy, col in assays.items():
                rec[f"{gene}_{asy}_pct"] = round(100 * (pd.to_numeric(vd[col], errors="coerce").fillna(0).values[m] > 0).mean(), 1)
        prof.append(rec)
    return pd.DataFrame(prof).sort_values("n_cells", ascending=False), \
        pd.DataFrame({"cell": list(cell_names), "clone": [f"C{c}" for c in labels]})


# --------------------------------------------------------------------------- #
# 2c. fixed N-gene window genotype scan (find a high-purity, sub-clonal CNV marker)
# --------------------------------------------------------------------------- #
def scan_genotype_windows(residual: pd.DataFrame, gene_order: pd.DataFrame,
                          mut_mask, wt_mask, ref_mask=None, window: int = 3,
                          ref_percentile: float = 95.0, min_mut_pos: int = 15):
    """Slide a fixed N-gene window across the genome over the pyInferCNV/inferCNV
    residual and correlate CNV-positive status to a binary genotype (e.g. believable
    GoT-LR-matching SRSF2 MUT vs WT).

    For each window x direction (gain/loss), CNV+ is called from the REFERENCE
    (healthy) distribution: gain if score > ``ref_percentile`` of the reference,
    loss if score < (100 - ref_percentile). Reports, per window:
      pct_MUT_CNVpos       % of all MUT cells that are CNV+   (denominator = all MUT)
      pct_WT_CNVpos        % of all WT cells that are CNV+    (denominator = all WT)
      purity_MUTfrac       MUT / (MUT+WT) among genotyped CNV+ cells
      frac_allcells_CNVpos fraction of query cells CNV+       (sub-clonality)
      OR, p                Fisher (MUT enriched in CNV+ vs WT)

    residual : genes x cells DataFrame; gene_order indexed by gene (chr/start/stop).
    mut_mask, wt_mask : per-cell booleans (genotype). ref_mask : per-cell boolean of
    the healthy reference cells used for the threshold (None => use all cells/proxy).
    Returns a DataFrame (one row per window x direction with >= min_mut_pos MUT CNV+).
    """
    genes = [g for g in residual.index if g in gene_order.index]
    R = residual.loc[genes].values.astype(np.float32)
    go = gene_order.loc[genes]; gchr = go["chr"].values
    gstart = go["start"].values.astype(float); gstop = go["stop"].values.astype(float)
    gname = np.asarray(genes)
    mut_mask = np.asarray(mut_mask, bool); wt_mask = np.asarray(wt_mask, bool)
    has_ref = ref_mask is not None and np.asarray(ref_mask, bool).any()
    ref_mask = np.asarray(ref_mask, bool) if has_ref else np.zeros(R.shape[1], bool)
    query_mask = ~ref_mask if has_ref else np.ones(R.shape[1], bool)
    nMUT, nWT, nQ = int(mut_mask.sum()), int(wt_mask.sum()), int(query_mask.sum())
    rows = []
    for ch in pd.unique(gchr):
        gi = np.where(gchr == ch)[0]; gi = gi[np.argsort(gstart[gi])]
        if len(gi) < window:
            continue
        pref = np.vstack([np.zeros((1, R.shape[1]), np.float32), np.cumsum(R[gi, :], axis=0)])
        for s in range(0, len(gi) - window + 1):
            e = s + window; sc = (pref[e] - pref[s]) / window
            ref = sc[ref_mask] if has_ref else sc
            hi = np.percentile(ref, ref_percentile); lo = np.percentile(ref, 100 - ref_percentile)
            for direction, pos in [("gain", sc > hi), ("loss", sc < lo)]:
                nmp = int((pos & mut_mask).sum()); nwp = int((pos & wt_mask).sum())
                if nmp < min_mut_pos:
                    continue
                npos = int((pos & query_mask).sum())
                pur = nmp / (nmp + nwp) if (nmp + nwp) else np.nan
                orr, p = stats.fisher_exact([[nmp, nMUT - nmp], [nwp, nWT - nwp]], alternative="greater")
                rows.append({"chr": ch, "start_Mb": round(gstart[gi[s]] / 1e6, 2),
                             "stop_Mb": round(gstop[gi[e - 1]] / 1e6, 2), "direction": direction,
                             "genes": ";".join(gname[gi[s:e]]),
                             "pct_MUT_CNVpos": round(100 * nmp / max(nMUT, 1), 1),
                             "pct_WT_CNVpos": round(100 * nwp / max(nWT, 1), 1),
                             "purity_MUTfrac": round(pur, 3) if pur == pur else np.nan,
                             "frac_allcells_CNVpos": round(npos / max(nQ, 1), 3),
                             "n_MUT_pos": nmp, "n_WT_pos": nwp,
                             "OR": round(float(orr), 2), "p": float(f"{p:.1e}")})
    df = pd.DataFrame(rows)
    return df.sort_values("p") if len(df) else df


def select_mutation_specific_windows(scan_df: pd.DataFrame, min_purity: float = 0.90,
                                     max_frac_cells: float = 0.50, require_wt_lower: bool = True):
    """Filter scan_genotype_windows() output for a high-purity, sub-clonal,
    WT-clean CNV: purity_MUTfrac >= min_purity, frac_allcells_CNVpos <= max_frac_cells,
    and (optionally) pct_WT_CNVpos < pct_MUT_CNVpos."""
    if not len(scan_df):
        return scan_df
    m = (scan_df["purity_MUTfrac"] >= min_purity) & (scan_df["frac_allcells_CNVpos"] <= max_frac_cells)
    if require_wt_lower:
        m &= scan_df["pct_WT_CNVpos"] < scan_df["pct_MUT_CNVpos"]
    return scan_df[m].sort_values(["purity_MUTfrac", "pct_MUT_CNVpos"], ascending=False)


# --------------------------------------------------------------------------- #
# 3. cross-validate assays and pick the most-associated locus per gene
# --------------------------------------------------------------------------- #
def _merge_significant_windows(pdf: pd.DataFrame, sig_p: float, merge_gap: int):
    """Merge overlapping/adjacent significant windows (per chromosome) into distinct
    loci (peaks). pdf = primary-assay rows for one gene. Returns list of locus dicts
    {chr, g0, g1, peak_g0, peak_g1, r, p, direction} (one per merged region)."""
    sig = pdf[pdf["p"] < sig_p]
    loci = []
    for ch, sub in sig.groupby("chr"):
        sub = sub.sort_values("g0")
        cur = None
        for _, w in sub.iterrows():
            if cur is None:
                cur = dict(chr=ch, g0=int(w.g0), g1=int(w.g1), peak=w)
            elif int(w.g0) <= cur["g1"] + merge_gap:
                cur["g1"] = max(cur["g1"], int(w.g1))
                if w.p < cur["peak"].p:
                    cur["peak"] = w
            else:
                loci.append(cur); cur = dict(chr=ch, g0=int(w.g0), g1=int(w.g1), peak=w)
        if cur is not None:
            loci.append(cur)
    out = []
    for L in loci:
        pk = L["peak"]
        out.append({"chr": L["chr"], "g0": L["g0"], "g1": L["g1"],
                    "peak_g0": int(pk.g0), "peak_g1": int(pk.g1), "r": float(pk.r), "p": float(pk.p),
                    "direction": "down(loss)" if pk.r < 0 else "up(gain)"})
    return sorted(out, key=lambda d: d["p"])


def call_associated_loci(assoc: pd.DataFrame, params: Variant2CloneParams = Variant2CloneParams(),
                         sig_p: float = 1e-3, merge_gap: int = 10, max_loci_per_gene: int = 6):
    """Call ALL distinct CNV loci associated with each mutation (associations are
    combinatorial -- a mutation can mark several CNV positions). Selection is by the
    PRIMARY (LR) assay; significant windows are merged into peaks; the VALIDATION
    (GoT) assay, when present, confirms each locus.

    Returns (loci_df with one row per (gene, locus), best_loci dict {gene: [loci...]}).
    """
    rows, best = [], {}
    for gene, gdf in assoc.groupby("gene"):
        assays = list(gdf["assay"].unique())
        primary = _match_assay(assays, params.primary_assay) or assays[0]
        validation = next((_match_assay(assays, v) for v in params.validation_assays
                           if _match_assay(assays, v) and _match_assay(assays, v) != primary), None)
        pdf = gdf[gdf.assay == primary]
        vdf = gdf[gdf.assay == validation] if validation is not None else None
        cross = np.nan
        if vdf is not None:
            mp = pdf.set_index(["g0", "g1"])["r"]; mv = vdf.set_index(["g0", "g1"])["r"]
            common = mp.index.intersection(mv.index)
            if len(common) > 2:
                cross = float(np.corrcoef(mp.loc[common].values, mv.loc[common].values)[0, 1])
        loci = _merge_significant_windows(pdf, sig_p, merge_gap)[:max_loci_per_gene]
        best[gene] = []
        for i, L in enumerate(loci):
            validated, rv, pv = None, np.nan, np.nan
            if vdf is not None:
                vr = vdf[(vdf.g0 == L["peak_g0"]) & (vdf.g1 == L["peak_g1"])]
                if len(vr):
                    rv, pv = float(vr["r"].iloc[0]), float(vr["p"].iloc[0])
                    validated = bool((np.sign(L["r"]) == np.sign(rv)) and (pv < params.concordance_p))
            grow = pdf[(pdf.g0 == L["peak_g0"]) & (pdf.g1 == L["peak_g1"])].iloc[0]
            rows.append({"gene": gene, "locus": f"{gene}_L{i}", "chr": L["chr"],
                         "start_Mb": grow["start_Mb"], "stop_Mb": grow["stop_Mb"],
                         "size_genes": int(grow["size_genes"]), "direction": L["direction"],
                         "primary_assay": primary, "r_primary": round(L["r"], 3),
                         "p_primary": float(f"{L['p']:.1e}"), "validation_assay": validation,
                         "r_validation": (round(rv, 3) if rv == rv else np.nan),
                         "p_validation": (float(f"{pv:.1e}") if pv == pv else np.nan),
                         "validated": validated, "assay_r_concordance": (round(cross, 3) if cross == cross else np.nan),
                         "genes": grow["genes"]})
            best[gene].append({"locus": f"{gene}_L{i}", "chr": L["chr"], "g0": L["g0"], "g1": L["g1"],
                               "peak_g0": L["peak_g0"], "peak_g1": L["peak_g1"], "direction": L["direction"],
                               "r": L["r"], "validated": validated})
    loci_df = pd.DataFrame(rows).sort_values(["gene", "p_primary"]) if rows else pd.DataFrame()
    return loci_df, best


# backward-compatible alias (single best locus per gene)
def most_associated_loci(assoc, params: Variant2CloneParams = Variant2CloneParams()):
    loci_df, best = call_associated_loci(assoc, params)
    top = {g: (L[0] if L else None) for g, L in best.items()}
    return loci_df, top


# --------------------------------------------------------------------------- #
# 4. divide cells into discrete CNV clones using the validated loci as a prior
# --------------------------------------------------------------------------- #
def assign_clones(residual: pd.DataFrame, gene_order: pd.DataFrame, best_loci: dict,
                  variant_df: pd.DataFrame, cell_names, mut_specs=None,
                  params: Variant2CloneParams = Variant2CloneParams()):
    """COMBINATORIAL clones: cluster cells on the CNV scores of ALL called loci
    (across every mutation), so clones are defined by the *combination* of CNV
    positions (e.g. del5q+chr19q-loss vs del5q-only). Each clone is annotated with
    per-assay mutation frequencies and a derived ``cnv_signature``.

    best_loci : {gene: [locus dicts]} from call_associated_loci.
    Returns (clone_profile, per_cell, locus_names).
    """
    from sklearn.cluster import KMeans
    flat = [(g, L) for g, loci in best_loci.items() for L in loci]
    if not flat:
        per = pd.DataFrame({"cell": list(cell_names), "clone": "C0"})
        return pd.DataFrame([{"clone": "C0", "n_cells": len(cell_names), "cnv_signature": "(no loci)"}]), per, []
    feats, names, dirs, chrs = [], [], [], []
    for g, L in flat:
        gi = list(residual.index[L["g0"]:L["g1"] + 1])
        feats.append(residual.loc[gi].mean(axis=0).values)
        names.append(L["locus"]); dirs.append(L["direction"]); chrs.append(L["chr"])
    F = np.column_stack(feats)                                   # cells x loci (raw FC)
    sign = np.array([1.0 if d.startswith("up") else -1.0 for d in dirs])
    Fz = ((F * sign) - (F * sign).mean(0)) / ((F * sign).std(0) + 1e-9)
    k = params.n_clones or min(8, max(3, len(flat) + 1))
    labels = KMeans(n_clusters=k, random_state=params.seed, n_init=10).fit_predict(Fz)

    if mut_specs is None:
        mut_specs = infer_mutation_specs(variant_df.columns)
    vd = variant_df.reindex(list(cell_names)).fillna(0)
    overall = F.mean(0)
    prof = []
    for c in range(k):
        m = labels == c
        rec = {"clone": f"C{c}", "n_cells": int(m.sum())}
        for gene, assays in mut_specs.items():
            for asy, col in assays.items():
                rec[f"{gene}_{asy}_pct"] = round(100 * (pd.to_numeric(vd[col], errors="coerce").fillna(0).values[m] > 0).mean(), 1)
        toks = []
        for j, (g, L) in enumerate(flat):
            mu = float(F[m, j].mean()); rec[f"{L['locus']}_FC"] = round(mu, 3)
            delta = mu - overall[j]
            if L["direction"].startswith("down") and delta < -0.01:
                toks.append(f"{L['chr']}:loss")
            elif L["direction"].startswith("up") and delta > 0.01:
                toks.append(f"{L['chr']}:gain")
        rec["cnv_signature"] = ";".join(dict.fromkeys(toks)) or "(baseline)"
        prof.append(rec)
    clone_profile = pd.DataFrame(prof).sort_values("n_cells", ascending=False)
    per_cell = pd.DataFrame({"cell": list(cell_names), "clone": [f"C{c}" for c in labels]})
    return clone_profile, per_cell, names


# --------------------------------------------------------------------------- #
# orchestrator
# --------------------------------------------------------------------------- #
def run_variant2clone(residual: pd.DataFrame, gene_order: pd.DataFrame, variant_df: pd.DataFrame,
                      cell_types=None, mut_specs=None, outdir: Optional[str] = None,
                      make_plots: bool = True, params: Variant2CloneParams = Variant2CloneParams()):
    """End-to-end: intervals -> per-assay association -> cross-validated loci ->
    clones (+ optional plots). Returns a results dict; writes CSVs/plots if outdir."""
    cell_names = list(residual.columns)
    if mut_specs is None:
        mut_specs = infer_mutation_specs(variant_df.columns)
    scores, intervals = breakdown_intervals(residual, gene_order, params)
    assoc = associate_variants(scores, intervals, variant_df, cell_names, mut_specs, cell_types, params)
    loci_df, best = call_associated_loci(assoc, params)
    clone_profile, per_cell, locus_feats = assign_clones(residual, gene_order, best, variant_df,
                                                         cell_names, mut_specs, params)
    res = {"intervals": intervals, "scores": scores, "associations": assoc,
           "loci": loci_df, "best_loci": best, "clone_profile": clone_profile,
           "per_cell_clone": per_cell, "clone_features": locus_feats, "mut_specs": mut_specs}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        assoc.drop(columns=["g0", "g1"]).to_csv(os.path.join(outdir, "variant_interval_associations.csv"), index=False)
        loci_df.to_csv(os.path.join(outdir, "associated_loci.csv"), index=False)
        clone_profile.to_csv(os.path.join(outdir, "cnv_clones.csv"), index=False)
        per_cell.to_csv(os.path.join(outdir, "per_cell_clone.csv"), index=False)
        if make_plots:
            from . import plotting
            plotting.plot_all(res, residual, gene_order, variant_df, outdir, params)
    return res
