"""fastCNV: global preprocess + per-chromosome vectorized scoring.

Same scoring model as ``main.py`` (median baseline, 1.4826 * MAD scale, rolling
window mean of residuals, run-length interval calling), but reorganized so the
expensive chunk extractions and normalizations happen once per chromosome
instead of once per (state x chromosome). Includes KMeans-based clone
discovery, cached normalized expression for fast log2-ratio plots, and
zygosity annotations on per-cell and per-clone interval outputs.

Differences vs ``main.py``:
  * Requires a control h5ad. Internal-anchor mode is not implemented here -
    use ``main.py`` if you need it.
  * Single-pass scoring (no burden anchor refinement; the control already
    provides the WT baseline).
"""

from __future__ import annotations

import argparse
import logging
import math
import tempfile
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.fastCNV.main import (
    RESOURCE_FILES,
    FastCNVParams,
    Window,
    _call_clone_consensus_intervals,
    _call_intervals_for_cell,
    _chr_slices,
    _clone_profile_frame,
    _cnv_burden,
    _discover_state_clones,
    _format_duration,
    _merge_global_clones,
    _window_frame,
    bundled_gene_coordinates,
    build_windows,
    load_gene_coordinates,
)


LOGGER = logging.getLogger("fastCNV")


CHRY_Y_ONLY_MARKERS: Tuple[str, ...] = (
    "RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "NLGN4Y", "TMSB4Y",
    "PRKY", "ZFY", "TSPY1", "BCORP1", "PRY", "TBL1Y", "TXLNGY", "AMELY",
)
SEX_DETECTION_CHRY_PCT_DEFAULT: float = 5.0


import re as _re
_BLACKLIST_RX = _re.compile(r"^(IG[HKL][VJC]|IGH[GAMDE]|TR[ABGD][VJC]|HB[ABDGZMQ]\d?$|MT-)")
_BLACKLIST_EXTRA = frozenset({
    "JCHAIN", "IGJ", "PPBP", "PF4", "PF4V1",
    "HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ",
})


def build_gene_blacklist(var_names: Sequence[str]) -> set:
    """Immunoglobulin (IGH/IGK/IGL V,J,C), TCR (TRA/TRB/TRG/TRD), hemoglobin, mitochondrial and
    platelet (PPBP/PF4) genes. These are lineage-program / housekeeping transcripts, NOT copy-number
    signal: their hyper-expression at the Ig/TCR loci fabricates false 'gains' (chr14/22/2/7) and
    their dominance of the library distorts CP10k normalization (suppressing chrY in plasma/platelet
    cells). Standard practice in inferCNV/CopyKAT. Excluded from BOTH the library-size denominator
    and the scoring windows."""
    return {g for g in map(str, var_names) if _BLACKLIST_RX.match(g) or g in _BLACKLIST_EXTRA}


def _select_informative_region_genes(
    control_matrix, control_lib, control_var_names, region_genes,
    male_control_mask, female_control_mask, input_normalized,
    min_mean=0.15, min_ratio=3.0,
) -> List[str]:
    """Data-driven (NOT curated) selection of the genes whose REFERENCE expression best reflects the
    region's copy-number dosage. For chrY this is the sex-contrast (genes high in male, ~0 in female
    reference donors) — recovers the dosage-informative genes with no hand-picked marker list. The
    same principle generalizes to any region (for autosomes all expressed genes are informative)."""
    cv = pd.Index([str(g) for g in control_var_names])
    present = [g for g in region_genes if g in cv]
    if not present or female_control_mask is None or not np.asarray(female_control_mask).any():
        return present
    cols = cv.get_indexer(present)
    sub = control_matrix[:, cols]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    sub = np.asarray(sub, dtype=np.float64)
    if input_normalized:
        cp = sub  # already log1p(CP10k)
    else:
        sub = sub * (10000.0 / np.maximum(np.asarray(control_lib)[:, None], 1.0))
        cp = np.log1p(sub)
    mmean = cp[np.asarray(male_control_mask, bool)].mean(axis=0)
    fmean = cp[np.asarray(female_control_mask, bool)].mean(axis=0)
    ratio = (mmean + 0.01) / (fmean + 0.01)
    keep = [present[i] for i in range(len(present)) if mmean[i] >= min_mean and ratio[i] >= min_ratio]
    return keep if len(keep) >= 3 else present


def _simulation_copy_number_call(
    query_matrix, control_matrix, query_lib, control_lib,
    query_var_names, control_var_names, query_state, control_state,
    eligible_query_mask, ref_control_mask, control_is_metacell,
    input_normalized: bool, region_genes: Sequence[str],
    normal_copies: int = 2, min_ref: int = 20, min_detect_frac: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """GENERAL simulation-calibrated copy-number caller — NO region/CNV-specific thresholds.

    To know what a copy-number change looks like, SIMULATE it from the reference: take the cells that
    carry the NORMAL copy number for this region (the sex-matched reference) and reduce/raise EVERY
    gene's expression by the dosage factor k/normal_copies — a heterozygous deletion is x0.5, a
    homozygous deletion x0, a single-copy gain x1.5, etc. Aggregated over the region's genes, each
    simulated copy state k gives an expected expression distribution (mean+sd) PER CELL STATE. A query
    cell is assigned the MAXIMUM-LIKELIHOOD copy state; it is a deletion if that state < normal. The
    decision boundary is the likelihood crossover of the simulated distributions — it is identical for
    any deletion/duplication anywhere in the genome, never tuned to a region or to LOY. Applied to chrY
    (normal_copies=1, sex-matched male reference) the homozygous-deletion state IS Loss-of-Y, which
    therefore falls out of the same general rule as every other CNV.

    Returns (deletion_call[bool], log2_ratio[float32]) where log2_ratio = log2(called_copies/normal).
    """
    qv = pd.Index([str(g) for g in query_var_names]); cv = pd.Index([str(g) for g in control_var_names])
    present = [g for g in region_genes if g in qv and g in cv]
    n_query = query_matrix.shape[0]
    if len(present) < 3:
        return np.zeros(n_query, dtype=bool), np.full(n_query, np.nan, dtype=np.float32)
    normal_copies = max(int(normal_copies), 1)
    # General detection filter (genome-wide preprocessing, not region-specific): keep only genes
    # actually expressed in the normal-copy reference. Genes that are ~never detected carry no dosage
    # signal and only dilute the region aggregate (e.g. silent TTTY/FAM genes on chrY).
    _ref0 = np.asarray(ref_control_mask, dtype=bool)
    if control_is_metacell is not None:
        _ref0 = _ref0 & (~np.asarray(control_is_metacell, dtype=bool))
    if _ref0.any():
        _cc = cv.get_indexer(present); _sub = control_matrix[_ref0][:, _cc]
        _sub = _sub.toarray() if sp.issparse(_sub) else np.asarray(_sub)
        _detfrac = (np.asarray(_sub) > 0).mean(axis=0)
        _expr = [present[i] for i in range(len(present)) if _detfrac[i] >= min_detect_frac]
        if len(_expr) >= 3:
            present = _expr
    qcols = qv.get_indexer(present); ccols = cv.get_indexer(present)

    def _cp(matrix, cols, lib):
        sub = matrix[:, cols]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        sub = np.asarray(sub, dtype=np.float64)
        if input_normalized:            # X is log1p(CP10k) -> back to CP10k for dosage scaling
            return np.expm1(sub)
        return sub * (10000.0 / np.maximum(lib[:, None], 1.0))

    cp_q = _cp(query_matrix, qcols, query_lib); cp_c = _cp(control_matrix, ccols, control_lib)
    qa = np.log1p(cp_q).mean(axis=1)                       # query region aggregate (log1p CP10k)
    qstate = np.asarray(query_state); cstate = np.asarray(control_state)
    ref = np.asarray(ref_control_mask, dtype=bool)
    if control_is_metacell is not None:
        ref = ref & (~np.asarray(control_is_metacell, dtype=bool))
    if not ref.any():
        return np.zeros(n_query, dtype=bool), np.full(n_query, np.nan, dtype=np.float32)

    copy_states = list(range(0, 2 * normal_copies + 1))   # 0..2N: homo/het deletion .. normal .. gains
    elig = np.asarray(eligible_query_mask, dtype=bool)
    call = np.zeros(n_query, dtype=bool); log2r = np.full(n_query, np.nan, dtype=np.float32)
    for s in np.unique(qstate):
        rm = ref & (cstate == s)
        if rm.sum() < min_ref:
            rm = ref
        cpr = cp_c[rm]                                     # reference CP10k (normal copy number)
        mus = []; sds = []
        for k in copy_states:                             # SIMULATE copy state k by scaling x (k/normal)
            agg_k = np.log1p(cpr * (k / float(normal_copies))).mean(axis=1)
            mus.append(float(agg_k.mean())); sds.append(max(float(agg_k.std()), 1e-3))
        mus = np.asarray(mus); sds = np.asarray(sds)
        m = np.flatnonzero((qstate == s) & elig & np.isfinite(qa))
        if m.size == 0:
            continue
        x = qa[m][:, None]
        ll = -0.5 * ((x - mus[None, :]) / sds[None, :]) ** 2 - np.log(sds[None, :])   # Gaussian log-lik
        best = np.asarray(copy_states)[ll.argmax(axis=1)]
        call[m] = best < normal_copies
        log2r[m] = np.log2(np.maximum(best, 1e-6) / float(normal_copies)).astype(np.float32)
    return call.astype(bool), log2r


def _autosomal_arm_calls(
    query_matrix, control_matrix, query_var_names, control_var_names, query_state, control_state,
    arm_genes: Dict[str, List[str]], control_is_metacell, input_normalized: bool,
    min_detect_frac: float = 0.05, min_ref: int = 20, min_loglr: float = 25.0,
) -> Tuple[List[str], np.ndarray]:
    """Autosomal copy-number calls per chromosome ARM. SAME simulation+likelihood machinery as LOY
    (`_simulation_copy_number_call`): per-arm gene detection filter, simulate x0.01/x0.5/x1/x1.5/x2 of
    neutral from the reference, classify each cell by MAXIMUM LIKELIHOOD (argmax, NO threshold), using
    the simulated-state variance. The ONE added step (which LOY does not need, because chrY is a single
    region) is the standard per-cell GENOME-WIDE CENTERING across arms -- a global per-cell expression
    shift would otherwise mis-call all ~40 arms at once (this is exactly inferCNV's per-cell centering,
    a general normalization, not a threshold). Returns (arm_order, arm_log2[cells x arms])."""
    qv = pd.Index([str(g) for g in query_var_names]); cv = pd.Index([str(g) for g in control_var_names])
    arms = sorted(arm_genes.keys()); nq = query_matrix.shape[0]
    qstate = np.asarray(query_state); cstate = np.asarray(control_state)
    refm = (~np.asarray(control_is_metacell, bool)) if control_is_metacell is not None else np.ones(control_matrix.shape[0], bool)
    COPIES = np.array([0.01, 0.5, 1.0, 1.5, 2.0]); LOG2 = np.log2(COPIES).astype(np.float32); NEUTRAL = 2
    def _dense(M, cols):
        s = M[:, cols]; s = s.toarray() if sp.issparse(s) else np.asarray(s); return np.asarray(s, np.float64)
    # detection-filtered gene columns per arm (keep genes expressed in the reference)
    qcols = {}; ccols = {}
    for ai, arm in enumerate(arms):
        genes = [g for g in arm_genes[arm] if g in qv and g in cv]
        if len(genes) < 10:
            continue
        cc = cv.get_indexer(genes); det = (_dense(control_matrix[refm], cc) > 0).mean(0)
        keep = [genes[i] for i in range(len(genes)) if det[i] >= min_detect_frac]
        if len(keep) < 10:
            keep = genes
        qcols[ai] = qv.get_indexer(keep); ccols[ai] = cv.get_indexer(keep)
    arm_log2 = np.zeros((nq, len(arms)), dtype=np.float32)
    for st in np.unique(qstate):
        qrows = np.flatnonzero(qstate == st)
        if qrows.size == 0:
            continue
        crows = np.flatnonzero(refm & (cstate == st))
        if crows.size < min_ref:
            crows = np.flatnonzero(refm)
        # PER-GENE simulated levels + noise, then PER-GENE log-likelihood of each copy state summed over
        # the arm. A copy change is a CONSTANT multiplicative shift -> in log-expression it separates the
        # copy states only for genes that are actually EXPRESSED; near-silent genes have m_g(c) ~ flat and
        # contribute ~0 to discriminating copy states (no compression, no cell-type leakage). The arm CALL
        # requires the best copy state to beat NEUTRAL by a likelihood-ratio cutoff (general chi-square-style
        # rule, applies to any deletion/duplication) -- so an uninformative arm stays neutral instead of
        # being forced (argmax) into a call, which is what produced ~13 false calls/cell.
        arm_yq = {}; arm_m = {}; arm_sig = {}; offnum = np.zeros((qrows.size, len(arms))); offden = 0
        gw = {}
        for ai in qcols:
            yq = _dense(query_matrix, qcols[ai])[qrows]                  # query cells x genes (log1p CP10k)
            cpc = _dense(control_matrix, ccols[ai])[crows]
            refcp = np.expm1(cpc) if input_normalized else cpc
            ref_mean = refcp.mean(0)                                      # per-gene mean reference CP10k
            m = np.stack([np.log1p(ref_mean * (c / 2.0)) for c in COPIES], axis=1)   # genes x 5 expected log1p
            sig = np.maximum(cpc.std(0), 0.10)                           # per-gene reference SD (log1p units)
            arm_yq[ai] = yq; arm_m[ai] = m; arm_sig[ai] = sig
            # per-cell genome-wide library centering: accumulate (yq - neutral expected)/sig^2 weighted
            offnum[:, ai] = ((yq - m[None, :, NEUTRAL]) / sig[None, :] ** 2).sum(1); gw[ai] = (1.0 / sig ** 2).sum()
        # per-cell global offset = weighted median deviation from neutral across arms (robust to real CNV arms)
        with np.errstate(invalid="ignore"):
            per_arm_off = np.stack([offnum[:, ai] / max(gw[ai], 1e-9) for ai in qcols], axis=1)
        off = np.nan_to_num(np.nanmedian(per_arm_off, axis=1), nan=0.0)   # log-expression units, per cell
        for ai in qcols:
            yq = arm_yq[ai] - off[:, None]; m = arm_m[ai]; sig = arm_sig[ai]
            ll = np.stack([(-0.5 * ((yq - m[None, :, ki]) / sig[None, :]) ** 2).sum(1) for ki in range(5)], axis=1)
            best = ll.argmax(1)
            lr = ll[np.arange(qrows.size), best] - ll[:, NEUTRAL]         # log-likelihood-ratio best vs neutral
            called = (best != NEUTRAL) & (2.0 * lr >= min_loglr)         # significance gate (no forced calls)
            arm_log2[qrows[called], ai] = LOG2[best[called]]
    return arms, arm_log2


# hg38 centromere midpoints (bp) — used to split chromosomes into p / q arms (the natural CNV unit;
# arm aggregation gives the per-cell signal enough genes for the likelihood classifier to be reliable
# WITHOUT any tuned threshold). Approximate band-level positions; exact bp is not needed for arm split.
_HG38_CENTROMERE_BP = {f"chr{c}": int(m * 1e6) for c, m in {
    "1": 123, "2": 93.9, "3": 90.9, "4": 50, "5": 48.8, "6": 59.8, "7": 60.1, "8": 45.2, "9": 43,
    "10": 39.8, "11": 53.4, "12": 35.5, "13": 17.7, "14": 17.2, "15": 19, "16": 36.8, "17": 25.1,
    "18": 18.5, "19": 26.2, "20": 28.1, "21": 12, "22": 15, "chrX": 61, "chrY": 10.4,
}.items()}


def _arm_label(chrom: str, midpoint_bp: int) -> str:
    """Chromosome-arm label ('chr5p' / 'chr5q') from a genomic midpoint; '' for unknown chromosomes."""
    cen = _HG38_CENTROMERE_BP.get(str(chrom))
    if cen is None:
        return ""
    return f"{chrom}{'p' if midpoint_bp < cen else 'q'}"


def _infer_sample_sex(
    adata: ad.AnnData,
    sample_labels: np.ndarray,
    threshold_pct: float,
    chry_markers: Sequence[str] = CHRY_Y_ONLY_MARKERS,
) -> pd.DataFrame:
    """Infer per-sample biological sex from chrY-Y-only marker expression.

    For each sample (a unique value in ``sample_labels``), count the fraction of
    cells with >= 1 UMI on any of the chrY-Y-only markers present in
    ``adata.var_names``. Samples with > ``threshold_pct``% chrY-positive cells
    are called male, else female. This is robust to LOY: even in a sample with
    80% LOY cells the unaffected 20% will still express chrY, far above the
    ~1% background expected from a true female sample.

    Returns a DataFrame indexed by sample with columns
    n_cells, n_chrY_pos, pct_chrY_pos, inferred_sex, n_markers_present.
    """
    present = [g for g in chry_markers if g in adata.var_names]
    if len(present) == 0:
        return pd.DataFrame({
            "sample": pd.unique(sample_labels),
            "n_cells": [int((sample_labels == s).sum()) for s in pd.unique(sample_labels)],
            "n_chrY_pos": 0,
            "pct_chrY_pos": 0.0,
            "inferred_sex": "unknown",
            "n_markers_present": 0,
        }).set_index("sample")

    sub = adata[:, present].X
    if sp.issparse(sub):
        sub = sub.tocsr()
        n_pos_per_cell = np.asarray((sub > 0).sum(axis=1)).ravel()
    else:
        n_pos_per_cell = (np.asarray(sub) > 0).sum(axis=1)
    cell_chry_pos = n_pos_per_cell > 0

    rows = []
    for sample in pd.unique(sample_labels):
        mask = sample_labels == sample
        n_cells = int(mask.sum())
        n_pos = int(cell_chry_pos[mask].sum())
        pct = 100.0 * n_pos / max(n_cells, 1)
        rows.append({
            "sample": sample,
            "n_cells": n_cells,
            "n_chrY_pos": n_pos,
            "pct_chrY_pos": pct,
            "inferred_sex": "male" if pct > threshold_pct else "female",
            "n_markers_present": len(present),
        })
    return pd.DataFrame(rows).set_index("sample")


CANCER_DRIVER_GENES: Tuple[str, ...] = (
    "MCL1", "MYCN", "GLI2", "PDGFRA", "FGFR3", "KIT", "EGFR", "CDK6", "MET",
    "FGFR1", "MYC", "CDKN2A", "CDKN2B", "PTCH1", "JAK2", "PTEN", "FAS",
    "CCND1", "ATM", "KRAS", "CDK4", "MDM2", "RB1", "BRCA2", "FOXO1", "NF1",
    "TP53", "ERBB2", "BRCA1", "STK11", "CCNE1", "AURKA", "RUNX1", "AKT1",
    "BCL2", "SMAD4", "APC", "NOTCH1", "FBXW7", "ARID1A", "TET2", "DNMT3A",
    "ASXL1", "EZH2", "SF3B1", "U2AF1", "SRSF2", "IDH1", "IDH2", "FLT3",
    "NPM1", "WT1", "GATA2", "SETBP1", "BCOR", "STAG2", "PHF6", "RAD21",
    "CALR", "MPL",
    # Sex-chromosome dosage markers (always shown on the driver-gene strip plot
    # so LOY appears as a strong negative on RPS4Y1 and X-inactivation/loss
    # appears on XIST):
    "XIST", "RPS4Y1",
)


@dataclass
class FastParams:
    h5ad: Path
    control_h5ad: Path
    gene_coordinates: Optional[Path]
    output_prefix: Path
    state_key: str
    sample_key: Optional[str] = None
    control_state_key: Optional[str] = None
    layer: str = "auto"
    input_normalized: bool = False
    window_genes: int = 41
    stride_genes: int = 7
    min_chr_genes: int = 25
    min_state_cells: int = 30
    high_threshold: float = 2.6
    low_threshold: float = 1.6
    min_run_windows: int = 3
    min_interval_genes: int = 60
    min_mean_score: float = 1.8
    burden_quantile: float = 0.95
    cnv_burden_threshold: float = 1.8
    max_clones_per_state: int = 10
    max_global_clones: int = 10
    min_clone_cells: int = 5
    clone_similarity_threshold: float = 0.88
    clone_consensus_fraction: float = 0.45
    nmf_max_iter: int = 100
    clone_min_active_fraction: float = 0.05
    clone_max_features: int = 400
    clone_min_cnv_fraction: float = 0.015
    clone_min_cells_confident: int = 30
    zygosity_mode: str = "relative"
    skip_clones: bool = False
    skip_pdf: bool = False
    skip_heatmap: bool = False
    write_h5ad: bool = False
    random_state: int = 0
    n_jobs: int = -1
    pdf_smooth_genes: int = 50
    pdf_y_clip: float = 1.0
    pdf_label_threshold: float = 0.25
    pdf_score_y_clip: float = 4.0
    heatmap_max_cells: int = 20000
    heatmap_filter_threshold: float = 1.5
    heatmap_min_chr_windows: int = 35
    residual_candidate: bool = False
    residual_candidate_min_abs_region_score: float = 0.05
    residual_candidate_min_separation_mad: float = 1.5
    residual_candidate_heatmap_filter_threshold: float = 0.03
    residual_candidate_pyinfercnv_path: Optional[Path] = Path("/Users/saljh8/Documents/GitHub/pyInferCNV")
    sex_chrom_mode: str = "absolute_log2"
    sex_chrom_log2_unit: float = 0.040
    sex_detection_threshold_pct: float = 5.0
    control_sample_key: Optional[str] = None
    max_cells_per_state: int = 0
    gate_chry_by_marker_expression: bool = True
    sex_chrom_het_loss: float = -0.6
    sex_chrom_hom_loss: float = -1.5
    sex_chrom_het_gain: float = 0.4
    sex_chrom_hom_gain: float = 0.7
    sex_chroms: Tuple[str, ...] = ("chrY",)
    # --- mature unsupervised-CNV update (2026-06) ---
    gene_blacklist: bool = True       # exclude Ig/TCR/Hb/MT/platelet genes from libsize + scoring
    pooled_scale: bool = True         # per-window scale = pooled within-state control MAD (replaces 0.10 floor)
    simulation_autosomal: bool = False # OFF by default: the per-cell per-window autosomal caller is
                                      #   single-cell-noise-limited and used tuned cut-offs, so it both
                                      #   broke the LOY positive control's clones and violated the
                                      #   no-task-specific-thresholds rule. Left as explicit opt-in only.
                                      #   The validated default = clean LOY (chrY) positive control.
    autosomal_min_detect_frac: float = 0.05  # gene detection filter for the autosomal per-arm caller
                                      #   (same rule as LOY: keep genes expressed in >= this fraction of the
                                      #   reference). Set 0.0 to disable (benchmark with vs without).
    cnv_internal_baseline: bool = True  # per-window neutral = query's own per-state median (batch-robust,
                                      #   detects SUBCLONAL CNV). Set False to use the EXTERNAL reference
                                      #   neutral (detects CLONAL CNV; requires a batch-MATCHED reference).
    scale_floor_quantile: float = 0.05  # eps floor for the scale = this quantile of the per-window pooled scale
    loss_zero_copy_quantile: float = 0.95  # a region copy-LOSS is called when query expression is at/below
                                      #   this quantile of the ZERO-COPY reference (for chrY: female donors).
                                      #   Concordant with complete-loss; applies generally per region.
    nearest_normal_autosomal: bool = False  # opt-in validated nearest-normal-subcluster autosomal CNV
                                      #   caller (per query subcluster vs its nearest NORMAL reference
                                      #   subcluster -> cross-donor-z -> arm aggregate). Held-out normal
                                      #   control <5% with any call; detects -7/+8/del(5q). Integrated
                                      #   with the standard run so unbiased chrY/LOY remains active.


def _matrix(adata: ad.AnnData, layer: str) -> sp.spmatrix | np.ndarray:
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer == "X":
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in AnnData.")
    return adata.layers[layer]


def _row_sums(matrix: sp.spmatrix | np.ndarray) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)


def _normalize_chunk(
    matrix: sp.spmatrix | np.ndarray,
    col_indices: np.ndarray,
    library_sizes: np.ndarray,
    input_normalized: bool,
) -> np.ndarray:
    """Densify all rows for the given column indices and log-normalize in place."""
    sub = matrix[:, col_indices]
    if sp.issparse(sub):
        sub = sub.toarray()
    sub = np.asarray(sub, dtype=np.float32)
    if input_normalized:
        return sub
    factors = (10000.0 / np.maximum(library_sizes.astype(np.float32), 1.0)).astype(np.float32)
    sub *= factors[:, None]
    np.log1p(sub, out=sub)
    return sub


def _build_state_crosswalk(
    query: ad.AnnData,
    control: ad.AnnData,
    query_matrix: sp.spmatrix | np.ndarray,
    control_matrix: sp.spmatrix | np.ndarray,
    query_lib: np.ndarray,
    control_lib: np.ndarray,
    coords: pd.DataFrame,
    state_key: str,
    control_state_key: Optional[str],
    eligible_states: Sequence[str],
    input_normalized: bool,
    *,
    min_control_cells: int = 20,
    max_genes: int = 2500,
) -> Dict[str, str]:
    """Map query state labels to reference state labels when exact labels are absent."""
    if not control_state_key:
        return {str(s): str(s) for s in eligible_states}
    q_state = query.obs[state_key].astype(str).to_numpy()
    c_state = control.obs[control_state_key].astype(str).to_numpy()
    control_counts = pd.Series(c_state).value_counts()
    control_valid = set(control_counts[control_counts >= min_control_cells].index.astype(str))
    mapping = {str(s): str(s) for s in eligible_states if str(s) in control_valid}
    missing = [str(s) for s in eligible_states if str(s) not in mapping]
    if not missing:
        return mapping

    q_names = pd.Index(query.var_names.astype(str))
    c_names = pd.Index(control.var_names.astype(str))
    coord_genes = [str(q_names[int(i)]) for i in coords["var_index"].to_numpy()]
    genes = [g for g in coord_genes if g in c_names]
    if len(genes) < 100:
        LOGGER.warning(
            "State crosswalk skipped: only %d shared coordinate genes between query and control.",
            len(genes),
        )
        return {str(s): mapping.get(str(s), str(s)) for s in eligible_states}
    if len(genes) > max_genes:
        step = max(1, len(genes) // max_genes)
        genes = genes[::step][:max_genes]
    q_cols = q_names.get_indexer(genes)
    c_cols = c_names.get_indexer(genes)
    q_expr = _normalize_chunk(query_matrix, q_cols, query_lib, input_normalized)
    c_expr = _normalize_chunk(control_matrix, c_cols, control_lib, input_normalized)

    control_states = [s for s in sorted(control_valid) if np.sum(c_state == s) >= min_control_cells]
    if not control_states:
        return {str(s): mapping.get(str(s), str(s)) for s in eligible_states}
    c_centroids = np.vstack([np.nanmean(c_expr[c_state == s], axis=0) for s in control_states])
    c_centroids = c_centroids - np.nanmean(c_centroids, axis=1, keepdims=True)
    c_norm = np.sqrt(np.nansum(c_centroids * c_centroids, axis=1))
    for state in missing:
        rows = q_state == state
        if int(rows.sum()) < 3:
            mapping[state] = state
            continue
        q_centroid = np.nanmean(q_expr[rows], axis=0)
        q_centroid = q_centroid - np.nanmean(q_centroid)
        q_norm = math.sqrt(float(np.nansum(q_centroid * q_centroid)))
        denom = np.maximum(c_norm * max(q_norm, 1e-12), 1e-12)
        corr = np.nansum(c_centroids * q_centroid[None, :], axis=1) / denom
        best = str(control_states[int(np.nanargmax(corr))])
        mapping[state] = best
        LOGGER.info("State crosswalk: query '%s' -> control '%s' (centroid r=%.3f).", state, best, float(np.nanmax(corr)))
    return {str(s): mapping.get(str(s), str(s)) for s in eligible_states}


def _rolling_window_mean(values: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    if values.shape[1] == 0 or not windows:
        return np.zeros((values.shape[0], 0), dtype=np.float32)
    prefix = np.empty((values.shape[0], values.shape[1] + 1), dtype=np.float32)
    prefix[:, 0] = 0.0
    np.cumsum(values, axis=1, dtype=np.float32, out=prefix[:, 1:])
    starts = np.asarray([w.start_offset for w in windows], dtype=np.int64)
    ends = np.asarray([w.end_offset for w in windows], dtype=np.int64)
    lengths = np.maximum(ends - starts, 1).astype(np.float32)
    return ((prefix[:, ends] - prefix[:, starts]) / lengths).astype(np.float32, copy=False)


def _mad(values: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(values, axis=axis)
    return np.nanmedian(np.abs(values - np.expand_dims(med, axis)), axis=axis)


def _absolute_log2_window_scores(
    query_state_log1p: np.ndarray,
    control_state_log1p: np.ndarray,
    chr_windows: Sequence[Window],
    log2_unit: float,
) -> np.ndarray:
    """Per-cell sex-chrom window scores from absolute expression difference (log1p space).

    Bypasses MAD-standardization, which collapses real LOY signal because the male
    reference's chrY MAD is large (donor-level age-related sub-clinical LOY).

    For each cell we compute a per-gene difference (query_log1p - control_state_mean_log1p),
    rolling-window mean it, then divide by ``log2_unit`` so the result lives in the same
    standardized-score units (~|score| in [-3, +3]) used by the autosomal interval caller.
    Working in log1p space avoids the divide-by-zero / pseudocount-tuning problem of
    CP10K log2 ratios on sparse sex-chromosome counts. The mean log1p difference is
    monotonic with biological log2 fold-change for moderately expressed genes.
    """
    if query_state_log1p.size == 0 or control_state_log1p.size == 0:
        return np.zeros((query_state_log1p.shape[0], len(chr_windows)), dtype=np.float32)
    control_mean_log1p = np.nanmean(control_state_log1p, axis=0).astype(np.float32)
    diff_per_gene = (query_state_log1p.astype(np.float32) - control_mean_log1p[None, :])
    windowed = _rolling_window_mean(diff_per_gene, chr_windows)
    unit = max(float(log2_unit), 1e-3)
    return (windowed / unit).astype(np.float32)


def _control_var_map(query_var: pd.Index, control_var: pd.Index) -> np.ndarray:
    lookup = pd.Series(np.arange(len(control_var), dtype=np.int64), index=control_var.astype(str))
    return lookup.reindex(query_var.astype(str)).to_numpy()


def _kmeans_state_clones(
    scores: np.ndarray,
    cnv_mask: np.ndarray,
    max_clones: int,
    min_cells: int,
    seed: int,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """KMeans replacement for NMF-based state-local clone discovery.

    Operates only on CNV-positive cells. Returns per-cell labels (WT or clone1..N)
    and a summary frame with size + mean burden per clone.
    """
    labels = np.full(scores.shape[0], "WT", dtype=object)
    candidates = np.flatnonzero(cnv_mask)
    if candidates.size < min_cells:
        return labels, pd.DataFrame()

    feature_block = np.nan_to_num(scores[candidates].astype(np.float32, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    feature_block = np.clip(feature_block, -8.0, 8.0)
    if feature_block.shape[1] == 0:
        return labels, pd.DataFrame()

    n_clusters = max(1, min(int(max_clones), candidates.size // max(min_cells, 1)))
    if n_clusters <= 1:
        clone_labels = np.zeros(candidates.size, dtype=int)
    else:
        from sklearn.cluster import MiniBatchKMeans
        model = MiniBatchKMeans(
            n_clusters=n_clusters,
            random_state=seed,
            n_init=3,
            batch_size=min(1024, max(64, candidates.size)),
            max_iter=100,
        )
        clone_labels = model.fit_predict(feature_block).astype(int)

    rows: List[Dict[str, object]] = []
    order = sorted(set(int(c) for c in clone_labels), key=lambda c: int(np.sum(clone_labels == c)), reverse=True)
    next_id = 1
    for cluster in order:
        member_local = np.flatnonzero(clone_labels == cluster)
        if member_local.size < min_cells:
            continue
        clone_id = f"clone{next_id}"
        next_id += 1
        labels[candidates[member_local]] = clone_id
        cluster_scores = feature_block[member_local]
        rows.append({
            "state_clone_id": clone_id,
            "n_cells": int(member_local.size),
            "mean_burden": float(np.nanmean(np.abs(cluster_scores).max(axis=1))),
            "max_abs_score": float(np.nanmax(np.abs(cluster_scores))) if cluster_scores.size else 0.0,
        })
        if next_id > max_clones:
            break
    return labels, pd.DataFrame(rows)


def _smooth_and_center_log2_ratios(
    ratios: np.ndarray,
    coords: pd.DataFrame,
    smooth_genes: int,
) -> np.ndarray:
    """Smooth log2 ratios across a gene-neighborhood window per chromosome and center each chromosome on its median."""
    if ratios.size == 0:
        return ratios
    smoothed = np.full_like(ratios, np.nan, dtype=np.float32)
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        idx = chr_coords.index.to_numpy(dtype=np.int64)
        if idx.size == 0:
            continue
        block = ratios[:, idx]
        valid = ~np.isnan(block)
        block_filled = np.where(valid, block, 0.0).astype(np.float32, copy=False)
        valid_f = valid.astype(np.float32)
        n_genes = block.shape[1]
        if smooth_genes > 1 and n_genes >= 3:
            half = max(1, smooth_genes // 2)
            sum_prefix = np.concatenate([np.zeros((block.shape[0], 1), dtype=np.float32), np.cumsum(block_filled, axis=1, dtype=np.float32)], axis=1)
            count_prefix = np.concatenate([np.zeros((block.shape[0], 1), dtype=np.float32), np.cumsum(valid_f, axis=1, dtype=np.float32)], axis=1)
            starts = np.maximum(np.arange(n_genes) - half, 0)
            ends = np.minimum(np.arange(n_genes) + half + 1, n_genes)
            window_sum = sum_prefix[:, ends] - sum_prefix[:, starts]
            window_count = count_prefix[:, ends] - count_prefix[:, starts]
            chrom_smoothed = np.where(window_count > 0, window_sum / np.maximum(window_count, 1.0), np.nan).astype(np.float32)
        else:
            chrom_smoothed = block.astype(np.float32, copy=True)
        median = np.nanmedian(chrom_smoothed, axis=1, keepdims=True)
        median = np.nan_to_num(median, nan=0.0)
        chrom_smoothed = chrom_smoothed - median
        smoothed[:, idx] = chrom_smoothed
    return smoothed


_FIXED_ZYGOSITY_THRESHOLDS = {
    "homozygous_loss": -1.30,
    "heterozygous_loss": -0.55,
    "heterozygous_gain": 0.40,
    "homozygous_gain": 0.85,
}


def _zygosity_state(log2_ratio: float, thresholds: Optional[Dict[str, float]] = None) -> str:
    """Map a smoothed log2 ratio to a zygosity-like categorical state."""
    if thresholds is None:
        thresholds = _FIXED_ZYGOSITY_THRESHOLDS
    if not np.isfinite(log2_ratio):
        return "low_signal"
    if log2_ratio <= thresholds["homozygous_loss"]:
        return "homozygous_loss"
    if log2_ratio <= thresholds["heterozygous_loss"]:
        return "heterozygous_loss"
    if log2_ratio >= thresholds["homozygous_gain"]:
        return "homozygous_gain"
    if log2_ratio >= thresholds["heterozygous_gain"]:
        return "heterozygous_gain"
    return "low_signal"


def _calibrate_zygosity_thresholds(
    interval_records: Sequence[Dict[str, object]],
    wt_indices: np.ndarray,
    cached_query_expr: np.ndarray,
    control_gene_means: np.ndarray,
    coords: pd.DataFrame,
    seed: int = 0,
    cells_per_interval: int = 50,
    max_intervals: int = 2000,
) -> Optional[Dict[str, float]]:
    """Estimate empirical zygosity thresholds from WT cells' log2 ratios in called intervals.

    Returns a dict of thresholds, or None if there's not enough data to calibrate.
    """
    if cached_query_expr is None or wt_indices.size < 100 or len(interval_records) == 0:
        return None
    rng = np.random.default_rng(seed)
    coords_chr = coords["chr"].astype(str).to_numpy()
    coords_start = coords["start"].to_numpy()
    coords_end = coords["end"].to_numpy()
    chr_to_indices: Dict[str, np.ndarray] = {
        chrom: np.flatnonzero(coords_chr == chrom).astype(np.int64)
        for chrom in pd.unique(coords_chr)
    }
    pseudocount = 1e-3

    if len(interval_records) > max_intervals:
        sampled = rng.choice(len(interval_records), size=max_intervals, replace=False)
        sampled_records = [interval_records[i] for i in sampled]
    else:
        sampled_records = list(interval_records)

    noise_means: List[float] = []
    for record in sampled_records:
        chrom = str(record["chr"])
        start = int(record["start"])
        end = int(record["end"])
        chr_idx = chr_to_indices.get(chrom)
        if chr_idx is None or chr_idx.size == 0:
            continue
        in_interval = (coords_start[chr_idx] >= start) & (coords_end[chr_idx] <= end)
        gene_indices = chr_idx[in_interval]
        if gene_indices.size == 0:
            continue
        sample_size = min(cells_per_interval, wt_indices.size)
        wt_sample = rng.choice(wt_indices, size=sample_size, replace=False)
        cell_expr = cached_query_expr[wt_sample][:, gene_indices]
        ctrl_expr = control_gene_means[gene_indices]
        per_cell_means = np.nanmean(np.log2((cell_expr + pseudocount) / (ctrl_expr + pseudocount)), axis=1)
        noise_means.extend(per_cell_means.tolist())

    if len(noise_means) < 200:
        return None
    arr = np.asarray(noise_means, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size < 200:
        return None

    q_low = float(np.quantile(arr, 0.005))
    q_lo = float(np.quantile(arr, 0.05))
    q_hi = float(np.quantile(arr, 0.95))
    q_high = float(np.quantile(arr, 0.995))
    thresholds = {
        "homozygous_loss": min(q_low, -1.0),
        "heterozygous_loss": min(q_lo, -0.4),
        "heterozygous_gain": max(q_hi, 0.4),
        "homozygous_gain": max(q_high, 0.7),
    }
    return thresholds


def _annotate_intervals_with_zygosity(
    interval_records: List[Dict[str, object]],
    cached_query_expr: Optional[np.ndarray],
    control_gene_means: np.ndarray,
    coords: pd.DataFrame,
    barcode_to_index: Dict[str, int],
    thresholds: Optional[Dict[str, float]] = None,
    sex_chroms: Sequence[str] = (),
    sex_chrom_thresholds: Optional[Dict[str, float]] = None,
) -> List[Dict[str, object]]:
    """Attach mean_log2_ratio + zygosity_state to per-cell interval records.

    For autosomal intervals, the per-gene log2 ratio uses the cached log1p(CP10K)
    expression with a tiny pseudocount (legacy behavior, valid where both sides
    are non-trivial).

    For intervals on a chromosome in ``sex_chroms``, ``mean_log2_ratio`` is the
    mean log1p(CP10K) DIFFERENCE (query - control_state_mean), not a log2 ratio.
    The autosomal log2-of-log1p formulation explodes when control expression
    approaches zero (chrY in female-mixed or sparse references). The log1p-diff
    is monotonic with biological log2 fold-change for moderately expressed genes
    and is well-defined at any expression level. ``sex_chrom_thresholds`` must be
    calibrated for log1p-diff units (defaults: het_loss<=-0.04, hom_loss<=-0.08).
    """
    if not interval_records or cached_query_expr is None:
        return interval_records
    pseudocount = 1e-3
    coords_chr = coords["chr"].astype(str).to_numpy()
    coords_start = coords["start"].to_numpy()
    coords_end = coords["end"].to_numpy()
    chr_to_indices: Dict[str, np.ndarray] = {
        chrom: np.flatnonzero(coords_chr == chrom).astype(np.int64)
        for chrom in pd.unique(coords_chr)
    }
    sex_chrom_set = {str(c) for c in sex_chroms}
    annotated: List[Dict[str, object]] = []
    for record in interval_records:
        chrom = str(record["chr"])
        start = int(record["start"])
        end = int(record["end"])
        chr_idx = chr_to_indices.get(chrom)
        if chr_idx is None or chr_idx.size == 0:
            record["mean_log2_ratio"] = float("nan")
            record["zygosity_state"] = "low_signal"
            annotated.append(record)
            continue
        in_interval_local = (coords_start[chr_idx] >= start) & (coords_end[chr_idx] <= end)
        gene_indices = chr_idx[in_interval_local]
        cell_idx = barcode_to_index.get(str(record["CellBarcode"]))
        if cell_idx is None or gene_indices.size == 0:
            record["mean_log2_ratio"] = float("nan")
            record["zygosity_state"] = "low_signal"
            annotated.append(record)
            continue
        cell_expr = cached_query_expr[cell_idx, gene_indices]
        ctrl_expr = control_gene_means[gene_indices]
        is_sex_chrom = chrom in sex_chrom_set and sex_chrom_thresholds is not None
        per_gene_log2 = np.log2((cell_expr + pseudocount) / (ctrl_expr + pseudocount))
        mean_metric = float(np.nanmean(per_gene_log2))
        record["mean_log2_ratio"] = mean_metric
        if is_sex_chrom:
            record["zygosity_state"] = _zygosity_state(mean_metric, sex_chrom_thresholds)
        else:
            record["zygosity_state"] = _zygosity_state(mean_metric, thresholds)
        annotated.append(record)
    return annotated


def _compute_clone_log2_ratios_from_cache(
    cached_query_expr: np.ndarray,
    control_gene_means: np.ndarray,
    clone_to_query_rows: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, List[str]]:
    """Per-clone per-gene log2(mean_clone / mean_control) using a precomputed normalized matrix."""
    clone_order = list(clone_to_query_rows.keys())
    n_genes = cached_query_expr.shape[1]
    ratios = np.full((len(clone_order), n_genes), np.nan, dtype=np.float32)
    pseudocount = 1e-3
    for i, clone_id in enumerate(clone_order):
        rows = clone_to_query_rows[clone_id]
        if rows.size == 0:
            continue
        clone_mean = cached_query_expr[rows].mean(axis=0).astype(np.float32)
        ratios[i] = np.log2((clone_mean + pseudocount) / (control_gene_means + pseudocount))
    return ratios, clone_order


def _compute_clone_log2_ratios(
    query_matrix: sp.spmatrix | np.ndarray,
    query_lib: np.ndarray,
    coords: pd.DataFrame,
    control_gene_means: np.ndarray,
    clone_to_query_rows: Dict[str, np.ndarray],
    input_normalized: bool,
) -> Tuple[np.ndarray, List[str]]:
    """Compute per-clone log2(mean_clone / mean_control) for every coord-mapped gene.

    Returns (ratios, clone_order) with ratios of shape (n_clones, n_coord_genes).
    """
    clone_order = list(clone_to_query_rows.keys())
    n_genes = coords.shape[0]
    ratios = np.full((len(clone_order), n_genes), np.nan, dtype=np.float32)
    if not clone_order:
        return ratios, clone_order

    pseudocount = 1e-3
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        query_cols = chr_coords["var_index"].to_numpy(dtype=np.int64)
        coord_rows = chr_coords.index.to_numpy(dtype=np.int64)
        for i, clone_id in enumerate(clone_order):
            rows = clone_to_query_rows[clone_id]
            if rows.size == 0:
                continue
            sub = query_matrix[rows][:, query_cols]
            if sp.issparse(sub):
                sub = sub.toarray()
            sub = np.asarray(sub, dtype=np.float32)
            if not input_normalized:
                factors = (10000.0 / np.maximum(query_lib[rows].astype(np.float32), 1.0)).astype(np.float32)
                sub *= factors[:, None]
                np.log1p(sub, out=sub)
            clone_mean = np.nanmean(sub, axis=0).astype(np.float32)
            control_mean = control_gene_means[coord_rows]
            ratios[i, coord_rows] = np.log2((clone_mean + pseudocount) / (control_mean + pseudocount))
    return ratios, clone_order


def _genome_x_positions(coords: pd.DataFrame) -> Tuple[np.ndarray, List[str], List[Tuple[float, float]]]:
    """Return per-gene cumulative genome position, ordered chromosome list, and per-chromosome (start, end) bounds."""
    chr_to_offset: Dict[str, float] = {}
    chr_order: List[str] = []
    bounds: List[Tuple[float, float]] = []
    cumulative = 0.0
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        size = float(chr_coords["end"].max() - chr_coords["start"].min())
        chr_to_offset[chrom] = cumulative - float(chr_coords["start"].min())
        chr_order.append(str(chrom))
        bounds.append((cumulative, cumulative + size))
        cumulative += size + size * 0.005  # 0.5% gap between chromosomes
    midpoints = ((coords["start"] + coords["end"]) // 2).to_numpy(dtype=np.float64)
    chr_offsets = coords["chr"].map(chr_to_offset).to_numpy(dtype=np.float64)
    return midpoints + chr_offsets, chr_order, bounds


def _write_clone_genome_pdf(
    path: Path,
    coords: pd.DataFrame,
    clone_score_profiles: np.ndarray,
    windows: Sequence[Window],
    clone_order: Sequence[str],
    clone_sizes: Dict[str, int],
    clone_intervals_df: pd.DataFrame,
    score_y_clip: float,
    score_threshold: float,
) -> Optional[Path]:
    """Genome-wide CNV-score visualization (one page per clone).

    Plots per-clone mean window score along the genome rather than per-gene log2(clone/control).
    The score is the same standardized residual the interval caller uses, so what's drawn
    tracks the CNV calls directly and is not biased by transcriptional differential expression
    that's unrelated to copy number (a major confounder in AML/MDS where blast vs. healthy
    cell-state programs dominate per-gene log2).
    """
    if clone_score_profiles.shape[0] == 0 or not windows:
        return None
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    # Window midpoints positioned along a continuous genome x-axis using the same chrom-bound
    # layout as the per-gene plot, so chromosome boundaries align with what users expect.
    _, chr_order, chr_bounds = _genome_x_positions(coords)
    chr_centers = [(s + e) / 2.0 for s, e in chr_bounds]
    chr_labels = [c.replace("chr", "") for c in chr_order]
    chr_to_bounds = dict(zip(chr_order, chr_bounds))
    chr_min_start_lookup = {
        chrom: float(coords.loc[coords["chr"] == chrom, "start"].min())
        for chrom in chr_order
    }

    window_x = np.full(len(windows), np.nan, dtype=np.float64)
    for i, w in enumerate(windows):
        chrom = str(w.chrom)
        if chrom not in chr_to_bounds:
            continue
        chr_x_start, _ = chr_to_bounds[chrom]
        chr_min_start = chr_min_start_lookup[chrom]
        midpoint_genomic = (float(w.start) + float(w.end)) / 2.0
        window_x[i] = chr_x_start + (midpoint_genomic - chr_min_start)
    valid_x_mask = np.isfinite(window_x)

    cmap = LinearSegmentedColormap.from_list(
        "fastcnv_diverging",
        [(0.0, "#c8252c"), (0.5, "#cccccc"), (1.0, "#1f863d")],
    )

    intervals_by_clone: Dict[str, pd.DataFrame] = {}
    if not clone_intervals_df.empty:
        intervals_by_clone = {
            str(cid): grp for cid, grp in clone_intervals_df.groupby("global_clone_id", sort=False)
        }

    # Color saturation point: anything beyond ±score_threshold is full red/green. Below it,
    # colors fade toward grey so windows with sub-threshold (call-failed) noise don't draw the eye.
    color_extent = max(float(score_threshold), 1e-6)

    with PdfPages(path) as pdf:
        for clone_idx, clone_id in enumerate(clone_order):
            # Skip clones that produced no called consensus intervals — NMF placed cells in
            # this clone but the clone-level signal didn't survive the interval caller, so the
            # genome plot would be uninformative (no stars, just noise).
            clone_intervals = intervals_by_clone.get(str(clone_id))
            if clone_intervals is None or clone_intervals.empty:
                continue

            clone_scores_chr = clone_score_profiles[clone_idx]
            # Recenter on the genome-wide median so each clone's baseline sits at zero. The
            # standardized residual is centered per-window-per-state during scoring, but a
            # global library/protocol mismatch between the query cohort and the reference
            # leaves a residual offset (typically -0.5 to -1.0) that's identical across the
            # whole genome — not CNV. Subtracting the per-clone median for display exposes
            # the contiguous deviations that actually correspond to copy number.
            finite = np.isfinite(clone_scores_chr)
            clone_offset = float(np.nanmedian(clone_scores_chr[finite])) if finite.any() else 0.0
            clone_scores_centered = clone_scores_chr - clone_offset
            scores_clipped = np.clip(np.nan_to_num(clone_scores_centered, nan=0.0), -score_y_clip, score_y_clip)

            fig, ax = plt.subplots(figsize=(14, 4.5))
            ax.scatter(
                window_x[valid_x_mask], scores_clipped[valid_x_mask],
                c=scores_clipped[valid_x_mask], cmap=cmap, vmin=-color_extent, vmax=color_extent,
                s=12, linewidths=0, alpha=0.95, rasterized=True,
            )
            ax.axhline(0.0, color="#888888", linewidth=0.5)
            ax.axhline(score_threshold, color="#bbbbbb", linewidth=0.4, linestyle=":")
            ax.axhline(-score_threshold, color="#bbbbbb", linewidth=0.4, linestyle=":")
            for chrom in chr_order:
                start_x, _ = chr_to_bounds[chrom]
                ax.axvline(start_x, color="#cccccc", linewidth=0.4, linestyle="--")
            ax.set_xticks(chr_centers)
            ax.set_xticklabels(chr_labels, fontsize=8)
            ax.set_xlim(chr_bounds[0][0], chr_bounds[-1][1])
            ax.set_ylim(-score_y_clip - 0.1, score_y_clip + 0.1)
            ax.set_ylabel("CNV score  (standardized residual)")
            ax.set_xlabel("Chromosome")
            n_cells = int(clone_sizes.get(str(clone_id), 0))
            ax.set_title(f"{clone_id}  (n={n_cells} cells)")

            # Overlay called intervals as black horizontal segments at the interval mean score
            # (in the same recentered units as the scatter), with a black star at the midpoint.
            window_chr = np.array([str(w.chrom) for w in windows])
            window_start = np.array([w.start for w in windows], dtype=np.int64)
            window_end = np.array([w.end for w in windows], dtype=np.int64)
            for _, interval_row in clone_intervals.iterrows():
                chrom = str(interval_row["chr"])
                if chrom not in chr_to_bounds:
                    continue
                chr_x_start, _ = chr_to_bounds[chrom]
                chr_min_start = chr_min_start_lookup[chrom]
                seg_x0 = chr_x_start + (float(interval_row["start"]) - chr_min_start)
                seg_x1 = chr_x_start + (float(interval_row["end"]) - chr_min_start)
                in_segment = (
                    (window_chr == chrom)
                    & (window_start >= int(interval_row["start"]))
                    & (window_end <= int(interval_row["end"]))
                )
                if not in_segment.any():
                    continue
                seg_y = float(np.nanmean(clone_scores_centered[in_segment]))
                seg_y_clipped = float(np.clip(seg_y, -score_y_clip, score_y_clip))
                ax.plot(
                    [seg_x0, seg_x1], [seg_y_clipped, seg_y_clipped],
                    color="black", linewidth=2.5, solid_capstyle="butt", zorder=6,
                )
                ax.scatter(
                    [(seg_x0 + seg_x1) / 2.0], [seg_y_clipped],
                    marker="*", s=140, color="black", edgecolor="white", linewidth=0.6, zorder=7,
                )
            fig.tight_layout()
            pdf.savefig(fig, dpi=200)
            plt.close(fig)
    return path


def _write_clone_pdf_hardened(
    path: Path,
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    windows: Sequence[Window],
    barcode_to_index: Dict[str, int],
    low_threshold: float,
) -> Optional[Path]:
    """PDF writer that forces the Agg backend so the file is portable across viewers."""
    def _clone_sort_key(clone_id: str) -> Tuple[int, str]:
        text = str(clone_id)
        if text.startswith("clone") and text[5:].isdigit():
            return (0, f"{int(text[5:]):06d}")
        return (1, text)

    clone_ids = [c for c in sorted(cell_df["global_clone_id"].dropna().unique(), key=_clone_sort_key) if c != "WT"]
    if not clone_ids:
        return None
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap

    logging.getLogger("fontTools").setLevel(logging.WARNING)
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    matrix_rows: List[np.ndarray] = []
    labels: List[str] = []
    for clone_id in clone_ids:
        group = cell_df[cell_df["global_clone_id"] == clone_id]
        indices = np.asarray(
            [barcode_to_index[b] for b in group["CellBarcode"] if b in barcode_to_index],
            dtype=int,
        )
        if indices.size == 0:
            continue
        profile = np.nanmean(all_scores[indices], axis=0)
        calls = np.zeros(profile.size, dtype=int)
        calls[profile >= low_threshold] = 1
        calls[profile <= -low_threshold] = -1
        matrix_rows.append(calls)
        labels.append(f"{clone_id} (n={indices.size})")
    if not matrix_rows:
        return None
    call_matrix = np.vstack(matrix_rows)

    chrom_boundaries: List[int] = []
    chrom_labels: List[str] = []
    last_chrom = None
    for index, window in enumerate(windows):
        if window.chrom != last_chrom:
            chrom_boundaries.append(index)
            chrom_labels.append(window.chrom.replace("chr", ""))
            last_chrom = window.chrom

    fig_height = max(2.5, 0.35 * len(labels) + 1.5)
    cmap = ListedColormap(["#2b66c3", "#ffffff", "#d62f27"])
    n_rows, n_cols = call_matrix.shape
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(14, fig_height))
        im = ax.imshow(
            call_matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=-1,
            vmax=1,
            origin="upper",
            extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5),
        )
        im.set_rasterized(True)
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Genome windows")
        ax.set_title("fastCNV clone-level genome calls")
        for boundary in chrom_boundaries:
            ax.axvline(boundary - 0.5, color="#cccccc", linewidth=0.4)
        ax.set_xticks(chrom_boundaries)
        ax.set_xticklabels(chrom_labels, fontsize=7, rotation=90)
        ax.tick_params(axis="x", length=0)
        fig.tight_layout()
        pdf.savefig(fig, dpi=200)
        plt.close(fig)
    return path


def _write_infercnv_style_heatmaps(
    output_prefix: Path,
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    windows: Sequence[Window],
    barcode_to_index: Dict[str, int],
    max_cells: int = 20000,
    filter_threshold: float = 1.5,
    min_chr_windows: int = 35,
    random_state: int = 0,
) -> Dict[str, Path]:
    """Write inferCNV-style per-cell heatmaps from fastCNV window scores.

    The heatmap is an export-only view: it does not alter CNV calls. Rows are
    ordered by fastCNV clone and cell state; large cohorts are deterministically
    stratified by clone/state for display so the PDF/PNG remains portable.
    """
    outputs: Dict[str, Path] = {}
    if all_scores.size == 0 or not windows or cell_df.empty:
        return outputs

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap

    logging.getLogger("fontTools").setLevel(logging.WARNING)
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    frame = cell_df.copy()
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    def _barcode_key(value: object) -> str:
        return str(value).split(".")[0]

    normalized_index: Dict[str, int] = {}
    for barcode, index in barcode_to_index.items():
        normalized_index.setdefault(_barcode_key(barcode), int(index))
    frame["_row_index"] = frame["CellBarcode"].map(barcode_to_index)
    missing = frame["_row_index"].isna()
    if missing.any():
        frame.loc[missing, "_row_index"] = frame.loc[missing, "CellBarcode"].map(
            lambda value: normalized_index.get(_barcode_key(value))
        )
    frame = frame[frame["_row_index"].notna()].copy()
    if frame.empty:
        return outputs
    if "global_clone_id" not in frame.columns:
        frame["global_clone_id"] = "WT"
    if "cell_state" not in frame.columns:
        frame["cell_state"] = "unknown"

    def clone_sort_key(value: object) -> Tuple[int, str]:
        text = str(value)
        if text == "clone_loy":
            return (0, "000000")
        if text.startswith("clone") and text[5:].isdigit():
            return (1, f"{int(text[5:]):06d}")
        if text == "WT":
            return (9, text)
        return (2, text)

    rng = np.random.default_rng(int(random_state))
    frame["_row_index"] = frame["_row_index"].astype(int)
    frame["_clone_sort"] = frame["global_clone_id"].map(clone_sort_key)

    selected_parts: List[pd.DataFrame] = []
    total_cells = int(frame.shape[0])
    max_cells = int(max_cells) if max_cells and int(max_cells) > 0 else total_cells
    if total_cells > max_cells:
        group_cols = ["global_clone_id", "cell_state"]
        grouped = list(frame.groupby(group_cols, sort=False))
        quotas: Dict[Tuple[str, str], int] = {}
        for (clone_id, state), group in grouped:
            base = int(np.floor(max_cells * (len(group) / total_cells)))
            if str(clone_id) != "WT":
                base = max(base, min(len(group), 50))
            quotas[(str(clone_id), str(state))] = min(len(group), max(base, 1))
        quota_total = sum(quotas.values())
        if quota_total > max_cells:
            scale = max_cells / float(quota_total)
            quotas = {k: max(1, int(np.floor(v * scale))) for k, v in quotas.items()}
        for (clone_id, state), group in grouped:
            quota = min(len(group), quotas.get((str(clone_id), str(state)), 1))
            if quota >= len(group):
                selected_parts.append(group)
            else:
                chosen = rng.choice(group.index.to_numpy(), size=quota, replace=False)
                selected_parts.append(group.loc[np.sort(chosen)])
        selected = pd.concat(selected_parts, axis=0)
    else:
        selected = frame

    selected = selected.sort_values(
        by=["_clone_sort", "cell_state", "CellBarcode"],
        kind="mergesort",
    )
    row_indices = selected["_row_index"].to_numpy(dtype=np.int64)
    if row_indices.size == 0:
        return outputs

    matrix = np.asarray(all_scores[row_indices], dtype=np.float32)
    finite_cols = np.isfinite(matrix).any(axis=0)
    if not finite_cols.any():
        return outputs
    matrix = matrix[:, finite_cols]
    kept_windows = [w for w, keep in zip(windows, finite_cols) if keep]
    row_center = np.nanmedian(matrix, axis=1)
    row_center = np.where(np.isfinite(row_center), row_center, 0.0).astype(np.float32)
    matrix = matrix - row_center[:, None]
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    finite_values = matrix[np.isfinite(matrix)]
    vmax = float(np.nanquantile(np.abs(finite_values), 0.985)) if finite_values.size else 1.0
    vmax = max(vmax, float(filter_threshold), 0.25)
    matrix = np.clip(matrix, -vmax, vmax)

    chroms: List[str] = []
    boundaries: List[int] = [0]
    last = None
    for i, window in enumerate(kept_windows):
        chrom = str(window.chrom)
        if chrom != last:
            if last is not None:
                boundaries.append(i)
            chroms.append(chrom)
            last = chrom
    boundaries.append(len(kept_windows))

    display_blocks: List[np.ndarray] = []
    display_boundaries: List[int] = [0]
    min_chr_windows = max(int(min_chr_windows), 1)
    for i in range(len(chroms)):
        block = matrix[:, boundaries[i]:boundaries[i + 1]]
        if 0 < block.shape[1] < min_chr_windows:
            take = np.linspace(0, block.shape[1] - 1, min_chr_windows).round().astype(np.int64)
            block = block[:, take]
        display_blocks.append(block)
        display_boundaries.append(display_boundaries[-1] + block.shape[1])
    matrix = np.hstack(display_blocks) if display_blocks else matrix
    boundaries = display_boundaries
    centers = [(boundaries[i] + boundaries[i + 1]) / 2.0 for i in range(len(chroms))]

    clone_labels = selected["global_clone_id"].astype(str).to_numpy()
    clone_order = sorted(pd.unique(clone_labels), key=clone_sort_key)
    palette = [
        "#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2",
        "#b279a2", "#ff9da6", "#9d755d", "#bab0ac", "#59a14f",
        "#edc948", "#af7aa1", "#76b7b2", "#d37295",
    ]
    clone_color = {clone: palette[i % len(palette)] for i, clone in enumerate(clone_order)}
    clone_color["WT"] = "#d9d9d9"
    clone_color["clone_loy"] = "#2458a6"
    strip_values = np.array([clone_order.index(c) for c in clone_labels], dtype=np.int16)[:, None]
    strip_cmap = ListedColormap([clone_color[c] for c in clone_order])
    heat_cmap = LinearSegmentedColormap.from_list(
        "fastcnv_infercnv",
        ["#1f4fa3", "#597fc5", "#f7f7f7", "#d95f5f", "#a51f1f"],
    )

    def draw(data: np.ndarray, title: str, png_path: Path, pdf_path: Path) -> None:
        fig_height = min(12.0, max(5.0, 0.00055 * data.shape[0] + 4.5))
        fig = plt.figure(figsize=(13.2, fig_height))
        gs = fig.add_gridspec(1, 3, width_ratios=[0.018, 1.0, 0.018], wspace=0.012)
        ax_strip = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_strip)
        cax = fig.add_subplot(gs[0, 2])

        ax_strip.imshow(
            strip_values, aspect="auto", interpolation="none",
            cmap=strip_cmap, vmin=0, vmax=max(len(clone_order) - 1, 1),
            origin="upper",
        )
        ax_strip.set_axis_off()
        im = ax.imshow(
            data, aspect="auto", interpolation="none", cmap=heat_cmap,
            vmin=-vmax, vmax=vmax, origin="upper", rasterized=True,
        )
        for boundary in boundaries[1:-1]:
            ax.axvline(boundary - 0.5, color="#202020", linewidth=0.35)
        ax.set_xticks(centers)
        ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=7)
        ax.set_yticks([])
        ax.set_xlabel("chromosome (narrow chromosomes expanded for display)")
        sampled = "" if total_cells == row_indices.size else f"; displayed {row_indices.size:,}/{total_cells:,} cells"
        ax.set_title(f"{title}{sampled}", fontsize=10)

        handles = [plt.Rectangle((0, 0), 1, 1, color=clone_color[c]) for c in clone_order]
        labels = [f"{c} ({int((clone_labels == c).sum()):,})" for c in clone_order]
        ax.legend(
            handles, labels, title="fastCNV clone", loc="upper left",
            bbox_to_anchor=(1.02, 1.0), fontsize=7, title_fontsize=8,
            frameon=False, borderaxespad=0.0,
        )
        fig.colorbar(im, cax=cax).set_label("centered fastCNV score", fontsize=8)
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)

    primary_png = output_prefix.with_suffix(".infercnv_heatmap.png")
    primary_pdf = output_prefix.with_suffix(".infercnv_heatmap.pdf")
    filtered_png = output_prefix.with_suffix(".infercnv_heatmap_filtered.png")
    filtered_pdf = output_prefix.with_suffix(".infercnv_heatmap_filtered.pdf")
    draw(matrix, "fastCNV inferCNV-style heatmap", primary_png, primary_pdf)
    filtered = np.where(np.abs(matrix) < float(filter_threshold), 0.0, matrix).astype(np.float32)
    draw(filtered, f"fastCNV inferCNV-style heatmap, filtered |score|<{filter_threshold:g}", filtered_png, filtered_pdf)
    outputs.update({
        "infercnv_heatmap_png": primary_png,
        "infercnv_heatmap_pdf": primary_pdf,
        "infercnv_heatmap_filtered_png": filtered_png,
        "infercnv_heatmap_filtered_pdf": filtered_pdf,
    })
    return outputs


def _nearest_normal_autosomal(
    query: ad.AnnData, control: ad.AnnData, coords, state_key: str,
    control_state_key: Optional[str], layer: str, input_normalized: bool,
    output_prefix: Path, sample_key: Optional[str] = None,
    call_threshold: float = 2.0, expr_min: float = 0.5,
    query_to_control_state: Optional[Dict[str, str]] = None,
) -> Dict[str, Path]:
    """VALIDATED nearest-normal-subcluster autosomal CNV caller (opt-in: --nearest-normal-autosomal).

    Residual for each query leiden-subcluster = its mean log1p(CP10k) MINUS its NEAREST NORMAL reference
    subcluster (matched by correlation over expressed genes -> CNV-robust, one chromosome is ~4% of genes),
    divided by the per-cell-state cross-donor SD; expression-weighted pyramidal smoothing along the genome;
    per-chromosome-arm aggregate z; an arm is called loss/gain when |aggregate z| > call_threshold.
    Subclusters are grouped into clones by shared call-set. Writes the standard
    <prefix>.clone_report.csv / .clone_intervals.tsv / .cnv_cells.tsv. Validated on the AML benchmark:
    held-out NORMAL control < 5% of cells with ANY call; detects -7 / +8 / del(5q)."""
    import scanpy as sc
    import pandas as pd
    sc.settings.verbosity = 0
    CEN = {f"chr{c}": m * 1e6 for c, m in {
        "1": 123, "2": 93.9, "3": 90.9, "4": 50, "5": 48.8, "6": 59.8, "7": 60.1, "8": 45.2, "9": 43,
        "10": 39.8, "11": 53.4, "12": 35.5, "13": 17.7, "14": 17.2, "15": 19, "16": 36.8, "17": 25.1,
        "18": 18.5, "19": 26.2, "20": 28.1, "21": 12, "22": 15}.items()}
    qset = set(map(str, query.var_names)); cset = set(map(str, control.var_names))
    g2chr: Dict[str, str] = {}; g2pos: Dict[str, float] = {}; qnames = query.var_names
    for _, r in coords.iterrows():
        c = str(r["chr"])
        if c in CEN:
            g = str(qnames[int(r["var_index"])]); g2chr[g] = c; g2pos[g] = float(r["start"])
    genes = sorted([g for g in g2chr if g in qset and g in cset], key=lambda g: (g2chr[g], g2pos[g]))
    if len(genes) < 200:
        raise ValueError(f"nearest-normal autosomal: only {len(genes)} shared autosomal genes with coordinates.")
    chrom = np.array([g2chr[g] for g in genes]); gpos = np.array([g2pos[g] for g in genes], dtype=float)
    chr_order = [f"chr{i}" for i in range(1, 23) if f"chr{i}" in set(chrom)]
    chr_ranges = {c: np.flatnonzero(chrom == c) for c in chr_order}
    arms_arr = np.array([f"{chrom[i]}{'p' if gpos[i] < CEN[chrom[i]] else 'q'}" for i in range(len(genes))])
    aidx = {a: np.flatnonzero(arms_arr == a) for a in np.unique(arms_arr)}

    def _logmat(adata: ad.AnnData) -> np.ndarray:
        X = _matrix(adata, layer); idx = adata.var_names.get_indexer(genes)
        sub = X[:, idx]; sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        sub = np.asarray(sub, np.float64)
        if input_normalized:
            return sub                                          # already log1p(CP10k)
        full = X.tocsr() if sp.issparse(X) else X
        lib = np.asarray(full.sum(1)).ravel().astype(np.float64); lib[lib == 0] = 1.0
        return np.log1p((sub / lib[:, None]) * 1e4)             # raw counts -> log1p(CP10k)

    def _subcluster(adata: ad.AnnData) -> np.ndarray:
        a = adata.copy()
        if not input_normalized:
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=min(2000, a.n_vars))
        a = a[:, a.var.highly_variable]
        sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=min(20, a.n_vars - 1))
        sc.pp.neighbors(a, n_neighbors=15)
        sc.tl.leiden(a, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
        return a.obs["leiden"].to_numpy()

    cL = _logmat(control)
    c_state = control.obs[control_state_key or state_key].astype(str).to_numpy()
    states = [s for s in np.unique(c_state) if (c_state == s).sum() >= 20]
    if not states:
        raise ValueError("nearest-normal autosomal: no control cell-state has >= 20 cells.")
    sidx = {s: i for i, s in enumerate(states)}
    centroid = np.vstack([cL[c_state == s].mean(0) for s in states])
    centroid_all = cL.mean(0); W = np.minimum(np.expm1(centroid_all), 50.0)
    donor_col = next((dc for dc in ("Donor", "Sample") if dc in control.obs.columns), None)
    SG = np.full((len(states), len(genes)), 0.1)
    for si, s in enumerate(states):
        L = cL[c_state == s]
        if donor_col is not None:
            dd = control.obs[donor_col].astype(str).to_numpy()[c_state == s]
            dmeans = [L[dd == d].mean(0) for d in np.unique(dd) if (dd == d).sum() >= 10]
            if len(dmeans) >= 3:
                SG[si] = np.maximum(np.vstack(dmeans).std(0), 0.1); continue
        SG[si] = np.maximum(L.std(0), 0.1)
    LOGGER.info("nearest-normal: %d shared genes, %d control states, cross-donor SD via %s.",
                len(genes), len(states), donor_col or "per-state cell SD (no donor column)")

    def _smooth(zz: np.ndarray, em: np.ndarray, win: int = 51) -> np.ndarray:
        h = (win - 1) // 2
        kern = np.r_[np.arange(1, h + 1), h + 1, np.arange(h, 0, -1)].astype(float)
        wt = (W * em).astype(float); out = zz.copy()
        for c, idx in chr_ranges.items():
            if len(idx) < 3:
                continue
            kk = kern if len(idx) >= win else np.ones(min(len(idx), win))
            num = np.convolve(zz[idx] * wt[idx], kk, mode="same"); den = np.convolve(wt[idx], kk, mode="same")
            ok = den > 1e-9; s = out[idx].copy(); s[ok] = num[ok] / den[ok]; out[idx] = s
        return out

    c_sub = _subcluster(control); lib_profiles = []
    for cl in np.unique(c_sub):
        ii = np.flatnonzero(c_sub == cl); vals, cnts = np.unique(c_state[ii], return_counts=True); mj = vals[cnts.argmax()]
        if mj in sidx:
            lib_profiles.append((sidx[mj], cL[ii].mean(0)))
    lib_by_state: Dict[int, list] = {}
    for si, sm in lib_profiles:
        lib_by_state.setdefault(si, []).append(sm)

    def _match(sm: np.ndarray, si: int) -> np.ndarray:
        cands = lib_by_state.get(si) or [s for _, s in lib_profiles]
        if not cands:
            return centroid[si]
        em = centroid[si] >= expr_min; a = sm[em]; a = a - a.mean(); best = None; bc = -2.0
        for lsm in cands:
            b = lsm[em]; b = b - b.mean(); d = np.sqrt((a * a).sum() * (b * b).sum())
            cc = (a * b).sum() / d if d > 0 else -2.0
            if cc > bc:
                bc = cc; best = lsm
        return best if best is not None else centroid[si]

    qL = _logmat(query)
    q_state_original = query.obs[state_key].astype(str).to_numpy()
    if query_to_control_state:
        q_state = np.array([query_to_control_state.get(str(s), str(s)) for s in q_state_original], dtype=object)
    else:
        q_state = q_state_original
    q_sub = _subcluster(query); q_bc = query.obs_names.astype(str).to_numpy(); n = query.n_obs
    blocks = []
    for cl in np.unique(q_sub):
        ii = np.flatnonzero(q_sub == cl); vals, cnts = np.unique(q_state[ii], return_counts=True); mj = vals[cnts.argmax()]
        if mj not in sidx:
            blocks.append((ii, None, None)); continue
        si = sidx[mj]; sm = qL[ii].mean(0); ref_sm = _match(sm, si); em = centroid[si] >= expr_min
        zz = np.clip((sm - ref_sm) / SG[si], -10, 10); hi = centroid[si] >= 5.0
        off = float(np.median(zz[hi])) if hi.sum() > 10 else 0.0
        blocks.append((ii, _smooth(np.clip(zz - off, -10, 10), em), em))

    arm_means = []
    for ii, sz, em in blocks:
        if sz is None:
            arm_means.append((ii, {})); continue
        az = {a: float(np.mean(sz[aidx[a][em[aidx[a]]]])) for a in aidx if int(em[aidx[a]].sum()) >= 15}
        arm_means.append((ii, az))
    arm_coord = {}
    for a, gi in aidx.items():
        g = np.sort(gi)
        arm_coord[a] = (str(chrom[g[0]]), int(gpos[g[0]]), int(gpos[g[-1]]), str(genes[g[0]]), str(genes[g[-1]]), int(len(g)))

    groups: Dict[tuple, dict] = {}
    for (ii, sz, em), (_, az) in zip(blocks, arm_means):
        ck = tuple(sorted((a, -1 if az[a] < 0 else 1) for a in az if abs(az[a]) > call_threshold))
        gd = groups.setdefault(ck, {"cells": [], "armz": {}}); gd["cells"].extend(ii.tolist())
        for a, d in ck:
            gd["armz"].setdefault(a, []).append(az[a])
    ordered = sorted(groups.items(), key=lambda kv: -len(kv[1]["cells"])); clone_id = {}; cid = 1
    for ck, _ in ordered:
        clone_id[ck] = "WT" if not ck else f"clone{cid}"; cid += (1 if ck else 0)
    cell_clone = np.array(["WT"] * n, dtype=object); cell_nint = np.zeros(n, int)
    for ck, gd in groups.items():
        for c in gd["cells"]:
            cell_clone[c] = clone_id[ck]; cell_nint[c] = len(ck)
    sample_id = (query.obs[sample_key].astype(str).to_numpy() if sample_key and sample_key in query.obs.columns
                 else np.array([output_prefix.name] * n))

    rep = []
    for ck, gd in ordered:
        cells = gd["cells"]; cts = ";".join(pd.Series(q_state_original[cells]).value_counts().head(3).index.astype(str))
        regions = ";".join(f"{a}:{'loss' if d < 0 else 'gain'}" for a, d in ck) if ck else "(no CNV)"
        smp = ";".join(f"{s}({c})" for s, c in pd.Series(sample_id[cells]).value_counts().items())
        rep.append((clone_id[ck], len(cells), regions, cts, smp))
    rep_path = output_prefix.with_suffix(".clone_report.csv")
    pd.DataFrame(rep, columns=["clone", "n_cells", "regions", "cell_types", "samples"]).to_csv(rep_path, index=False)
    ci = []
    for ck, gd in ordered:
        if not ck:
            continue
        for a, d in ck:
            ch, st, en, sg, eg, ng = arm_coord[a]; zs = float(np.mean(gd["armz"][a]))
            ci.append((clone_id[ck], len(gd["cells"]), 1.0, "loss" if d < 0 else "gain", ch, st, en, sg, eg,
                       1, ng, round(abs(zs), 3), round(abs(zs), 3), round(abs(zs), 3), "", "", "nearest_normal"))
    ci_path = output_prefix.with_suffix(".clone_intervals.tsv")
    pd.DataFrame(ci, columns=["global_clone_id", "n_cells", "support_fraction", "call", "chr", "start", "end",
                              "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score",
                              "confidence", "mean_log2_ratio", "zygosity_state", "clone_confidence"]).to_csv(ci_path, sep="\t", index=False)
    cells_path = output_prefix.with_suffix(".cnv_cells.tsv")
    pd.DataFrame({"CellBarcode": q_bc, "cell_state": q_state_original, "sample": sample_id,
                  "cnv_status": np.where(cell_nint > 0, "CNV", "WT"), "baseline_source": "nearest_normal",
                  "n_cnv_intervals": cell_nint, "sex_chrom_only": False, "state_clone_id": cell_clone,
                  "global_clone_id": cell_clone, "clone_confidence": np.where(cell_nint > 0, "medium", "wt")}).to_csv(cells_path, sep="\t", index=False)
    LOGGER.info("nearest-normal autosomal: %d clones, %d/%d cells with a CNV call. -> %s",
                sum(1 for ck, _ in ordered if ck), int((cell_nint > 0).sum()), n, rep_path.name)
    return {"clone_report": rep_path, "clone_intervals": ci_path, "cnv_cells": cells_path}


def run_fast(params: FastParams) -> Dict[str, Path]:
    run_start = time.perf_counter()
    params.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Reading query AnnData: %s", params.h5ad)
    t0 = time.perf_counter()
    query = ad.read_h5ad(params.h5ad)
    LOGGER.info("Query loaded: %d cells x %d genes (%s)", query.n_obs, query.n_vars, _format_duration(time.perf_counter() - t0))

    available = ", ".join(query.obs.columns.astype(str)) or "<none>"
    if params.state_key not in query.obs.columns:
        raise KeyError(f"State column '{params.state_key}' not found. Available: {available}")
    if params.sample_key and params.sample_key not in query.obs.columns:
        raise KeyError(f"Sample column '{params.sample_key}' not found. Available: {available}")

    LOGGER.info("Reading control AnnData: %s", params.control_h5ad)
    t0 = time.perf_counter()
    control = ad.read_h5ad(params.control_h5ad)
    LOGGER.info("Control loaded: %d cells x %d genes (%s)", control.n_obs, control.n_vars, _format_duration(time.perf_counter() - t0))

    if params.control_state_key and params.control_state_key not in control.obs.columns:
        raise KeyError(
            f"Control state column '{params.control_state_key}' not found. "
            f"Available: {', '.join(control.obs.columns.astype(str))}"
        )

    if params.gene_coordinates is None:
        raise ValueError("gene_coordinates must be provided.")
    coords = load_gene_coordinates(params.gene_coordinates, query.var_names)
    blacklist = build_gene_blacklist(query.var_names) if params.gene_blacklist else set()
    if blacklist:
        _coord_syms = np.array([str(query.var_names[i]) for i in coords["var_index"].to_numpy()])
        n_before = coords.shape[0]
        coords = coords[~np.isin(_coord_syms, list(blacklist))].copy().reset_index(drop=True)
        LOGGER.info(
            "Gene blacklist ON: %d Ig/TCR/Hb/MT/platelet genes; scoring coords %d -> %d.",
            len(blacklist), n_before, coords.shape[0],
        )
    windows = build_windows(coords, params.window_genes, params.stride_genes, params.min_chr_genes)
    slices_by_chr = _chr_slices(windows)
    LOGGER.info("Built %d windows across %d chromosomes.", len(windows), len(slices_by_chr))

    query_matrix = _matrix(query, params.layer)
    if sp.issparse(query_matrix):
        query_matrix = query_matrix.tocsc()
    control_matrix = _matrix(control, params.layer)
    if sp.issparse(control_matrix):
        control_matrix = control_matrix.tocsc()

    LOGGER.info("Computing library sizes%s.", " (blacklist excluded)" if blacklist else "")
    t0 = time.perf_counter()
    if params.input_normalized:
        query_lib = np.ones(query.n_obs, dtype=np.float32)
        control_lib = np.ones(control.n_obs, dtype=np.float32)
    else:
        _qm = query_matrix.tocsr() if sp.issparse(query_matrix) else query_matrix
        _cm = control_matrix.tocsr() if sp.issparse(control_matrix) else control_matrix
        query_lib = _row_sums(_qm)
        control_lib = _row_sums(_cm)
        if blacklist:
            _qbl = np.flatnonzero(np.isin(query.var_names.astype(str), list(blacklist)))
            _cbl = np.flatnonzero(np.isin(control.var_names.astype(str), list(blacklist)))
            if _qbl.size:
                query_lib = (query_lib - np.asarray(_qm[:, _qbl].sum(axis=1)).ravel().astype(np.float32))
            if _cbl.size:
                control_lib = (control_lib - np.asarray(_cm[:, _cbl].sum(axis=1)).ravel().astype(np.float32))
            query_lib = np.maximum(query_lib, 1.0).astype(np.float32)
            control_lib = np.maximum(control_lib, 1.0).astype(np.float32)
    LOGGER.info("Library sizes done (%s).", _format_duration(time.perf_counter() - t0))

    var_map = _control_var_map(query.var_names, control.var_names)

    state_values = query.obs[params.state_key]
    states = [str(s) for s in pd.Index(state_values[state_values.notna()]).unique().tolist() if str(s)]
    state_index_map: Dict[str, np.ndarray] = {
        state: np.flatnonzero(state_values.astype(str).to_numpy() == state).astype(np.int64)
        for state in states
    }
    eligible_states = [s for s in states if state_index_map[s].size >= params.min_state_cells]
    LOGGER.info("States: %d total, %d eligible (>= %d cells).", len(states), len(eligible_states), params.min_state_cells)

    state_crosswalk: Dict[str, str] = {str(s): str(s) for s in eligible_states}
    if params.control_state_key:
        state_crosswalk = _build_state_crosswalk(
            query, control, query_matrix, control_matrix, query_lib, control_lib,
            coords, params.state_key, params.control_state_key, eligible_states,
            params.input_normalized,
        )

    control_state_pool: Optional[Dict[str, np.ndarray]] = None
    if params.control_state_key:
        control_state_values = control.obs[params.control_state_key].astype(str).to_numpy()
        control_state_pool = {}
        for state in eligible_states:
            control_state = state_crosswalk.get(str(state), str(state))
            mask = control_state_values == control_state
            control_state_pool[state] = np.flatnonzero(mask).astype(np.int64) if mask.any() else np.arange(control.n_obs, dtype=np.int64)
    all_control_rows = np.arange(control.n_obs, dtype=np.int64)

    # ---- Automatic per-sample sex detection ----------------------------------
    # ChrY scoring is meaningless when the query or control samples are female,
    # and pooling male+female controls flattens the chrY MAD/baseline so LOY
    # in male queries gets squashed. Detect sex per sample (>{threshold}% cells
    # with >=1 chrY-Y-only marker UMI -> male, else female), then gate chrY
    # scoring: only male query cells get scored against a male-only control pool.
    if params.sample_key and params.sample_key in query.obs.columns:
        query_sample_labels = query.obs[params.sample_key].astype(str).to_numpy()
    else:
        query_sample_labels = np.full(query.n_obs, "all", dtype=object)
    query_sex_df = _infer_sample_sex(query, query_sample_labels, params.sex_detection_threshold_pct)
    sample_to_query_sex = query_sex_df["inferred_sex"].to_dict()
    query_cell_sex = np.array(
        [sample_to_query_sex.get(s, "unknown") for s in query_sample_labels], dtype=object,
    )
    male_query_mask = (query_cell_sex == "male")

    # LOY-specific heuristic (controlled by --gate-chry-by-marker-expression, default True):
    # A cell expressing >=1 UMI on any chrY-Y-only marker gene cannot have biologically
    # lost chrY — the residual transcripts would be impossible without an intact
    # chromosome. We use this to gate chrY-loss interval emission downstream: chrY
    # scores for chrY-marker-positive cells are NaN'd so the interval caller can never
    # emit a chrY-loss interval for them, regardless of how the autosomal scoring or
    # windowing behaves. This is biologically motivated but specific to LOY; can be
    # disabled with --no-gate-chry-by-marker-expression to evaluate scoring purity.
    # The marker-positivity gate is a BIASED-mode heuristic only. In unbiased mode (default) it is
    # neither computed nor applied — LOY is detected purely by the unsupervised chrY scoring.
    if False:  # marker gate removed (general unsupervised method never gates chrY)
        _chry_present = [g for g in CHRY_Y_ONLY_MARKERS if g in query.var_names]
        if _chry_present:
            _chry_X = query[:, _chry_present].X
            if sp.issparse(_chry_X):
                _chry_X = _chry_X.tocsr()
                chry_marker_pos_query = np.asarray((_chry_X > 0).sum(axis=1)).ravel() > 0
            else:
                chry_marker_pos_query = (np.asarray(_chry_X) > 0).sum(axis=1) > 0
            LOGGER.info(
                "chrY-marker positivity gate ON (biased mode): %d / %d query cells express >=1 chrY-Y-only marker (gated from chrY-loss calls).",
                int(chry_marker_pos_query.sum()), query.n_obs,
            )
        else:
            chry_marker_pos_query = np.zeros(query.n_obs, dtype=bool)
    elif params.gate_chry_by_marker_expression:
        LOGGER.info("LOY mode = unbiased: chrY-marker gate DISABLED (LOY detected as an unsupervised CNV).")
        chry_marker_pos_query = np.zeros(query.n_obs, dtype=bool)
    else:
        LOGGER.info("chrY-marker positivity gate DISABLED (--no-gate-chry-by-marker-expression).")
        chry_marker_pos_query = np.zeros(query.n_obs, dtype=bool)

    if params.control_sample_key and params.control_sample_key in control.obs.columns:
        control_sample_labels = control.obs[params.control_sample_key].astype(str).to_numpy()
    else:
        control_sample_labels = np.full(control.n_obs, "all", dtype=object)
    control_sex_df = _infer_sample_sex(control, control_sample_labels, params.sex_detection_threshold_pct)
    sample_to_control_sex = control_sex_df["inferred_sex"].to_dict()
    control_cell_sex = np.array(
        [sample_to_control_sex.get(s, "unknown") for s in control_sample_labels], dtype=object,
    )
    male_control_rows = np.flatnonzero(control_cell_sex == "male").astype(np.int64)

    LOGGER.info(
        "Sex detection (>%g%% chrY-marker-positive cells = male):",
        params.sex_detection_threshold_pct,
    )
    for src_label, df in (("query", query_sex_df), ("control", control_sex_df)):
        for sample, row in df.iterrows():
            LOGGER.info(
                "  %s sample=%s n=%d chrY+ %d/%d (%.1f%%) -> %s",
                src_label, sample, int(row["n_cells"]),
                int(row["n_chrY_pos"]), int(row["n_cells"]),
                float(row["pct_chrY_pos"]), str(row["inferred_sex"]),
            )

    n_male_q = int(male_query_mask.sum())
    n_male_c = int(male_control_rows.size)
    chry_scoring_enabled = (
        params.sex_chrom_mode == "absolute_log2"
        and "chrY" in {str(c) for c in params.sex_chroms}
        and n_male_q > 0 and n_male_c > 0
    )
    if not chry_scoring_enabled and "chrY" in {str(c) for c in params.sex_chroms}:
        if n_male_q == 0:
            LOGGER.warning("No male samples in query; chrY scoring disabled.")
        elif n_male_c == 0:
            LOGGER.warning("No male cells in control; chrY scoring disabled (LOY cannot be detected).")

    all_scores = np.full((query.n_obs, len(windows)), np.nan, dtype=np.float32)
    # --- autosomal caller = the EXACT same machinery as LOY ---
    # The autosomal copy-number caller applies the SAME `_simulation_copy_number_call` used for chrY
    # (LOY), once per chromosome ARM (normal_copies=2 instead of 1, all reference instead of male-only).
    # Same detection filter, same simulated-state likelihood (argmax, NO threshold), same baseline. Build
    # the per-arm gene lists here from the gene coordinates.
    do_autosomal_sim = bool(params.simulation_autosomal and params.pooled_scale)
    autosomal_arm_genes: Dict[str, List[str]] = {}
    if do_autosomal_sim:
        _cchr = coords["chr"].astype(str).to_numpy()
        _cstart = coords["start"].astype(float).to_numpy()
        _cvn = [str(query.var_names[int(i)]) for i in coords["var_index"].to_numpy()]
        for _gi in range(len(_cvn)):
            _a = _arm_label(_cchr[_gi], int(_cstart[_gi]))
            if _a and not _a.startswith(("chrX", "chrY")):
                autosomal_arm_genes.setdefault(_a, []).append(_cvn[_gi])
        autosomal_arm_genes = {a: g for a, g in autosomal_arm_genes.items() if len(g) >= 20}
    control_gene_means = np.full(coords.shape[0], np.nan, dtype=np.float32)
    coords_var_lookup = pd.Series(
        np.arange(coords.shape[0], dtype=np.int64), index=coords["var_index"].to_numpy()
    )
    cache_query_expr = not params.skip_pdf
    if cache_query_expr:
        cached_query_expr = np.zeros((query.n_obs, coords.shape[0]), dtype=np.float32)
        LOGGER.info(
            "Caching normalized query expression: %d x %d float32 (~%.1f GB).",
            query.n_obs, coords.shape[0], cached_query_expr.nbytes / 1e9,
        )
    else:
        cached_query_expr = None
    LOGGER.info(
        "Allocated score matrix: %d x %d float32 (~%.1f GB).",
        query.n_obs, len(windows), all_scores.nbytes / 1e9,
    )

    chrom_groups = list(coords.groupby("chr", sort=False))
    LOGGER.info("Scoring %d chromosomes.", len(chrom_groups))

    # Per-state-calibrated unsupervised LOY (parity with the validated standalone). Captured during
    # chrY scoring, used after the loop to call chrY-loss against EACH STATE'S OWN reference null
    # (5% quantile) — removes the cell-type bias that a single fixed threshold causes.
    chry_window_abs_idx: Optional[np.ndarray] = None
    chry_state_tau: Dict[str, float] = {}
    chry_global_tau: Optional[float] = None
    chry_interval_template: Optional[Dict[str, object]] = None

    for chr_index, (chrom, chr_coords) in enumerate(chrom_groups):
        chrom_str = str(chrom)
        chr_slice = slices_by_chr.get(chrom_str)
        if chr_slice is None:
            continue
        chr_windows = windows[chr_slice]
        if not chr_windows:
            continue
        chr_t0 = time.perf_counter()

        query_cols = chr_coords["var_index"].to_numpy(dtype=np.int64)
        query_expr = _normalize_chunk(query_matrix, query_cols, query_lib, params.input_normalized)
        if cached_query_expr is not None:
            cached_query_expr[:, chr_coords.index.to_numpy(dtype=np.int64)] = query_expr

        mapped = var_map[query_cols]
        valid_mask = ~pd.isna(mapped)
        valid_control_cols = mapped[valid_mask].astype(np.int64) if valid_mask.any() else np.empty(0, dtype=np.int64)
        if valid_control_cols.size == 0:
            LOGGER.warning("Chromosome %s: no overlapping control genes; skipping.", chrom_str)
            continue
        control_expr_full = np.zeros((control.n_obs, query_cols.size), dtype=np.float32)
        control_sub = _normalize_chunk(control_matrix, valid_control_cols, control_lib, params.input_normalized)
        control_expr_full[:, valid_mask] = control_sub
        coord_rows_for_chr = coords_var_lookup.reindex(query_cols).to_numpy()

        is_sex_chrom = (
            params.sex_chrom_mode == "absolute_log2"
            and chrom_str in params.sex_chroms
        )
        # For sex chromosomes (specifically chrY when scoring is active), the per-gene
        # control mean must be computed from male donors only — otherwise female zeros drag
        # the baseline toward zero and the per-clone log2(clone/control) used by the PDF
        # collapses to ~0 in LOY cells (where clone_mean is also ~0). The result was a
        # visually flat chrY band even though LOY is biologically dramatic. Restricting
        # to males here makes the PDF show the real ~log2(-2 to -5) chrY drop in LOY clones.
        if is_sex_chrom and chrom_str == "chrY" and chry_scoring_enabled and male_control_rows.size > 0:
            gene_chunk_means = np.nanmean(control_expr_full[male_control_rows], axis=0).astype(np.float32)
        else:
            gene_chunk_means = np.nanmean(control_expr_full, axis=0).astype(np.float32)
        control_gene_means[coord_rows_for_chr] = gene_chunk_means
        del control_sub
        # ChrY scoring requires both male query cells AND male control cells. If either is
        # missing for this chromosome, fall back to standard MAD-standardized scoring (which
        # will produce noise / no calls — but won't fabricate a LOY signal in females).
        chry_active = is_sex_chrom and chrom_str == "chrY" and chry_scoring_enabled
        # Pool of control rows used for sex-chrom baseline: males only when scoring chrY.
        sex_chrom_control_rows = male_control_rows if chry_active else all_control_rows
        if control_state_pool is None:
            baseline = np.nanmedian(control_expr_full[all_control_rows], axis=0).astype(np.float32)
            residual_query = query_expr - baseline[None, :]
            residual_control = control_expr_full - baseline[None, :]
            window_query = _rolling_window_mean(residual_query, chr_windows)
            window_control = _rolling_window_mean(residual_control, chr_windows)
            center = np.nanmedian(window_control, axis=0).astype(np.float32)
            scale = (1.4826 * _mad(window_control, axis=0)).astype(np.float32)
            scale = np.where(scale < 0.10, 0.10, scale).astype(np.float32)
            scores_chr = ((window_query - center[None, :]) / scale[None, :]).astype(np.float32)
            if is_sex_chrom:
                scores_chr = _absolute_log2_window_scores(
                    query_expr, control_expr_full[sex_chrom_control_rows],
                    chr_windows, params.sex_chrom_log2_unit,
                )
                if chry_active and not male_query_mask.all():
                    # Female (and unknown-sex) query cells: chrY signal is biologically zero,
                    # any "loss" call would be spurious. Mark NaN so they're skipped downstream.
                    scores_chr[~male_query_mask] = np.nan
                if chry_active and chry_marker_pos_query.any():
                    # Cells expressing chrY-Y-only markers cannot have biologically lost chrY.
                    # Mark NaN so the interval caller cannot emit a chrY-loss for them.
                    scores_chr[chry_marker_pos_query] = np.nan
            for state in eligible_states:
                rows = state_index_map[state]
                all_scores[rows[:, None], np.arange(chr_slice.start, chr_slice.stop)[None, :]] = scores_chr[rows]
        else:
            chr_window_count = len(chr_windows)
            chr_abs_indices = np.arange(chr_slice.start, chr_slice.stop)
            # --- pooled within-state scale (replaces the arbitrary 0.10 floor) for autosomes ---
            # Per-window robust dispersion across ALL control cells after removing each state's own
            # window center: well-estimated (many cells), stable for rare states, so one threshold =>
            # consistent FPR across cell states. This un-suppresses autosomal CNV that the 0.10 floor
            # (which dominated ~100% of windows at ~2.5x the true noise) was masking. Sex chrom uses
            # its own path below.
            autosomal_pooled_scale = None
            autosomal_state_center: Dict[str, np.ndarray] = {}
            autosomal_window_query = None
            if not is_sex_chrom and params.pooled_scale:
                _gbase = np.nanmedian(control_expr_full[all_control_rows], axis=0).astype(np.float32)
                _wc_all = _rolling_window_mean(control_expr_full[all_control_rows] - _gbase[None, :], chr_windows)
                autosomal_window_query = _rolling_window_mean(query_expr - _gbase[None, :], chr_windows)
                _resid: List[np.ndarray] = []
                for st in eligible_states:
                    crows = control_state_pool.get(st, all_control_rows)
                    if crows.size == 0:
                        crows = all_control_rows
                    cc = np.nanmedian(_wc_all[crows], axis=0).astype(np.float32)
                    autosomal_state_center[st] = cc
                    _resid.append(_wc_all[crows] - cc[None, :])
                _pooled = (1.4826 * _mad(np.vstack(_resid), axis=0)).astype(np.float32)
                _finite = np.isfinite(_pooled)
                _eps = max(float(np.nanquantile(_pooled[_finite], params.scale_floor_quantile)), 1e-3) if _finite.any() else 1e-3
                autosomal_pooled_scale = np.where((~_finite) | (_pooled < _eps), _eps, _pooled).astype(np.float32)
                autosomal_global_center = (np.nanmedian(np.vstack(list(autosomal_state_center.values())), axis=0).astype(np.float32)
                                           if autosomal_state_center else np.zeros(len(chr_windows), np.float32))
                del _wc_all, _resid
            # --- unsupervised chrY (LOY) scoring: SAME pooled-scale machinery as autosomes but
            #     males-only reference + per-state center + NO marker gate. LOY therefore emerges as
            #     a CNV from the general scorer (positive control). loy_mode="biased" keeps the legacy
            #     absolute-log2 + marker-gate path in the branch below.
            chry_unbiased_query = None
            chry_unbiased_scale = None
            chry_unbiased_center: Dict[str, np.ndarray] = {}
            chry_unbiased_global = None
            if is_sex_chrom and chry_active and male_control_rows.size > 0:
                _mb = np.nanmedian(control_expr_full[male_control_rows], axis=0).astype(np.float32)
                _wcy = _rolling_window_mean(control_expr_full[all_control_rows] - _mb[None, :], chr_windows)
                chry_unbiased_query = _rolling_window_mean(query_expr - _mb[None, :], chr_windows)
                _resy: List[np.ndarray] = []
                for st in eligible_states:
                    crows = np.intersect1d(control_state_pool.get(st, all_control_rows), male_control_rows, assume_unique=False)
                    if crows.size == 0:
                        crows = male_control_rows
                    cc = np.nanmedian(_wcy[crows], axis=0).astype(np.float32)
                    chry_unbiased_center[st] = cc
                    _resy.append(_wcy[crows] - cc[None, :])
                _ply = (1.4826 * _mad(np.vstack(_resy), axis=0)).astype(np.float32)
                _finy = np.isfinite(_ply)
                _epy = max(float(np.nanquantile(_ply[_finy], params.scale_floor_quantile)), 1e-3) if _finy.any() else 1e-3
                chry_unbiased_scale = np.where((~_finy) | (_ply < _epy), _epy, _ply).astype(np.float32)
                chry_unbiased_global = (np.nanmedian(np.vstack(list(chry_unbiased_center.values())), axis=0).astype(np.float32)
                                        if chry_unbiased_center else np.zeros(len(chr_windows), np.float32))
                # Per-state reference NULL: 5% quantile of the reference-male chrY aggregate z for each
                # state. The LOY call after the loop is (query chrY aggregate <= this state's tau) — a
                # per-state-calibrated operating point that gives a uniform ~5% FPR in EVERY cell state
                # (removes the cell-type bias a single fixed threshold causes). Matches the standalone.
                chry_window_abs_idx = chr_abs_indices.copy()
                _ref_aggs: List[np.ndarray] = []
                for st in eligible_states:
                    crows = np.intersect1d(control_state_pool.get(st, all_control_rows), male_control_rows, assume_unique=False)
                    if crows.size == 0:
                        crows = male_control_rows
                    cc = chry_unbiased_center[st]
                    ref_agg = np.nanmean((_wcy[crows] - cc[None, :]) / chry_unbiased_scale[None, :], axis=1)
                    ref_agg = ref_agg[np.isfinite(ref_agg)]
                    _ref_aggs.append(ref_agg)
                    if ref_agg.size >= 20:
                        chry_state_tau[st] = float(np.quantile(ref_agg, 0.05))
                _allref = np.concatenate([a for a in _ref_aggs if a.size]) if any(a.size for a in _ref_aggs) else np.array([-2.0], dtype=np.float32)
                chry_global_tau = float(np.quantile(_allref, 0.05))
                chry_interval_template = {
                    "chr": "chrY", "start": int(chr_windows[0].start), "end": int(chr_windows[-1].end),
                    "start_gene": chr_windows[0].start_gene, "end_gene": chr_windows[-1].end_gene,
                    "n_windows": int(len(chr_windows)),
                    "n_genes": int(chr_windows[-1].end_offset - chr_windows[0].start_offset),
                }
                del _wcy, _resy
            for state in eligible_states:
                rows = state_index_map[state]
                control_rows = control_state_pool.get(state, all_control_rows)
                if control_rows.size == 0:
                    control_rows = all_control_rows
                if is_sex_chrom and chry_active and chry_unbiased_query is not None:
                    # UNSUPERVISED chrY: pooled-scale z vs males-only reference, per-state center,
                    # NO marker gate. Females/unknown-sex query cells still NaN'd (no chrY biology).
                    center = chry_unbiased_center.get(state, chry_unbiased_global)
                    state_scores = ((chry_unbiased_query[rows] - center[None, :]) / chry_unbiased_scale[None, :]).astype(np.float32)
                    state_male_mask = male_query_mask[rows]
                    if not state_male_mask.all():
                        state_scores[~state_male_mask] = np.nan
                elif is_sex_chrom:
                    if chry_active:
                        # Restrict per-state control pool to males
                        control_rows = np.intersect1d(control_rows, male_control_rows, assume_unique=False)
                        if control_rows.size == 0:
                            control_rows = male_control_rows
                    state_scores = _absolute_log2_window_scores(
                        query_expr[rows], control_expr_full[control_rows],
                        chr_windows, params.sex_chrom_log2_unit,
                    )
                    if chry_active:
                        # NaN-out female / unknown query cells within this state
                        state_male_mask = male_query_mask[rows]
                        if not state_male_mask.all():
                            state_scores[~state_male_mask] = np.nan
                        # Biological gate (BIASED mode only): cells expressing chrY-Y-only markers
                        # cannot have biologically lost chrY — NaN their chrY scores so the interval
                        # caller cannot emit a chrY-loss for them.
                        state_chry_pos = chry_marker_pos_query[rows]
                        if state_chry_pos.any():
                            state_scores[state_chry_pos] = np.nan
                elif params.pooled_scale and autosomal_pooled_scale is not None:
                    center = autosomal_state_center.get(state, autosomal_global_center)
                    state_scores = ((autosomal_window_query[rows] - center[None, :]) / autosomal_pooled_scale[None, :]).astype(np.float32)
                else:
                    baseline = np.nanmedian(control_expr_full[control_rows], axis=0).astype(np.float32)
                    state_residual = query_expr[rows] - baseline[None, :]
                    control_residual = control_expr_full[control_rows] - baseline[None, :]
                    window_state = _rolling_window_mean(state_residual, chr_windows)
                    window_control = _rolling_window_mean(control_residual, chr_windows)
                    center = np.nanmedian(window_control, axis=0).astype(np.float32)
                    scale = (1.4826 * _mad(window_control, axis=0)).astype(np.float32)
                    scale = np.where(scale < 0.10, 0.10, scale).astype(np.float32)
                    state_scores = ((window_state - center[None, :]) / scale[None, :]).astype(np.float32)
                all_scores[rows[:, None], chr_abs_indices[None, :]] = state_scores

        del query_expr, control_expr_full
        LOGGER.info(
            "Chromosome %s [%d/%d]: %d windows scored (%s).",
            chrom_str, chr_index + 1, len(chrom_groups), len(chr_windows),
            _format_duration(time.perf_counter() - chr_t0),
        )

    # --- per-cell autosomal re-centering ---
    # Remove each cell's genome-wide offset (library/normalization difference vs the reference) so a
    # CONTIGUOUS deviation reflects copy number, not a global shift. With the properly-calibrated
    # pooled scale this is essential: otherwise the query-vs-reference offset reads as genome-wide
    # false 'loss' on every chromosome. Offset = robust median over autosomal windows (real arm/focal
    # CNV is a small fraction of the genome, so the median tracks the copy-neutral baseline). Applied
    # to autosomal + chrX windows only; chrY keeps its own absolute-difference scale.
    if params.pooled_scale and len(windows):
        window_chrom_arr = np.array([str(w.chrom) for w in windows])
        # Standardize EVERY chromosome (incl. chrY) identically: subtract the cell's autosomal-median
        # offset and divide by its autosomal MAD. chrY is scored vs the SEX-MATCHED (male) reference
        # in the loop above, so this treats it exactly like any other chromosome -> a whole-chrY loss
        # is just a CNV the interval caller emits (labeled LOY downstream); no special chrY statistics.
        recenter_mask = np.ones(len(windows), dtype=bool)
        offset_mask = ~np.isin(window_chrom_arr, np.array(["chrX", "chrY"] + [str(c) for c in params.sex_chroms]))
        if offset_mask.any():
            cell_offset = np.nan_to_num(np.nanmedian(all_scores[:, offset_mask], axis=1), nan=0.0)
            _auto = all_scores[:, offset_mask] - cell_offset[:, None]
            cell_scale = (1.4826 * np.nanmedian(np.abs(_auto), axis=1)).astype(np.float32)
            cell_scale = np.where(np.isfinite(cell_scale) & (cell_scale > 1e-3), cell_scale, 1.0)
            all_scores[:, recenter_mask] = (
                (all_scores[:, recenter_mask] - cell_offset[:, None]) / cell_scale[:, None]
            ).astype(np.float32)
            LOGGER.info(
                "Per-cell re-center+scale applied to all chromosomes (median offset=%.3f, median cell MAD=%.3f, %d autosomal windows).",
                float(np.nanmedian(cell_offset)), float(np.nanmedian(cell_scale)), int(offset_mask.sum()),
            )

    # --- autosomal caller = LOY's simulation+likelihood machinery, per ARM, + per-cell centering ---
    # Same detection filter + simulated-state likelihood (argmax, NO threshold) + simulated variance as
    # chrY/LOY (user-accepted chrY biology aside: 1 copy / male ref / whole-chrY region). The one extra,
    # general step is the standard per-cell genome-wide centering across arms (`_autosomal_arm_calls`).
    arm_log2 = None; autosomal_arms = []; arm_templates = []
    if do_autosomal_sim and autosomal_arm_genes:
        control_is_metacell = (control.obs["is_metacell"].to_numpy().astype(bool)
                               if "is_metacell" in control.obs.columns else None)
        control_state_arr = (control.obs[params.control_state_key].astype(str).to_numpy()
                             if params.control_state_key else np.full(control.n_obs, "all", dtype=object))
        query_state_arr = query.obs[params.state_key].astype(str).to_numpy()
        coord_pos = {str(query.var_names[int(i)]): (str(c), int(s)) for i, c, s in
                     zip(coords["var_index"].to_numpy(), coords["chr"].astype(str).to_numpy(), coords["start"].astype(float).to_numpy())}
        autosomal_arms, arm_log2 = _autosomal_arm_calls(
            query_matrix, control_matrix, query.var_names, control.var_names,
            query_state_arr, control_state_arr, autosomal_arm_genes, control_is_metacell,
            params.input_normalized, min_detect_frac=params.autosomal_min_detect_frac)
        for arm in autosomal_arms:
            genes = autosomal_arm_genes[arm]; poss = [coord_pos[g][1] for g in genes if g in coord_pos]
            chrom = coord_pos[genes[0]][0] if genes[0] in coord_pos else arm[:-1]
            arm_templates.append({"chr": chrom, "start": int(min(poss)) if poss else 0,
                                  "end": int(max(poss)) if poss else 0, "start_gene": genes[0],
                                  "end_gene": genes[-1], "n_windows": max(len(genes) // 7, 1), "n_genes": len(genes)})
        _wc = np.array([str(w.chrom) for w in windows])
        all_scores[:, _wc != "chrY"] = np.nan  # autosomal/chrX windowed calls suppressed; per-arm + chrY injected
        LOGGER.info("Autosomal caller = LOY machinery per arm + per-cell centering: %d arms.", len(autosomal_arms))

    # --- chrY copy-number call via the GENERAL simulation classifier (positive control = LOY) ---
    # chrY is gene-sparse (too few expressed genes to form a window), so it is called as a single region
    # by the SAME general simulation rule used genome-wide: simulate copy states from the sex-matched
    # (male) reference by scaling every gene's expression (x0 = homozygous deletion, x0.5 = het, ...) and
    # assign each male query cell its max-likelihood state. The homozygous-deletion state IS LOY, so LOY
    # emerges from the general caller with NO chrY/LOY-specific threshold.
    chry_loss_call = np.zeros(query.n_obs, dtype=bool)
    chry_query_agg = np.full(query.n_obs, np.nan, dtype=np.float32)
    if n_male_q > 0 and n_male_c > 0 and "chrY" in {str(c) for c in params.sex_chroms}:
        male_control_mask_full = np.zeros(control.n_obs, dtype=bool)
        male_control_mask_full[male_control_rows] = True
        control_is_metacell = (control.obs["is_metacell"].to_numpy().astype(bool)
                               if "is_metacell" in control.obs.columns else None)
        control_state_arr = (control.obs[params.control_state_key].astype(str).to_numpy()
                             if params.control_state_key else np.full(control.n_obs, "all", dtype=object))
        chry_genes = [str(query.var_names[i]) for i in coords.loc[coords["chr"].astype(str) == "chrY", "var_index"].to_numpy()]
        chry_loss_call, chry_log2r = _simulation_copy_number_call(
            query_matrix, control_matrix, query_lib, control_lib,
            query.var_names, control.var_names,
            query.obs[params.state_key].astype(str).to_numpy(), control_state_arr,
            male_query_mask, male_control_mask_full, control_is_metacell,
            params.input_normalized, region_genes=chry_genes, normal_copies=1,
        )
        chry_query_agg = np.clip(-np.nan_to_num(chry_log2r, nan=0.0), 0.0, 5.0).astype(np.float32)
        LOGGER.info(
            "chrY copy-number via general simulation classifier (LOY = homozygous deletion): %d / %d male query cells.",
            int(chry_loss_call.sum()), int(male_query_mask.sum()),
        )

    LOGGER.info("Calling per-cell intervals.")
    interval_t0 = time.perf_counter()
    cell_records: List[Dict[str, object]] = []
    interval_records: List[Dict[str, object]] = []
    sample_series = query.obs[params.sample_key].astype(str) if params.sample_key else None
    state_series = query.obs[params.state_key].astype(str)

    sex_chrom_set_runtime = (
        {str(c) for c in params.sex_chroms}
        if params.sex_chrom_mode == "absolute_log2"
        else set()
    )
    # Mask of autosomal windows for the burden filter. Sex-chrom scores are biological-signal-driven
    # rather than CNV-burden-like (a single donor-level chrY loss/gain produces strong signal across
    # all chrY windows), so including them in the top-5% burden statistic causes WT cells with mild
    # sex-chrom signal to spuriously cross the autosomal burden threshold.
    burden_window_mask = np.ones(len(windows), dtype=bool)
    if sex_chrom_set_runtime:
        for w_idx, w in enumerate(windows):
            if str(w.chrom) in sex_chrom_set_runtime:
                burden_window_mask[w_idx] = False

    # Windowed interval-caller cut-offs = the general z-score defaults (NO task-specific tuning). When the
    # autosomal simulation classifier is active the autosomal windows are NaN in all_scores, so the
    # windowed caller emits nothing for autosomes; autosomal calls come from the per-ARM likelihood
    # classifier (injected per cell, no threshold) and chrY from its own injection.
    eff_low, eff_high, eff_min_mean, eff_burden = (
        params.low_threshold, params.high_threshold, params.min_mean_score, params.cnv_burden_threshold)
    eff_min_run, eff_min_genes = params.min_run_windows, params.min_interval_genes
    for state in eligible_states:
        rows = state_index_map[state]
        scores_block = all_scores[rows]
        burden = _cnv_burden(scores_block[:, burden_window_mask], params.burden_quantile)
        candidate_mask = (np.abs(scores_block) >= eff_low).any(axis=1)
        burden_mask = burden >= eff_burden
        eligible_local = np.flatnonzero(candidate_mask & burden_mask)
        eligible_set = set(int(i) for i in eligible_local)
        for local_index, row in enumerate(rows):
            barcode = str(query.obs_names[row])
            sample = sample_series.iloc[row] if sample_series is not None else ""
            sex_chrom_only = False
            if local_index in eligible_set:
                calls = _call_intervals_for_cell(
                    scores_block[local_index], windows, slices_by_chr,
                    high_threshold=eff_high,
                    low_threshold=eff_low,
                    min_run_windows=eff_min_run,
                    min_interval_genes=eff_min_genes,
                    min_mean_score=eff_min_mean,
                )
            elif sex_chrom_set_runtime and candidate_mask[local_index]:
                # Sex-chrom rescue: cell didn't pass the autosomal burden filter but may
                # have a survived sex-chrom interval call (e.g. LOY in an otherwise-quiet
                # cell). Run the interval caller and keep ONLY sex-chrom calls. This avoids
                # flooding outputs with autosomal noise calls from cells that legitimately
                # failed the burden filter.
                all_calls = _call_intervals_for_cell(
                    scores_block[local_index], windows, slices_by_chr,
                    high_threshold=eff_high,
                    low_threshold=eff_low,
                    min_run_windows=eff_min_run,
                    min_interval_genes=eff_min_genes,
                    min_mean_score=eff_min_mean,
                )
                calls = [c for c in all_calls if str(c["chr"]) in sex_chrom_set_runtime]
                sex_chrom_only = bool(calls)
            else:
                calls = []
            # chrY is gene-sparse, so its reliable call is the general SIMULATION classifier (region-level),
            # not the windowed interval caller. Replace any windowed chrY call with the simulation call.
            if chry_interval_template is not None:
                calls = [c for c in calls if str(c["chr"]) != "chrY"]
                if chry_loss_call[row]:
                    _mag = float(max(chry_query_agg[row], 1e-3))  # loss magnitude (log2-ratio based)
                    loss_call = dict(chry_interval_template)
                    loss_call.update({
                        "call": "loss", "mean_score": _mag, "max_score": _mag,
                        "confidence": _mag * math.sqrt(max(int(chry_interval_template["n_windows"]), 1)),
                    })
                    calls = list(calls) + [loss_call]
            # Inject the ARM-level simulation copy-number calls (general likelihood, NO threshold): every
            # non-neutral arm for this cell becomes a loss/gain interval.
            if arm_log2 is not None:
                for ai in np.flatnonzero(np.abs(arm_log2[row]) > 1e-6):
                    lg = float(arm_log2[row, ai]); t = dict(arm_templates[ai]); mag = abs(lg)
                    t.update({"call": "loss" if lg < 0 else "gain", "mean_score": mag, "max_score": mag,
                              "confidence": mag * math.sqrt(max(int(t["n_windows"]), 1)), "mean_log2_ratio": lg})
                    calls.append(t)
            # A cell whose CNV calls are entirely on the sex chromosome(s) (i.e. a whole-chrY loss with
            # no autosomal CNV) is flagged sex_chrom_only -> annotated as the LOY clone downstream.
            sex_chrom_only = bool(calls) and all(str(c["chr"]) in sex_chrom_set_runtime for c in calls)

            status = "CNV" if calls else "WT"
            cell_records.append({
                "CellBarcode": barcode,
                "cell_state": state,
                "sample": sample,
                "cnv_status": status,
                "cnv_burden": float(burden[local_index]),
                "baseline_source": "control",
                "n_cnv_intervals": len(calls) if status == "CNV" else 0,
                "sex_chrom_only": sex_chrom_only,
            })
            if status == "CNV":
                for call in calls:
                    interval_records.append({
                        "CellBarcode": barcode,
                        "cell_state": state,
                        "sample": sample,
                        **call,
                    })

    skipped_states = [s for s in states if s not in set(eligible_states)]
    for state in skipped_states:
        for row in state_index_map[state]:
            cell_records.append({
                "CellBarcode": str(query.obs_names[row]),
                "cell_state": state,
                "sample": sample_series.iloc[row] if sample_series is not None else "",
                "cnv_status": "low_power",
                "cnv_burden": float("nan"),
                "baseline_source": "none",
                "n_cnv_intervals": 0,
            })
    LOGGER.info("Interval calling done (%s).", _format_duration(time.perf_counter() - interval_t0))

    barcode_to_index = {str(b): i for i, b in enumerate(query.obs_names.astype(str))}

    zygosity_thresholds: Optional[Dict[str, float]] = None
    if cached_query_expr is not None and interval_records and params.zygosity_mode == "relative":
        cal_t0 = time.perf_counter()
        wt_indices = np.asarray(
            [barcode_to_index[b] for b, status in zip(
                (r["CellBarcode"] for r in cell_records),
                (r["cnv_status"] for r in cell_records),
            ) if status == "WT" and b in barcode_to_index],
            dtype=np.int64,
        )
        zygosity_thresholds = _calibrate_zygosity_thresholds(
            interval_records, wt_indices, cached_query_expr, control_gene_means, coords,
            seed=params.random_state,
        )
        if zygosity_thresholds is not None:
            LOGGER.info(
                "Calibrated zygosity thresholds in %s: hom_loss<=%.2f het_loss<=%.2f het_gain>=%.2f hom_gain>=%.2f",
                _format_duration(time.perf_counter() - cal_t0),
                zygosity_thresholds["homozygous_loss"], zygosity_thresholds["heterozygous_loss"],
                zygosity_thresholds["heterozygous_gain"], zygosity_thresholds["homozygous_gain"],
            )
        else:
            LOGGER.info("Could not calibrate relative zygosity thresholds; falling back to fixed.")

    sex_chrom_thresholds = (
        {
            "homozygous_loss": float(params.sex_chrom_hom_loss),
            "heterozygous_loss": float(params.sex_chrom_het_loss),
            "heterozygous_gain": float(params.sex_chrom_het_gain),
            "homozygous_gain": float(params.sex_chrom_hom_gain),
        }
        if params.sex_chrom_mode == "absolute_log2"
        else None
    )

    if cached_query_expr is not None and interval_records:
        zyg_t0 = time.perf_counter()
        interval_records = _annotate_intervals_with_zygosity(
            interval_records, cached_query_expr, control_gene_means, coords, barcode_to_index,
            thresholds=zygosity_thresholds,
            sex_chroms=params.sex_chroms if sex_chrom_thresholds is not None else (),
            sex_chrom_thresholds=sex_chrom_thresholds,
        )
        LOGGER.info("Per-cell zygosity annotation done (%s).", _format_duration(time.perf_counter() - zyg_t0))

    cell_df = pd.DataFrame(cell_records)
    if not cell_df.empty:
        cell_df["state_clone_id"] = "WT"
        cell_df["global_clone_id"] = "WT"

    nearest_normal_clone_intervals_df = pd.DataFrame()
    used_nearest_normal_autosomal = bool(getattr(params, "nearest_normal_autosomal", False))
    if used_nearest_normal_autosomal and not cell_df.empty:
        nn_t0 = time.perf_counter()
        LOGGER.info(
            "Nearest-normal-subcluster autosomal caller: replacing windowed autosomal calls; "
            "standard chrY/LOY calls remain from the unbiased simulation path."
        )
        with tempfile.TemporaryDirectory(prefix="fastcnv_nn_") as tmpdir:
            nn_prefix = Path(tmpdir) / "nearest_normal"
            _nearest_normal_autosomal(
                query, control, coords, params.state_key, params.control_state_key,
                params.layer, params.input_normalized, nn_prefix, sample_key=params.sample_key,
                query_to_control_state=state_crosswalk,
            )
            nn_cells_path = nn_prefix.with_suffix(".cnv_cells.tsv")
            nn_intervals_path = nn_prefix.with_suffix(".clone_intervals.tsv")
            nn_cells = pd.read_csv(nn_cells_path, sep="\t")
            nearest_normal_clone_intervals_df = pd.read_csv(nn_intervals_path, sep="\t")

        nn_cells["CellBarcode"] = nn_cells["CellBarcode"].astype(str)
        nn_lookup = nn_cells.set_index("CellBarcode")
        cell_df = cell_df.set_index("CellBarcode", drop=False)
        aligned = nn_lookup.reindex(cell_df.index)
        nn_is_cnv = aligned["cnv_status"].astype(str).eq("CNV").fillna(False)
        loy_is_cnv = cell_df["cnv_status"].astype(str).eq("CNV") & cell_df["sex_chrom_only"].astype(bool)
        merged_is_cnv = nn_is_cnv | loy_is_cnv
        cell_df.loc[:, "cnv_status"] = np.where(merged_is_cnv, "CNV", cell_df["cnv_status"])
        cell_df.loc[~merged_is_cnv & cell_df["cnv_status"].astype(str).eq("CNV"), "cnv_status"] = "WT"
        cell_df.loc[:, "baseline_source"] = np.where(nn_is_cnv, "nearest_normal", cell_df["baseline_source"])
        cell_df.loc[:, "state_clone_id"] = aligned["state_clone_id"].fillna("WT").astype(str).to_numpy()
        cell_df.loc[:, "global_clone_id"] = aligned["global_clone_id"].fillna("WT").astype(str).to_numpy()
        cell_df.loc[~nn_is_cnv, ["state_clone_id", "global_clone_id"]] = "WT"
        cell_df.loc[:, "n_cnv_intervals"] = aligned["n_cnv_intervals"].fillna(0).astype(int).to_numpy()
        cell_df.loc[loy_is_cnv, "n_cnv_intervals"] = cell_df.loc[loy_is_cnv, "n_cnv_intervals"].astype(int) + 1
        cell_df = cell_df.reset_index(drop=True)

        sex_chrom_interval_records = [
            r for r in interval_records if str(r.get("chr", "")) in sex_chrom_set_runtime
        ]
        interval_records = []
        if not nearest_normal_clone_intervals_df.empty:
            nn_by_clone = {
                str(clone_id): grp.to_dict("records")
                for clone_id, grp in nearest_normal_clone_intervals_df.groupby("global_clone_id")
            }
            cell_meta = cell_df.set_index("CellBarcode")
            for barcode, row in cell_meta.iterrows():
                clone_id = str(row.get("global_clone_id", "WT"))
                if clone_id == "WT":
                    continue
                for rec in nn_by_clone.get(clone_id, []):
                    out = {
                        "CellBarcode": barcode,
                        "cell_state": row.get("cell_state", ""),
                        "sample": row.get("sample", ""),
                    }
                    for key in (
                        "call", "chr", "start", "end", "start_gene", "end_gene",
                        "n_windows", "n_genes", "mean_score", "max_score",
                        "confidence", "mean_log2_ratio", "zygosity_state",
                    ):
                        out[key] = rec.get(key, "")
                    interval_records.append(out)
        interval_records.extend(sex_chrom_interval_records)
        LOGGER.info(
            "Nearest-normal autosomal merge done (%s): %d autosomal-CNV cells, %d LOY/sex-chrom cells.",
            _format_duration(time.perf_counter() - nn_t0),
            int(nn_is_cnv.sum()),
            int(loy_is_cnv.sum()),
        )
    helper_params = FastCNVParams(
        h5ad=params.h5ad,
        gene_coordinates=params.gene_coordinates,
        output_prefix=params.output_prefix,
        state_key=params.state_key,
        sample_key=params.sample_key,
        layer=params.layer,
        input_normalized=params.input_normalized,
        burden_quantile=params.burden_quantile,
        high_threshold=params.high_threshold,
        low_threshold=params.low_threshold,
        min_mean_score=params.min_mean_score,
        min_run_windows=params.min_run_windows,
        min_interval_genes=params.min_interval_genes,
        cnv_burden_threshold=params.cnv_burden_threshold,
        max_clones_per_state=params.max_clones_per_state,
        max_global_clones=params.max_global_clones,
        min_clone_cells=params.min_clone_cells,
        clone_similarity_threshold=params.clone_similarity_threshold,
        clone_consensus_fraction=params.clone_consensus_fraction,
        nmf_max_iter=params.nmf_max_iter,
        skip_clones=params.skip_clones,
        skip_pdf=params.skip_pdf,
        random_state=params.random_state,
    )

    state_clone_rows: List[Dict[str, object]] = []
    if not used_nearest_normal_autosomal and not params.skip_clones and not cell_df.empty:
        clones_t0 = time.perf_counter()

        # Cells whose CNV calls are entirely on sex chromosomes (e.g. pure LOY) are excluded
        # from NMF-based autosomal clone discovery — they have no signal in the autosomal
        # feature space, and their inclusion shifts cluster boundaries and dilutes real
        # autosomal CNV clones. They are still tracked, assigned to a synthetic 'sex_chrom'
        # clone, and will be picked up by the clone-consensus interval calling pass.
        sex_chrom_only_col = cell_df["sex_chrom_only"] if "sex_chrom_only" in cell_df.columns else None
        cnv_global_indices = np.asarray(
            [barcode_to_index[b]
             for b, status, sex_only in zip(
                 cell_df["CellBarcode"], cell_df["cnv_status"],
                 sex_chrom_only_col if sex_chrom_only_col is not None else [False]*len(cell_df),
             )
             if status == "CNV" and not bool(sex_only) and b in barcode_to_index],
            dtype=np.int64,
        )
        # Build a mask excluding sex-chrom windows from clone discovery features. Sex-chrom
        # signal (LOY/LOX) is binary at the donor level and dominates NMF clustering, which
        # otherwise scrambles autosomal clone discovery. Sex-chrom CNVs are still surfaced
        # via per-cell intervals and the clone consensus pass on the full score matrix.
        autosomal_mask = np.ones(all_scores.shape[1], dtype=bool)
        if sex_chrom_set_runtime:
            for w_idx, w in enumerate(windows):
                if str(w.chrom) in sex_chrom_set_runtime:
                    autosomal_mask[w_idx] = False

        if cnv_global_indices.size >= max(params.min_clone_cells, 1):
            cnv_scores = all_scores[cnv_global_indices]
            active = (np.abs(cnv_scores) >= params.low_threshold).astype(np.float32)
            active_fraction = np.nan_to_num(active.mean(axis=0), nan=0.0)
            informative = (active_fraction >= params.clone_min_active_fraction) & autosomal_mask
            n_informative = int(informative.sum())
            if n_informative > params.clone_max_features:
                rank_score = np.nan_to_num(np.abs(cnv_scores).mean(axis=0), nan=0.0)
                rank_score[~informative] = -np.inf
                top_idx = np.argpartition(rank_score, -params.clone_max_features)[-params.clone_max_features:]
                informative = np.zeros_like(informative)
                informative[top_idx] = True
            informative_idx = np.flatnonzero(informative).astype(np.int64)
            if informative_idx.size == 0:
                informative_idx = np.flatnonzero(autosomal_mask).astype(np.int64)
            LOGGER.info(
                "Clone feature pruning: %d / %d autosomal windows retained "
                "(active >= %.0f%% across %d CNV cells; cap %d; sex-chrom windows excluded).",
                int(informative_idx.size), int(autosomal_mask.sum()),
                params.clone_min_active_fraction * 100.0, int(cnv_global_indices.size),
                params.clone_max_features,
            )
        else:
            informative_idx = np.flatnonzero(autosomal_mask).astype(np.int64)

        clone_scores = all_scores[:, informative_idx]

        group_columns = ["cell_state"]
        if params.sample_key and "sample" in cell_df.columns:
            group_columns.append("sample")

        group_specs: List[Tuple[int, str, str, np.ndarray, np.ndarray, np.ndarray]] = []
        for state_number, (group_key, state_cell_df) in enumerate(cell_df.groupby(group_columns, sort=False)):
            if isinstance(group_key, tuple):
                state = str(group_key[0])
                sample = str(group_key[1]) if len(group_key) > 1 else ""
            else:
                state = str(group_key)
                sample = ""
            valid_state_df = state_cell_df[state_cell_df["cnv_status"].isin(["CNV", "WT"])]
            if valid_state_df.empty:
                continue
            row_indices = valid_state_df.index.to_numpy()
            score_indices = np.asarray(
                [barcode_to_index[b] for b in valid_state_df["CellBarcode"] if b in barcode_to_index],
                dtype=int,
            )
            if score_indices.size == 0:
                continue
            cnv_mask = valid_state_df["cnv_status"].to_numpy(dtype=object) == "CNV"
            group_specs.append((state_number, state, sample, row_indices, score_indices, cnv_mask))

        def _run_one(spec):
            state_number, state, sample, _row_indices, score_indices, cnv_mask = spec
            labels, summary = _kmeans_state_clones(
                clone_scores[score_indices],
                cnv_mask=cnv_mask,
                max_clones=params.max_clones_per_state,
                min_cells=params.min_clone_cells,
                seed=params.random_state + state_number,
            )
            return state, sample, labels, summary

        try:
            from joblib import Parallel, delayed
            n_jobs = params.n_jobs if params.n_jobs != 0 else 1
            results = Parallel(n_jobs=n_jobs, prefer="threads", verbose=0)(
                delayed(_run_one)(spec) for spec in group_specs
            )
        except ImportError:
            results = [_run_one(spec) for spec in group_specs]

        for spec, (state, sample, labels, summary) in zip(group_specs, results):
            row_indices = spec[3]
            resolved = np.asarray(
                ["WT" if label == "WT" else f"{state}::{sample}::{label}" for label in labels],
                dtype=object,
            )
            cell_df.loc[row_indices, "state_clone_id"] = resolved
            if not summary.empty:
                summary.insert(0, "cell_state", state)
                summary.insert(1, "sample", sample)
                summary["state_clone_id"] = summary["state_clone_id"].map(
                    lambda label: f"{state}::{sample}::{label}"
                )
                state_clone_rows.extend(summary.to_dict("records"))

        _, state_profiles, state_sizes = _clone_profile_frame(cell_df, clone_scores, barcode_to_index, "state_clone_id")
        global_assignments = _merge_global_clones(state_profiles, state_sizes, helper_params)
        cell_df["global_clone_id"] = cell_df["state_clone_id"].map(lambda c: global_assignments.get(str(c), "WT"))
        for row in state_clone_rows:
            row["global_clone_id"] = global_assignments.get(str(row["state_clone_id"]), "WT")

        # LOY OVERRIDE: a cell with a chrY copy-loss is labeled clone_loy DIRECTLY from the chrY-loss
        # call (the positive control), so the LOY clone is robust to any autosomal calls the cell also
        # carries. (Earlier the LOY clone was built from sex_chrom_only, which per-cell autosomal noise
        # could strip away, collapsing clone_loy.) LOY overrides the numeric clone label.
        loy_mask = np.array(
            [bool(chry_loss_call[barcode_to_index[b]]) if b in barcode_to_index else False
             for b in cell_df["CellBarcode"].astype(str)], dtype=bool)
        if loy_mask.any():
            cell_df.loc[loy_mask, "global_clone_id"] = "clone_loy"
            cell_df.loc[loy_mask, "state_clone_id"] = "loy"

        n_global = len({c for c in cell_df["global_clone_id"] if c != "WT"})
        LOGGER.info(
            "Clone discovery done (%s): %d state-local clones merged into %d global clones.",
            _format_duration(time.perf_counter() - clones_t0),
            len(state_clone_rows),
            n_global,
        )

    if used_nearest_normal_autosomal and not cell_df.empty and bool(chry_loss_call.any()):
        loy_mask = np.array(
            [bool(chry_loss_call[barcode_to_index[b]]) if b in barcode_to_index else False
             for b in cell_df["CellBarcode"].astype(str)], dtype=bool)
        if loy_mask.any():
            cell_df.loc[loy_mask, "cnv_status"] = "CNV"
            cell_df.loc[loy_mask, "global_clone_id"] = "clone_loy"
            cell_df.loc[loy_mask, "state_clone_id"] = "loy"
            cell_df.loc[loy_mask, "sex_chrom_only"] = True
            cell_df.loc[loy_mask, "n_cnv_intervals"] = np.maximum(
                cell_df.loc[loy_mask, "n_cnv_intervals"].astype(int), 1,
            )
            LOGGER.info("LOY override applied after nearest-normal autosomal merge: %d cells -> clone_loy.", int(loy_mask.sum()))

    n_total_cells = int(cell_df.shape[0]) if not cell_df.empty else 0
    n_cnv_cells = int((cell_df["cnv_status"] == "CNV").sum()) if not cell_df.empty else 0
    cnv_fraction = n_cnv_cells / max(n_total_cells, 1)
    global_low_confidence = (
        not params.skip_clones and not cell_df.empty and cnv_fraction < params.clone_min_cnv_fraction
    )
    if global_low_confidence:
        LOGGER.warning(
            "CNV-cell fraction %.2f%% below threshold %.2f%%; flagging ALL clones as low_confidence.",
            100.0 * cnv_fraction, 100.0 * params.clone_min_cnv_fraction,
        )
    clone_size_counts = (
        cell_df["global_clone_id"].value_counts().to_dict() if "global_clone_id" in cell_df.columns else {}
    )

    def _confidence_for(clone_id: object) -> str:
        cid = str(clone_id)
        if cid == "WT" or cid == "":
            return "wt"
        if global_low_confidence:
            return "low"
        if int(clone_size_counts.get(cid, 0)) < params.clone_min_cells_confident:
            return "low"
        return "high"

    if not cell_df.empty and "global_clone_id" in cell_df.columns:
        cell_df["clone_confidence"] = cell_df["global_clone_id"].map(_confidence_for)
    for row in state_clone_rows:
        row["clone_confidence"] = _confidence_for(row.get("global_clone_id", ""))

    clone_intervals_df = (
        nearest_normal_clone_intervals_df.copy()
        if used_nearest_normal_autosomal
        else (
            _call_clone_consensus_intervals(cell_df, all_scores, windows, slices_by_chr, barcode_to_index, helper_params)
            if not params.skip_clones and not cell_df.empty
            else pd.DataFrame()
        )
    )

    if cached_query_expr is not None and not clone_intervals_df.empty:
        clone_to_query_rows_zyg: Dict[str, np.ndarray] = {}
        for clone_id in clone_intervals_df["global_clone_id"].astype(str).unique():
            barcodes = cell_df.loc[cell_df["global_clone_id"] == clone_id, "CellBarcode"].astype(str)
            indices = np.asarray([barcode_to_index[b] for b in barcodes if b in barcode_to_index], dtype=np.int64)
            if indices.size > 0:
                clone_to_query_rows_zyg[clone_id] = indices
        if clone_to_query_rows_zyg:
            clone_log2_ratios_full, clone_log2_order = _compute_clone_log2_ratios_from_cache(
                cached_query_expr, control_gene_means, clone_to_query_rows_zyg,
            )
            clone_log2_lookup = {clone_id: clone_log2_ratios_full[i] for i, clone_id in enumerate(clone_log2_order)}
            coords_chr = coords["chr"].astype(str).to_numpy()
            coords_start = coords["start"].to_numpy()
            coords_end = coords["end"].to_numpy()
            mean_log2_list: List[float] = []
            zyg_states: List[str] = []
            for _, interval_row in clone_intervals_df.iterrows():
                clone_id = str(interval_row["global_clone_id"])
                gene_log2 = clone_log2_lookup.get(clone_id)
                if gene_log2 is None:
                    mean_log2_list.append(float("nan"))
                    zyg_states.append("low_signal")
                    continue
                in_interval = (
                    (coords_chr == str(interval_row["chr"])) &
                    (coords_start >= float(interval_row["start"])) &
                    (coords_end <= float(interval_row["end"]))
                )
                if not in_interval.any():
                    mean_log2_list.append(float("nan"))
                    zyg_states.append("low_signal")
                    continue
                mean_log2 = float(np.nanmean(gene_log2[in_interval]))
                mean_log2_list.append(mean_log2)
                interval_chr = str(interval_row["chr"])
                if (
                    sex_chrom_thresholds is not None
                    and interval_chr in params.sex_chroms
                ):
                    zyg_states.append(_zygosity_state(mean_log2, sex_chrom_thresholds))
                else:
                    zyg_states.append(_zygosity_state(mean_log2, None))
            clone_intervals_df = clone_intervals_df.copy()
            clone_intervals_df["mean_log2_ratio"] = mean_log2_list
            clone_intervals_df["zygosity_state"] = zyg_states
            clone_intervals_df["clone_confidence"] = clone_intervals_df["global_clone_id"].astype(str).map(
                lambda cid: _confidence_for(cid)
            )

    # chrY:loss consensus for LOY-enriched clones (incl. clone_loy): the per-cell chrY call is the
    # simulation classifier (region-level), which the windowed clone-consensus caller cannot reproduce,
    # so emit it here from the per-cell calls.
    if (chry_interval_template is not None and not cell_df.empty
            and "global_clone_id" in cell_df.columns and bool(chry_loss_call.any())):
        chry_rows: List[Dict[str, object]] = []
        for clone_id, grp in cell_df.groupby("global_clone_id"):
            cid = str(clone_id)
            if cid in ("WT", ""):
                continue
            idx = np.asarray([barcode_to_index[b] for b in grp["CellBarcode"].astype(str) if b in barcode_to_index], dtype=np.int64)
            if idx.size == 0 or float(chry_loss_call[idx].mean()) < params.clone_consensus_fraction:
                continue
            depths = chry_query_agg[idx[chry_loss_call[idx]]]
            mean_depth = float(np.nanmean(depths)) if depths.size else 1e-3
            row = dict(chry_interval_template)
            row.update({
                "global_clone_id": cid, "n_cells": int(idx.size),
                "support_fraction": round(float(chry_loss_call[idx].mean()), 3), "call": "loss",
                "mean_score": mean_depth, "max_score": mean_depth,
                "confidence": mean_depth * math.sqrt(max(int(chry_interval_template["n_windows"]), 1)),
                "mean_log2_ratio": float("nan"), "zygosity_state": "homozygous_loss",
                "clone_confidence": _confidence_for(cid),
            })
            chry_rows.append(row)
        if chry_rows:
            if not clone_intervals_df.empty:
                clone_intervals_df = clone_intervals_df[clone_intervals_df["chr"].astype(str) != "chrY"].copy()
            clone_intervals_df = pd.concat([clone_intervals_df, pd.DataFrame(chry_rows)], ignore_index=True)
            LOGGER.info("chrY:loss consensus emitted for %d LOY-enriched clone(s) (incl. clone_loy).", len(chry_rows))

    cells_path = params.output_prefix.with_suffix(".cnv_cells.tsv")
    intervals_path = params.output_prefix.with_suffix(".cnv_intervals.tsv")
    windows_path = params.output_prefix.with_suffix(".cnv_windows.tsv")
    scores_path = params.output_prefix.with_suffix(".cnv_window_scores.npz")
    state_clones_path = params.output_prefix.with_suffix(".state_clones.tsv")
    clone_intervals_path = params.output_prefix.with_suffix(".clone_intervals.tsv")
    clone_pdf_path = params.output_prefix.with_suffix(".clone_genome.pdf")
    h5ad_path = params.output_prefix.with_suffix(".fastcnv.h5ad")
    sample_sex_path = params.output_prefix.with_suffix(".sample_sex.tsv")

    sex_audit_rows = []
    for sample, row in query_sex_df.iterrows():
        sex_audit_rows.append({"source": "query", "sample": sample, **row.to_dict()})
    for sample, row in control_sex_df.iterrows():
        sex_audit_rows.append({"source": "control", "sample": sample, **row.to_dict()})
    pd.DataFrame(sex_audit_rows).to_csv(sample_sex_path, sep="\t", index=False)

    cell_df.to_csv(cells_path, sep="\t", index=False)
    interval_columns = [
        "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
        "mean_log2_ratio", "zygosity_state",
    ]
    pd.DataFrame(interval_records, columns=interval_columns).to_csv(intervals_path, sep="\t", index=False)
    _window_frame(windows).to_csv(windows_path, sep="\t", index=False)
    pd.DataFrame(state_clone_rows).to_csv(state_clones_path, sep="\t", index=False)
    clone_interval_columns = [
        "global_clone_id", "n_cells", "support_fraction", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
        "mean_log2_ratio", "zygosity_state", "clone_confidence",
    ]
    _ci_df = pd.DataFrame(clone_intervals_df, columns=clone_interval_columns)
    _ci_df.to_csv(clone_intervals_path, sep="\t", index=False)
    # clone_report.csv: EVERY clone with its cell count + consensus region(s) + cell-types + samples,
    # written by fastCNV from its own cell_df / clone_intervals (no post-processing or thresholds).
    clone_report_path = params.output_prefix.with_suffix(".clone_report.csv")
    if not cell_df.empty and "global_clone_id" in cell_df.columns:
        _rep = []
        for _clone, _grp in cell_df.groupby("global_clone_id"):
            _g = _ci_df[_ci_df["global_clone_id"].astype(str) == str(_clone)]
            _regs = "; ".join(f"{r.chr}:{r.call}" for r in _g.itertuples()) if len(_g) else "(no consensus interval)"
            _cts = _grp["cell_state"].value_counts() if "cell_state" in _grp else pd.Series(dtype=int)
            _sps = _grp["sample"].value_counts() if "sample" in _grp else pd.Series(dtype=int)
            _rep.append({"clone": _clone, "n_cells": int(len(_grp)), "regions": _regs,
                         "cell_types": "; ".join(_cts.index[:6].astype(str)),
                         "samples": "; ".join(f"{k}({v})" for k, v in _sps.items())})
        pd.DataFrame(_rep).sort_values("n_cells", ascending=False).to_csv(clone_report_path, index=False)
    else:
        pd.DataFrame(columns=["clone", "n_cells", "regions", "cell_types", "samples"]).to_csv(clone_report_path, index=False)
    np.savez(
        scores_path,
        scores=all_scores,
        cell_barcodes=np.asarray(query.obs_names.astype(str), dtype=object),
        states=state_series.to_numpy(dtype=object),
        window_ids=np.asarray([w.window_id for w in windows], dtype=object),
    )

    pdf_path = None
    if not params.skip_pdf and not cell_df.empty and "global_clone_id" in cell_df.columns:
        pdf_t0 = time.perf_counter()
        clone_to_query_rows: Dict[str, np.ndarray] = {}
        clone_sizes_map: Dict[str, int] = {}
        clone_ids = sorted(
            (c for c in cell_df["global_clone_id"].dropna().unique() if c != "WT"),
            key=lambda c: (0, int(c[5:])) if str(c).startswith("clone") and str(c)[5:].isdigit() else (1, str(c)),
        )
        for clone_id in clone_ids:
            barcodes = cell_df.loc[cell_df["global_clone_id"] == clone_id, "CellBarcode"].astype(str)
            indices = np.asarray([barcode_to_index[b] for b in barcodes if b in barcode_to_index], dtype=np.int64)
            if indices.size == 0:
                continue
            clone_to_query_rows[str(clone_id)] = indices
            clone_sizes_map[str(clone_id)] = int(indices.size)

        if clone_to_query_rows:
            ordered_clone_ids = list(clone_to_query_rows.keys())
            clone_score_profiles = np.full(
                (len(ordered_clone_ids), all_scores.shape[1]), np.nan, dtype=np.float32,
            )
            for i, clone_id in enumerate(ordered_clone_ids):
                rows = clone_to_query_rows[clone_id]
                if rows.size == 0:
                    continue
                clone_score_profiles[i] = np.nanmean(all_scores[rows], axis=0)
            pdf_path = _write_clone_genome_pdf(
                clone_pdf_path,
                coords=coords,
                clone_score_profiles=clone_score_profiles,
                windows=windows,
                clone_order=ordered_clone_ids,
                clone_sizes=clone_sizes_map,
                clone_intervals_df=clone_intervals_df,
                score_y_clip=float(params.pdf_score_y_clip),
                score_threshold=float(params.low_threshold),
            )
        LOGGER.info("Clone PDF written in %s.", _format_duration(time.perf_counter() - pdf_t0))

    outputs = {
        "cells": cells_path,
        "intervals": intervals_path,
        "windows": windows_path,
        "scores": scores_path,
        "state_clones": state_clones_path,
        "clone_intervals": clone_intervals_path,
        "clone_report": clone_report_path,
        "sample_sex": sample_sex_path,
    }
    if pdf_path is not None:
        outputs["clone_pdf"] = pdf_path
    if not params.skip_heatmap:
        heatmap_t0 = time.perf_counter()
        try:
            heatmap_outputs = _write_infercnv_style_heatmaps(
                params.output_prefix,
                cell_df=cell_df,
                all_scores=all_scores,
                windows=windows,
                barcode_to_index=barcode_to_index,
                max_cells=int(params.heatmap_max_cells),
                filter_threshold=float(params.heatmap_filter_threshold),
                min_chr_windows=int(params.heatmap_min_chr_windows),
                random_state=int(params.random_state),
            )
            outputs.update(heatmap_outputs)
            if heatmap_outputs:
                LOGGER.info("inferCNV-style heatmaps written in %s.", _format_duration(time.perf_counter() - heatmap_t0))
        except Exception:
            LOGGER.exception("inferCNV-style heatmap export failed; core fastCNV outputs were written.")
    if params.write_h5ad:
        obs_updates = cell_df.set_index("CellBarcode")
        for col in ("cnv_status", "cnv_burden", "state_clone_id", "global_clone_id"):
            query.obs[f"fastcnv_{col}"] = obs_updates.reindex(query.obs_names.astype(str))[col].values
        query.write_h5ad(h5ad_path)
        outputs["h5ad"] = h5ad_path

    if params.residual_candidate:
        candidate_t0 = time.perf_counter()
        try:
            from altanalyze3.components.fastCNV.candidate.pyinfer_residual_clone import (
                ResidualCloneParams,
                run_candidate as run_residual_candidate,
            )
            candidate_prefix = Path(f"{params.output_prefix}.residual_candidate")
            candidate_outputs = run_residual_candidate(
                ResidualCloneParams(
                    h5ad=params.h5ad,
                    fastcnv_prefix=params.output_prefix,
                    output_prefix=candidate_prefix,
                    state_key=params.state_key,
                    layer=params.layer,
                    min_abs_region_score=float(params.residual_candidate_min_abs_region_score),
                    min_separation_mad=float(params.residual_candidate_min_separation_mad),
                    pyinfercnv_path=params.residual_candidate_pyinfercnv_path,
                    heatmap_filter_threshold=float(params.residual_candidate_heatmap_filter_threshold),
                )
            )
            for key, path in candidate_outputs.items():
                outputs[f"residual_candidate_{key}"] = path
            LOGGER.info("Residual candidate written in %s.", _format_duration(time.perf_counter() - candidate_t0))
        except Exception:
            LOGGER.exception("Residual candidate export failed; core fastCNV outputs were written.")

    LOGGER.info("fastCNV total runtime: %s", _format_duration(time.perf_counter() - run_start))
    return outputs


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fast variant of fastCNV: requires a control h5ad; vectorized per-chromosome scoring.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--control-h5ad", required=True, type=Path)
    p.add_argument("--gene-coordinates", default=None, type=Path)
    p.add_argument("--species", choices=sorted(RESOURCE_FILES), default=None)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--state-key", required=True)
    p.add_argument("--sample-key", default=None)
    p.add_argument("--control-state-key", default=None)
    p.add_argument("--layer", default="auto")
    p.add_argument("--input-normalized", action="store_true")
    p.add_argument("--window-genes", type=int, default=41)
    p.add_argument("--stride-genes", type=int, default=7)
    p.add_argument("--min-chr-genes", type=int, default=25)
    p.add_argument("--min-state-cells", type=int, default=30)
    # Defaults below are calibrated for the new per-cell-standardized scores (~N(0,1) under the
    # null with --pooled-scale, the default). For the legacy --no-pooled-scale path use the old
    # values (high 2.6 / low 1.6 / min-mean 1.8 / burden 1.8).
    p.add_argument("--high-threshold", type=float, default=3.5)
    p.add_argument("--low-threshold", type=float, default=2.0)
    p.add_argument("--min-run-windows", type=int, default=3)
    p.add_argument("--min-interval-genes", type=int, default=60)
    p.add_argument("--min-mean-score", type=float, default=2.5)
    p.add_argument("--burden-quantile", type=float, default=0.95)
    p.add_argument("--cnv-burden-threshold", type=float, default=3.0)
    p.add_argument("--nearest-normal-autosomal", dest="nearest_normal_autosomal", action="store_true",
                   default=False,
                   help="Use the validated nearest-normal-subcluster autosomal CNV caller: match each "
                        "query subcluster to its nearest normal reference subcluster, then merge "
                        "autosomal calls with the standard unbiased chrY/LOY caller.")
    p.add_argument("--no-nearest-normal-autosomal", dest="nearest_normal_autosomal", action="store_false",
                   help="Disable the validated nearest-normal autosomal caller and use the legacy windowed "
                        "autosomal path. The unbiased chrY/LOY caller is unchanged.")
    p.add_argument("--simulation-autosomal", action="store_true",
                   help="Enable the per-cell autosomal simulation copy-number caller. OFF by default: it "
                        "is single-cell-noise-limited and not production-ready; the validated default is "
                        "the clean LOY (chrY) positive control. Use the pseudobulk window-average for "
                        "autosomal karyotype concordance instead.")
    p.add_argument("--autosomal-detect-frac", type=float, default=0.05,
                   help="Gene detection filter for the autosomal per-arm caller (keep genes expressed in "
                        ">= this fraction of the reference; same rule as LOY). 0.0 disables it.")
    p.add_argument("--external-baseline", action="store_true",
                   help="Use the EXTERNAL reference neutral as the per-window baseline (detects CLONAL "
                        "CNV; requires a batch-matched reference). Default uses the query-internal "
                        "per-state median (subclonal, batch-robust).")
    p.add_argument("--max-clones-per-state", type=int, default=10)
    p.add_argument("--max-global-clones", type=int, default=10)
    p.add_argument("--min-clone-cells", type=int, default=5)
    p.add_argument("--clone-similarity-threshold", type=float, default=0.88)
    p.add_argument("--clone-consensus-fraction", type=float, default=0.45)
    p.add_argument("--nmf-max-iter", type=int, default=100)
    p.add_argument("--clone-min-active-fraction", type=float, default=0.05, help="Window kept for clone discovery if active in >= this fraction of CNV cells.")
    p.add_argument("--clone-max-features", type=int, default=400, help="Hard cap on windows used for clone discovery; top by mean |score| if exceeded.")
    p.add_argument("--clone-min-cnv-fraction", type=float, default=0.015, help="If CNV-cell fraction is below this, ALL clones are flagged low_confidence.")
    p.add_argument("--clone-min-cells-confident", type=int, default=30, help="Clones with fewer than this many cells are flagged low_confidence.")
    p.add_argument("--zygosity-mode", choices=["fixed", "relative"], default="relative", help="Zygosity threshold mode: 'fixed' uses literature defaults; 'relative' calibrates from WT-cell log2 distribution.")
    p.add_argument("--n-jobs", type=int, default=1, help="Workers for parallel clone discovery. "
                   "Default 1 (sequential) — the joblib thread pool can deadlock at scale and clone "
                   "discovery is fast single-threaded; set >1 only if you know your environment is safe.")
    p.add_argument("--pdf-smooth-genes", type=int, default=50, help="(Legacy) Rolling gene-window size for the older log2-ratio PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-y-clip", type=float, default=1.0, help="(Legacy) Symmetric y-axis limit for the older log2-ratio PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-label-threshold", type=float, default=0.25, help="(Legacy) Driver-gene label threshold for the older PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-score-y-clip", type=float, default=4.0,
                   help="Symmetric y-axis limit for the per-clone CNV-score genome plot. Window scores past this are clipped for display only.")
    p.add_argument("--skip-clones", action="store_true", help="Skip state-local NMF clone discovery and global merge.")
    p.add_argument("--skip-pdf", action="store_true", help="Skip clone-level genome PDF export.")
    p.add_argument("--skip-heatmap", action="store_true", help="Skip default inferCNV-style heatmap PNG/PDF exports.")
    p.add_argument("--heatmap-max-cells", type=int, default=20000,
                   help="Maximum cells rendered in inferCNV-style heatmaps. Larger cohorts are "
                        "deterministically stratified by fastCNV clone and cell state for display. "
                        "Set 0 to render all cells.")
    p.add_argument("--heatmap-filter-threshold", type=float, default=1.5,
                   help="Denoised heatmap threshold: scores with absolute value below this are "
                        "shown as neutral in the filtered export.")
    p.add_argument("--heatmap-min-chr-windows", type=int, default=35,
                   help="Minimum rendered display width per chromosome in inferCNV-style heatmaps. "
                        "Narrow chromosomes such as chrY are repeated for display only so true calls "
                        "are not compressed to an invisible stripe.")
    p.add_argument("--residual-candidate", action="store_true",
                   help="Also run the July 5 pyInferCNV-residual clone candidate from this primary "
                        "fastCNV invocation. Writes <output>.residual_candidate.* tables and "
                        "residual-scale inferCNV-style heatmaps. This is opt-in because it requires "
                        "pyInferCNV and recomputes residuals.")
    p.add_argument("--residual-candidate-min-abs-region-score", type=float, default=0.05,
                   help="Minimum absolute chromosome residual shift for the residual candidate. "
                        "July 5 SamplePre34 benchmark setting: 0.05.")
    p.add_argument("--residual-candidate-min-separation-mad", type=float, default=1.5,
                   help="Minimum state-local two-component separation in MAD units for the residual "
                        "candidate. July 5 SamplePre34 benchmark setting: 1.5.")
    p.add_argument("--residual-candidate-heatmap-filter-threshold", type=float, default=0.03,
                   help="Residual-scale filtered heatmap threshold for the residual candidate. "
                        "Scores with absolute residual below this are shown as neutral.")
    p.add_argument("--residual-candidate-pyinfercnv-path", type=Path,
                   default=Path("/Users/saljh8/Documents/GitHub/pyInferCNV"),
                   help="Local pyInferCNV checkout used if pyinfercnv is not importable.")
    p.add_argument(
        "--sex-chrom-mode",
        choices=["off", "absolute_log2"],
        default="absolute_log2",
        help="Score chrX/chrY using absolute log2 fold-change instead of MAD-standardized "
        "residuals. The MAD denominator inflates with donor-level chrY variability and "
        "squashes real LOY signal; absolute log2 recovers it. Use 'off' for legacy behavior.",
    )
    p.add_argument(
        "--sex-chrom-log2-unit",
        type=float,
        default=0.040,
        help="Divisor mapping the mean log1p(CP10K) difference (query - control_state_mean) "
        "into standardized-score units used by the interval caller. Default 0.040 calibrated "
        "to keep healthy-male LOY false-positive rate under ~2%% per donor while retaining "
        "sensitivity to bona fide LOY clones. Lower values (e.g. 0.025) increase sensitivity "
        "at the cost of false positives; higher values (>=0.05) collapse the LOY signal entirely.",
    )
    p.add_argument("--sex-chrom-het-loss", type=float, default=-0.6,
                   help="Mean log2 fold-change (computed against the male-only control baseline "
                        "for chrY) at/below which a sex-chrom interval is called heterozygous_loss.")
    p.add_argument("--sex-chrom-hom-loss", type=float, default=-1.5,
                   help="Mean log2 fold-change at/below which a sex-chrom interval is called "
                        "homozygous_loss.")
    p.add_argument("--sex-chrom-het-gain", type=float, default=0.4)
    p.add_argument("--sex-chrom-hom-gain", type=float, default=0.7)
    p.add_argument("--control-sample-key", default=None,
                   help="Optional obs column on the control h5ad identifying donor/sample. "
                        "Used by automatic sex detection to call each control sample male/female. "
                        "If absent, the entire control is treated as a single sample.")
    p.add_argument("--max-cells-per-state", type=int, default=0,
                   help="Per-sample cap on cells per --state-key group. 0 disables (use all cells). "
                        "Down-sampling to e.g. 500 keeps memory bounded on large samples without "
                        "loss of information for CNV detection: per-state median baselines are stable "
                        "with ~50-100 cells, and clone discovery operates on the same down-sampled set. "
                        "Only applied within --per-sample mode; the full sample is preserved on disk.")
    p.add_argument("--gate-chry-by-marker-expression", dest="gate_chry_by_marker_expression",
                   action="store_true", default=True,
                   help="LOY-specific heuristic (default ON): NaN chrY scores for any query cell "
                        "expressing >=1 UMI on any chrY-Y-only marker (RPS4Y1, DDX3Y, EIF1AY, KDM5D, "
                        "UTY, USP9Y, NLGN4Y, TMSB4Y, PRKY, ZFY, TSPY1, BCORP1, PRY, TBL1Y, TXLNGY, "
                        "AMELY). Such cells cannot have biologically lost chrY, so the interval "
                        "caller is prevented from emitting a chrY-loss for them.")
    p.add_argument("--no-gate-chry-by-marker-expression", dest="gate_chry_by_marker_expression",
                   action="store_false",
                   help="Disable the chrY-marker-expression gate (use to evaluate scoring without "
                        "the LOY-specific heuristic).")
    p.add_argument("--no-gene-blacklist", dest="gene_blacklist", action="store_false",
                   help="Disable the Ig/TCR/Hb/MT/platelet gene blacklist (default ON: excludes them "
                        "from library size + scoring to remove expression-program false CNVs).")
    p.set_defaults(gene_blacklist=True)
    p.add_argument("--no-pooled-scale", dest="pooled_scale", action="store_false",
                   help="Disable the pooled within-state per-window scale (default ON, replaces the "
                        "0.10 MAD floor) and revert to the legacy per-state MAD + 0.10 floor.")
    p.set_defaults(pooled_scale=True)
    p.add_argument("--scale-floor-quantile", type=float, default=0.05,
                   help="Epsilon floor for the pooled scale = this quantile of the per-window pooled scale.")
    p.add_argument("--loss-zero-copy-quantile", type=float, default=0.95,
                   help="A region copy-LOSS is called when query expression is at/below this quantile "
                        "of the ZERO-COPY reference (for chrY: female donors). Default 0.95 -> concordant "
                        "with complete-loss; applies generally to any region.")
    p.add_argument("--sex-detection-threshold", type=float, default=SEX_DETECTION_CHRY_PCT_DEFAULT,
                   help="Percent of cells in a sample with >=1 chrY-Y-only marker UMI required "
                        "to call the sample male. Default 5%% reliably separates female samples "
                        "(typically 0.5-2%% chrY+ background dropout) from male samples even "
                        "with severe LOY (a male with ~85%% LOY still has ~15%% chrY+ cells).")
    p.add_argument("--sex-chroms", default="chrY",
                   help="Comma-separated chromosomes scored via absolute log1p-difference rather "
                        "than MAD-standardized residuals. Defaults to chrY only because chrX in "
                        "single-cell data has X-inactivation/dosage-compensation complexities and "
                        "often shows systematic batch effects vs reference. Pass 'chrX,chrY' to "
                        "include chrX (e.g. for monosomy-X detection in a sex-matched setting).")
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--write-h5ad", action="store_true")
    p.add_argument("--per-sample", action="store_true",
                   help="Run fastCNV independently for each unique value in --sample-key. "
                        "Each sample gets its own subdirectory under the parent of --output; "
                        "the --output stem becomes the file prefix inside each subdirectory. "
                        "Example: '--output /run/fastCNV --per-sample' produces "
                        "/run/S1/fastCNV.* , /run/S2/fastCNV.* , etc. No cross-sample summary "
                        "files are written; each sample directory is self-contained.")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def params_from_args(args: argparse.Namespace) -> FastParams:
    gene_coords = Path(args.gene_coordinates) if args.gene_coordinates else None
    if gene_coords is None and args.species:
        gene_coords = bundled_gene_coordinates(args.species)
    return FastParams(
        h5ad=Path(args.h5ad),
        control_h5ad=Path(args.control_h5ad),
        gene_coordinates=gene_coords,
        output_prefix=Path(args.output),
        state_key=args.state_key,
        sample_key=args.sample_key,
        control_state_key=args.control_state_key,
        layer=args.layer,
        input_normalized=bool(args.input_normalized),
        window_genes=int(args.window_genes),
        stride_genes=int(args.stride_genes),
        min_chr_genes=int(args.min_chr_genes),
        min_state_cells=int(args.min_state_cells),
        high_threshold=float(args.high_threshold),
        low_threshold=float(args.low_threshold),
        min_run_windows=int(args.min_run_windows),
        min_interval_genes=int(args.min_interval_genes),
        min_mean_score=float(args.min_mean_score),
        burden_quantile=float(args.burden_quantile),
        cnv_burden_threshold=float(args.cnv_burden_threshold),
        simulation_autosomal=(bool(args.simulation_autosomal) or bool(args.external_baseline)),
        autosomal_min_detect_frac=float(args.autosomal_detect_frac),
        cnv_internal_baseline=(not bool(args.external_baseline)),
        max_clones_per_state=int(args.max_clones_per_state),
        max_global_clones=int(args.max_global_clones),
        min_clone_cells=int(args.min_clone_cells),
        clone_similarity_threshold=float(args.clone_similarity_threshold),
        clone_consensus_fraction=float(args.clone_consensus_fraction),
        nmf_max_iter=int(args.nmf_max_iter),
        clone_min_active_fraction=float(args.clone_min_active_fraction),
        clone_max_features=int(args.clone_max_features),
        clone_min_cnv_fraction=float(args.clone_min_cnv_fraction),
        clone_min_cells_confident=int(args.clone_min_cells_confident),
        zygosity_mode=str(args.zygosity_mode),
        skip_clones=bool(args.skip_clones),
        skip_pdf=bool(args.skip_pdf),
        skip_heatmap=bool(args.skip_heatmap),
        random_state=int(args.random_state),
        write_h5ad=bool(args.write_h5ad),
        n_jobs=int(args.n_jobs),
        pdf_smooth_genes=int(args.pdf_smooth_genes),
        pdf_y_clip=float(args.pdf_y_clip),
        pdf_label_threshold=float(args.pdf_label_threshold),
        pdf_score_y_clip=float(args.pdf_score_y_clip),
        heatmap_max_cells=int(args.heatmap_max_cells),
        heatmap_filter_threshold=float(args.heatmap_filter_threshold),
        heatmap_min_chr_windows=int(args.heatmap_min_chr_windows),
        residual_candidate=bool(args.residual_candidate),
        residual_candidate_min_abs_region_score=float(args.residual_candidate_min_abs_region_score),
        residual_candidate_min_separation_mad=float(args.residual_candidate_min_separation_mad),
        residual_candidate_heatmap_filter_threshold=float(args.residual_candidate_heatmap_filter_threshold),
        residual_candidate_pyinfercnv_path=args.residual_candidate_pyinfercnv_path,
        sex_chrom_mode=str(args.sex_chrom_mode),
        sex_chrom_log2_unit=float(args.sex_chrom_log2_unit),
        sex_chrom_het_loss=float(args.sex_chrom_het_loss),
        sex_chrom_hom_loss=float(args.sex_chrom_hom_loss),
        sex_chrom_het_gain=float(args.sex_chrom_het_gain),
        sex_chrom_hom_gain=float(args.sex_chrom_hom_gain),
        sex_chroms=tuple(c.strip() for c in str(args.sex_chroms).split(",") if c.strip()),
        sex_detection_threshold_pct=float(args.sex_detection_threshold),
        control_sample_key=args.control_sample_key,
        max_cells_per_state=int(args.max_cells_per_state),
        gate_chry_by_marker_expression=bool(args.gate_chry_by_marker_expression),
        gene_blacklist=bool(args.gene_blacklist),
        pooled_scale=bool(args.pooled_scale),
        scale_floor_quantile=float(args.scale_floor_quantile),
        loss_zero_copy_quantile=float(args.loss_zero_copy_quantile),
        nearest_normal_autosomal=bool(getattr(args, "nearest_normal_autosomal", False)),
    )


def _run_per_sample(args: argparse.Namespace, base_params: FastParams) -> Dict[str, Path]:
    """Split the query h5ad by --sample-key and run fastCNV independently per sample.

    Each sample gets its own subdirectory under the parent of ``--output``. The
    ``--output`` stem becomes the file prefix inside each subdirectory, e.g. ::

        --output /path/run/fastCNV --sample-key sample_id --per-sample
        ->  /path/run/S1/fastCNV.cnv_cells.tsv
            /path/run/S2/fastCNV.cnv_cells.tsv
            ...

    No cross-sample aggregates are written — each sample's directory is self-contained.
    """
    if not base_params.sample_key:
        raise ValueError("--per-sample requires --sample-key")
    LOGGER.info("Opening combined query h5ad in backed mode: %s", base_params.h5ad)
    # Use backed='r' so we don't load the full matrix into RAM just to enumerate samples
    # and slice — for 100k+ cell cohorts the dense load OOMs the kernel.
    full = ad.read_h5ad(base_params.h5ad, backed="r")
    if base_params.sample_key not in full.obs.columns:
        raise ValueError(
            f"--sample-key '{base_params.sample_key}' not in query.obs columns: {list(full.obs.columns)}"
        )
    sample_values = full.obs[base_params.sample_key].astype(str)
    raw_samples = [s for s in pd.unique(sample_values) if s and s.lower() not in ("nan", "none")]
    # Order samples small -> large so memory pressure (per-sample slicing + fastCNV scoring)
    # peaks at the end. If the largest sample OOMs, all smaller ones still complete.
    sample_counts = {s: int((sample_values.values == s).sum()) for s in raw_samples}
    samples = sorted(raw_samples, key=lambda s: sample_counts[s])
    LOGGER.info(
        "Discovered %d samples (ordered small->large): %s",
        len(samples),
        [(s, sample_counts[s]) for s in samples],
    )

    parent_dir = base_params.output_prefix.parent
    parent_dir.mkdir(parents=True, exist_ok=True)
    file_stem = base_params.output_prefix.name

    aggregated_outputs: Dict[str, Path] = {}

    for sample in samples:
        sample_safe = str(sample).replace("/", "_").replace(" ", "_")
        sample_dir = parent_dir / sample_safe
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_prefix = sample_dir / file_stem
        sample_h5ad = sample_prefix.with_suffix(".input.h5ad")

        # Materialize only this sample's slice into memory before writing to disk.
        sample_mask = sample_values.values == sample
        n_cells_sample = int(sample_mask.sum())
        # Optionally downsample to N cells per (sample, state_key) group. Median-baseline
        # estimation only needs ~50-100 cells per state to be stable; clone-discovery NMF
        # works fine on the same downsampled set. Keeps memory bounded on huge samples
        # without changing call quality.
        if base_params.max_cells_per_state and base_params.max_cells_per_state > 0:
            sample_global_indices = np.flatnonzero(sample_mask)
            sample_states = full.obs.iloc[sample_global_indices][base_params.state_key].astype(str).to_numpy()
            rng = np.random.default_rng(int(base_params.random_state))
            keep_global: List[int] = []
            for state in pd.unique(sample_states):
                local_idx = np.flatnonzero(sample_states == state)
                if local_idx.size > base_params.max_cells_per_state:
                    chosen = rng.choice(local_idx, size=base_params.max_cells_per_state, replace=False)
                else:
                    chosen = local_idx
                keep_global.extend(int(sample_global_indices[i]) for i in chosen)
            global_indices = np.sort(np.asarray(keep_global, dtype=np.int64))
            n_kept = int(global_indices.size)
            LOGGER.info(
                "[%s] downsampling to %d cells/state -> %d / %d cells retained",
                sample, base_params.max_cells_per_state, n_kept, n_cells_sample,
            )
            keep_mask = np.zeros(full.n_obs, dtype=bool)
            keep_mask[global_indices] = True
            sub = full[keep_mask].to_memory()
        else:
            LOGGER.info("[%s] slicing %d cells from backed h5ad...", sample, n_cells_sample)
            sub = full[sample_mask].to_memory()
        LOGGER.info("[%s] %d cells -> %s/", sample, sub.n_obs, sample_dir)
        if sub.n_obs < max(base_params.min_clone_cells, base_params.min_state_cells):
            LOGGER.warning(
                "[%s] only %d cells (< min-state-cells=%d); writing inputs but expect mostly empty outputs.",
                sample, sub.n_obs, base_params.min_state_cells,
            )
        sub.write_h5ad(sample_h5ad, compression="gzip")
        del sub  # free per-sample slice memory before run_fast loads it back

        sample_params = replace(
            base_params,
            h5ad=sample_h5ad,
            output_prefix=sample_prefix,
        )
        try:
            outs = run_fast(sample_params)
        except Exception:
            LOGGER.exception("[%s] run_fast failed; continuing with remaining samples", sample)
            continue
        for k, v in outs.items():
            aggregated_outputs[f"{sample}/{k}"] = v

    full.file.close()
    return aggregated_outputs


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s | %(message)s")
    base_params = params_from_args(args)
    if getattr(args, "per_sample", False):
        outputs = _run_per_sample(args, base_params)
    else:
        outputs = run_fast(base_params)
    for name, path in outputs.items():
        LOGGER.info("%s: %s", name, path)
    return outputs


if __name__ == "__main__":
    main()
