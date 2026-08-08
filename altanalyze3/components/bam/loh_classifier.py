#!/usr/bin/env python3
"""altanalyze3 LOH classifier -- multi-timepoint, driver-anchored single-cell loss-of-heterozygosity.

Faithful generalization of the validated MDS5801 RUNX1/chr21 classifier
(Xuan, 2026-03-05_RUNX1_LOH_classifier/S5_binary_cv.py, AUC-ROC 0.977) to ANY chromosome and ANY
driver. The method needs >=2 timepoints per patient; the driver's per-timepoint bulk VAF trajectory
both selects the informative SNVs and defines the LOH vs HET training labels.

Method (unchanged from the validated pipeline):
  1. Restrict candidate germline het-SNVs to the driver's chromosome.
  2. Depth filter: keep SNVs with >= min_reads_per_group pooled reads in EVERY timepoint group.
  3. Spearman select: keep SNVs whose per-group bulk VAF tracks the driver's per-group VAF trajectory
     (|rho| >= spearman_thresh). +rho = retained allele, -rho = lost allele.
  4. Features: per-cell ALT and REF counts at the selected SNVs (depth-aware, 0 = uncovered).
  5. Labels: cells carrying the anchor mutation (a founding clonal marker), labeled LOH in high-driver-VAF
     groups and HET in het-VAF groups.
  6. Classifier: gradient-boosted trees (sklearn HistGradientBoosting) + L1 logistic baseline, 5-fold CV.
  7. Apply to all cells -> per-cell P(LOH) and the per-group LOH trajectory.

The driver trajectory, LOH groups and HET groups default to data-derived values but can be set
explicitly to reproduce a published run exactly.

Run: python -m altanalyze3.components.bam.loh_classifier --matrix-dir DIR --driver RUNX1p.W279* \
        --chrom chr21 --anchor SRSF2p.P95R --output OUT ...
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, confusion_matrix
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline


# ── data loading ──────────────────────────────────────────────────────────────
def load_matrices(matrix_dir):
    """Read combined_alt.txt, combined_ref.txt (cells x features) and cell_metadata.tsv."""
    alt = pd.read_csv(os.path.join(matrix_dir, "combined_alt.txt"), sep="\t", index_col=0)
    ref = pd.read_csv(os.path.join(matrix_dir, "combined_ref.txt"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(matrix_dir, "cell_metadata.tsv"), sep="\t", index_col=0)
    meta = meta.loc[alt.index]
    return alt, ref, meta


def _bulk_vaf_per_group(alt, ref, feature, groups_of_cell, group_order):
    """Pooled (bulk) VAF and pooled depth of one feature within each group, in group_order."""
    vaf, dep = [], []
    for g in group_order:
        cells = groups_of_cell.index[groups_of_cell == g]
        a = float(alt.loc[cells, feature].sum())
        r = float(ref.loc[cells, feature].sum())
        d = a + r
        vaf.append(a / d if d > 0 else np.nan)
        dep.append(d)
    return np.array(vaf), np.array(dep)


# ── SNV selection (depth filter -> Spearman vs driver trajectory) ──────────────
def select_informative_snvs(alt, ref, groups_of_cell, group_order, chrom, driver_traj,
                            min_reads_per_group=100, spearman_thresh=0.70):
    """Return (selected_snvs, snv_table). SNVs restricted to `chrom`; kept if pooled depth
    >= min_reads_per_group in EVERY group AND |Spearman rho(per-group VAF, driver_traj)| >= thresh.
    Vectorized: one groupby-sum builds the group x SNV pooled ALT/REF tables."""
    snv_cols = [c for c in alt.columns if c.startswith(chrom + ":")]
    A = alt[snv_cols].groupby(groups_of_cell).sum().reindex(group_order)      # group x SNV pooled ALT
    R = ref[snv_cols].groupby(groups_of_cell).sum().reindex(group_order)      # group x SNV pooled REF
    D = A + R
    V = A / D.replace(0, np.nan)                                              # group x SNV bulk VAF
    depth_pass = (D >= min_reads_per_group).all(axis=0)                       # per SNV
    traj = np.asarray(driver_traj, dtype=float)
    rho = pd.Series(np.nan, index=snv_cols)
    for snv in snv_cols:
        v = V[snv].values
        if depth_pass[snv] and not np.any(np.isnan(v)):
            rho[snv] = spearmanr(v, traj)[0]
    selected = [s for s in snv_cols if depth_pass[s] and not np.isnan(rho[s]) and abs(rho[s]) >= spearman_thresh]
    sel_set = set(selected)
    tbl = {"snv_id": snv_cols, "spearman_rho": rho.values,
           "depth_pass": depth_pass.values, "selected": [s in sel_set for s in snv_cols]}
    for g in group_order:
        tbl[f"vaf_{g}"] = V.loc[g].values
        tbl[f"dep_{g}"] = D.loc[g].values
    return selected, pd.DataFrame(tbl)


# ── labels (anchor clone x group) ──────────────────────────────────────────────
def derive_loh_het_groups(driver_traj, group_order, loh_vaf_min=0.75, het_vaf=(0.35, 0.65)):
    """Data-derived default: a group is LOH if driver VAF >= loh_vaf_min, HET if within het_vaf band."""
    loh = [g for g, v in zip(group_order, driver_traj) if not np.isnan(v) and v >= loh_vaf_min]
    het = [g for g, v in zip(group_order, driver_traj) if not np.isnan(v) and het_vaf[0] <= v <= het_vaf[1]]
    return loh, het


def build_labels(alt, groups_of_cell, anchor, loh_groups, het_groups, anchor_min_alt=1):
    """LOH=anchor-mutant cell in a loh_group; HET=anchor-mutant cell in a het_group; else unlabeled."""
    anchor_mut = alt[anchor] >= anchor_min_alt
    labels = pd.Series("unlabeled", index=alt.index)
    labels[anchor_mut & groups_of_cell.isin(loh_groups)] = "LOH"
    labels[anchor_mut & groups_of_cell.isin(het_groups)] = "HET"
    return labels


# ── features + classifier ──────────────────────────────────────────────────────
def build_features(alt, ref, selected_snvs):
    data = {}
    for snv in selected_snvs:
        data[snv + "__ALT"] = alt[snv].values.astype(float)
        data[snv + "__REF"] = ref[snv].values.astype(float)
    return pd.DataFrame(data, index=alt.index)


def _model(kind):
    if kind == "lr":
        return Pipeline([("scale", StandardScaler()),
                         ("lr", LogisticRegression(penalty="l1", solver="liblinear", C=1.0,
                                                   random_state=42, max_iter=1000))])
    # gradient-boosted trees: sklearn analogue of the validated XGBoost model
    return HistGradientBoostingClassifier(max_depth=5, learning_rate=0.05, max_iter=300,
                                          l2_regularization=1.0, random_state=42)


def crossval(X, y, kind="gbm", n_splits=5):
    """Stratified n-fold out-of-fold probabilities + AUC-ROC / AUC-PR.

    Defensive against degenerate label vectors: returns (oof, nan, nan) when y is empty or
    has a single class (AUC/AP are undefined there), and clamps the fold count to the smallest
    class so StratifiedKFold never raises. Callers should guard degenerate y upstream; this keeps
    the function from crashing an otherwise-completing run."""
    y = np.asarray(y)
    oof = np.full(len(y), np.nan, dtype=float)
    classes, counts = np.unique(y, return_counts=True)
    if len(y) == 0 or len(classes) < 2:
        return oof, float("nan"), float("nan")
    k = int(min(n_splits, counts.min()))
    if k < 2:
        return oof, float("nan"), float("nan")
    cv = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
    oof = np.zeros(len(y), dtype=float)
    for tr, va in cv.split(X, y):
        m = _model(kind)
        m.fit(X[tr], y[tr])
        oof[va] = m.predict_proba(X[va])[:, 1]
    return oof, roc_auc_score(y, oof), average_precision_score(y, oof)


def _write_status(out_dir, status, reason, driver, chrom, anchor, group_order, driver_traj,
                  loh_groups, het_groups, n_loh, n_het, n_selected):
    """Write loh_status.txt so a patient with no trainable LOH signal COMPLETES with a clear,
    machine-readable reason instead of crashing the LSF task."""
    os.makedirs(out_dir, exist_ok=True)
    traj = {g: (None if np.isnan(v) else round(float(v), 4)) for g, v in zip(group_order, driver_traj)}
    with open(os.path.join(out_dir, "loh_status.txt"), "w") as fh:
        fh.write("status\t%s\n" % status)
        fh.write("reason\t%s\n" % reason)
        fh.write("driver\t%s\nchrom\t%s\nanchor\t%s\n" % (driver, chrom, anchor))
        fh.write("driver_trajectory\t%s\n" % traj)
        fh.write("loh_groups\t%s\nhet_groups\t%s\n" % (loh_groups, het_groups))
        fh.write("n_selected_snvs\t%s\n" % ("" if n_selected is None else n_selected))
        fh.write("n_labeled_LOH\t%s\nn_labeled_HET\t%s\n" % (
            "" if n_loh is None else n_loh, "" if n_het is None else n_het))
    print("[loh_status] %s: %s" % (status, reason))


# ── orchestration ──────────────────────────────────────────────────────────────
def run_loh_classifier(matrix_dir, driver, chrom, anchor, out_dir, group_col="stage",
                       group_order=None, loh_groups=None, het_groups=None,
                       min_reads_per_group=100, spearman_thresh=0.70, model="gbm",
                       driver_traj_override=None):
    os.makedirs(out_dir, exist_ok=True)
    alt, ref, meta = load_matrices(matrix_dir)
    goc = meta[group_col]
    if group_order is None:
        group_order = list(pd.unique(goc))
    for col in (driver, anchor):
        if col not in alt.columns:
            raise SystemExit(f"column not in matrices: {col}")

    # 1. driver trajectory (per-group bulk VAF) — the selection + label anchor
    if driver_traj_override is not None:
        driver_traj = np.array([driver_traj_override[g] for g in group_order], dtype=float)
    else:
        driver_traj, _ = _bulk_vaf_per_group(alt, ref, driver, goc, group_order)
    print("driver %s trajectory: %s" % (driver, dict(zip(group_order, np.round(driver_traj, 3)))))

    # 2. LOH/HET groups
    if loh_groups is None or het_groups is None:
        dl, dh = derive_loh_het_groups(driver_traj, group_order)
        loh_groups = loh_groups or dl
        het_groups = het_groups or dh
    print("LOH groups: %s | HET groups: %s" % (loh_groups, het_groups))

    # 3. SNV selection — restricted to timepoints where the driver is COVERED. A timepoint with no
    # driver coverage (VAF = nan) carries no trajectory signal, yet the "depth in EVERY group" filter
    # would let it veto every SNV (this is why 5801's cohort run selected 0 SNVs while its bespoke run
    # gave AUC 0.977 — the sAML timepoint had no coverage). A trajectory needs >= 2 covered timepoints.
    finite = ~np.isnan(driver_traj)
    sel_groups = [g for g, f in zip(group_order, finite) if f]
    sel_traj = driver_traj[finite]
    if len(sel_groups) < 2:
        _write_status(out_dir, "insufficient_timepoints",
                      "only %d timepoint(s) have driver coverage; a trajectory needs >= 2" % len(sel_groups),
                      driver, chrom, anchor, group_order, driver_traj, loh_groups, het_groups, None, None, 0)
        print("INSUFFICIENT TIMEPOINTS (%d covered) -- wrote loh_status.txt, skipping classifier" % len(sel_groups))
        return None, None
    selected, snv_tbl = select_informative_snvs(alt, ref, goc, sel_groups, chrom, sel_traj,
                                                min_reads_per_group, spearman_thresh)
    snv_tbl.sort_values("spearman_rho", ascending=False).to_csv(
        os.path.join(out_dir, "snv_selection.csv"), index=False)
    npos = int((snv_tbl.loc[snv_tbl.selected, "spearman_rho"] > 0).sum())
    nneg = int((snv_tbl.loc[snv_tbl.selected, "spearman_rho"] < 0).sum())
    print("selected %d %s SNVs (retained +rho=%d, lost -rho=%d) at dp>=%d, |rho|>=%.2f"
          % (len(selected), chrom, npos, nneg, min_reads_per_group, spearman_thresh))
    if not selected:
        _write_status(out_dir, "no_informative_snvs",
                      "no %s SNV passed depth >= %d in every covered timepoint with |rho| >= %.2f "
                      "vs the driver trajectory" % (chrom, min_reads_per_group, spearman_thresh),
                      driver, chrom, anchor, group_order, driver_traj, loh_groups, het_groups, 0, 0, 0)
        print("NO INFORMATIVE SNVs -- wrote loh_status.txt, skipping classifier")
        return None, None

    # 4. features + labels
    X_all = build_features(alt, ref, selected)
    labels = build_labels(alt, goc, anchor, loh_groups, het_groups)
    lab = labels.isin(["LOH", "HET"])
    Xl = X_all[lab].values
    yl = (labels[lab] == "LOH").astype(int).values
    print("labeled cells: LOH=%d HET=%d (anchor=%s)" % (int(yl.sum()), int((yl == 0).sum()), anchor))
    n_loh, n_het = int(yl.sum()), int((yl == 0).sum())
    if n_loh < 1 or n_het < 1:
        reason = ("no anchor-mutant cell falls in both an LOH and a HET timepoint" if len(yl) == 0
                  else ("only one training class present (LOH=%d, HET=%d): the driver does not "
                        "contrast LOH vs HET across timepoints" % (n_loh, n_het)))
        _write_status(out_dir, "not_trainable", reason, driver, chrom, anchor, group_order,
                      driver_traj, loh_groups, het_groups, n_loh, n_het, len(selected))
        print("NOT TRAINABLE: %s -- wrote loh_status.txt, skipping classifier" % reason)
        return None, None

    # 5. cross-validated metrics (both models)
    metrics = []
    oof_main = None
    for kind in ("gbm", "lr"):
        oof, auc, ap = crossval(Xl, yl, kind=kind)
        metrics.append({"model": kind, "AUC_ROC": auc, "AUC_PR": ap})
        print("  %-3s  AUC-ROC=%.4f  AUC-PR=%.4f" % (kind, auc, ap))
        if kind == model:
            oof_main = oof
    pd.DataFrame(metrics).to_csv(os.path.join(out_dir, "model_metrics.csv"), index=False)

    # 6. fit final model on all labeled, predict every cell
    final = _model(model)
    final.fit(Xl, yl)
    prob_all = final.predict_proba(X_all.values)[:, 1]
    pred = pd.DataFrame({"cell": X_all.index, group_col: goc.values,
                         "label": labels.values, "loh_prob": prob_all}).set_index("cell")
    pred.to_csv(os.path.join(out_dir, "cell_loh_probabilities.csv"))

    # mutant-retained LOH: in predicted-LOH cells, is the driver's mutant allele the ONLY one expressed
    # (VAF -> 1.0, WT allele lost)? This is the clinically relevant readout (e.g. RUNX1 in 5801).
    _driver_retention_report(alt, ref, pred, [d for d in (driver, anchor) if d in alt.columns], out_dir)

    # 7. per-group trajectory + confusion (out-of-fold, threshold 0.5)
    traj = []
    for g in group_order:
        v = pred.loc[pred[group_col] == g, "loh_prob"].values
        traj.append({"group": g, "driver_vaf": float(driver_traj[group_order.index(g)]),
                     "mean_P_LOH": float(np.mean(v)), "median_P_LOH": float(np.median(v)),
                     "n_cells": int(len(v))})
    traj_df = pd.DataFrame(traj)
    traj_df.to_csv(os.path.join(out_dir, "per_group_trajectory.csv"), index=False)
    print("\nper-group P(LOH) trajectory:")
    print(traj_df.to_string(index=False))
    cm = confusion_matrix(yl, (oof_main >= 0.5).astype(int))
    print("\nconfusion (rows true HET/LOH, cols pred):\n", cm)
    _plot(out_dir, driver, chrom, group_order, traj_df, yl, oof_main, metrics, loh_groups)
    return traj_df, metrics


def _driver_retention_report(alt, ref, pred, drivers, out_dir):
    """For each driver: pooled mutant-allele VAF in predicted-LOH cells vs predicted-HET cells.
    LOH-cell VAF -> 1.0 means the WT allele was lost and only the mutant driver is expressed."""
    loh = pred.index[pred["loh_prob"] >= 0.5]
    het = pred.index[pred["loh_prob"] < 0.5]
    rows = []
    for d in drivers:
        for cls, cells in (("LOH", loh), ("HET", het)):
            a = float(alt.loc[cells, d].sum()); r = float(ref.loc[cells, d].sum()); tot = a + r
            rows.append({"driver": d, "cell_class": cls, "n_cells": len(cells),
                         "mut_umis": int(a), "wt_umis": int(r),
                         "mutant_vaf": (a / tot) if tot else np.nan})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "driver_loh_retention.csv"), index=False)
    print("\nmutant-retained LOH report (driver mutant-allele VAF, predicted-LOH vs predicted-HET cells):")
    for d in drivers:
        vl = next((x["mutant_vaf"] for x in rows if x["driver"] == d and x["cell_class"] == "LOH"), np.nan)
        vh = next((x["mutant_vaf"] for x in rows if x["driver"] == d and x["cell_class"] == "HET"), np.nan)
        tag = "  <- MUTANT-RETAINED LOH (WT allele lost, mutant-only expression)" if (vl == vl and vl > 0.9) else ""
        print(f"  {d:22s} LOH-cell VAF={vl:.2f}  HET-cell VAF={vh:.2f}{tag}")


# ── figure (vector, Arial, hex colors, 2-point lines) ──────────────────────────
def _plot(out_dir, driver, chrom, group_order, traj_df, yl, oof, metrics, loh_groups):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    auc = next(m["AUC_ROC"] for m in metrics if m["model"] == "gbm")
    fpr, tpr, _ = roc_curve(yl, oof)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.2))
    a0.plot([0, 1], [0, 1], ls="--", lw=0.8, color="#999999")
    a0.plot(fpr, tpr, lw=2, color="#d62728", label="GBM  AUC=%.3f" % auc)
    a0.set_xlabel("False positive rate"); a0.set_ylabel("True positive rate")
    a0.set_title("%s LOH vs HET (5-fold CV)" % driver, fontsize=10, fontweight="bold")
    a0.legend(fontsize=8); a0.set_xlim(0, 1); a0.set_ylim(0, 1.02)
    x = np.arange(len(group_order))
    for xi, g in enumerate(group_order):
        if g in loh_groups:
            a1.axvspan(xi - 0.45, xi + 0.45, color="#ffe0e0", alpha=0.5)
    a1.plot([-0.5, len(group_order) - 0.5], [0.5, 0.5], ls="--", lw=0.8, color="#999999")
    a1.plot(x, traj_df["mean_P_LOH"].values, lw=2, color="#d62728", marker="o", ms=6, label="mean P(LOH)")
    a1.plot(x, traj_df["driver_vaf"].values, lw=2, color="#1f77b4", marker="s", ms=5, label="driver bulk VAF")
    a1.set_xticks(x); a1.set_xticklabels(group_order, rotation=30, ha="right", fontsize=8)
    a1.set_ylabel("P(LOH) / VAF"); a1.set_ylim(-0.05, 1.1)
    a1.set_title("%s per-timepoint LOH trajectory" % chrom, fontsize=10, fontweight="bold")
    a1.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "loh_classifier.pdf"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-dir", required=True, help="dir with combined_alt.txt/combined_ref.txt/cell_metadata.tsv")
    ap.add_argument("--driver", required=True, help="driver mutation column (defines the VAF trajectory)")
    ap.add_argument("--chrom", required=True, help="chromosome of the driver (SNV features restricted to it)")
    ap.add_argument("--anchor", required=True, help="clonal-marker mutation column for labeling (e.g. the founding driver)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--group-col", default="stage")
    ap.add_argument("--group-order", nargs="*", default=None, help="timepoint order (default: as encountered)")
    ap.add_argument("--loh-groups", nargs="*", default=None, help="override LOH timepoints (default: driver VAF>=0.75)")
    ap.add_argument("--het-groups", nargs="*", default=None, help="override HET timepoints (default: 0.35<=driver VAF<=0.65)")
    ap.add_argument("--min-reads", dest="min_reads", type=int, default=100, help="min pooled reads/group per SNV")
    ap.add_argument("--spearman", type=float, default=0.70, help="min |Spearman rho| vs driver trajectory")
    ap.add_argument("--model", default="gbm", choices=["gbm", "lr"])
    ap.add_argument("--driver-traj", dest="driver_traj", default=None,
                    help="comma-separated per-group driver VAF (in --group-order) to override the "
                         "data-derived trajectory (reproduce a curated published run exactly)")
    args = ap.parse_args()
    traj_override = None
    if args.driver_traj:
        vals = [float(x) for x in args.driver_traj.split(",")]
        if args.group_order is None or len(vals) != len(args.group_order):
            ap.error("--driver-traj needs one value per --group-order entry")
        traj_override = dict(zip(args.group_order, vals))
    run_loh_classifier(args.matrix_dir, args.driver, args.chrom, args.anchor, args.output,
                       group_col=args.group_col, group_order=args.group_order,
                       loh_groups=args.loh_groups, het_groups=args.het_groups,
                       min_reads_per_group=args.min_reads, spearman_thresh=args.spearman, model=args.model,
                       driver_traj_override=traj_override)


if __name__ == "__main__":
    main()
