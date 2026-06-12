#!/usr/bin/env python3
"""SATAY-UDON core engine — donor-level covariate enrichment of UDON clusters.

Canonical API (matches the original pyudon `satay_udon.py`, extended):
  make_mean_based_binary(df, cols, direction)             numeric -> binary at column mean (opt-in)
  load_metadata(path, sample_col, mean_binarize, ...)     standard --metadata table -> (binary cov df, donor map)
  fishers_clinical_feats(meta, key, udon_clusters, ...)   Fisher's exact per (cell type x cluster)
  cmh_clinical_feats(meta, key, batch_key, udon_clusters) Cochran-Mantel-Haenszel stratified by batch
  fdr_correction(p_val_matrix, alpha, method)             Benjamini-Hochberg over a matrix
  satay_udon(meta, keys, batch_key, udon_clusters, ...)   run every covariate -> {key: {'p_val','OR','fdr'}}

Difference vs pyudon (per project decision): pseudobulks are collapsed to DISTINCT DONORS before
counting (a donor is feature-positive if ANY of its pseudobulks is positive, and in a cluster if
ANY of its pseudobulks is) and >= MIN_UNITS distinct positive donors are required to test — a
generalization of pyudon's per-pseudobulk `n_samples` floor for datasets with multiple samples per
donor (avoids pseudoreplication). Pass donor_of=None to fall back to per-sample units (== pyudon's
per-pseudobulk counting within a cell type). Pure Python (numpy/scipy/statsmodels). No R.

Metadata format spec: SATAY_METADATA_FORMAT.md.
"""
import os, sys, argparse, numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
try:
    from statsmodels.stats.contingency_tables import StratifiedTable   # CMH
except Exception:
    StratifiedTable = None

MIN_UNITS = 3            # >= this many distinct positive donors (or samples) required to test
MIN_PB = 8               # a cell type needs >= this many pseudobulks (with a defined covariate) to test
DROP_VALUES = {"nan", "Unknown", "None", "n.a.", "NA", "", "na", "NaN", "<NA>"}


# --------------------------------------------------------------------------- #
# binarization + standard metadata loading
# --------------------------------------------------------------------------- #
def make_mean_based_binary(df, cols_to_convert, direction="greater"):
    """Convert numeric columns to binary at the column mean (from pyudon). direction='greater'
    -> 1 if value >= mean; 'lesser' -> 1 if value <= mean (use for vars interesting when low)."""
    numeric = df[cols_to_convert].apply(pd.to_numeric, errors="coerce")
    means = numeric.mean(axis=0, skipna=True)
    out = (numeric <= means) if direction == "lesser" else (numeric >= means)
    return out.astype(int)


def load_metadata(path, sample_col="Sample", mean_binarize=None, direction="greater"):
    """Read + validate the standard SATAY-UDON metadata table (SATAY_METADATA_FORMAT.md).
    Returns (binary covariate DataFrame indexed by Sample, donor_of dict Sample->Donor_ID or None).
    Binary {0,1} columns kept; categorical -> one-hot per level (undefined samples = NaN, not 0);
    numeric columns are an ERROR unless named in `mean_binarize` (then make_mean_based_binary)."""
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str)
    if sample_col not in df.columns:
        raise ValueError(f"[metadata] required column '{sample_col}' not found. Columns: {list(df.columns)}")
    df[sample_col] = df[sample_col].astype(str).str.strip()
    dups = sorted(df[sample_col][df[sample_col].duplicated()].unique())
    if dups:
        raise ValueError(f"[metadata] duplicate {sample_col} values not allowed (one row/sample): {dups[:8]}")
    df = df.set_index(sample_col)
    donor_of = df["Donor_ID"].astype(str).to_dict() if "Donor_ID" in df.columns else None
    mean_binarize = set(mean_binarize or [])
    reserved = {"Donor_ID", "Study", "Dataset"}
    parts, bad = [], []
    for col in [c for c in df.columns if c not in reserved]:
        s = df[col].astype(str).str.strip().where(lambda x: ~x.isin(DROP_VALUES), other=np.nan)
        vals = set(s.dropna().unique())
        if not vals:
            continue
        if vals <= {"0", "1", "0.0", "1.0"}:                              # already binary
            parts.append(s.map(lambda v: np.nan if pd.isna(v) else float(v)).rename(col))
        elif col in mean_binarize:                                        # opt-in numeric -> binary
            parts.append(make_mean_based_binary(pd.DataFrame({col: s}), [col], direction)[col].rename(col))
        elif pd.to_numeric(s.dropna(), errors="coerce").notna().mean() > 0.9 and s.nunique() > 12:
            bad.append(f"{col} (continuous numeric; add to --mean-binarize or pre-bin it)")
        else:                                                             # categorical -> one-hot
            oneh = pd.get_dummies(s).astype(float)
            oneh.columns = [f"{col}={lev}" for lev in oneh.columns]
            oneh[s.isna().values] = np.nan                               # undefined samples -> NaN, not 0
            parts.append(oneh)
    if bad:
        raise ValueError("[metadata] columns neither binary {0,1} nor categorical:\n  " + "\n  ".join(bad))
    if not parts:
        raise ValueError("[metadata] no usable covariate columns after filtering.")
    return pd.concat(parts, axis=1), donor_of


# --------------------------------------------------------------------------- #
# donor-level contingency counting (shared by Fisher + CMH)
# --------------------------------------------------------------------------- #
def _split_ids(udon_clusters):
    uc = udon_clusters.copy()
    uc.columns = ["cluster"] if uc.shape[1] == 1 else uc.columns
    uc["cell_type"] = [str(i).split("__", 1)[0] for i in uc.index]
    uc["sample"] = [str(i).split("__", 1)[1] if "__" in str(i) else str(i) for i in uc.index]
    return uc


def _donor_counts(T, clusters, cl_idx):
    """T: rows for one cell type (and optionally one batch); columns 'cluster','unit','cm' (0/1).
    Returns {cluster: (a,b,c,d)} donor-level: a=donors in-cluster & feature+, b=in-cluster & feature-,
    c=out-cluster & feature+, d=out-cluster & feature-. Donor feature+ if ANY pseudobulk cm==1;
    donor in a cluster if ANY pseudobulk is in it."""
    units = T["unit"].values
    uniq, uinv = np.unique(units, return_inverse=True)
    member = np.zeros((len(uniq), len(clusters)), bool)
    member[uinv, [cl_idx[c] for c in T["cluster"].values]] = True
    donor_f = np.zeros(len(uniq), int)
    np.maximum.at(donor_f, uinv, T["cm"].values.astype(int))
    Fp, ntot, fpos = int(donor_f.sum()), len(uniq), donor_f == 1
    out = {}
    for c in clusters:
        inc = member[:, cl_idx[c]]
        a = int((inc & fpos).sum()); b = int((inc & ~fpos).sum())
        cc = Fp - a; dd = (ntot - (a + b)) - cc
        out[c] = (a, b, cc, dd)
    return out


def _prep(udon_clusters, clinical_metadata_df, clinical_measure_key, donor_of):
    uc = _split_ids(udon_clusters)
    cm = clinical_metadata_df[clinical_measure_key].dropna().astype(float)
    uc = uc[uc["sample"].isin(cm.index)].copy()
    uc["cm"] = uc["sample"].map(cm)
    uc["unit"] = [(donor_of or {}).get(str(s), str(s)) for s in uc["sample"]]
    cell_types = sorted(uc["cell_type"].unique())
    clusters = sorted(uc["cluster"].unique(), key=lambda x: (len(str(x)), str(x)))
    return uc, cell_types, clusters, {c: j for j, c in enumerate(clusters)}


# --------------------------------------------------------------------------- #
# Fisher's exact (donor-level), per (cell type x cluster)
# --------------------------------------------------------------------------- #
def fishers_clinical_feats(clinical_metadata_df, clinical_measure_key, udon_clusters,
                           donor_of=None, n_units=MIN_UNITS, min_pb=MIN_PB, alternative="greater"):
    """Donor-level Fisher's exact for one binary covariate. udon_clusters: index 'celltype__Sample',
    one column = cluster. clinical_metadata_df: index Sample, the covariate column is 0/1 (NaN=undefined).
    Returns {'p_val': df, 'OR': df} (cell types x clusters); cells with < n_units positive donors -> NaN."""
    uc, cell_types, clusters, cl_idx = _prep(udon_clusters, clinical_metadata_df, clinical_measure_key, donor_of)
    p_mat = pd.DataFrame(np.nan, index=cell_types, columns=clusters)
    or_mat = pd.DataFrame(np.nan, index=cell_types, columns=clusters)
    for ct in cell_types:
        T = uc[uc["cell_type"] == ct]
        if len(T) < min_pb:
            continue
        for c, (a, b, cc, dd) in _donor_counts(T, clusters, cl_idx).items():
            try:
                orr, p = fisher_exact([[a, b], [cc, dd]], alternative=alternative)
            except Exception:
                orr, p = np.nan, np.nan
            or_mat.at[ct, c] = orr
            if a >= n_units:                                # >= n_units distinct positive donors
                p_mat.at[ct, c] = p
    return {"p_val": p_mat, "OR": or_mat}


# --------------------------------------------------------------------------- #
# Cochran-Mantel-Haenszel (batch-stratified), per (cell type x cluster) -- ported + FIXED from pyudon
# --------------------------------------------------------------------------- #
def cmh_clinical_feats(clinical_metadata_df, clinical_measure_key, batch_key, udon_clusters,
                       batch_of=None, donor_of=None, n_units=MIN_UNITS, min_pb=MIN_PB):
    """Batch-stratified enrichment: per (cell type x cluster), build a donor-level 2x2 within EACH
    batch and combine with the CMH test (controls for batch). The batch label per Sample comes from
    `batch_of` (a Sample->batch dict) if given, else from the `batch_key` column of the metadata.
    Returns a p-value matrix (cell types x clusters)."""
    if StratifiedTable is None:
        raise ImportError("cmh_clinical_feats requires statsmodels.stats.contingency_tables.StratifiedTable")
    uc, cell_types, clusters, cl_idx = _prep(udon_clusters, clinical_metadata_df, clinical_measure_key, donor_of)
    if batch_of is not None:
        b_of = pd.Series({str(k): str(v) for k, v in dict(batch_of).items()})
    elif batch_key is not None and batch_key in clinical_metadata_df.columns:
        b_of = clinical_metadata_df[batch_key].astype(str)
    else:
        raise ValueError(f"[cmh] need batch_of (Sample->batch) or a batch_key column present in the "
                         f"metadata; got batch_key={batch_key!r}")
    uc["batch"] = uc["sample"].map(b_of)
    p_mat = pd.DataFrame(np.nan, index=cell_types, columns=clusters)
    for ct in cell_types:
        T = uc[uc["cell_type"] == ct]
        if len(T) < min_pb:
            continue
        per_cluster_tables = {c: [] for c in clusters}
        per_cluster_pos = {c: 0 for c in clusters}
        for bv in sorted(x for x in T["batch"].dropna().unique() if str(x) not in DROP_VALUES):
            Tb = T[T["batch"] == bv]
            if Tb["unit"].nunique() < 2:
                continue
            for c, (a, b, cc, dd) in _donor_counts(Tb, clusters, cl_idx).items():
                per_cluster_tables[c].append(np.array([[a, b], [cc, dd]]))
                per_cluster_pos[c] += a
        for c in clusters:
            tabs = [t for t in per_cluster_tables[c] if t.sum() > 0]
            if len(tabs) < 2 or per_cluster_pos[c] < n_units:     # need >=2 informative strata + support
                continue
            try:
                st = StratifiedTable([t.astype(float) for t in tabs])
                p_mat.at[ct, c] = float(st.test_null_odds().pvalue)
            except Exception:
                pass
    return p_mat


# --------------------------------------------------------------------------- #
# FDR + wrapper
# --------------------------------------------------------------------------- #
def fdr_correction(p_val_matrix, alpha=0.05, method="fdr_bh"):
    """Benjamini-Hochberg over all non-NaN entries of a p-value matrix (from pyudon)."""
    flat = p_val_matrix.values.flatten()
    mask = ~np.isnan(flat)
    out = np.full_like(flat, np.nan, dtype=float)
    if mask.sum():
        out[mask] = multipletests(flat[mask], alpha=alpha, method=method)[1]
    return pd.DataFrame(out.reshape(p_val_matrix.shape), index=p_val_matrix.index, columns=p_val_matrix.columns)


def satay_udon(clinical_metadata_df, clinical_measure_keys, udon_clusters, batch_key=None, batch_of=None,
               donor_of=None, n_units=MIN_UNITS, min_pb=MIN_PB, fdr_alpha=0.05):
    """Run SATAY-UDON over every covariate in `clinical_measure_keys`. If a batch is given (batch_key
    column or batch_of map), uses the batch-stratified CMH test; otherwise Fisher's exact.
    Returns {key: {'p_val','OR'(Fisher only),'fdr'}}."""
    use_cmh = batch_key is not None or batch_of is not None
    results = {}
    for key in clinical_measure_keys:
        if use_cmh:
            pv = cmh_clinical_feats(clinical_metadata_df, key, batch_key, udon_clusters,
                                    batch_of=batch_of, donor_of=donor_of, n_units=n_units, min_pb=min_pb)
            res = {"p_val": pv, "OR": None}
        else:
            res = fishers_clinical_feats(clinical_metadata_df, key, udon_clusters,
                                         donor_of=donor_of, n_units=n_units, min_pb=min_pb)
        res["fdr"] = fdr_correction(res["p_val"], alpha=fdr_alpha)
        results[key] = res
    return results


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _results_to_long(results):
    recs = []
    for key, res in results.items():
        pv, fdr, orr = res["p_val"], res["fdr"], res.get("OR")
        for ct in pv.index:
            for cl in pv.columns:
                p = pv.at[ct, cl]
                if pd.isna(p):
                    continue
                recs.append({"covariate": key, "celltype": ct, "cluster": cl, "p": p,
                             "fdr": fdr.at[ct, cl],
                             "odds_ratio": (orr.at[ct, cl] if orr is not None else np.nan)})
    return pd.DataFrame(recs)


def _byfield_mats(long, p_thr=0.1):
    clusters = sorted(long["cluster"].astype(str).unique(), key=lambda x: (len(x), x))
    out = {}
    for cov in sorted(long["covariate"].unique()):
        sub = long[long["covariate"] == cov]
        piv = sub.pivot_table(index="celltype", columns="cluster", values="p", aggfunc="min").reindex(columns=clusters)
        M = piv.fillna(1.0).values
        keep = [i for i in range(len(piv.index)) if (M[i] < p_thr).any()]
        if keep:
            out[cov] = (list(np.array(piv.index)[keep]), M[keep], clusters)
    return out


def render_satay(long, udon_clusters, outdir, donor_of=None, study_of=None, title="SATAY-UDON"):
    """Render the standard SATAY-UDON figures from the enrichment long table:
    satay_byfield (metadata-field-grouped) + donor/cell-type composition heatmaps (Adobe-safe)."""
    os.makedirs(outdir, exist_ok=True)
    try:
        from satay_heatmap_byfield import plot_byfield
        mats = _byfield_mats(long)
        if mats:
            plot_byfield(mats, f"{title} by metadata field", os.path.join(outdir, "satay_byfield"))
        else:
            print("[satay] no field-level enrichments to plot (satay_byfield skipped)")
    except Exception as e:
        print(f"[satay] satay_byfield skipped: {e}")
    try:
        from udon_binary_heatmaps import make_udon_binary_heatmaps
        make_udon_binary_heatmaps(udon_clusters, outdir, donor_of=donor_of, study_of=study_of)
    except Exception as e:
        print(f"[satay] composition heatmaps skipped: {e}")


def main():
    ap = argparse.ArgumentParser(description="SATAY-UDON: donor-level covariate enrichment of UDON clusters")
    ap.add_argument("--metadata", required=True, help="standard metadata table (SATAY_METADATA_FORMAT.md)")
    ap.add_argument("--udon-clusters", required=True, help="udon_clusters.txt (index celltype__Sample, col cluster)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample-col", default="Sample")
    ap.add_argument("--batch-key", default=None, help="metadata column to stratify by (uses CMH instead of Fisher)")
    ap.add_argument("--mean-binarize", default=None, help="comma list of numeric covariate columns to binarize at their mean")
    ap.add_argument("--min-donors", type=int, default=MIN_UNITS)
    ap.add_argument("--min-pb", type=int, default=MIN_PB)
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--heatmaps", action="store_true",
                    help="also render satay_byfield + donor/cell-type composition heatmaps")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    mb = args.mean_binarize.split(",") if args.mean_binarize else None
    binary, donor_of = load_metadata(args.metadata, sample_col=args.sample_col, mean_binarize=mb)
    print(f"[satay] metadata: {binary.shape[1]} binary covariates x {binary.shape[0]} samples; "
          f"donor-level={'yes' if donor_of else 'no (per-sample)'}")
    cl = pd.read_csv(args.udon_clusters, sep="\t", index_col=0)
    cl.columns = ["cluster"]
    print(f"[satay] udon_clusters: {len(cl)} pseudobulks, {cl['cluster'].nunique()} clusters")
    batch_of = None
    if args.batch_key:                                   # batch column lives in the raw metadata (reserved)
        _sep0 = "," if str(args.metadata).lower().endswith(".csv") else "\t"
        _rawb = pd.read_csv(args.metadata, sep=_sep0, dtype=str).drop_duplicates(args.sample_col).set_index(args.sample_col)
        if args.batch_key not in _rawb.columns:
            sys.exit(f"[satay] --batch-key '{args.batch_key}' is not a column in {args.metadata}")
        batch_of = _rawb[args.batch_key].astype(str).to_dict()
        print(f"[satay] batch-stratified (CMH) by '{args.batch_key}': {len(set(batch_of.values()))} batches")
    results = satay_udon(binary, list(binary.columns), cl, batch_key=args.batch_key, batch_of=batch_of,
                         donor_of=donor_of, n_units=args.min_donors, min_pb=args.min_pb, fdr_alpha=args.fdr_alpha)
    long = _results_to_long(results)
    long.to_csv(os.path.join(args.outdir, "satay_cluster_enrichment.tsv"), sep="\t", index=False)
    sig = long[long["fdr"] < args.fdr_alpha] if len(long) else long
    print(f"[satay] {len(long)} testable (celltype x cluster x covariate); {len(sig)} significant (FDR<{args.fdr_alpha})")
    print(f"[satay] wrote {os.path.join(args.outdir, 'satay_cluster_enrichment.tsv')}")
    if args.heatmaps and len(long):
        _sep = "," if str(args.metadata).lower().endswith(".csv") else "\t"
        _raw = pd.read_csv(args.metadata, sep=_sep, dtype=str).set_index(args.sample_col)
        study_of = _raw["Study"].astype(str).to_dict() if "Study" in _raw.columns else None
        render_satay(long, cl, args.outdir, donor_of=donor_of, study_of=study_of,
                     title=f"SATAY-UDON ({'CMH/batch' if args.batch_key else 'Fisher'})")
        print(f"[satay] wrote heatmaps to {args.outdir}")


if __name__ == "__main__":
    main()
