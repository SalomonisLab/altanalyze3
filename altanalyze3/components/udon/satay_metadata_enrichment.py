#!/usr/bin/env python3
"""
[LEGACY] SATAY-UDON metadata enrichment. The CANONICAL, pyudon-API-compatible engine is now
`satay_udon_core.py` (make_mean_based_binary / fishers_clinical_feats / cmh_clinical_feats /
fdr_correction / satay_udon + the standard --metadata loader), driven by `run_workflow.py`. This
module is retained for the multi-run cross-comparison orchestration; new work should use the core.

SATAY-UDON metadata enrichment: test each UDON cluster for enrichment of harmonized
AML metadata, run identically across three UDON analyses to contrast sensitivity and
specificity, and to compare study-level vs batch-corrected merges.

Feature categories (considered SEPARATELY):
  mutation : 53 genotype/fusion columns (Mutation_Matrix)
  clinical : disease_category, timepoint, is_pediatric, sex, FAB, ELN_risk, WHO_classification
  celltype : the cell population of each pseudobulk (Hs-BM-titrated-reference-centroid)

Test: per (cluster, feature), pseudobulk-level Fisher's exact (a cluster's pseudobulks
that are feature+ vs the rest), BH-FDR within each run x category. Background for a
mutation/clinical feature = pseudobulks whose Sample has that field defined.

Runs:
  1. per_study  (each of 8 studies analyzed within itself -> study-level, batch-free)
  2. matched_sex_platform  (merged, sex+platform matched -> batch-corrected merge)
  3. udon_core  (original merged run -> uncorrected merge)

Pure Python (pandas/scipy). No R.
"""
import os, glob, numpy as np, pandas as pd
from scipy.stats import fisher_exact

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"
XLSX = os.path.join(PB, "UDON", "AML_harmonized_metadata.xlsx")   # updated/re-curated mutations
OUT = os.path.join(PB, "UDON", "satay_metadata"); os.makedirs(OUT, exist_ok=True)
FDR_THR = 0.05
CLIN_CATS = ["disease_category", "timepoint", "is_pediatric", "sex", "FAB", "ELN_risk", "WHO_classification"]

# harmonize redundant/synonymous metadata values (data cleaning, applied in build_features)
NORMALIZE = {
    "ELN_risk": {"FAV": "Favorable", "LR": "Favorable", "Favorable": "Favorable",
                 "INT": "Intermediate", "IR": "Intermediate", "Intermediate": "Intermediate",
                 "ADV": "Adverse", "HR": "Adverse", "Adverse": "Adverse", "n.a.": None},
    "FAB": {"M5A": "M5a", "M5B": "M5b"},
}
DROP_VALUES = {"nan", "Unknown", "None", "n.a.", "NA", "", "na"}

# --- distinct-unit (donor/sample) support: an enrichment must be backed by >= MIN_UNITS
# distinct DONORS (or samples), never by many pseudobulks of one donor (pseudoreplication) ---
UNIT_MODE = "donor"        # "donor" (default, ideal) | "sample"
MIN_UNITS = 3              # >= this many distinct donors/samples required for an enrichment
DONOR_OF = {}              # Sample -> Donor_ID (populated by build_features)


def _units(samples):
    if UNIT_MODE == "donor":
        return np.array([DONOR_OF.get(str(s), str(s)) for s in samples])
    return np.asarray([str(s) for s in samples])


def unit_fisher(samples, in_cluster, fpos):
    """Collapse pseudobulks to distinct units (donor/sample) -> one vote per donor, then
    one-sided Fisher. A unit is in-cluster / feature+ if ANY of its pseudobulks is.
    Returns (p, odds, a, cluster_units, feature_units); p=1 if a < MIN_UNITS."""
    g = (pd.DataFrame({"u": _units(samples), "inc": np.asarray(in_cluster, int),
                       "f": np.asarray(fpos, int)})
         .groupby("u").agg(inc=("inc", "max"), f=("f", "max")))
    a = int(((g["inc"] == 1) & (g["f"] == 1)).sum()); b = int(((g["inc"] == 1) & (g["f"] == 0)).sum())
    c = int(((g["inc"] == 0) & (g["f"] == 1)).sum()); d = int(((g["inc"] == 0) & (g["f"] == 0)).sum())
    if a < MIN_UNITS or (a + c) < MIN_UNITS:
        return 1.0, 0.0, a, a + b, a + c
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return p, orr, a, a + b, a + c


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); r = np.empty(n); r[o] = np.arange(1, n + 1)
    q = p * n / r
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


MIN_LEVEL = 5          # a categorical level must have >= this many samples to be testable


def _features_from_table(df, sample_col="Sample"):
    """Parse + STRICTLY validate the standard SATAY-UDON metadata table (see
    SATAY_METADATA_FORMAT.md). Raises ValueError on any spec violation rather than coercing."""
    global DONOR_OF, UNIT_MODE
    if sample_col not in df.columns:
        raise ValueError(f"[metadata] required column '{sample_col}' not found. Columns: {list(df.columns)}")
    df = df.copy(); df[sample_col] = df[sample_col].astype(str).str.strip()
    dup = sorted(df[sample_col][df[sample_col].duplicated()].unique())
    if dup:
        raise ValueError(f"[metadata] duplicate {sample_col} values not allowed (one row per sample): {dup[:8]}")
    df = df.set_index(sample_col)
    if "Donor_ID" in df.columns:
        DONOR_OF = dict(zip(df.index, df["Donor_ID"].astype(str)))
    if not DONOR_OF:
        UNIT_MODE = "sample"          # no Donor_ID -> distinct samples
    cov_cols = [c for c in df.columns if c not in ("Donor_ID", "Study", "Dataset")]
    if not cov_cols:
        raise ValueError("[metadata] no covariate columns (only Sample/Donor_ID/Study present).")
    feat, defined, bad = {}, {}, []
    for col in cov_cols:
        present = df[col].astype(str).str.strip()
        present = present[~present.isin(DROP_VALUES)]
        vals = set(present.values)
        if not vals:
            continue
        if vals <= {"0", "1", "0.0", "1.0"}:                       # binary covariate -> presence test
            dfn = set(present.index)
            feat[("binary", col)] = {s: int(float(present[s]) == 1) for s in dfn}
            defined[("binary", col)] = dfn
            continue
        numeric_frac = pd.to_numeric(present, errors="coerce").notna().mean()
        if numeric_frac > 0.9 and present.nunique() > 12:          # continuous numeric -> must pre-bin
            bad.append(f"{col} (continuous: {present.nunique()} numeric values -> pre-bin into binary/categorical)")
            continue
        dfn = set(present.index)                                    # categorical -> per-level tests
        for lev, n in present.value_counts().items():
            if n < MIN_LEVEL:
                continue
            feat[("categorical", f"{col}={lev}")] = {s: int(present.get(s) == lev) for s in dfn}
            defined[("categorical", f"{col}={lev}")] = dfn
    if bad:
        raise ValueError("[metadata] columns that are neither binary {0,1} nor categorical:\n  " +
                         "\n  ".join(bad))
    if not feat:
        raise ValueError("[metadata] no testable covariates (all columns empty or below MIN_LEVEL).")
    return feat, defined


def build_features(metadata_path=None):
    """SATAY-UDON covariate features. PRIMARY input = a standard tab-delimited metadata table
    (--metadata / UDON_METADATA env; see _features_from_table for the spec). The legacy 2-tab AML
    xlsx is only a deprecated fallback when no standard file is provided."""
    global DONOR_OF, UNIT_MODE
    metadata_path = metadata_path or os.environ.get("UDON_METADATA")
    if metadata_path:
        sep = "," if str(metadata_path).lower().endswith(".csv") else "\t"
        return _features_from_table(pd.read_csv(metadata_path, sep=sep, dtype=str))
    import sys as _sys
    print("[satay] WARNING: no --metadata / UDON_METADATA given -- falling back to the LEGACY AML "
          "xlsx. Provide a standard metadata file for reproducible runs.", file=_sys.stderr)
    mm = pd.read_excel(XLSX, sheet_name="Mutation_Matrix")
    clin = pd.read_excel(XLSX, sheet_name="Clinical_Metadata")
    if "Donor_ID" in clin.columns:
        DONOR_OF = dict(zip(clin["Sample"].astype(str), clin["Donor_ID"].astype(str)))
    if not DONOR_OF:
        UNIT_MODE = "sample"          # donor IDs not supplied -> fall back to distinct samples
    mut_cols = [c for c in mm.columns if c not in ("Dataset", "Sample")]
    feat = {}           # (category, feature) -> {sample: 0/1}
    defined = {}        # (category, feature) -> set(samples with the field defined)
    # mutations: NaN/0 -> not mutated; sample defined if it's in the mutation matrix
    mm_samples = set(mm["Sample"].astype(str))
    for g in mut_cols:
        col = mm.set_index(mm["Sample"].astype(str))[g]
        feat[("mutation", g)] = {s: int(pd.notna(v) and float(v) == 1) for s, v in col.items()}
        defined[("mutation", g)] = mm_samples
    # clinical categoricals -> binary indicators per non-trivial level (harmonized)
    clin = clin.set_index(clin["Sample"].astype(str))
    clin = clin[~clin.index.duplicated(keep="first")]                  # dedupe duplicate samples
    for cat in CLIN_CATS:
        if cat not in clin.columns: continue
        nmap = NORMALIZE.get(cat, {})
        s = clin[cat].astype(str).map(lambda v: nmap.get(v, v))        # harmonize synonyms
        dfn = set(s[~s.isin(DROP_VALUES) & s.notna()].index)
        s = s[s.index.isin(dfn)]
        for lev, n in s.value_counts().items():
            if lev in DROP_VALUES or n < 5: continue                   # skip empties / tiny levels
            feat[("clinical", f"{cat}={lev}")] = {smp: int(s.get(smp) == lev) for smp in dfn}
            defined[("clinical", f"{cat}={lev}")] = dfn
    return feat, defined


def load_run(path):
    cl = pd.read_csv(path, sep="\t", index_col=0); cl.columns = ["cluster"]
    # ID = celltype__Sample. The Sample can itself contain "__" (WashU/EBI/Milan),
    # but a cell type never does -> split on the FIRST "__".
    cl["celltype"] = [i.split("__", 1)[0] for i in cl.index]
    cl["Sample"] = [i.split("__", 1)[1] if "__" in i else i for i in cl.index]
    # normalize raw UDON cluster ids to U# (per-study clusters are bare integers)
    cl["cluster"] = cl["cluster"].astype(str).apply(lambda x: x if x[:1].isalpha() else f"U{x}")
    return cl


def enrich(cl, feat, defined, run_label):
    rows = []
    clusters = sorted(cl["cluster"].unique())
    # mutation + clinical (+ demographic): collapse to distinct donors/samples (>= MIN_UNITS required)
    for k, smap in feat.items():
        cat, fname = (k if len(k) == 2 else ("demographic", k[0]))   # age key is a 1-tuple
        dfn = defined[k]
        sub = cl[cl["Sample"].isin(dfn)]
        if len(sub) < 10: continue
        fpos = sub["Sample"].map(smap).fillna(0).astype(int).values
        samples = sub["Sample"].values; clusvals = sub["cluster"].values
        if len(set(_units(samples[fpos == 1]))) < MIN_UNITS: continue   # < MIN_UNITS positive donors total
        for c in clusters:
            p, orr, a, cn, fn = unit_fisher(samples, (clusvals == c), fpos)
            if a < MIN_UNITS: continue
            rows.append([run_label, cat, fname, c, a, cn, fn, len(set(_units(samples))), orr, p, UNIT_MODE])
    # celltype: also donor-collapsed
    cts = cl["celltype"].value_counts(); cts = cts[cts >= 10].index
    samples_all = cl["Sample"].values; clus_all = cl["cluster"].values
    nunit_all = len(set(_units(samples_all)))
    for ct in cts:
        ctpos = (cl["celltype"] == ct).astype(int).values
        for c in clusters:
            p, orr, a, cn, fn = unit_fisher(samples_all, (clus_all == c), ctpos)
            if a < MIN_UNITS: continue
            rows.append([run_label, "celltype", ct, c, a, cn, fn, nunit_all, orr, p, UNIT_MODE])
    df = pd.DataFrame(rows, columns=["run", "category", "feature", "cluster", "overlap",
                                     "cluster_n", "feature_n", "bg_n", "odds_ratio", "p", "unit"])
    # NOTE: 'overlap'=distinct positive units (donors) in the cluster (>= MIN_UNITS required);
    # 'feature_n'/'bg_n' are distinct-unit counts; 'unit' column records donor vs sample.
    # FDR within run x category
    df["fdr"] = np.nan
    for cat, idx in df.groupby("category").groups.items():
        df.loc[idx, "fdr"] = bh_fdr(df.loc[idx, "p"].values)
    return df


def main():
    feat, defined = build_features()
    print(f"features: {sum(1 for k in feat if k[0]=='mutation')} mutation, "
          f"{sum(1 for k in feat if k[0]=='clinical')} clinical levels")

    runs = {}
    # 1) per-study (each study separately)
    perstudy = []
    for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt"))):
        st = os.path.basename(os.path.dirname(d))
        cl = load_run(d)
        e = enrich(cl, feat, defined, f"per_study:{st}")
        perstudy.append(e)
    runs["per_study"] = pd.concat(perstudy, ignore_index=True)
    # 2) matched_sex_platform (batch-corrected merge)
    runs["matched_sex_platform"] = enrich(load_run(os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt")),
                                          feat, defined, "matched_sex_platform")
    # 3) udon_core (uncorrected merge)
    runs["udon_core"] = enrich(load_run(os.path.join(PB, "UDON/udon_core/udon_clusters.txt")),
                               feat, defined, "udon_core")
    # 4) study_aware_final (final batch-corrected programs P#)
    fa = os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv")
    if os.path.exists(fa):
        f = pd.read_csv(fa, sep="\t")
        clf = pd.DataFrame(index=f["pseudobulk"])
        clf["cluster"] = f["final_program"].astype(str).values
        clf["Sample"] = f["Sample"].astype(str).values
        clf["celltype"] = [str(p).split("__", 1)[0] for p in f["pseudobulk"]]
        runs["study_aware_final"] = enrich(clf, feat, defined, "study_aware_final")

    alldf = pd.concat(runs.values(), ignore_index=True)
    alldf.to_csv(os.path.join(OUT, "all_enrichments.tsv"), sep="\t", index=False)
    sig = alldf[alldf["fdr"] < FDR_THR].copy()
    sig.sort_values(["run", "category", "fdr"]).to_csv(os.path.join(OUT, "significant_enrichments.tsv"), sep="\t", index=False)

    # ---- contrast: sensitivity (# sig) + specificity (features/cluster, clusters/feature) ----
    print(f"\n=== SENSITIVITY: significant (FDR<{FDR_THR}) enrichments per analysis x category ===")
    sens = sig.groupby(["run", "category"]).size().unstack(fill_value=0)
    # for per_study, also count distinct features (across studies) and # studies
    print(sens.to_string())
    rep = open(os.path.join(OUT, "SATAY_CONTRAST.txt"), "w")
    rep.write("SENSITIVITY (# FDR<0.05 enrichments) per run x category\n" + sens.to_string() + "\n\n")

    for run in ["per_study", "matched_sex_platform", "udon_core"]:
        s = sig[sig["run"].str.startswith(run.split(":")[0])] if run == "per_study" else sig[sig["run"] == run]
        if run == "per_study":
            s = sig[sig["run"].str.startswith("per_study")]
        print(f"\n=== {run} ===")
        for cat in ["mutation", "clinical", "celltype"]:
            sc = s[s["category"] == cat]
            feats = sc["feature"].nunique(); ncl = sc["cluster"].nunique()
            print(f"  {cat:9s}: {len(sc)} sig | {feats} distinct features | {ncl} clusters involved")
            rep.write(f"{run} {cat}: {len(sc)} sig, {feats} features, {ncl} clusters\n")
            # specificity: features enriched in only 1 cluster (specific) vs many
            if len(sc):
                per_feat = sc.groupby("feature")["cluster"].nunique()
                spec = int((per_feat == 1).sum()); promisc = int((per_feat >= 3).sum())
                print(f"             specificity: {spec} feature(s) in exactly 1 cluster, {promisc} in >=3 clusters")
                top = sc.sort_values("fdr").head(6)
                for _, r in top.iterrows():
                    print(f"             {r['feature']:28s} cl={r['cluster']:4s} OR={r['odds_ratio']:.1f} "
                          f"fdr={r['fdr']:.1e} ({r['overlap']}/{r['cluster_n']})")

    # ---- which mutations/clinical are recovered in which runs (feature x run) ----
    for cat in ["mutation", "clinical", "celltype"]:
        piv = (sig[sig["category"] == cat]
               .assign(runtag=lambda x: np.where(x["run"].str.startswith("per_study"), "per_study", x["run"]))
               .groupby(["feature", "runtag"]).size().unstack(fill_value=0))
        piv.to_csv(os.path.join(OUT, f"feature_by_run_{cat}.tsv"), sep="\t")
    rep.close()

    # ---- DEFAULT EXPORT: both SATAY-UDON heatmap layouts for every result ----
    print("\n[default export] rendering SATAY-UDON heatmaps (cluster-grouped) ...")
    import satay_heatmap, satay_heatmap_byfield
    satay_heatmap.generate_all(feat, defined)
    print("[default export] rendering SATAY-UDON heatmaps (metadata-field-grouped) ...")
    satay_heatmap_byfield.generate_all(feat, defined)
    print("\nwrote:", OUT)


if __name__ == "__main__":
    main()
