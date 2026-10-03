#!/usr/bin/env python3
"""Permute the normalization and assignment choices, scored on MarkerFinder concordance.

ARI is recorded but never used to pick a winner: it rewards agreeing with FlowSOM's particular
20-channel partition, and on this data it ranked the most collapsed method first.
"""
import argparse, itertools, os, sys, time
# xgboost ships its own OpenMP runtime. On macOS it clashes with the one numpy and sklearn
# have already loaded, and the process dies with no traceback, losing every row not yet
# written. n_jobs=1 inside the model is not enough: the limit has to be in the environment
# BEFORE numpy loads. Setting it here makes the module safe on its own, rather than relying
# on the caller remembering to export it. On 2026-10-02 a run without it lost 7 of 12 rows.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd, anndata as ad
from sklearn.neighbors import KNeighborsClassifier

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")
from altanalyze3.components.rna2flow.io import read_flowjo_rds
from altanalyze3.components.rna2flow.crosswalk import build_crosswalk
from altanalyze3.components.rna2flow.normalize import scale_feature, reference_spline, kde_quantile_map
from altanalyze3.components.rna2flow.evaluate import score_transfer
from altanalyze3.components.rna2flow.marker_concordance import concordance

FLOW_RDS = ("/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/Ungated_J8DW/"
            "FlowSOM_RData_J8DW.RDS")
FLOW_LAB = ("/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/Ungated_J8DW/"
            "FlowSOM_Results_J8DW.csv")


def _dense(X):
    """np.asarray on a scipy sparse matrix returns a 0-d object array, which silently
    transposes everything downstream. Densify explicitly."""
    import scipy.sparse as sp
    return np.asarray(X.toarray() if sp.issparse(X) else X, dtype=np.float64)


def log(m): print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def prep(C, F, clip, joint):
    """joint=True scales both platforms on one shared percentile range, so a marker that is
    genuinely dimmer on one platform stays dimmer instead of being stretched to match."""
    lo, hi = clip
    if not joint:
        Cs = np.column_stack([scale_feature(C[:, j], lo, hi) for j in range(C.shape[1])])
        Fs = np.column_stack([scale_feature(F[:, j], lo, hi) for j in range(F.shape[1])])
        return Cs, Fs
    Cs, Fs = np.empty_like(C), np.empty_like(F)
    for j in range(C.shape[1]):
        both = np.concatenate([C[:, j], F[:, j]])
        a, b = np.percentile(both, lo), np.percentile(both, hi)
        rng = (b - a) or 1.0
        Cs[:, j] = np.clip((C[:, j] - a) / rng, 0, 1)
        Fs[:, j] = np.clip((F[:, j] - a) / rng, 0, 1)
    return Cs, Fs


def kde_transfer(C, F, labels, clip, joint, n_ref, ties, k, seed=0):
    Cs, Fs = prep(C, F, clip, joint)
    rng = np.random.default_rng(seed)
    Cm = np.empty_like(Cs)
    for j in range(Cs.shape[1]):
        ref = Fs[:, j]
        if n_ref and ref.size > n_ref:
            ref = ref[rng.choice(ref.size, n_ref, replace=False)]
        Cm[:, j] = kde_quantile_map(Cs[:, j], None, spline=reference_spline(ref), ties=ties)
    return KNeighborsClassifier(n_neighbors=k, n_jobs=-1).fit(Cm, labels).predict(Fs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adt-h5ad", required=True)
    ap.add_argument("--label-col", default="celltype")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="sweep")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fl = read_flowjo_rds(FLOW_RDS)
    flowsom = pd.read_csv(FLOW_LAB).iloc[:, 1].astype(str).values
    adt = ad.read_h5ad(a.adt_h5ad)
    labels = adt.obs[a.label_col].astype(str).values
    pairs, _, _ = build_crosswalk(fl.antibodies, list(adt.var_names), extra_aliases={"CD8": "CD8a"})
    # Dedup on BOTH sides. One flow marker can match several CITE features (two clones of
    # one protein), and a duplicated flow name gives MarkerFinder duplicate feature names,
    # whose row lookup then returns several rows and breaks the rank comparison.
    pairs = (pairs.drop_duplicates(subset=["cite"], keep="first")
                  .drop_duplicates(subset=["flow"], keep="first").reset_index(drop=True))
    fi = {c: i for i, c in enumerate(fl.channels)}
    cols = [fi[[c for c in fl.channels
                if (c.split("::", 1)[1].strip() if "::" in c else c.strip()) == f][0]]
            for f in pairs["flow"]]
    F = np.asarray(fl.X[:, cols], dtype=np.float64)
    ci = {c: i for i, c in enumerate(adt.var_names)}
    C = _dense(adt.X)[:, [ci[c] for c in pairs["cite"]]]
    keep = labels != "unassigned"
    C, L = C[keep], labels[keep]
    log("%d shared markers, %d CITE cells, %d labels" % (len(pairs), len(C), len(set(L))))
    if len(pairs) < 3:
        # Running the grid anyway produces 36 identical failures and a TSV of error strings,
        # which reads like a result. Stop at the real cause instead.
        raise SystemExit(
            "only %d shared marker(s) between the flow panel and %s.\n"
            "  --adt-h5ad must hold the ADT PANEL (one column per antibody) AND the label\n"
            "  column. An h5ad carrying only the label column has no markers to match.\n"
            "  flow antibodies: %s\n"
            "  CITE features (first 10 of %d): %s"
            % (len(pairs), a.adt_h5ad, ", ".join(map(str, fl.antibodies[:8])),
               adt.shape[1], ", ".join(map(str, list(adt.var_names)[:10]))))

    grid = list(itertools.product([(1, 99), (0.5, 99.5), (5, 95)], [False, True],
                                  [20000, 0], ["first"], [5, 15, 50]))
    rows, tsv = [], os.path.join(a.out, "sweep_%s.tsv" % a.tag)
    for clip, joint, n_ref, ties, k in grid:
        t = time.time()
        try:
            pred = kde_transfer(C, F, L, clip, joint, n_ref, ties, k)
            s = score_transfer(pred, flowsom)
            _, cs = concordance(C, L, F, pred, pairs); s.update(cs)
            s.update(clip_lo=clip[0], clip_hi=clip[1], joint_scaling=joint,
                     kde_reference_n=n_ref or len(F), ties=ties, k=k,
                     seconds=round(time.time() - t, 1))
            rows.append(s)
            log("clip=%-10s joint=%-5s ref=%-6s k=%-3d -> pass %2d/%-2d rho=%.3f ARI=%.3f top=%.2f"
                % (str(clip), joint, n_ref or "all", k, s["populations_passing"],
                   s["populations_compared"], s["median_rank_spearman"], s["ARI"],
                   s["largest_label_fraction"]))
        except Exception as e:
            log("clip=%s joint=%s ref=%s k=%d FAILED %s" % (clip, joint, n_ref, k, str(e)[:120]))
            rows.append(dict(clip_lo=clip[0], clip_hi=clip[1], joint_scaling=joint,
                             kde_reference_n=n_ref, k=k, error=str(e)[:200]))
        pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
    log("wrote %s" % tsv)


if __name__ == "__main__":
    main()
