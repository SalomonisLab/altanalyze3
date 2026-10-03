#!/usr/bin/env python3
"""Generalized CITE-seq -> flow benchmark: any dataset, any label sets, any method set.

Primary criterion is MarkerFinder concordance (populations passing), not ARI.
"""
import argparse, json, os, sys, time
# xgboost ships its own OpenMP runtime. On macOS it clashes with the one numpy and sklearn
# have already loaded, and the process dies with no traceback, losing every row not yet
# written. n_jobs=1 inside the model is not enough: the limit has to be in the environment
# BEFORE numpy loads. Setting it here makes the module safe on its own, rather than relying
# on the caller remembering to export it. On 2026-10-02 a run without it lost 7 of 12 rows.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")
from altanalyze3.components.rna2flow.io import read_flow, align_event_labels
from altanalyze3.components.rna2flow.crosswalk import build_crosswalk
from altanalyze3.components.rna2flow.transfer import METHODS
from altanalyze3.components.rna2flow.evaluate import score_transfer, adversarial_controls
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


def to_annot(n):
    """cellHarmony <bc>-1.<lib> -> <lib>_<bc>-1 ; anything else passes through."""
    return "%s_%s" % (n.rsplit(".", 1)[1], n.rsplit(".", 1)[0]) if "." in n else n


def load_labels(adt, spec):
    """spec: name=path:obscol  (h5ad) or name=obs:col (from the ADT object itself)."""
    out = {}
    for s in spec:
        name, src = s.split("=", 1)
        if src.startswith("obs:"):
            out[name] = adt.obs[src[4:]].astype(str).values
        else:
            path, col = src.rsplit(":", 1)
            o = ad.read_h5ad(path, backed="r")
            vals = o.obs[col].astype(str).values
            # Barcode conventions differ between objects. Try the raw names and the
            # cellHarmony <bc>.<lib> <-> <lib>_<bc> conversion, and keep whichever actually
            # overlaps. Guessing one and silently joining 0 cells is the failure this avoids.
            best, best_n = None, -1
            for nm, idx in (("raw", pd.Index([str(n) for n in o.obs_names])),
                            ("converted", pd.Index([to_annot(str(n)) for n in o.obs_names]))):
                ser = pd.Series(vals, index=idx)
                ser = ser[~ser.index.duplicated()]
                n = len(set(ser.index) & set(adt.obs_names))
                log("    %s join via %-9s -> %d of %d cells" % (name, nm, n, adt.n_obs))
                if n > best_n:
                    best, best_n = ser, n
            if best_n == 0:
                raise ValueError("label set %r joins 0 cells under either barcode convention" % name)
            out[name] = best.reindex(adt.obs_names).fillna("unassigned").values
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adt-h5ad", required=True)
    ap.add_argument('--flow-input', default=FLOW_RDS, help='FCS, FlowJo RDS, or event × channel CSV')
    ap.add_argument('--flow-format', choices=['fcs','rds','csv'])
    ap.add_argument('--flow-labels', default=FLOW_LAB)
    ap.add_argument('--flowsom-column')
    ap.add_argument('--event-id-column', default='EventNumberDP')
    ap.add_argument("--label", action="append", required=True, help="name=path:obscol | name=obs:col")
    ap.add_argument("--out", required=True)
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--tag", default="run")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fl = read_flow(a.flow_input, a.flow_format)
    flowsom = align_event_labels(fl, pd.read_csv(a.flow_labels), a.flowsom_column, a.event_id_column)
    log("flow: %d events x %d channels, %d metaclusters" % (*fl.X.shape, len(set(flowsom))))

    adt = ad.read_h5ad(a.adt_h5ad)
    log("CITE: %s from %s" % (adt.shape, os.path.basename(a.adt_h5ad)))
    label_sets = load_labels(adt, a.label)
    for k, v in label_sets.items():
        log("  %-18s %d labels, %d unassigned" % (k, len(set(v)), int((v == "unassigned").sum())))

    pairs, unf, _ = build_crosswalk(fl.antibodies, list(adt.var_names),
                                    extra_aliases={"CD8": "CD8a"})
    # Dedup on BOTH sides. One flow marker can match several CITE features (two clones of
    # one protein), and a duplicated flow name gives MarkerFinder duplicate feature names,
    # whose row lookup then returns several rows and breaks the rank comparison.
    pairs = (pairs.drop_duplicates(subset=["cite"], keep="first")
                  .drop_duplicates(subset=["flow"], keep="first").reset_index(drop=True))
    log("%d shared markers; unmatched flow: %s" % (len(pairs), ", ".join(unf)))
    pairs.to_csv(os.path.join(a.out, "crosswalk_%s.tsv" % a.tag), sep="\t", index=False)
    if len(pairs) < 3:
        log("FEWER THAN 3 SHARED MARKERS; nothing to transfer"); return

    fi = {c: i for i, c in enumerate(fl.channels)}
    cols = []
    for f in pairs["flow"]:
        hit = [c for c in fl.channels if (c.split("::", 1)[1].strip() if "::" in c else c.strip()) == f]
        cols.append(fi[hit[0]])
    Fx = np.asarray(fl.X[:, cols], dtype=np.float64)
    ci = {c: i for i, c in enumerate(adt.var_names)}
    Cx = _dense(adt[:, list(pairs['cite'])].X)

    rows, tsv = [], os.path.join(a.out, "transfer_benchmark_%s.tsv" % a.tag)
    for lname, labels in label_sets.items():
        keep = labels != "unassigned"
        C, L = Cx[keep], labels[keep]
        log("== %s: %d CITE cells, %d labels" % (lname, keep.sum(), len(set(L))))
        for mname in a.methods.split(","):
            fn = METHODS[mname]; t = time.time()
            try:
                pred = fn(C, Fx, L)
                s = score_transfer(pred, flowsom)
                s.update(adversarial_controls(fn, C, Fx, L, flowsom))
                cdf, csum = concordance(C, L, Fx, pred, pairs); s.update(csum)
                cdf.to_csv(os.path.join(a.out, "concordance_%s_%s_%s.tsv" % (a.tag, lname, mname)),
                           sep="\t", index=False)
                s.update(label_set=lname, method=mname, seconds=round(time.time() - t, 1),
                         n_markers=len(pairs), n_cite=int(keep.sum()), dataset=a.tag)
                rows.append(s)
                log("   %-16s ARI=%.4f | pass %d/%d rho=%.3f %s | top=%.2f %.0fs"
                    % (mname, s["ARI"], s["populations_passing"], s["populations_compared"],
                       s["median_rank_spearman"], "MEETS" if s["meets_target"] else "",
                       s["largest_label_fraction"], s["seconds"]))
                np.save(os.path.join(a.out, "pred_%s_%s_%s.npy" % (a.tag, lname, mname)), pred.astype(str))
            except Exception as e:
                log("   %-16s FAILED: %s" % (mname, str(e)[:150]))
                rows.append(dict(dataset=a.tag, label_set=lname, method=mname, error=str(e)[:300]))
            pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
    log("wrote %s" % tsv)


if __name__ == "__main__":
    main()
