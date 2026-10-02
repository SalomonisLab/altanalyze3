#!/usr/bin/env python3
"""Assemble an rna2flow viewer bundle from a benchmark run.

Copies the benchmark table, the per-method predictions and the FlowSOM labels into one
directory, and precomputes the per-marker distribution overlays the viewer draws.
"""
import argparse, json, os, shutil, sys
import numpy as np, pandas as pd

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")


def hist(x, bins=60, lo=None, hi=None):
    x = np.asarray(x, dtype=float)
    lo = np.percentile(x, 0.5) if lo is None else lo
    hi = np.percentile(x, 99.5) if hi is None else hi
    c, e = np.histogram(np.clip(x, lo, hi), bins=bins, range=(lo, hi))
    return {"centers": ((e[:-1] + e[1:]) / 2).round(5).tolist(),
            "density": (c / max(c.sum(), 1)).round(6).tolist()}


def main():
    ap = argparse.ArgumentParser(description="Build an rna2flow viewer bundle.")
    ap.add_argument("--run-dir", required=True, help="directory written by run_chinese.py")
    ap.add_argument("--flowsom-csv", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    for f in os.listdir(a.run_dir):
        if f.endswith((".tsv", ".npy")):
            shutil.copy2(os.path.join(a.run_dir, f), os.path.join(a.out, f))
    lab = pd.read_csv(a.flowsom_csv)
    np.save(os.path.join(a.out, "flowsom_labels.npy"), lab.iloc[:, 1].astype(str).values)

    cw = os.path.join(a.out, "crosswalk_chinese.tsv")
    marker_json = {}
    if os.path.exists(cw):
        from altanalyze3.components.rna2flow.io import read_flowjo_rds
        from altanalyze3.components.rna2flow.normalize import scale_feature, kde_quantile_map
        import anndata as ad
        pairs = pd.read_csv(cw, sep="\t")
        fl = read_flowjo_rds("/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/"
                             "Ungated_J8DW/FlowSOM_RData_J8DW.RDS")
        adt = ad.read_h5ad("/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/Chinese/"
                           "PRJCA039526_thymus_ADT_dsb_annotated.h5ad")
        fi = {(c.split("::", 1)[1].strip() if "::" in c else c.strip()): i
              for i, c in enumerate(fl.channels)}
        ci = {c: i for i, c in enumerate(adt.var_names)}
        rng = np.random.default_rng(0)
        for _, r in pairs.iterrows():
            f = scale_feature(fl.X[:, fi[r["flow"]]], 1, 99)
            c = scale_feature(np.asarray(adt.X)[:, ci[r["cite"]]], 1, 99)
            ref = f[rng.choice(f.size, min(20000, f.size), replace=False)]
            marker_json[r["flow"]] = {"cite_raw": hist(c, lo=0, hi=1),
                                      "cite_mapped": hist(kde_quantile_map(c, ref), lo=0, hi=1),
                                      "flow": hist(f, lo=0, hi=1), "cite_feature": r["cite"]}
    json.dump(marker_json, open(os.path.join(a.out, "marker_distributions.json"), "w"))
    print("bundle written to %s (%d markers)" % (a.out, len(marker_json)))


if __name__ == "__main__":
    main()
