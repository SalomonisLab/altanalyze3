#!/usr/bin/env python3
"""rna2flow viewer: a Flask app over a precomputed bundle, in the shape of isv_web.

It serves what a reader needs to judge a CITE-seq -> flow transfer:
  * the benchmark table, every method x label set with its adversarial controls
  * a contingency heatmap, transferred label vs FlowSOM metacluster
  * per-marker distribution overlays, CITE before and after mapping against flow

Nothing is computed here. `precompute.py` writes the bundle; this only serves it.
"""
from __future__ import annotations

import argparse, json, os
import numpy as np
import pandas as pd
from flask import Flask, jsonify, render_template, request

app = Flask(__name__)
BUNDLE = {"root": None}


def _load(name):
    p = os.path.join(BUNDLE["root"], name)
    if not os.path.exists(p):
        return None
    if p.endswith(".tsv"):
        return pd.read_csv(p, sep="\t")
    if p.endswith(".json"):
        return json.load(open(p))
    if p.endswith(".npy"):
        return np.load(p, allow_pickle=True)
    return None


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/benchmark")
def api_benchmark():
    df = _load("transfer_benchmark.tsv")
    if df is None:
        return jsonify({"error": "no benchmark in bundle"}), 404
    return jsonify(df.fillna("").to_dict(orient="records"))


@app.route("/api/crosswalk")
def api_crosswalk():
    df = _load("crosswalk_chinese.tsv")
    return jsonify([] if df is None else df.to_dict(orient="records"))


@app.route("/api/contingency")
def api_contingency():
    """Transferred label x FlowSOM metacluster, row-normalized."""
    ls, m = request.args.get("label_set"), request.args.get("method")
    pred = _load("pred_%s_%s.npy" % (ls, m))
    meta = _load("flowsom_labels.npy")
    if pred is None or meta is None:
        return jsonify({"error": "missing pred_%s_%s.npy" % (ls, m)}), 404
    ct = pd.crosstab(pd.Series(pred.astype(str), name="transferred"),
                     pd.Series(meta.astype(str), name="flowsom"))
    frac = ct.div(ct.sum(axis=1).replace(0, 1), axis=0)
    return jsonify({"rows": list(ct.index), "cols": list(ct.columns),
                    "counts": ct.values.tolist(), "fraction": frac.round(4).values.tolist()})


@app.route("/api/marker")
def api_marker():
    """Distributions for one shared marker: CITE raw, CITE mapped, flow."""
    d = _load("marker_distributions.json")
    if d is None:
        return jsonify({"error": "no marker distributions in bundle"}), 404
    mk = request.args.get("marker")
    return jsonify(d.get(mk, d) if mk else d)


def main():
    ap = argparse.ArgumentParser(description="Serve an rna2flow concordance bundle.")
    ap.add_argument("--bundle", required=True, help="directory written by precompute.py")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8077)
    a = ap.parse_args()
    BUNDLE["root"] = os.path.abspath(a.bundle)
    print("serving %s on http://%s:%d" % (BUNDLE["root"], a.host, a.port))
    app.run(host=a.host, port=a.port, debug=False)


if __name__ == "__main__":
    main()
