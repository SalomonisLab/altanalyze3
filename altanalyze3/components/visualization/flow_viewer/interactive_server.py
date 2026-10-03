#!/usr/bin/env python3
"""Interactive rna2flow viewer: UMAPs, label overlays, marker colouring, biaxial gating.

Binary arrays stream once to the browser; panning, colouring and gate drawing then happen
client-side. Gate statistics are computed server-side because they cross-tabulate a polygon
against every label set, which is the comparison the viewer exists for.
"""
from __future__ import annotations

import argparse, json, os
import numpy as np
from flask import Flask, Response, jsonify, render_template, request

app = Flask(__name__)
B = {"root": None, "man": None, "chan": None, "lab": {}, "emb": {}}


def man():
    return B["man"]


def channels():
    if B["chan"] is None:
        m = man()
        B["chan"] = np.fromfile(os.path.join(B["root"], "arrays", m["channels_file"]),
                                dtype=np.float32).reshape(m["n_events"], len(m["channels"]))
    return B["chan"]


def labels(name):
    if name not in B["lab"]:
        spec = man()["labels"][name]
        B["lab"][name] = np.fromfile(os.path.join(B["root"], "arrays", spec["file"]), np.int16)
    return B["lab"][name]


def embedding(name):
    if name not in B["emb"]:
        spec = man()["embeddings"][name]
        B["emb"][name] = np.fromfile(os.path.join(B["root"], "arrays", spec["file"]),
                                     np.float32).reshape(-1, 2)
    return B["emb"][name]


@app.route("/")
def index():
    return render_template("interactive.html")


@app.route("/api/manifest")
def api_manifest():
    m = dict(man())
    m["labels"] = {k: {kk: vv for kk, vv in v.items() if kk != "file"}
                   for k, v in m["labels"].items()}
    m["embeddings"] = {k: {kk: vv for kk, vv in v.items() if kk != "file"}
                       for k, v in m["embeddings"].items()}
    return jsonify(m)


@app.route("/api/embedding/<name>")
def api_embedding(name):
    return Response(embedding(name).tobytes(), mimetype="application/octet-stream")


@app.route("/api/labels/<name>")
def api_labels(name):
    return Response(labels(name).tobytes(), mimetype="application/octet-stream")


@app.route("/api/channel/<path:name>")
def api_channel(name):
    m = man()
    if name not in m["channels"]:
        return jsonify({"error": "unknown channel %r" % name}), 404
    return Response(np.ascontiguousarray(
        channels()[:, m["channels"].index(name)]).tobytes(), mimetype="application/octet-stream")


def _inside(poly, X, Y):
    """Even-odd ray casting, vectorized over all events."""
    p = np.asarray(poly, dtype=np.float64)
    inside = np.zeros(len(X), dtype=bool)
    j = len(p) - 1
    for i in range(len(p)):
        xi, yi = p[i]; xj, yj = p[j]
        cross = ((yi > Y) != (yj > Y)) & (X < (xj - xi) * (Y - yi) / (yj - yi + 1e-30) + xi)
        inside ^= cross
        j = i
    return inside


@app.route("/api/gate", methods=["POST"])
def api_gate():
    """A polygon on any two axes -> which events fall in, and how every label set divides them.

    This is the virtual flow experiment: draw a gate, see its composition under FlowSOM and
    under each transferred annotation at once.
    """
    q = request.get_json(force=True)
    poly = q["polygon"]
    space = q.get("space", "channel")
    if space == "embedding":
        E = embedding(q["embedding"]); X, Y = E[:, 0], E[:, 1]
    else:
        m = man(); C = channels()
        X = C[:, m["channels"].index(q["x"])]
        Y = C[:, m["channels"].index(q["y"])]
    sel = _inside(poly, X.astype(np.float64), Y.astype(np.float64))
    n = int(sel.sum())
    out = {"n_gated": n, "n_total": int(len(X)),
           "percent": round(100.0 * n / max(len(X), 1), 3), "composition": {}}
    if n:
        for ls in q.get("label_sets", []):
            lev = man()["labels"][ls]["levels"]
            c = labels(ls)[sel]
            vals, cnt = np.unique(c, return_counts=True)
            order = np.argsort(-cnt)
            out["composition"][ls] = [
                {"label": lev[vals[i]] if 0 <= vals[i] < len(lev) else "NA",
                 "n": int(cnt[i]), "percent": round(100.0 * cnt[i] / n, 2)}
                for i in order[:15]]
    return jsonify(out)


@app.route("/api/concordance", methods=["POST"])
def api_concordance():
    """Cross-tabulate two label sets over all events, or within a gate."""
    q = request.get_json(force=True)
    a, b = labels(q["a"]), labels(q["b"])
    la, lb = man()["labels"][q["a"]]["levels"], man()["labels"][q["b"]]["levels"]
    sel = slice(None)
    if q.get("polygon"):
        m = man(); C = channels()
        X = C[:, m["channels"].index(q["x"])]; Y = C[:, m["channels"].index(q["y"])]
        sel = _inside(q["polygon"], X.astype(np.float64), Y.astype(np.float64))
    M = np.zeros((len(la), len(lb)), dtype=np.int64)
    np.add.at(M, (a[sel].astype(int), b[sel].astype(int)), 1)
    return jsonify({"rows": la, "cols": lb, "counts": M.tolist()})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8078)
    a = ap.parse_args()
    B["root"] = os.path.abspath(a.bundle)
    B["man"] = json.load(open(os.path.join(B["root"], "manifest.json")))
    print("rna2flow interactive viewer: http://%s:%d" % (a.host, a.port))
    app.run(host=a.host, port=a.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
