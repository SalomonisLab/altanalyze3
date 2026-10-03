#!/usr/bin/env python3
"""Register MarkerFinder-derived virtual gates in a bundle, in the manual gates' schema.

A virtual gate is a 2-channel box that rna2flow.gating derived for one transferred
population: MarkerFinder ranks the channels, coordinate ascent fits the box. Writing it in
the same shape as the FlowJo tree means the viewer draws both with one code path, and the
comparison is between two objects of the same kind rather than between a gate and a
clustering.

Each virtual gate is also SCORED here against the manual tree it is meant to rival, so the
viewer never shows a derived gate without the F1 of its closest manual population.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2flow.gate_compare import pairwise_f1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--gates", required=True, help="gates_*.json from rna2flow.gating")
    ap.add_argument("--name", required=True)
    ap.add_argument("--space", default="flow")
    ap.add_argument("--compare-to", default="FlowJo_manual")
    a = ap.parse_args()
    root = os.path.abspath(a.bundle); arr = os.path.join(root, "arrays")
    man = json.load(open(os.path.join(root, "manifest.json")))
    sp = man["spaces"][a.space]
    n = int(sp["n"]); feats = list(sp["features"])
    X = np.memmap(os.path.join(arr, sp["features_file"]), np.float32, "r").reshape(n, len(feats))

    gates = json.load(open(a.gates))
    ref = {}
    if a.compare_to in man.get("gatesets", {}):
        for nd in man["gatesets"][a.compare_to]["nodes"]:
            m = np.fromfile(os.path.join(arr, nd["mask_file"]), np.uint8).astype(bool)
            if m.sum() >= 20:
                ref[nd["path"]] = m

    nodes, masks = [], {}
    for i, g in enumerate(gates):
        xi, yi = feats.index(g["x_channel"]), feats.index(g["y_channel"])
        x0, x1, y0, y1 = g["gate"]
        m = ((X[:, xi] >= x0) & (X[:, xi] <= x1) & (X[:, yi] >= y0) & (X[:, yi] <= y1))
        f = "virtual_%s_%03d.bool" % (a.name, i)
        m.astype(np.uint8).tofile(os.path.join(arr, f))
        masks[g["population"]] = m
        nodes.append(dict(
            path=g["population"], name=g["population"], depth=0, gate="rectangle",
            channels=[g["x_channel"], g["y_channel"]], raw_channels=[],
            dims="%s|%s" % (g["x_channel"], g["y_channel"]),
            bounds_display=[[x0, x1], [y0, y1]],
            vertices_display=g.get("polygon"),
            n=int(m.sum()), pct_of_parent=round(100.0 * m.sum() / n, 3),
            pct_of_total=round(100.0 * m.sum() / n, 3), workspace_count=None,
            unevaluable=None, mask_file=f,
            f1_vs_own_population=g.get("f1"), f1_initial=g.get("f1_initial"),
            precision=g.get("precision"), recall=g.get("recall"),
            rho_x=g.get("rho_x"), rho_y=g.get("rho_y"),
            n_population=g.get("n_population")))

    if ref:
        F = pairwise_f1(ref, masks)
        for nd in nodes:
            s = F[nd["path"]]
            j = int(np.argmax(s.values))
            nd["best_manual_match"] = str(s.index[j])
            nd["f1_vs_manual"] = round(float(s.values[j]), 4)
        rows = [{"derived": nd["path"], "n": nd["n"],
                 "f1_vs_own_transferred_population": nd["f1_vs_own_population"],
                 "best_manual_match": nd["best_manual_match"],
                 "f1_vs_manual": nd["f1_vs_manual"]} for nd in nodes]
        os.makedirs(os.path.join(root, "comparison"), exist_ok=True)
        out = os.path.join(root, "comparison", "virtual_vs_manual_%s.tsv" % a.name)
        pd.DataFrame(rows).sort_values("f1_vs_manual", ascending=False).to_csv(
            out, sep="\t", index=False)
        print("scored %d virtual gates against %d manual populations -> %s"
              % (len(nodes), len(ref), out))
        print("  median F1 vs manual: %.3f | >=0.5: %d of %d"
              % (float(np.median([nd["f1_vs_manual"] for nd in nodes])),
                 sum(1 for nd in nodes if nd["f1_vs_manual"] >= 0.5), len(nodes)))

    man.setdefault("gatesets", {})[a.name] = {
        "space": a.space, "kind": "virtual", "source": os.path.abspath(a.gates),
        "n_events": n, "nodes": nodes}
    json.dump(man, open(os.path.join(root, "manifest.json"), "w"), indent=1)
    print("gate sets now: %s" % ", ".join(man["gatesets"]))


if __name__ == "__main__":
    main()
