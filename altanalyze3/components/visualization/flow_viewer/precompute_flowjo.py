#!/usr/bin/env python3
"""Load the authors' FlowJo gating strategy into a flow_viewer bundle.

Writes three things beside the existing arrays:

  flowjo_gates.json   the gate TREE as the authors drew it: every population, its parent,
                      its gate type, its two channels, its vertices in RAW space and the
                      SAME vertices mapped through this bundle's logicle, plus the count
                      this run recomputed and the count the workspace recorded.
  lab_FlowJo_leaf     int16 codes: each event labelled by the DEEPEST manual population
                      that contains it, so the manual strategy is a label set like any other.
  one boolean array   per population, so a population can be shown or used as a parent
                      without re-running the tree in the browser.

Gates are applied in RAW space, which is what their coordinates mean. The display mapping
exists only so an outline lands on its own events on screen.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")
from altanalyze3.components.rna2flow.flowjo import (
    parse_workspace, apply_tree, population_table, build_lookup)
from altanalyze3.components.rna2flow.io import read_fcs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--workspace", required=True, help=".wsp or .wspt holding the gate tree")
    ap.add_argument("--fcs", required=True)
    ap.add_argument("--space", default="flow")
    a = ap.parse_args()
    root = os.path.abspath(a.bundle); arr = os.path.join(root, "arrays")
    man = json.load(open(os.path.join(root, "manifest.json")))
    sp = man["spaces"][a.space]

    pops, transforms = parse_workspace(a.workspace)
    print("workspace: %s" % os.path.basename(a.workspace))
    print("  populations %d | gates %s" % (
        len(pops), {k: sum(1 for p in pops if p.gate and p.gate.kind == k)
                    for k in ("polygon", "rectangle", "ellipsoid")}))

    raw = read_fcs(a.fcs)
    reg, dup = build_lookup(raw.channels, raw.antibodies)
    need = set()
    for p in pops:
        if p.gate:
            need.update(p.gate.dims)
    missing = sorted(need - set(reg))
    print("  events %d | gate channels %d | unresolved %s | ambiguous %s"
          % (raw.X.shape[0], len(need), missing or "none", dup or "none"))

    applied = apply_tree(pops, {k: raw.X[:, i] for k, i in reg.items()})
    rows = population_table(pops, applied, raw.X.shape[0])
    uneval = [r for r in rows if r["unevaluable"]]
    print("  unevaluable %d of %d" % (len(uneval), len(rows)))

    # ---- display mapping: the viewer's own logicle, per channel ----------------------
    from altanalyze3.components.rna2flow.transform import read_fcs_anndata, logicle
    adv = read_fcs_anndata(a.fcs).var
    need_cols = ["LOGICLE_T", "LOGICLE_W", "LOGICLE_M", "LOGICLE_A"]
    missing_cols = [c for c in need_cols if c not in adv.columns]
    if missing_cols:
        # Falling back to raw coordinates here would write plausible-looking numbers in the
        # wrong space, and the gate would be drawn off-screen with no error. Refuse instead.
        raise SystemExit("FCS var lacks %s; cannot map gates into display space"
                         % ", ".join(missing_cols))
    det2ab = {str(d).split(" :: ")[0].strip(): ab for d, ab in zip(raw.channels, raw.antibodies)}
    ab2det = {v: k for k, v in det2ab.items()}

    def disp(chan, vals):
        """Map raw coordinates on `chan` into the bundle's display space (logicle)."""
        key = chan if chan in adv.index else ab2det.get(chan)
        if key is None or key not in adv.index:
            raise SystemExit("channel %r absent from the FCS var table" % chan)
        r = adv.loc[key]
        out = logicle(np.asarray(vals, float), T=float(r["LOGICLE_T"]),
                      W=float(r["LOGICLE_W"]), M=float(r["LOGICLE_M"]),
                      A=float(r["LOGICLE_A"]))
        out = np.asarray(out, float)
        if not np.isfinite(out[np.isfinite(np.asarray(vals, float))]).all():
            raise SystemExit("logicle produced a non-finite coordinate on %r" % chan)
        return [float(v) for v in out]

    tree, codes_leaf = [], np.full(raw.X.shape[0], -1, np.int16)
    leaf_levels = []
    for r, p in zip(rows, pops):
        g = p.gate
        node = dict(r)
        node["channels"] = [det2ab.get(c, c) for c in (g.dims if g else [])]
        node["raw_channels"] = list(g.dims) if g else []
        if g and g.kind == "polygon":
            xs = [v[0] for v in g.vertices]; ys = [v[1] for v in g.vertices]
            node["vertices_raw"] = [[float(x), float(y)] for x, y in zip(xs, ys)]
            node["vertices_display"] = [[x, y] for x, y in
                                        zip(disp(g.dims[0], xs), disp(g.dims[1], ys))]
        elif g and g.kind == "rectangle":
            node["bounds_raw"] = [[float(lo), float(hi)] for lo, hi in g.bounds]
            node["bounds_display"] = [disp(c, [lo, hi]) for c, (lo, hi) in
                                      zip(g.dims, g.bounds)]
        tree.append(node)
        f = "flowjo_%03d.bool" % len(tree)
        applied[p.path]["mask"].astype(np.uint8).tofile(os.path.join(arr, f))
        node["mask_file"] = f

    order = sorted(range(len(pops)), key=lambda i: -rows[i]["depth"])
    for i in order:                       # deepest wins, so a leaf is not overwritten
        m = applied[pops[i].path]["mask"]
        tgt = m & (codes_leaf < 0)
        if tgt.any():
            codes_leaf[tgt] = len(leaf_levels)
            leaf_levels.append(pops[i].path)
    codes_leaf.tofile(os.path.join(arr, "lab_FlowJo_leaf.i16"))
    sp["labels"]["FlowJo_manual_gates"] = {
        "file": "lab_FlowJo_leaf.i16", "levels": leaf_levels, "source": "FlowJo workspace"}

    man.setdefault("gatesets", {})["FlowJo_manual"] = {
        "space": a.space, "workspace": os.path.abspath(a.workspace),
        "fcs": os.path.abspath(a.fcs), "n_events": int(raw.X.shape[0]),
        "transforms_in_workspace": len(transforms), "nodes": tree}
    # ---- invariant: a gate must intersect the data range of the channels it names -----
    sp_feats = list(sp["features"])
    Xd = np.memmap(os.path.join(arr, sp["features_file"]), np.float32, "r").reshape(
        int(sp["n"]), len(sp_feats))
    off = []
    for node in tree:
        box = node.get("bounds_display") or (
            [[min(p[0] for p in node["vertices_display"]),
              max(p[0] for p in node["vertices_display"])],
             [min(p[1] for p in node["vertices_display"]),
              max(p[1] for p in node["vertices_display"])]]
            if node.get("vertices_display") else None)
        if not box:
            continue
        for ch, (lo, hi) in zip(node["channels"], box):
            if ch not in sp_feats:
                continue
            v = Xd[:, sp_feats.index(ch)]
            if hi < float(v.min()) or lo > float(v.max()):
                off.append("%s on %s: gate [%.3g, %.3g] vs data [%.3g, %.3g]"
                           % (node["name"], ch, lo, hi, v.min(), v.max()))
    print("  gates outside their channel's data range: %d of %d%s"
          % (len(off), len(tree), ("" if not off else "\n    " + "\n    ".join(off[:6]))))

    json.dump(man, open(os.path.join(root, "manifest.json"), "w"), indent=1)
    pd.DataFrame(rows).to_csv(os.path.join(root, "flowjo_population_table.tsv"),
                              sep="\t", index=False)
    print("  leaf label set: %d populations covering %d of %d events (%.1f%%)"
          % (len(leaf_levels), int((codes_leaf >= 0).sum()), len(codes_leaf),
             100.0 * (codes_leaf >= 0).mean()))
    print("wrote %s/flowjo_population_table.tsv" % root)


if __name__ == "__main__":
    main()
