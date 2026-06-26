#!/usr/bin/env python3
"""Standalone test of the ISV-web network query: reproduce a regulator ego-network for two contexts and
render a side-by-side radial preview (sanity-check the data layer that the D3 tab will consume)."""
import sys, os, json
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/isv_web")
import network as N
import numpy as np
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42; matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
GR = LinearSegmentedColormap.from_list("gr", ["#f7f7f7", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"])
ISO_COL = ["#E45756", "#4C78A8", "#54A24B", "#B279A2", "#EECA3B"]

ART = "/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM/isvweb_network.json"
net = N.load_network(ART)
genesrc = {e["src_gene"] for e in net["edges"]}
print("regulators (source genes) available:", len(genesrc), "| MAX:", "MAX" in genesrc, "DDIT3:", "DDIT3" in genesrc, "ATF2:", "ATF2" in genesrc)

GENE = "MAX" if "MAX" in genesrc else sorted(genesrc)[0]
CTXA, CTXB = "HSC-1||young", "MEP-Eryth-2||young"
res = N.query_network(net, CTXA, CTXB, threshold=1.0, genes=[GENE], edge_type=None)
for side in ("A", "B"):
    s = res[side]; print(f"[{s['ctx']}] {GENE}: {len(s['edges'])} edges, {len(s['nodes'])} nodes")
    for e in s["edges"][:4]:
        print(f"    {e['src_iso']} -[{e['type']}]-> {e['tgt_gene']}  src_cpm={e['src_cpm']} tgt_cpm={e['tgt_cpm']}")

# ---- side-by-side radial preview ----
def draw(ax, side, title):
    eds = side["edges"]
    if not eds: ax.text(0.5, 0.5, "no edges", ha="center"); ax.axis("off"); ax.set_title(title); return
    isos = sorted({e["src_iso"] for e in eds}); icol = {iso: ISO_COL[i % len(ISO_COL)] for i, iso in enumerate(isos)}
    tgts = sorted({e["tgt_gene"] for e in eds})
    nco = {n["id"]: n for n in side["nodes"]}
    vmax = max([n["cpm"] for n in side["nodes"]] + [1.0])
    # sources along the bottom, targets on an upper arc
    sx = {iso: (0.5 + (i - (len(isos) - 1) / 2) * 0.16) for i, iso in enumerate(isos)}
    import math
    tpos = {}
    for i, t in enumerate(tgts):
        ang = math.pi * (0.12 + 0.76 * (i / max(1, len(tgts) - 1)))
        tpos[t] = (0.5 + 0.42 * math.cos(ang), 0.30 + 0.62 * math.sin(ang))
    for e in eds:
        x0, y0 = sx[e["src_iso"]], 0.10; x1, y1 = tpos[e["tgt_gene"]]
        ax.plot([x0, x1], [y0, y1], color=icol[e["src_iso"]], lw=1.3, alpha=0.7, zorder=1, solid_capstyle="round")
        if e["type"] == "PDI":
            ax.annotate("", xy=(x1, y1), xytext=(x0 + (x1 - x0) * 0.82, y0 + (y1 - y0) * 0.82),
                        arrowprops=dict(arrowstyle="-|>", color=icol[e["src_iso"]], lw=1.2), zorder=2)
    for t in tgts:
        x, y = tpos[t]; cpm = nco[t]["cpm"]
        ax.scatter([x], [y], s=240, c=[GR(cpm / vmax)], edgecolors="#333", linewidths=0.5, zorder=3)
        ax.text(x, y + 0.045, t, ha="center", va="bottom", fontsize=7)
    for iso in isos:
        x = sx[iso]; cpm = nco[iso]["cpm"]
        ax.scatter([x], [0.10], s=420, c=[GR(cpm / vmax)], marker="D", edgecolors="#000", linewidths=0.8, zorder=4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title(title, fontsize=11)
    ax.legend(handles=[plt.Line2D([0], [0], color=icol[i], lw=2, label=i) for i in isos], fontsize=6, loc="lower center", ncol=min(3, len(isos)), frameon=False)

fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))
draw(axes[0], res["A"], f"{GENE} — {CTXA}"); draw(axes[1], res["B"], f"{GENE} — {CTXB}")
fig.suptitle(f"{GENE} interaction network (expression-supported, >1 CPM both endpoints) — {res['mode']} mode", fontsize=12)
out = "/tmp/isvweb_network_preview.pdf"; fig.savefig(out, bbox_inches="tight"); print("wrote", out)
