"""Single-cell evidence map for one variant: one column per cell, one row per inference marker.

Columns are grouped in blocks (for example: cells with a mutant read, cells with reads and no mutant read, cells with no
read that a marker infers as mutant). Each row is one marker of the variant (a passenger variant, a copy-number clone).
A cell-by-row tile is 'present', 'absent' (the cell was assessed and lacks the marker) or 'not assessed' (no read over
the marker). The reader sees, cell by cell, whether a marker stays inside the read-mutant block and out of the
read-wild-type block, and which cells each marker adds.

Input: a table with one row per cell, a block column, and one column per track holding 1 (present), 0 (absent) or
-1 (not assessed).  Illustrator-ready vector PDF: live Arial text, one compound path per tile colour, no clipping.

CLI:  python -m altanalyze3.components.visualization.variant_cell_evidence_map --cells CELLS.tsv --tracks TRACKS.tsv
          --out OUT.pdf [--title T] [--legend-file LEGEND.txt]
  TRACKS.tsv columns: key (column of CELLS.tsv), label, note (free text shown at the right, may be empty)
"""
from __future__ import annotations

import argparse
import textwrap

import numpy as np
import pandas as pd

C_PRESENT, C_ABSENT, C_NOT_ASSESSED, C_FRAME, C_TEXT = "#B2182B", "#D0D0D0", "#FFFFFF", "#4D4D4D", "#000000"
BLOCK_HEAD = ["#B2182B", "#4D4D4D", "#F08C7D", "#9AC2E5", "#8C6BB1"]


def _tiles(ax, x0, w, y0, h, colour):
    """One compound path of axis-aligned rectangles (x0, w arrays; fixed y0, h)."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    n = len(x0)
    if n == 0:
        return
    v = np.zeros((n, 5, 2)); c = np.tile(np.array([Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY], dtype=np.uint8), n)
    v[:, 0] = np.c_[x0, np.full(n, y0)]; v[:, 1] = np.c_[x0 + w, np.full(n, y0)]
    v[:, 2] = np.c_[x0 + w, np.full(n, y0 + h)]; v[:, 3] = np.c_[x0, np.full(n, y0 + h)]; v[:, 4] = v[:, 0]
    p = PathPatch(Path(v.reshape(-1, 2), c), facecolor=colour, edgecolor="none", linewidth=0)
    p.set_clip_on(False)
    ax.add_patch(p)


def block_counts(cells: pd.DataFrame, tracks, block_col="block", block_order=None) -> pd.DataFrame:
    """present / assessed / n per track and block (the numbers printed on the figure)."""
    order = block_order or list(pd.unique(cells[block_col]))
    out = []
    for t in tracks:
        for b in order:
            v = cells.loc[cells[block_col] == b, t["key"]].to_numpy()
            out.append(dict(track=t["label"].replace("\n", " "), block=b, n_cells=len(v), present=int((v == 1).sum()),
                            absent=int((v == 0).sum()), not_assessed=int((v == -1).sum())))
    return pd.DataFrame(out)


def plot_variant_cell_evidence_map(cells: pd.DataFrame, tracks, out_pdf, title="", legend="", block_col="block",
                                   block_order=None, min_block_share=0.16, figsize_width=13.0,
                                   block_short=None):
    """Draw the map. `tracks`: list of {key, label, note}. Returns the counts table of block_counts()."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    order = block_order or list(pd.unique(cells[block_col]))
    keys = [t["key"] for t in tracks]
    bad = set(np.unique(cells[keys].to_numpy())) - {-1, 0, 1}
    if bad:
        raise ValueError(f"track values must be 1, 0 or -1; found {sorted(bad)}")
    if len(order) > len(BLOCK_HEAD):
        raise ValueError("at most %d blocks" % len(BLOCK_HEAD))
    counts = block_counts(cells, tracks, block_col, order)
    n_by = [int((cells[block_col] == b).sum()) for b in order]
    share = np.sqrt(np.maximum(n_by, 1)).astype(float); share = np.maximum(share / share.sum(), min_block_share); share = share / share.sum()
    gap = 0.012; span = 1.0 - gap * (len(order) - 1)
    x_start = np.r_[0, np.cumsum(share * span + gap)[:-1]]

    nt = len(tracks); row_h = 0.62
    short = block_short or {}
    leg_lines = textwrap.wrap(legend, 245) if legend else []
    fig_h = 1.5 + row_h * nt + 0.115 * len(leg_lines) + 0.5
    fig = plt.figure(figsize=(figsize_width, fig_h))
    top_in, bot_in = 1.0, 0.115 * len(leg_lines) + 0.75
    ax = fig.add_axes([0.2, bot_in / fig_h, 0.5, (fig_h - top_in - bot_in) / fig_h])
    ax.set_xlim(0, 1); ax.set_ylim(nt, -1.0); ax.axis("off")

    for bi, b in enumerate(order):
        sub = cells[cells[block_col] == b]
        if len(sub):
            m = sub[keys].to_numpy()
            # most markers present first, then by the track order
            o = np.lexsort([-(m[:, j] == 1).astype(int) for j in range(nt - 1, -1, -1)] + [-(m == 1).sum(1)])
            m = m[o]
            w = share[bi] * span / len(sub); xs = x_start[bi] + np.arange(len(sub)) * w
            for j in range(nt):
                for val, col in ((1, C_PRESENT), (0, C_ABSENT)):          # not assessed: left empty inside the frame
                    _tiles(ax, xs[m[:, j] == val], w, j + 0.12, 0.76, col)
        x0, x1 = x_start[bi], x_start[bi] + share[bi] * span
        _tiles(ax, np.array([x0]), x1 - x0, -0.62, 0.3, BLOCK_HEAD[bi])
        for j in range(nt):                                   # frame of each row of the block: four 2-point lines
            for xa, ya, xb, yb in ((x0, j + 0.12, x1, j + 0.12), (x0, j + 0.88, x1, j + 0.88), (x0, j + 0.12, x0, j + 0.88), (x1, j + 0.12, x1, j + 0.88)):
                ax.plot([xa, xb], [ya, yb], color=C_FRAME, linewidth=0.4, clip_on=False)
        ax.text((x0 + x1) / 2, -0.72, "%s\n(n = %d cells)" % (b, n_by[bi]), ha="center", va="bottom", fontsize=7, color=C_TEXT)
    for j, t in enumerate(tracks):
        ax.text(-0.012, j + 0.5, t["label"], ha="right", va="center", fontsize=7, color=C_TEXT)
        c = counts[counts.track == t["label"].replace("\n", " ")]
        lines = []
        for r in c.itertuples():
            s = "%s: %d of %d" % (short.get(r.block, r.block.replace("\n", " ")), r.present, r.n_cells)
            if r.not_assessed:
                s += " (%d of %d with reads over the marker)" % (r.present, r.present + r.absent)
            lines.append(s)
        if t.get("note"):
            lines.append(str(t["note"]))
        ax.text(1.012, j + 0.5, "\n".join(lines), ha="left", va="center", fontsize=6, color=C_TEXT)

    from matplotlib.patches import Rectangle
    hs = [Rectangle((0, 0), 1, 1, facecolor=C_PRESENT, edgecolor="none"), Rectangle((0, 0), 1, 1, facecolor=C_ABSENT, edgecolor="none"),
          Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=C_FRAME, linewidth=0.4)]
    lg = fig.legend(hs, ["marker present in the cell", "cell assessed, marker absent", "no read over the marker in the cell"],
                    loc="lower left", bbox_to_anchor=(0.2, (bot_in - 0.42) / fig_h), ncol=3, frameon=False, fontsize=7)
    fig.text(0.01, 1 - 0.28 / fig_h, title, fontsize=8, ha="left", va="top", color=C_TEXT)
    if leg_lines:
        fig.text(0.01, 0.06 / fig_h, "\n".join(leg_lines), fontsize=6, ha="left", va="bottom", color=C_TEXT)
    fig.patch.set_visible(False)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)
    return counts


def main(argv=None):
    ap = argparse.ArgumentParser(description="Single-cell evidence map for one variant.")
    ap.add_argument("--cells", required=True); ap.add_argument("--tracks", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=""); ap.add_argument("--legend-file", default=None); ap.add_argument("--block-col", default="block")
    a = ap.parse_args(argv)
    cells = pd.read_csv(a.cells, sep="\t"); tr = pd.read_csv(a.tracks, sep="\t").fillna("")
    legend = open(a.legend_file).read().strip() if a.legend_file else ""
    c = plot_variant_cell_evidence_map(cells, tr.to_dict("records"), a.out, a.title, legend, a.block_col)
    c.to_csv(a.out.rsplit(".", 1)[0] + "__counts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
