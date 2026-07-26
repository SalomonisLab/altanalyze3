#!/usr/bin/env python3
"""ISV-web-faithful isoform structure-track renderer (matplotlib), as a reusable iso2function component.

Reuses the AltAnalyze3 ``isoform_structure_view`` engine (read_gene_model / gene_model.db,
build_gene_maps, map_coord, strip_terminal_coords, build_isoform_segments) so geometry matches ISV-web
exactly, then applies ISV-web's frontend ``static/app.js`` ``mergeSegments()`` step (coalesce
adjacent/overlapping exons by block label -> dedupe duplicate exon tokens and merge a block's sub-exons;
keep retained introns separate) so the drawn structure matches what ISV-web renders -- in particular for
novel collapsed isoforms whose raw tokens carry terminal coords (e.g. ``E34.2_1887588|I34.1_1891712|
E34.2|I34.1``). CDS region from ``gff-output/coding_regions.txt`` is drawn TALL, UTR SHORT, with
strand-aware translation begin/end markers; NMD isoforms hatched.

Fast + faithful: the gene-model exon lookup comes from ISV's ``gene_model.db`` index (per-gene, identical
coords); isoform token structures come from the TSV ``transcript_associations.txt`` (the ISV
``transcripts.db`` index drops novel per-sample isoforms, so it is deliberately NOT used here);
``intron_scale`` defaults to 0.2 (ISV's ``build_gene_maps`` default) and is overridable.

Usage:
    R = StructRenderer(gene_model="…/Hs_Ensembl_exon.txt", gff="…/gff-output", intron_scale=0.2)
    R.load_genes(["ENSG…"])
    xlim = R.shared_xlim([("ENSG…", iso) for iso in isoforms])   # align stacked tracks
    R.draw(ax, "ENSG…", iso, color="#377EB8", xlim=xlim)         # one isoform per axes
"""
import os
import sys
import subprocess
import tempfile

import numpy as np
import matplotlib.patches as mpatches

# Import the validated ISV engine the same proven way isoform_track.py does: put components/ (so the
# engine's ``from long_read import …`` resolves) and the visualization dir on sys.path, then import the
# module by name. This avoids re-implementing any of the engine.
_HERE = os.path.dirname(os.path.abspath(__file__))                     # …/iso2function/plotting
_COMPONENTS = os.path.dirname(os.path.dirname(_HERE))                  # …/altanalyze3/components
_ENGINE_DIR = os.path.join(_COMPONENTS, "visualization")
if _COMPONENTS not in sys.path:
    sys.path.insert(0, _COMPONENTS)
if _ENGINE_DIR not in sys.path:
    sys.path.insert(0, _ENGINE_DIR)
import isoform_structure_view as isv  # noqa: E402


class StructRenderer:
    """Render isoform structure tracks for a dataset's gene model + gff-output.

    Parameters
    ----------
    gene_model : path to the Ensembl exon coords TSV (``Hs_Ensembl_exon.txt``). ISV builds/uses a
        ``gene_model.db`` next to it for fast per-gene exon lookup.
    gff : path to the dataset's ``gff-output`` dir (holds ``transcript_associations.txt`` +
        ``coding_regions.txt``).
    intron_scale : float, default 0.2 (ISV's ``build_gene_maps`` default). ``INTRON_SCALE`` env overrides
        when ``intron_scale`` is left None.
    """

    def __init__(self, gene_model, gff, intron_scale=None):
        self.intron_scale = float(os.environ.get("INTRON_SCALE", 0.2 if intron_scale is None else intron_scale))
        self.gff = gff
        self.gene_model = gene_model
        self._gmdb = None                                      # ISV gene_model.db (fast per-gene exon lookup)
        try:
            self._gmdb = isv._ensure_gene_model_db(gene_model)
        except Exception as e:                                 # pragma: no cover - index build failure
            print(f"[struct] gene_model.db unavailable ({e}); using full text model")
        if self._gmdb:
            self.gs, self.el = {}, {}                          # loaded lazily per gene from sqlite
        else:
            self.gs, self.el = isv.read_gene_model(gene_model)
        self.gmaps, self.tokens, self.cds, self._loaded = {}, {}, {}, set()

    def load_genes(self, ensgs):
        """Load exon lookup (gene_model.db), isoform token structures + coding regions for these genes."""
        ensgs = [e for e in set(ensgs) if e not in self._loaded]
        if not ensgs:
            return
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write("\n".join(ensgs) + "\n"); pat = fh.name
        if self._gmdb:                                         # per-gene exon lookup from the ISV index (fast)
            for g in ensgs:
                try:
                    gseg, el = isv.read_gene_model_from_sqlite(self._gmdb, g)
                    self.gs.update(gseg); self.el.update(el)
                except Exception as e:
                    print(f"[struct] gene-model load warn {g}: {e}")
        # isoform tokens from the TSV -- the ISV transcripts.db index DROPS novel per-sample isoforms, so
        # it cannot be used; the TSV has both known + novel.
        try:
            out = subprocess.run(["grep", "-F", "-f", pat, f"{self.gff}/transcript_associations.txt"],
                                 capture_output=True, text=True, timeout=1800).stdout
            for ln in out.splitlines():
                p = ln.split("\t")
                if len(p) >= 4 and p[0] in ensgs:
                    self.tokens[p[3]] = p[2].split("|")
        except Exception as e:
            print(f"[struct] token load warn: {e}")
        # coding_regions.txt: tid gene chrom strand cds_start cds_end gmin gmax prot_len tx_len nmd
        try:
            out = subprocess.run(["grep", "-F", "-f", pat, f"{self.gff}/coding_regions.txt"],
                                 capture_output=True, text=True, timeout=900).stdout
            for ln in out.splitlines():
                p = ln.split("\t")
                if len(p) >= 11 and p[1] in ensgs:
                    try:
                        self.cds[p[0]] = dict(strand=p[3], cds_start=int(p[4]), cds_end=int(p[5]),
                                              gmin=int(p[6]), gmax=int(p[7]),
                                              prot=int(float(p[8])) if p[8] else None, nmd=p[10].strip())
                    except Exception:
                        pass
        except Exception as e:
            print(f"[struct] cds load warn: {e}")
        os.remove(pat); self._loaded |= set(ensgs)

    def gene_map(self, gene):
        if gene not in self.gmaps and gene in self.gs:
            self.gmaps[gene] = isv.build_gene_maps({gene: self.gs[gene]}, self.intron_scale)[gene]
        return self.gmaps.get(gene)

    def gene_strand(self, gene):
        if not hasattr(self, "_gstrand"):
            self._gstrand = {}
        if gene not in self._gstrand:
            s = "+"
            for k, v in self.el.items():
                if k[0] == gene and v.get("strand") in ("+", "-"):
                    s = v["strand"]; break
            self._gstrand[gene] = s
        return self._gstrand[gene]

    def _segs(self, gene, iso):
        tok = self.tokens.get(iso)
        if tok is None:
            return None
        tok = isv.strip_terminal_coords(tok, gene)             # ISV-faithful token normalization
        return isv.build_isoform_segments(tok, self.el, gene)

    def _cds(self, iso):
        return self.cds.get(iso) or self.cds.get(iso.split(".", 1)[0])   # ENST direct, novel by numeric id

    @staticmethod
    def _merge_segments(segs):
        """Port of ISV-web static/app.js mergeSegments(): coalesce adjacent/overlapping exons by block
        label (dedupe duplicate exon tokens, merge a block's sub-exons), keep retained introns separate.
        Returns [(gstart, gend, 'E'|'I'), ...] in genomic coordinates."""
        _E = sorted((min(s["start"], s["end"]), max(s["start"], s["end"]),
                     str(s.get("label") or "").split(".")[0])
                    for s in segs if str(s.get("type", "E")).upper().startswith("E"))
        merged = []
        for a, b, lab in _E:
            if merged and ((a <= merged[-1][1] + 2 and lab == merged[-1][2]) or a <= merged[-1][1]):
                merged[-1] = [merged[-1][0], max(merged[-1][1], b), merged[-1][2] or lab]
            else:
                merged.append([a, b, lab])
        introns = [(min(s["start"], s["end"]), max(s["start"], s["end"]))
                   for s in segs if str(s.get("type", "")).upper().startswith("I")]
        return [(a, b, "E") for a, b, _ in merged] + [(a, b, "I") for a, b in introns]

    def extent(self, gene, iso):
        gm = self.gene_map(gene); segs = self._segs(gene, iso)
        if not gm or not segs:
            return None
        xs = [isv.map_coord(gm, c) for a, b, _ in self._merge_segments(segs) for c in (a, b)]
        return (min(xs), max(xs))

    def shared_xlim(self, pairs, pad=0.04):
        ext = [self.extent(g, i) for g, i in pairs]; ext = [e for e in ext if e]
        if not ext:
            return None
        lo = min(e[0] for e in ext); hi = max(e[1] for e in ext); d = (hi - lo) * pad or 1
        return (lo - d, hi + d)

    def draw(self, ax, gene, iso, color="#444444", xlim=None, label=None, lw_back=0.9, label_fontsize=7):
        """Draw one isoform structure track filling ``ax`` (ax ylim fixed to [0,1])."""
        gm = self.gene_map(gene); segs = self._segs(gene, iso)
        ax.set_ylim(0, 1); ax.set_yticks([])
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.set_xticks([]); ax.patch.set_visible(False)
        lab = label if label is not None else iso
        if not gm or not segs:
            ax.text(0.0, 0.5, f"{lab}  (structure n/a)", fontsize=label_fontsize, va="center", ha="left",
                    transform=ax.transAxes, color=color)
            return
        D = [(min(isv.map_coord(gm, a), isv.map_coord(gm, b)),
              max(isv.map_coord(gm, a), isv.map_coord(gm, b)), t) for a, b, t in self._merge_segments(segs)]
        tx0 = min(d[0] for d in D); tx1 = max(d[1] for d in D)
        if xlim is None:
            pad = (tx1 - tx0) * 0.02 or 1; xlim = (tx0 - pad, tx1 + pad)
        ax.set_xlim(*xlim)
        c = self._cds(iso); cmin = cmax = None; nmd = False
        if c:
            cmin = min(isv.map_coord(gm, c["gmin"]), isv.map_coord(gm, c["gmax"]))
            cmax = max(isv.map_coord(gm, c["gmin"]), isv.map_coord(gm, c["gmax"]))
            nmd = ("NMD" in c["nmd"] and "Not" not in c["nmd"])
        CY = 0.40; UT_H, CD_H, IR_H = 0.24, 0.50, 0.16
        strand = c["strand"] if c else self.gene_strand(gene)
        ax.plot([tx0, tx1], [CY, CY], color="#999999", lw=lw_back, zorder=1, solid_capstyle="butt", clip_on=False)
        coding = []
        for d0, d1, typ in D:
            if cmin is None:
                h = IR_H if typ == "I" else UT_H; al = 0.35 if typ == "I" else 1.0
                ax.add_patch(mpatches.Rectangle((d0, CY - h / 2), d1 - d0, h, facecolor=color, edgecolor="none",
                                                alpha=al, zorder=(2 if typ == "I" else 3), clip_on=False))
                continue
            for a, b, iscd in [(d0, min(d1, cmin), False), (max(d0, cmin), min(d1, cmax), True), (max(d0, cmax), d1, False)]:
                if b <= a:
                    continue
                if iscd:
                    coding.append((a, b))
                if typ == "I":
                    if iscd:
                        pat = mpatches.Rectangle((a, CY - CD_H / 2), b - a, CD_H, facecolor=color, edgecolor=color,
                                                 lw=0.4, alpha=0.55, zorder=3, clip_on=False); pat.set_hatch("....")
                    else:
                        pat = mpatches.Rectangle((a, CY - IR_H / 2), b - a, IR_H, facecolor=color, edgecolor="none",
                                                 alpha=0.35, zorder=2, clip_on=False)
                else:
                    h = CD_H if iscd else UT_H
                    pat = mpatches.Rectangle((a, CY - h / 2), b - a, h, facecolor=color, edgecolor="none",
                                             zorder=(4 if iscd else 3), clip_on=False)
                    if iscd and nmd:
                        pat.set_hatch("////"); pat.set_edgecolor("#ffffff")
                ax.add_patch(pat)
        if coding:
            lo = min(x[0] for x in coding); hi = max(x[1] for x in coding)
            beg, end = ((lo, hi) if strand == "+" else (hi, lo)); my = CY + CD_H / 2 + 0.18
            for mx, mc in [(beg, "#1a9850"), (end, "#d73027")]:
                ax.plot([mx, mx], [CY - CD_H / 2, my], color=mc, lw=0.8, zorder=6, clip_on=False)
                ax.plot([mx], [my], marker="v", ms=5, color=mc, zorder=7, clip_on=False)
        arr = "▶" if strand == "+" else "◀"
        aa = f"{c['prot']}aa" if (c and c.get("prot") is not None) else "no CDS"
        ax.text(0.0, 0.97, f"{lab} {arr} {aa}{'  (NMD)' if nmd else ''}", fontsize=label_fontsize, va="top",
                ha="left", transform=ax.transAxes, color=color, zorder=8)
