#!/usr/bin/env python3
"""Workflow summary tables + bar charts for the long-read single-cell isoform pipeline.

Writes FOUR TSV/chart pairs into a ``summary/`` folder beside the run outputs. Every number is
read from the live collapse results (passed in) or from the artifacts already on disk (per-sample
junction h5ads, transcript_associations, FINAL_isoform_catalog, protein_summary) -- nothing is
guessed or re-derived by re-running analysis.

  1. per_sample_counts.tsv      reads (molecules), junctions detected, unique junctions per sample
  2. overall_counts.tsv         total reads, total unique structures, collapsed (kept) isoforms,
                                collapse_method (wta/em), removed-low-detection, + per-sample
                                final-isoform / read-conservation breakdown
  3. known_vs_novel.tsv         Ensembl(ENST)-known vs novel isoform counts
  4. protein_class.tsv          NMD vs truncated(<25% of the gene's longest isoform) vs protein-coding

Charts are vector PDFs with editable text (Type-42 fonts, Arial) so they drop straight into
Illustrator for figures.
"""
import os
import csv

# Charts are OPTIONAL: the TSV tables must always be written even where matplotlib is unavailable
# (e.g. a minimal cluster python env). Import lazily; if it fails, skip charts and keep the tables.
_HAVE_MPL = True
try:
    import matplotlib
    matplotlib.use("Agg")
    # Manuscript-ready vector text: keep fonts as editable text (Type 42) in Arial.
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "Arial"
    import matplotlib.pyplot as plt
except Exception:
    _HAVE_MPL = False


def _write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _bar(path, labels, series, title, ylabel, legend=None, rotate=90):
    """Vector bar chart. series: list of (name, values) for grouped bars, or a single value list.
    No-op if matplotlib is unavailable (the TSV table is always written regardless)."""
    if not _HAVE_MPL:
        return
    import numpy as np
    fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(labels) + 2), 4.2))
    x = np.arange(len(labels))
    if isinstance(series[0], tuple):
        n = len(series)
        w = 0.8 / n
        for i, (nm, vals) in enumerate(series):
            ax.bar(x + (i - (n - 1) / 2) * w, vals, w, label=nm)
        if legend is not False:
            ax.legend(frameon=False, fontsize=8)
    else:
        ax.bar(x, series, 0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rotate, ha="center" if rotate == 90 else "right", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _read_junction_count(h5ad_path):
    """Unique junctions for a sample = var count of its <lib>-junction.h5ad (no full load)."""
    try:
        import h5py
        with h5py.File(h5ad_path, "r", locking=False) as f:
            for p in ("var/_index", "var/index"):
                if p in f:
                    return int(f[p].shape[0])
            g = f["var"]
            return int(g[g.attrs.get("_index", "_index")].shape[0])
    except Exception:
        return None


def _count_reads_from_ta(ta_path):
    """Reads (molecules) for a sample = non-UNK rows of its transcript_associations.txt."""
    n = 0
    try:
        with open(ta_path) as f:
            for line in f:
                p = line.split("\t", 1)
                if p and "UNK" not in p[0]:
                    n += 1
    except Exception:
        return None
    return n


def write_summary(outdir, *, collapse_method, min_total, n_structures, catalog, h5ad_results,
                  sample_junction_h5ads=None, sample_ta_paths=None, protein_summary_path=None,
                  truncated_frac=0.25, log=print):
    """Write the four summary TSVs + bar charts into ``<outdir>/summary/``.

    Parameters mirror what run_pipeline holds at the end:
      n_structures        : total unique isoform structures pre-filter (n_final).
      catalog             : kept collapsed isoforms (list of rows; row[5] truthy => known/ENST).
      h5ad_results        : list of (name, out_path, n_final_iso, raw_reads, final_reads, rss, secs).
      sample_junction_h5ads / sample_ta_paths : {sample -> path} for per-sample junctions / reads.
      protein_summary_path: gff-output/protein_summary.txt for the protein-class breakdown.
    """
    sdir = os.path.join(outdir, "summary")
    os.makedirs(sdir, exist_ok=True)
    sample_junction_h5ads = sample_junction_h5ads or {}
    sample_ta_paths = sample_ta_paths or {}

    # ---- 1. per-sample: reads (molecules), junctions detected, unique junctions ----
    per_rows = []
    samples = [r[0] for r in h5ad_results] if h5ad_results else sorted(
        set(sample_ta_paths) | set(sample_junction_h5ads))
    for name in samples:
        reads = _count_reads_from_ta(sample_ta_paths.get(name, "")) if name in sample_ta_paths else None
        nj = _read_junction_count(sample_junction_h5ads.get(name, "")) if name in sample_junction_h5ads else None
        per_rows.append([name, reads if reads is not None else "NA",
                         nj if nj is not None else "NA"])
    _write_tsv(os.path.join(sdir, "per_sample_counts.tsv"),
               ["sample", "reads_molecules", "unique_junctions"], per_rows)
    labels = [r[0] for r in per_rows]
    reads_v = [r[1] if isinstance(r[1], int) else 0 for r in per_rows]
    junc_v = [r[2] if isinstance(r[2], int) else 0 for r in per_rows]
    _bar(os.path.join(sdir, "per_sample_counts.pdf"), labels,
         [("reads (molecules)", reads_v), ("unique junctions", junc_v)],
         "Per-sample reads and junctions", "count")
    log(f"[summary] wrote per_sample_counts.tsv ({len(per_rows)} samples)")

    # ---- 2. overall counts (+ per-sample final-isoform / read conservation) ----
    n_kept = len(catalog)
    n_removed = n_structures - n_kept
    tot_raw = sum(r[3] for r in h5ad_results) if h5ad_results else 0
    tot_final = sum(r[4] for r in h5ad_results) if h5ad_results else 0
    overall = [
        ["collapse_method", collapse_method],
        ["min_total_filter", min_total],
        ["total_reads_all_samples", tot_raw],
        ["total_reads_retained_final", tot_final],
        ["total_unique_isoform_structures", n_structures],
        ["collapsed_isoforms_kept", n_kept],
        ["isoforms_removed_low_detection", n_removed],
    ]
    _write_tsv(os.path.join(sdir, "overall_counts.tsv"), ["metric", "value"], overall)
    # overall bar (the count metrics, not method/min_total)
    obar = [("total reads", tot_raw), ("reads retained", tot_final),
            ("unique structures", n_structures), ("collapsed kept", n_kept),
            ("removed (<min_total)", n_removed)]
    _bar(os.path.join(sdir, "overall_counts.pdf"), [k for k, _ in obar],
         [v for _, v in obar], f"Overall isoform counts (method={collapse_method}, min_total={min_total})",
         "count", rotate=30)
    # per-sample breakdown sidecar TSV (final isoforms + conservation), referenced by overall
    if h5ad_results:
        psb = [[r[0], r[2], r[3], r[4],
                f"{(100.0 * r[4] / r[3]):.1f}" if r[3] else "NA"] for r in h5ad_results]
        _write_tsv(os.path.join(sdir, "overall_per_sample_breakdown.tsv"),
                   ["sample", "final_isoforms", "raw_reads", "final_reads", "pct_reads_retained"], psb)
    log(f"[summary] wrote overall_counts.tsv (method={collapse_method}, kept={n_kept:,}, "
        f"removed={n_removed:,})")

    # ---- 3. known (Ensembl/ENST) vs novel ----
    n_known = sum(1 for row in catalog if len(row) > 5 and row[5])
    n_novel = n_kept - n_known
    _write_tsv(os.path.join(sdir, "known_vs_novel.tsv"),
               ["class", "isoforms"], [["known_ENST", n_known], ["novel", n_novel]])
    _bar(os.path.join(sdir, "known_vs_novel.pdf"), ["known (ENST)", "novel"],
         [n_known, n_novel], "Known (Ensembl) vs novel isoforms", "isoforms", rotate=0)
    log(f"[summary] wrote known_vs_novel.tsv (known={n_known:,} novel={n_novel:,})")

    # ---- 4. protein class: NMD vs truncated(<25% longest) vs protein-coding ----
    if protein_summary_path and os.path.exists(protein_summary_path):
        n_nmd = n_trunc = n_coding = n_unknown = 0
        with open(protein_summary_path) as f:
            next(f, None)  # header
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                nmd = p[3].strip()
                try:
                    length = int(p[2]); longest = int(p[5])
                except ValueError:
                    n_unknown += 1
                    continue
                # NMD status values are exactly 'Not-NMD' / 'Potential-NMD' / 'NMD'. Must NOT use
                # `"NMD" in nmd` -- that matches 'Not-NMD' too. NMD class = anything but 'Not-NMD'.
                is_nmd = (nmd != "Not-NMD") and ("NMD" in nmd)
                if is_nmd:
                    n_nmd += 1
                elif longest > 0 and length < truncated_frac * longest:
                    n_trunc += 1
                else:
                    n_coding += 1
        rows = [["NMD", n_nmd],
                [f"truncated_lt_{int(truncated_frac*100)}pct_longest", n_trunc],
                ["protein_coding", n_coding]]
        if n_unknown:
            rows.append(["unclassified", n_unknown])
        _write_tsv(os.path.join(sdir, "protein_class.tsv"), ["class", "isoforms"], rows)
        _bar(os.path.join(sdir, "protein_class.pdf"), [r[0] for r in rows], [r[1] for r in rows],
             "Isoform protein class", "isoforms", rotate=20)
        log(f"[summary] wrote protein_class.tsv (NMD={n_nmd:,} truncated={n_trunc:,} "
            f"coding={n_coding:,})")
    else:
        log("[summary] protein_summary.txt not found; skipping protein_class table")

    return sdir
