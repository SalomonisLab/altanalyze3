#!/usr/bin/env python3
"""Validate a run_workflow.py output directory: required files present + non-empty, PDFs render
(no image-decode failures) with embedded fonts, enrichment tables non-trivial, logs/timeline complete.
Usage: python validate_workflow_outputs.py <workflow_outdir>"""
import os, sys, pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
problems = []


def chk(cond, msg):
    print(("  [OK ] " if cond else "  [FAIL] ") + msg)
    if not cond:
        problems.append(msg)
    return cond


def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def pdf_ok(path):
    """returns (valid, image_decode_failures, embedded_font_count)"""
    try:
        import pypdf
        r = pypdf.PdfReader(path)
        bad = 0
        for pg in r.pages:
            for im in pg.images:
                try:
                    _ = im.data
                except Exception:
                    bad += 1
        d = open(path, "rb").read()
        emb = d.count(b"/Type0") + d.count(b"/TrueType")
        return True, bad, emb
    except Exception as e:
        return False, -1, -1


print(f"=== validating workflow outputs: {OUT} ===")

print("\n[logs + timeline]")
tl_path = os.path.join(OUT, "logs", "timeline.tsv")
chk(nonempty(os.path.join(OUT, "logs", "workflow.log")), "logs/workflow.log present")
chk(nonempty(tl_path), "logs/timeline.tsv present")
chk(nonempty(os.path.join(OUT, "WORKFLOW_SUMMARY.txt")), "WORKFLOW_SUMMARY.txt present")
chk(nonempty(os.path.join(OUT, "command.txt")), "command.txt present (reproducible)")
if nonempty(tl_path):
    tl = pd.read_csv(tl_path, sep="\t")
    chk((tl["status"] == "OK").all(), f"all {len(tl)} steps status=OK ({list(tl['status'])})")
    chk("seconds" in tl.columns and tl["seconds"].notna().all(), "every step has a timing")
    print("   timeline:\n" + tl.to_string(index=False).replace("\n", "\n   "))

print("\n[UDON core outputs]")
udon = os.path.join(OUT, "udon")
chk(nonempty(os.path.join(udon, "udon_core", "udon_clusters.txt")), "udon_clusters.txt non-empty")
chk(nonempty(os.path.join(udon, "goelite", "GOElite_UDON.tsv")), "goelite/GOElite_UDON.tsv non-empty")
for f in ("marker_heatmap.pdf", "donor_cluster_heatmap.pdf", "celltype_cluster_heatmap.pdf"):
    p = os.path.join(udon, f)
    if chk(nonempty(p), f"{f} non-empty"):
        v, bad, emb = pdf_ok(p)
        chk(v and bad == 0, f"{f} renders (image-decode-failures={bad}, Adobe-safe)")
        chk(emb > 0, f"{f} has embedded fonts (Type0/TrueType count={emb})")

print("\n[SATAY Fisher outputs]")
sf = os.path.join(OUT, "satay_fisher")
ep = os.path.join(sf, "satay_cluster_enrichment.tsv")
if chk(nonempty(ep), "satay_cluster_enrichment.tsv non-empty"):
    e = pd.read_csv(ep, sep="\t")
    chk(len(e) > 0 and (e["fdr"] < 0.05).any(), f"has significant enrichments ({int((e['fdr']<0.05).sum())} FDR<0.05)")
for f in ("satay_byfield.pdf", "donor_cluster_heatmap.pdf", "celltype_cluster_heatmap.pdf"):
    p = os.path.join(sf, f)
    if chk(nonempty(p), f"{f} non-empty"):
        v, bad, emb = pdf_ok(p)
        chk(v and bad == 0, f"{f} renders (image-decode-failures={bad})")

print("\n[SATAY CMH (batch-stratified) outputs]")
cp = os.path.join(OUT, "satay_cmh", "satay_cluster_enrichment.tsv")
if chk(nonempty(cp), "satay_cmh/satay_cluster_enrichment.tsv non-empty"):
    c = pd.read_csv(cp, sep="\t")
    chk(len(c) > 0, f"CMH produced {len(c)} testable, {int((c['fdr']<0.05).sum())} significant (FDR<0.05)")

print("\n=== " + ("ALL CHECKS PASSED" if not problems else f"{len(problems)} PROBLEM(S)") + " ===")
sys.exit(0 if not problems else 1)
