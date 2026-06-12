#!/usr/bin/env python3
"""Study-aware UDON — the COMPLETE, HARDWIRED, logged + timed workflow.

ONE command runs the entire batch-corrected analysis into UDON/study_aware[<_tag>]/ with a logs/
folder (per-step logs + timeline.tsv + timing) and a WORKFLOW_SUMMARY.txt. NO manual steps, NO
hand-run scripts, NO AI-generated artifacts.

Pipeline:
  1 study_aware_udon      per-study UDON (matched sex+platform controls, protein-coding gene filter)
  2 study_aware_genelist  conserved vs unique programs (top-50 marker gene-list overlap)
  3 study_aware_integrate SVA integration of conserved centroids -> final_program_assignments.tsv
  4 final_heatmap         MarkerFinder heatmap of the final programs WITH hardwired gene+GO callouts
                          (computes per-program top gene + runs GO-Elite itself) + per-program GO-Elite
                          (goelite/) + donor/cell-type composition heatmaps
  5 satay_fisher          SATAY-UDON donor-level Fisher of the final programs + byfield/composition heatmaps
  6 satay_cmh             SATAY-UDON batch-stratified CMH (controls for --batch-key)
The satay_*/satay_cluster_enrichment.tsv tables ARE the per-program enrichment summary.

Reproducible (the exact command is saved to <SA>/command.txt):
  python run_study_aware.py --metadata satay_metadata.tsv [--mean-binarize age,blast_pct]
      [--batch-key Study] [--tag NAME] [--n-run 5] [--no-gene-filter]
"""
import os, sys, time, argparse, datetime, subprocess
import pandas as pd

UDON_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, UDON_DIR)
import pseudobulk_protocol as P

PB = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk"


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def run_step(name, cmd, logs, env, timeline, master):
    log_path = os.path.join(logs, f"{name}.log")
    t0 = time.time(); start = _now()
    print(f"[{start}] >>> {name}", flush=True); master.write(f"[{start}] >>> {name}\n"); master.flush()
    with open(log_path, "w") as fh:
        fh.write(f"# step: {name}\n# start: {start}\n# cmd: {' '.join(cmd)}\n\n"); fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=UDON_DIR).returncode
        dur = time.time() - t0
        fh.write(f"\n# end: {_now()}  ({dur:.1f}s)  exit={rc}\n")
    status = "OK" if rc == 0 else f"FAIL(exit={rc})"
    timeline.append({"step": name, "start": start, "seconds": round(dur, 1), "status": status,
                     "log": os.path.basename(log_path)})
    print(f"[{_now()}] <<< {name}: {status} in {dur:.1f}s  (log: logs/{name}.log)", flush=True)
    master.write(f"[{_now()}] <<< {name}: {status} in {dur:.1f}s\n"); master.flush()
    return rc == 0


def main():
    ap = argparse.ArgumentParser(description="Study-aware UDON master workflow (hardwired, logged + timed)")
    ap.add_argument("--metadata", required=True, help="SATAY covariate table (SATAY_METADATA_FORMAT.md)")
    ap.add_argument("--mean-binarize", default=None, help="comma list of numeric covariates to binarize at the mean")
    ap.add_argument("--batch-key", default="Study", help="metadata column for the batch-stratified CMH step")
    ap.add_argument("--tag", default="", help="output-folder suffix -> UDON/study_aware_<tag>/")
    ap.add_argument("--n-run", type=int, default=5, help="NMF restarts (5=stability default; 1=fast)")
    ap.add_argument("--no-gene-filter", action="store_true")
    args = ap.parse_args()

    sfx = ("_" + args.tag.replace("/", "-")) if args.tag else ""
    SA = os.path.join(PB, "UDON", "study_aware" + sfx)
    os.makedirs(SA, exist_ok=True)
    logs = os.path.join(SA, "logs"); os.makedirs(logs, exist_ok=True)
    with open(os.path.join(SA, "command.txt"), "w") as fh:
        fh.write("python " + " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]) + "\n")

    env = os.environ.copy()
    if args.tag:
        env["UDON_TAG"] = args.tag
    env["UDON_NMF_RUNS"] = str(args.n_run)
    env["UDON_GENE_FILTER"] = "0" if args.no_gene_filter else "1"
    env["UDON_METADATA"] = os.path.abspath(args.metadata)   # so final_heatmap can render the covariate heatmap
    if args.mean_binarize:
        env["UDON_MEAN_BINARIZE"] = args.mean_binarize      # same numeric->binary list the covariate heatmap needs
    PY = sys.executable
    timeline = []; t_all = time.time()
    master = open(os.path.join(logs, "workflow.log"), "w")
    master.write(f"# study-aware UDON workflow\n# output: {SA}\n# started: {_now()}\n\n")

    # ---- steps 1-5: the study-aware integration + final heatmap (callouts hardwired in final_heatmap) ----
    # final_heatmap is fully self-contained (computes its own gene + GO callouts + per-program GO-Elite),
    # so the whole critical path is: per-study UDON -> conserved/unique -> SVA integrate -> final heatmap.
    # (The SATAY enrichment tables produced below ARE the per-program enrichment summary; the legacy
    # annotate_summary step is retired -- its callout/GO-Elite role is now hardwired into final_heatmap.)
    seq = [("01_study_aware_udon", "study_aware_udon.py", True),
           ("02_genelist", "study_aware_genelist.py", True),
           ("03_integrate", "study_aware_integrate.py", True),
           ("04_final_heatmap", "final_heatmap.py", True)]
    for name, script, critical in seq:
        ok = run_step(name, [PY, os.path.join(UDON_DIR, script)], logs, env, timeline, master)
        if not ok and critical:
            if name == "03_integrate" and not os.path.exists(os.path.join(SA, "final_program_assignments.tsv")):
                master.write("!! integration produced no final programs (see integration_report.txt); stopping.\n")
            else:
                master.write(f"!! critical step {name} failed; stopping.\n")
            break

    fa = os.path.join(SA, "final_program_assignments.tsv")
    # ---- steps 6-7: SATAY-UDON via the canonical engine on the final programs ----
    if os.path.exists(fa):
        clusters = os.path.join(SA, "final_program_clusters.txt")
        _f = pd.read_csv(fa, sep="\t")
        _f[["pseudobulk", "final_program"]].set_index("pseudobulk").to_csv(clusters, sep="\t")
        common = [PY, os.path.join(UDON_DIR, "satay_udon_core.py"), "--metadata", args.metadata,
                  "--udon-clusters", clusters]
        if args.mean_binarize:
            common += ["--mean-binarize", args.mean_binarize]
        run_step("06_satay_fisher", common + ["--outdir", os.path.join(SA, "satay_fisher"), "--heatmaps"],
                 logs, env, timeline, master)
        run_step("07_satay_cmh", common + ["--outdir", os.path.join(SA, "satay_cmh"), "--batch-key", args.batch_key],
                 logs, env, timeline, master)

    # ---- timeline + summary + built-in validation ----
    tl = pd.DataFrame(timeline)
    tl.to_csv(os.path.join(logs, "timeline.tsv"), sep="\t", index=False)
    total = time.time() - t_all

    # built-in checks (no external script): final heatmap callouts + SATAY enrichments
    checks = []
    hp = os.path.join(SA, "final_program_heatmap.pdf")
    if os.path.exists(hp):
        try:
            import pypdf
            txt = pypdf.PdfReader(hp).pages[0].extract_text()
            mk = pd.read_csv(os.path.join(SA, "final_program_markers.txt"), sep="\t")
            col = "pearson_r" if "pearson_r" in mk.columns else mk.columns[2]
            genes = {str(c): str(g.sort_values(col, ascending=False).iloc[0]["marker"])
                     for c, g in mk.groupby("top_cluster")}
            ng = sum(1 for v in genes.values() if v in txt)
            sel = os.path.join(SA, "goelite", "GOElite_UDON_selected.tsv")
            ngo = 0
            if os.path.exists(sel):
                s = pd.read_csv(sel, sep="\t")
                gos = {str(c): str(g.sort_values("fdr").iloc[0]["term_name"]) for c, g in s.groupby("cluster")}
                ngo = sum(any(w in txt for w in str(v).split()) for v in gos.values())
            checks.append(f"final_program_heatmap callouts: {ng}/{len(genes)} gene + {ngo} GO-term (hardwired)")
        except Exception as e:
            checks.append(f"heatmap callout check skipped: {e}")
    for lab, sub in [("satay_fisher", "satay_fisher"), ("satay_cmh", "satay_cmh")]:
        ep = os.path.join(SA, sub, "satay_cluster_enrichment.tsv")
        if os.path.exists(ep):
            e = pd.read_csv(ep, sep="\t")
            checks.append(f"{lab}: {len(e)} testable, {int((e['fdr'] < 0.05).sum())} significant (FDR<0.05)")

    with open(os.path.join(SA, "WORKFLOW_SUMMARY.txt"), "w") as fh:
        fh.write(f"study-aware UDON workflow\noutput: {SA}\nfinished: {_now()}\n")
        fh.write(f"total wall-clock: {total:.1f}s ({total/60:.1f} min)\n\n=== step timeline ===\n")
        fh.write(tl.to_string(index=False) + "\n\n=== built-in validation ===\n")
        for c in checks:
            fh.write("  " + c + "\n")
    master.write(f"\n# finished: {_now()}  total {total:.1f}s\n"); master.close()
    print(f"\n[study-aware] DONE in {total/60:.1f} min -> {SA}")
    for c in checks:
        print("  " + c)
    print(f"  summary: {SA}/WORKFLOW_SUMMARY.txt | timeline: {logs}/timeline.tsv")


if __name__ == "__main__":
    main()
