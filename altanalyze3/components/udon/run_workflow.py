#!/usr/bin/env python3
"""Master UDON + SATAY-UDON workflow with built-in per-step LOGGING + TIMING.

Everything runs into a single --outdir; per-step logs and a timeline go to <outdir>/logs/.
Steps:
  01_udon_core      run_udon_pseudobulk (matched sex+platform controls, protein-coding gene filter,
                    NMF n_run=5) -> clusters, markers, marker heatmap, GO-Elite, donor/celltype heatmaps
  02_satay_fisher   satay_udon_core (donor-level Fisher) on the clusters + standard --metadata
                    -> enrichment table + satay_byfield + composition heatmaps
  03_satay_cmh      satay_udon_core --batch-key <col> (Cochran-Mantel-Haenszel, batch-stratified)
                    -> batch-controlled enrichment table

Outputs:  <outdir>/udon/, <outdir>/satay_fisher/, <outdir>/satay_cmh/
Logs:     <outdir>/logs/<step>.log , <outdir>/logs/timeline.tsv , <outdir>/WORKFLOW_SUMMARY.txt

Reproducible (the exact command is saved to <outdir>/command.txt):
  python run_workflow.py --outdir DIR --pseudobulk H5AD --counts H5AD \
      --sample-metadata SAMPLE.tsv --metadata COVARIATES.tsv [--cell-type X] [--mean-binarize a,b] [--batch-key Study]
"""
import os, sys, time, json, argparse, datetime, subprocess

UDON_DIR = os.path.dirname(os.path.abspath(__file__))


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def run_step(name, cmd, logs_dir, env, timeline, master):
    log_path = os.path.join(logs_dir, f"{name}.log")
    t0 = time.time(); start = _now()
    line = f"[{start}] >>> {name}"
    print(line, flush=True); master.write(line + "\n"); master.flush()
    with open(log_path, "w") as fh:
        fh.write(f"# step:  {name}\n# start: {start}\n# cmd:   {' '.join(cmd)}\n\n"); fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=UDON_DIR).returncode
        dur = time.time() - t0
        fh.write(f"\n# end: {_now()}  ({dur:.1f}s)  exit={rc}\n")
    status = "OK" if rc == 0 else f"FAIL(exit={rc})"
    timeline.append({"step": name, "start": start, "seconds": round(dur, 1), "status": status,
                     "log": os.path.basename(log_path)})
    line = f"[{_now()}] <<< {name}: {status} in {dur:.1f}s  (log: logs/{name}.log)"
    print(line, flush=True); master.write(line + "\n"); master.flush()
    return rc == 0


def main():
    ap = argparse.ArgumentParser(description="UDON + SATAY-UDON master workflow (logged + timed)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--counts", required=True, help="raw-count pseudobulk h5ad (sex calling)")
    ap.add_argument("--sample-metadata", required=True, help="UDON sample metadata (Sample, chemistry, Dataset)")
    ap.add_argument("--metadata", required=True, help="SATAY covariate table (SATAY_METADATA_FORMAT.md)")
    ap.add_argument("--cell-type", default=None, help="restrict the whole analysis to one cell type")
    ap.add_argument("--mean-binarize", default=None, help="comma list of numeric covariates to binarize at the mean")
    ap.add_argument("--batch-key", default="Study", help="metadata column for the batch-stratified CMH step")
    ap.add_argument("--species", default="Hs")
    ap.add_argument("--no-gene-filter", action="store_true")
    ap.add_argument("--fast", action="store_true", help="vectorized feature selection (NMF n_run unchanged)")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    logs = os.path.join(outdir, "logs"); os.makedirs(logs, exist_ok=True)
    with open(os.path.join(outdir, "command.txt"), "w") as fh:
        fh.write("python " + " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]) + "\n")
    env = os.environ.copy()
    env["UDON_METADATA"] = os.path.abspath(args.metadata)   # so the UDON step renders the covariate heatmap
    if args.mean_binarize:
        env["UDON_MEAN_BINARIZE"] = args.mean_binarize      # same numeric->binary list the covariate heatmap needs
    PY = sys.executable
    timeline = []
    t_all = time.time()
    master = open(os.path.join(logs, "workflow.log"), "w")
    master.write(f"# UDON + SATAY-UDON workflow\n# outdir: {outdir}\n# started: {_now()}\n\n")

    udon_out = os.path.join(outdir, "udon")
    cmd1 = [PY, "run_udon_pseudobulk.py", "--pseudobulk", args.pseudobulk, "--counts-pseudobulk", args.counts,
            "--sample-metadata", args.sample_metadata, "--output-dir", udon_out,
            "--species", args.species, "--control-mode", "matched"]
    if args.cell_type:
        cmd1 += ["--cell-type", args.cell_type]
    if args.no_gene_filter:
        cmd1 += ["--no-gene-filter"]
    if args.fast:
        cmd1 += ["--fast"]
    ok1 = run_step("01_udon_core", cmd1, logs, env, timeline, master)

    clusters = os.path.join(udon_out, "udon_core", "udon_clusters.txt")
    if ok1 and os.path.exists(clusters):
        sat_f = os.path.join(outdir, "satay_fisher")
        cmd2 = [PY, "satay_udon_core.py", "--metadata", args.metadata, "--udon-clusters", clusters,
                "--outdir", sat_f, "--heatmaps"]
        if args.mean_binarize:
            cmd2 += ["--mean-binarize", args.mean_binarize]
        run_step("02_satay_fisher", cmd2, logs, env, timeline, master)

        sat_c = os.path.join(outdir, "satay_cmh")
        cmd3 = [PY, "satay_udon_core.py", "--metadata", args.metadata, "--udon-clusters", clusters,
                "--outdir", sat_c, "--batch-key", args.batch_key]
        if args.mean_binarize:
            cmd3 += ["--mean-binarize", args.mean_binarize]
        run_step("03_satay_cmh", cmd3, logs, env, timeline, master)
    else:
        master.write("!! UDON core failed or produced no clusters; SATAY steps skipped.\n")
        print("!! UDON core failed; SATAY steps skipped.")

    # timeline + summary
    import pandas as pd
    tl = pd.DataFrame(timeline)
    tl.to_csv(os.path.join(logs, "timeline.tsv"), sep="\t", index=False)
    total = time.time() - t_all
    expected = {
        "udon/udon_core/udon_clusters.txt": os.path.join(udon_out, "udon_core", "udon_clusters.txt"),
        "udon/marker_heatmap.pdf": os.path.join(udon_out, "marker_heatmap.pdf"),
        "udon/goelite/GOElite_UDON.tsv": os.path.join(udon_out, "goelite", "GOElite_UDON.tsv"),
        "udon/donor_cluster_heatmap.pdf": os.path.join(udon_out, "donor_cluster_heatmap.pdf"),
        "satay_fisher/satay_cluster_enrichment.tsv": os.path.join(outdir, "satay_fisher", "satay_cluster_enrichment.tsv"),
        "satay_fisher/satay_byfield.pdf": os.path.join(outdir, "satay_fisher", "satay_byfield.pdf"),
        "satay_cmh/satay_cluster_enrichment.tsv": os.path.join(outdir, "satay_cmh", "satay_cluster_enrichment.tsv"),
    }
    with open(os.path.join(outdir, "WORKFLOW_SUMMARY.txt"), "w") as fh:
        fh.write(f"UDON + SATAY-UDON workflow\noutdir: {outdir}\nfinished: {_now()}\n")
        fh.write(f"total wall-clock: {total:.1f}s ({total/60:.1f} min)\n\n=== step timeline ===\n")
        fh.write(tl.to_string(index=False) + "\n\n=== expected outputs ===\n")
        for rel, path in expected.items():
            exists = os.path.exists(path)
            size = (os.path.getsize(path) if exists else 0)
            fh.write(f"  [{'OK ' if exists and size > 0 else 'MISSING'}] {rel}  ({size} bytes)\n")
    master.write(f"\n# finished: {_now()}  total {total:.1f}s\n"); master.close()

    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write(f"""# UDON + SATAY-UDON workflow run

Generated by `run_workflow.py` on {_now()} (total {total/60:.1f} min). The exact command is in
`command.txt`; per-step logs + timing are in `logs/` (`logs/timeline.tsv`, `logs/<step>.log`).

## Steps
| # | step | what | output |
|---|------|------|--------|
| 1 | UDON core | matched sex+platform controls, protein-coding gene filter (drop RPL/RPS, XIST/TSIX, keep Y), NMF (n_run=5), MarkerFinder, GO-Elite | `udon/` |
| 2 | SATAY Fisher | donor-level Fisher per (cell type x cluster) of the `--metadata` covariates | `satay_fisher/` |
| 3 | SATAY CMH | batch-stratified Cochran-Mantel-Haenszel (controls for `--batch-key`) | `satay_cmh/` |

## Key outputs
- `udon/udon_core/udon_clusters.txt` — UDON cluster per pseudobulk (`celltype__Sample`)
- `udon/marker_heatmap.pdf`, `udon/goelite/GOElite_UDON.tsv`
- `udon/donor_cluster_heatmap.pdf`, `udon/celltype_cluster_heatmap.pdf`, `udon/donor_specificity_per_cluster.tsv`
- `satay_fisher/satay_cluster_enrichment.tsv` (covariate, celltype, cluster, p, fdr, odds_ratio) + `satay_byfield.pdf`
- `satay_cmh/satay_cluster_enrichment.tsv` (batch-controlled p/fdr)

## Reproduce / validate
```
python run_workflow.py ...            # see command.txt
python validate_workflow_outputs.py {os.path.basename(outdir)}
```
Metadata format: `SATAY_METADATA_FORMAT.md`. SATAY engine validated against the original pyudon
(`validate_satay_core.py`).
""")
    print(f"\n[workflow] DONE in {total:.1f}s ({total/60:.1f} min). "
          f"Summary: {outdir}/WORKFLOW_SUMMARY.txt ; timeline: {logs}/timeline.tsv ; README: {outdir}/README.md")


if __name__ == "__main__":
    main()
