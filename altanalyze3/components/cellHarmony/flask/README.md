# cellHarmony-lite Flask demo

This directory contains a lightweight Flask application that wraps the existing
`cellHarmony_lite`/`approximate_umap` workflows behind a modern web
experience. The current implementation focuses on scaffolding (upload/QC/results
views, Plotly visualisations, downloadable artifacts, etc.) and ships with a
synthetic pipeline so the UI can be exercised without long-running compute. Hook
it up to the production pipelines by replacing the placeholder logic in
`tasks.py` with calls to `cellHarmony_lite` and `approximate_umap`.

## Quick start

```bash
cd altanalyze3/altanalyze3/components/cellHarmony/flask
python3 -m venv .venv
source .venv/bin/activate
pip install flask plotly
flask --app app run --debug
```

Open http://127.0.0.1:5000 in your browser. The demo lets you:

1. Upload up to seven H5/H5AD files (each requires a unique sample name).
2. Pick a species/reference combination from `reference_config.json`.
3. Define QC thresholds and start the alignment with a progress/log view.
4. Inspect placeholder UMAP and expression plots powered by Plotly and download
   generated TSV/log artifacts.

## Integrating the real pipeline

The placeholder steps inside `tasks.JobRunner._generate_placeholder_outputs`
should be replaced with calls to the actual tooling. A typical integration looks
like:

1. Invoke `cellHarmony_lite.py` via `subprocess.run` (or by importing the
   module) with `--h5dir`, `--refdir`, `--outdir`, and all fast-mode flags
   (SoupX disabled, no metacells, no unsupervised clustering, no h5ad export).
2. After cellHarmony-lite finishes, call
   `altanalyze3.components.visualization.approximate_umap` with the generated
   TSV outputs to place the query cells onto the reference UMAP.
3. Persist the marker genes into the AnnData object and export the TSVs/logs to
   `outputs/` so the Flask endpoints can serve them.

Remember to update `reference_config.json` with the real species/reference paths
available on your cluster.

## File layout

```
flask/
  app.py                # WSGI entrypoint
  config.py             # Default app settings
  job_manager.py        # Filesystem-backed job registry
  tasks.py              # Background executor (placeholder pipeline)
  routes.py             # Flask blueprints (HTML + JSON APIs)
  reference_config.json # Species/reference metadata used by the UI
  templates/            # Jinja templates (Bootstrap layout)
  static/               # Plotly-enabled JS/CSS bundle
  jobs/                 # (gitignored) local workspace for uploads/outputs
```

## Testing

The app ships with a self-contained fake pipeline so you can validate the web
flow without heavy dependencies. Simply upload a few small files (any content)
and observe the status/progress updates and Plotly charts after the job
completes. Once integrated with the production workflow, keep an eye on the
log view to monitor each stage (`cellHarmony-lite`, marker detection, approximate
UMAP, etc.).
