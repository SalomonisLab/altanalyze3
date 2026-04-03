# cellHarmony-lite FastAPI portal

This directory provides a browser-testable FastAPI front end for the existing
`cellHarmony_lite` and `approximate_umap` pipeline. It is intended as the
production-oriented successor to the older Flask demo: same core workflow,
cleaner API layer, easier deployment behind a reverse proxy, and a UI that can
be tested immediately in a web browser.

Current web interface capabilities include:

- upload one or more `.h5` or `.h5ad` single-cell files as a job
- preserve compatible `.obs` metadata from uploaded `.h5ad` files
- filter the displayed cells in the UMAP and expression viewers by selected
  `.obs` fields and values
- run differential analysis using a user-selected cell-state field from `.obs`
- define numerator and denominator values from a separate `.obs` grouping field
- optionally switch differential comparison type between `cells` and
  `pseudobulk` when enough grouped values are available

## Related docs

- [cellHarmony-lite documentation](../../../../docs/cellHarmony.md)
- [cellHarmony-differential documentation](../../../../docs/cellHarmony_differential.md)
- [How to Use cellHarmony web](./HOW_TO_USE.md)

## Quick start

```bash
cd altanalyze3/altanalyze3/components/cellHarmony/webapp
python3 -m venv .venv
source .venv/bin/activate
pip install fastapi uvicorn python-multipart jinja2 plotly
uvicorn altanalyze3.components.cellHarmony.webapp.app:app --reload
```

Open `http://127.0.0.1:8000`.

## Docker for a cloud VM

This web app can run as a single container with persisted job storage. The
Docker assets live in this directory:

- [Dockerfile](./Dockerfile)
- [docker-compose.yml](./docker-compose.yml)
- [requirements.docker.txt](./requirements.docker.txt)

From `altanalyze3/altanalyze3/components/cellHarmony/webapp`:

```bash
docker compose up --build
```

Then open `http://<vm-ip>:8000`.

Notes:

- Job data is persisted to `cellHarmony/webapp/jobs` on the host via a bind mount.
- If you serve the app under a subpath such as `/cell-harmony`, set
  `CELLHARMONY_ROOT_PATH=/cell-harmony` so the frontend and API requests use
  that prefix.
- The container sets `CELLHARMONY_JOB_STORAGE=/srv/cellharmony/jobs`.
- The default reference registry is baked into the image at:
  `/app/altanalyze3/components/cellHarmony/flask/reference_config.json`
- If you deploy behind Nginx or Caddy, keep the container on port `8000` and
  terminate TLS at the reverse proxy.

Useful commands:

```bash
docker compose up --build -d
docker compose logs -f
docker compose down
```

## API docs

- Browser API: FastAPI serves interactive OpenAPI docs at `http://127.0.0.1:8000/docs`
- Dedicated approximate UMAP endpoint notes: [approximate_umap_api.md](./approximate_umap_api.md)

## Configuration

Environment variables:

- `CELLHARMONY_ROOT_PATH`
- `CELLHARMONY_JOB_STORAGE`
- `CELLHARMONY_REFERENCE_REGISTRY`
- `CELLHARMONY_MAX_FILES`
- `CELLHARMONY_JOB_WORKERS`

By default the app reuses the existing reference registry at:

`altanalyze3/components/cellHarmony/flask/reference_config.json`

## Daily retention cleanup

The web app stores each job under `cellHarmony/webapp/jobs/<job_id>` with a
`job.json` metadata file. A cleanup helper is included at:

`cellHarmony/webapp/cleanup_jobs.py`

It deletes only jobs that are safe to purge:

- top-level terminal status only: `completed`, `failed`, `cancelled`, `canceled`
- never deletes `uploaded`, `queued`, or `processing` jobs
- never deletes a job whose nested differential status is still `queued` or `processing`
- uses `updated_at` from `job.json` and falls back to `created_at`

Recommended first pass:

```bash
cd altanalyze3/altanalyze3/components
python3 cellHarmony/webapp/cleanup_jobs.py --retain-days 7 --keep-latest 5 --dry-run
```

Then enable real deletion:

```bash
cd altanalyze3/altanalyze3/components
python3 cellHarmony/webapp/cleanup_jobs.py --retain-days 7 --keep-latest 5
```

For macOS, a launchd template is provided at:

`cellHarmony/webapp/launchd/org.cellharmony.jobs-retention.plist`

Update the paths in that plist for your machine, then install it:

```bash
cp cellHarmony/webapp/launchd/org.cellharmony.jobs-retention.plist ~/Library/LaunchAgents/
launchctl unload ~/Library/LaunchAgents/org.cellharmony.jobs-retention.plist 2>/dev/null || true
launchctl load ~/Library/LaunchAgents/org.cellharmony.jobs-retention.plist
launchctl start org.cellharmony.jobs-retention
```

Check the logs at:

- `/tmp/cellharmony_jobs_retention.log`
- `/tmp/cellharmony_jobs_retention.err`

## Notes

- The portal reuses the current filesystem-backed `JobStore` and threaded
  `JobRunner`, so it is immediately testable without extra infrastructure.
- For a high-scale deployment, keep this FastAPI layer and swap the job runner
  to Redis/Celery or Redis/RQ plus a real persistent metadata store.
