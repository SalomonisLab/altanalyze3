# cellHarmony-lite FastAPI portal

This directory provides a browser-testable FastAPI front end for the existing
`cellHarmony_lite` and `approximate_umap` pipeline. It is intended as the
production-oriented successor to the older Flask demo: same core workflow,
cleaner API layer, easier deployment behind a reverse proxy, and a UI that can
be tested immediately in a web browser.

## Quick start

```bash
cd altanalyze3/altanalyze3/components/cellHarmony/webapp
python3 -m venv .venv
source .venv/bin/activate
pip install fastapi uvicorn python-multipart jinja2 plotly
uvicorn altanalyze3.components.cellHarmony.webapp.app:app --reload
```

Open `http://127.0.0.1:8000`.

## API docs

- Browser API: FastAPI serves interactive OpenAPI docs at `http://127.0.0.1:8000/docs`
- Dedicated approximate UMAP endpoint notes: [approximate_umap_api.md](./approximate_umap_api.md)

## Configuration

Environment variables:

- `CELLHARMONY_JOB_STORAGE`
- `CELLHARMONY_REFERENCE_REGISTRY`
- `CELLHARMONY_MAX_FILES`
- `CELLHARMONY_JOB_WORKERS`

By default the app reuses the existing reference registry at:

`altanalyze3/components/cellHarmony/flask/reference_config.json`

## Notes

- The portal reuses the current filesystem-backed `JobStore` and threaded
  `JobRunner`, so it is immediately testable without extra infrastructure.
- For a high-scale deployment, keep this FastAPI layer and swap the job runner
  to Redis/Celery or Redis/RQ plus a real persistent metadata store.
