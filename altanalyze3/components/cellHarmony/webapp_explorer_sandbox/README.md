# cellHarmony Explorer Sandbox

This package is an evaluation-only alternative frontend for `cellHarmony web`.

It reuses the existing backend routes and workflow logic from:

- `altanalyze3.components.cellHarmony.webapp`

but serves a different template and static bundle with:

- persistent left control rail
- tabbed explorer workspace
- separate activity/log tab

## Run locally

From the repo root:

```bash
uvicorn altanalyze3.components.cellHarmony.webapp_explorer_sandbox.app:app --reload --port 8001
```

Then open:

```text
http://127.0.0.1:8001
```

The sandbox stores its jobs separately from the production web app.
