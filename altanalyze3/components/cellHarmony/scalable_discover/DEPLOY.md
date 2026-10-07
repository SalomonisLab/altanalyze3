# Deploying scALABLE-discover

scALABLE-discover is a separate container from scALABLE-web. It clusters uploaded
single-cell data with ICGS3 and serves scALABLE's Explore and Chat views. It needs no
reference atlas and no S3 data. The image holds everything it reads.

| Item | Value |
| --- | --- |
| Repository | `SalomonisLab/altanalyze3`, commit `d0a7da4` or later |
| Build file | `altanalyze3/components/cellHarmony/scalable_discover/Dockerfile` |
| Compose file | `altanalyze3/components/cellHarmony/scalable_discover/docker-compose.discover.yml` |
| Container, image | `scalable-discover` |
| Host port | `127.0.0.1:8007`, container port 8000 |
| Public path | `/scalable-discover/` |
| Persistent volume | `./jobs` beside the compose file, mounted at `/srv/scalable-discover/jobs` |
| Memory | `mem_limit: 30g`, the same as scALABLE-web; 2 analysis workers |

## 1. Build and start

    git clone --depth 1 https://github.com/SalomonisLab/altanalyze3.git
    cd altanalyze3/altanalyze3/components/cellHarmony/scalable_discover
    docker network create lungmap_default 2>/dev/null || true
    docker compose -f docker-compose.discover.yml up -d --build

The compose file joins the `lungmap_default` network, as scALABLE-web's does. The command
`docker network create` does nothing when the network already exists.

Put `./jobs` on real disk. Uploads and results land there. The build ignores that
directory, so a rebuild never copies old jobs into the image.

Rebuild the Discover image for the disk-import handoff and QC/UMAP status fixes. No new
environment variable is required. Multiple H5AD uploads require AnnData 0.12 or newer,
already specified in `requirements.docker.txt`. The 1 GiB import trigger refers to
uncompressed arrays, separately from the 10,000-cell Explore-bundle threshold.

If port 8007 is taken on the host, change `127.0.0.1:8007:8000` in the compose file and use
the new port below.

## 2. Check the container

    docker inspect --format '{{.State.Health.Status}}' scalable-discover   # healthy
    curl -s http://127.0.0.1:8007/openapi.json | jq '.paths | length'       # 41
    curl -s -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8007/        # 200

The health check needs about 40 s after start. A count other than 41 means a stale image:
rebuild with `--build`.

## 3. Proxy

Route the public path to the container and pass the prefix through, as for `/scalable/`:

    ProxyPass        /scalable-discover/ http://127.0.0.1:8007/scalable-discover/
    ProxyPassReverse /scalable-discover/ http://127.0.0.1:8007/scalable-discover/

The app carries the same prefix in `SCALABLE_DISCOVER_ROOT_PATH`. A proxy that strips the
prefix breaks every page link. Check:

    curl -s https://<site origin>/scalable-discover/openapi.json | jq '.paths | length'   # 41

## 4. Chat

Chat calls the LungMAP.net assistant. The compose file sets

    CELLHARMONY_ASSISTANT_URL=http://site:8001/lungmap.net/api/assistant/viewer-intent

which reaches the site container on `lungmap_default`, as scALABLE-web does. If the site runs
on another host, set this variable to that host's
`/lungmap.net/api/assistant/viewer-intent` URL. Without the assistant, Chat answers HTTP 503
to questions it cannot parse itself; the other tabs work.

## 5. First job

Open `https://<site origin>/scalable-discover/`, choose Human, upload one Cell Ranger `.h5`
file and click `Save QC and run`. The status line shows each stage: QC, then ICGS3 steps 1
to 10, then results. On a laptop, 3,126 cells took 44 s and 8,485 cells took 134 s. When it
finishes, Explore opens on the predicted cell states.

A first job on the server is the first run of this image. Report the job's
`jobs/<job id>/logs/pipeline.log` if it fails.

## 6. Operation

| Task | Command |
| --- | --- |
| Update | `git pull`, then `docker compose -f docker-compose.discover.yml up -d --build` |
| Restart | `docker compose -f docker-compose.discover.yml restart` |
| Server log | `docker logs --tail 200 scalable-discover` |
| One job's log | `jobs/<job id>/logs/pipeline.log` |
| Purge old jobs | `docker exec scalable-discover python /app/altanalyze3/components/cellHarmony/webapp/cleanup_jobs.py --job-root /srv/scalable-discover/jobs --retain-days 7 --keep-latest 5` |

Discover uses scALABLE-web's retention policy: completed, failed and cancelled jobs
last updated more than 8 hours ago are deleted, including uploads, outputs and logs.
Unfinished jobs and jobs with a queued/processing differential are retained. Cleanup
runs on each new upload, at server startup and every hour while the server runs.
The hourly sweep runs outside the request event loop; cleanup errors are logged and
retried at the next sweep. No cron job or additional dependency is required.
Rebuild/restart the Discover service to activate the startup/hourly sweeps.

The manual purge command uses its explicit seven-day policy; add `--dry-run` to
list what it would remove. It does not override the automatic eight-hour policy.

### Normalized inputs and accelerated UMAP update

Rebuild the Discover image to include the normalized-input handoff fix, the UMAP fitting
control and `clustering/{umap_fit,umap_input,accelerated_umap}.py`. No additional dependency
or environment variable is required. Discover defaults to a 15-neighbor correlation
fit on the complete final MarkerFinder panel, targeting 30,000 landmarks and at least
200 cells per cluster (all cells from smaller clusters). Additional places are allocated
proportionally and sampled within clusters; the budget expands if needed for minimum
coverage. The Run tab also offers the original full feature
fitting. ICGS3's API and CLI defaults remain full fitting. Existing completed jobs retain
their saved coordinates; rerunning a job uses its saved choice, or accelerated fitting if
it predates this option. The failed normalized-input jobs can be rerun from their uploads.
The saved UI/API `landmark` choice now uses feature landmarks directly. Rebuild the image
and restart Discover to apply the profile and sampling update. No new dependency,
environment setting or result migration is needed; existing coordinates remain saved.

Include the BioMarkers identifier fix in `clustering/ICGS.py` and the packaged
`clustering/biomarkers/{Hs,Mm}/Ensembl-BioMarkers.txt.gz` catalogs. Ensembl-indexed uploads
now use the catalog's Ensembl column, with version suffixes normalized for lookup;
symbol-indexed uploads continue to use symbols. Unknown prediction labels use lowercase
`c18`. Existing jobs retain their coordinates and raw categories, with abbreviated UMAP
and filter-menu text; rerun to recompute biological predictions and gene-symbol aliases.
Include `cellHarmony/webapp/feature_lookup.py` and its shared app integration. Explore,
gene-set views and Chat accept the uploaded symbol aliases and Ensembl IDs without
renaming expression features. The browser offers one readable suggestion per feature;
all primary IDs remain accepted. Ambiguous symbols require an exact feature ID.

## 7. Settings

The compose file sets these. Change them there, then rerun `up -d`.

| Variable | Compose value | Effect |
| --- | --- | --- |
| `SCALABLE_DISCOVER_ROOT_PATH` | `/scalable-discover` | public path prefix |
| `SCALABLE_DISCOVER_JOB_STORAGE` | `/srv/scalable-discover/jobs` | job folders inside the container |
| `SCALABLE_DISCOVER_EXPORTS` | `minimal` | `minimal`: final h5ad and MarkerFinder files; `full`: every ICGS3 file and archive |
| `CELLHARMONY_JOB_WORKERS` | `2` | analyses that run at once |
| `CELLHARMONY_WORKER_MEMORY_LIMIT_GIB` | `15` | hold new starts while a worker uses this much |
| `CELLHARMONY_TOTAL_MEMORY_LIMIT_GIB` | `27` | hold new starts at this container use |
| `CELLHARMONY_BUNDLE_MIN_CELLS` | `10000` | build a disk-backed store at or above this size |
| `CELLHARMONY_ASSISTANT_URL` | see section 4 | Chat intent router |

## 8. Status of this image

The image has not been built yet; your first build is its first. Everything else was
tested locally on 2026-10-05 and 2026-10-06:

- the app on two human and one mouse dataset;
- its Python requirements in a clean environment;
- every Explore view in a browser.

The README in this directory lists those checks.
