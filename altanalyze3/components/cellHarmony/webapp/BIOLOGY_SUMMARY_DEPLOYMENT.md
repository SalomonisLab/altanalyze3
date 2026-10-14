# scALABLE deployment — Michal

Deploy the latest scALABLE code and rebuild the web image. Keep the existing jobs
volume, references, imputation models and 30 GiB container limit.

## DrugMatch library

The complete library is already in S3. Download it to the server and attach it
with the supplied Compose override. From `cellHarmony/webapp`:

```bash
# Choose a persistent directory on the server.
export SCALABLE_DRUGSTORE_HOST=/srv/scalable-resources/drugstore_v1
aws s3 sync s3://lungmap-breath-data/deploy/lungmap/full/cellref2/drugstore_v1/ "$SCALABLE_DRUGSTORE_HOST/"
docker compose -f docker-compose.scalable.yml -f docker-compose.drugmatch.yml up -d --build
```

The override mounts the library read-only and automatically sets
`CELLREF2_DRUGSTORE`. If deploying with a different container configuration,
mount the same library and set that variable to its path inside the container.

## AI summary service

Run the included `biology_model.py` with `llama-cpp-python` and the
Qwen2.5-7B-Instruct Q4_K_M GGUF weights in a separate internal service named
`scalable-summary`, listening on port 8086. Configure scALABLE with:

```text
SCALABLE_SUMMARY_URL=http://scalable-summary:8086/v1/chat/completions
SCALABLE_SUMMARY_MODEL=Qwen2.5-7B-Instruct-Q4_K_M
SCALABLE_SUMMARY_PROTOCOL=openai
```

This is the model/runtime tested locally. The model service and weights must be
provisioned separately; the web Compose file does not create them. The LungMAP
intent router cannot supply summary text. Keep model memory outside the 30 GiB
analysis container.

For the eight-CPU model service, pass `--threads 8` when starting `biology_model.py`.
Rebuild/restart the web app for summary progress and persistent response caching;
restart the model service to expose inference timing and thread/GPU settings.
The initial summary now reads saved evidence without model inference or literature
lookups. Expanding an interpretation requests only that result, with a loading
indicator and cached responses. Literature checks and AI prioritization are also
loaded on demand. Rebuild the web app and viewer to enable this behavior.

Rebuild the scALABLE viewer from `f6ecf18` as well: it imports the shared
cellHarmony web app. Keep its existing viewer bundles and assets.

Allow outbound HTTPS for interpretation lookups to ToppGene, PubChem, openFDA,
ClinicalTrials.gov and NCBI.

## Verify after deployment

- Open an existing completed job and run a DrugMatch query from Chat.
- Generate an AI summary and check that its result links open.
- Run a GRN-edge differential comparison and confirm that cell-state progress updates appear.

No database migration, reference/model retraining or proxy-timeout change is required.
