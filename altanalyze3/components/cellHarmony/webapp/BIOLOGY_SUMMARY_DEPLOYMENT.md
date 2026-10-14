# AI biological summary deployment

## Expanded summary and native Discover interpretation (2026-10-11)

Also package `biology_sections.py`, `biology_literature.py`, `drug_chat.py` and
the complete `cellHarmony/discover_resources` Python package, including
`source_manifest.json`. The source manifest records the authoritative LungMAP
Discover file hashes; ten core DrugMatch selection/scoring/pagination definitions
have pinned AST parity checks. No external LungMAP checkout is imported at runtime.

Mount the existing source-signature resource as a read-only volume and set:

```yaml
CELLREF2_DRUGSTORE: /resources/cellref2/drugstore_v1
```

This is **the built Discover `drugstore_v1` bundle**, including its manifest,
partition metadata, gene axes and effect arrays. The tested bundle was
`drugstore_v1_20260908_083141` (about 5.6 GB on disk), with 323,661 eligible
treatment signatures. Do not substitute the fastComm ligand/receptor resource:
that resource has a different purpose. The DrugMatch data volume is external to
the Git source package. An optional `SCALABLE_INTERPRET_CACHE` sets the writable
native enrichment-cache directory; otherwise it lives under `webapp/jobs`.

### Mounting DrugMatch on the scALABLE server

The Python image contains the scorer, not the approximately 5.6 GB source library.
Rebuilding Git code alone does not install it. The verified production source is
`/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/cellref2/drugstore_v1/`.
Transfer that **entire directory**, including all four partitions, to the server;
do not use the separate `lungmap-data-test` or `lungmap-data-s3test` fixtures.

From `cellHarmony/webapp`, set `SCALABLE_DRUGSTORE_HOST` to its existing absolute
directory on the server, then add the supplied read-only mount override:

```bash
docker compose -f docker-compose.scalable.yml -f docker-compose.drugmatch.yml config --quiet
docker compose -f docker-compose.scalable.yml -f docker-compose.drugmatch.yml up -d --build
docker compose -f docker-compose.scalable.yml -f docker-compose.drugmatch.yml exec scalable-web python -m altanalyze3.components.cellHarmony.discover_resources.verify_drugstore
```

The override refuses to create an empty host directory for a missing source.
The check uses the native reader and validates the complete declared signature
roster, array dimensions, metadata offsets, source p-values and measured-gene
masks without loading effect matrices into RAM. The expected production version
is `drugstore_v1_20260908_083141`: 333,562 source signatures across four partitions,
of which 323,661 are treatment signatures eligible for the native drug screen.
For another service layout, apply the same read-only bind and
`CELLREF2_DRUGSTORE` environment value to that service; the optional override
targets the repository's `scalable-web` service. A configured path alone is
insufficient unless the complete source library exists inside the container.

The complete production bundle is already published at:

```text
s3://lungmap-breath-data/deploy/lungmap/full/cellref2/drugstore_v1/
```

Verified on 2026-10-13: all 25 objects (6,045,138,558 bytes) matched freshly
computed local SHA-256 hashes against S3 `sha256` object metadata and byte sizes.
No re-upload was required. On the deployment host, sync that existing prefix
into the host directory before applying the mount override:

```bash
aws s3 sync s3://lungmap-breath-data/deploy/lungmap/full/cellref2/drugstore_v1/ "${SCALABLE_DRUGSTORE_HOST:?Set the absolute destination directory}/"
```

Keep AWS credentials out of Compose files and Git; use the host's normal AWS
credential chain or instance role for the download.

`/chat-drugs` and `/chat-interpret` return 202 while building and 200 when ready,
using the bounded summary worker/cache. Compound matching uses local partition
memmaps. Interpret requires outbound HTTPS to ToppGene, PubChem, openFDA,
ClinicalTrials.gov and NCBI E-utilities. PubMed uses the native rate limiter.
Unavailable retrieval is exposed rather than reported as no evidence/trials.
There is no requirement to increase Apache ProxyTimeout.

Methods appear first; sections and per-state/model details are expandable.
Per-state regulatory synthesis binds identical comparison arms and state layers.
GO-Elite list enrichment, hub-target ToppGene enrichment, and RNA/lipid/metabolite
diagram overlays are separate evidence types. Both-mode reference-gap screening
retains the complete cluster ledger and reports purity, unassigned-cell
concentration, GO-Elite names and supervised assignments. Its numeric thresholds
are explicit descriptive criteria, not new significance tests.

On saved job `537ad55735014adda24cb98623882809` (7,445 cells, 46 states,
six completed modality comparisons), native GRN review took 6.78 s. DrugMatch
scored all 323,661 eligible signatures in 21.88 s with all supported genes and
2.71 s with measured genes; maximum process RSS was about 1.04 GiB. These timings
do not establish performance for a million-cell analysis. Live checks also
retrieved complete hub-target ToppGene mappings, saved GO-Elite interpretation,
native pathway interpretation and exact-compound disease-specific trial results.

**Unresolved provenance:** several saved model records lack release-specific
training-dataset citations and independent validation/confidence reports. Methods
show this gap and retain the executed bundle/version; current registry metadata
is not retroactively substituted. The missing-source question remains pending.
Literature title/abstract searches retrieve citation candidates, with documented
state aliases and broader disease-context searches. Hits do not certify an exact
species/state/direction claim, and zero hits do not establish novelty. The
“Verified known biology” section therefore does not falsely certify these hits;
claim-level verification remains incomplete.

Audited observations can be stored under `reviewed_observations` in the biology
resource. A record must identify species, state, disease, feature, modality and
direction (`up`/`down`), have `status: reviewed`, and retain `reviewed_by`,
`reviewed_at`, `claim`, `evidence_excerpt` and an HTTPS `citation_url`.
Only an exact identity and direction match can enter “Verified known biology”.
No audited findings have been fabricated or approval status inferred for this
deployment. This gate is tested against wrong genes, species, states and direction.

Commit/package `webapp/biology_summary.py`, `webapp/biology_synthesis.py`,
`webapp/biology_questions.py`, `webapp/biology_model.py`,
`webapp/biology_interpretation.py`,
`webapp/static/ai-summary.gif`, `webapp/static/ai-summary.png`, the updated Chat UI and background queue,
and **`cellHarmony/resources/biology/lungmap.json`** (included in setuptools package
data). Rebuild the scALABLE web image. No reference/model retraining, database
migration or proxy timeout increase is required. Existing jobs are supported.

The existing LungMAP `/api/assistant/viewer-intent` endpoint classifies questions;
it does **not** generate biological summaries. A runnable local completion sidecar
is included. In a separate model-service environment with `llama-cpp-python>=0.3.34`,
mount the installed GGUF weights and run:

```sh
python /app/altanalyze3/components/cellHarmony/webapp/biology_model.py \
  --model /models/qwen2.5-7b-instruct-q4_k_m-00001-of-00002.gguf \
  --model-name Qwen2.5-7B-Instruct-Q4_K_M --host 0.0.0.0 --port 8086
```

Both shards must be present for a split GGUF model. The service does not download
weights. The optional package extra is `biology-summary`. CPU execution is the
portable default; `--gpu-layers -1` enables available GPU offload. Only expose the
service on the internal application network. `/health` reports the configured model
and loaded state. Native inference runs in a child process, is serialized and unloads
after five minutes idle; a native abort or timeout restarts that child on the next call.
It cannot terminate scALABLE or an analysis worker.

Configure the scALABLE container:

```yaml
SCALABLE_SUMMARY_URL: http://local-model:8086/v1/chat/completions
SCALABLE_SUMMARY_MODEL: Qwen2.5-7B-Instruct-Q4_K_M
SCALABLE_SUMMARY_PROTOCOL: openai
```

An Ollama service can use `http://local-model:11434/api/chat` and protocol `ollama`.
An optional `SCALABLE_SUMMARY_API_KEY` is sent only as an Authorization header and
is never included in the report/cache identity. No model URL is user-controlled.
Choose a locally hosted instruction model with sufficient context for the evidence
digest. The web application does not download or load model weights. Allocate the
model service's memory separately from the 30 GiB analysis container.

The framework is code-owned. An explicit physiological question catalogue covers
host defense, surfactant, lipid handling, matrix remodeling, barriers, mucociliary
clearance, development/repair and cellular stress. These are question categories,
not substitute gene sets. Native directional GO-Elite enrichment defines which
questions can actually be evaluated. Exact saved RNA calls are recovered for
retained pathway genes and regulator targets, rather than allowing ranking truncation
to hide their support. Every join uses the original state layer and comparison arms.
Full overlap-gene lists are retained.

The ledger records eligible contrasts, original units/parameters, missing modality
comparisons, replication/confounder limits, RNA/pathway anchors, target-specific GRN
support, predicted follow-ups, contradictions and distinguishing experiments. The
priority hierarchy uses supported regulator–target–pathway chains, named RNA support,
the explicit physiological catalogue and saved pathway FDR. It is an editorial
hierarchy, not a new significance score. Counting dependent predicted modalities
cannot raise confidence. All observed states and native results remain represented.
A negative/absent summary question is not evidence of no biological effect.

Code owns the biological question ledger, comparisons and evidence links. The
first model call selects experimental priorities from the defined alternatives.
A separate model phase proposes biological systems interpretations above the
GRN, gene-set enrichment and multimodal pathway sections (and communication when
supported). It receives native hub–target relationships, target RNA directions,
TF/edge predictions, exact directional pathway overlaps and distinct molecular
nodes from the displayed diagrams. Rendering coordinates and repeated diagram
glyphs are omitted from the prompt; saved sources and all per-state results remain
in the report. Introductions use bounded leading results, with direct result links.

Code binds each section to its supplied sources and states the observed results.
The model returns a biological implication framed as a testable question and
suggestions for alternatives and experiments. Only the implication is promoted
to the report: alternatives and functional tests are code-defined. Raw model
suggestions remain in the audit record, since remeasuring RNA does not validate
an RNA-derived molecular prediction. It cannot rewrite the
observations or declare untested contrasts significant. Source/state/feature and
comparison guards reject borrowed evidence, invented quantitative claims and known
unsupported signed-regulation or signaling claims. Biological consequences are
hypotheses, not validation. These checks cannot prove that every interpretation is
biologically correct. Provider failures do not echo credentials or discard results.
Each call allows 240 seconds, with 2,200 tokens for priorities. The interpretation
phase makes one request per leading native source, capped at 900 tokens each;
isolated state/source contexts prevent findings being borrowed by another state.
Code adds functional follow-up checks from the existing physiological question
catalogue, keyed to saved enriched processes; these are proposed assays, not new
analyses or evidence. The two roles
are labelled separately in the UI. Summary polling stays
asynchronous; the local model is cached between calls and unloaded after inactivity.

Software reproducibility retains its JSON in both the browser and methods export.
Other methods use manuscript prose and rendered tables; raw timing/scale logs stay
in the ZIP as original logs, with machine-readable evidence in run_manifest.json.
New execution snapshots capture the package version and, when available, the Git
commit/date and modified-checkout status. Historical runs are not backfilled with
current checkout metadata. Component names are attributed to altanalyze3.

An actual installed-model audit rejected free-prose responses that inferred signed
regulation, surface-protein function or measured signaling from predictions. The
closed-choice contract remains in place for the hypothesis ledger; section prose
is separately checked against native result sources. On the 7,445-cell six-comparison
job, a real closed-choice call completed in 22.5 seconds with all six decisions.
The evidence collector took about 0.75 seconds and generated 12 directional questions,
retaining all 46 observed states. This is a tested candidate framework, not a validated
universal mechanistic discovery method.

Local tests used the already installed Qwen2.5-7B-Instruct Q4_K_M and Metal offload.
Model memory is additional to the analysis container budget; allocate it separately.
The CPU deployment's performance has not been benchmarked. No native analysis,
imputation model, sample roster or baseline feature panel is changed by this summary.

Bundle-backed viewers without a retained source H5AD report missing observation
coverage explicitly, while retaining available saved marker and comparison evidence.

The browser posts `/api/jobs/{job}/chat-summary`; 202 responses mean it is preparing,
200 carries the report, 409 means the analysis has not completed. One background
summary worker deduplicates requests and keeps four cached/in-flight entries per
server process; plots have an independent queue. Cache entries expire after an hour
and job/source/model configuration changes invalidate results. Across multiple web
processes, each has its own queue; configure the local model service's concurrency
limit to prevent overlapping inference allocations.

Validation covers both provider protocols with controlled responses, every supported
molecular modality, both state layers, correct saved comparison routing, unavailable
or invalid model output, and million-cell metadata-only extraction. A real installed model also passed the closed-choice contract. Neither contract tests
nor that call establish the biological accuracy of the code-defined hypotheses.
