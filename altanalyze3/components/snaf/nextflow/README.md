# SNAF portable Nextflow workflow

Implements BAM → reference indexing → parallel junction/intron counts → aggregation
and annotation → ≥20-read filtering → supplied/inferred HLA → SNAF predictions,
with optional MS evidence interpretation by the independent pyNeoQuant library.
Existing tumor-specificity (`ts`), matrix-based prediction (`full`) and surface
antigen (`surface`) modes remain available.

Nextflow ≥25.04.8. Run `nextflow run main.nf -profile test,local` for a complete
small BAM-to-count-matrix test. The test data and reference are included and no
external models are downloaded. Prediction/model acceptance is a separate test.

See [deployment instructions](../../../deployment/README.md) for installation,
samplesheets, CLI integration, containers, Galaxy and testing.

Parameter definitions: `nextflow_schema.json`. Trace/report/timeline/DAG artifacts
are written beneath `--outdir`. Container names are local build tags until the
separate images are published. pyNeoQuant is optional and runs in its own image.

This workflow is prepared for further nf-core review; community acceptance and
registry publication have not been claimed.
