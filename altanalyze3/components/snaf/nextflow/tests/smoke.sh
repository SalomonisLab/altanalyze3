#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
work=${SNAF_TEST_WORK:-$(mktemp -d)}
mkdir -p "$work"
python3 "$here/make_fixture.py" --outdir "$work/data"
cd "$work"
nf=${NEXTFLOW_BIN:-nextflow}
"$nf" run "$here/../main.nf" -ansi-log false \
  --mode bam --input "$work/data/samples.csv" --gtf "$work/data/reference.gtf" \
  --stop_after_counts --outdir "$work/counts" "$@"
python3 - "$work" <<'PY'
import csv,sys
from pathlib import Path
p=Path(sys.argv[1]);rows=list(csv.DictReader((p/'counts/filter_counts/filtered.tsv').open(),delimiter='\t'))
assert len(rows)==1,rows
assert rows[0]['S1']=='20' and rows[0]['S2']=='19',rows
assert rows[0]['UID'].startswith('ENSG00000000001:E1.1-E2.1'),rows
assert (p/'counts/index_reference/reference/gene_model.bed.gz').exists(), 'shared reference was deleted'
introns=list(csv.reader((p/'counts/count_introns/S1_intcounts.bed').open(),delimiter='\t'))
assert len(introns)==2 and all(int(row[4])==3 for row in introns), introns
PY
"$nf" run "$here/../main.nf" -ansi-log false -stub-run \
  --mode bam --input "$work/data/samples.csv" --gtf "$work/data/reference.gtf" \
  --db_dir "$work/data/snaf_db" --infer_hla --hla "$work/data/hla.tsv" \
  --outdir "$work/full_stub" "$@"
printf 'Validation artifacts: %s\n' "$work"
