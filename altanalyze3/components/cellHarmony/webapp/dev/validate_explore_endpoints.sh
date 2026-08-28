#!/bin/bash
# Sweep every Explore endpoint of a running cellHarmony-web server and print the
# status code and the wall time of each. Give the job id as the first argument.
#   ./validate_explore_endpoints.sh 8677866414b84a07a0b66a56fb0d1321 [http://127.0.0.1:8000]
J=${1:?give a job id}
BASE=${2:-http://127.0.0.1:8000}
OUT=$(mktemp)
printf "%-58s %s\n" "ENDPOINT" "CODE / TIME"
for q in "status" "genes" "plot-variables" "display-filters" \
         "umap?modality=rna" "umap?modality=rna&color_by=cell_type" "umap?modality=rna&coords=X_umap" \
         "umap?modality=rna&x_field=pct_counts_mt&y_field=pct_counts_mt" \
         "expression?gene=Mpo&modality=rna" \
         "dotplot" "dotplot?genes=Mpo,Elane" "combplot?min_cells=5" "combplot?genes=Mpo&min_cells=5" \
         "chat-examples" "marker/heatmap.tsv?modality=rna"; do
  printf "%-58s %s\n" "$q" "$(curl -s -o "$OUT" -w '%{http_code} %{time_total}s' --max-time 300 "$BASE/api/jobs/$J/$q")"
done
rm -f "$OUT"
