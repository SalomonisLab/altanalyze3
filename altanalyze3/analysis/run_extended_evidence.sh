#!/usr/bin/env bash
set -euo pipefail
OUT=/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927/extended
mkdir -p "$OUT"
/opt/homebrew/opt/python@3.11/bin/python3.11 /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/analysis/extended_evidence_worker.py > "$OUT/worker.log" 2>&1
