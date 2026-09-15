#!/usr/bin/env python3
"""Rename CellRef2.0 cell states to short_name and stretch the viewer X axis.

Two edits, applied in place to the three cellHarmony CellRef2.0 reference files:

1. Every cell-state label becomes its `short_name` from the CellRef2 v7 name map.
   The map key is `v7_name`, which covers 86 of 86 labels in both files.
2. `UMAP1` is multiplied by X_SCALE so the viewer spreads the cells wider on X.
   `UMAP2` does not change.

Every file is copied to a `.longname-backup-<STAMP>` sibling before the write.
The script raises on any row-count, column-count or coverage mismatch.
"""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

REF_DIR = Path(
    "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/"
    "cellHarmony/flask/references/Human/Lung/CellRef2.0"
)
MAP_TSV = Path(
    "/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/"
    "external_metadata/cellref2_v7_name_map.tsv"
)
STATES_TSV = REF_DIR / "Hs-CellRef2.txt"
CLUSTERS_TSV = REF_DIR / "Hs-CellRef2_clusters.tsv"
UMAP_TSV = REF_DIR / "Hs-CellRef2_umap.tsv"

MAP_KEY = "v7_name"
MAP_VALUE = "short_name"
X_SCALE = 1.5
STAMP = "20260904"

EXPECTED_STATE_ROWS = 5139
EXPECTED_STATE_FIELDS = 87
EXPECTED_CELL_ROWS = 90003
EXPECTED_CELL_FIELDS = 3
EXPECTED_LABELS = 86


def log(message: str) -> None:
    print(message, flush=True)


def load_name_map() -> dict[str, str]:
    with MAP_TSV.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mapping = {row[MAP_KEY]: row[MAP_VALUE] for row in rows}
    if len(mapping) != len(rows):
        raise ValueError(f"{MAP_KEY} is not unique in {MAP_TSV}")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f"{MAP_VALUE} collides in {MAP_TSV}")
    log(f"[map] {MAP_TSV}: {len(rows)} rows, {MAP_KEY} -> {MAP_VALUE}, 0 collisions")
    return mapping


def read_lines(path: Path) -> list[str]:
    text = path.read_text()
    if not text.endswith("\n"):
        raise ValueError(f"{path} has no trailing newline")
    return text.splitlines()


def backup(path: Path) -> Path:
    target = path.with_name(f"{path.name}.longname-backup-{STAMP}")
    if target.exists():
        raise FileExistsError(f"backup already exists: {target}")
    shutil.copy2(path, target)
    log(f"[backup] {path} -> {target}")
    return target


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def rename_states(mapping: dict[str, str]) -> None:
    lines = read_lines(STATES_TSV)
    if len(lines) != EXPECTED_STATE_ROWS:
        raise ValueError(f"{STATES_TSV}: {len(lines)} rows, expected {EXPECTED_STATE_ROWS}")
    header = lines[0].split("\t")
    if len(header) != EXPECTED_STATE_FIELDS:
        raise ValueError(f"{STATES_TSV}: {len(header)} header fields, expected {EXPECTED_STATE_FIELDS}")
    old = header[1:]
    missing = [name for name in old if name not in mapping]
    if missing:
        raise KeyError(f"{STATES_TSV}: {len(missing)} unmapped columns, first: {missing[:5]}")
    new = [mapping[name] for name in old]
    if len(set(new)) != len(new):
        raise ValueError(f"{STATES_TSV}: short_name collision in the header")
    backup(STATES_TSV)
    lines[0] = "\t".join([header[0]] + new)
    write_lines(STATES_TSV, lines)
    log(f"[states] {STATES_TSV}: renamed {len(new)} of {len(old)} columns, {len(lines)} rows written")


def rename_clusters(mapping: dict[str, str]) -> None:
    lines = read_lines(CLUSTERS_TSV)
    if len(lines) != EXPECTED_CELL_ROWS:
        raise ValueError(f"{CLUSTERS_TSV}: {len(lines)} rows, expected {EXPECTED_CELL_ROWS}")
    header = lines[0].split("\t")
    if header != ["barcode", "cell_state", "Population"]:
        raise ValueError(f"{CLUSTERS_TSV}: unexpected header {header}")
    out = [lines[0]]
    seen_old: set[str] = set()
    seen_new: set[str] = set()
    for line_no, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != EXPECTED_CELL_FIELDS:
            raise ValueError(f"{CLUSTERS_TSV}:{line_no}: {len(fields)} fields")
        for name in fields[1:]:
            if name not in mapping:
                raise KeyError(f"{CLUSTERS_TSV}:{line_no}: unmapped label {name!r}")
            seen_old.add(name)
        renamed = [mapping[name] for name in fields[1:]]
        seen_new.update(renamed)
        out.append("\t".join([fields[0]] + renamed))
    if len(out) != len(lines):
        raise ValueError(f"{CLUSTERS_TSV}: {len(out)} rows out, {len(lines)} rows in")
    if len(seen_old) != EXPECTED_LABELS or len(seen_new) != EXPECTED_LABELS:
        raise ValueError(
            f"{CLUSTERS_TSV}: {len(seen_old)} labels in, {len(seen_new)} labels out, "
            f"expected {EXPECTED_LABELS}"
        )
    backup(CLUSTERS_TSV)
    write_lines(CLUSTERS_TSV, out)
    log(
        f"[clusters] {CLUSTERS_TSV}: {len(out) - 1} cells renamed, "
        f"{len(seen_old)} labels in -> {len(seen_new)} labels out"
    )


def rescale_umap() -> None:
    lines = read_lines(UMAP_TSV)
    if len(lines) != EXPECTED_CELL_ROWS:
        raise ValueError(f"{UMAP_TSV}: {len(lines)} rows, expected {EXPECTED_CELL_ROWS}")
    header = lines[0].split("\t")
    if header != ["barcode", "UMAP1", "UMAP2"]:
        raise ValueError(f"{UMAP_TSV}: unexpected header {header}")
    out = [lines[0]]
    x_in: list[float] = []
    x_out: list[float] = []
    for line_no, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != EXPECTED_CELL_FIELDS:
            raise ValueError(f"{UMAP_TSV}:{line_no}: {len(fields)} fields")
        old_x = float(fields[1])
        new_x = old_x * X_SCALE
        x_in.append(old_x)
        x_out.append(new_x)
        out.append("\t".join([fields[0], repr(new_x), fields[2]]))
    if len(out) != len(lines):
        raise ValueError(f"{UMAP_TSV}: {len(out)} rows out, {len(lines)} rows in")
    range_in = max(x_in) - min(x_in)
    range_out = max(x_out) - min(x_out)
    ratio = range_out / range_in
    if abs(ratio - X_SCALE) > 1e-9:
        raise ValueError(f"{UMAP_TSV}: X range ratio {ratio}, expected {X_SCALE}")
    backup(UMAP_TSV)
    write_lines(UMAP_TSV, out)
    log(
        f"[umap] {UMAP_TSV}: UMAP1 x {X_SCALE}; "
        f"range {range_in:.4f} -> {range_out:.4f} (ratio {ratio:.6f}); "
        f"min {min(x_in):.4f} -> {min(x_out):.4f}, max {max(x_in):.4f} -> {max(x_out):.4f}; "
        f"UMAP2 unchanged"
    )


def main() -> int:
    mapping = load_name_map()
    rename_states(mapping)
    rename_clusters(mapping)
    rescale_umap()
    log("[done] all three files rewritten")
    return 0


if __name__ == "__main__":
    sys.exit(main())
