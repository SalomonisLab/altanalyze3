"""Small, dependency-free file contracts used by CLI, Galaxy and Nextflow."""
import csv
import hashlib
import json
from pathlib import Path


def read_table(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'Missing or duplicate columns: {path}')
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f'Ragged table: {path}')
        return rows, reader.fieldnames


def write_table(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def manifest(path, **values):
    Path(path).write_text(json.dumps({'schema_version': '1.0', **values}, indent=2, sort_keys=True)+'\n')
