"""Lossless column transport for large plots; legacy row payloads remain available."""
from __future__ import annotations


def compact_plot_payload(payload):
    result = dict(payload)
    for name in ("query", "reference", "umap", "scatter"):
        rows = result.get(name)
        if not rows or not isinstance(rows, list):
            continue
        fields = list(rows[0])
        # Only homogeneous point records: do not silently discard extra fields.
        if any(row.keys() != rows[0].keys() for row in rows):
            continue
        columns, dictionaries = {}, {}
        for field in fields:
            values = [row[field] for row in rows]
            if field in ("population", "sample"):
                labels = list(dict.fromkeys(values))
                codes = {label: i for i, label in enumerate(labels)}
                columns[field] = [codes[value] for value in values]
                dictionaries[field] = labels
            else:
                columns[field] = values
        result[name] = {"encoding": "columns-v1", "length": len(rows),
                        "columns": columns, "dictionaries": dictionaries}
    return result
