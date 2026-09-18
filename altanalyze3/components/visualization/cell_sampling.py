"""Reproducible sampling of individual observations for shared plot displays."""
from hashlib import blake2b

import numpy as np


CELL_SAMPLE_LIMITS = (5, 10, 20, 50)


def sample_cell_indices(cell_ids, sample_labels, *, limit=10, group_labels=None):
    """Return indices of up to ``limit`` cells per sample (and optional group).

    Hash ranking makes selections reproducible and nested: the five-cell selection
    is contained in the ten-cell selection. Return indices in input order to keep
    the plot's existing group/cell ordering. Zero means all observations.
    """
    if limit not in (0, *CELL_SAMPLE_LIMITS):
        raise ValueError("Cells per sample must be 5, 10, 20, 50, or 0 for all cells.")
    ids = list(map(str, cell_ids))
    samples = list(map(str, sample_labels))
    groups = [""] * len(ids) if group_labels is None else list(map(str, group_labels))
    if len(ids) != len(samples) or len(ids) != len(groups):
        raise ValueError("Cell identifiers and sampling annotations must align.")
    if not limit:
        return np.arange(len(ids), dtype=np.int64)
    strata = {}
    for index, (cell, sample, group) in enumerate(zip(ids, samples, groups)):
        rank = blake2b(cell.encode("utf-8"), digest_size=16, person=b"scALABLE-cells").digest()
        strata.setdefault((sample, group), []).append((rank, cell, index))
    selected = [index for rows in strata.values() for _, _, index in sorted(rows)[:limit]]
    return np.asarray(sorted(selected), dtype=np.int64)
