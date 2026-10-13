"""Deterministic capped cell selection, using native scores and exact barcodes.

This module changes selection only, not normalization, model inference or tests.
Reference scores come from cellHarmony's recorded assignment artifacts and
unsupervised scores from ICGS3's saved winning SVM decisions; absent
or conflicting scores are never reconstructed, imputed or replaced by random ranks.
"""
from pathlib import Path
import json

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem


def select_cells(obs, population_col, sample_col, selected_samples, limit=500, scores=None):
    """Keep highest-scoring cells per state/sample; break score ties by barcode."""
    if not 1 <= int(limit) <= 500:
        raise ValueError("scALABLE differential cell limit must be between 1 and 500.")
    limit = int(limit)
    if not obs.index.is_unique or not obs.index.astype(str).is_unique:
        raise ValueError("Differential selection requires unique cell barcodes; no cells were discarded.")
    positions = np.flatnonzero(obs[sample_col].astype(str).isin([str(s) for s in selected_samples]))
    annotations = obs.iloc[positions][[population_col, sample_col]].astype(str)
    groups = annotations.groupby([population_col, sample_col], sort=True, observed=True).indices.values()
    if scores is None and 'AlignmentScore' in obs:
        scores = obs['AlignmentScore']
    if isinstance(scores, pd.Series):
        if not scores.index.astype(str).is_unique:
            raise ValueError("Differential score barcodes must be unique.")
        scores = scores.copy(deep=False)
        scores.index = scores.index.astype(str)
        scores = scores.reindex(obs.index.astype(str))
    values = (pd.to_numeric(pd.Series(np.asarray(scores)), errors='coerce').to_numpy(dtype=float, na_value=np.nan)
              if scores is not None else np.full(len(obs), np.nan))
    if len(values) != len(obs):
        raise ValueError("Differential scores must match the complete observation roster.")
    barcodes = obs.index.astype(str).to_numpy()
    chosen = []
    for indices in groups:
        rows = positions[indices]
        if len(rows) > limit:
            if not np.isfinite(values[rows]).all():
                missing = barcodes[rows[~np.isfinite(values[rows])]]
                raise ValueError("Native scores are missing or nonfinite for capped differential cells: "
                                 + ', '.join(missing[:5]) + ". No random selection or score substitution was performed.")
            rows = rows[np.lexsort((barcodes[rows], -values[rows]))[:limit]]
        chosen.extend(rows)
    chosen = np.asarray(chosen, dtype=np.int64)
    # Canonical analysis order also stabilizes floating-point aggregation when
    # the same input cells are stored in a different row order.
    return chosen[np.argsort(barcodes[chosen], kind='stable')]


def needs_cap(obs, population_col, sample_col, selected_samples, limit=500):
    selected = obs[sample_col].astype(str).isin([str(s) for s in selected_samples])
    return bool((obs.loc[selected].groupby([population_col, sample_col], observed=True).size() > int(limit)).any())


def reference_scores(meta, barcodes):
    """Resolve the native reference scores for all supplied cell IDs, including Both.

    Both's top-level assignments are the merged labels, so read the supervised
    BranchStore metadata for its original assignment artifact instead. Below-cutoff
    scores are retained in the explicitly registered unaligned-cells H5AD.
    """
    source = meta
    if meta.get('analysis_mode') == 'both':
        combined = Path(meta['artifacts']['combined_h5ad'])
        branch = combined.parent / 'branches' / 'supervised' / meta['job_id'] / 'job.json'
        if not branch.is_file():
            raise ValueError(f"Native supervised score provenance is missing: {branch}")
        source = json.loads(branch.read_text())
    artifacts = source.get('artifacts') or {}
    assignment_path = artifacts.get('assignments')
    if not assignment_path or not Path(assignment_path).is_file():
        raise ValueError("The native cellHarmony assignment artifact is required for score-ranked selection.")
    assignments = pd.read_csv(assignment_path, sep='\t', usecols=['CellBarcode', 'AlignmentScore'],
                              dtype={'CellBarcode': str})
    scores = assignments.set_index('CellBarcode')['AlignmentScore']
    origins = [str(assignment_path)]
    below = artifacts.get('unaligned_cells_h5ad')
    if below and Path(below).is_file():
        with h5py.File(below, 'r') as f:
            obs = f['obs']
            index = pd.Index(read_elem(obs[obs.attrs.get('_index', '_index')]).astype(str))
            values = read_elem(obs['AlignmentScore'])
            scores = pd.concat([scores, pd.Series(values, index=index)])
        origins.append(str(below))
    if not scores.index.is_unique:
        raise ValueError("Native reference score artifacts contain duplicate cell barcodes; no records were discarded.")
    return scores.reindex(pd.Index(barcodes).astype(str)), origins


def web_assignment_scores(meta, barcodes):
    """scALABLE web policy only; standalone cellHarmony does not call this.

    Both retains reference-score selection even when grouping by ICGS3 labels.
    Read unsupervised scores from the authoritative RNA output, so all modality
    sidecars select the same cells without loading the expression matrix.
    """
    mode = meta.get('analysis_mode') or (meta.get('qc') or {}).get('analysis_mode', 'supervised')
    if mode != 'unsupervised':
        scores, sources = reference_scores(meta, barcodes)
        return scores, sources, 'highest_alignment_score', 'AlignmentScore'
    path = (meta.get('artifacts') or {}).get('combined_h5ad')
    if not path or not Path(path).is_file():
        raise ValueError("The native ICGS3 combined H5AD is required for SVM-ranked selection.")
    with h5py.File(path, 'r') as f:
        obs = f['obs']
        # scALABLE's documented rename of ICGS3_SVM_score is SVM_score.
        columns = [key for key in ('SVM_score', 'ICGS3_SVM_score') if key in obs]
        if not columns:
            raise ValueError("Native ICGS3 winning SVM scores are missing; no alternative scores or random ranks were used.")
        index = pd.Index(read_elem(obs[obs.attrs.get('_index', '_index')]).astype(str))
        if not index.is_unique:
            raise ValueError("Native ICGS3 score barcodes must be unique; no records were discarded.")
        scores = pd.Series(read_elem(obs[columns[0]]), index=index)
        if len(columns) == 2 and not scores.equals(pd.Series(read_elem(obs[columns[1]]), index=index)):
            raise ValueError("Saved ICGS3 SVM score aliases conflict; no score substitution was performed.")
    return scores.reindex(pd.Index(barcodes).astype(str)), [str(path)], 'highest_icgs3_svm_score', columns[0]
