"""Resolve existing input storage without replacing observations or fitting models."""
import json
import re
from types import SimpleNamespace
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem, sparse_dataset


def embedding_label(key):
    name = str(key)
    if 'tsne' in name.lower() or 't_sne' in name.lower(): return 't-SNE'
    if 'umap' in name.lower(): return 'UMAP'
    return name.removeprefix('X_').replace('_', ' ')


def choose_embedding(keys, default=None):
    keys = list(keys)
    if default in keys: return default
    return next((k for k in keys if 'umap' in k.lower()), keys[0] if keys else None)


def _probe(node, cells, genes):
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    matrix = node if isinstance(node, h5py.Dataset) else sparse_dataset(node)
    rows = np.unique(np.linspace(0, cells - 1, min(cells, 128), dtype=int)) if cells else np.array([], dtype=int)
    # This is a header-time diagnostic. Execution rechecks the native classifier
    # using its established sample policy before count-space operations.
    values = matrix[rows, :] if len(rows) else np.empty((0, genes))
    return infer_expression_scale(values)['verdict']


def mitochondrial_mask(var):
    mask = pd.Index(var.index.astype(str)).str.upper().str.startswith('MT-')
    for key in ('gene_symbols', 'gene_symbol', 'feature_name', 'Gene', 'Symbol'):
        if key in var:
            mask |= var[key].astype(object).fillna('').astype(str).str.upper().str.startswith('MT-').to_numpy()
    return np.asarray(mask, bool)


MISSING_LAYOUT = ('UMAP or t-SNE coordinates missing. Please select Supervised or Unsupervised '
                  'to assign coordinates before proceeding.')


def embedding_candidates(obs, obsm):
    """Recognize named, finite supplied layouts; never reinterpret PCA as UMAP."""
    candidates = []
    for key, values in obsm.items():
        if not any(term in key.lower() for term in ('umap', 'tsne', 't_sne')):
            continue
        if isinstance(values, h5py.Group):
            values = read_elem(values)
        if isinstance(values, pd.DataFrame):
            values = values.to_numpy()
        shape = values.shape
        if len(shape) == 2 and shape[0] == len(obs) and shape[1] >= 2:
            coordinates = np.asarray(values[:, :2])
            if np.issubdtype(coordinates.dtype, np.number) and np.isfinite(coordinates).all():
                candidates.append({'key': key, 'kind': 'obsm', 'label': embedding_label(key) + ' · obsm[' + key + ']'})
    # Existing tools also export named coordinates directly into .obs. Pair
    # exact axis suffixes, retaining every source observation unchanged.
    groups = {}
    for key in obs:
        match = re.match(r'^(.*(?:umap|t[_-]?sne).*?)[_ .-]?([012xy])$', str(key), re.I)
        if match and pd.api.types.is_numeric_dtype(obs[key]):
            groups.setdefault(match[1].lower().rstrip('_ .-'), {})[match[2].lower()] = key
    for prefix, axes in groups.items():
        pair = next(([axes[x], axes[y]] for x, y in [('x', 'y'), ('0', '1'), ('1', '2')]
                     if x in axes and y in axes), None)
        if pair and np.isfinite(obs[pair].to_numpy(dtype=float)).all():
            key = 'obs:' + json.dumps(pair, separators=(',', ':'))
            candidates.append({'key': key, 'kind': 'obs', 'columns': pair,
                               'label': embedding_label(prefix) + ' · obs[' + ', '.join(pair) + ']'})
    return candidates


def inspect_storage(handle, obs, var):
    note = None
    candidates = [('layers/counts', handle.get('layers/counts'))]
    candidates += [('layers/' + key, node) for key, node in handle.get('layers', {}).items() if key != 'counts']
    raw = handle.get('raw')
    if raw is not None and 'var' in raw:
        raw_names = read_elem(raw['var']).index
        if raw_names.is_unique and len(raw_names) == len(var) and var.index.isin(raw_names).all():
            candidates.append(('raw/X', raw.get('X')))
        else:
            note = 'Original data use a different gene panel; a matching counts matrix is needed for count-based QC.'
    candidates.append(('X', handle['X']))
    count_candidates = [{'key': name, 'label': name.replace('/', '.'), 'scale': _probe(node, len(obs), len(var))}
                        for name, node in candidates if node is not None]
    counts = [c['key'] for c in count_candidates if c['scale'] == 'counts']
    layouts = embedding_candidates(obs, handle.get('obsm', {}))
    embeddings = [entry['key'] for entry in layouts]
    default = read_elem(handle['uns/default_embedding']) if 'uns/default_embedding' in handle else None
    selected = choose_embedding(embeddings, default)
    return {'count_source': counts[0] if counts else None, 'count_candidates': count_candidates,
            'count_note': note, 'embeddings': embeddings, 'embedding_candidates': layouts,
            'embedding_key': selected, 'embedding_label': embedding_label(selected) if selected else None,
            'mitochondrial_genes': int(mitochondrial_mask(var).sum())}


def select_counts(matrix, source='auto'):
    """Validate the explicit source with the native classifier before using it.

    All source matrices remain available, including a pre-existing counts layer
    when a different counts source is requested.
    """
    automatic = not source or source == 'auto'
    if automatic:
        from altanalyze3.components.clustering.ICGS import infer_expression_scale
        workspace = getattr(matrix, '_matrix_workspace', None)
        def is_counts(values):
            probe = SimpleNamespace(X=values, layers={}, n_obs=matrix.n_obs, n_vars=matrix.n_vars)
            if workspace:
                probe = workspace.scale_probe(probe)
            return infer_expression_scale(probe.X)['verdict'] == 'counts'
        sources = ['counts'] if 'counts' in matrix.layers else []
        sources += [key for key in matrix.layers if key != 'counts']
        source = next(('layers/' + key for key in sources if is_counts(matrix.layers[key])), None)
        if source is None:
            # .raw is admitted only with the complete matching gene roster.
            raw = matrix.raw
            compatible = raw is not None and raw.var_names.is_unique and len(raw.var_names) == matrix.n_vars and matrix.var_names.isin(raw.var_names).all()
            source = 'raw/X' if compatible and is_counts(raw.X) else 'X'
    if source == 'X':
        values = matrix.X
    elif source.startswith('layers/') and source[7:] in matrix.layers:
        values = matrix.layers[source[7:]]
    elif source == 'raw/X' and matrix.raw is not None:
        raw = matrix.raw
        if not raw.var_names.is_unique or len(raw.var_names) != matrix.n_vars or not matrix.var_names.isin(raw.var_names).all():
            raise ValueError('Selected raw counts must contain the complete matching gene panel; no features were removed.')
        workspace = getattr(matrix, '_matrix_workspace', None)
        order = raw.var_names.get_indexer(matrix.var_names)
        values = raw.X if raw.var_names.equals(matrix.var_names) else (workspace.select(raw.X, cols=order) if workspace else raw.X[:, order])
    else:
        raise ValueError('Selected counts source is unavailable in the uploaded input.')
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    workspace = getattr(matrix, '_matrix_workspace', None)
    probe = SimpleNamespace(X=values, layers={}, n_obs=matrix.n_obs, n_vars=matrix.n_vars)
    if workspace:
        probe = workspace.scale_probe(probe)
    if infer_expression_scale(probe.X, name=source)['verdict'] != 'counts':
        if automatic and source == 'X':
            return None  # normalized-only input is checked against requested QC below
        raise ValueError('Selected counts source is not count-scale expression. Choose a genuine counts matrix.')
    if source != 'layers/counts':
        if 'counts' in matrix.layers:
            backup = 'source_counts'
            while backup in matrix.layers:
                backup = 'source_' + backup
            matrix.layers[backup] = matrix.layers['counts']
            matrix.uns['scalable_original_counts_layer'] = backup
        workspace = getattr(matrix, '_matrix_workspace', None)
        matrix.layers['counts'] = ((workspace.copy_matrix(values) if workspace else values.copy())
                                   if values is matrix.X else values)
    matrix.uns['scalable_count_source'] = source
    print(f'[qc] Use Existing: selected {source} for count-based operations; original source arrays retained.')
    return source


def attach_raw_counts(matrix):
    """Use verified .raw counts only when the complete RNA gene roster matches."""
    if 'counts' in matrix.layers or matrix.raw is None: return None
    raw = matrix.raw
    if not raw.var_names.is_unique or len(raw.var_names) != matrix.n_vars or not matrix.var_names.isin(raw.var_names).all():
        return None
    workspace = getattr(matrix, '_matrix_workspace', None)
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    probe = workspace.scale_probe(SimpleNamespace(X=raw.X, layers={}, n_obs=matrix.n_obs, n_vars=raw.n_vars)) if workspace else raw
    if infer_expression_scale(probe.X, name='raw.X')['verdict'] != 'counts': return None
    order = raw.var_names.get_indexer(matrix.var_names)
    counts = raw.X if raw.var_names.equals(matrix.var_names) else (workspace.select(raw.X, cols=order) if workspace else raw.X[:, order])
    # No new dense copy: QC reads this matrix; normalization leaves logged X
    # intact, while ambient correction's native hook makes its own work arrays.
    matrix.layers['counts'] = counts
    matrix.uns['scalable_count_source'] = 'raw/X'
    print('[qc] Use Existing: verified raw.X counts with the complete matching gene panel; source raw data retained.')
    return 'raw/X'
