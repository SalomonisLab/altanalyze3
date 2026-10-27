"""Use supplied labels/coordinates with the standard scALABLE downstream hooks."""
import numpy as np
import pandas as pd


def prepare_annotations(matrix, meta, store, job_id):
    from .existing_annotations import inspect_annotations, upload_paths, validate_selection
    report = inspect_annotations(upload_paths(meta, store.uploads_dir(job_id)))
    selected = (meta.get('qc') or {}).get('existing_state_column') or report['selected']
    candidate = validate_selection(report, selected)
    if selected not in matrix.obs:
        raise ValueError(f'Supplied annotation column {selected!r} was not retained through QC.')
    if not matrix.obs_names.is_unique or not matrix.var_names.is_unique:
        raise ValueError('Use Existing requires unique barcode and feature identities; none have been discarded.')
    if not matrix.n_obs or not matrix.n_vars:
        raise ValueError('Requested QC retained no cells or features. Review the QC thresholds.')
    from ..scalable_discover.pipeline import attach_serving_gene_symbols
    attach_serving_gene_symbols(matrix)
    # Keep source values and dtypes intact, including missing annotations. The
    # serving column represents missing labels explicitly rather than dropping rows.
    key = 'scalable_existing_state'
    while key in matrix.obs:
        key = 'source_' + key
    source = matrix.obs[selected]
    labels = source.astype(object).where(source.notna(), 'Unassigned').astype(str)
    labels = labels.where(labels.str.strip().ne(''), 'Unassigned')
    matrix.obs[key] = pd.Categorical(labels)
    record = dict(meta.get('existing_annotations') or {}, selected=selected,
                  effective_column=key, source='supplied_obs', qc_cells=matrix.n_obs,
                  missing_cells=int((labels == 'Unassigned').sum()))
    matrix.uns['scalable_existing_analysis'] = {k: v for k, v in record.items()
        if k in {'selected', 'effective_column', 'source', 'qc_cells', 'missing_cells', 'method', 'model', 'reason'}}
    store.update_job(job_id, existing_annotations=record)
    store.append_log(job_id, f'[existing] Using supplied obs[{selected!r}] ({candidate["categories"]} groups); '
                     'supervised alignment and ICGS3 are skipped. Original annotations are retained.')
    if candidate['kind'] == 'cluster':
        store.append_log(job_id, '[existing] Supplied labels are cluster identifiers; no biological cell-type labels have been inferred.')
    return key


def supplied_embedding_result(matrix, cluster_key, store, job_id):
    from altanalyze3.components.visualization.approximate_umap import ApproximateUMAPResult
    key = next((k for k in ['X_umap'] + list(matrix.obsm) if k in matrix.obsm and 'umap' in k.lower()
                and len(matrix.obsm[k].shape) == 2 and matrix.obsm[k].shape[1] >= 2), None)
    if key:
        coordinates = np.asarray(matrix.obsm[key])[:, :2]
        if 'X_umap' not in matrix.obsm:
            matrix.obsm['X_umap'] = coordinates.copy()
        store.append_log(job_id, f'[existing] Preserving supplied embedding {key}; no UMAP is fitted or projected.')
    else:
        coordinates = np.full((matrix.n_obs, 2), np.nan)
        store.append_log(job_id, '[existing] WARNING: no supplied UMAP; use non-UMAP Explore plots. No coordinates were invented.')
    store.update_job(job_id, progress=82, message='Use Existing: supplied annotations and coordinates retained; preparing downstream analyses.')
    frame = pd.DataFrame(coordinates, index=matrix.obs_names, columns=['umap_0', 'umap_1'])
    return ApproximateUMAPResult(query_adata=matrix, coordinates=frame, reference_choices={},
                                 cluster_key=cluster_key, cluster_order=list(dict.fromkeys(matrix.obs[cluster_key].astype(str))))


def reuse_embedded_modalities(matrix, outputs_dir, selected, cluster_key, meta, compression):
    """Expose embedded scALABLE predictions without rerunning their models.

    External sidecar paths cannot establish availability in a new upload. Values
    and their recorded scales/provenance are copied, never reconstructed.
    """
    import anndata as ad
    from ..flask import pipeline as P
    artifacts, markers, warnings = {}, {}, []
    imported = matrix.uns.get('imputed_modalities') or {}
    for modality, info in imported.items():
        if modality in selected or (modality == 'grn_tf' and 'grn' in selected):
            continue
        key, names_key = info.get('obsm_key'), info.get('feature_names_key')
        if not key or key not in matrix.obsm or not names_key or names_key not in matrix.uns:
            warnings.append(f'{modality}: prior predictions are not embedded in the uploaded H5AD; external sidecar paths were not followed.')
            continue
        values, features = matrix.obsm[key], [str(v) for v in matrix.uns[names_key]]
        if values.shape != (matrix.n_obs, len(features)) or not pd.Index(features).is_unique:
            raise ValueError(f'Embedded {modality} prediction identities are inconsistent; no outputs were discarded.')
        if modality not in P._MODALITY_DEFINITIONS:
            warnings.append(f'{modality}: embedded values are retained in the H5AD; this modality has no registered Explore definition.')
            continue
        prediction = ad.AnnData(X=values, obs=matrix.obs.copy(), var=pd.DataFrame(index=features))
        if 'X_umap' in matrix.obsm:
            prediction.obsm['X_umap'] = np.asarray(matrix.obsm['X_umap']).copy()
        summary = dict(info.get('summary') or {})
        prediction.uns['prediction_summary'] = summary
        prediction.uns['modality'] = modality
        for field in ('expression_scale', 'log_base', 'log_pseudocount', 'target_encoding_complete', 'abundance_units', 'BH_policy'):
            if field in summary:
                prediction.uns[field] = summary[field]
        path = outputs_dir / f'existing_{modality}.h5ad'
        P.approx_mod.ensure_h5ad_compat_for_write(prediction)
        prediction.write_h5ad(path, compression=compression)
        artifacts[modality] = {'h5ad': str(path), 'imported_existing': True}
        if P._MODALITY_DEFINITIONS[modality].get('supports_marker_heatmap', True):
            markers[modality] = P._emit_modality_marker_heatmap(modality, prediction, outputs_dir, cluster_key, meta)
    return artifacts, markers, warnings
