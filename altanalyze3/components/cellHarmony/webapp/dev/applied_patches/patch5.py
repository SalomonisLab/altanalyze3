P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:160]}"
    src = src.replace(old, new)

# --- 1. the cache remembers which obsm embeddings the h5ad carries -----------
sub('''            "umap_x": umap_x,
            "umap_y": umap_y,''',
    '''            "umap_x": umap_x,
            "umap_y": umap_y,
            "obsm_keys": _obsm_embedding_keys(adata),''')

# --- 2. helpers: which embeddings, which colour columns, and how to read one -
anchor = "def _build_umap_payload("
helpers = '''def _obsm_embedding_keys(adata) -> List[str]:
    """The obsm entries that can be drawn as a 2-D map, in the h5ad's own order.

    An entry needs at least two columns; the first two are plotted. PCA and
    scVI latent spaces qualify as much as a UMAP does, so nothing is filtered on
    the name - a dataset that stores `X_umap_harmony`, `X_tsne` or `X_scvi` gets
    all of them.
    """
    keys = []
    for key in getattr(adata, "obsm", {}) or {}:
        try:
            matrix = adata.obsm[key]
            if getattr(matrix, "ndim", 0) == 2 and matrix.shape[1] >= 2:
                keys.append(str(key))
        except Exception:  # noqa: BLE001 - an unreadable entry is simply not offered
            continue
    return keys


def _umap_coordinate_options(cache: Dict[str, Any]) -> List[Dict[str, str]]:
    """The coordinate sets a UMAP panel may be drawn on.

    The first entry is the cellHarmony projection, which is what the panel has
    always drawn: the coordinates the alignment wrote, read from the job's
    `umap_coordinates` artifact.
    """
    options = [{"key": "", "label": "cellHarmony UMAP"}]
    for key in cache.get("obsm_keys", []) or []:
        options.append({"key": str(key), "label": str(key)})
    return options


def _coordinates_for_key(cache: Dict[str, Any], coords_key: str = "") -> tuple:
    """(x, y, resolved key). An unknown key falls back to the cellHarmony one."""
    key = str(coords_key or "").strip()
    if not key:
        return cache["umap_x"], cache["umap_y"], ""
    adata = cache["adata"]
    obsm = getattr(adata, "obsm", {}) or {}
    if key not in obsm:
        return cache["umap_x"], cache["umap_y"], ""
    matrix = np.asarray(obsm[key])
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        return cache["umap_x"], cache["umap_y"], ""
    return (np.asarray(matrix[:, 0], dtype=float),
            np.asarray(matrix[:, 1], dtype=float), key)


def _labels_for_color_by(cache: Dict[str, Any], color_by: str = "") -> tuple:
    """(per-cell labels, resolved column). Empty means the cellHarmony states."""
    column = str(color_by or "").strip()
    if not column or column == str(cache["cluster_key"]):
        return cache["populations"], ""
    adata = cache["adata"]
    if column not in adata.obs.columns:
        return cache["populations"], ""
    values = (adata.obs[column].astype(str).str.strip()
              .replace({"nan": "", "None": ""}).to_numpy(dtype=str))
    return values, column


'''
assert src.count(anchor) == 1
src = src.replace(anchor, helpers + anchor)

# --- 3. the payload honours both choices -------------------------------------
sub('''    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict[str, List[Dict]]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    obs_names = cache_entry["obs_names"]
    populations = cache_entry["populations"]
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x = cache_entry["umap_x"]
    umap_y = cache_entry["umap_y"]''',
    '''    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
    color_by: str = "",
    coords_key: str = "",
) -> Dict[str, List[Dict]]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    obs_names = cache_entry["obs_names"]
    populations, resolved_color_by = _labels_for_color_by(cache_entry, color_by)
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x, umap_y, resolved_coords = _coordinates_for_key(cache_entry, coords_key)''')

sub('''    reference_points = []
    ref_adata = _load_reference_adata(app, meta)
    if ref_adata is not None:''',
    '''    # The reference atlas is drawn behind the query only while both choices are
    # the default. Another obs column has no counterpart in the reference, and a
    # second embedding is a different coordinate space, so an overlay there would
    # put the reference cells in positions that mean nothing.
    reference_points = []
    ref_adata = None if (resolved_color_by or resolved_coords) else _load_reference_adata(app, meta)
    if ref_adata is not None:''')

sub('''    return {"reference": reference_points, "query": query_points, "sample_field": sample_field}''',
    '''    return {"reference": reference_points, "query": query_points,
            "sample_field": sample_field,
            "color_by": resolved_color_by,
            "color_label": resolved_color_by or str(cache_entry["cluster_key"]),
            "coords_key": resolved_coords,
            "coords_label": resolved_coords or "cellHarmony UMAP",
            "reference_hidden": bool(resolved_color_by or resolved_coords)}''')

# --- 4. the endpoint takes the two new parameters ----------------------------
sub('''        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            return JSONResponse(_build_umap_payload(app, meta, modality=modality, display_filters=display_filters))''',
    '''        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
        color_by: str = Query(""),
        coords: str = Query(""),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            return JSONResponse(_build_umap_payload(
                app, meta, modality=modality, display_filters=display_filters,
                color_by=color_by, coords_key=coords))''')

# --- 5. plot-variables also answers "what can I colour and plot on?" ---------
sub('''        cache = _get_expression_cache(app, store.get_job(job_id), modality=modality)
        return JSONResponse({"cluster_key": str(cache["cluster_key"]),
                             "variables": _groupable_columns(cache)})''',
    '''        cache = _get_expression_cache(app, store.get_job(job_id), modality=modality)
        variables = _groupable_columns(cache)
        return JSONResponse({"cluster_key": str(cache["cluster_key"]),
                             "variables": variables,
                             # The UMAP panel colours by any of the same columns
                             # and draws on any 2-D embedding the h5ad carries.
                             "color_variables": variables,
                             "coords": _umap_coordinate_options(cache)})''')

open(P, "w").write(src)
print("patch5 applied")
