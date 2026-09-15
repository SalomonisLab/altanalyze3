P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:160]}"
    src = src.replace(old, new)

# --- 1. which numeric obs columns can serve as an axis ------------------------
sub('''def _umap_coordinate_options(cache: Dict[str, Any]) -> List[Dict[str, str]]:''',
    '''def _numeric_obs_columns(cache: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The numeric obs columns a panel may plot on an axis.

    ShinyCell lets a reader put any numeric cell annotation on an axis, so a
    stored `UMAP_1`/`UMAP_2` pair, a pseudotime, a module score or a QC measure
    all work. Counts live in X, not in obs, so nothing here is expression. A
    boolean column is a category, not a measurement, and is left out. So is a
    column with no finite value, which would draw an empty panel.
    """
    adata = cache["adata"]
    out = []
    for name in adata.obs.columns:
        series = adata.obs[name]
        if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
            continue
        values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        if not finite.any():
            continue
        out.append({
            "field": str(name),
            "n_finite": int(finite.sum()),
            "n_missing": int(values.size - finite.sum()),
            "n_unique": int(np.unique(values[finite]).size),
            "min": float(values[finite].min()),
            "max": float(values[finite].max()),
        })
    return out


def _axis_values(cache: Dict[str, Any], field: str) -> Optional[np.ndarray]:
    """One numeric obs column as floats, or None when it cannot serve as an axis."""
    column = str(field or "").strip()
    if not column:
        return None
    adata = cache["adata"]
    if column not in adata.obs.columns:
        return None
    series = adata.obs[column]
    if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
        return None
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    return values if np.isfinite(values).any() else None


def _umap_coordinate_options(cache: Dict[str, Any]) -> List[Dict[str, str]]:''')

# --- 2. obs axes take precedence over an embedding ---------------------------
sub('''    populations, resolved_color_by = _labels_for_color_by(cache_entry, color_by)
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x, umap_y, resolved_coords = _coordinates_for_key(cache_entry, coords_key)''',
    '''    populations, resolved_color_by = _labels_for_color_by(cache_entry, color_by)
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x, umap_y, resolved_coords = _coordinates_for_key(cache_entry, coords_key)
    # A pair of numeric obs columns replaces the embedding outright. Both have to
    # resolve: one axis alone would silently mix a metadata value against a UMAP
    # coordinate, which reads as a map and is not one.
    axis_x = _axis_values(cache_entry, x_field)
    axis_y = _axis_values(cache_entry, y_field)
    resolved_x, resolved_y = "", ""
    if axis_x is not None and axis_y is not None:
        umap_x, umap_y = axis_x, axis_y
        resolved_x, resolved_y = str(x_field).strip(), str(y_field).strip()
        resolved_coords = ""
    axes_source = "obs" if resolved_x else ("obsm" if resolved_coords else "cellharmony")''')

sub('''    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
    color_by: str = "",
    coords_key: str = "",
) -> Dict[str, List[Dict]]:''',
    '''    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
    color_by: str = "",
    coords_key: str = "",
    x_field: str = "",
    y_field: str = "",
) -> Dict[str, List[Dict]]:''')

sub('''    reference_points = []
    ref_adata = None if (resolved_color_by or resolved_coords) else _load_reference_adata(app, meta)''',
    '''    reference_points = []
    ref_adata = (None if (resolved_color_by or resolved_coords or resolved_x)
                 else _load_reference_adata(app, meta))''')

sub('''    return {"reference": reference_points, "query": query_points,
            "sample_field": sample_field,
            "color_by": resolved_color_by,
            "color_label": resolved_color_by or str(cache_entry["cluster_key"]),
            "coords_key": resolved_coords,
            "coords_label": resolved_coords or "cellHarmony UMAP",
            "reference_hidden": bool(resolved_color_by or resolved_coords)}''',
    '''    # How many cells the panel could not place. A metadata axis is often
    # recorded for part of the dataset only, and a silently shorter plot would
    # read as a real absence of cells.
    n_kept = int(np.count_nonzero(np.asarray(display_mask, dtype=bool)))
    return {"reference": reference_points, "query": query_points,
            "sample_field": sample_field,
            "color_by": resolved_color_by,
            "color_label": resolved_color_by or str(cache_entry["cluster_key"]),
            "coords_key": resolved_coords,
            "coords_label": resolved_coords or "cellHarmony UMAP",
            "axes_source": axes_source,
            "x_field": resolved_x,
            "y_field": resolved_y,
            "x_label": resolved_x or "UMAP 1",
            "y_label": resolved_y or "UMAP 2",
            "n_cells_selected": n_kept,
            "n_points_drawn": len(query_points),
            "n_dropped_no_coordinate": max(0, n_kept - len(query_points)),
            "reference_hidden": bool(resolved_color_by or resolved_coords or resolved_x)}''')

# --- 3. the endpoint takes the two fields ------------------------------------
sub('''        color_by: str = Query(""),
        coords: str = Query(""),
    ):''',
    '''        color_by: str = Query(""),
        coords: str = Query(""),
        x_field: str = Query(""),
        y_field: str = Query(""),
    ):''')

sub('''            return JSONResponse(_build_umap_payload(
                app, meta, modality=modality, display_filters=display_filters,
                color_by=color_by, coords_key=coords))''',
    '''            return JSONResponse(_build_umap_payload(
                app, meta, modality=modality, display_filters=display_filters,
                color_by=color_by, coords_key=coords,
                x_field=x_field, y_field=y_field))''')

# --- 4. the lists the panel fills its menus from -----------------------------
sub('''                             # The UMAP panel colours by any of the same columns
                             # and draws on any 2-D embedding the h5ad carries.
                             "color_variables": variables,
                             "coords": _umap_coordinate_options(cache)})''',
    '''                             # The UMAP panel colours by any of the same columns,
                             # draws on any 2-D embedding the h5ad carries, and
                             # takes any pair of numeric obs columns as axes.
                             "color_variables": variables,
                             "coords": _umap_coordinate_options(cache),
                             "numeric_variables": _numeric_obs_columns(cache)})''')

open(P, "w").write(src)
print("patch12 applied")
