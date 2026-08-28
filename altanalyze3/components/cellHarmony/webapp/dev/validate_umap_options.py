"""Prove the coordinate and colour switches read obsm/obs, not the defaults.

The job h5ad on this machine stores one embedding whose values equal the
cellHarmony TSV, so the live payload cannot tell the two sources apart. This
builds an AnnData with two DIFFERENT embeddings and three annotation columns and
calls the helpers directly.
"""
import sys
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
import numpy as np, pandas as pd, anndata as ad
import importlib
# webapp/__init__.py rebinds the name `app` to the FastAPI instance, so the
# module itself has to be pulled out of sys.modules.
A = importlib.import_module("altanalyze3.components.cellHarmony.webapp.app")

n = 40
rng = np.random.default_rng(0)
obs = pd.DataFrame({
    "cellHarmony_state": pd.Categorical(["HSC"] * 20 + ["GMP"] * 20),
    "author_atlas":      pd.Categorical(["A"] * 10 + ["B"] * 10 + ["C"] * 20),
}, index=[f"cell{i}" for i in range(n)])
adata = ad.AnnData(X=rng.random((n, 5), dtype=np.float32), obs=obs)
adata.obsm["X_umap"] = np.column_stack([np.arange(n), np.arange(n) * 2.0])
adata.obsm["X_umap_harmony"] = np.column_stack([np.arange(n) * -1.0, np.arange(n) * -3.0])
adata.obsm["X_pca"] = rng.random((n, 10))
adata.obsm["too_narrow"] = np.arange(n).reshape(n, 1)

cache = {
    "adata": adata,
    "cluster_key": "cellHarmony_state",
    "populations": adata.obs["cellHarmony_state"].astype(str).to_numpy(),
    "umap_x": np.full(n, 99.0),          # the cellHarmony TSV coordinates
    "umap_y": np.full(n, -99.0),
    "obsm_keys": A._obsm_embedding_keys(adata),
    "var_names": adata.var_names.astype(str).to_numpy(),
}

print("obsm keys offered:", cache["obsm_keys"])
assert cache["obsm_keys"] == ["X_umap", "X_umap_harmony", "X_pca"], cache["obsm_keys"]
print("  'too_narrow' (1 column) excluded: yes")

print("\ncoordinate options:", A._umap_coordinate_options(cache))

for key, want_x0, want_y0, want_key in [
    ("",                99.0,  -99.0, ""),
    ("X_umap",           0.0,    0.0, "X_umap"),
    ("X_umap_harmony",   0.0,   -0.0, "X_umap_harmony"),
    ("X_pca",  float(adata.obsm["X_pca"][0, 0]), float(adata.obsm["X_pca"][0, 1]), "X_pca"),
    ("missing_key",     99.0,  -99.0, ""),
    ("too_narrow",      99.0,  -99.0, ""),
]:
    x, y, resolved = A._coordinates_for_key(cache, key)
    ok = (abs(float(x[0]) - want_x0) < 1e-9 and abs(float(y[0]) - want_y0) < 1e-9
          and resolved == want_key)
    print(f"  coords={key!r:18} -> first cell ({float(x[0]):.3f}, {float(y[0]):.3f}) "
          f"resolved={resolved!r:16} {'OK' if ok else 'WRONG'}")
    assert ok
    if key == "X_umap_harmony":
        assert float(x[5]) == -5.0 and float(y[5]) == -15.0, (x[5], y[5])

print("\ncolour column:")
for column, want_first, want_key in [
    ("",                  "HSC", ""),
    ("cellHarmony_state", "HSC", ""),          # the default, reported as default
    ("author_atlas",      "A",   "author_atlas"),
    ("no_such_column",    "HSC", ""),
]:
    labels, resolved = A._labels_for_color_by(cache, column)
    ok = str(labels[0]) == want_first and resolved == want_key
    print(f"  color_by={column!r:20} -> first label {str(labels[0])!r:6} "
          f"resolved={resolved!r:14} n_levels={len(set(map(str, labels)))} {'OK' if ok else 'WRONG'}")
    assert ok

print("\nall coordinate and colour cases pass")


# ---------------------------------------------------------------------------
# Numeric obs columns as axes, and the whole payload around them.
#
# No h5ad on this machine holds two numeric obs columns or a second embedding,
# so the payload is exercised against a cache built here. `_get_expression_cache`
# and `_load_reference_adata` are replaced so the builder reads this cache and a
# small stand-in reference instead of touching disk.
# ---------------------------------------------------------------------------
print("\n=== numeric obs columns as axes ===")

obs2 = pd.DataFrame({
    "cellHarmony_state": pd.Categorical(["HSC"] * 20 + ["GMP"] * 20),
    "author_atlas":      pd.Categorical(["A"] * 10 + ["B"] * 10 + ["C"] * 20),
    "pseudotime":        np.linspace(0.0, 1.0, n),
    "percent_mt":        np.linspace(10.0, 0.0, n),
    "flagged":           np.array([True, False] * (n // 2)),      # boolean: not a measurement
    "all_missing":       np.full(n, np.nan),                      # nothing to draw
}, index=[f"cell{i}" for i in range(n)])
# One axis recorded for half the cells, which is the case the caption must report.
obs2["late_marker_score"] = np.where(np.arange(n) < n // 2, np.nan, np.arange(n, dtype=float))
adata2 = ad.AnnData(X=rng.random((n, 5), dtype=np.float32), obs=obs2)
adata2.obsm["X_umap"] = np.column_stack([np.arange(n), np.arange(n) * 2.0])

cache2 = {
    "adata": adata2,
    "cluster_key": "cellHarmony_state",
    "populations": adata2.obs["cellHarmony_state"].astype(str).to_numpy(),
    "obs_names": adata2.obs_names.astype(str).to_numpy(),
    "umap_x": np.full(n, 99.0),
    "umap_y": np.full(n, -99.0),
    "obsm_keys": A._obsm_embedding_keys(adata2),
    "var_names": adata2.var_names.astype(str).to_numpy(),
    "sample_field": "",
    "sample_labels": None,
}

offered = A._numeric_obs_columns(cache2)
names = [entry["field"] for entry in offered]
print("offered as axes:", names)
assert names == ["pseudotime", "percent_mt", "late_marker_score"], names
print("  boolean 'flagged' excluded: yes | all-NaN 'all_missing' excluded: yes")
print("  categorical columns excluded: yes")
late = next(e for e in offered if e["field"] == "late_marker_score")
print(f"  late_marker_score: {late['n_finite']} finite, {late['n_missing']} missing, "
      f"range {late['min']:.1f} .. {late['max']:.1f}")
assert late["n_finite"] == n // 2 and late["n_missing"] == n // 2

# A stand-in reference, so "is the reference drawn?" is a real test.
reference = ad.AnnData(X=rng.random((6, 5), dtype=np.float32),
                       obs=pd.DataFrame({"state": ["HSC"] * 6},
                                        index=[f"ref{i}" for i in range(6)]))
reference.obsm["X_umap"] = np.column_stack([np.arange(6, dtype=float), np.arange(6, dtype=float)])

A._get_expression_cache = lambda app, meta, modality="rna": cache2
A._load_reference_adata = lambda app, meta: reference
meta = {"job_id": "synthetic", "cluster_key": "cellHarmony_state", "reference_cluster_key": "state"}

def payload(**kwargs):
    return A._build_umap_payload(None, meta, modality="rna", display_filters=None, **kwargs)

print("\n  case                                 axes_source  ref pts  drawn  dropped  x[0]     y[0]")
cases = [
    ("default",                         {}),
    ("coords=X_umap",                   {"coords_key": "X_umap"}),
    ("x/y = pseudotime, percent_mt",     {"x_field": "pseudotime", "y_field": "percent_mt"}),
    ("y half missing",                   {"x_field": "pseudotime", "y_field": "late_marker_score"}),
    ("one axis only (must fall back)",   {"x_field": "pseudotime"}),
    ("boolean axis (must fall back)",    {"x_field": "flagged", "y_field": "percent_mt"}),
]
for label, kwargs in cases:
    out = payload(**kwargs)
    first = out["query"][0]
    print(f"  {label:36} {out['axes_source']:11}  {len(out['reference']):7}  "
          f"{out['n_points_drawn']:5}  {out['n_dropped_no_coordinate']:7}  "
          f"{first['x']:7.2f}  {first['y']:7.2f}")

# The assertions behind that table.
default = payload()
assert default["axes_source"] == "cellharmony" and len(default["reference"]) == 6
assert default["query"][0]["x"] == 99.0 and default["n_dropped_no_coordinate"] == 0

embedded = payload(coords_key="X_umap")
assert embedded["axes_source"] == "obsm" and embedded["reference"] == []
assert embedded["query"][0]["x"] == 0.0 and embedded["query"][5]["y"] == 10.0

axes = payload(x_field="pseudotime", y_field="percent_mt")
assert axes["axes_source"] == "obs" and axes["reference"] == []
assert axes["x_label"] == "pseudotime" and axes["y_label"] == "percent_mt"
assert abs(axes["query"][0]["x"] - 0.0) < 1e-9 and abs(axes["query"][0]["y"] - 10.0) < 1e-9
assert axes["n_points_drawn"] == n and axes["n_dropped_no_coordinate"] == 0

partial = payload(x_field="pseudotime", y_field="late_marker_score")
assert partial["n_points_drawn"] == n // 2, partial["n_points_drawn"]
assert partial["n_dropped_no_coordinate"] == n // 2
assert partial["n_cells_selected"] == n

for bad in ({"x_field": "pseudotime"}, {"x_field": "flagged", "y_field": "percent_mt"},
            {"x_field": "no_such_column", "y_field": "percent_mt"}):
    out = payload(**bad)
    assert out["axes_source"] == "cellharmony" and out["query"][0]["x"] == 99.0
    assert len(out["reference"]) == 6

print("\nall axis cases pass, including the half-missing column and every fallback")
