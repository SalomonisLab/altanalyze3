P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

old = '''    if cache_path and cache_path.exists():
        with np.load(cache_path, allow_pickle=False) as npz_data:
            matrix = np.asarray(npz_data["matrix"], dtype=np.float32)
            row_ids = np.asarray(npz_data["row_ids"], dtype=str)
            col_ids = np.asarray(npz_data["col_ids"], dtype=str)
            if "col_barcodes" in npz_data:
                col_barcodes = np.asarray(npz_data["col_barcodes"], dtype=str)
            else:
                col_barcodes = np.asarray(
                    [value.split(":", 1)[1] if ":" in value else value for value in col_ids],
                    dtype=str,
                )
        entry = {'''

new = '''    if cache_path and cache_path.exists():
        # The fold-matrix cache moved from .npz to .h5ad, at
        # visualization/marker_heatmap_h5ad.py:73. np.load cannot open an HDF5
        # file, so every MarkerHeatmap request raised
        # "Cannot load file containing pickled data when allow_pickle=False".
        # `_read_heatmap_cache` is the writer's own reader and it opens both
        # formats, so this calls it rather than parsing the cache a second way.
        # The import sits here because that module loads scanpy, which would add
        # seconds to app start-up for a view many jobs never open.
        from altanalyze3.components.visualization.marker_heatmap_h5ad import _read_heatmap_cache

        heatmap_df, row_clusters, column_clusters, _, _ = _read_heatmap_cache(cache_path)
        matrix = np.asarray(heatmap_df.to_numpy(), dtype=np.float32)
        col_barcodes = heatmap_df.columns.astype(str).to_numpy()
        genes = heatmap_df.index.astype(str).to_numpy()
        # The ids the .npz held were "cluster:gene" and "cluster:barcode". They
        # are rebuilt here so every reader downstream sees what it saw before.
        row_ids = np.asarray([f"{c}:{g}" for c, g in zip(row_clusters, genes)], dtype=str)
        col_ids = np.asarray([f"{c}:{b}" for c, b in zip(column_clusters, col_barcodes)], dtype=str)
        entry = {'''

assert src.count(old) == 1, src.count(old)
open(P, "w").write(src.replace(old, new))
print("patch11 applied")
