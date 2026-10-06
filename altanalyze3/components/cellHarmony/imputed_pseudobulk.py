"""Mean existing predictions per sample/state without panel-total normalization."""
from __future__ import annotations

import copy

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .disk_differential import root_and_rows, read_rows
from .imputed_scale import IMPUTED_MODALITIES, prediction_encoding, inverse_predictions, encode_predictions


def aggregate_imputed_predictions(adata, *, population_col, sample_col, covariate_col, min_cells):
    scale, base, pseudocount = prediction_encoding(adata)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Imputed aggregation requires unique cell and feature identities; do not silently merge them.")
    obs = adata.obs
    grouping = obs[[population_col, sample_col]].astype(str)
    groups = grouping.groupby([population_col, sample_col], sort=False, observed=True).indices
    found = root_and_rows(adata)
    if found:
        root, positions = found
        positions = np.arange(root.n_obs) if isinstance(positions, slice) else positions
    means, records, names, skipped = [], [], [], []
    for (population, sample), rows in groups.items():
        if len(rows) < min_cells:
            skipped.append({"population": population, "sample": sample, "n_cells": len(rows),
                            "reason": "below the established pseudobulk minimum cell count"})
            continue
        labels = obs.iloc[rows][covariate_col].astype(str).unique()
        if len(labels) != 1:
            raise ValueError(f"{population}/{sample}: conflicting comparison assignments; resolve them before aggregation.")
        total = np.zeros(adata.n_vars, dtype=np.float64)
        for start in range(0, len(rows), 128):
            selected = rows[start:start + 128]
            block = read_rows(root.X, positions[selected]) if found else adata.X[selected]
            block = block.toarray() if sparse.issparse(block) else np.asarray(block)
            block = np.asarray(block, dtype=np.float64)
            if not np.isfinite(block).all():
                raise ValueError("Nonfinite predictions; stop instead of filling or dropping features/cells.")
            linear = inverse_predictions(block, (scale, base, pseudocount))
            total += linear.sum(axis=0)
        means.append(total / len(rows))
        records.append({population_col: population, sample_col: sample, "Sample": sample,
                        "Population": population, covariate_col: labels[0], "n_cells": len(rows)})
        names.append(f"{population}|{sample}")
    if not means:
        raise ValueError("No sample/state groups meet the existing pseudobulk minimum cell count.")
    if len(set(names)) != len(names):
        raise ValueError("Ambiguous sample/state group names; do not combine these groups.")
    linear_means = np.asarray(means)
    values = encode_predictions(linear_means, (scale, base, pseudocount))
    uns = copy.deepcopy(dict(adata.uns))
    for key in ("broadcast_profiles", "broadcast_profile_codes"):
        uns.pop(key, None)
    uns.update(pseudobulk_method="pseudobulk", pseudobulk_statistic="arithmetic_mean_of_existing_linear_predictions",
               lipid_panel_total_normalized=False)
    out = ad.AnnData(X=values.astype(np.float32),
                     obs=pd.DataFrame(records, index=pd.Index(names, name="GroupID")),
                     var=adata.var.copy(), uns=uns)
    out.layers["counts"] = linear_means.astype(np.float32)
    if list(out.var_names) != list(adata.var_names):
        raise ValueError("The imputed feature panel changed during aggregation.")
    audit = {"source_cells": adata.n_obs, "aggregated_cells": sum(r["n_cells"] for r in records),
             "groups": len(records), "features": adata.n_vars, "skipped_groups": skipped,
             "statistic": uns["pseudobulk_statistic"], "feature_panel_total_normalized": False,
             "model_recomputed": False, "expression_scale": scale,
             "log_base": base, "log_pseudocount": pseudocount}
    return out, audit
