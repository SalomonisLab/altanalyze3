"""Does the lung model do better on linear or on log1p ADT targets?

The bone marrow bundle trains on log1p ADT values and the cellHarmony-web
registry declares ``expression_scale: log1p`` for it, so the lung bundle
follows. This script measures whether that choice costs anything on the lung
panel: same architecture, same split, same features, targets either
``log1p(denoised)`` or the raw linear denoised values. Correlation metrics are
computed on each model's own target scale and, for a like-for-like comparison,
also after mapping both predictions back to the linear scale.

Reads the panel matrix from the cache ``train_lung --panel-cache`` writes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import data as lung_data
from .model import HumanLungRna2AdtModel
from .train_lung import donor_of, score_per_adt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-cache", type=Path, required=True)
    parser.add_argument("--rna-h5ad", type=Path, default=lung_data.RNA_H5AD_DEFAULT)
    parser.add_argument("--adt-txt", type=Path, default=lung_data.ADT_TXT_DEFAULT)
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--max-train-cells", type=int, default=40000)
    parser.add_argument("--max-test-cells", type=int, default=15000)
    parser.add_argument("--extra-train-cells", type=int, default=150000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    cached = np.load(args.panel_cache, allow_pickle=True)
    X_all = cached["X"]
    union = [str(g) for g in cached["genes"]]
    cells = [str(c) for c in cached["obs_names"]]
    print(f"[cache] X {X_all.shape} over {len(union)} genes", flush=True)

    adt = lung_data.load_adt(args.adt_txt, log1p=False)          # linear
    index = lung_data.align_cells(args.rna_h5ad, adt)
    if [str(c) for c in index.obs_names] != cells:
        raise SystemExit("cache cell order does not match the current alignment")
    Y_linear = adt.values[index.adt_rows]
    Y_log1p = np.log1p(Y_linear.astype(np.float64)).astype(np.float32)
    adt_names = [str(v) for v in adt.var_names]
    libraries = lung_data.rna_obs_column(args.rna_h5ad, "Library")[index.rna_rows]
    donors = np.array([donor_of(v) for v in libraries], dtype=object)

    n_cells = X_all.shape[0]

    def cell_split(rng, n_train: int):
        # permutation, not a sorted sample - see train_lung for why
        n_keep = min(n_cells, n_train + args.max_test_cells)
        selected = rng.permutation(n_cells)[:n_keep]
        cut = min(n_train, n_keep - 1)
        return np.sort(selected[:cut]), np.sort(selected[cut:])

    def donor_split(rng, n_train: int):
        unique = np.array(sorted(set(donors)))
        held = rng.choice(unique, 10, replace=False)
        is_held = np.isin(donors, held)
        train_pool, test_pool = np.where(~is_held)[0], np.where(is_held)[0]
        train = (np.sort(rng.choice(train_pool, n_train, replace=False))
                 if train_pool.size > n_train else train_pool)
        test = (np.sort(rng.choice(test_pool, args.max_test_cells, replace=False))
                if test_pool.size > args.max_test_cells else test_pool)
        return train, test

    runs = [
        ("log1p", "cells", args.max_train_cells),
        ("linear", "cells", args.max_train_cells),
        ("log1p", "donor", args.max_train_cells),
        ("linear", "donor", args.max_train_cells),
        ("log1p", "cells", args.extra_train_cells),
        ("log1p", "donor", args.extra_train_cells),
    ]

    rows = []
    for target, split_kind, n_train in runs:
        rng = np.random.default_rng(args.seed)   # identical split for every config
        train_idx, test_idx = (cell_split(rng, n_train) if split_kind == "cells"
                               else donor_split(rng, n_train))
        Y = Y_log1p if target == "log1p" else Y_linear
        model = HumanLungRna2AdtModel()
        model.fit(X_all[train_idx], Y[train_idx], rna_genes=union, adt_names=adt_names)
        predicted = model.predict(X_all[test_idx])
        own = score_per_adt(Y[test_idx], predicted)[0]
        # map both back to linear so the two targets are comparable
        pred_linear = predicted if target == "linear" else np.expm1(predicted.astype(np.float64))
        linear = score_per_adt(Y_linear[test_idx], np.asarray(pred_linear, dtype=np.float32))[0]
        label = f"{target}|{split_kind}|{len(train_idx)}train"
        print(f"[{label}] test={len(test_idx)} "
              f"mean_pearson_own_scale={np.nanmean(own):.3f} "
              f"median={np.nanmedian(own):.3f} "
              f"mean_pearson_linear_scale={np.nanmean(linear):.3f}", flush=True)
        for j, name in enumerate(adt_names):
            rows.append({"config": label, "target_scale": target, "split": split_kind,
                         "n_train": len(train_idx), "n_test": len(test_idx),
                         "adt_raw": name,
                         "pearson_own_scale": own[j],
                         "pearson_linear_scale": linear[j]})

    frame = pd.DataFrame(rows)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"\nwrote {args.out_tsv}")
    summary = frame.groupby("config").agg(
        mean_own=("pearson_own_scale", "mean"), median_own=("pearson_own_scale", "median"),
        mean_linear=("pearson_linear_scale", "mean"),
        median_linear=("pearson_linear_scale", "median"))
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
