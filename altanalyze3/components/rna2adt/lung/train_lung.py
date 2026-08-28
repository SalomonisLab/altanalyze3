"""Train the human lung rna2adt bundle: per-protein ElasticNet on whitelist union.

Same architecture and hyper-parameters as the bone marrow production bundle
(``Rna2LipidArchPanelPerProtein`` with ``feature_source="whitelist_union"``,
alpha 0.01, l1_ratio 0.5) and the mouse port (``MouseRna2AdtModel``):

* input  - per-cell RNA on the panel-union of whitelist genes, z-scored
* output - one ElasticNet head per ADT, inverse-z-scored back to the
           log1p TotalVI-denoised scale the ADT file carries
* heads are fit independently; L1 inside each head selects its own features

Two evaluation splits are reported:

* ``cells``  - random cell holdout, the split the mouse and bone marrow
               READMEs report. This is the default and the split the shipped
               bundle is trained on.
* ``donor``  - whole donors held out. Cells from one donor share ambient RNA,
               antibody staining batch and clonal composition, so a random
               cell holdout leaks; the donor split is the honest estimate of
               what the bundle does on a new dataset.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import data as lung_data
from .adt_rna_map import strip_prefix
from .model import HumanLungRna2AdtModel


def load_whitelist(path: Path) -> Dict[str, List[str]]:
    frame = pd.read_csv(path, sep="\t")
    return {str(row["adt_raw"]): [g.strip() for g in str(row["feature_genes"]).split(",") if g.strip()]
            for _, row in frame.iterrows()}


def panel_union(whitelist: Dict[str, List[str]], adt_names: Sequence[str],
                atlas_genes: Sequence[str]) -> Tuple[List[str], int]:
    present = set(atlas_genes)
    union: List[str] = []
    seen: set = set()
    n_no_entry = 0
    for adt in adt_names:
        features = whitelist.get(str(adt)) or []
        if not features:
            n_no_entry += 1
        for gene in features:
            if gene in present and gene not in seen:
                seen.add(gene)
                union.append(gene)
    return union, n_no_entry


def score_per_adt(Y_true: np.ndarray, Y_pred: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    pearsons, spearmans = [], []
    for j in range(Y_true.shape[1]):
        yt, yp = Y_true[:, j], Y_pred[:, j]
        if yt.std() < 1e-9 or yp.std() < 1e-9:
            pearsons.append(np.nan)
            spearmans.append(np.nan)
            continue
        pearsons.append(float(np.corrcoef(yt, yp)[0, 1]))
        try:
            spearmans.append(float(spearmanr(yt, yp).correlation))
        except Exception:
            spearmans.append(np.nan)
    return np.asarray(pearsons), np.asarray(spearmans)


def donor_of(library: str) -> str:
    return str(library).split("_", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna-h5ad", type=Path, default=lung_data.RNA_H5AD_DEFAULT)
    parser.add_argument("--adt-txt", type=Path, default=lung_data.ADT_TXT_DEFAULT)
    parser.add_argument("--whitelist", type=Path,
                        default=Path(__file__).parent / "configs" / "empirical_whitelist.tsv")
    parser.add_argument("--bundle-out", type=Path,
                        default=Path(__file__).parent / "rna2adt_hs_lung_bundle.pkl")
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--params-out", type=Path, default=None)
    parser.add_argument("--holdout", choices=("cells", "donor"), default="cells")
    parser.add_argument("--n-holdout-donors", type=int, default=10)
    parser.add_argument("--library-col", default="Library")
    parser.add_argument("--max-train-cells", type=int, default=40000)
    parser.add_argument("--max-test-cells", type=int, default=15000)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--l1-ratio", type=float, default=0.5)
    parser.add_argument("--max-iter", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-save-bundle", action="store_true")
    parser.add_argument("--panel-cache", type=Path, default=None,
                        help="npz holding the extracted (cells x panel genes) matrix. Written "
                             "on a miss, reused only when the gene list and the cell list match "
                             "exactly.")
    args = parser.parse_args()

    started = time.time()

    def step(message: str) -> None:
        print(f"[{time.time() - started:7.1f}s] {message}", flush=True)

    step(f"loading ADT {args.adt_txt}")
    adt = lung_data.load_adt(args.adt_txt)
    adt_names = [str(v) for v in adt.var_names]
    step(f"{adt.obs_names.size} ADT rows, {len(adt_names)} ADTs")

    step(f"aligning against {args.rna_h5ad}")
    index = lung_data.align_cells(args.rna_h5ad, adt)
    alignment = lung_data.report_alignment(index)
    step(f"matched {alignment['n_matched_cells']} cells "
         f"(RNA retained {100 * alignment['rna_retained_fraction']:.2f}%, "
         f"ADT retained {100 * alignment['adt_retained_fraction']:.2f}%)")

    atlas_genes = [str(g) for g in lung_data.rna_var_names(args.rna_h5ad)]
    whitelist = load_whitelist(args.whitelist)
    union, n_no_entry = panel_union(whitelist, adt_names, atlas_genes)
    step(f"{len(union)} panel-union genes (ADTs with no whitelist entry: {n_no_entry})")
    if not union:
        raise SystemExit("panel union is empty; check the whitelist path")

    gene_position = {g: i for i, g in enumerate(atlas_genes)}
    gene_indices = np.asarray([gene_position[g] for g in union], dtype=np.int64)

    X_all = None
    if args.panel_cache and args.panel_cache.exists():
        cached = np.load(args.panel_cache, allow_pickle=True)
        cached_genes = [str(g) for g in cached["genes"]]
        cached_cells = [str(c) for c in cached["obs_names"]]
        wanted_cells = [str(c) for c in index.obs_names]
        if cached_genes == list(union) and cached_cells == wanted_cells:
            X_all = cached["X"]
            step(f"reused panel cache {args.panel_cache} {X_all.shape}")
        else:
            step(f"panel cache {args.panel_cache} does not match this gene set or cell set; "
                 f"re-streaming")
    if X_all is None:
        step(f"streaming RNA for {index.rna_rows.size} cells x {len(union)} panel genes")
        X_all = lung_data.read_rna_rows(args.rna_h5ad, index.rna_rows,
                                        gene_indices=gene_indices, dense=True, verbose=False)
        if args.panel_cache:
            args.panel_cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez(args.panel_cache, X=X_all,
                     genes=np.array(union, dtype=object),
                     obs_names=np.array([str(c) for c in index.obs_names], dtype=object))
            step(f"wrote panel cache {args.panel_cache}")
    Y_all = adt.values[index.adt_rows]
    if X_all.shape[0] != Y_all.shape[0]:
        raise ValueError(f"cell mismatch: X {X_all.shape[0]} vs Y {Y_all.shape[0]}")
    step(f"X {X_all.shape} Y {Y_all.shape}")

    libraries = lung_data.rna_obs_column(args.rna_h5ad, args.library_col)[index.rna_rows]
    donors = np.array([donor_of(v) for v in libraries], dtype=object)

    rng = np.random.default_rng(args.seed)
    n_cells = X_all.shape[0]
    if args.holdout == "cells":
        # Draw the split from a PERMUTATION, not from a sorted sample. Sorting the
        # sample and then slicing head/tail makes the test set the highest-indexed
        # cells, which in this atlas means the last libraries in file order - a
        # positional holdout wearing a random holdout's label. The mouse trainer
        # (mouse/train_mouse.py) still slices a sorted sample; that difference is
        # deliberate and is recorded in this package's README.
        n_keep = min(n_cells, args.max_train_cells + args.max_test_cells)
        selected = rng.permutation(n_cells)[:n_keep]
        n_train = min(args.max_train_cells, n_keep - 1)
        train_idx = np.sort(selected[:n_train])
        test_idx = np.sort(selected[n_train:])
        split_label = f"cell_holdout[{len(train_idx)}train/{len(test_idx)}test]"
        held_out_donors: List[str] = []
    else:
        unique_donors = np.array(sorted(set(donors)))
        held = rng.choice(unique_donors, min(args.n_holdout_donors, unique_donors.size),
                          replace=False)
        held_out_donors = sorted(str(d) for d in held)
        is_held = np.isin(donors, held)
        train_pool = np.where(~is_held)[0]
        test_pool = np.where(is_held)[0]
        if train_pool.size > args.max_train_cells:
            train_idx = np.sort(rng.choice(train_pool, args.max_train_cells, replace=False))
        else:
            train_idx = train_pool
        if test_pool.size > args.max_test_cells:
            test_idx = np.sort(rng.choice(test_pool, args.max_test_cells, replace=False))
        else:
            test_idx = test_pool
        split_label = (f"donor_holdout[{len(held_out_donors)} donors: "
                       f"{','.join(held_out_donors)}]")
    step(f"split {split_label}: train={train_idx.size} test={test_idx.size}")
    if set(donors[train_idx]) & set(donors[test_idx]) and args.holdout == "donor":
        raise ValueError("donor leak between train and test")

    X_train, Y_train = X_all[train_idx], Y_all[train_idx]
    X_test, Y_test = X_all[test_idx], Y_all[test_idx]

    step("fitting model")
    model = HumanLungRna2AdtModel(alpha=args.alpha, l1_ratio=args.l1_ratio,
                                  max_iter=args.max_iter)
    model.fit(X_train, Y_train, rna_genes=union, adt_names=adt_names)

    step(f"evaluating on {split_label}")
    predicted = model.predict(X_test)
    pearsons, spearmans = score_per_adt(Y_test, predicted)
    print(f"\n=== {split_label} (n_test={Y_test.shape[0]}, n_adts={Y_test.shape[1]}) ===")
    print(f"mean Pearson:  {np.nanmean(pearsons):.3f}   median: {np.nanmedian(pearsons):.3f}")
    print(f"mean Spearman: {np.nanmean(spearmans):.3f}  median: {np.nanmedian(spearmans):.3f}")
    print(f"valid ADTs: {int(np.isfinite(pearsons).sum())} / {len(pearsons)}")

    metrics_path = args.metrics_out or (
        args.bundle_out.parent / f"rna2adt_hs_lung_per_adt_metrics_{args.holdout}.tsv")
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "adt_raw": adt_names,
        "adt_clean": [strip_prefix(n) for n in adt_names],
        "pearson": pearsons,
        "spearman": spearmans,
        "split": split_label,
    }).to_csv(metrics_path, sep="\t", index=False)
    print(f"  per-ADT metrics: {metrics_path}")

    metadata = {
        "approach": HumanLungRna2AdtModel.head_kind,
        "feature_source": "whitelist_union",
        "label": "rna2adt_hs_lung_bundle",
        "species": "human",
        "tissue": "lung",
        "reference": "COVID-TotalVI human lung CITE-seq",
        "adt_name_format": "Hu.<marker>, catalogue suffix stripped",
        "target_scaling": {"mode": "none"},
        "target_transform": "log1p(TotalVI-denoised ADT), natural log",
        "expression_scale": "log1p",
        "log_base": "e",
        "rna_input_scale": "CP10k+log1p (natural log), the atlas X layer",
        "alpha": args.alpha,
        "l1_ratio": args.l1_ratio,
        "max_iter": args.max_iter,
        "seed": args.seed,
        "split": split_label,
        "held_out_donors": held_out_donors,
        "n_train_cells": int(X_train.shape[0]),
        "n_test_cells": int(X_test.shape[0]),
        "n_panel_union": int(len(union)),
        "n_adts": int(len(adt_names)),
        "alignment": alignment,
        "rna_h5ad": str(args.rna_h5ad),
        "adt_txt": str(args.adt_txt),
        "whitelist": str(args.whitelist),
        "evaluation_summary": {
            args.holdout: {
                "mean_pearson": float(np.nanmean(pearsons)),
                "median_pearson": float(np.nanmedian(pearsons)),
                "mean_spearman": float(np.nanmean(spearmans)),
                "median_spearman": float(np.nanmedian(spearmans)),
                "n_valid_adts": int(np.isfinite(pearsons).sum()),
            }
        },
    }

    if args.params_out:
        args.params_out.parent.mkdir(parents=True, exist_ok=True)
        with args.params_out.open("w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, default=str)
        print(f"  parameters: {args.params_out}")

    if not args.no_save_bundle:
        step("saving bundle")
        bundle = {
            "model": model,
            "scaler_x": None,      # handled inside the model
            "scaler_y": None,
            "X_columns": union,
            "Y_columns": adt_names,
            "metadata": metadata,
        }
        args.bundle_out.parent.mkdir(parents=True, exist_ok=True)
        with args.bundle_out.open("wb") as handle:
            pickle.dump(bundle, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"  saved {args.bundle_out}")

    step("done")


if __name__ == "__main__":
    main()
