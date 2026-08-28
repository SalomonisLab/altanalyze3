"""End-to-end check of the shipped human lung bundle.

Loads the bundle through the real inference entry point
(``components.rna2adt.api.load_bundle``), predicts on held-out lung cells
straight from the query .h5ad, and scores the predictions against the measured
TotalVI-denoised ADT. Also drives the cellHarmony-web builder
(``cellHarmony.flask.pipeline._build_imputed_adt_adata``) with the registry
entry so the wiring is exercised, not just the model.

    python -m altanalyze3.components.rna2adt.lung.validate_bundle --n-cells 5000
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from ..api import load_bundle
from . import data as lung_data
from .train_lung import score_per_adt


BUNDLE_DEFAULT = Path(__file__).parent / "rna2adt_hs_lung_bundle.pkl"
REGISTRY_DEFAULT = (Path(__file__).resolve().parents[2]
                    / "cellHarmony" / "flask" / "reference_config.json")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, default=BUNDLE_DEFAULT)
    parser.add_argument("--rna-h5ad", type=Path, default=lung_data.RNA_H5AD_DEFAULT)
    parser.add_argument("--adt-txt", type=Path, default=lung_data.ADT_TXT_DEFAULT)
    parser.add_argument("--registry", type=Path, default=REGISTRY_DEFAULT)
    parser.add_argument("--reference-id", default="hs_lung_hlca_reference")
    parser.add_argument("--n-cells", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--out-tsv", type=Path, default=None)
    args = parser.parse_args()

    bundle = load_bundle(args.bundle)
    info = bundle.model_info()
    print("[bundle] " + json.dumps({k: v for k, v in info.items() if k != "metadata"}, indent=1))
    held = list((bundle.metadata or {}).get("held_out_donors") or [])
    print(f"[bundle] held_out_donors declared: {held or 'none'}")

    adt = lung_data.load_adt(args.adt_txt)
    index = lung_data.align_cells(args.rna_h5ad, adt)
    libraries = lung_data.rna_obs_column(args.rna_h5ad, "Library")[index.rna_rows]
    donors = np.array([str(v).split("_", 1)[0] for v in libraries], dtype=object)

    eligible = np.where(np.isin(donors, held))[0] if held else np.arange(index.rna_rows.size)
    scope = ("held-out donors only" if held else
             "all matched cells - NOT a clean holdout, the shipped bundle used a random "
             "cell split so ~13% of these cells were in its training set; the clean number "
             "is the bundle's own cell_holdout metric")
    rng = np.random.default_rng(args.seed)
    pick = np.sort(rng.choice(eligible, min(args.n_cells, eligible.size), replace=False))
    print(f"[scope] scoring on {pick.size} cells from {scope}; "
          f"{len(set(donors[pick]))} donors")

    rna_rows = index.rna_rows[pick]
    genes = [str(g) for g in lung_data.rna_var_names(args.rna_h5ad)]
    matrix = lung_data.read_rna_rows(args.rna_h5ad, rna_rows, verbose=False)
    query = ad.AnnData(X=matrix,
                       obs=pd.DataFrame(index=[str(n) for n in index.obs_names[pick]]),
                       var=pd.DataFrame(index=genes))
    print(f"[query] AnnData {query.shape} built straight from the .h5ad X layer")

    result = bundle.predict_from_adata(query)
    print("[predict] summary: " + json.dumps(result.summary, indent=1, default=str))
    predicted = result.predictions
    measured = pd.DataFrame(adt.values[index.adt_rows[pick]],
                            index=query.obs_names,
                            columns=[str(v) for v in adt.var_names])
    shared = [c for c in predicted.columns if c in measured.columns]
    if len(shared) != predicted.shape[1]:
        print(f"FAIL {predicted.shape[1] - len(shared)} predicted ADTs are absent from the "
              f"measured table")
        return 1
    pearsons, spearmans = score_per_adt(measured[shared].to_numpy(dtype=np.float32),
                                        predicted[shared].to_numpy(dtype=np.float32))
    frame = pd.DataFrame({"adt_raw": shared,
                          "pearson": pearsons,
                          "spearman": spearmans}).sort_values("pearson", ascending=False)
    print(f"\n=== bundle vs measured ADT ({pick.size} cells, {len(shared)} ADTs, {scope}) ===")
    print(f"mean Pearson:  {np.nanmean(pearsons):.3f}   median: {np.nanmedian(pearsons):.3f}")
    print(f"mean Spearman: {np.nanmean(spearmans):.3f}  median: {np.nanmedian(spearmans):.3f}")
    print(f"valid ADTs: {int(np.isfinite(pearsons).sum())} / {len(pearsons)}")
    print(frame.to_string(index=False))
    if args.out_tsv:
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.out_tsv, sep="\t", index=False)
        print(f"  wrote {args.out_tsv}")

    # --- cellHarmony-web wiring ---
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from altanalyze3.components.cellHarmony.flask.pipeline import (  # noqa: E402
        _build_imputed_adt_adata, _lookup_reference)
    entry = _lookup_reference("human", args.reference_id, args.registry)
    config = (entry.get("impute_config") or {}).get("adt")
    if not config:
        print(f"FAIL registry entry {args.reference_id} carries no impute_config.adt")
        return 1
    print(f"\n[registry] {args.reference_id} impute_modalities={entry.get('impute_modalities')}")
    print(f"[registry] adt config={json.dumps(config)}")
    if not Path(str(config["bundle_path"])).exists():
        print(f"FAIL bundle_path does not resolve: {config['bundle_path']}")
        return 1
    imputed, summary = _build_imputed_adt_adata(query, entry)
    print(f"[pipeline] imputed AnnData {imputed.shape} "
          f"expression_scale={imputed.uns.get('expression_scale')} "
          f"log_base={summary.get('log_base')} "
          f"clipped_negatives={summary.get('clipped_negative_values')}")
    print(f"[pipeline] var_names sample: {list(imputed.var_names[:6])}")
    linear = imputed.layers["counts"]
    print(f"[pipeline] counts layer (linear) min={linear.min():.3f} max={linear.max():.1f}")
    if imputed.shape[0] != query.n_obs or imputed.shape[1] != len(shared):
        print("FAIL imputed AnnData shape does not match the query")
        return 1
    print("\nPASS bundle loads through api.load_bundle, predicts from a raw query .h5ad, "
          "and the cellHarmony-web ADT builder runs on the registry entry")
    return 0


if __name__ == "__main__":
    sys.exit(main())
