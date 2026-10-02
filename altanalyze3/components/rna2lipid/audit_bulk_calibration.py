"""Audit preprocessing and construct explicitly provisional bulk-referenced proxies.

Run with /usr/bin/python3 (pandas/numpy). A proxy is a per-lipid multiplicative
calibration in linear space. It preserves folds; it does not establish physical
units, validate cross-study assay comparability, or validate an imputation model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def normalize_mode(name):
    # Preserve structural annotations, isomer/adduct labels and ionization mode.
    name = str(name).strip()
    if name.endswith("_POS"):
        return name[:-4] + "_P"
    if name.endswith("_NEG"):
        return name[:-4] + "_N"
    return name


def calibrate(native, bulk):
    if not native.columns.is_unique or not bulk.columns.is_unique:
        raise ValueError("Duplicate feature IDs require explicit resolution.")
    shared = native.columns.intersection(bulk.columns)
    pmx = native.loc[native.index.str.endswith("_PMX"), shared]
    if pmx.empty or len(shared) == 0:
        raise ValueError("No shared features or PMX references.")
    reference = pmx.median()
    baseline = bulk.loc[:, shared].median()
    valid = reference.notna() & baseline.notna()
    shared = shared[valid]
    offsets = baseline.loc[shared] - reference.loc[shared]
    anchored = native.loc[:, shared].add(offsets, axis=1)
    recovered = anchored.subtract(offsets, axis=1)
    error = float(np.nanmax(abs(recovered.values - native.loc[:, shared].values)))
    if error > 1e-10:
        raise ValueError("Calibration failed fold-preserving round trip.")
    parameters = pd.DataFrame({
        "bulk_median_log2": baseline.loc[shared],
        "pmx_median_native_log2": reference.loc[shared],
        "log2_offset": offsets,
        "bulk_n": bulk.loc[:, shared].notna().sum(),
        "pmx_n": pmx.loc[:, shared].notna().sum(),
    })
    return anchored, parameters, error


def audit(source, native_path, output_dir):
    source, output = Path(source), Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    files = {
        "bulk_raw": source / "unprocessed_bulk.csv",
        "cell_raw": source / "Unprocessed_cell_type_clair_unnormalized.csv",
        "bulk_final": source / "Final_bulk_Bulk_lipids_cleaned_normalized_median_527_log210.csv",
        "cell_final": source / "Final_cell_cell_lipids_cleaned_norm_median_286_log210.csv",
    }
    bulk = pd.read_csv(files["bulk_raw"], index_col=0).T
    cell = pd.read_csv(files["cell_raw"], index_col=0).T
    bulk_final = pd.read_csv(files["bulk_final"], index_col=0)
    cell_final = pd.read_csv(files["cell_final"], index_col=0)
    native = pd.read_csv(native_path, index_col=0)
    bulk.columns = [normalize_mode(c) for c in bulk.columns]

    profiles = cell.index.intersection(cell_final.index)
    lipids = cell.columns.intersection(cell_final.columns)
    raw, final = cell.loc[profiles, lipids], cell_final.loc[profiles, lipids]
    # Reproduce the supplied PDF's stated sorted-cell steps for diagnosis only.
    cleaned = cell.where(cell >= 0)
    cleaned = cleaned.T.fillna(cleaned.median(axis=1)).T.fillna(0)
    stated = np.log2(1 + 10 * cleaned)
    correlations = raw.corrwith(final).dropna()

    # Final bulk IDs strip mode: compare only features unambiguous after stripping.
    base_bulk = bulk.copy()
    base_bulk.columns = base_bulk.columns.str.replace(r"_[PN]$", "", regex=True)
    unambiguous = base_bulk.columns[~base_bulk.columns.duplicated(keep=False)]
    bulk_lipids = unambiguous.intersection(bulk_final.columns)
    bulk_profiles = base_bulk.index.intersection(bulk_final.index)
    expected_bulk = np.log2(1 + 10 * np.exp2(base_bulk.loc[bulk_profiles, bulk_lipids]))
    bulk_error = bulk_final.loc[bulk_profiles, bulk_lipids] - expected_bulk
    discrepancies = bulk_error.stack().abs().sort_values(ascending=False)
    discrepancies.rename("absolute_log2_error").to_csv(output / "bulk_formula_discrepancies.csv")

    anchored, parameters, error = calibrate(native, bulk)
    anchored.index.name = "sample"
    anchored.to_csv(output / "PROVISIONAL_bulk_referenced_proxy_log2.csv")
    np.exp2(anchored).to_csv(output / "PROVISIONAL_bulk_referenced_proxy_linear.csv")
    anchored.loc[:, ~anchored.isna().any()].to_csv(output / "PROVISIONAL_complete_proxy_log2.csv")
    parameters.index.name = "lipid"
    parameters.to_csv(output / "calibration_parameters.csv")
    bulk.index.name = "sample"
    bulk.to_csv(output / "bulk_native_log2_mode_preserved.csv")

    # Test native CSV versus the supplied cell CSV without losing ion-mode identity.
    native_to_rna_style = native.copy()
    native_to_rna_style.index = ["D" + x.split("_")[1].zfill(3) + "_" + x.split("_")[2]
                                for x in native.index]
    matches = []
    for lipid in native.columns:
        base = lipid[:-2] if lipid.endswith(("_P", "_N")) else lipid
        if base in cell.columns:
            rows = native_to_rna_style.index.intersection(cell.index)
            diff = (native_to_rna_style.loc[rows, lipid] - cell.loc[rows, base]).dropna()
            if len(diff):
                matches.append({"feature": lipid, "n": len(diff), "max_abs_error": float(abs(diff).max())})
    pd.DataFrame(matches).to_csv(output / "native_to_supplied_cell_checks.csv", index=False)

    report = {
        "inputs_sha256": {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in files.items()},
        "bulk_raw_shape_samples_features": list(bulk.shape),
        "cell_raw_shape_samples_features": list(cell.shape),
        "negative_sorted_values": int((cell < 0).sum().sum()),
        "sorted_range": [float(cell.min().min()), float(cell.max().max())],
        "stated_protocol_sorted_output_range": [float(stated.min().min()), float(stated.max().max())],
        "provided_final_sorted_range": [float(cell_final.min().min()), float(cell_final.max().max())],
        "sorted_final_comparison_shape": list(raw.shape),
        "sorted_raw_final_per_lipid_correlation_quantiles": correlations.quantile([0, .25, .5, .75, 1]).to_dict(),
        "sorted_raw_final_negative_correlations": int((correlations < 0).sum()),
        "direction_reversal_example": {
            "lipid": "CE(18:2)", "contrast": "D019_END minus D019_EPI",
            "native_log2_difference": float(cell.loc["D019_END", "CE(18:2)"] - cell.loc["D019_EPI", "CE(18:2)"]),
            "provided_final_difference": float(cell_final.loc["D019_END", "CE(18:2)"] - cell_final.loc["D019_EPI", "CE(18:2)"]),
        },
        "bulk_protocol_comparison_values": int(bulk_error.size),
        "bulk_protocol_values_within_0_101_log2": int((bulk_error.abs() <= .101).sum().sum()),
        "bulk_protocol_max_absolute_error": float(bulk_error.abs().max().max()),
        "bulk_profiles_missing_from_final": bulk.index.difference(bulk_final.index).tolist(),
        "calibration": {
            "status": "PROVISIONAL_REFERENCE_PROXY_NOT_VALIDATED_ABUNDANCE",
            "formula": "proxy_log2 = native_sorted_log2 - median_PMX_native_log2 + median_bulk_log2",
            "shared_mode_preserved_features": int(anchored.shape[1]),
            "complete_shared_features": int((~anchored.isna().any()).sum()),
            "round_trip_max_error": error,
            "linear_range": [float(np.exp2(anchored).min().min()), float(np.exp2(anchored).max().max())],
            "interpretation": "Positive bulk-referenced arbitrary units. PMX is an unsorted-cell reference, not a known physical match to intact tissue.",
            "limitations": [
                "The protocol documents bulk as log2; original bulk measurement and normalization metadata are not independently established.",
                "Same transformation does not establish same assay scale, cell/protein/tissue denominator or matrix response.",
                "Per-lipid anchoring preserves within-lipid folds exactly and cannot fix attenuated or reversed training effects.",
                "Cross-lipid quantities, concentrations and cell-type contributions to tissue totals are not identifiable here.",
                "No RNA model was retrained or validated and no production bundle was changed.",
                "Shared D labels do not establish specimen pairing across studies.",
                "Only literal structural annotation matches with matching ionization mode were used; unmatched features need annotation review.",
            ],
        },
    }
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--native-targets", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.source_dir, args.native_targets, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
