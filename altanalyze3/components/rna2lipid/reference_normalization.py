"""Compare normalization methods under an exact differential-preservation gate.

Run with /usr/bin/python3. All fitting, calibration and alpha selection in
validation exclude the held-out donor. Raw cohort MS is not a prerequisite.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import pickle
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

from prepare_lungmap_targets import extract
from audit_bulk_calibration import normalize_mode, audit as audit_legacy_preprocessing

METHODS = ["no_alignment", "global_pmx_shift", "lipid_pmx_shift",
           "shared_donor_pmx_shift", "sample_median", "quantile", "population_matching"]
ELIGIBLE = ["global_pmx_shift", "lipid_pmx_shift", "shared_donor_pmx_shift"]
ALPHAS = np.array([1., 10., 100., 1000., 10000.])


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def donor_id(sample):
    return str(sample).split("_")[0]


def rna_id(sample):
    if str(sample).startswith("Sample_"):
        _, number, population = str(sample).split("_")
        return f"D{int(number):03d}_{population}"
    match = re.fullmatch(r"D(\d+)_(\w+)", str(sample))
    if match:
        return f"D{int(match[1]):03d}_{match[2]}"
    return str(sample)


def load_data(source, out):
    source, out = Path(source), Path(out)
    extract(source / "LMEX0000003692_data_table.xlsx", out / "source_targets")
    cell = pd.read_csv(out / "source_targets/targets_native_log2.csv", index_col=0)
    native_sample_median_error = float(cell.median(axis=1).abs().max())
    cell.index = [rna_id(s) for s in cell.index]
    bulk = pd.read_csv(source / "Preprocessing Process/unprocessed_bulk.csv", index_col=0).T
    bulk.columns = [normalize_mode(c) for c in bulk.columns]
    crna = pd.read_csv(source / "newrna_cell_clair_filtered_symbol.csv", index_col=0).T
    brna = pd.read_csv(source / "feature_blankreduiction.csv", index_col=0)
    duplicate_genes = crna.columns[crna.columns.duplicated(keep=False)]
    pd.Series(duplicate_genes).value_counts().rename_axis("gene").rename("source_rows").to_csv(out / "collapsed_RNA_gene_symbols.csv")
    # Matches the established predictor's gene-symbol aggregation on supplied scale.
    crna = crna.T.groupby(level=0, sort=True).mean().T
    for label, df in [("cell lipids", cell), ("bulk lipids", bulk), ("cell RNA", crna), ("bulk RNA", brna)]:
        if not df.index.is_unique or not df.columns.is_unique:
            raise ValueError(f"Duplicate IDs in {label}; resolve explicitly.")
        if not np.isfinite(df.to_numpy(dtype=float)[~df.isna().to_numpy()]).all():
            raise ValueError(f"Nonfinite measurements in {label}.")
    genes = sorted(crna.columns.intersection(brna.columns))
    lipids = sorted(cell.columns.intersection(bulk.columns))
    # Complete observed targets: never fabricate lipid measurements by imputation.
    keep = cell[lipids].notna().all() & bulk[lipids].notna().all()
    dropped = list(keep.index[~keep])
    lipids = list(keep.index[keep])
    cell, bulk = cell[lipids], bulk[lipids]
    ci, bi = cell.index.intersection(crna.index), bulk.index.intersection(brna.index)
    X = pd.concat([crna.loc[ci, genes], brna.loc[bi, genes]])
    Y = pd.concat([cell.loc[ci], bulk.loc[bi]])
    if not X.index.is_unique or X.isna().any().any():
        raise ValueError("Ambiguous/missing matched RNA.")
    meta = pd.DataFrame(index=X.index)
    meta["donor"] = [donor_id(s) for s in X.index]
    meta["dataset"] = ["sorted" if s in ci else "bulk" for s in X.index]
    meta["population"] = [s.split("_")[-1] if s in ci else "BULK" for s in X.index]
    X.index.name = Y.index.name = meta.index.name = "sample"
    X.to_csv(out / "paired_RNA.csv")
    Y.to_csv(out / "native_lipid_targets.csv")
    meta.to_csv(out / "sample_metadata.csv")
    pd.DataFrame({"lipid": dropped, "reason": "missing observed sorted or bulk measurement"}).to_csv(out / "excluded_lipids.csv", index=False)
    hashes = {}
    for path in [source / "LMEX0000003692_data_table.xlsx", source / "newrna_cell_clair_filtered_symbol.csv",
                 source / "feature_blankreduiction.csv", source / "Preprocessing Process/unprocessed_bulk.csv"]:
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"input_sha256": hashes, "paired_profiles": len(X), "sorted_profiles": len(ci),
                "bulk_profiles": len(bi), "donors": int(meta.donor.nunique()), "genes": len(genes),
                "mode_matched_lipids": len(keep), "complete_lipids": len(lipids),
                "excluded_incomplete_lipids": len(dropped),
                "duplicate_RNA_gene_rows_collapsed": int(len(duplicate_genes) - len(set(duplicate_genes))),
                "original_sorted_sample_medians_max_absolute": native_sample_median_error,
                "bulk_profiles_without_RNA": bulk.index.difference(brna.index).tolist(),
                "sorted_profiles_without_RNA": cell.index.difference(crna.index).tolist(),
                "RNA_representation": "Exact supplied numeric representation; no inferred log conversion or library normalization.",
                "feature_policy": "Exact structural annotation and ion mode. No mode merging, negative clipping or target imputation.",
                "identifier_policy": "Sample_N_POP -> DNNN_POP, corroborated by provided original cell CSV value matches; RNA .1 profiles are excluded.",
                "validation_policy": "A held-out donor is excluded from BOTH modalities, RNA fitting, calibration, and alpha selection."}
    dump(out / "manifest.json", manifest)
    return X, Y, meta, cell, bulk, crna, manifest


def fit_offsets(method, cell, bulk, exclude=()):
    cs = cell.loc[~cell.index.map(donor_id).isin(exclude)]
    bs = bulk.loc[~bulk.index.isin(exclude)]
    pmx = cs.loc[cs.index.str.endswith("_PMX")]
    if pmx.empty or bs.empty:
        raise ValueError("Insufficient independent calibration references.")
    if method == "no_alignment":
        return pd.Series(0., index=cell.columns)
    offset = bs.median() - pmx.median()
    if method == "global_pmx_shift":
        return pd.Series(float(offset.median()), index=cell.columns)
    if method == "shared_donor_pmx_shift":
        pairs = []
        for s in pmx.index:
            if donor_id(s) in bs.index:
                pairs.append(bs.loc[donor_id(s)] - pmx.loc[s])
        if not pairs:
            raise ValueError("No shared-donor calibration pairs.")
        return pd.DataFrame(pairs).median()
    return offset


def transform(Y, meta, offsets, method):
    result = Y.copy()
    result.loc[meta.dataset == "sorted"] += offsets
    if method == "sample_median":
        result = result.sub(result.median(axis=1), axis=0)
    elif method == "population_matching":
        target = result.loc[meta.dataset == "bulk"].median()
        for population in sorted(meta.loc[meta.dataset == "sorted", "population"].unique()):
            rows = meta.population == population
            result.loc[rows] += target - result.loc[rows].median()
    elif method == "quantile":
        values = result.to_numpy()
        reference = np.sort(values, axis=1).mean(axis=0)
        # Average rank for ties; interpolate the empirical shared distribution.
        ranks = result.rank(axis=1, method="average").to_numpy() - 1
        result.iloc[:, :] = np.interp(ranks, np.arange(len(reference)), reference)
    return result


def preservation_gate(native, combined, meta, tolerance=1e-9):
    records, max_error = [], 0.
    for dataset in ["sorted", "bulk"]:
        ids = meta.index[meta.dataset == dataset]
        before, after = native.loc[ids].values, combined.loc[ids].values
        i, j = np.triu_indices(len(ids), 1)
        error = float(np.max(np.abs((after[i] - after[j]) - (before[i] - before[j]))))
        max_error = max(max_error, error)
        records.append({"dataset": dataset, "sample_pairs": len(i), "lipids": native.shape[1],
                        "differentials_checked": len(i) * native.shape[1], "max_abs_log2fc_error": error})
    return {"pass": max_error <= tolerance, "tolerance": tolerance, "max_abs_log2fc_error": max_error,
            "checks": records}


def population_differentials(native, combined, meta, method):
    rows = []
    populations = sorted(meta.loc[meta.dataset == "sorted", "population"].unique())
    for case, ctrl in combinations(populations, 2):
        a = native.loc[meta.population == case].copy()
        b = native.loc[meta.population == ctrl].copy()
        c = combined.loc[a.index].copy()
        d = combined.loc[b.index].copy()
        for frame in [a, b, c, d]:
            frame.index = frame.index.map(donor_id)
        donors = sorted(a.index.intersection(b.index))
        fc0, fc1 = (a.loc[donors] - b.loc[donors]).mean(), (c.loc[donors] - d.loc[donors]).mean()
        p0 = ttest_rel(a.loc[donors], b.loc[donors], axis=0).pvalue
        p1 = ttest_rel(c.loc[donors], d.loc[donors], axis=0).pvalue
        fdr0, fdr1 = multipletests(p0, method="fdr_bh")[1], multipletests(p1, method="fdr_bh")[1]
        for j, lipid in enumerate(native.columns):
            rows.append([method, case, ctrl, lipid, len(donors), fc0[lipid], fc1[lipid],
                         fc1[lipid] - fc0[lipid], p0[j], p1[j], fdr0[j], fdr1[j]])
    return pd.DataFrame(rows, columns=["method", "case", "control", "lipid", "donors", "native_log2fc",
                                      "combined_log2fc", "difference", "native_pvalue", "combined_pvalue",
                                      "native_FDR", "combined_FDR"])


def verify_reference_artifacts(out):
    """Recheck saved references without refitting any predictive model."""
    out = Path(out)
    native = pd.read_csv(out / "native_lipid_targets.csv", index_col=0)
    meta = pd.read_csv(out / "sample_metadata.csv", index_col=0)
    checks, tables = {}, []
    for method in METHODS:
        combined = pd.read_csv(out / f"reference_{method}.csv", index_col=0)
        gate = preservation_gate(native, combined, meta)
        table = population_differentials(native, combined, meta, method)
        p_error = float(abs(table.native_pvalue - table.combined_pvalue).max())
        fdr_error = float(abs(table.native_FDR - table.combined_FDR).max())
        before = (table.native_FDR < .05) & (abs(table.native_log2fc) >= np.log2(1.5))
        after = (table.combined_FDR < .05) & (abs(table.combined_log2fc) >= np.log2(1.5))
        gate.update({"population_tests_checked": len(table), "max_pvalue_error": p_error,
                     "max_FDR_error": fdr_error, "significant_tests_before": int(before.sum()),
                     "significant_tests_after": int(after.sum()), "changed_significance_calls": int((before != after).sum()),
                     "eligible_for_training": method in ELIGIBLE and gate["pass"]})
        if method in ELIGIBLE and (not gate["pass"] or p_error > 1e-9 or fdr_error > 1e-9 or (before != after).any()):
            raise ValueError("Saved combined reference fails differential/FDR preservation: " + method)
        checks[method] = gate
        tables.append(table)
    pd.concat(tables).to_csv(out / "population_differential_checks.csv", index=False)
    dump(out / "differential_preservation.json", checks)
    selected = json.loads((out / "workflow_result.json").read_text())["selected_method"]
    final = pd.read_csv(out / "combined_training_reference_log2.csv", index_col=0)
    expected = pd.read_csv(out / f"reference_{selected}.csv", index_col=0)
    np.testing.assert_allclose(final, expected, atol=1e-12, rtol=0)
    linear = pd.read_csv(out / "combined_training_reference_linear.csv", index_col=0)
    if not (linear.values > 0).all() or not np.isfinite(linear.values).all():
        raise ValueError("Reference abundances must be strictly positive and finite.")
    np.testing.assert_allclose(np.log2(linear.values), final.values, atol=1e-12, rtol=0)
    report = json.loads((out / "workflow_result.json").read_text())
    report["mandatory_differential_gate"] = checks[selected]
    dump(out / "workflow_result.json", report)
    return checks


def calibration_bridge_checks(cell, bulk, out):
    """Assess shared donor identifiers as candidate bridges, leaving each out."""
    rows = []
    shared = sorted({donor_id(s) for s in cell.index if s.endswith("_PMX")} & set(bulk.index))
    for held in shared:
        sample = held + "_PMX"
        for method in ELIGIBLE:
            offsets = fit_offsets(method, cell, bulk, {held})
            residual = cell.loc[sample] + offsets - bulk.loc[held]
            for lipid in cell.columns:
                rows.append([method, held, lipid, float(residual[lipid])])
    table = pd.DataFrame(rows, columns=["method", "held_out_shared_ID", "lipid", "PMX_minus_bulk_log2"])
    table.to_csv(Path(out) / "calibration_bridge_heldout_residuals.csv", index=False)
    report = {"shared_donor_identifiers": shared,
              "interpretation": "Candidate PMX/tissue bridge consistency; matching donor IDs do not prove matching specimen amounts or assay response.",
              "methods": {method: {"RMSE_log2": float(np.sqrt(np.mean(frame.PMX_minus_bulk_log2 ** 2))),
                                    "median_abs_error_log2": float(frame.PMX_minus_bulk_log2.abs().median())}
                          for method, frame in table.groupby("method")}}
    dump(Path(out) / "calibration_bridge_checks.json", report)
    return report


def loss_scale(y, meta):
    residual = y - y.groupby(meta.dataset).transform("mean")
    return np.maximum(residual.std(ddof=0).values, .2)


def fast_predict(xtrain, ytrain, xtest, alphas):
    sx = StandardScaler().fit(xtrain)
    a, b = sx.transform(xtrain), sx.transform(xtest)
    mean = ytrain.mean(axis=0)
    eigen, q = np.linalg.eigh(a @ a.T)
    eigen = np.maximum(eigen, 0)
    left = b @ a.T @ q
    right = q.T @ (ytrain - mean)
    return [mean + left @ (right / (eigen[:, None] + alpha)) for alpha in alphas]


def tune(X, Y, meta, cell, bulk, excluded):
    donors = sorted(meta.donor.unique())
    losses = {method: [] for method in ELIGIBLE}
    for held in donors:
        train, test = meta.donor != held, meta.donor == held
        scale = loss_scale(Y.loc[train], meta.loc[train])
        for method in ELIGIBLE:
            offsets = fit_offsets(method, cell, bulk, set(excluded) | {held})
            y = transform(Y, meta, offsets, method)
            predictions = fast_predict(X.loc[train].values, y.loc[train].values, X.loc[test].values, ALPHAS)
            losses[method].append([float(np.mean(((p - y.loc[test].values) / scale) ** 2)) for p in predictions])
    choices = {}
    for method in ELIGIBLE:
        score = np.mean(losses[method], axis=0)
        j = int(np.argmin(score))
        choices[method] = {"alpha": float(ALPHAS[j]), "inner_donor_loss": float(score[j]),
                           "alpha_losses": [float(x) for x in score]}
    return choices


def select_method(choices, cell, bulk, excluded):
    """Near-ties in predictive loss are resolved by held-out bridge consistency.

    The 1% practical equivalence margin is fixed for the workflow. Every bridge
    used here excludes the outer test donor; it is tuning evidence, not a new
    independent validation claim.
    """
    best_loss = min(c["inner_donor_loss"] for c in choices.values())
    contenders = [m for m in ELIGIBLE if choices[m]["inner_donor_loss"] <= best_loss * 1.01]
    shared = sorted(({donor_id(s) for s in cell.index if s.endswith("_PMX")} & set(bulk.index)) - set(excluded))
    for method in ELIGIBLE:
        residuals = []
        for donor in shared:
            offsets = fit_offsets(method, cell, bulk, set(excluded) | {donor})
            residuals.append((cell.loc[donor + "_PMX"] + offsets - bulk.loc[donor]).values)
        choices[method]["bridge_tuning_RMSE_log2"] = float(np.sqrt(np.mean(np.square(residuals)))) if residuals else None
        choices[method]["within_1pct_of_best_predictive_loss"] = method in contenders
    return min(contenders, key=lambda m: (choices[m]["bridge_tuning_RMSE_log2"]
                                          if choices[m]["bridge_tuning_RMSE_log2"] is not None else float("inf"),
                                          choices[m]["inner_donor_loss"]))


def safe_corr(a, b):
    if np.std(a) < 1e-10 or np.std(b) < 1e-10:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def effect_metrics(truth, predicted, meta):
    within_type_truth, within_type_pred, between_type_truth, between_type_pred = [], [], [], []
    rows = meta.loc[meta.dataset == "sorted"]
    for population in sorted(rows.population.unique()):
        ids = rows.index[rows.population == population]
        a, b = truth.loc[ids].values, predicted.loc[ids].values
        i, j = np.triu_indices(len(ids), 1)
        within_type_truth.append(a[i] - a[j])
        within_type_pred.append(b[i] - b[j])
    for donor in sorted(rows.donor.unique()):
        ids = rows.index[rows.donor == donor]
        a, b = truth.loc[ids].values, predicted.loc[ids].values
        i, j = np.triu_indices(len(ids), 1)
        between_type_truth.append(a[i] - a[j])
        between_type_pred.append(b[i] - b[j])
    report = {}
    arrays = {}
    for label, true_list, pred_list in [("within_population_donor_contrasts", within_type_truth, within_type_pred),
                                       ("within_donor_population_contrasts", between_type_truth, between_type_pred)]:
        a, b = np.concatenate(true_list), np.concatenate(pred_list)
        mask = np.abs(a) >= .5
        report[label] = {"comparisons": int(a.size), "correlation": safe_corr(a.ravel(), b.ravel()),
                         "RMSE_log2fc": float(np.sqrt(np.mean((a - b) ** 2))),
                         "direction_agreement_for_abs_truth_ge_0_5": float(np.mean(np.sign(a[mask]) == np.sign(b[mask]))) if mask.any() else None}
        arrays[label] = (a, b)
    return report, arrays


def evaluate(X, Y, meta, cell, bulk, out):
    full_offsets = {m: fit_offsets(m, cell, bulk) for m in ELIGIBLE}
    truth = {m: transform(Y, meta, full_offsets[m], m) for m in ELIGIBLE}
    predictions = {m: pd.DataFrame(np.nan, index=Y.index, columns=Y.columns) for m in ELIGIBLE}
    selected = pd.DataFrame(np.nan, index=Y.index, columns=Y.columns)
    baseline = pd.DataFrame(np.nan, index=Y.index, columns=Y.columns)
    # One common reporting gauge, independent of the fold's fitted calibration.
    common = truth["lipid_pmx_shift"]
    records = []
    for k, held in enumerate(sorted(meta.donor.unique()), 1):
        train, test = meta.donor != held, meta.donor == held
        choices = tune(X.loc[train], Y.loc[train], meta.loc[train], cell, bulk, {held})
        best = select_method(choices, cell, bulk, {held})
        record = {"held_out_donor": held, "selected_method": best, "choices": choices}
        for method in ELIGIBLE:
            fold_offsets = fit_offsets(method, cell, bulk, {held})
            yy = transform(Y, meta, fold_offsets, method)
            pred = fast_predict(X.loc[train].values, yy.loc[train].values, X.loc[test].values,
                                [choices[method]["alpha"]])[0]
            frame = pd.DataFrame(pred, index=Y.index[test], columns=Y.columns)
            sorted_ids = meta.index[test & (meta.dataset == "sorted")]
            frame.loc[sorted_ids] += full_offsets[method] - fold_offsets
            predictions[method].loc[frame.index] = frame
            if method == best:
                aligned = frame.copy()
                aligned.loc[sorted_ids] += full_offsets["lipid_pmx_shift"] - full_offsets[method]
                selected.loc[aligned.index] = aligned
        for population in meta.loc[test, "population"].unique():
            ids = meta.index[test & (meta.population == population)]
            reference = common.loc[train & (meta.population == population)].mean()
            baseline.loc[ids] = reference.values
        records.append(record)
        dump(out / "outer_fold_choices.json", records)
        print(f"Donor validation {k}/{meta.donor.nunique()}: {held}, inner-selected {best}", flush=True)
    summary = {}
    for method, pred in {**predictions, "nested_selected": selected, "population_mean_baseline": baseline}.items():
        target = truth.get(method, common)
        pred.to_csv(out / f"OOF_{method}.csv")
        effects, arrays = effect_metrics(target, pred, meta)
        mask = meta.dataset == "sorted"
        effects["sorted_RMSE_log2"] = float(np.sqrt(np.mean((target.loc[mask].values - pred.loc[mask].values) ** 2)))
        effects["bulk_RMSE_log2"] = float(np.sqrt(np.mean((target.loc[~mask].values - pred.loc[~mask].values) ** 2)))
        summary[method] = effects
    selected_effects, arrays = effect_metrics(common, selected, meta)
    sorted_rows = meta.dataset == "sorted"
    sq = (common.loc[sorted_rows] - selected.loc[sorted_rows]) ** 2
    bsq = (common.loc[sorted_rows] - baseline.loc[sorted_rows]) ** 2
    table = pd.DataFrame(index=Y.columns)
    table["sorted_RMSE_log2"] = np.sqrt(sq.mean())
    table["skill_vs_population_mean"] = 1 - sq.mean() / bsq.mean()
    for label, (a, b) in arrays.items():
        table[label + "_RMSE_log2fc"] = np.sqrt(np.mean((a - b) ** 2, axis=0))
        mask = abs(a) >= .5
        correct = (np.sign(a) == np.sign(b)) & mask
        n = mask.sum(axis=0)
        table[label + "_n_abs_truth_ge_0_5"] = n
        table[label + "_direction_agreement"] = np.divide(correct.sum(axis=0), n, out=np.full(len(n), np.nan), where=n > 0)
    # Fixed provisional gates, not disease validation or a selected-by-test model list.
    table["internal_support"] = ((table.skill_vs_population_mean > 0) & (table.sorted_RMSE_log2 <= 1.) &
                                 (table.within_population_donor_contrasts_direction_agreement >= .7))
    table.index.name = "lipid"
    table.to_csv(out / "per_lipid_validation.csv")
    dump(out / "donor_validation.json", summary)
    return summary, table


def fit_bundle(X, Y, method, offsets, alpha, out, supported_lipids=()):
    sx, sy = StandardScaler().fit(X), StandardScaler().fit(Y)
    model = KernelRidge(alpha=alpha, kernel="linear").fit(sx.transform(X), sy.transform(Y))
    bundle = {"model": model, "scaler_x": sx, "scaler_y": sy, "X_columns": list(X.columns),
              "Y_columns": list(Y.columns), "metadata": {
                  "model_name": "ReferenceCalibratedLinearKernelRidge", "target_scaling": {"mode": "standard"},
                  "normalization_method": method, "expression_scale": "log2", "inverse": "2**prediction; no pseudocount or clipping",
                  "RNA_input_scale": "Exact representation of supplied paired_RNA.csv",
                  "calibration_interpretation": "Bulk-referenced normalized lipid abundance; physical units not established",
                  "donor_disjoint_validation": True, "production_default_changed": False,
                  "internally_supported_lipids": list(supported_lipids),
                  "prediction_policy": "CLI defaults to internally supported lipids. Full output requires --include-unsupported. Support is provisional, not disease validation.",
                  "calibration_offsets": offsets.to_dict(), "alpha": alpha}}
    with Path(out).open("wb") as handle:
        pickle.dump(bundle, handle)
    return model, sx, sy


def sorted_only_comparison(X, Y, meta, cell, bulk, out):
    """Donor-held-out comparator establishes whether combining data helps."""
    keep = meta.dataset == "sorted"
    x, y, md = X.loc[keep], Y.loc[keep], meta.loc[keep]
    result = pd.DataFrame(np.nan, index=y.index, columns=y.columns)
    choices = []
    for held in sorted(md.donor.unique()):
        outer_train, outer_test = md.donor != held, md.donor == held
        errors = []
        for inner in sorted(md.loc[outer_train, "donor"].unique()):
            train = outer_train & (md.donor != inner)
            test = outer_train & (md.donor == inner)
            scale = loss_scale(y.loc[train], md.loc[train])
            predictions = fast_predict(x.loc[train].values, y.loc[train].values, x.loc[test].values, ALPHAS)
            errors.append([float(np.mean(((p - y.loc[test].values) / scale) ** 2)) for p in predictions])
        alpha = float(ALPHAS[int(np.argmin(np.mean(errors, axis=0)))])
        result.loc[outer_test] = fast_predict(x.loc[outer_train].values, y.loc[outer_train].values,
                                             x.loc[outer_test].values, [alpha])[0]
        choices.append({"held_out_donor": held, "alpha": alpha})
    common_offsets = fit_offsets("lipid_pmx_shift", cell, bulk)
    anchored = result + common_offsets
    anchored.to_csv(Path(out) / "OOF_sorted_only.csv")
    effects, _ = effect_metrics(y + common_offsets, anchored, md)
    effects["sorted_RMSE_log2"] = float(np.sqrt(np.mean((y.values - result.values) ** 2)))
    effects["fold_choices"] = choices
    dump(Path(out) / "sorted_only_comparison.json", effects)
    return effects


def external_validation(X, Y, meta, cell, bulk, crna, reference_path, out):
    reference = pd.read_csv(reference_path, index_col=0)
    reference.index = [rna_id(s) for s in reference.index]
    columns = {}
    for name in reference.columns:
        lipid, adduct = name.split("|", 1)
        if adduct.endswith("+"):
            canonical = lipid + "_P"
        elif adduct.endswith("-"):
            canonical = lipid + "_N"
        else:
            continue
        if canonical in columns:
            raise ValueError("Ambiguous external feature annotation.")
        columns[canonical] = name
    shared = [c for c in Y.columns if c in columns]
    matched = reference.index.intersection(crna.index)
    ext_donors = {donor_id(s) for s in reference.index}
    train = ~meta.donor.isin(ext_donors)
    # Bulk profiles for these donors must also be excluded from training/calibration.
    choices = tune(X.loc[train], Y.loc[train], meta.loc[train], cell, bulk, ext_donors)
    method = select_method(choices, cell, bulk, ext_donors)
    offsets = fit_offsets(method, cell, bulk, ext_donors)
    yy = transform(Y, meta, offsets, method)
    pred = pd.DataFrame(fast_predict(X.loc[train].values, yy.loc[train].values,
                                    crna.loc[matched, X.columns].values, [choices[method]["alpha"]])[0],
                        index=matched, columns=Y.columns)
    rows = []
    for donor in sorted({donor_id(s) for s in matched}):
        pmx = donor + "_PMX"
        if pmx not in matched:
            continue
        for sample in matched:
            if donor_id(sample) != donor or sample == pmx:
                continue
            for lipid in shared:
                true = reference.loc[sample, columns[lipid]] - reference.loc[pmx, columns[lipid]]
                predicted = pred.loc[sample, lipid] - pred.loc[pmx, lipid]
                rows.append([donor, sample, pmx, lipid, true, predicted])
    contrasts = pd.DataFrame(rows, columns=["donor", "sample", "reference", "lipid", "measured_log2fc", "predicted_log2fc"])
    contrasts.to_csv(out / "MSV000081973_heldout_contrasts.csv", index=False)
    pred.to_csv(out / "MSV000081973_heldout_predictions.csv")
    a, b = contrasts.measured_log2fc.values, contrasts.predicted_log2fc.values
    mask = abs(a) >= .5
    summary = {"excluded_donors_from_both_modalities_and_calibration": sorted(ext_donors),
               "matched_profiles": len(matched), "unmatched_profiles": reference.index.difference(crna.index).tolist(),
               "shared_exact_mode_features": len(shared), "contrasts": len(contrasts), "selected_method": method,
               "RMSE_log2fc": float(np.sqrt(np.mean((a - b) ** 2))), "correlation": safe_corr(a, b),
               "direction_agreement_abs_truth_ge_0_5": float(np.mean(np.sign(a[mask]) == np.sign(b[mask]))),
               "evidence": "Processed independent sorted profiles; same-donor bulk samples excluded. Not raw-MS or physical concentration validation."}
    dump(out / "MSV000081973_validation.json", summary)
    return summary


def run(args):
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    X, Y, meta, cell, bulk, crna, manifest = load_data(args.source_dir, out)
    legacy = audit_legacy_preprocessing(Path(args.source_dir) / "Preprocessing Process",
                                      out / "source_targets/targets_native_log2.csv", out / "legacy_preprocessing_audit")
    gates, differential_tables = {}, []
    for method in METHODS:
        offsets = fit_offsets(method, cell, bulk)
        combined = transform(Y, meta, offsets, method)
        gate = preservation_gate(Y, combined, meta)
        gate["eligible_for_training"] = method in ELIGIBLE and gate["pass"]
        gates[method] = gate
        combined.to_csv(out / f"reference_{method}.csv")
        differential_tables.append(population_differentials(Y, combined, meta, method))
    dump(out / "differential_preservation.json", gates)
    pd.concat(differential_tables).to_csv(out / "population_differential_checks.csv", index=False)
    if not all(gates[m]["eligible_for_training"] for m in ELIGIBLE):
        raise ValueError("A proposed normalization failed the mandatory differential gate.")
    print("Differential gate complete; starting nested donor validation", flush=True)
    summary, lipid_metrics = evaluate(X, Y, meta, cell, bulk, out)
    sorted_only = sorted_only_comparison(X, Y, meta, cell, bulk, out)
    choices = tune(X, Y, meta, cell, bulk, set())
    method = select_method(choices, cell, bulk, set())
    offsets = fit_offsets(method, cell, bulk)
    yy = transform(Y, meta, offsets, method)
    final_gate = preservation_gate(Y, yy, meta)
    assert final_gate["pass"]
    yy.to_csv(out / "combined_training_reference_log2.csv")
    np.exp2(yy).to_csv(out / "combined_training_reference_linear.csv")
    offsets.rename("sorted_log2_offset").to_csv(out / "selected_calibration_offsets.csv")
    supported = lipid_metrics.index[lipid_metrics.internal_support].tolist()
    model, sx, sy = fit_bundle(X, yy, method, offsets, choices[method]["alpha"], out / "reference_calibrated_bundle.pkl", supported)
    # End-to-end serialization/API check with the supported existing predictor.
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from altanalyze3.components.rna2lipid.api import load_bundle
    loaded = load_bundle(out / "reference_calibrated_bundle.pkl")
    actual = loaded.predict_from_dataframe(X).predictions
    expected = sy.inverse_transform(model.predict(sx.transform(X)))
    serialization_error = float(np.max(abs(actual.values - expected)))
    if serialization_error > 1e-10:
        raise ValueError("Bundle API round trip failed.")
    actual.to_csv(out / "training_predictions_log2.csv")
    np.exp2(actual).to_csv(out / "training_predictions_linear.csv")
    crna = crna.loc[:, X.columns]
    extrapolation = loaded.predict_from_dataframe(crna).predictions
    extrapolation.to_csv(out / "supplied_sorted_RNA_predictions_log2.csv")
    np.exp2(extrapolation).to_csv(out / "supplied_sorted_RNA_predictions_linear.csv")
    extrapolation[supported].to_csv(out / "internally_supported_predictions_log2.csv")
    np.exp2(extrapolation[supported]).to_csv(out / "internally_supported_predictions_linear.csv")
    external = external_validation(X, Y, meta, cell, bulk, crna, args.external_reference, out)
    bridge = calibration_bridge_checks(cell, bulk, out)
    report = {"selected_method": method, "final_alpha": choices[method]["alpha"], "method_selection": choices,
              "method_selection_policy": "Lowest nested predictive loss; contenders within 1% are ranked by training-only leave-one-shared-ID-out bridge RMSE.",
              "mandatory_differential_gate": final_gate, "API_serialization_max_error": serialization_error,
              "nested_donor_validation": summary["nested_selected"], "external_validation": external,
              "sorted_only_comparison": sorted_only,
              "calibration_bridge_checks": bridge,
              "internally_supported_lipids": int(lipid_metrics.internal_support.sum()), "total_lipids": len(yy.columns),
              "production_default_changed": False,
              "source_chain_audit": {"sorted_direction_reversal": legacy["direction_reversal_example"],
                                     "sorted_raw_final_negative_correlations": legacy["sorted_raw_final_negative_correlations"],
                                     "bulk_formula_max_error": legacy["bulk_protocol_max_absolute_error"]},
              "important": ["Training-reference fold preservation and prediction accuracy are separate checks.",
                            "Cross-dataset absolute folds are not asserted from assay alignment alone.",
                            "Full-training predictions are resubstitution, not validation.",
                            "Internal support thresholds: positive skill against population mean, RMSE<=1 log2, within-population direction agreement>=0.7 for measured abs fold>=0.5.",
                            "External validation has only two RNA-matched donors and cannot establish disease transfer."]}
    dump(out / "workflow_result.json", report)
    verify_reference_artifacts(out)
    print(json.dumps(report, indent=2), flush=True)


def predict(args):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from altanalyze3.components.rna2lipid.api import load_bundle
    bundle = load_bundle(args.bundle)
    frame = pd.read_csv(args.input, index_col=0)
    if args.transpose:
        frame = frame.T
    if not frame.index.is_unique or not np.isfinite(frame.values).all():
        raise ValueError("Input requires finite values and unique sample IDs.")
    frame.columns = frame.columns.astype(str).str.strip()
    frame = frame.T.groupby(level=0, sort=True).mean().T
    missing = set(bundle.input_genes) - set(frame.columns)
    if missing:
        raise ValueError(f"Input lacks {len(missing)} required genes; no silent zero filling is permitted.")
    result = bundle.predict_from_dataframe(frame).predictions
    if not args.include_unsupported:
        supported = bundle.metadata.get("internally_supported_lipids", [])
        if not supported:
            raise ValueError("No internally supported lipid list in bundle; inspect validation before exporting all outputs.")
        result = result.loc[:, supported]
    result.to_csv(args.output_log2)
    linear = np.exp2(result)
    if not np.isfinite(linear.values).all() or not (linear.values > 0).all():
        raise ValueError("Exponentiation overflow/underflow; input may be outside training scope.")
    linear.to_csv(args.output_linear)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--source-dir", required=True)
    r.add_argument("--output-dir", required=True)
    r.add_argument("--external-reference", required=True)
    p = sub.add_parser("predict")
    p.add_argument("--bundle", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--transpose", action="store_true")
    p.add_argument("--output-log2", required=True)
    p.add_argument("--output-linear", required=True)
    p.add_argument("--include-unsupported", action="store_true", help="Explicitly export exploratory lipid outputs that fail internal support criteria.")
    v = sub.add_parser("verify")
    v.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.command == "run":
            run(args)
        elif args.command == "predict":
            predict(args)
        else:
            print(json.dumps(verify_reference_artifacts(args.output_dir), indent=2))


if __name__ == "__main__":
    main()
