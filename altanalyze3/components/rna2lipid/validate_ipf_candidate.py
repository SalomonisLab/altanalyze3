"""Fit isolated all-reference candidates and compare corrected lung IPF cohorts.

No IPF lipid outcome is used for fitting, calibration or model selection.
Native transformed lipid tissues are aggregated by donor before primary tests.
"""
from __future__ import annotations

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('validate_ipf_candidate.py')


import argparse
import csv
import gzip
import hashlib
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from scipy.stats import mannwhitneyu, spearmanr, t as student_t
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

from reference_normalization import fit_bundle, preservation_gate
from ipf_utils import metadata, map_expression, donor_means, welch_pvalues, fold_table
from ipf_lung_targets import load_lung, compare


def run(source, reference, raw_reference, supplied_reference, out, gene_annotation, bootstraps=2000):
    reject_retired_workflow('validate_ipf_candidate.py:run')
    source, reference, raw_reference, out = map(Path, [source, reference, raw_reference, out])
    out.mkdir(parents=True, exist_ok=True)
    evidence = out / "source_evidence"
    meta = metadata(source / "GSE213001_series_matrix.txt.gz")
    meta.to_csv(out / "RNA_sample_metadata.csv")
    log_path = source / "GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-logRPKMs-normalised.csv.gz"
    counts_path = source / "GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz"
    published = pd.read_csv(log_path, index_col=0)
    counts = pd.read_csv(counts_path, index_col=0)
    annotation = pd.read_csv(gene_annotation, sep="\t")
    if annotation.gene_id.duplicated().any():
        raise ValueError("Ambiguous gene annotation")
    mapping = annotation.set_index("gene_id").gene.to_dict()
    RNA = map_expression(published, mapping, True)
    raw = map_expression(counts, mapping, False)
    effective_lib = meta["lib.size"].astype(float) * meta["norm.factors"].astype(float)
    # Explicit sensitivity representation; not an inferred recreation of RPKM.
    logcpm = np.log2((raw + .5).div(effective_lib.loc[raw.index] + 1., axis=0) * 1e6)
    X = pd.read_csv(reference / "paired_RNA.csv", index_col=0)
    native = pd.read_csv(reference / "native_lipid_targets.csv", index_col=0)
    training_meta = pd.read_csv(reference / "sample_metadata.csv", index_col=0)
    genes = X.columns.intersection(RNA.columns).intersection(logcpm.columns).tolist()
    pd.Series(X.columns.difference(genes), name="missing_IPF_gene").to_csv(out / "excluded_RNA_genes.csv", index=False)
    X = X[genes]
    X.to_csv(out / "candidate_training_RNA.csv")
    RNA[genes].to_csv(out / "IPF_publisher_logRPKM_model_genes.csv")
    donor_groups = meta.groupby("donorid").diseasegroup.first()
    if (meta.groupby("donorid").diseasegroup.nunique() != 1).any():
        raise ValueError("Inconsistent donor disease label")
    selected = meta.diseasegroup.isin(["IPF", "NDC"]) & meta.lunglocation.isin(["Apex", "Base"])
    meta.loc[selected].to_csv(out / "primary_RNA_sample_metadata.csv")
    # Full RNA disease contrasts, independent of the model feature subset.
    for region in [None, "Apex", "Base"]:
        donor_log = donor_means(RNA, meta, region)
        donor_linear = donor_means(np.exp2(RNA), meta, region)
        groups = donor_groups.loc[donor_log.index]
        a, b = donor_log.loc[groups.eq("IPF")], donor_log.loc[groups.eq("NDC")]
        pvals, testable = welch_pvalues(a, b)
        arithmetic = np.log2(donor_linear.loc[groups.eq("IPF")].mean() / donor_linear.loc[groups.eq("NDC")].mean())
        pd.DataFrame({"geometric_log2FC": a.mean() - b.mean(), "arithmetic_log2FC": arithmetic,
                      "pvalue": pvals, "FDR_all_mapped_RNA_genes": multipletests(pvals, method="fdr_bh")[1],
                      "nonzero_within_group_variance": testable,
                      "candidate_input_gene": RNA.columns.isin(genes)}).rename_axis("gene").to_csv(
                          out / f"RNA_{region or 'donor_balanced'}_IPF_vs_control_all_genes.csv")
    rna_soup = BeautifulSoup((evidence / "rna_paper.xml").read_text(), "xml")
    rna_doi = rna_soup.find("article-id", attrs={"pub-id-type": "doi"}).text
    rng = np.random.default_rng(20261001)
    supplied_path, measured_stats, annotation_table, source_audit = load_lung(source, out, rng, bootstraps)
    results, gates, snapshots = {}, {}, {}
    supplied_reference = Path(supplied_reference)
    remapped = pd.read_csv(supplied_reference / "combined_sorted_and_supplied_bulk_MS1_log2.csv", index_col=0)
    remapped.index = remapped.index.str.replace("_BULK", "", regex=False)
    remapped_native = pd.read_csv(supplied_reference / "supplied_bulk_inverted_native_log2.csv", index_col=0)
    source_native = native[remapped.columns].copy()
    common_bulk = source_native.index.intersection(remapped_native.index)
    source_native.loc[common_bulk] = remapped_native.loc[common_bulk]
    missing_reference_RNA = remapped.index.difference(X.index)
    pd.Series(missing_reference_RNA, name="reference_profile_without_matched_RNA").to_csv(out / "remapped_reference_profiles_without_RNA.csv", index=False)
    pd.Series(X.index.difference(remapped.index), name="earlier_reference_profile_absent_from_final_bulk").to_csv(out / "earlier_reference_profiles_absent_from_final_bulk.csv", index=False)
    targets = {
        "supplied_bulk_MS1_47": remapped.loc[remapped.index.intersection(X.index)],
        "MS1_47": pd.read_csv(raw_reference / "combined_reference_MS1_log2.csv", index_col=0),
        "bulk_gauge_219_exploratory": pd.read_csv(reference / "combined_training_reference_log2.csv", index_col=0)}
    training_records = {}
    for label, Y in targets.items():
        directory = out / label
        directory.mkdir(exist_ok=True)
        ids = X.index.intersection(Y.index)
        xx, md = X.loc[ids], training_meta.loc[ids]
        Y = Y.loc[ids]
        xx.to_csv(directory / "candidate_training_RNA.csv")
        md.to_csv(directory / "candidate_training_sample_metadata.csv")
        training_records[label] = {"paired_samples": len(xx), "sorted_profiles": int(md.dataset.eq("sorted").sum()),
                                  "bulk_profiles": int(md.dataset.eq("bulk").sum()), "donors": int(md.donor.nunique()),
                                  "genes": len(genes), "heldout_samples": 0, "alpha": 100}
        Y.to_csv(directory / "candidate_training_lipids_log2.csv")
        check_source = source_native if label == "supplied_bulk_MS1_47" else native
        gate = preservation_gate(check_source.loc[ids, Y.columns], Y, md)
        if not gate["pass"]:
            raise ValueError("Candidate reference changed measured differentials")
        gates[label] = gate
        model, sx, sy = fit_bundle(xx, Y, "all_reference_candidate", pd.Series(0., index=Y.columns), 100., directory / "candidate_bundle.pkl", [])
        fit_prediction = pd.DataFrame(sy.inverse_transform(model.predict(sx.transform(xx))), index=xx.index, columns=Y.columns)
        fit_prediction.to_csv(directory / "training_predictions_log2.csv")
        np.exp2(fit_prediction).to_csv(directory / "training_predictions_linear.csv")
        training_records[label]["in_sample_RMSE_log2_abundance_not_validation"] = float(np.sqrt(np.mean((fit_prediction - Y).to_numpy() ** 2)))
        path = directory / "candidate_bundle.pkl"
        with path.open("rb") as handle:
            bundle = pickle.load(handle)
        bundle["metadata"].update({"candidate_only": True, "production_default_changed": False,
            "donor_disjoint_validation": False, "training_samples": len(xx), "training_donors": int(md.donor.nunique()),
            "heldout_samples": 0, "hyperparameter_policy": "alpha=100 frozen from prior reference work; no IPF lipid-informed tuning",
            "internally_supported_lipids": [], "RNA_input_scale": "Publisher-normalized log2 RPKM; per-gene healthy-control additive alignment",
            "calibration_interpretation": "Reference-normalized MS1 signal, not molar concentration" if label != "bulk_gauge_219_exploratory" else "Relative bulk gauge; additional features lack raw-MS replication support",
            "external_validation": "Unmatched whole-lung RNA and lung lipid cohorts; see IPF validation report",
            "prediction_policy": "Experimental all-feature output for validation; no IPF-validated lipid support inferred"})
        with path.open("wb") as handle:
            pickle.dump(bundle, handle)
        for representation, expression in [("publisher_logRPKM", RNA), ("TMM_logCPM_sensitivity", logcpm)]:
            expression = expression.loc[meta.index, genes]
            # Only healthy RNA is used to align gene baselines. A single offset
            # per gene applies identically to ALL groups. Slopes are never fitted.
            healthy = donor_means(expression, meta).loc[donor_groups.eq("NDC")]
            offset = xx.loc[md.dataset.eq("bulk")].mean() - healthy.mean()
            bundle.setdefault("external_RNA_calibration", {})[representation] = {
                "additive_gene_offsets": offset.to_dict(), "RNA_study": "GSE213001",
                "policy": "Add same gene offsets to every sample; never recompute by disease group"}
            aligned = expression + offset
            offset.to_csv(directory / f"{representation}_RNA_offsets.csv", header=["common_added_log2_offset"])
            prediction = pd.DataFrame(sy.inverse_transform(model.predict(sx.transform(aligned))), index=expression.index, columns=Y.columns)
            direct = pd.DataFrame(sy.inverse_transform(model.predict(sx.transform(expression))), index=expression.index, columns=Y.columns)
            # Additive RNA baseline shifts must change only per-lipid intercepts.
            translated = prediction - direct
            translation_error = float(np.max(abs((translated - translated.iloc[0]).to_numpy())))
            if translation_error > 1e-8:
                raise ValueError("RNA baseline alignment altered predicted log differentials")
            prediction.to_csv(directory / f"{representation}_all_samples_predictions_log2.csv")
            np.exp2(prediction).to_csv(directory / f"{representation}_all_samples_predictions_linear.csv")
            z = sx.transform(aligned.loc[selected])
            results[f"{label}/{representation}/RNA_QC"] = {
                "baseline_alignment_prediction_translation_error": translation_error,
                "fraction_gene_sample_values_abs_training_z_gt_3": float(np.mean(abs(z) > 3)),
                "median_profile_fraction_abs_training_z_gt_3": float(np.median(np.mean(abs(z) > 3, axis=1)))}
            for region in [None, "Apex", "Base"]:
                contrast = region or "donor_balanced"
                linear_donors = donor_means(np.exp2(prediction), meta, region)
                log_donors = donor_means(prediction, meta, region)
                groups = donor_groups.loc[linear_donors.index]
                stats = fold_table(linear_donors, log_donors, groups, rng, bootstraps)
                linear_donors.to_csv(directory / f"{representation}_{contrast}_donor_predictions_linear.csv")
                stats.to_csv(directory / f"{representation}_{contrast}_predicted_IPF_vs_control.csv")
                comparison, match_audit, numeric_summary = compare(measured_stats, annotation_table, stats)
                match_audit.to_csv(directory / "lipid_match_audit.csv", index=False)
                comparison.to_csv(directory / f"{representation}_{contrast}_measured_vs_predicted.csv", index=False)
                results[f"{label}/{representation}/{contrast}"] = {
                    "RNA_IPF_donors": int(groups.eq("IPF").sum()), "RNA_control_donors": int(groups.eq("NDC").sum()),
                    **numeric_summary}
                if representation == "publisher_logRPKM" and region is None:
                    snapshots[label] = comparison
        with path.open("wb") as handle:
            pickle.dump(bundle, handle)
    # External phenotype sanity checks use ALL available mapped genes, independent
    # of which genes were eligible for model fitting.
    marker_logs = donor_means(RNA, meta)
    marker_linear = donor_means(np.exp2(RNA), meta)
    markers = [g for g in ["MMP7", "COL1A1", "COL3A1", "ACTA2", "SFTPC", "SFTPA2"] if g in RNA.columns]
    fold_table(marker_linear[markers], marker_logs[markers], donor_groups.loc[marker_logs.index], rng, bootstraps).to_csv(out / "RNA_phenotype_sanity_checks.csv")
    report = {"candidate_only": True, "official_model_replaced": False,
              "primary_candidate": "supplied_bulk_MS1_47", "training_by_candidate": training_records,
              "RNA": {"accession": "GSE213001", "paper_DOI": rna_doi,
                      "IPF_donors": 20, "control_donors": 14, "whole_lung_profiles_primary": int(selected.sum())},
              "measured_lipids": source_audit, "comparisons": results,
              "normalization_preservation_gates": gates,
              "fold_definition": "Primary: ratio of arithmetic donor mean linear signals, conditional on log2 source scale. Geometric and supplied source log-effects are saved separately.",
              "superseded_results": "../../delete/IPF_plasma_20261001",
              "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [supplied_path, log_path, counts_path,
                  source / "GSE213001_series_matrix.txt.gz", Path(gene_annotation), reference / "paired_RNA.csv",
                  reference / "combined_training_reference_log2.csv", raw_reference / "combined_reference_MS1_log2.csv",
                  supplied_reference / "combined_sorted_and_supplied_bulk_MS1_log2.csv"]}}
    (out / "IPF_candidate_validation_result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    from render_ipf_report import render_report
    render_report(out)
    return report

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="/Users/saljh8/Downloads/Lipidomics/IPF")
    parser.add_argument("--reference", default="components/rna2lipid/artifacts/reference_normalization_20261001")
    parser.add_argument("--raw-reference", default="components/rna2lipid/artifacts/MS1_reference_normalization_20261001")
    parser.add_argument("--supplied-reference", default="components/rna2lipid/artifacts/supplied_bulk_MS1_remapping_20261001")
    parser.add_argument("--out", default="components/rna2lipid/artifacts/IPF_candidate_validation_20261001")
    parser.add_argument("--gene-annotation", default="components/fastCNV/resources/Hs_Ensembl_GRCh38_genes.tsv")
    parser.add_argument("--bootstraps", type=int, default=2000)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        result = run(args.source, args.reference, args.raw_reference, args.supplied_reference, args.out, args.gene_annotation, args.bootstraps)
    print(json.dumps({"training": result["training_by_candidate"], "primary": result["comparisons"]["supplied_bulk_MS1_47/publisher_logRPKM/donor_balanced"]}, indent=2))
