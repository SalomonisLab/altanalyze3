"""Release, scale, differential and mean-aggregation regression checks."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests

from altanalyze3.components.rna2lipid.api import load_bundle, DEFAULT_BUNDLE_PATH, PREVIOUS_BUNDLE_PATH
from altanalyze3.components.rna2lipid.release import release_manifest, file_sha256, exact_identifiers
from altanalyze3.components.cellHarmony.flask.pipeline import _build_imputed_lipid_adata, _lookup_reference
from altanalyze3.components.cellHarmony.flask.disk_predictions import PredictionWriter
from altanalyze3.components.cellHarmony.flask.pipeline import _read_differential_h5ad
from altanalyze3.components.cellHarmony.cellHarmony_differential import (
    _moderated_t_test, _rank_genes_scanpy, _compute_pseudobulk_log2fc)
from altanalyze3.components.cellHarmony.imputed_pseudobulk import aggregate_imputed_predictions

HERE = Path(__file__).resolve().parents[1]


class NativeLog2DeploymentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = release_manifest()
        cls.bundle = load_bundle()
        cls.training = pd.read_csv(HERE / "artifacts/LungMAP_full202_native_log2_20261005/candidate_training_RNA.csv", index_col=0)

    def lipid_data(self):
        rng = np.random.default_rng(516)
        values = rng.normal(0, .25, size=(12, 202)) - 5
        values[:6, :4] += 1
        obs = pd.DataFrame({"condition": pd.Categorical(["case"] * 6 + ["control"] * 6),
                            "state": "macrophage", "sample": [f"s{i // 2}" for i in range(12)]},
                           index=[f"cell{i}" for i in range(12)])
        return ad.AnnData(values, obs=obs, var=pd.DataFrame(index=self.manifest["Y_columns"]),
                          uns={"expression_scale": "native_relative_log2", "modality": "lipids"})

    def test_release_is_exact_validated_model_and_preserves_rollback(self):
        self.assertEqual(file_sha256(DEFAULT_BUNDLE_PATH), self.manifest["bundle_sha256"])
        self.assertEqual(file_sha256(PREVIOUS_BUNDLE_PATH), self.manifest["previous_bundle_sha256"])
        self.assertEqual((len(self.bundle.input_genes), len(self.bundle.output_lipids)), (1303, 202))
        self.assertFalse(self.bundle.metadata["candidate_only"])
        reference = load_bundle(HERE / "artifacts/LungMAP_full202_native_log2_20261005/candidate_bundle.pkl")
        actual = self.bundle.predict_from_dataframe(self.training).predictions
        expected = reference.predict_from_dataframe(self.training).predictions
        np.testing.assert_array_equal(actual, expected)
        self.assertTrue((actual.to_numpy() < 0).any())

    def test_all_stored_estimators_match_optimized_inference(self):
        original = self.bundle._lipidwise_linear_parameters
        try:
            optimized = self.bundle.predict_from_dataframe(self.training).predictions
            self.bundle._lipidwise_linear_parameters = lambda: None
            individual = self.bundle.predict_from_dataframe(self.training).predictions
            np.testing.assert_allclose(optimized, individual, atol=1e-12, rtol=1e-12)
        finally:
            self.bundle._lipidwise_linear_parameters = original

    def test_all_five_lung_references_use_the_promoted_default(self):
        registry = HERE.parent / "cellHarmony/flask/reference_config.json"
        refs = json.loads(registry.read_text())["species"][0]["references"]
        found = []
        for reference in refs:
            if "lipids" in reference.get("impute_config", {}):
                entry = _lookup_reference("human", reference["id"], registry)
                self.assertEqual(Path(entry["impute_config"]["lipids"]["bundle_path"]), DEFAULT_BUNDLE_PATH)
                found.append(reference["id"])
        self.assertEqual(len(found), 5)

    def test_reference_cannot_relabel_a_legacy_bundle_as_native_log2(self):
        legacy = SimpleNamespace(metadata={})
        with patch("altanalyze3.components.cellHarmony.flask.pipeline.load_rna2lipid_bundle", return_value=legacy):
            with self.assertRaisesRegex(ValueError, "target encoding disagree"):
                _build_imputed_lipid_adata(ad.AnnData(np.ones((2, 2))),
                    {"impute_config": {"lipids": {"expression_scale": "native_relative_log2"}}})

    def test_export_retains_negatives_and_uses_exact_exp2(self):
        expression = self.training.iloc[:8]
        query = ad.AnnData(expression.to_numpy(), obs=pd.DataFrame(index=expression.index),
                           var=pd.DataFrame(index=expression.columns), uns={"log1p": {"base": 2}})
        values, summary = _build_imputed_lipid_adata(query)
        expected = self.bundle.predict_from_dataframe(expression).predictions.to_numpy(dtype=np.float32)
        np.testing.assert_array_equal(values.X, expected)
        np.testing.assert_array_equal(values.layers["counts"], np.exp2(expected.astype(np.float64)).astype(np.float32))
        self.assertEqual(summary["clipped_negative_values"], 0)
        self.assertGreater(summary["negative_log2_values_preserved"], 0)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "lipids.h5ad"
            values.write_h5ad(path)
            saved = ad.read_h5ad(path)
            np.testing.assert_array_equal(saved.X, expected)
            self.assertEqual(saved.uns["expression_scale"], "native_relative_log2")

    def test_declared_natural_log_RNA_matches_log2_RNA_without_mutation(self):
        expression = self.training.iloc[:8]
        logged = expression.to_numpy() * np.log(2)
        query = ad.AnnData(sparse.csr_matrix(logged), obs=pd.DataFrame(index=expression.index),
                           var=pd.DataFrame(index=expression.columns), uns={"log1p": {"base": None}})
        before = query.X.copy()
        values, summary = _build_imputed_lipid_adata(query)
        expected = self.bundle.predict_from_dataframe(expression).predictions.to_numpy(dtype=np.float32)
        np.testing.assert_allclose(values.X, expected, atol=1e-6, rtol=1e-6)
        self.assertEqual((query.X != before).nnz, 0)
        self.assertAlmostEqual(summary["RNA_log_base_conversion_factor"], 1 / np.log(2))

    def test_negative_targets_have_full_panel_BH_for_both_tests(self):
        values = self.lipid_data()
        moderated, tested = _moderated_t_test(values, "condition", "case", "control", "macrophage")
        self.assertEqual(tested, 202)
        np.testing.assert_allclose(moderated.fdr, multipletests(moderated.pval, method="fdr_bh")[1])
        names, fdr, fold, pval = _rank_genes_scanpy(values, "condition", "case", "control", "wilcoxon")
        self.assertEqual(len(names), 202)
        np.testing.assert_allclose(fdr, multipletests(pval, method="fdr_bh")[1])
        expected = _compute_pseudobulk_log2fc(values, "condition", "case", "control")
        np.testing.assert_allclose(fold, expected.log2fc.reindex(names))

    def test_raw_tests_match_unmodified_reference_and_BH_is_baseline_invariant(self):
        values = self.lipid_data()
        unmarked = values.copy()
        unmarked.uns = {"log1p": {"base": 2}}
        # Unmarked negative values fail historical independent filtering, so
        # compare raw Wilcoxon tests, and independently reconstruct the moderated t.
        _, _, _, raw_old = _rank_genes_scanpy(unmarked, "condition", "case", "control", "wilcoxon")
        _, fdr, folds, raw_new = _rank_genes_scanpy(values, "condition", "case", "control", "wilcoxon")
        np.testing.assert_array_equal(raw_old.sort_index(), raw_new.sort_index())
        shifted = values.copy(); shifted.X += 20
        _, other_fdr, other_fold, other_raw = _rank_genes_scanpy(shifted, "condition", "case", "control", "wilcoxon")
        np.testing.assert_array_equal(fdr.sort_index(), other_fdr.sort_index())
        np.testing.assert_array_equal(raw_new.sort_index(), other_raw.sort_index())
        np.testing.assert_allclose(folds.sort_index(), other_fold.sort_index(), atol=1e-12)
        actual, _ = _moderated_t_test(values, "condition", "case", "control", "macrophage")
        shifted_actual, _ = _moderated_t_test(shifted, "condition", "case", "control", "macrophage")
        np.testing.assert_allclose(actual.pval, shifted_actual.pval, atol=1e-12)
        np.testing.assert_allclose(actual.fdr, shifted_actual.fdr, atol=1e-12)

    def test_disk_backed_tests_match_in_memory(self):
        values = self.lipid_data()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "source.h5ad"; values.write_h5ad(path)
            disk = _read_differential_h5ad(path, disk_backed=True)
            try:
                for fn in [_rank_genes_scanpy, _moderated_t_test]:
                    args = ("condition", "case", "control", "wilcoxon" if fn == _rank_genes_scanpy else "macrophage")
                    a, b = fn(values, *args), fn(disk, *args)
                    if fn == _rank_genes_scanpy:
                        for x, y in zip(a[1:], b[1:]):
                            np.testing.assert_allclose(x.sort_index(), y.sort_index(), atol=1e-12)
                    else:
                        pd.testing.assert_frame_equal(a[0], b[0], atol=1e-12)
            finally:
                disk._analysis_h5_handle.close()

    def test_missing_or_wrong_lipid_identity_blocks_DE(self):
        for values in [self.lipid_data()[:, :-1].copy(), self.lipid_data()]:
            if values.n_vars == 202:
                values.var_names = list(values.var_names[:-1]) + ["wrong_lipid"]
            with self.assertRaisesRegex(ValueError, "complete approved release"):
                _moderated_t_test(values, "condition", "case", "control", "macrophage")

    def test_common_aggregation_is_mean_linear_not_mean_log_or_panel_normalization(self):
        for modality, scale, x in [("lipids", "native_relative_log2", [-3., 1.]),
                                   ("adt", "log1p", [0., np.log(4.)]),
                                   ("grn", "linear", [-2., 4.]), ("grn_tf", "linear", [-2., 4.])]:
            values = ad.AnnData(np.array(x)[:, None],
                obs=pd.DataFrame({"state": "s", "sample": "donor", "condition": "case"}, index=["c1", "c2"]),
                var=pd.DataFrame(index=["feature"]),
                uns={"modality": modality, "expression_scale": scale,
                     **({"log1p": {"base": None}} if scale == "log1p" else {})})
            out, audit = aggregate_imputed_predictions(values, population_col="state", sample_col="sample", covariate_col="condition", min_cells=2)
            expected = np.mean(x) if scale == "linear" else (np.exp2(x).mean() if scale == "native_relative_log2" else np.expm1(x).mean())
            np.testing.assert_allclose(out.layers["counts"], [[expected]], rtol=1e-7)
            self.assertFalse(audit["model_recomputed"])
            self.assertFalse(audit["feature_panel_total_normalized"])
            self.assertEqual(audit["aggregated_cells"], 2)
            self.assertEqual(list(out.var_names), list(values.var_names))

    def test_ambiguous_target_scale_and_conflicting_samples_are_not_guessed(self):
        values = self.lipid_data()
        values.uns.update(modality="metabolite", expression_scale="log2")
        with self.assertRaisesRegex(ValueError, "preprocessing protocol"):
            aggregate_imputed_predictions(values, population_col="state", sample_col="sample", covariate_col="condition", min_cells=2)
        values = self.lipid_data(); values.obs.loc["cell0", "condition"] = "control"
        with self.assertRaisesRegex(ValueError, "conflicting"):
            aggregate_imputed_predictions(values, population_col="state", sample_col="sample", covariate_col="condition", min_cells=2)

    def test_streamed_native_export_matches_memory_and_rejects_nonfinite(self):
        values = self.lipid_data()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "streamed.h5ad"
            writer = PredictionWriter(path, values.obs, values.var, dict(values.uns), {}, expression_scale="native_relative_log2")
            for start in range(0, len(values), 3):
                writer.append(values.X[start:start + 3])
            writer.close()
            saved = ad.read_h5ad(path)
            np.testing.assert_array_equal(saved.X, values.X.astype(np.float32))
            np.testing.assert_array_equal(saved.layers["counts"], np.exp2(saved.X.astype(np.float64)).astype(np.float32))
            bad = PredictionWriter(Path(folder) / "bad.h5ad", values.obs, values.var, {}, {}, expression_scale="native_relative_log2")
            try:
                with self.assertRaisesRegex(ValueError, "Nonfinite"):
                    bad.append(np.full((3, 202), np.nan))
            finally:
                bad.abort()


if __name__ == "__main__":
    unittest.main()
