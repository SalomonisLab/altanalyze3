"""Inference API for the AML rna2lipid component: impute lipid abundance from
bulk/pseudobulk RNA.

This is a SEPARATE, AML-specific model and code path that does not touch the lung/IPF
``components.rna2lipid`` multitask elastic-net model. Model: one per-lipid ridge
(120 genes/target, α=100) trained on the CPTAC AML cohort (87 RNA+lipidomics cases).
"""
from __future__ import annotations

from pathlib import Path

from ._impute import PerTargetImputeBundle, PredictionResult, cp10k_log1p  # noqa: F401

DEFAULT_BUNDLE_PATH = Path(__file__).with_name("artifacts") / "rna2lipid_aml_bundle.pkl.gz"


class Rna2LipidAmlBundle(PerTargetImputeBundle):
    """AML RNA -> lipid imputation bundle."""


def load_bundle(bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> Rna2LipidAmlBundle:
    return Rna2LipidAmlBundle.load(bundle_path)
