"""Inference API for the AML rna2metabolite component: impute metabolite abundance from
bulk/pseudobulk RNA, mirroring ``components.rna2grn.api`` / ``components.rna2adt.api`` so
the cellHarmony-web pipeline can register a metabolite-imputation modality on the same
plumbing.

Model: one per-metabolite ridge (1000 genes/target, α=100) trained on the CPTAC AML
cohort (84 RNA+metabolomics cases). This is an AML-specific classifier and is independent
of the (non-AML) lung ``rna2lipid`` models.
"""
from __future__ import annotations

from pathlib import Path

from ._impute import PerTargetImputeBundle, PredictionResult, cp10k_log1p  # noqa: F401

DEFAULT_BUNDLE_PATH = Path(__file__).with_name("artifacts") / "rna2metabolite_aml_bundle.pkl.gz"


class Rna2MetaboliteBundle(PerTargetImputeBundle):
    """AML RNA -> metabolite imputation bundle."""


def load_bundle(bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> Rna2MetaboliteBundle:
    return Rna2MetaboliteBundle.load(bundle_path)
