"""Inference API for the ``rna2psi`` component: impute splicing-event PSI from bulk RNA
counts, mirroring ``components.rna2metabolite.api`` so the cellHarmony-web pipeline can
register a PSI-imputation modality on the same plumbing.

Default model (v2, 2026-07-09): one per-event bootstrap-LASSO stability-selected Ridge
(<=5 predictor genes/event from all ~39,771 genes), over 7,630 mutation-signature events,
trained on all 437 Leucegene AML samples (AltAnalyze2 MultiPath-PSI). Feature genes are
matched to the input by Ensembl gene ID. See README.md / METHODS.md for full detail.
The prior v1 model (protein-coding, ElasticNet/Ridge, 18,168 events) is in ``legacy/``.
"""
from __future__ import annotations

from pathlib import Path

from ._impute import PerEventPSIBundle, PredictionResult, cp10k_log1p  # noqa: F401

DEFAULT_BUNDLE_PATH = Path(__file__).with_name("artifacts") / "rna2psi_leucegene_focused.pkl.gz"
RECLASSIFIER_PATH = Path(__file__).with_name("artifacts") / "rna2psi_reclassifiers.pkl.gz"


class Rna2PsiBundle(PerEventPSIBundle):
    """AML RNA -> splicing PSI imputation bundle."""


def load_bundle(bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> Rna2PsiBundle:
    return Rna2PsiBundle.load(bundle_path)
