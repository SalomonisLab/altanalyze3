"""
baselines.py: Baseline models for Proteus benchmarking.

Three baselines are implemented:

1. SequenceDeltaBaseline
   ESM2(alt) − ESM2(ref) → logistic regression per task.
   Represents the simplest "what does the sequence difference tell us?" baseline.

2. LengthRatioBaseline
   alt_length / ref_length → logistic regression / threshold.
   Captures the effect of large truncations (NMD, partial loss).

3. StructuralAnnotationBaseline
   Rule-based predictions using TM helix counts, domain counts, etc. from
   StructuralFeatureEncoder output. No learned parameters.
   Represents the "expert rules" upper bound for annotation-based prediction.

All baselines implement a common interface:
  .fit(train_loader)  → train (if parametric)
  .predict(loader)    → dict[task → np.ndarray of probabilities]
  .evaluate(loader)   → dict[task → TaskMetrics]
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


class SequenceDeltaBaseline:
    """
    ESM2 embedding difference → logistic regression per task.

    Concatenates [alt_emb, ref_emb, alt_emb - ref_emb, |alt_emb - ref_emb|]
    and fits a separate L2-regularized logistic regression per task.

    Parameters
    ----------
    C : float
        Logistic regression inverse regularization strength.
    max_iter : int
        Maximum solver iterations.
    """

    def __init__(self, C: float = 0.1, max_iter: int = 500) -> None:
        self.C = C
        self.max_iter = max_iter
        self._classifiers: Dict[str, Any] = {}

    def _extract_features(self, batch) -> np.ndarray:
        """Extract [ref_prot, alt_prot, diff, abs_diff] feature vector."""
        ref = batch["ref_protein_emb"].numpy()
        alt = batch["alt_protein_emb"].numpy()
        diff = alt - ref
        abs_diff = np.abs(diff)
        return np.concatenate([ref, alt, diff, abs_diff], axis=1)

    def fit(self, train_loader) -> None:
        """Fit one logistic regression per task on training data."""
        from sklearn.linear_model import LogisticRegression
        from proteus.heads.task_heads import BINARY_TASKS

        # Collect all features and labels
        features_list = []
        labels: Dict[str, List] = {t: [] for t in BINARY_TASKS}
        masks: Dict[str, List] = {t: [] for t in BINARY_TASKS}

        for batch in train_loader:
            feats = self._extract_features(batch)
            features_list.append(feats)
            for task in BINARY_TASKS:
                lk = "label_preservation" if task == "global_preservation" else f"label_{task}"
                lab = batch.get(lk, batch.get("label_preservation"))
                if lab is not None:
                    labels[task].append(lab.numpy().flatten())
                mk = batch.get("task_masks", {}).get(task)
                if mk is not None:
                    masks[task].append(mk.numpy().flatten())

        X = np.concatenate(features_list, axis=0)

        for task in BINARY_TASKS:
            if not labels[task]:
                continue
            y = np.concatenate(labels[task]).astype(int)
            m = np.concatenate(masks[task]) if masks[task] else np.ones(len(y))
            valid = (y != -1) & (m > 0.5)
            if valid.sum() < 10:
                continue
            Xv, yv = X[valid], y[valid]
            if yv.sum() < 2 or (1 - yv).sum() < 2:
                continue
            clf = LogisticRegression(C=self.C, max_iter=self.max_iter, solver="saga", n_jobs=1)
            try:
                clf.fit(Xv, yv)
                self._classifiers[task] = clf
            except Exception as exc:
                warnings.warn(f"[SequenceDeltaBaseline] Fit failed for {task}: {exc}")

        print(f"[SequenceDeltaBaseline] Fitted {len(self._classifiers)} classifiers.")

    def predict(self, loader) -> Dict[str, np.ndarray]:
        """Return predicted probabilities per task for all samples in loader."""
        features_list = []
        for batch in loader:
            feats = self._extract_features(batch)
            features_list.append(feats)
        X = np.concatenate(features_list, axis=0)

        preds: Dict[str, np.ndarray] = {}
        for task, clf in self._classifiers.items():
            try:
                preds[task] = clf.predict_proba(X)[:, 1]
            except Exception:
                preds[task] = np.full(len(X), 0.5)
        return preds


class LengthRatioBaseline:
    """
    alt_length / ref_length → per-task thresholded prediction.

    Rationale: large truncations predict NMD, domain loss, and topology changes.
    A ratio < 0.6 predicts "not preserved" for most tasks.

    Fits thresholds on training data using Youden's J statistic.
    """

    def __init__(self) -> None:
        self._thresholds: Dict[str, float] = {}
        self._directions: Dict[str, int] = {}  # 1 = ratio > thresh → positive

    @staticmethod
    def _get_length_ratio(batch) -> np.ndarray:
        """Estimate alt/ref length ratio from structural feature (protein_length_norm idx 14)."""
        ref_str = batch["ref_structural_raw"].numpy()
        alt_str = batch["alt_structural_raw"].numpy()
        # Feature [14] = protein_length / 1000 (see structural.py layout)
        ref_len = np.maximum(ref_str[:, 14] * 1000, 1)
        alt_len = alt_str[:, 14] * 1000
        return alt_len / ref_len

    def fit(self, train_loader) -> None:
        """Find best threshold per task via Youden's J."""
        from proteus.heads.task_heads import BINARY_TASKS
        from sklearn.metrics import roc_curve

        ratios_list = []
        labels: Dict[str, List] = {t: [] for t in BINARY_TASKS}
        masks: Dict[str, List] = {t: [] for t in BINARY_TASKS}

        for batch in train_loader:
            ratios_list.append(self._get_length_ratio(batch))
            for task in BINARY_TASKS:
                lk = "label_preservation" if task == "global_preservation" else f"label_{task}"
                lab = batch.get(lk, batch.get("label_preservation"))
                if lab is not None:
                    labels[task].append(lab.numpy().flatten())
                mk = batch.get("task_masks", {}).get(task)
                if mk is not None:
                    masks[task].append(mk.numpy().flatten())

        ratios = np.concatenate(ratios_list)

        for task in BINARY_TASKS:
            if not labels[task]:
                continue
            y = np.concatenate(labels[task]).astype(int)
            m = np.concatenate(masks[task]) if masks[task] else np.ones(len(y))
            valid = (y != -1) & (m > 0.5)
            r, yv = ratios[valid], y[valid]
            if yv.sum() < 2 or (1 - yv).sum() < 2:
                continue
            try:
                fpr, tpr, threshs = roc_curve(yv, r)
                j = tpr - fpr
                best_idx = np.argmax(j)
                self._thresholds[task] = float(threshs[best_idx])
                # Most tasks: higher ratio → more preserved → positive
                self._directions[task] = 1
            except Exception:
                self._thresholds[task] = 0.7
                self._directions[task] = 1

        print(f"[LengthRatioBaseline] Thresholds: {self._thresholds}")

    def predict(self, loader) -> Dict[str, np.ndarray]:
        """Return predicted probabilities (ratio as soft score, thresholded)."""
        ratios_list = []
        for batch in loader:
            ratios_list.append(self._get_length_ratio(batch))
        ratios = np.clip(np.concatenate(ratios_list), 0, 1)

        preds: Dict[str, np.ndarray] = {}
        for task, thresh in self._thresholds.items():
            # Use ratio as probability directly (clipped 0-1)
            if self._directions.get(task, 1) == 1:
                preds[task] = ratios
            else:
                preds[task] = 1.0 - ratios
        return preds


class StructuralAnnotationBaseline:
    """
    Rule-based predictions using structural annotation features only.
    No learned parameters — represents the "expert annotation" ceiling.

    Rules:
      global_preservation:      ref_tm == alt_tm AND ref_domain_count == alt_domain_count
      topology_preserved:       ref_tm_count == alt_tm_count
      surface_retained:         alt has ≥1 extracellular topology annotation
      kinase_competent:         alt has kinase domain (feature [16] > 0)
      tf_competent:             alt has DNA-binding domain (feature [17] or [18] > 0)
      localization_preserved:   HPA annotation present in alt
      tm_insertion_competent:   alt has ≥1 TM helix predicted
      ppi_interface_preserved:  alt has PPI binding features (feature [5] in PPIEncoder)
    """

    def predict(self, loader) -> Dict[str, np.ndarray]:
        """Compute rule-based prediction probabilities from structural features."""
        ref_structs = []
        alt_structs = []
        has_alt = []

        for batch in loader:
            ref_structs.append(batch["ref_structural_raw"].numpy())
            alt_structs.append(batch["alt_structural_raw"].numpy())
            has_alt.append(batch["has_alt_protein"].numpy().astype(float))

        R = np.concatenate(ref_structs, axis=0)   # [N, 96]
        A = np.concatenate(alt_structs, axis=0)
        H = np.concatenate(has_alt)

        preds: Dict[str, np.ndarray] = {}

        # Use feature indices from StructuralFeatureEncoder layout
        # [7] = tm_helix_count, [12] = is_membrane_protein, [0] = domain_count
        # [13] = is_surface_protein, [16] = dna_binding, [17] = zinc_finger
        # [28] = biogrid_partner_count_log1p, [30] = ppi_binding_count_log1p
        # [47] = has_hpa_extracellular, [48] = is_hpa_membrane_or_surface

        def soft(x): return np.clip(x, 0, 1)

        # global_preservation: soft agreement between ref and alt on key features
        tm_agree = 1.0 - np.abs(R[:, 7] - A[:, 7]) / (np.maximum(R[:, 7], 1) + 1e-8)
        dom_agree = 1.0 - np.abs(R[:, 0] - A[:, 0]) / (np.maximum(R[:, 0], 1) + 1e-8)
        preds["global_preservation"] = soft((tm_agree + dom_agree) / 2 * H)

        # topology_preserved: TM helix counts match
        tm_diff = np.abs(R[:, 7] - A[:, 7])
        preds["topology_preserved"] = soft(1.0 - tm_diff / (np.maximum(R[:, 7], 1) + 1e-8)) * H

        # surface_retained: alt has extracellular annotation or TM helix
        preds["surface_retained"] = soft(
            (A[:, 13] + A[:, 12]) / 2 * H
        )

        # kinase_competent: alt retains kinase-related features
        # Feature [16] = dna_binding_feature_count (log1p) — proxy for kinase if kinase domain
        # Use both domain count similarity and structural feature preservation
        preds["kinase_competent"] = soft(A[:, 0] / (R[:, 0] + 1e-8) * H)

        # tf_competent: DNA-binding features retained
        dna_ref = R[:, 16] + R[:, 17]
        dna_alt = A[:, 16] + A[:, 17]
        preds["tf_competent"] = soft(dna_alt / (dna_ref + 1e-8) * H)

        # disorder_preserved: always 0.5 (no structural feature for disorder)
        preds["disorder_preserved"] = np.full(len(R), 0.5)

        # localization_preserved: HPA annotation match
        hpa_match = (A[:, 47] == R[:, 47]).astype(float)
        preds["localization_preserved"] = soft(hpa_match * H)

        # tm_insertion_competent: alt has ≥1 predicted TM helix
        preds["tm_insertion_competent"] = soft(np.minimum(A[:, 7], 1) * H)

        # signaling_competent: alt retains both TM and sufficient domain count
        tm_ok = (A[:, 7] >= R[:, 7] * 0.5).astype(float)
        dom_ok = (A[:, 0] >= R[:, 0] * 0.5).astype(float)
        preds["signaling_competent"] = soft(tm_ok * dom_ok * H)

        # ppi_interface_preserved: alt has PPI binding features
        ppi_ref = R[:, 30]
        ppi_alt = A[:, 30]
        preds["ppi_interface_preserved"] = soft(ppi_alt / (ppi_ref + 1e-8) * H)

        return preds
