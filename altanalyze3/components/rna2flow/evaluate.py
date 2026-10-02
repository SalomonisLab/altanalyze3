"""Scoring and adversarial controls for a CITE-seq -> flow transfer.

A high ARI against FlowSOM proves nothing on its own: a method that assigns one label to every
event, or that is learning batch rather than biology, can still look respectable. Every run
here therefore carries its controls.
"""
from __future__ import annotations

import numpy as np
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

__all__ = ["score_transfer", "adversarial_controls"]


def score_transfer(predicted, flowsom):
    p = np.asarray(predicted).astype(str)
    f = np.asarray(flowsom).astype(str)
    _, counts = np.unique(p, return_counts=True)
    return {
        "ARI": float(adjusted_rand_score(f, p)),
        "AMI": float(adjusted_mutual_info_score(f, p)),
        "n_labels_used": int(len(counts)),
        "largest_label_fraction": float(counts.max() / counts.sum()),
        "entropy_bits": float(-(counts / counts.sum() * np.log2(counts / counts.sum())).sum()),
    }


def adversarial_controls(fn, cite, flow, labels, flowsom, seed=0, n_dropout=3):
    """Four questions asked of every method, each one able to invalidate its headline ARI.

    shuffled_labels   train on permuted labels. ARI must collapse toward 0. If it does not,
                      the score is coming from cluster-size structure, not from the markers.
    held_out_cite     split CITE in two, transfer one half onto the other. This is the method's
                      ceiling when there is no platform gap at all.
    marker_dropout    drop one shared marker at a time. A score that depends on a single
                      channel is not a transfer, it is a gate.
    label_collapse    the fraction of flow events taking the single most common label.
    """
    rng = np.random.default_rng(seed)
    out = {}

    perm = rng.permutation(len(labels))
    out["shuffled_labels_ARI"] = score_transfer(fn(cite, flow, np.asarray(labels)[perm]), flowsom)["ARI"]

    idx = rng.permutation(len(cite))
    a, b = idx[: len(idx) // 2], idx[len(idx) // 2:]
    pred_b = fn(cite[a], cite[b], np.asarray(labels)[a])
    out["held_out_cite_accuracy"] = float((pred_b == np.asarray(labels)[b]).mean())

    drops = []
    for j in rng.choice(cite.shape[1], min(n_dropout, cite.shape[1]), replace=False):
        keep = [c for c in range(cite.shape[1]) if c != j]
        drops.append(score_transfer(fn(cite[:, keep], flow[:, keep], labels), flowsom)["ARI"])
    out["marker_dropout_ARI_min"] = float(min(drops))
    out["marker_dropout_ARI_max"] = float(max(drops))
    return out
