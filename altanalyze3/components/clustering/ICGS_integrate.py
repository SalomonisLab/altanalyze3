"""Integrate the cluster results of several independent ICGS3 runs.

Each ICGS3 run reduces one dataset to clusters, markers and centroid folds. This module compares
those clusters ACROSS datasets and merges the ones that describe the same state, without moving a
single cell and without a batch-correction step.

Why the comparison acts on clusters rather than cells: an embedding method must decide which
distance is batch and which is biology. Clusters already carry markers, so a cluster-to-cluster
comparison never needs that decision.

Memory: the whole comparison runs on the centroid matrices, which hold genes x clusters, a few
thousand by a few dozen. Cells are read only for the optional rebuild, in row blocks, one dataset
at a time.

Usage:
    python -m altanalyze3.components.clustering.ICGS_integrate \\
        --run Adams=/path/ICGS3_Adams_corr03_ne3 \\
        --run Basil2022=/path/ICGS3_Basil2022 \\
        --output-dir /path/ICGS3_integration \\
        --cells-per-cluster 100
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy.stats import hypergeom


def _log(message: str) -> None:
    print(f"[ICGS-integrate] {message}", flush=True)


# ----------------------------------------------------------------------------------- inputs
@dataclass
class RunInputs:
    """One completed ICGS3 run, read from the files it already wrote."""

    name: str
    path: str
    centroids: pd.DataFrame            # genes x clusters, folds
    markers: Dict[str, List[str]]      # cluster -> ordered marker list
    marker_r: Dict[str, Dict[str, float]]
    names: Dict[str, str]              # cluster -> GO-Elite cell_type_prediction
    name_fdr: Dict[str, float]
    cells: pd.DataFrame                # barcode, cluster, svm score
    cluster_sizes: Dict[str, int] = field(default_factory=dict)
    unique_markers: Optional[pd.DataFrame] = None   # marker, top_cluster, pearson_r
    annotations: Optional[pd.DataFrame] = None      # barcode x reference annotation columns

    @property
    def clusters(self) -> List[str]:
        return list(self.centroids.columns)


def unique_marker_table(path: str) -> pd.DataFrame:
    """Per-cluster unique markers with correlation, from the ICGS3 all-correlations file.

    Each gene appears exactly once in that file: ICGS3 assigns every gene to the single cluster
    it correlates with best. Genes shared by a whole lineage are therefore already competed away
    within a study, which is what makes an overlap between two studies informative. Comparing
    recomputed per-gene correlations instead reintroduces those lineage genes and made the median
    cluster share 30 of its top 50 genes with some reference state.
    """
    f = os.path.join(path, "MarkerFinder", "icgs3_markers_all_correlations.tsv")
    if not os.path.exists(f):
        f = os.path.join(path, "MarkerFinder", "icgs3_markers.tsv")
    d = pd.read_csv(f, sep="\t")
    d["top_cluster"] = d["top_cluster"].astype(str)
    d["marker"] = d["marker"].astype(str)
    return d[["marker", "top_cluster", "pearson_r"]]


def load_run(name: str, path: str, marker_top_n: int = 60) -> RunInputs:
    mf = os.path.join(path, "MarkerFinder")
    cen_path = os.path.join(mf, "icgs3_marker_heatmap_fold_matrix.centroids.tsv")
    mk_path = os.path.join(mf, "icgs3_markers.tsv")
    cell_path = os.path.join(path, "icgs3_cell_barcode_clusters.tsv")
    go_path = os.path.join(path, "GO-Elite", "icgs3_cell_state_predictions.tsv")
    for required in (cen_path, mk_path, cell_path):
        if not os.path.exists(required):
            raise FileNotFoundError(f"{name}: missing {required}")

    centroids = pd.read_csv(cen_path, sep="\t", index_col=0)
    centroids = centroids.astype(np.float32)
    centroids = centroids.loc[~centroids.index.duplicated(keep="first")]

    mk = pd.read_csv(mk_path, sep="\t")
    markers: Dict[str, List[str]] = {}
    marker_r: Dict[str, Dict[str, float]] = {}
    for cluster, grp in mk.groupby("top_cluster"):
        grp = grp.sort_values("pearson_r", ascending=False)
        markers[str(cluster)] = grp["marker"].astype(str).head(marker_top_n).tolist()
        marker_r[str(cluster)] = dict(zip(grp["marker"].astype(str), grp["pearson_r"].astype(float)))

    names: Dict[str, str] = {}
    name_fdr: Dict[str, float] = {}
    if os.path.exists(go_path):
        go = pd.read_csv(go_path, sep="\t")
        if {"cluster", "cell_type_prediction"}.issubset(go.columns):
            go = go.sort_values("fdr")
            for _, row in go.iterrows():
                c = str(row["cluster"])
                if c not in names:
                    names[c] = str(row["cell_type_prediction"])
                    name_fdr[c] = float(row.get("fdr", np.nan))

    full = pd.read_csv(cell_path, sep="\t", low_memory=False)
    usecols = ["barcode", "ICGS3_cluster", "ICGS3_SVM_score"]
    cells = full[[c for c in usecols if c in full.columns]].copy()
    cells["ICGS3_cluster"] = cells["ICGS3_cluster"].astype(str)
    sizes = cells["ICGS3_cluster"].value_counts().to_dict()

    # Every remaining column is a candidate reference annotation. The technical columns ICGS3
    # writes are named here explicitly, so a new reference column is picked up without a code
    # change. Library and sample are experimental metadata, not cell types.
    technical = {"barcode", "Library", "sample", "ICGS3_cluster", "ICGS3_original_NMF_cluster",
                 "ICGS3_SVM_score", "ICGS3_SVM_margin", "ICGS3_cell_state_prediction"}
    ann_cols = [c for c in full.columns if c not in technical]
    annotations = None
    if ann_cols:
        annotations = full[["barcode"] + ann_cols].copy()
        annotations["barcode"] = annotations["barcode"].astype(str)

    unique_tbl = unique_marker_table(path)
    run = RunInputs(name=name, path=path, centroids=centroids, markers=markers,
                    marker_r=marker_r, names=names, name_fdr=name_fdr, cells=cells,
                    cluster_sizes=sizes, unique_markers=unique_tbl, annotations=annotations)
    _log(f"{name}: {centroids.shape[1]} clusters, {centroids.shape[0]} centroid genes, "
         f"{cells.shape[0]} cells, {len(names)} GO-Elite names, "
         f"reference annotations: {ann_cols if ann_cols else 'none'}")
    return run


# --------------------------------------------------------------------------- representative cells
def select_representative_cells(run: RunInputs, n_per_cluster: int) -> pd.DataFrame:
    """Top n cells per cluster by SVM score. A high score sits far from every other centroid."""
    cells = run.cells
    if "ICGS3_SVM_score" in cells.columns:
        ordered = cells.sort_values("ICGS3_SVM_score", ascending=False)
    else:
        ordered = cells
    picked = ordered.groupby("ICGS3_cluster", sort=False).head(int(n_per_cluster)).copy()
    picked["dataset"] = run.name
    short = picked["ICGS3_cluster"].value_counts()
    thin = [c for c, n in short.items() if n < int(n_per_cluster)]
    if thin:
        _log(f"{run.name}: {len(thin)} of {len(run.clusters)} clusters hold fewer than "
             f"{n_per_cluster} cells (smallest {int(short.min())})")
    return picked


# ------------------------------------------------------------------------------------- scoring
def _bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


@dataclass
class HState:
    """One harmonized state, grown across datasets."""
    sid: int
    members: List[str] = field(default_factory=list)          # "dataset|cluster"
    centroid: Optional[pd.Series] = None                      # genes -> mean fold
    marker_weight: Dict[str, float] = field(default_factory=dict)   # gene -> mean 1/rank
    n_contrib: int = 0
    support_datasets: List[str] = field(default_factory=list)   # datasets redundant with it

    @property
    def datasets(self) -> set:
        return {m.split("|")[0] for m in self.members}

    def ranked_markers(self, top_n: int) -> List[str]:
        return [g for g, _ in sorted(self.marker_weight.items(), key=lambda kv: -kv[1])][:top_n]

    def absorb(self, member: str, centroid: pd.Series, markers: Sequence[str]) -> None:
        self.members.append(member)
        w = {g: 1.0 / (i + 1) for i, g in enumerate(markers)}
        n = self.n_contrib
        for g in set(self.marker_weight) | set(w):
            # running mean of 1/rank across members: a gene ranked high in every dataset that
            # contributed keeps a high weight, a gene seen once is diluted
            self.marker_weight[g] = (self.marker_weight.get(g, 0.0) * n + w.get(g, 0.0)) / (n + 1)
        if self.centroid is None:
            self.centroid = centroid.astype(np.float32)
        else:
            joined = pd.concat([self.centroid, centroid.astype(np.float32)], axis=1)
            self.centroid = joined.mean(axis=1, skipna=True).astype(np.float32)
        self.n_contrib += 1


def enrichment_against_states(query_genes: Sequence[str], db: Dict[str, List[str]],
                              universe: int, *, overlap_fraction: float = 0.25,
                              min_overlap_floor: int = 3) -> pd.DataFrame:
    """Hypergeometric enrichment of a query marker set against each state's marker database.

    The genes each cluster owns vary from 1 to 449 across the 214 input clusters, median 41 at
    r > 0.25. A fixed minimum overlap is therefore incoherent: 10 genes cannot be reached by a
    cluster holding 11, and is trivial for one holding 449. The requirement is instead a fraction
    of the smaller of the two sets being compared, with an absolute floor so that very small sets
    still need real agreement. Expected overlap by chance is 0.31 genes for a 50-gene query
    against a median set, so the hypergeometric FDR remains the primary control.
    """
    q = set(query_genes)
    rows = []
    for state, genes_list in db.items():
        target = set(genes_list)
        k = len(q & target)
        if k == 0:
            continue
        p = hypergeom.sf(k - 1, universe, len(q), len(target))
        required = max(int(min_overlap_floor),
                       int(np.ceil(float(overlap_fraction) * min(len(q), len(target)))))
        rows.append({"state": state, "overlap": k, "query_size": len(q),
                     "db_size": len(target), "p_value": float(p),
                     "required_overlap": required, "meets_overlap": bool(k >= required)})
    df = pd.DataFrame(rows)
    if df.shape[0]:
        df["fdr"] = _bh_adjust(df["p_value"].to_numpy())
        df = df.sort_values(["overlap", "p_value"], ascending=[False, True])
    return df



def integrate_hierarchically(runs: List[RunInputs], *, marker_top_n: int, outdir: str,
                             pool: Optional[sp.csr_matrix] = None,
                             genes: Optional[pd.Index] = None,
                             barcodes: Optional[List[str]] = None,
                             node_of: Optional[List[str]] = None,
                             survival_rho: float = 0.3,
                             survival_min_markers: int = 1,
                             nomination_query_top: int = 60,
                             nomination_min_overlap: int = 10,
                             nomination_fdr: float = 0.05,
                             nomination_specificity: float = 0.0,
                             damage_floor: int = 1,
                             nomination_overlap_fraction: float = 0.17,
                             survival_ref_cells: int = 60,
                             seed_dataset: Optional[str] = None,
                             dataset_order: Optional[Sequence[str]] = None,
                             nomination_identity_top: int = 10,
                             nomination_identity_min: int = 5,
                             redundancy_top_n: int = 0,
                             redundancy_min_fraction: float = 0.30,
                             redundancy_centroid_r: float = 0.0,
                             compare_to_excluded: bool = False,
                             force_admit: Optional[set] = None,
                             suppress_cluster: Optional[set] = None,
                             legacy_damage_test: bool = False
                             ) -> Tuple[List[HState], pd.DataFrame, pd.DataFrame]:
    """Build a reference of non-redundant cell states across datasets.

    The dataset with the most clusters initialises the reference. Every cluster of every later
    dataset is tested once against that reference by gene set enrichment. A cluster whose marker
    genes are enriched in an existing state is redundant and is NOT included; the state it
    matched is left exactly as it was, with its own cells and its own markers. A cluster enriched
    in no existing state is added as a new state.

    Nothing is merged. A state therefore always holds the cells of the single dataset that first
    contributed it. Cross-study support is recorded separately, as the list of datasets whose
    clusters proved redundant with that state.

    Entry order decides which dataset is never tested and which must justify every cluster, so
    both overrides change the result. `dataset_order` names the full order and wins outright.
    `seed_dataset` names only the dataset that initialises the reference; the rest follow in
    descending cluster count. Passing neither keeps the default, descending cluster count.
    """
    by_name = {r.name: r for r in runs}
    names = resolve_dataset_order(list(by_name), {n: len(r.clusters) for n, r in by_name.items()},
                                  seed_dataset=seed_dataset, dataset_order=dataset_order)
    order = [by_name[n] for n in names]
    if dataset_order:
        why = "set explicitly by --dataset-order"
    elif seed_dataset:
        why = f"seed forced to {seed_dataset} by --seed-dataset, rest by descending cluster count"
    else:
        why = "most clusters first"
    _log(f"dataset order, {why}: " +
         ", ".join(f"{r.name}({len(r.clusters)})" for r in order))
    if len(order) > 1 and len(order[0].clusters) < max(len(r.clusters) for r in order[1:]):
        _log(f"NOTE seed {order[0].name} holds {len(order[0].clusters)} clusters, fewer than "
             f"{max(len(r.clusters) for r in order[1:])} in a later dataset. The seed is never "
             f"tested, so a coarser seed admits later clusters the finer partition would have "
             f"held.")

    states: List[HState] = []
    audit_rows: List[dict] = []
    support: Dict[int, List[str]] = {}
    next_sid = [0]

    def new_state() -> HState:
        next_sid[0] += 1
        return HState(sid=next_sid[0])

    # Mean expression of every input cluster over the shared gene space, taken from the same
    # pooled matrix the final centroids come from. Used for the centroid-correlation condition.
    node_centroid: Dict[str, np.ndarray] = {}
    if pool is not None and node_of is not None:
        idx_by_node: Dict[str, List[int]] = {}
        for i, n in enumerate(node_of):
            idx_by_node.setdefault(n, []).append(i)
        for n, ix in idx_by_node.items():
            node_centroid[n] = np.asarray(pool[ix].mean(axis=0)).ravel().astype(np.float64)
        _log(f"centroids over the shared gene space for {len(node_centroid)} input clusters")

    def _centroid_r(node: str, target_nodes: Sequence[str]) -> float:
        a = node_centroid.get(node)
        if a is None or not target_nodes:
            return float("nan")
        cols = [node_centroid[t] for t in target_nodes if t in node_centroid]
        if not cols:
            return float("nan")
        b = np.mean(np.vstack(cols), axis=0)
        ac, bc = a - a.mean(), b - b.mean()
        den = np.sqrt((ac ** 2).sum()) * np.sqrt((bc ** 2).sum())
        return float((ac * bc).sum() / den) if den > 0 else float("nan")

    evidence_rows: List[dict] = []
    excluded_nodes: List[str] = []          # dataset|cluster of every set-aside cluster
    force_admit = set(force_admit or ())
    suppress_cluster = set(suppress_cluster or ())
    if force_admit or suppress_cluster:
        _log(f"sub-cluster replacement active: {len(force_admit)} clusters forced in, "
             f"{len(suppress_cluster)} suppressed")
    both = force_admit & suppress_cluster
    if both:
        raise ValueError(
            f"{len(both)} cluster(s) appear in BOTH --force-admit and --suppress-cluster, so the "
            f"request contradicts itself: {sorted(both)[:8]}")
    all_nodes = {f"{r.name}|{c}" for r in order for c in r.clusters}
    unmatched = (force_admit | suppress_cluster) - all_nodes
    if unmatched:
        raise ValueError(
            f"{len(unmatched)} --force-admit/--suppress-cluster token(s) match no "
            f"dataset|cluster in this run, so the request would be silently ignored: "
            f"{sorted(unmatched)[:8]}")
    seed_forced = {n for n in force_admit if n.split("|", 1)[0] == order[0].name}
    if seed_forced:
        _log(f"NOTE --force-admit names {len(seed_forced)} cluster(s) of the seed dataset "
             f"{order[0].name}; every seed cluster already becomes a state, so those are no-ops")
    seed = order[0]
    for c in seed.clusters:
        if f"{seed.name}|{c}" in suppress_cluster:
            audit_rows.append({"step": seed.name, "cluster": c,
                               "action": "suppressed_replaced_by_finer_dataset",
                               "state": "", "overlap_r025": np.nan, "fdr_r025": np.nan,
                               "overlap_r030": np.nan, "fdr_r030": np.nan, "matched_state": ""})
            continue
        st = new_state()
        st.absorb(f"{seed.name}|{c}", seed.centroids[c], seed.markers.get(c, []))
        states.append(st)
        audit_rows.append({"step": seed.name, "cluster": c, "action": "seed_state",
                           "state": f"S{st.sid}", "overlap_r025": np.nan, "fdr_r025": np.nan,
                           "overlap_r030": np.nan, "fdr_r030": np.nan, "matched_state": ""})
        evidence_rows.append({"order": len(evidence_rows) + 1, "dataset": seed.name,
                              "cluster": c, "source": f"{seed.name}|{c}",
                              "decision": "retained", "state": f"S{st.sid}",
                              "matched_target": "", "target_kind": "", "target_name": "",
                              "overlap": np.nan, "required_overlap": np.nan, "fdr": np.nan,
                              "identity_shared": np.nan, "identity_required": np.nan,
                              "head_overlap": np.nan, "head_required": np.nan,
                              "centroid_r": np.nan, "centroid_r_required": np.nan,
                              "failed_test": "", "reason": "seed dataset, admitted without test"})
    for c in seed.clusters:
        if f"{seed.name}|{c}" in suppress_cluster:
            evidence_rows.append({"order": len(evidence_rows) + 1, "dataset": seed.name,
                                  "cluster": c, "source": f"{seed.name}|{c}",
                                  "decision": "excluded", "state": "", "matched_target": "",
                                  "target_kind": "", "target_name": "", "overlap": np.nan,
                                  "required_overlap": np.nan, "fdr": np.nan,
                                  "identity_shared": np.nan, "identity_required": np.nan,
                                  "head_overlap": np.nan, "head_required": np.nan,
                                  "centroid_r": np.nan, "centroid_r_required": np.nan,
                                  "failed_test": "qualification",
                                  "reason": "removed before the comparison; not used as a "
                                            "redundancy target"})
    _log(f"seed {seed.name}: {len(states)} states")

    for run in order[1:]:
        # Label every reference cell by its state, and every cell of the incoming dataset by its
        # own cluster. One MarkerFinder run over that combined set asks whether each incoming
        # cluster can hold unique markers while competing against the whole current reference.
        sid_of_node = {}
        for st in states:
            for m in st.members:
                sid_of_node[m] = f"S{st.sid}"
        cand_nodes = {f"{run.name}|{c}": c for c in run.clusters}
        labels_all = []
        for n in node_of:
            if n in sid_of_node:
                labels_all.append(sid_of_node[n])
            elif n in cand_nodes:
                labels_all.append(f"CAND|{cand_nodes[n]}")
            else:
                labels_all.append("")
        keep = [i for i, l in enumerate(labels_all) if l]
        sub_labels = [labels_all[i] for i in keep]

        # enrichment statistics are still recorded, as evidence about which reference state each
        # excluded cluster resembles, but they no longer decide inclusion
        # Gate 1 databases and queries both come from ICGS3's unique-marker assignment, so a
        # gene counted for one cluster is counted for no other within that study. Every state is
        # a single source cluster, because nothing is merged, so a state's marker set is simply
        # its source cluster's.
        by_run = {r.name: r for r in runs}
        db_a: Dict[str, List[str]] = {}
        db_b: Dict[str, List[str]] = {}
        gene_set = set(genes)
        for st in states:
            src = st.members[0]
            d_name, d_cl = src.split("|", 1)
            tbl = by_run[d_name].unique_markers
            sub = tbl.loc[tbl["top_cluster"] == d_cl]
            sub = sub.loc[sub["marker"].isin(gene_set)]
            db_a[f"S{st.sid}"] = list(sub.loc[sub["pearson_r"] > 0.25]
                                      .sort_values("pearson_r", ascending=False)
                                      .head(200)["marker"])
            db_b[f"S{st.sid}"] = list(sub.loc[sub["pearson_r"] > 0.30]
                                      .sort_values("pearson_r", ascending=False)
                                      .head(100)["marker"])
        target_nodes: Dict[str, List[str]] = {f"S{st.sid}": list(st.members) for st in states}
        target_name: Dict[str, str] = {f"S{st.sid}": st.members[0] for st in states}
        # Clusters already set aside stay in the comparison. A cluster redundant with an excluded
        # cluster describes the same population as something the reference already rejected, so
        # admitting it would reintroduce that population under a different dataset's label.
        if compare_to_excluded:
            for node in excluded_nodes:
                d_name, d_cl = node.split("|", 1)
                tbl = by_run[d_name].unique_markers
                sub = tbl.loc[(tbl["top_cluster"] == d_cl) & (tbl["marker"].isin(gene_set))]
                key = f"X|{node}"
                db_a[key] = list(sub.loc[sub["pearson_r"] > 0.25]
                                 .sort_values("pearson_r", ascending=False).head(200)["marker"])
                db_b[key] = list(sub.loc[sub["pearson_r"] > 0.30]
                                 .sort_values("pearson_r", ascending=False).head(100)["marker"])
                target_nodes[key] = [node]
                target_name[key] = node
            _log(f"{run.name}: comparison set holds {len(states)} retained states and "
                 f"{len(excluded_nodes)} excluded clusters")
        _sz = [len(v) for v in db_a.values()] or [0]
        _log(f"{run.name}: reference marker databases from unique markers, "
             f"median {int(np.median(_sz))} genes per state at r>0.25")
        cand_tbl = run.unique_markers
        universe = int(len(genes))

        # Gate 1, redundancy. A cluster whose markers are specifically enriched in one existing
        # reference state describes a population the reference already holds, and is excluded
        # before survival is considered. Both the query and the databases come from ICGS3's
        # unique-marker assignment, where each gene belongs to exactly one cluster within a
        # study, so genes shared by a whole lineage are already competed away.
        provisional: List[str] = []
        withdrawn: Dict[str, str] = {}
        enrich_hit: Dict[str, str] = {}
        def _ev(c, decision, best=None, db=None, shared_ident=np.nan, head_ov=np.nan,
                cr=np.nan, failed="", reason=""):
            key = str(best["state"]) if best is not None else ""
            evidence_rows.append({
                "order": len(evidence_rows) + 1, "dataset": run.name, "cluster": c,
                "source": f"{run.name}|{c}", "decision": decision, "state": "",
                "matched_target": key,
                "target_kind": ("excluded" if key.startswith("X|") else
                                ("state" if key else "")),
                "target_name": target_name.get(key, ""),
                "overlap": int(best["overlap"]) if best is not None else np.nan,
                "required_overlap": int(best["required_overlap"]) if best is not None else np.nan,
                "fdr": float(best["fdr"]) if best is not None else np.nan,
                "identity_shared": shared_ident,
                "identity_required": int(nomination_identity_min),
                "head_overlap": head_ov,
                "head_required": (int(np.ceil(float(redundancy_min_fraction)
                                              * int(redundancy_top_n)))
                                  if int(redundancy_top_n) > 0 else np.nan),
                "centroid_r": cr,
                "centroid_r_required": (float(redundancy_centroid_r)
                                        if float(redundancy_centroid_r) > 0 else np.nan),
                "failed_test": failed, "reason": reason})

        for c in run.clusters:
            if f"{run.name}|{c}" in suppress_cluster:
                continue
            sub = cand_tbl.loc[(cand_tbl["top_cluster"] == c)
                               & (cand_tbl["marker"].isin(gene_set))]
            q = list(sub.sort_values("pearson_r", ascending=False)
                     .head(int(nomination_query_top))["marker"])
            if not q:
                provisional.append(c)
                _ev(c, "retained", failed="no_query_markers",
                    reason="no unique markers inside the shared gene space")
                continue
            ea = enrichment_against_states(q, db_a, universe,
                                           overlap_fraction=nomination_overlap_fraction,
                                           min_overlap_floor=nomination_min_overlap)
            eb = enrichment_against_states(q, db_b, universe,
                                           overlap_fraction=nomination_overlap_fraction,
                                           min_overlap_floor=nomination_min_overlap)

            def _verdict(e):
                if e.shape[0] == 0:
                    return None, 0.0
                top = e.iloc[0]
                if not bool(top["meets_overlap"]) or float(top["fdr"]) > float(nomination_fdr):
                    return None, 0.0
                second = int(e.iloc[1]["overlap"]) if e.shape[0] > 1 else 0
                margin = (int(top["overlap"]) - second) / max(int(top["overlap"]), 1)
                return top, margin

            ta, ma = _verdict(ea)
            tb, mb = _verdict(eb)
            pa = ta is not None and ma >= float(nomination_specificity)
            pb = tb is not None and mb >= float(nomination_specificity)
            if pa or pb:
                best = ta if (pa and (not pb or int(ta["overlap"]) >= int(tb["overlap"]))) else tb
                # An overlap of any n genes is not evidence of one population. Activation genes,
                # cell-cycle genes and lineage-shared genes are held in common by populations
                # that are not the same cell type, and on 7 lung datasets that merged HLCA
                # haematopoietic stem cells into megakaryocytes on 10 shared genes while PRSS57
                # (r = 0.82), SPINK2, AVP, CRHBP and CYTL1 were absent from the target, and
                # merged MAIT cells and conventional CD4 T cells into regulatory T cells on
                # pan-T costimulatory genes.
                # Redundancy therefore also requires the candidate's OWN identity markers, its
                # highest-correlation genes, to be present in the target state. Measured on those
                # datasets, every correct call carried 8 to 10 of the candidate's top 10 markers
                # in the target and every incorrect call carried 1 to 4.
                ident = list(q[:int(nomination_identity_top)])
                db = db_a if (best is ta) else db_b
                shared_ident = len(set(ident) & set(db.get(str(best["state"]), [])))
                if shared_ident < int(nomination_identity_min):
                    provisional.append(c)
                    _ev(c, "retained", best, db, shared_ident,
                        cr=_centroid_r(f"{run.name}|{c}",
                                       target_nodes.get(str(best["state"]), [])),
                        failed="identity",
                        reason=f"shares {shared_ident} of its top "
                               f"{int(nomination_identity_top)} markers, below "
                               f"{int(nomination_identity_min)}")
                    _log(f"{run.name}: {c} overlaps {best['state']} by "
                         f"{int(best['overlap'])} genes but shares only {shared_ident} of its "
                         f"top {int(nomination_identity_top)} markers, below "
                         f"{int(nomination_identity_min)}; kept as a candidate")
                    continue
                # Second redundancy test, on the head of each marker list only. Enrichment
                # over 60 query genes can be carried by mid-ranked genes a lineage shares. Two
                # clusters that describe one population also agree on their strongest genes.
                # The test therefore takes the top `redundancy_top_n` markers of the candidate
                # and of the target and requires `redundancy_min_fraction` of that many genes in
                # common. A candidate below the requirement is kept and goes to the survival
                # gate rather than being excluded.
                head_ov = np.nan
                if int(redundancy_top_n) > 0:
                    n_top = int(redundancy_top_n)
                    need = int(np.ceil(float(redundancy_min_fraction) * n_top))
                    head_q = list(q[:n_top])
                    head_t = list(db.get(str(best["state"]), []))[:n_top]
                    head_ov = len(set(head_q) & set(head_t))
                    if head_ov < need:
                        provisional.append(c)
                        _ev(c, "retained", best, db, shared_ident, head_ov,
                            cr=_centroid_r(f"{run.name}|{c}",
                                           target_nodes.get(str(best["state"]), [])),
                            failed="head_overlap",
                            reason=f"shares {head_ov} of the top {n_top} markers of each, "
                                   f"below {need}")
                        _log(f"{run.name}: {c} enriched in {best['state']} but shares only "
                             f"{head_ov} of the top {n_top} markers of each "
                             f"({len(head_q)} and {len(head_t)} available), below {need} "
                             f"({100 * float(redundancy_min_fraction):.0f}%); kept as a candidate")
                        continue
                # Centroid condition. Gene overlap alone can be carried by a shared programme,
                # so redundancy also requires the two mean expression profiles to agree over the
                # whole shared gene space.
                cr = _centroid_r(f"{run.name}|{c}", target_nodes.get(str(best["state"]), []))
                if float(redundancy_centroid_r) > 0:
                    if not np.isfinite(cr) or cr < float(redundancy_centroid_r):
                        provisional.append(c)
                        _ev(c, "retained", best, db, shared_ident, head_ov, cr,
                            failed="centroid_r",
                            reason=f"centroid correlation {cr:.3f} with {best['state']}, below "
                                   f"{float(redundancy_centroid_r):.2f}")
                        _log(f"{run.name}: {c} enriched in {best['state']} but its centroid "
                             f"correlates at {cr:.3f}, below "
                             f"{float(redundancy_centroid_r):.2f}; kept as a candidate")
                        continue
                excluded_nodes.append(f"{run.name}|{c}")
                _ev(c, "excluded", best, db, shared_ident, head_ov, cr, failed="",
                    reason=("redundant with an excluded cluster"
                            if str(best["state"]).startswith("X|")
                            else "redundant with a retained state"))
                enrich_hit[c] = str(best["state"])
                withdrawn[c] = (f"markers specifically enriched in reference state "
                                f"{best['state']} ({int(best['overlap'])} of "
                                f"{int(best['required_overlap'])} required genes, FDR "
                                f"{float(best['fdr']):.2e}; {shared_ident} of its top "
                                f"{int(nomination_identity_top)} markers shared"
                                + (f"; {head_ov} of the top {int(redundancy_top_n)} markers "
                                   f"of each shared" if int(redundancy_top_n) > 0 else "")
                                + ")")
            else:
                provisional.append(c)
                near = ta if ta is not None else tb
                if near is None and ea.shape[0]:
                    near = ea.iloc[0]
                _ev(c, "retained", near, db_a, np.nan, np.nan,
                    cr=(_centroid_r(f"{run.name}|{c}",
                                    target_nodes.get(str(near["state"]), []))
                        if near is not None else np.nan),
                    failed="enrichment",
                    reason="no comparison target reached the required overlap at the set FDR")
        _log(f"{run.name}: gate 1, {len(run.clusters) - len(provisional)} of "
             f"{len(run.clusters)} clusters excluded as redundant; "
             f"{len(provisional)} proceed to the survival test")

        # Gate 2, survival, evaluated ONE CANDIDATE AT A TIME.
        #
        # Testing all candidates together cannot attribute damage. Admitting twenty candidates
        # spreads a reference state's markers across many of them, and the rule then withdraws
        # whichever took the most, which was measured at one gene: the Migratory DC clusters of
        # Adams and Natri-2024, holding 19 and 26 unique markers of their own, were withdrawn for
        # supposedly erasing state S30, a fibroblast state (CTHRC1, FAP, LUM, MXRA5) sharing only
        # CCL19 with them. Adding one candidate at a time isolates its actual effect.
        #
        # Reference cells are subsampled for this test because only marker structure matters.
        from altanalyze3.components.udon.markerFinder import marker_finder_wrapper

        rng = np.random.default_rng(0)
        ref_by_state: Dict[str, List[int]] = {}
        for i in keep:
            if not labels_all[i].startswith("CAND|"):
                ref_by_state.setdefault(labels_all[i], []).append(i)
        ref_rows: List[int] = []
        for st_label, idxs in ref_by_state.items():
            take = idxs if len(idxs) <= survival_ref_cells else \
                list(rng.choice(idxs, size=survival_ref_cells, replace=False))
            ref_rows.extend(int(x) for x in take)
        ref_rows.sort()
        _log(f"{run.name}: gate 2 reference panel, {len(ref_by_state)} states, "
             f"{len(ref_rows)} cells at up to {survival_ref_cells} per state")

        df_ref = pd.DataFrame(pool[ref_rows].toarray().astype(np.float32),
                              index=[barcodes[i] for i in ref_rows], columns=genes)
        g_ref = pd.DataFrame({"cluster": [labels_all[i] for i in ref_rows]}, index=df_ref.index)
        _, m0, _ = marker_finder_wrapper(input_df=df_ref, groups=g_ref, top_n=marker_top_n,
                                         rho_threshold=survival_rho,
                                         marker_finder_rho=survival_rho,
                                         min_markers_per_cluster=survival_min_markers)
        baseline = m0.groupby("top_cluster").size() if m0.shape[0] else pd.Series(dtype=int)
        n_below = sum(1 for st in states
                      if int(baseline.get(f"S{st.sid}", 0)) < int(survival_min_markers))
        _log(f"{run.name}: baseline, {len(states) - n_below} of {len(states)} reference states "
             f"hold >= {survival_min_markers} unique markers before any candidate is added")

        cells_by_cand: Dict[str, List[int]] = {}
        for i in keep:
            if labels_all[i].startswith("CAND|"):
                cells_by_cand.setdefault(labels_all[i].split("|", 1)[1], []).append(i)

        measured: Dict[str, int] = {}
        admitted_list: List[str] = []
        # --force-admit overrides gate 1, which judges REDUNDANCY. It must not override gate 2,
        # which protects existing states from having their markers erased, and which also
        # requires a candidate to hold its own markers. Gate 2 iterates `provisional`, so a
        # gate-1-rejected cluster would never reach it. Add the forced clusters here.
        gate2_candidates = list(provisional)
        forced_here = [c for c in run.clusters
                       if f"{run.name}|{c}" in force_admit and c not in set(provisional)]
        if forced_here:
            gate2_candidates.extend(forced_here)
            _log(f"{run.name}: {len(forced_here)} forced cluster(s) skip gate 1 and are sent "
                 f"through gate 2: {forced_here[:6]}")
        for c in gate2_candidates:
            rows_c = cells_by_cand.get(c, [])
            if len(rows_c) < 5:
                withdrawn[c] = "too few cells to test"
                continue
            idx = ref_rows + rows_c
            df_c = pd.DataFrame(pool[idx].toarray().astype(np.float32),
                                index=[barcodes[i] for i in idx], columns=genes)
            g_c = pd.DataFrame({"cluster": [labels_all[i] for i in idx]}, index=df_c.index)
            _, mt, _ = marker_finder_wrapper(input_df=df_c, groups=g_c, top_n=marker_top_n,
                                             rho_threshold=survival_rho,
                                             marker_finder_rho=survival_rho,
                                             min_markers_per_cluster=survival_min_markers)
            cnt = mt.groupby("top_cluster").size() if mt.shape[0] else pd.Series(dtype=int)
            own = int(cnt.get(f"CAND|{c}", 0))
            measured[c] = own
            if own < int(survival_min_markers):
                withdrawn[c] = (f"holds {own} unique markers against the reference, "
                                f"{survival_min_markers} required")
                continue
            # The damage test protects an existing state from losing its markers to a new one.
            # Applied without regard to the relative evidence it inverts: a state holding a
            # single weak marker vetoes a candidate holding many strong ones, because marker
            # assignment is winner-take-all and recomputed over all states, so ANY admission can
            # take that last gene. On 7 lung datasets this excluded PedDev Schwann cells, which
            # own 60 markers up to r = 0.94 (PLP1, S100B, MPZ, SOX10) and share no gene with the
            # state they were said to damage, TGEN MyoFB, whose best marker reaches r = 0.44 and
            # which ended the run with no unique markers of its own.
            # A state may therefore veto a candidate only when the state is at least as well
            # supported as the candidate. Better-supported evidence wins; the weaker state is
            # left for the survival gate, which is the step that exists to remove it.
            hurt = [st for st in states
                    if int(baseline.get(f"S{st.sid}", 0)) >= int(survival_min_markers)
                    and int(cnt.get(f"S{st.sid}", 0)) < int(damage_floor)
                    and (legacy_damage_test
                         or int(baseline.get(f"S{st.sid}", 0)) >= own)]
            if not legacy_damage_test:
                waived = [st for st in states
                          if int(baseline.get(f"S{st.sid}", 0)) >= int(survival_min_markers)
                          and int(cnt.get(f"S{st.sid}", 0)) < int(damage_floor)
                          and int(baseline.get(f"S{st.sid}", 0)) < own]
                for st in waived:
                    _log(f"{run.name}: {c} holds {own} unique markers and would reduce "
                         f"S{st.sid} from {int(baseline.get(f'S{st.sid}', 0))} to "
                         f"{int(cnt.get(f'S{st.sid}', 0))}; the better supported cluster is "
                         f"admitted and S{st.sid} is left to the survival gate")
            if hurt:
                st = hurt[0]
                withdrawn[c] = (f"admission alone drops reference state S{st.sid} from "
                                f"{int(baseline.get(f'S{st.sid}', 0))} to "
                                f"{int(cnt.get(f'S{st.sid}', 0))} unique markers, and S{st.sid} "
                                f"is at least as well supported as this cluster "
                                f"({int(baseline.get(f'S{st.sid}', 0))} vs {own})")
                continue
            admitted_list.append(c)
        provisional = admitted_list
        _log(f"{run.name}: gate 2, {len(admitted_list)} of {len(cells_by_cand)} candidates hold "
             f"their own markers without erasing a reference state")

        excluded, admitted = 0, 0
        for c in run.clusters:
            n_unique = int(measured.get(c, 0))
            sub = cand_tbl.loc[(cand_tbl["top_cluster"] == c)
                               & (cand_tbl["marker"].isin(gene_set))]
            query = list(sub.sort_values("pearson_r", ascending=False)
                         .head(int(nomination_query_top))["marker"])
            e_a = enrichment_against_states(query, db_a, universe,
                                            overlap_fraction=nomination_overlap_fraction,
                                            min_overlap_floor=nomination_min_overlap)
            e_b = enrichment_against_states(query, db_b, universe,
                                            overlap_fraction=nomination_overlap_fraction,
                                            min_overlap_floor=nomination_min_overlap)
            best_a = e_a.iloc[0] if e_a.shape[0] else None
            best_b = e_b.iloc[0] if e_b.shape[0] else None
            row = {"step": run.name, "cluster": c, "unique_markers": n_unique,
                   "required_unique_markers": int(survival_min_markers),
                   "rho": float(survival_rho),
                   "withdrawn_reason": withdrawn.get(c, ""),
                   "nearest_state_r025": str(best_a["state"]) if best_a is not None else "",
                   "overlap_r025": int(best_a["overlap"]) if best_a is not None else 0,
                   "fdr_r025": float(best_a["fdr"]) if best_a is not None else 1.0,
                   "nearest_state_r030": str(best_b["state"]) if best_b is not None else "",
                   "overlap_r030": int(best_b["overlap"]) if best_b is not None else 0,
                   "fdr_r030": float(best_b["fdr"]) if best_b is not None else 1.0}
            node = f"{run.name}|{c}"
            if node in suppress_cluster:
                row["action"] = "suppressed_replaced_by_finer_dataset"
                row["state"] = ""
                audit_rows.append(row)
                evidence_rows.append({
                    "order": len(evidence_rows) + 1, "dataset": run.name, "cluster": c,
                    "source": node, "decision": "excluded", "state": "", "matched_target": "",
                    "target_kind": "", "target_name": "", "overlap": np.nan,
                    "required_overlap": np.nan, "fdr": np.nan, "identity_shared": np.nan,
                    "identity_required": np.nan, "head_overlap": np.nan,
                    "head_required": np.nan, "centroid_r": np.nan,
                    "centroid_r_required": np.nan, "failed_test": "qualification",
                    "reason": "removed before the comparison; not used as a redundancy target"})
                excluded += 1
                continue
            row["forced"] = bool(node in force_admit)
            if c in set(provisional):
                st = new_state()
                st.absorb(f"{run.name}|{c}", run.centroids[c], run.markers.get(c, []))
                states.append(st)
                row["action"] = ("added_new_state_forced" if node in force_admit
                                 else "added_new_state")
                row["state"] = f"S{st.sid}"
                for ev in evidence_rows:
                    if ev["source"] == node and ev["decision"] == "retained":
                        ev["state"] = f"S{st.sid}"
                admitted += 1
            else:
                row["action"] = ("excluded_redundant" if c in enrich_hit
                                 else "excluded_not_admitted")
                row["state"] = ""
                if c not in enrich_hit:
                    excluded_nodes.append(node)
                    for ev in evidence_rows:
                        if ev["source"] == node and ev["decision"] == "retained":
                            ev["decision"] = "excluded"
                            ev["failed_test"] = "survival"
                            ev["reason"] = ("passed the redundancy test but held too few "
                                            "unique markers to be admitted")
                # Only a cluster that PASSED gate 1 evidences that another dataset holds the
                # same population. A cluster withdrawn for any other reason, "too few cells to
                # test" among them, was never compared, and `best_a` is an unfiltered nearest
                # hit that ignores both the FDR threshold and the overlap floor. Counting it
                # would put fabricated replication into `supporting_datasets`, which is the
                # field a reader cites as cross-dataset validation.
                target = enrich_hit.get(c, "")
                if target.startswith("S"):
                    support.setdefault(int(target[1:]), []).append(run.name)
                excluded += 1
            audit_rows.append(row)
        n_gate1 = sum(1 for c in run.clusters
                      if c in enrich_hit and f"{run.name}|{c}" not in force_admit)
        _log(f"+ {run.name}: {len(run.clusters)} clusters tested; {n_gate1} excluded as "
             f"redundant, {excluded - n_gate1} excluded by the survival test, "
             f"{admitted} added as new states; reference now {len(states)}")

    for st in states:
        st.support_datasets = sorted(set(support.get(st.sid, [])))

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(os.path.join(outdir, "hierarchical_audit.tsv"), sep="\t", index=False)
    ev = pd.DataFrame(evidence_rows)
    if ev.shape[0]:
        ev.to_csv(os.path.join(outdir, "state_retention_evidence.tsv"), sep="\t", index=False)
        n_ret = int((ev["decision"] == "retained").sum())
        _log(f"wrote {os.path.join(outdir, 'state_retention_evidence.tsv')}: {ev.shape[0]} "
             f"input clusters, {n_ret} retained, {ev.shape[0] - n_ret} excluded, each with the "
             f"comparison target and every statistic that decided it")
    nomination = audit.loc[audit["action"] != "seed_state"].copy()
    nomination.to_csv(os.path.join(outdir, "nomination_decisions.tsv"), sep="\t", index=False)
    return states, audit, nomination


# ----------------------------------------------------------------------- cell state annotation
def reference_cell_type_enrichment(cell_state: Sequence[str], cell_dataset: Sequence[str],
                                   cell_barcode: Sequence[str], runs: List[RunInputs],
                                   *, columns: Optional[Sequence[str]] = None
                                   ) -> pd.DataFrame:
    """Source 1. Enrich each state for reference cell types carried by its own cells.

    The universe is every final cell that carries the reference column, so the denominator is the
    cell types actually present in this dataset rather than the whole reference vocabulary. For
    state S and reference type T the module computes

        k = cells of S labelled T          n = cells of S carrying the column
        K = all final cells labelled T     N = all final cells carrying the column
        p = hypergeom.sf(k - 1, N, n, K)

    and adjusts p across every state-by-type pair within one column, by Benjamini-Hochberg.

    `accuracy` is k / n, the fraction of the state's own annotated cells carrying T. Precision
    rather than recall: a state is well named when its cells agree, not when it captures every
    cell of that type, because one type may legitimately split across several states.

    Each dataset is scored only on the columns it carries. Two studies annotated with different
    vocabularies never enter the same test, because their label strings are not comparable.
    """
    lut: Dict[str, Dict[str, Dict[str, str]]] = {}
    for r in runs:
        if r.annotations is None:
            continue
        for col in r.annotations.columns:
            if col == "barcode":
                continue
            if columns and col not in set(columns):
                continue
            lut.setdefault(col, {})[r.name] = dict(
                zip(r.annotations["barcode"].astype(str),
                    r.annotations[col].astype(str)))
    if not lut:
        return pd.DataFrame()

    rows = []
    for col, per_run in lut.items():
        labels = []
        for st, ds, bc in zip(cell_state, cell_dataset, cell_barcode):
            raw = bc.split("|", 1)[1] if "|" in bc else bc
            v = per_run.get(ds, {}).get(raw)
            labels.append(v if v not in (None, "", "nan", "NA") else None)
        frame = pd.DataFrame({"state": list(cell_state), "label": labels}).dropna()
        if frame.empty:
            continue
        N = int(frame.shape[0])
        K_of = frame["label"].value_counts().to_dict()
        col_rows = []
        for state, grp in frame.groupby("state"):
            n = int(grp.shape[0])
            for label, k in grp["label"].value_counts().items():
                k = int(k); K = int(K_of[label])
                p = float(hypergeom.sf(k - 1, N, n, K))
                col_rows.append({"state": str(state), "source": "reference_cells",
                                 "reference": col, "annotation": str(label),
                                 "evidence": k, "state_size": n, "reference_size": K,
                                 "universe": N, "accuracy": float(k) / float(n),
                                 "p_value": p})
        if not col_rows:
            continue
        cf = pd.DataFrame(col_rows)
        cf["fdr"] = _bh_adjust(cf["p_value"].to_numpy())
        rows.append(cf)
        _log(f"annotation source 1: reference '{col}' scored {cf['state'].nunique()} states "
             f"over {N} annotated cells and {len(K_of)} cell types present")
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def marker_gene_set_enrichment(markers: pd.DataFrame, background: Sequence[str], *,
                               species: str = "Hs",
                               biomarker_file: Optional[str] = None,
                               outdir: str) -> pd.DataFrame:
    """Source 2. Enrich each state's marker genes against the GO-Elite BioMarkers reference.

    The module calls ICGS3's own `biomarker_enrichment`, so the integration and the source runs
    share one implementation and one reference file. The `background` argument is the gene
    denominator: `biomarker_enrichment` intersects every BioMarkers term with it before testing,
    so the universe is the genes this dataset measured, not the whole genome. Passing the union
    feature space therefore does what GO-Elite calls a restricted background.

    Term labels pass through ICGS3's `clean_biomarker_prediction_labels`, the legacy RNASeq.py
    rules: drop the parenthetical citation, strip Adult/Fetal/Embryonic/Term, infer the dominant
    tissue across clusters and remove it when it is not the whole term.

    `accuracy` is overlap / query size, the fraction of the state's markers inside the term.
    """
    from altanalyze3.components.clustering.ICGS import (
        biomarker_enrichment, clean_biomarker_prediction_labels, _default_biomarker_file)

    sp = {"human": "Hs", "mouse": "Mm", "Hs": "Hs", "Mm": "Mm"}.get(str(species), str(species))
    path = biomarker_file or _default_biomarker_file(sp)
    if not path or not os.path.exists(path):
        _log(f"annotation source 2 skipped: no BioMarkers file for species '{sp}'"
             f"{f' at {biomarker_file}' if biomarker_file else ''}")
        return pd.DataFrame()
    if markers.empty:
        _log("annotation source 2 skipped: no markers")
        return pd.DataFrame()

    class _Cfg:
        pass
    cfg = _Cfg()
    cfg.biomarker_file = path
    cfg.species = sp
    _log(f"annotation source 2: BioMarkers {path}, gene denominator restricted to "
         f"{len(set(background))} genes measured in this integration")
    biomarker_enrichment(markers, background, cfg, outdir)
    enr_path = os.path.join(outdir, "GO-Elite", "icgs3_biomarker_enrichment.tsv")
    if not os.path.exists(enr_path):
        _log("annotation source 2: enrichment produced no result")
        return pd.DataFrame()
    enr = pd.read_csv(enr_path, sep="\t")
    if enr.empty:
        return pd.DataFrame()
    # ICGS3 cleans one row per cluster, its top term, and infers the dominant tissue from that
    # row set. `clean_biomarker_prediction_labels` records the tissue in a dict keyed by cluster,
    # so the LAST qualifying row per cluster casts that cluster's vote. Every enriched term needs
    # a cleaned name here, because the accuracy rule can adopt a term that is not the top hit.
    # Ordering each cluster's top term last therefore reproduces ICGS3's vote exactly while still
    # cleaning every row through the same rules. Section 20.3 of ICGS_INTEGRATE_METHODS.md records
    # the check that the top-term names match ICGS3's own output.
    enr = enr.sort_values(["cluster", "fdr", "p_value", "term_name"],
                          ascending=[True, True, True, True])
    rank = enr.groupby("cluster").cumcount()
    order = pd.Series(np.where(rank == 0, 1, 0), index=enr.index)   # top term sorts last
    enr = enr.assign(_top_last=order).sort_values(["cluster", "_top_last"],
                                                  kind="mergesort").drop(columns="_top_last")
    named = clean_biomarker_prediction_labels(
        enr[["cluster", "term_name", "fdr", "overlap"]].copy())
    # clean_biomarker_prediction_labels appends _c<cluster>; the state id is already the cluster
    clean = {}
    for _, r in named.iterrows():
        lab = str(r["cell_type_prediction"])
        suffix = f"_c{r['cluster']}"
        clean[(str(r["cluster"]), str(r["term_name"]))] = (
            lab[: -len(suffix)] if lab.endswith(suffix) else lab)
    out = pd.DataFrame({
        "state": enr["cluster"].astype(str),
        "source": "marker_gene_sets",
        "reference": f"GO-Elite BioMarkers ({sp})",
        "annotation": [clean.get((str(c), str(tn)), str(tn))
                       for c, tn in zip(enr["cluster"], enr["term_name"])],
        "evidence": enr["overlap"].astype(int),
        "state_size": enr["query_size"].astype(int),
        "reference_size": enr["term_size"].astype(int),
        "universe": len(set(background)),
        "accuracy": enr["overlap"].astype(float) / enr["query_size"].astype(float),
        "p_value": enr["p_value"].astype(float),
        "fdr": enr["fdr"].astype(float)})
    _log(f"annotation source 2: {out.shape[0]} enriched terms across "
         f"{out['state'].nunique()} states")
    return out


def assign_state_annotations(candidates: pd.DataFrame, states: Sequence[str], *,
                             min_evidence: int = 4, fdr: float = 0.05) -> pd.DataFrame:
    """Pick one annotation per state: highest accuracy among candidates that clear both bars.

    A candidate must carry at least `min_evidence` matching cells or genes and reach FDR <= `fdr`.
    Ties break on the smaller FDR, then the larger evidence, then the annotation name, so the
    result never depends on row order.

    The two sources are not on one scale. Cell-label accuracy is a purity, which routinely reaches
    0.7 to 1.0, while marker-overlap accuracy is a fraction of 60 markers inside one term, which
    rarely passes 0.3. Source 1 therefore wins whenever it clears both bars. The behaviour is
    intended: a reference label carried by the state's own cells is direct evidence, and a marker
    set overlap is inferred evidence. The per-source columns record what source 2 would have said,
    so a disagreement stays visible.
    """
    cols = ["state", "annotation", "source", "reference", "evidence", "state_size",
            "reference_size", "universe", "accuracy", "p_value", "fdr"]
    if candidates.empty:
        return pd.DataFrame({"state": list(states)}).assign(
            annotation="", annotation_source="", annotation_reference="",
            annotation_evidence=0, annotation_accuracy=np.nan, annotation_fdr=np.nan)
    ok = candidates.loc[(candidates["evidence"] >= int(min_evidence))
                        & (candidates["fdr"] <= float(fdr))].copy()
    _log(f"annotation: {ok.shape[0]} of {candidates.shape[0]} candidates clear "
         f">= {min_evidence} matching cells or genes at FDR <= {fdr}")
    rows = []
    for st in states:
        sub = ok.loc[ok["state"].astype(str) == str(st)]
        if sub.empty:
            rows.append({"state": str(st), "annotation": "", "annotation_source": "",
                         "annotation_reference": "", "annotation_evidence": 0,
                         "annotation_accuracy": np.nan, "annotation_fdr": np.nan})
            continue
        sub = sub.sort_values(["accuracy", "fdr", "evidence", "annotation"],
                              ascending=[False, True, False, True])
        best = sub.iloc[0]
        rows.append({"state": str(st), "annotation": str(best["annotation"]),
                     "annotation_source": str(best["source"]),
                     "annotation_reference": str(best["reference"]),
                     "annotation_evidence": int(best["evidence"]),
                     "annotation_accuracy": round(float(best["accuracy"]), 4),
                     "annotation_fdr": float(best["fdr"])})
    return pd.DataFrame(rows)


def resolve_dataset_order(names: Sequence[str], cluster_counts: Dict[str, int], *,
                          seed_dataset: Optional[str] = None,
                          dataset_order: Optional[Sequence[str]] = None) -> List[str]:
    """Return the entry order, after checking that every name the user gave exists.

    Called once before any file is read, so a typo fails in a second rather than after the
    pooled matrix is built. `integrate_hierarchically` repeats the checks, because it is also
    callable directly.
    """
    known = set(names)
    if dataset_order:
        wanted = [str(n).strip() for n in dataset_order if str(n).strip()]
        missing = [n for n in wanted if n not in known]
        if missing:
            raise ValueError(f"--dataset-order names datasets that were not given with --run: "
                             f"{missing}. Available: {sorted(known)}")
        if len(wanted) != len(set(wanted)):
            dup = sorted({n for n in wanted if wanted.count(n) > 1})
            raise ValueError(f"--dataset-order repeats a dataset: {dup}")
        omitted = sorted(known - set(wanted))
        if omitted:
            raise ValueError(f"--dataset-order must list every dataset given with --run. "
                             f"Missing: {omitted}")
        return wanted
    order = sorted(known, key=lambda n: cluster_counts.get(n, 0), reverse=True)
    if seed_dataset:
        if seed_dataset not in known:
            raise ValueError(f"--seed-dataset '{seed_dataset}' was not given with --run. "
                             f"Available: {sorted(known)}")
        order = [seed_dataset] + [n for n in order if n != seed_dataset]
    return order


def states_to_frame(states: List[HState], runs: List[RunInputs]) -> pd.DataFrame:
    name_of, fdr_of, size_of = {}, {}, {}
    for r in runs:
        for c in r.clusters:
            key = f"{r.name}|{c}"
            name_of[key] = r.names.get(c, "")
            fdr_of[key] = r.name_fdr.get(c, np.nan)
            size_of[key] = r.cluster_sizes.get(c, 0)
    rows = []
    for st in states:
        labels = [name_of.get(m, "") for m in st.members if name_of.get(m, "")]
        vote = pd.Series(labels).value_counts() if labels else pd.Series(dtype=int)
        rows.append({
            "state": f"S{st.sid}",
            "n_members": len(st.members),
            "source_dataset": ",".join(sorted(st.datasets)),
            "n_supporting_datasets": 1 + len(getattr(st, "support_datasets", [])),
            "supporting_datasets": ",".join(getattr(st, "support_datasets", [])),
            "members": ",".join(st.members),
            "n_cells": int(sum(size_of.get(m, 0) for m in st.members)),
            "name": vote.index[0] if vote.size else "",
            "name_agreement": round(float(vote.iloc[0] / vote.sum()), 3) if vote.size else np.nan,
            "competing_names": ";".join(f"{k}({v})" for k, v in vote.items()) if vote.size else "",
            "best_name_fdr": _safe_min_fdr([fdr_of.get(m, np.nan) for m in st.members]),
            "top_markers": ",".join(st.ranked_markers(15)),
            "kind": ("supported_by_multiple_datasets"
                     if getattr(st, "support_datasets", []) else "single_dataset_only"),
        })
    return pd.DataFrame(rows).sort_values(["n_supporting_datasets", "n_cells"],
                                          ascending=[False, False])



# ------------------------------------------------- union feature space, redundancy, HOPACH
def union_feature_space(runs: List[RunInputs]) -> pd.Index:
    """Every gene any dataset selected as informative, not the intersection.

    Each ICGS3 run reports only the genes ITS MarkerFinder chose, so a gene informative in one
    dataset is simply absent from another dataset's centroid file. Comparing on the intersection
    would discard exactly the genes that separate datasets. Scoring uses the shared genes of each
    pair, but the FINAL matrix is built on the union and filled from the h5ad, so no state is
    described by a gene list truncated to another dataset's choices.
    """
    genes = pd.Index([])
    for r in runs:
        genes = genes.union(r.centroids.index)
    return genes


# A dataset read over less than this share of the shared gene space enters the comparison as
# mostly zeros, so its centroids are not comparable with the rest.
MIN_UNION_GENE_COVERAGE = 0.90


def _union_gene_share(h5_path: str, genes: pd.Index) -> float:
    """Fraction of `genes` an h5ad carries. Reads var only, so it costs almost nothing."""
    if not h5_path or not os.path.exists(h5_path):
        return 0.0
    import h5py

    try:
        with h5py.File(h5_path, "r") as handle:
            var = handle["var"]
            key = var.attrs.get("_index", "_index")
            key = key.decode() if isinstance(key, bytes) else str(key)
            names = [x.decode() if isinstance(x, bytes) else x for x in var[key][:]]
    except (OSError, KeyError):
        return 0.0
    if len(genes) == 0:
        return 0.0
    return len(pd.Index(names).intersection(genes)) / float(len(genes))


def _existing_path(path: str) -> Optional[str]:
    """Return `path`, or the local mount that holds it, or None.

    A run submitted on the cluster records a /data path. The same share mounts on a Mac under
    /Volumes, and the user's home under /data/saljh8 is /Users/saljh8 locally.
    """
    if not path:
        return None
    if os.path.exists(path):
        return path
    for old, new in (("/data/saljh8", "/Users/saljh8"),
                     ("/data/", "/Volumes/")):
        if path.startswith(old):
            stem = path[len(old):]
            for root in (new, new.rstrip("/") + "-1/"):
                candidate = os.path.join(root, stem) if root.endswith("/") else root + stem
                candidate = candidate.replace("//", "/")
                if os.path.exists(candidate):
                    return candidate
            # /Volumes mounts a second copy as <share>-1 when the first is still attached.
            head = stem.split("/", 1)
            if len(head) == 2:
                for suffix in ("", "-1", "-2"):
                    candidate = os.path.join(new, head[0] + suffix, head[1])
                    if os.path.exists(candidate):
                        return candidate
    return None


def source_h5ad_for_run(run_path: str) -> Optional[str]:
    """The h5ad an ICGS3 run started from, read off its own icgs3_config.json.

    The run's icgs3_result.h5ad is not always a full matrix: a streamed load writes only the
    selected variable features and drops layers['counts']. The source file always carries every
    gene and the counts layer, so the pooled comparison reads it instead.
    """
    config = os.path.join(run_path, "icgs3_config.json")
    if not os.path.exists(config):
        return None
    try:
        with open(config) as handle:
            paths = json.load(handle).get("input_paths") or []
    except (OSError, ValueError):
        return None
    for path in paths:
        resolved = _existing_path(str(path))
        if resolved:
            return resolved
    return None


def _rows_from_h5ad(path: str, want: Sequence[str], genes: pd.Index, *, with_counts: bool
                    ) -> Tuple[sp.csr_matrix, Optional[sp.csr_matrix], List[str], int]:
    """Read the `want` cells of one h5ad into the shared `genes` space, with h5py.

    Reads consecutive wanted rows as one slice. One h5py read per cell would decompress the
    containing chunk once per cell, which costs minutes on a large gzip-chunked file.
    """
    import h5py

    with h5py.File(path, "r", rdcc_nbytes=256 * 1024 * 1024, rdcc_nslots=100003) as handle:
        names = np.array([x.decode() if isinstance(x, bytes) else x
                          for x in handle["obs"]["_index"][:]])
        var = np.array([x.decode() if isinstance(x, bytes) else x
                        for x in handle["var"]["_index"][:]])
        position = {name: i for i, name in enumerate(names)}
        kept = [b for b in want if b in position]
        rows = np.array([position[b] for b in kept], dtype=np.int64)
        if rows.size == 0:
            empty = sp.csr_matrix((0, len(genes)), dtype=np.float32)
            return empty, (empty if with_counts else None), [], 0
        # Read in file order, then restore the order of `kept` with `inverse`.
        order = np.argsort(rows)
        sorted_rows = rows[order]
        inverse = np.empty(len(order), dtype=np.int64)
        inverse[order] = np.arange(len(order))

        cols = pd.Index(var)
        take = cols.intersection(genes)
        pos = genes.get_indexer(take)
        col_map = np.full(len(var), -1, dtype=np.int64)
        col_map[cols.get_indexer(take)] = pos

        def read_block(group) -> sp.csr_matrix:
            indptr = group["indptr"][:].astype(np.int64)
            data_src, idx_src = group["data"], group["indices"]
            starts, ends = indptr[sorted_rows], indptr[sorted_rows + 1]
            lengths = ends - starts
            total = int(lengths.sum())
            data = np.empty(total, dtype=np.float32)
            indices = np.empty(total, dtype=np.int64)
            offset = j = 0
            while j < len(sorted_rows):
                k = j
                while k + 1 < len(sorted_rows) and sorted_rows[k + 1] == sorted_rows[k] + 1:
                    k += 1
                span_start, span_end = int(starts[j]), int(ends[k])
                span = span_end - span_start
                data[offset:offset + span] = data_src[span_start:span_end]
                indices[offset:offset + span] = idx_src[span_start:span_end]
                offset += span
                j = k + 1
            out_indptr = np.zeros(len(sorted_rows) + 1, dtype=np.int64)
            out_indptr[1:] = np.cumsum(lengths)
            mapped = col_map[indices]
            keep_mask = mapped >= 0
            if not keep_mask.all():
                per_row = np.add.reduceat(keep_mask.astype(np.int64), out_indptr[:-1]) \
                    if total else np.zeros(len(sorted_rows), dtype=np.int64)
                per_row = np.where(np.diff(out_indptr) > 0, per_row, 0)
                out_indptr = np.zeros(len(sorted_rows) + 1, dtype=np.int64)
                out_indptr[1:] = np.cumsum(per_row)
                data, mapped = data[keep_mask], mapped[keep_mask]
            matrix = sp.csr_matrix((data, mapped, out_indptr),
                                   shape=(len(sorted_rows), len(genes)))
            return matrix[inverse]

        block = read_block(handle["X"])
        counts = None
        if with_counts and "layers" in handle and "counts" in handle["layers"]:
            counts = read_block(handle["layers/counts"])
    return block, counts, kept, len(take)


def expression_for_representative_cells(runs: List[RunInputs], reps: pd.DataFrame,
                                        genes: pd.Index, *, with_counts: bool = True,
                                        cells_from: str = "auto"
                                        ) -> Tuple[sp.csr_matrix, List[str], List[str], Optional[sp.csr_matrix]]:
    """Pull the representative cells of every dataset into one cells x union-genes matrix.

    Every dataset is read from the SAME kind of file, because a centroid built on one dataset's
    truncated gene list is not comparable with a centroid built on another's. `cells_from`:

      source  the h5ad the run started from. It carries every gene and layers['counts'].
      result  the run's own icgs3_result.h5ad. A streamed load writes only the selected
              variable features there and drops the counts layer, so coverage can collapse.
      auto    source for every run when every run resolves one, otherwise result.

    Reads each h5ad once and keeps only the selected rows. Memory stays at
    (cells kept) x (union genes).
    """
    chosen = {}
    for run in runs:
        result = os.path.join(run.path, "icgs3_result.h5ad")
        source = source_h5ad_for_run(run.path)
        if cells_from == "result":
            chosen[run.name] = (result, "result")
            continue
        if cells_from == "source":
            if not source:
                raise FileNotFoundError(
                    f"--cells-from source: {run.name} resolves no input h5ad from its "
                    f"icgs3_config.json")
            chosen[run.name] = (source, "source")
            continue
        # auto: keep the run's own result h5ad when it spans the shared gene space, and fall
        # back to the source only for the runs where it does not. A streamed ICGS3 load stores
        # only the selected variable features, which is what makes that check necessary.
        share = _union_gene_share(result, genes)
        if share >= MIN_UNION_GENE_COVERAGE:
            chosen[run.name] = (result, "result")
        elif source:
            _log(f"{run.name}: icgs3_result.h5ad spans {100 * share:.1f}% of the union gene "
                 f"space; reading its cells from the source h5ad instead ({source})")
            chosen[run.name] = (source, "source")
        else:
            raise FileNotFoundError(
                f"{run.name}: icgs3_result.h5ad spans {100 * share:.1f}% of the union gene "
                f"space and no source h5ad resolves from its icgs3_config.json")
    _log("cell source per dataset: " + ", ".join(f"{n}={k}" for n, (_, k) in sorted(chosen.items())))

    blocks, barcodes, datasets, count_blocks = [], [], [], []
    coverage = {}
    for run in runs:
        want = reps.loc[reps["dataset"] == run.name, "barcode"].astype(str).tolist()
        if not want:
            continue
        h5, kind = chosen[run.name]
        if not h5 or not os.path.exists(h5):
            _log(f"WARNING {run.name}: no {kind} h5ad; skipping its cells in the rebuild")
            continue
        block, cblock, keep, n_take = _rows_from_h5ad(h5, want, genes, with_counts=with_counts)
        if not keep:
            _log(f"WARNING {run.name}: none of its {len(want)} representative barcodes are in "
                 f"{h5}; skipping its cells")
            continue
        blocks.append(block)
        if with_counts:
            count_blocks.append(cblock if cblock is not None
                                else sp.csr_matrix((len(keep), len(genes)), dtype=np.float32))
        barcodes.extend([f"{run.name}|{b}" for b in keep])
        datasets.extend([run.name] * len(keep))
        coverage[run.name] = n_take / float(len(genes)) if len(genes) else 0.0
        _log(f"{run.name}: pulled {len(keep)} representative cells x {n_take} of "
             f"{len(genes)} union genes ({100 * coverage[run.name]:.1f}%)"
             f"{'' if cblock is not None else ', NO counts layer'}")
    if not blocks:
        raise RuntimeError("no representative expression could be read")
    # A dataset measured on a fraction of the shared gene space contributes near-empty
    # centroids, which silently distorts every comparison it enters. Fail instead.
    poor = {n: v for n, v in coverage.items() if v < MIN_UNION_GENE_COVERAGE}
    if poor:
        detail = ", ".join(f"{n} {100 * v:.1f}%" for n, v in sorted(poor.items()))
        raise RuntimeError(
            f"gene-space coverage below {100 * MIN_UNION_GENE_COVERAGE:.0f}% of the "
            f"{len(genes)}-gene union for: {detail}. Those datasets would enter the comparison "
            f"as mostly zeros. Re-run them so their h5ad holds the full gene space, or pass "
            f"--cells-from source so the original input h5ad is read.")
    counts = sp.vstack(count_blocks, format="csr") if (with_counts and count_blocks) else None
    return sp.vstack(blocks, format="csr"), barcodes, datasets, counts



def recompute_centroids_from_pool(matrix: sp.csr_matrix, genes: pd.Index,
                                  datasets: List[str], clusters: List[str]) -> Dict[str, pd.DataFrame]:
    """Per-dataset cluster centroids computed from the pooled union matrix.

    Every centroid now spans the SAME union gene space, so a gene missing from one dataset's
    MarkerFinder selection is still measured there. Comparing centroids built from each dataset's
    own truncated gene list would compare different feature spaces and give garbage.
    """
    key = pd.Series([f"{d}|{c}" for d, c in zip(datasets, clusters)])
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for label, idx in key.groupby(key).groups.items():
        rows = np.asarray(idx, dtype=int)
        block = matrix[rows]
        mean = np.asarray(block.mean(axis=0)).ravel().astype(np.float32)
        ds, cl = str(label).split("|", 1)
        out.setdefault(ds, {})[cl] = mean
    frames = {ds: pd.DataFrame(cols, index=genes).astype(np.float32) for ds, cols in out.items()}
    for ds, f in frames.items():
        _log(f"{ds}: recomputed {f.shape[1]} centroids over {f.shape[0]} union genes")
    return frames


def hopach_order_states(states: List["HState"], genes: pd.Index) -> Tuple[List[int], pd.DataFrame]:
    """Cluster and order the harmonized centroids with HOPACH before the final MarkerFinder.

    HOPACH gives the ordering the ICGS and AltAnalyze heatmaps have always used, so the final
    figure groups related states together instead of listing them by the accident of merge order.
    """
    from altanalyze3.components.clustering.hopach import hopach as run_hopach

    frame = pd.DataFrame({f"S{st.sid}": st.centroid for st in states}).reindex(genes).fillna(0.0)
    mat = frame.to_numpy(dtype=np.float64).T          # states x genes
    _log(f"HOPACH on {mat.shape[0]} state centroids x {mat.shape[1]} union genes")
    res = run_hopach(mat, d="cosangle", kmax=9, kmin=2, mincluster=2, random_state=0)
    order = [int(i) for i in np.asarray(res.order).ravel()]
    labels = np.asarray(res.clust.labels).ravel()
    table = pd.DataFrame({"state": [f"S{states[i].sid}" for i in order],
                          "hopach_position": range(1, len(order) + 1),
                          "hopach_cluster": [labels[i] for i in order]})
    _log(f"HOPACH produced {len(set(labels))} top-level groups over {len(order)} states")
    return order, table


def final_markerfinder(matrix: sp.csr_matrix, genes: pd.Index, barcodes: List[str],
                       labels: List[str], outdir: str, top_n: int = 60) -> pd.DataFrame:
    """Run MarkerFinder once on the pooled representative cells and the harmonized labels."""
    from altanalyze3.components.udon.markerFinder import marker_finder_wrapper

    # marker_finder expects observations as the index and features as the columns, so the
    # pooled matrix goes in as cells x genes, not genes x cells.
    df = pd.DataFrame(matrix.toarray().astype(np.float32), index=barcodes, columns=genes)
    groups = pd.DataFrame({"cluster": labels}, index=barcodes)
    _log(f"final MarkerFinder on {df.shape[0]} cells x {df.shape[1]} genes, "
         f"{groups['cluster'].nunique()} harmonized states "
         f"({df.memory_usage(deep=True).sum() / 1e9:.2f} GB dense)")
    markers_all, markers_top, heat = marker_finder_wrapper(
        input_df=df, groups=groups, top_n=top_n,
        rho_threshold=0.3, marker_finder_rho=0.3, min_markers_per_cluster=2)
    markers_all.to_csv(os.path.join(outdir, "harmonized_markers_all.tsv"), sep="\t", index=False)
    markers_top.to_csv(os.path.join(outdir, "harmonized_markers.tsv"), sep="\t", index=False)
    _log(f"final MarkerFinder: {markers_top.shape[0]} markers across "
         f"{markers_top['top_cluster'].nunique()} states")
    heat.to_csv(os.path.join(outdir, "harmonized_heatmap_matrix.tsv"), sep="\t")
    return markers_top



# ------------------------------------------------ survival test, re-MarkerFinder, deliverables
def markerfinder_survival(pool: sp.csr_matrix, genes: pd.Index, barcodes: List[str],
                          labels: List[str], *, rho: float, min_markers: int,
                          top_n: int, max_rounds: int = 10,
                          centroids: Optional[pd.DataFrame] = None,
                          datasets_of_state: Optional[Dict[str, set]] = None,
                          enforce: bool = True):
    """MarkerFinder, then the survival test, then MERGE the states that fail.

    unique_marker_finder gives each gene to exactly one cluster. Running it over 107 states at
    once starves the losers: a state whose genes were all claimed by a near neighbour ends with
    zero markers even though its cells are real. Deleting it would throw away cells and leave
    fewer states than a single input dataset already resolved.

    A state with no unique marker is, by the definition of this test, not distinguishable from
    the state that took its markers. So it is MERGED into that state when the datasets allow,
    and only dropped when no legal partner exists. Cells are never discarded silently.
    """
    from altanalyze3.components.udon.markerFinder import marker_finder_wrapper

    labels = list(labels)
    merged_log, dropped_log = [], []
    markers_all = markers_top = heat = None
    for rnd in range(1, max_rounds + 1):
        df = pd.DataFrame(pool.toarray().astype(np.float32), index=barcodes, columns=genes)
        groups = pd.DataFrame({"cluster": labels}, index=df.index)
        present = sorted(set(labels))
        markers_all, markers_top, heat = marker_finder_wrapper(
            input_df=df, groups=groups, top_n=top_n,
            rho_threshold=rho, marker_finder_rho=rho, min_markers_per_cluster=min_markers)
        counts = markers_top.groupby("top_cluster").size() if markers_top.shape[0] else pd.Series(dtype=int)
        survivors = [c for c in present if int(counts.get(c, 0)) >= int(min_markers)]
        failed = [c for c in present if c not in set(survivors)]
        _log(f"survival round {rnd}: {len(survivors)} of {len(present)} states own >= "
             f"{min_markers} unique markers at rho >= {rho}; {len(failed)} fail")
        if not failed or not survivors:
            break
        if not enforce:
            _log(f"   survival-mode=report: keeping all {len(present)} states; "
                 f"{len(failed)} recorded as low-marker")
            for f in failed:
                dropped_log.append({"round": rnd, "state": f,
                                    "unique_markers": int(counts.get(f, 0)),
                                    "reason": "low unique markers, kept under report mode"})
            break

        # each failing state joins its closest surviving state, if the datasets permit it
        remap: Dict[str, str] = {}
        if centroids is not None:
            cen = centroids
            surv_cols = [c for c in survivors if c in cen.columns]
            M = cen[surv_cols].to_numpy(dtype=np.float64)
            Mn = np.linalg.norm(M, axis=0); Mn[Mn == 0] = 1.0
            for f in failed:
                if f not in cen.columns:
                    continue
                v = cen[f].to_numpy(dtype=np.float64)
                nv = np.linalg.norm(v) or 1.0
                cos = (M / Mn).T @ (v / nv)
                order = np.argsort(-cos)
                for k in order:
                    target = surv_cols[int(k)]
                    if datasets_of_state is not None:
                        if datasets_of_state.get(f, set()) & datasets_of_state.get(target, set()):
                            continue          # would place one dataset twice
                    remap[f] = target
                    merged_log.append({"round": rnd, "state": f, "merged_into": target,
                                       "centroid_cosine": round(float(cos[int(k)]), 4),
                                       "unique_markers": int(counts.get(f, 0))})
                    if datasets_of_state is not None:
                        datasets_of_state[target] = datasets_of_state.get(target, set()) | \
                                                    datasets_of_state.get(f, set())
                    break
        for f in failed:
            if f not in remap:
                dropped_log.append({"round": rnd, "state": f,
                                    "unique_markers": int(counts.get(f, 0)),
                                    "reason": "no legal merge partner"})
        if not remap:
            labels = [l for l in labels]
            keep = [i for i, l in enumerate(labels) if l not in set(failed)]
            pool, barcodes, labels = pool[keep], [barcodes[i] for i in keep], [labels[i] for i in keep]
            continue
        labels = [remap.get(l, l) for l in labels]
        drop = {d["state"] for d in dropped_log}
        if drop:
            keep = [i for i, l in enumerate(labels) if l not in drop]
            pool, barcodes, labels = pool[keep], [barcodes[i] for i in keep], [labels[i] for i in keep]
        _log(f"   merged {len(remap)} starved states into their closest surviving state")
    return (markers_all, markers_top, heat, labels, barcodes, pool,
            pd.DataFrame(merged_log), pd.DataFrame(dropped_log))


def write_integration_deliverables(pool: sp.csr_matrix, genes: pd.Index, barcodes: List[str],
                                   labels: List[str], datasets: List[str], clusters: List[str],
                                   names: Dict[str, str], outdir: str, *, top_n: int,
                                   heatmap_cells: int,
                                   counts: Optional[sp.csr_matrix] = None) -> None:
    """Cell-level annotated MarkerFinder result, centroids and the canonical heatmap PDF."""
    import anndata as ad

    obs = pd.DataFrame({
        "harmonized_state": labels,
        "state_name": [names.get(l, "") for l in labels],
        "source_dataset": datasets,
        "source_cluster": clusters,
    }, index=pd.Index(barcodes, name="cell"))
    obs.to_csv(os.path.join(outdir, "harmonized_cell_annotations.tsv"), sep="\t")

    adata = ad.AnnData(X=pool.tocsr(), obs=obs,
                       var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    if counts is not None:
        adata.layers["counts"] = counts.tocsr()
    adata.write_h5ad(os.path.join(outdir, "harmonized_integrated.h5ad"))

    centroids = {}
    lab = pd.Series(labels)
    for state, idx in lab.groupby(lab).groups.items():
        rows = np.asarray(idx, dtype=int)
        centroids[str(state)] = np.asarray(pool[rows].mean(axis=0)).ravel().astype(np.float32)
    pd.DataFrame(centroids, index=genes).to_csv(
        os.path.join(outdir, "harmonized_final_centroids.tsv"), sep="\t")

    try:
        from altanalyze3.components.visualization.marker_heatmap_h5ad import (
            generate_marker_heatmap_from_adata)
        mf = os.path.join(outdir, "MarkerFinder")
        os.makedirs(mf, exist_ok=True)
        generate_marker_heatmap_from_adata(
            adata, cluster_key="harmonized_state",
            out=os.path.join(mf, "harmonized_marker_heatmap.pdf"),
            top_n=top_n,
            markers_tsv=os.path.join(mf, "harmonized_marker_heatmap_markers.tsv"),
            heatmap_tsv=os.path.join(mf, "harmonized_marker_heatmap_fold_matrix.tsv"),
            heatmap_cache=os.path.join(mf, "harmonized_marker_heatmap_fold_matrix.npz"),
            marker_method="markerfinder", cells_per_cluster=heatmap_cells, seed=0,
            species="human", covariate_columns=["source_dataset"],
            layer=("counts" if counts is not None else None),
            scale_data=bool(counts is not None))
        _log(f"heatmap PDF written to {mf}/harmonized_marker_heatmap.pdf")
    except Exception as exc:
        _log(f"WARNING heatmap PDF failed ({type(exc).__name__}: {exc}); "
             f"matrices and centroids were still written")

# ----------------------------------------------------------------------------------- assembly
def _safe_min_fdr(values: Sequence[float]) -> float:
    """Smallest FDR among members. Returns NaN when no member carries a GO-Elite name,
    instead of letting np.nanmin warn on an all-NaN list."""
    finite = [v for v in values if v == v]
    return float(min(finite)) if finite else float("nan")


def _pairs(spec: Optional[str]) -> set:
    """Parse a DS|CLUSTER comma list. A cluster name may hold a comma, so split on '|' first."""
    if not spec:
        return set()
    out, buf = set(), []
    for tok in str(spec).split(","):
        if "|" in tok:
            if buf:
                out.add(",".join(buf))
            buf = [tok]
        elif buf:
            buf.append(tok)
    if buf:
        out.add(",".join(buf))
    # A stray or trailing comma would otherwise survive as part of a cluster name and
    # match nothing. The caller then validates every token against the real clusters.
    return {x.strip().strip(",").strip() for x in out if x.strip().strip(",").strip()}


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=PATH",
                    help="One completed ICGS3 output directory. Repeat for each dataset.")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--seed-dataset", default=None, metavar="NAME",
                    help="Dataset that initialises the reference, named as in --run. The seed "
                         "is never tested and never loses a cluster, so it anchors the whole "
                         "result. Default: the dataset with the most clusters, because the "
                         "finest partition makes every later cluster face the most demanding "
                         "reference. Ignored when --dataset-order is given.")
    ap.add_argument("--dataset-order", default=None, metavar="A,B,C",
                    help="Explicit entry order, comma separated, naming every dataset given "
                         "with --run. Overrides --seed-dataset and the default sort. Entry "
                         "order changes the result: a dataset entering last faces the largest "
                         "reference and contributes fewest states.")
    ap.add_argument("--cells-per-cluster", type=int, default=200,
                    help="Representative cells taken per source cluster at input, ranked by "
                         "SVM score. Default 200.")
    ap.add_argument("--require-multi-dataset", action="store_true",
                    help="Keep only states seen in two or more datasets. Off by default, so a "
                         "state unique to one dataset survives and is reported as such.")
    ap.add_argument("--marker-top-n", type=int, default=60,
                    help="Markers MarkerFinder reports per cluster in every internal call. "
                         "Default 60.")
    ap.add_argument("--nomination-query-top", type=int, default=60,
                    help="Top markers of a candidate tested for enrichment. ICGS3 reports "
                         "60 markers per cluster by default, so the query matches that. "
                         "Default 60.")
    ap.add_argument("--nomination-min-overlap", type=int, default=10,
                    help="Shared genes required to call a candidate redundant. Against a "
                         "60-marker query, 10 shared genes is the working expectation. "
                         "Default 10.")
    ap.add_argument("--nomination-fdr", type=float, default=0.05,
                    help="BH-adjusted enrichment FDR a nominee must reach against an existing "
                         "state before it counts as that state. Default 0.05.")
    ap.add_argument("--nomination-specificity", type=float, default=0.0,
                    help="How much better a cluster must match its best reference state than "
                         "its second best, as a fraction of the best overlap. 0 disables the "
                         "condition and is the default: an ablation on four lung datasets "
                         "(2026-08-21) showed 0.30 held apart 10 cross-study copies of the same "
                         "population, 8 of them above the 98th percentile of all cross-state "
                         "centroid correlations, and split the ionocyte marker panel across two "
                         "states. Raise it only with evidence that a specific reference is "
                         "over-merging. Default 0.0.")
    ap.add_argument("--nomination-overlap-fraction", type=float, default=0.17,
                    help="Overlap a cluster must share with a reference state, as a fraction of "
                         "the smaller of the two marker sets, used when that set is smaller "
                         "than the query. 0.17 corresponds to 10 genes out of 60. Default 0.17.")
    ap.add_argument("--survival-ref-cells", type=int, default=60,
                    help="Reference cells per state used in the one-at-a-time survival test. "
                         "Only marker structure matters there, so subsampling keeps the test "
                         "affordable. Default 60.")
    ap.add_argument("--damage-floor", type=int, default=1,
                    help="Unique markers an existing reference state must retain when a new "
                         "cluster is admitted. Equal to --survival-min-markers protects the "
                         "reference completely; 1 forbids only erasure. Default 1.")
    ap.add_argument("--survival-mode", choices=["report", "enforce"], default="report",
                    help="report: the final global MarkerFinder records states with too few "
                         "unique markers but keeps them, because winner-take-all assignment "
                         "starves real states at this scale. enforce: remove them.")
    ap.add_argument("--survival-rho", type=float, default=0.3,
                    help="Pearson r a marker must reach in the survival test. Default 0.3.")
    ap.add_argument("--survival-min-markers", type=int, default=1,
                    help="Unique markers a candidate must own at --survival-rho to count as a "
                         "distinct population. The purpose of this step is to admit valid "
                         "populations rather than remove them, so a single marker above the "
                         "correlation threshold is sufficient. Default 1.")
    ap.add_argument("--final-cells-per-state", type=int, default=100,
                    help="Cells per surviving state in the final MarkerFinder, drawn from "
                         "the input pool. Default 100.")
    ap.add_argument("--no-rebuild", action="store_true",
                    help="Skip the pooled MarkerFinder and heatmap matrix.")
    ap.add_argument("--cells-from", default="auto", choices=["auto", "result", "source"],
                    help="Which h5ad supplies the representative cells. 'result' uses each run's "
                         "icgs3_result.h5ad, which a streamed ICGS3 load truncates to the "
                         "selected variable features. 'source' uses the input h5ad named in each "
                         "run's icgs3_config.json, which carries every gene and the counts "
                         "layer. 'auto', the default, keeps the result h5ad for every run that "
                         "spans at least 90%% of the shared gene space and falls back to the "
                         "source only for the runs that do not.")
    ap.add_argument("--species", default="Hs", choices=["Hs", "Mm", "human", "mouse"],
                    help="Species for the GO-Elite BioMarkers reference. Default Hs (human).")
    ap.add_argument("--biomarker-file", default=None,
                    help="GO-Elite Ensembl-BioMarkers.txt to use. Default: the file ICGS3 uses "
                         "for --species.")
    ap.add_argument("--annotation-columns", default=None, metavar="A,B",
                    help="Reference annotation columns of the source runs to test, comma "
                         "separated, for example HLCA. Default: every non-technical column each "
                         "run carries. Columns are never pooled across vocabularies.")
    ap.add_argument("--annotation-min-evidence", type=int, default=4,
                    help="Matching cells (source 1) or genes (source 2) an annotation needs "
                         "before the module will assign it. Default 4.")
    ap.add_argument("--annotation-fdr", type=float, default=0.05,
                    help="BH-adjusted FDR an annotation must reach. Default 0.05.")
    ap.add_argument("--nomination-identity-top", type=int, default=10,
                    help="How many of a candidate's highest-correlation markers define its "
                         "identity. Default 10.")
    ap.add_argument("--nomination-identity-min", type=int, default=5,
                    help="How many of those identity markers must also mark the reference state "
                         "before the candidate is called redundant with it. A shared overlap of "
                         "generic genes is not evidence of one population. 0 restores the "
                         "overlap-only behaviour. Default 5.")
    ap.add_argument("--redundancy-top-n", type=int, default=0,
                    help="Second redundancy test on the head of each marker list. A candidate "
                         "already called redundant is kept unless its top N markers and the "
                         "target state's top N markers share at least "
                         "--redundancy-min-fraction of N genes. 0 disables the test.")
    ap.add_argument("--redundancy-min-fraction", type=float, default=0.30,
                    help="Fraction of --redundancy-top-n that the two head marker lists must "
                         "share for the exclusion to stand. Default 0.30.")
    ap.add_argument("--redundancy-centroid-r", type=float, default=0.0,
                    help="Additional condition for redundancy: the candidate's mean expression "
                         "centroid must correlate with the matched target at least this much "
                         "over the shared gene space. 0 disables the condition.")
    ap.add_argument("--compare-to-excluded", action="store_true",
                    help="Test every candidate against the clusters already set aside as well "
                         "as against the retained states. A candidate redundant with an "
                         "excluded cluster joins the excluded set.")
    ap.add_argument("--legacy-damage-test", action="store_true",
                    help="Let any reference state veto a candidate that would take its last "
                         "unique marker, whatever the relative evidence. The default instead "
                         "requires the vetoing state to hold at least as many unique markers as "
                         "the candidate, because a state holding 1 weak marker would otherwise "
                         "exclude a distinct cell type holding 60 strong ones.")
    ap.add_argument("--force-admit", default=None, metavar="DS|CLUSTER,...",
                    help="Clusters that must become their own state even when gate 1 or gate 2 "
                         "would exclude them. Use with --suppress-cluster to let a finer "
                         "dataset REPLACE a coarser reference state. "
                         "ICGS_integrate_subcluster.py derives both lists from a first pass, "
                         "so the replacement runs through this entry point and produces the "
                         "full standard output set rather than a hand-edited state table.")
    ap.add_argument("--suppress-cluster", default=None, metavar="DS|CLUSTER,...",
                    help="Clusters that must never become a state. Names the coarse cluster a "
                         "finer dataset replaces.")
    ap.add_argument("--no-annotation", action="store_true",
                    help="Skip cell state annotation. Both sources run by default.")
    args = ap.parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)
    runs: List[RunInputs] = []
    for spec in args.run:
        if "=" not in spec:
            raise ValueError(f"--run expects NAME=PATH, got {spec}")
        name, path = spec.split("=", 1)
        runs.append(load_run(name, path, marker_top_n=args.marker_top_n))
    if len(runs) < 2:
        raise ValueError("give at least two runs")
    # fail on a bad --seed-dataset or --dataset-order now, before the pooled matrix is built
    _requested_order = resolve_dataset_order(
        [r.name for r in runs], {r.name: len(r.clusters) for r in runs},
        seed_dataset=args.seed_dataset,
        dataset_order=(args.dataset_order.split(",") if args.dataset_order else None))
    _log("dataset entry order: " + " -> ".join(_requested_order))

    reps = pd.concat([select_representative_cells(r, args.cells_per_cluster) for r in runs],
                     ignore_index=True)
    reps.to_csv(os.path.join(args.output_dir, "representative_cells.tsv"), sep="\t", index=False)
    _log(f"representative cells: {reps.shape[0]} across {reps['dataset'].nunique()} datasets")

    # STEP 1, before any comparison: pull those cells out of the downsampled h5ads and restrict
    # to the union of the MarkerFinder genes. Every later step reads this one matrix.
    genes = union_feature_space(runs)
    _log(f"union MarkerFinder gene space: {len(genes)} genes "
         f"({', '.join(f'{r.name} {r.centroids.shape[0]}' for r in runs)})")
    pool, pool_barcodes, pool_datasets, pool_counts = expression_for_representative_cells(
        runs, reps, genes, cells_from=args.cells_from)
    bc_to_cluster = dict(zip(reps["dataset"].astype(str) + "|" + reps["barcode"].astype(str),
                             reps["ICGS3_cluster"].astype(str)))
    pool_clusters = [bc_to_cluster.get(b, "") for b in pool_barcodes]
    keep = [i for i, c in enumerate(pool_clusters) if c]
    pool = pool[keep]
    if pool_counts is not None:
        pool_counts = pool_counts[keep]
    pool_barcodes = [pool_barcodes[i] for i in keep]
    pool_datasets = [pool_datasets[i] for i in keep]
    pool_clusters = [pool_clusters[i] for i in keep]
    _log(f"pooled matrix: {pool.shape[0]} cells x {pool.shape[1]} union genes, "
         f"{pool.nnz / (pool.shape[0] * pool.shape[1]):.1%} dense")

    # STEP 2: rebuild every cluster centroid from that shared matrix
    recomputed = recompute_centroids_from_pool(pool, genes, pool_datasets, pool_clusters)
    for r in runs:
        if r.name in recomputed:
            cols = [c for c in r.clusters if c in recomputed[r.name].columns]
            r.centroids = recomputed[r.name][cols]

    if True:
        node_of_all = [f"{d}|{c}" for d, c in zip(pool_datasets, pool_clusters)]
        states_list, audit, nomination = integrate_hierarchically(
            runs, marker_top_n=args.marker_top_n,
            outdir=args.output_dir, pool=pool, genes=genes, barcodes=pool_barcodes,
            node_of=node_of_all, survival_rho=args.survival_rho,
            survival_min_markers=args.survival_min_markers,
            nomination_query_top=args.nomination_query_top,
            nomination_min_overlap=args.nomination_min_overlap,
            nomination_fdr=args.nomination_fdr,
            nomination_specificity=args.nomination_specificity,
            damage_floor=args.damage_floor,
            nomination_overlap_fraction=args.nomination_overlap_fraction,
            survival_ref_cells=args.survival_ref_cells,
            seed_dataset=args.seed_dataset,
            dataset_order=(args.dataset_order.split(",") if args.dataset_order else None),
            nomination_identity_top=args.nomination_identity_top,
            nomination_identity_min=args.nomination_identity_min,
            redundancy_top_n=args.redundancy_top_n,
            redundancy_min_fraction=args.redundancy_min_fraction,
            redundancy_centroid_r=args.redundancy_centroid_r,
            compare_to_excluded=bool(args.compare_to_excluded),
            force_admit=_pairs(args.force_admit),
            suppress_cluster=_pairs(args.suppress_cluster),
            legacy_damage_test=bool(args.legacy_damage_test))
        # the order that actually ran, read off the audit rather than recomputed
        _resolved_order = list(dict.fromkeys(audit["step"].astype(str).tolist()))
        states = states_to_frame(states_list, runs)
        if args.require_multi_dataset:
            dropped = states.loc[states["n_supporting_datasets"] < 2]
            _log(f"--require-multi-dataset: dropping {dropped.shape[0]} single-dataset states "
                 f"holding {int(dropped['n_cells'].sum())} cells")
            states = states.loc[states["n_supporting_datasets"] >= 2].copy()
        states.to_csv(os.path.join(args.output_dir, "harmonized_states.tsv"), sep="\t", index=False)
        n_rep = int((states["kind"] == "supported_by_multiple_datasets").sum())
        n_uni = int((states["kind"] == "single_dataset_only").sum())
        _log(f"reference states: {states.shape[0]} total, {n_rep} corroborated by another "
             f"dataset, {n_uni} seen in one dataset only")
        cen = pd.DataFrame({f"S{st.sid}": st.centroid for st in states_list
                            if f"S{st.sid}" in set(states["state"])}).astype(np.float32)
        cen.to_csv(os.path.join(args.output_dir, "harmonized_centroids.tsv"), sep="\t")
        _log(f"harmonized centroids: {cen.shape[0]} genes x {cen.shape[1]} states")
        bad = states.loc[states["members"].map(
            lambda m: len({x.split("|")[0] for x in str(m).split(",")}) != len(str(m).split(",")))]
        _log(f"invariant check, one cluster per dataset per state: "
             f"{'PASS' if bad.empty else f'FAIL on {bad.shape[0]} states'}")

        # 2. HOPACH orders the harmonized centroids before anything downstream reads them
        keep_ids = set(states["state"])
        ordered_states = [st for st in states_list if f"S{st.sid}" in keep_ids]
        order, hop = hopach_order_states(ordered_states, genes)
        hop.to_csv(os.path.join(args.output_dir, "hopach_state_order.tsv"), sep="\t", index=False)
        states = states.merge(hop, on="state", how="left").sort_values("hopach_position")
        states.to_csv(os.path.join(args.output_dir, "harmonized_states.tsv"), sep="\t", index=False)
        cen = cen[[c for c in hop["state"] if c in cen.columns]]
        cen.to_csv(os.path.join(args.output_dir, "harmonized_centroids.tsv"), sep="\t")

        # STEPS 5-6: MarkerFinder, survival test at the required rho, drop and repeat
        if not args.no_rebuild:
            member_state = {}
            for _, row in states.iterrows():
                for m in str(row["members"]).split(","):
                    if m:
                        member_state[m] = row["state"]
            node_of = [f"{d}|{c}" for d, c in zip(pool_datasets, pool_clusters)]
            labels_all = [member_state.get(n, "") for n in node_of]
            rows0 = [i for i, l in enumerate(labels_all) if l]
            pool_k = pool[rows0]
            pool_counts_k = pool_counts[rows0] if pool_counts is not None else None
            bc_k = [pool_barcodes[i] for i in rows0]
            lab_k = [labels_all[i] for i in rows0]
            ds_k = [pool_datasets[i] for i in rows0]
            cl_k = [pool_clusters[i] for i in rows0]

            ds_of_state = {row["state"]: set(str(row["source_dataset"]).split(","))
                           for _, row in states.iterrows()}
            m_all, m_top, heat, lab_k, bc_k, pool_k, merged, dropped = markerfinder_survival(
                pool_k, genes, bc_k, lab_k, rho=args.survival_rho,
                min_markers=args.survival_min_markers, top_n=args.marker_top_n,
                centroids=cen, datasets_of_state=ds_of_state,
                max_rounds=(10 if args.survival_mode == "enforce" else 1),
                enforce=(args.survival_mode == "enforce"))
            merged.to_csv(os.path.join(args.output_dir, "states_merged_by_survival.tsv"),
                          sep="\t", index=False)
            dropped.to_csv(os.path.join(args.output_dir, "states_removed_by_survival.tsv"),
                           sep="\t", index=False)
            m_top.to_csv(os.path.join(args.output_dir, "harmonized_markers.tsv"),
                         sep="\t", index=False)
            m_all.to_csv(os.path.join(args.output_dir, "harmonized_markers_all.tsv"),
                         sep="\t", index=False)
            heat.to_csv(os.path.join(args.output_dir, "harmonized_heatmap_matrix.tsv"), sep="\t")
            # cells were re-subset inside the survival step, so rebuild the per-cell
            # dataset and source-cluster vectors from the surviving barcodes
            meta = dict(zip(pool_barcodes, zip(pool_datasets, pool_clusters)))
            ds_k = [meta[b][0] for b in bc_k]
            cl_k = [meta[b][1] for b in bc_k]
            score_lookup = dict(zip(
                reps["dataset"].astype(str) + "|" + reps["barcode"].astype(str),
                reps.get("ICGS3_SVM_score", pd.Series(0.0, index=reps.index))))
            survivors = sorted(set(lab_k))
            _log(f"survival test kept {len(survivors)} states: "
                 f"{merged.shape[0]} merged, {dropped.shape[0]} dropped, "
                 f"from {states.shape[0]} candidates")

            states = states.loc[states["state"].isin(survivors)].copy()
            states.to_csv(os.path.join(args.output_dir, "harmonized_states.tsv"),
                          sep="\t", index=False)
            cen = cen[[c for c in cen.columns if c in set(survivors)]]

            # STEP 7: HOPACH arranges the surviving states
            kept_objs = [st for st in states_list if f"S{st.sid}" in set(survivors)]
            if len(kept_objs) > 2:
                order2, hop2 = hopach_order_states(kept_objs, genes)
                hop2.to_csv(os.path.join(args.output_dir, "hopach_state_order.tsv"),
                            sep="\t", index=False)
                states = states.drop(columns=[c for c in ("hopach_position", "hopach_cluster")
                                              if c in states.columns]).merge(
                    hop2, on="state", how="left").sort_values("hopach_position")
                states.to_csv(os.path.join(args.output_dir, "harmonized_states.tsv"),
                              sep="\t", index=False)
                cen = cen[[c for c in hop2["state"] if c in cen.columns]]
            cen.to_csv(os.path.join(args.output_dir, "harmonized_centroids.tsv"), sep="\t")

            # STEP 8: final MarkerFinder on the top N cells per surviving state
            score_of = score_lookup
            per_state: Dict[str, List[int]] = {}
            for i in range(len(lab_k)):
                per_state.setdefault(lab_k[i], []).append(i)
            final_rows: List[int] = []
            for st, idxs in per_state.items():
                idxs.sort(key=lambda i: -float(score_of.get(bc_k[i], 0.0)))
                final_rows.extend(idxs[:int(args.final_cells_per_state)])
            final_rows.sort()
            _log(f"STEP 8: final MarkerFinder on top {args.final_cells_per_state} cells per state, "
                 f"{len(final_rows)} cells across {len(per_state)} states")
            name_of_state = dict(zip(states["state"].astype(str), states["name"].astype(str)))
            f_all, f_top, f_heat, _fl, _fb, _fp, _fm, _fd = markerfinder_survival(
                pool_k[final_rows], genes, [bc_k[i] for i in final_rows],
                [lab_k[i] for i in final_rows], rho=args.survival_rho,
                min_markers=args.survival_min_markers, top_n=args.marker_top_n, max_rounds=1)
            f_top.to_csv(os.path.join(args.output_dir, "final_markers.tsv"),
                         sep="\t", index=False)

            # STEP 9: annotate every state from both sources, then adopt the better result
            if not args.no_annotation:
                fin_states = [lab_k[i] for i in final_rows]
                fin_ds = [ds_k[i] for i in final_rows]
                fin_bc = [bc_k[i] for i in final_rows]
                want_cols = ([c.strip() for c in args.annotation_columns.split(",") if c.strip()]
                             if args.annotation_columns else None)
                src1 = reference_cell_type_enrichment(fin_states, fin_ds, fin_bc, runs,
                                                      columns=want_cols)
                src2 = marker_gene_set_enrichment(f_top, list(genes), species=args.species,
                                                  biomarker_file=args.biomarker_file,
                                                  outdir=args.output_dir)
                cand = pd.concat([x for x in (src1, src2) if not x.empty], ignore_index=True) \
                    if (not src1.empty or not src2.empty) else pd.DataFrame()
                if not cand.empty:
                    cand.sort_values(["state", "source", "fdr", "accuracy"],
                                     ascending=[True, True, True, False]).to_csv(
                        os.path.join(args.output_dir, "state_annotation_candidates.tsv"),
                        sep="\t", index=False)
                chosen = assign_state_annotations(
                    cand, sorted(set(fin_states)),
                    min_evidence=args.annotation_min_evidence, fdr=args.annotation_fdr)
                # the best call from each source separately, so a disagreement stays visible
                for tag, frame in (("reference_cells", src1), ("marker_gene_sets", src2)):
                    if frame.empty:
                        continue
                    b = assign_state_annotations(
                        frame, sorted(set(fin_states)),
                        min_evidence=args.annotation_min_evidence, fdr=args.annotation_fdr)
                    chosen = chosen.merge(
                        b[["state", "annotation", "annotation_accuracy", "annotation_evidence"]]
                        .rename(columns={"annotation": f"{tag}_annotation",
                                         "annotation_accuracy": f"{tag}_accuracy",
                                         "annotation_evidence": f"{tag}_evidence"}),
                        on="state", how="left")
                chosen.to_csv(os.path.join(args.output_dir, "state_annotations.tsv"),
                              sep="\t", index=False)
                n_named = int((chosen["annotation"].astype(str) != "").sum())
                by_src = chosen.loc[chosen["annotation"].astype(str) != "",
                                    "annotation_source"].value_counts().to_dict()
                _log(f"STEP 9 annotation: {n_named} of {chosen.shape[0]} states named "
                     f"(>= {args.annotation_min_evidence} matching cells or genes, "
                     f"FDR <= {args.annotation_fdr}); by source {by_src}")
                states = states.drop(columns=[c for c in chosen.columns
                                              if c != "state" and c in states.columns]) \
                    .merge(chosen, on="state", how="left")
                states.to_csv(os.path.join(args.output_dir, "harmonized_states.tsv"),
                              sep="\t", index=False)
                name_of_state = {s: (a if isinstance(a, str) and a else name_of_state.get(s, ""))
                                 for s, a in zip(chosen["state"].astype(str),
                                                 chosen["annotation"].astype(str))}

            write_integration_deliverables(
                pool_k[final_rows], genes, [bc_k[i] for i in final_rows],
                [lab_k[i] for i in final_rows], [ds_k[i] for i in final_rows],
                [cl_k[i] for i in final_rows], name_of_state, args.output_dir,
                top_n=args.marker_top_n, heatmap_cells=0,
                # pool_counts_k is already aligned with bc_k, so final_rows indexes it directly;
                # mapping through pool_barcodes indexed the full 39,396-cell pool instead
                counts=(pool_counts_k[final_rows] if pool_counts_k is not None else None))

        # `n_cells` came from the input draw (--cells-per-cluster) while the atlas exports
        # --final-cells-per-state cells per state, so the summary table reported up to twice the
        # cell support the per-cell table holds. Rewrite it from the file that was actually
        # written, and keep the draw size under its own name.
        # Every cluster that did not become a state, and the reference state it was judged
        # closest to. A reader renaming a state needs this: it says which cluster of which
        # dataset maps onto which reference state, and on what evidence.
        try:
            aud = audit.copy()
            st_tab = pd.read_csv(os.path.join(args.output_dir, "harmonized_states.tsv"), sep="\t")
            nm = dict(zip(st_tab["state"].astype(str), st_tab["name"].astype(str)))
            mem = dict(zip(st_tab["state"].astype(str), st_tab["members"].astype(str)))
            out = aud[aud["action"].astype(str).str.startswith("excluded")].copy()
            tgt = out["nearest_state_r025"].astype(str)
            out["reference_state"] = tgt
            out["reference_state_name"] = tgt.map(nm).fillna("")
            out["reference_state_members"] = tgt.map(mem).fillna("")
            out["excluded_cluster"] = (out["step"].astype(str) + "|" + out["cluster"].astype(str))
            keep = ["step", "cluster", "excluded_cluster", "action", "reference_state",
                    "reference_state_name", "reference_state_members", "overlap_r025",
                    "fdr_r025", "overlap_r030", "fdr_r030", "withdrawn_reason"]
            keep = [c for c in keep if c in out.columns]
            out = out[keep].rename(columns={"step": "dataset"})
            path = os.path.join(args.output_dir, "excluded_cluster_assignments.tsv")
            out.sort_values(["dataset", "cluster"]).to_csv(path, sep="\t", index=False)
            _log(f"wrote {path}: {out.shape[0]} clusters that did not become a state, each with "
                 f"the reference state it matched")
        except Exception as exc:
            _log(f"WARNING could not write excluded_cluster_assignments.tsv ({exc})")

        cells_file = os.path.join(args.output_dir, "harmonized_cell_annotations.tsv")
        if os.path.exists(cells_file):
            actual = pd.read_csv(cells_file, sep="\t")["harmonized_state"].astype(str).value_counts()
            path = os.path.join(args.output_dir, "harmonized_states.tsv")
            st = pd.read_csv(path, sep="\t")
            st["n_cells_representative"] = st["n_cells"]
            st["n_cells"] = st["state"].astype(str).map(actual).fillna(0).astype(int)
            st.to_csv(path, sep="\t", index=False)
            n_fixed = int((st["n_cells"] != st["n_cells_representative"]).sum())
            _log(f"n_cells rewritten from harmonized_cell_annotations.tsv for {n_fixed} of "
                 f"{st.shape[0]} states; the input draw is kept as n_cells_representative")

        summary = {
            "mode": "hierarchical",
            "runs": {r.name: {"clusters": len(r.clusters), "cells": int(r.cells.shape[0])} for r in runs},
            "seed": _resolved_order[0],
            "dataset_order": _resolved_order,
            "seed_dataset_override": args.seed_dataset,
            "dataset_order_override": args.dataset_order,
            "cells_per_cluster": args.cells_per_cluster,
            "require_multi_dataset": bool(args.require_multi_dataset),
            "nomination_query_top": args.nomination_query_top,
            "nomination_min_overlap": args.nomination_min_overlap,
            "nomination_overlap_fraction": args.nomination_overlap_fraction,
            "nomination_fdr": args.nomination_fdr,
            "nomination_specificity": args.nomination_specificity,
            "survival_rho": args.survival_rho,
            "survival_min_markers": args.survival_min_markers,
            "damage_floor": args.damage_floor,
            "survival_ref_cells": args.survival_ref_cells,
            "survival_mode": args.survival_mode,
            "final_cells_per_state": args.final_cells_per_state,
            "species": args.species,
            "annotation": (not args.no_annotation),
            "annotation_min_evidence": args.annotation_min_evidence,
            "annotation_fdr": args.annotation_fdr,
            "annotation_columns": args.annotation_columns,
            "cells_from": args.cells_from,
            "n_states": int(states.shape[0]),
            "n_states_repeated": n_rep,
            "n_states_dataset_unique": n_uni,
            "n_input_clusters": int(sum(len(r.clusters) for r in runs)),
            "invariant_one_cluster_per_dataset": bool(bad.empty),
        }
        with open(os.path.join(args.output_dir, "integration_summary.json"), "w") as handle:
            json.dump(summary, handle, indent=2)
        _log(f"wrote outputs to {args.output_dir}")
        return 0

if __name__ == "__main__":
    raise SystemExit(main())
