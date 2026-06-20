"""Per-junction isoform IMPORTANCE (Goals 1 & 2) -- a run-once, expression-based
alternative to the length-based isoform pick.

For each junction, of the isoforms that CONTAIN it, report a 0-1 importance (the
isoform's share of that junction's reads) and the top (representative) isoform. This
replaces the current length-based `Isoform_1|Length` / `Isoform_2|Length` choice with
an expression-weighted one, and adds the full 0-1 list.

Method (benchmark winner -- see code-dev importance_benchmark): the ROBUST per-pseudobulk
read share. For junction j and isoform i (containing j), in cluster c:
    importance_{i,j,c} = mean over pseudobulks p in c (with junction reads >= min_reads)
                         of   Y_{i,p} / sum_{k contains j} Y_{k,p}
Averaging per-pseudobulk shares (rather than pooling totals) avoids the bias that pools
introduce toward high-variance isoforms; the min-reads filter drops noisy sparse
pseudobulks. It is clique-free (a junction's read share needs only its containing
isoforms, NOT the competing/background junctions), so no exon database is required.

This is the STATIC association importance; it is NOT the differential change-driver
(that is isoform_attribute.py / method M4). PSI-weighted and PSI-correlation variants
were benchmarked and are worse for this read-share association, so plain read share is used.

Originals are never modified -- a new '*-importance.txt' is written.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _structure_junctions(tokens: List[str]) -> List[str]:
    return [f"{tokens[i]}-{tokens[i + 1]}" for i in range(len(tokens) - 1)]


def _strip_gene(j: str) -> str:
    return j.split(":", 1)[1] if ":" in j else j


def _strip_version(x: str) -> str:
    return x.split(".", 1)[0].strip()


def load_junction_isoforms(assoc_path: str, genes_needed: set) -> Dict[str, Dict[str, List[str]]]:
    """gene -> {junction_key -> [isoform_id,...]} for the requested genes."""
    out: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    with open(assoc_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] not in genes_needed:
                continue
            gene, structure, iso = p[0], p[2], p[3]
            toks = [t for t in structure.split("|") if t]
            for jn in _structure_junctions(toks):
                out[gene][jn].append(iso)
    return out


class IsoformMatrix:
    """Pseudobulk x isoform counts, grouped by cluster (obs name '<cluster>.<sample>...')."""

    def __init__(self, h5ad_path: str):
        import anndata as ad
        self.ad = ad.read_h5ad(h5ad_path)
        self.col = {self._iso_key(v): i for i, v in enumerate(self.ad.var_names)}
        self.cluster_rows: Dict[str, List[int]] = defaultdict(list)
        for r, name in enumerate(self.ad.obs_names):
            self.cluster_rows[self._cluster(name)].append(r)
        self._X = None

    @staticmethod
    def _iso_key(v: str) -> str:
        # var names are 'ENSG:isoform_id' -> key on 'gene:iso'
        return str(v)

    @staticmethod
    def _cluster(obs_name: str) -> str:
        # 'MPP-1.WF40_ND21-251_HSC_3k-isoform' -> 'MPP-1'
        return str(obs_name).split(".", 1)[0]

    def X(self) -> np.ndarray:
        if self._X is None:
            X = self.ad.X
            self._X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        return self._X

    def col_for(self, gene: str, iso: str) -> Optional[int]:
        return self.col.get(f"{gene}:{iso}")


def importance_for_junction(
    gene: str, iso_ids: List[str], rows: List[int], mat: IsoformMatrix,
    min_reads: float = 5.0,
) -> List[Tuple[str, float]]:
    """Robust per-pseudobulk read share of each containing isoform, restricted to the
    given pseudobulk rows (a cluster). Returns [(iso_id, importance)] summing to ~1."""
    cols, keep_iso = [], []
    for iso in iso_ids:
        c = mat.col_for(gene, iso)
        if c is not None:
            cols.append(c); keep_iso.append(iso)
    if not cols:
        return []
    Y = mat.X()[np.ix_(rows, cols)]            # (n_pb_c x n_iso_j)
    N = Y.sum(1)                               # junction reads per pseudobulk
    ok = N >= min_reads
    if ok.sum() == 0:                          # fallback: pool all reads
        tot = Y.sum(0)
        s = tot.sum()
        imp = tot / s if s > 0 else np.zeros_like(tot)
    else:
        shares = Y[ok] / N[ok][:, None]        # per-pseudobulk share
        imp = shares.mean(0)
        si = imp.sum()
        imp = imp / si if si > 0 else imp
    order = np.argsort(-imp)
    return [(keep_iso[i], float(imp[i])) for i in order]


def resolve_cluster(path: str, known_clusters) -> Optional[str]:
    """Resolve the cluster from an event filename by matching the matrix's known cluster
    names as a suffix of '<cond1>-<cond2>-<cluster>' (handles cluster/condition names that
    themselves contain dashes). Returns the longest matching known cluster, else None."""
    base = os.path.basename(path).replace("_stats-annotated.txt", "").replace("_stats.txt", "")
    matches = [c for c in known_clusters if base == c or base.endswith("-" + c)]
    return max(matches, key=len) if matches else None


def annotate_event_table(event_tsv: str, assoc_path: str, mat: IsoformMatrix,
                         cluster: str, top_k: int = 5, output_tsv: Optional[str] = None,
                         min_reads: float = 5.0) -> pd.DataFrame:
    events = pd.read_csv(event_tsv, sep="\t", dtype=str).fillna("")
    rows_c = mat.cluster_rows.get(cluster)
    if not rows_c:
        raise ValueError(f"cluster '{cluster}' not found in isoform matrix "
                         f"(have e.g. {list(mat.cluster_rows)[:5]})")
    # genes present in the events
    genes = set()
    for feat in events.get("Feature", []):
        if "|" in feat:
            for half in feat.split("|"):
                if ":" in half:
                    genes.add(half.split(":", 1)[0])
    jiso = load_junction_isoforms(assoc_path, genes)

    def annotate_side(feat_half: str) -> Tuple[str, str, str]:
        if ":" not in feat_half:
            return "", "", ""
        gene = feat_half.split(":", 1)[0]
        jkey = _strip_gene(feat_half.strip())
        iso_ids = jiso.get(gene, {}).get(jkey, [])
        if not iso_ids:
            return "", "", "no_isoforms_in_models"
        ranked = importance_for_junction(gene, iso_ids, rows_c, mat, min_reads=min_reads)
        if not ranked:
            return "", "", "not_expressed"
        top = ranked[0][0]
        lst = ";".join(f"{i}:{v:.3f}" for i, v in ranked[:top_k])
        return top, lst, ""

    top1, list1, top2, list2, status = [], [], [], [], []
    for _, row in events.iterrows():
        feat = row.get("Feature", "")
        if "|" in feat:
            j1, j2 = feat.split("|", 1)
            t1, l1, s1 = annotate_side(j1)
            t2, l2, s2 = annotate_side(j2)
        else:
            t1 = l1 = t2 = l2 = ""; s1 = s2 = "no_feature"
        top1.append(t1); list1.append(l1); top2.append(t2); list2.append(l2)
        status.append(s1 or s2)

    events["imp_top_isoform_1"] = top1
    events["imp_list_1"] = list1
    events["imp_top_isoform_2"] = top2
    events["imp_list_2"] = list2
    events["imp_status"] = status

    if output_tsv is None:
        base, ext = os.path.splitext(event_tsv)
        output_tsv = f"{base}-importance{ext or '.txt'}"
    events.to_csv(output_tsv, sep="\t", index=False)
    n_ok = sum(1 for t in top1 if t)
    print(f"[isoform_importance] cluster={cluster}: annotated {n_ok}/{len(events)} junction-1 sides; "
          f"wrote {output_tsv}")
    return events


def annotate_differentials(work_dir: str, isoform_h5ad: Optional[str] = None,
                           assoc: Optional[str] = None, top_k: int = 5,
                           min_reads: float = 5.0) -> List[str]:
    """Workflow hook: annotate every dPSI event table with per-junction isoform
    importance (run-once, expression-based; clique-free; no exon DB). Writes companion
    'importance/dPSI-*-events/*-importance.txt'. Best-effort, non-destructive."""
    import glob
    isoform_h5ad = isoform_h5ad or os.path.join(work_dir, "isoform_combined_pseudo_cluster_counts.h5ad")
    assoc = assoc or os.path.join(work_dir, "gff-output", "transcript_associations.txt")
    if not (os.path.exists(isoform_h5ad) and os.path.exists(assoc)):
        print("[isoform_importance] isoform matrix or transcript_associations missing; skipping")
        return []
    mat = IsoformMatrix(isoform_h5ad)
    written = []
    for ev_dir in ("dPSI-covariate-events", "dPSI-cluster-events"):
        evp = os.path.join(work_dir, ev_dir)
        if not os.path.isdir(evp):
            continue
        out_dir = os.path.join(work_dir, "importance", ev_dir)
        os.makedirs(out_dir, exist_ok=True)
        for ev_file in sorted(glob.glob(os.path.join(evp, "*_stats-annotated.txt"))):
            cluster = resolve_cluster(ev_file, mat.cluster_rows.keys())
            if cluster is None:
                continue
            base = os.path.basename(ev_file).replace("_stats-annotated.txt", "-importance.txt")
            out = os.path.join(out_dir, base)
            try:
                annotate_event_table(ev_file, assoc, mat, cluster, top_k=top_k,
                                     output_tsv=out, min_reads=min_reads)
                written.append(out)
            except Exception as e:
                print(f"[isoform_importance] skipped {os.path.basename(ev_file)}: {e}")
    print(f"[isoform_importance] importance complete: {len(written)} tables -> {work_dir}/importance/")
    return written


def _build_parser():
    p = argparse.ArgumentParser(description="Per-junction isoform importance (Goals 1 & 2).")
    p.add_argument("event_tsv", help="dPSI event table with a Feature column.")
    p.add_argument("isoform_h5ad", help="isoform_combined_pseudo_cluster_counts.h5ad (pseudobulk x isoform).")
    p.add_argument("transcript_associations", help="gff-output/transcript_associations.txt")
    p.add_argument("--cluster", default=None, help="cluster name (else parsed from filename).")
    p.add_argument("--min-reads", type=float, default=5.0)
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("-o", "--output", default=None)
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    mat = IsoformMatrix(args.isoform_h5ad)
    cluster = args.cluster or resolve_cluster(args.event_tsv, mat.cluster_rows.keys())
    if cluster is None:
        raise ValueError(f"could not resolve cluster from {args.event_tsv}; pass --cluster")
    annotate_event_table(args.event_tsv, args.transcript_associations, mat, cluster,
                         top_k=args.top_k, output_tsv=args.output, min_reads=args.min_reads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
