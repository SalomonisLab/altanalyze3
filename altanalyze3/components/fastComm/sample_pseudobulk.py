"""Per-sample fastComm scores as a (sample x receiver state) pseudobulk for cellHarmony_differential.

Input is the split_scores_long.tsv that run_benchmark (benchmark-h5ad) writes: one row per scored interaction per
split (sample). Output is an h5ad with one row per (sample, receiver state) and one column per incoming interaction
named "<sender>|<ligand>-><receptor>"; X holds the score column (fastcomm_score by default). Keying rows on the
receiver state gives the rows the same cell_type vocabulary as the RNA pseudobulks. An interaction a sample did not
report is a structural 0 (fastComm emitted no edge: a state below min_cells or an LR pair below threshold).
uns["pseudobulk_method"] = "pseudobulk" routes cellHarmony_differential to its eBayes pseudobulk test.

python -m altanalyze3.components.fastComm.cli sample-pseudobulk --split-scores split_scores_long.tsv \\
    --output cell_communication_pseudobulk.h5ad [--covariates covariates.tsv --covariates-output cov.txt]
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

REQUIRED = ("split", "sender_state", "receiver_state", "ligand", "receptor")


def build_sample_pseudobulk(split_scores_long: Path, output: Path, score: str = "fastcomm_score",
                            sample_col: str = "Sample", population_col: str = "cell_type",
                            log=print) -> Tuple[object, Dict[str, object]]:
    import anndata as ad
    from scipy import sparse

    split_scores_long, output = Path(split_scores_long), Path(output)
    t = pd.read_csv(split_scores_long, sep="\t", usecols=lambda c: c in REQUIRED + (score,), dtype={"split": str})
    missing = [c for c in REQUIRED + (score,) if c not in t.columns]
    if missing:
        raise ValueError(f"{split_scores_long}: missing columns {missing}")
    n_records = len(t)
    if n_records == 0:
        raise ValueError(f"{split_scores_long}: no scored interactions")
    for c in ("split", "receiver_state"):
        bad = t[c].astype(str).str.contains("|", regex=False)
        if bad.any():
            raise ValueError(f"{c} values contain '|', which the row key uses: {sorted(t.loc[bad, c].unique())[:5]}")
    vals = pd.to_numeric(t[score], errors="coerce")
    if vals.isna().any():
        raise ValueError(f"{int(vals.isna().sum())} non-numeric {score} values")
    t[score] = vals
    t["row"] = t["split"].astype(str) + "|" + t["receiver_state"].astype(str)
    t["feature"] = (t["sender_state"].astype(str) + "|" + t["ligand"].astype(str) + "->"
                    + t["receptor"].astype(str))
    t = t.groupby(["row", "feature"], as_index=False, sort=False)[score].max()
    collapsed = n_records - len(t)
    log(f"[pseudobulk] {n_records:,} scored records; duplicate (sample|receiver, interaction) pairs collapsed by max: "
        f"{collapsed:,} ({100.0 * collapsed / n_records:.3f}%)")
    rows = pd.Index(sorted(t.row.unique()))
    cols = pd.Index(sorted(t.feature.unique()))
    X = sparse.csr_matrix((t[score].to_numpy(np.float32), (rows.get_indexer(t.row), cols.get_indexer(t.feature))),
                          shape=(len(rows), len(cols)), dtype=np.float32)
    X.eliminate_zeros()
    if X.nnz != int((t[score] != 0).sum()):
        raise AssertionError(f"stored values {X.nnz} != nonzero records {int((t[score] != 0).sum())}")
    obs = pd.DataFrame(index=rows)
    obs[sample_col] = [r.split("|", 1)[0] for r in rows]
    obs[population_col] = [r.split("|", 1)[1] for r in rows]
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=cols))
    a.uns["pseudobulk_method"] = "pseudobulk"
    a.uns["modality"] = "cell_communication"
    a.uns["value"] = score
    a.uns["source"] = str(split_scores_long.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(output, compression="gzip")
    summary = dict(source=str(split_scores_long.resolve()), output=str(output.resolve()), score=score,
                   records=n_records, collapsed_duplicates=collapsed, rows=len(rows), samples=int(obs[sample_col].nunique()),
                   receiver_states=int(obs[population_col].nunique()), interactions=len(cols), nnz=int(X.nnz),
                   density=float(X.nnz / max(1, X.shape[0] * X.shape[1])))
    log(f"[pseudobulk] wrote {output}: {len(rows):,} rows ({summary['samples']} samples x "
        f"{summary['receiver_states']} receiver states present) x {len(cols):,} interactions; "
        f"density {100 * summary['density']:.3f}%")
    return a, summary


def build_candidate_edges(split_scores_long: Path, h5ad: Path, output: Path, state_key: str, split_key: str,
                          min_cells: int, log=print) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Candidate edges for a per-sample differential: the union of edges any split reported in a thresholded
    benchmark-h5ad run, restricted to senders present (>= min_cells cells) in EVERY split. A sender absent from a split
    has no measured ligand there, so its edges would enter the pseudobulk as a stand-in 0; those edges are written to
    <output stem>_excluded_senders.tsv, never tested. Rescore the candidates with benchmark-h5ad --candidate-edges."""
    import h5py

    keys = ["sender_state", "receiver_state", "ligand", "receptor"]
    t = pd.read_csv(split_scores_long, sep="\t", usecols=lambda c: c in keys + ["split"], dtype=str)
    union = t[keys].drop_duplicates().reset_index(drop=True)
    with h5py.File(h5ad, "r") as h:
        obs = h["obs"]

        def column(k):
            g = obs[k]
            if isinstance(g, h5py.Group) and "categories" in g:
                cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
                codes = g["codes"][:]
                if (codes < 0).any():
                    raise ValueError(f"{h5ad}: obs[{k!r}] has missing values")
                return cats[codes]
            return np.array([v.decode() if isinstance(v, bytes) else v for v in g[:]], dtype=object)

        splits, states = column(split_key), column(state_key)
    counts = pd.crosstab(pd.Series(splits, name=split_key), pd.Series(states, name=state_key))
    present = counts >= min_cells
    in_all = set(counts.columns[present.all(axis=0)])
    reported = set(t["split"].astype(str))
    if reported != set(counts.index.astype(str)):
        log(f"[candidates] splits in scores {len(reported)} vs h5ad {counts.shape[0]}; "
            f"absent from scores: {sorted(set(counts.index.astype(str)) - reported)}")
    keep = union["sender_state"].isin(in_all)
    cand = union[keep].sort_values(keys).reset_index(drop=True)
    excl = (union[~keep].groupby("sender_state").size().rename("edges").reset_index()
            .assign(splits_with_min_cells=lambda d: d.sender_state.map(present.sum(axis=0)).fillna(0).astype(int),
                    splits_total=counts.shape[0]))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    cand.to_csv(output, sep="\t", index=False)
    excl.to_csv(output.with_name(output.stem + "_excluded_senders.tsv"), sep="\t", index=False)
    summary = dict(source=str(Path(split_scores_long).resolve()), h5ad=str(Path(h5ad).resolve()), min_cells=min_cells,
                   union_edges=int(len(union)), candidate_edges=int(len(cand)), excluded_edges=int((~keep).sum()),
                   senders_in_all_splits=len(in_all), states_total=int(counts.shape[1]),
                   excluded_senders=excl.sender_state.tolist(), splits=int(counts.shape[0]))
    log(f"[candidates] {len(union):,} edges reported in >= 1 split; {len(cand):,} kept (sender in all "
        f"{counts.shape[0]} splits with >= {min_cells} cells: {len(in_all)} of {counts.shape[1]} states); "
        f"{int((~keep).sum()):,} excluded from {len(excl)} senders")
    return cand, summary


def write_covariates(obs_samples, covariates: Path, output: Path, sample_col: str = "Sample",
                     covariate_col: str = "Condition", log=print) -> pd.DataFrame:
    """Covariates for cellHarmony_differential (Sample, Condition), restricted to the samples in the pseudobulk.
    Every pseudobulk sample must have exactly one covariate value."""
    cov = pd.read_csv(covariates, sep="\t", dtype=str)
    for c in (sample_col, covariate_col):
        if c not in cov.columns:
            raise ValueError(f"{covariates}: missing column {c}")
    cov = cov[[sample_col, covariate_col]].apply(lambda s: s.str.strip()).drop_duplicates()
    if cov[sample_col].duplicated().any():
        raise ValueError(f"{covariates}: a sample has more than one {covariate_col}")
    samples = sorted(set(obs_samples))
    absent = sorted(set(samples) - set(cov[sample_col]))
    if absent:
        raise ValueError(f"pseudobulk samples without a covariate: {absent}")
    out = cov[cov[sample_col].isin(samples)].sort_values(sample_col)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)
    log(f"[covariates] wrote {output}: {out[covariate_col].value_counts().to_dict()}")
    return out


def main_from_args(args) -> int:
    a, summary = build_sample_pseudobulk(args.split_scores, args.output, score=args.score,
                                         sample_col=args.sample_col, population_col=args.population_col)
    if args.covariates:
        out = args.covariates_output or Path(args.output).with_name(Path(args.output).stem + "_covariates.txt")
        write_covariates(a.obs[args.sample_col], args.covariates, out, sample_col=args.sample_col,
                         covariate_col=args.covariate_col)
        summary["covariates_output"] = str(Path(out).resolve())
    summary_path = args.summary_json or Path(args.output).with_suffix(".summary.json")
    Path(summary_path).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0
