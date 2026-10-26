"""fastComm candidate mode: rescore a fixed edge set in every split with thresholds off (benchmark --candidate-edges)."""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from altanalyze3.components.fastComm.benchmark import FastCommBenchmarkParams, run_benchmark
from altanalyze3.components.fastComm.sample_pseudobulk import build_candidate_edges

KEYS = ["sender_state", "receiver_state", "ligand", "receptor"]


def _inputs(tmp_path):
    rng = np.random.default_rng(3)
    genes = ["LIGA", "RECA", "LIGB", "RECB"] + [f"G{i}" for i in range(40)]
    rows, obs = [], []
    level = {"D1": 1.0, "D2": 0.35, "D3": 1.2}
    for donor in ("D1", "D2", "D3"):
        for state in ("S1", "S2", "S3"):
            n = 2 if (donor == "D3" and state == "S3") else 12
            for k in range(n):
                x = rng.gamma(1.0, 1.0, len(genes))
                x[:4] = 0.0
                if state == "S1":
                    # D2 expresses LIGA in 1 of 12 sender cells: its LIGA->RECA edges fall below the 0.2 threshold
                    x[0] = (0.02 if k == 0 else 0.0) if donor == "D2" else 40 * level[donor]
                    x[2] = 30
                if state == "S2":
                    x[1] = 40 * level[donor]
                if state == "S3":
                    x[3] = 30
                    x[2] = 25
                cp = x / x.sum() * 1e4
                rows.append(np.log1p(cp))
                obs.append((donor, state))
    X = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
    o = pd.DataFrame(obs, columns=["Donor", "cell_state"], index=[f"c{i}" for i in range(len(obs))])
    h5 = tmp_path / "toy.h5ad"
    ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=genes)).write_h5ad(h5)
    lr = tmp_path / "lr.tsv"
    pd.DataFrame({"ligand": ["LIGA", "LIGB"], "receptor": ["RECA", "RECB"], "source": "CellChatDB",
                  "pathway": ["PA", "PB"], "evidence_weight": 1.0}).to_csv(lr, sep="\t", index=False)
    return h5, lr


def _params(h5, lr, out, **kw):
    base = dict(h5ad=h5, output_dir=out, state_key="cell_state", split_key="Donor", lr_table=lr, response_matrix=None,
                lr_sources=["CellChatDB"], min_cells=5, min_lr_expression_score=0.2,
                max_lr_candidates_per_state_pair=5, include_self_edges=False)
    base.update(kw)
    return FastCommBenchmarkParams(**base)


def test_candidate_rescoring_recovers_below_threshold_values(tmp_path):
    h5, lr = _inputs(tmp_path)
    run_benchmark(_params(h5, lr, tmp_path / "thr"))
    thr = pd.read_csv(tmp_path / "thr" / "split_scores_long.tsv", sep="\t")
    assert len(thr) > 0
    cand, summary = build_candidate_edges(tmp_path / "thr" / "split_scores_long.tsv", h5, tmp_path / "cand.tsv",
                                          "cell_state", "Donor", 5, log=lambda *a: None)
    assert "S3" not in set(cand.sender_state), "S3 has 2 cells in D3, so it is not in every split"
    excl = pd.read_csv(tmp_path / "cand_excluded_senders.tsv", sep="\t")
    assert set(excl.sender_state) <= {"S3"} and summary["senders_in_all_splits"] == 2

    run_benchmark(_params(h5, lr, tmp_path / "re", candidate_edges=tmp_path / "cand.tsv"))
    re = pd.read_csv(tmp_path / "re" / "split_scores_long.tsv", sep="\t")
    got = set(map(tuple, re[KEYS].drop_duplicates().to_numpy()))
    assert got == set(map(tuple, cand[KEYS].to_numpy())), "rescoring returns exactly the candidate edges"
    per_split = re.groupby("split").size().to_dict()
    # a receiver below min_cells in a split has no row there (D3/S3: 2 cells); every other candidate edge is scored
    expected = {"D1": len(cand), "D2": len(cand), "D3": int((cand.receiver_state != "S3").sum())}
    assert per_split == expected, per_split

    both = thr.merge(re, on=["split"] + KEYS, suffixes=("_thr", "_re"))
    assert len(both) == len(thr[thr.sender_state.isin(cand.sender_state)])
    assert np.allclose(both.lr_expression_score_thr, both.lr_expression_score_re, rtol=0, atol=1e-12)
    low = re.merge(thr[["split"] + KEYS], on=["split"] + KEYS, how="left", indicator=True)
    low = low[low["_merge"] == "left_only"]
    assert len(low) > 0 and (low.lr_expression_score < 0.2).all(), "below-threshold values are now measured"
    assert (re.lr_expression_score >= 0).all()
