"""
Unit tests for the metacell construction module.

These tests exercise the major execution paths offered by the CLI, including
boundary-restricted aggregation, graph-based clustering, and random metacell
generation. The tests produce a detailed log containing runtime statistics and
matrix dimensionality summaries for each step to support rapid validation.
"""

from __future__ import annotations

import json
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Tuple
from unittest import mock

import numpy as np
import pandas as pd
from tqdm import tqdm

import anndata as ad

import altanalyze3.components.metacells.main as mc_main
from altanalyze3.components.tests.metacell_synthetic_data import build_synthetic_metacell_adata
from altanalyze3.components.metacells.main import (
    MetacellParams,
    aggregate_expression,
    assemble_metacells,
    build_membership_df,
    build_metacell_obs,
    generate_random_groups,
    group_iterable,
    initial_labels,
    process_group,
    resolve_boundary_columns,
    subset_adata,
)

# Type alias for test return payloads
TestResult = Dict[str, object]
TestFunc = Callable[[ad.AnnData], TestResult]


def make_params(**overrides) -> MetacellParams:
    """Helper to instantiate MetacellParams with sensible defaults."""
    defaults = dict(
        target_size=6,
        min_size=3,
        max_size=12,
        preserve_small=False,
        aggregation="sum",
        graph_algorithm="kmeans",
        resolution=None,
        resolution_steps=4,
        resolution_tolerance=0.25,
        n_neighbors=12,
        neighbor_method="auto",
        neighbor_metric="euclidean",
        n_top_genes=50,
        n_pcs=10,
        pca_svd_solver="randomized",
        hvg_layer=None,
        expression_layer=None,
        use_raw=False,
        random_state=7,
        random_metacell_count=None,
        random_cells_per_metacell=5,
        random_sampling_with_replacement=False,
    )
    defaults.update(overrides)
    return MetacellParams(**defaults)


def test_subset_restriction(adata: ad.AnnData) -> TestResult:
    restricted = subset_adata(
        adata,
        sample_column="sample_id",
        cell_type_column="cell_type",
        restrict_samples=["D1"],
        restrict_cell_types=["B_cell"],
    )
    assert restricted.n_obs == 8, "Expected 8 cells after restricting to D1/B_cell"
    return {
        "success": True,
        "message": "Subset restriction retained expected cells.",
        "input_shape": (adata.n_obs, adata.n_vars),
        "output_shape": (restricted.n_obs, restricted.n_vars),
    }


def test_group_iterable(adata: ad.AnnData) -> TestResult:
    groups = list(group_iterable(adata, ["sample_id", "cell_type"]))
    assert len(groups) == 9, "Should produce 9 sample/cell-type groups."
    all_sizes = [group.n_obs for _, group in groups]
    assert all(size == 8 for size in all_sizes), "Each subgroup should contain 8 cells."
    return {
        "success": True,
        "message": "Grouping by sample/cell-type produced 9 balanced subsets.",
        "input_shape": (adata.n_obs, adata.n_vars),
        "output_shape": (sum(all_sizes), adata.n_vars),
    }


def test_generate_random_groups(_: ad.AnnData) -> TestResult:
    rng = np.random.default_rng(11)
    groups = generate_random_groups(
        n_cells=30,
        count=6,
        cells_per_group=5,
        with_replacement=False,
        rng=rng,
    )
    assert len(groups) == 6, "Random group generator should return requested number of groups."
    assert all(len(g) == 5 for g in groups), "Each random group must contain 5 indices."
    flat = np.concatenate(groups)
    assert flat.max() < 30 and flat.min() >= 0, "Generated indices must respect bounds."
    return {
        "success": True,
        "message": "Random aggregation generated the expected group structure.",
        "input_shape": (30, 0),
        "output_shape": (len(groups), len(groups[0])),
    }


def test_process_group_random(adata: ad.AnnData) -> TestResult:
    group = subset_adata(
        adata,
        sample_column="sample_id",
        cell_type_column="cell_type",
        restrict_samples=["D2"],
        restrict_cell_types=["T_cell"],
    )
    params = make_params(
        graph_algorithm="random",
        random_metacell_count=5,
        random_cells_per_metacell=4,
        random_sampling_with_replacement=False,
    )
    metacells, membership, next_idx = process_group(
        group=group,
        params=params,
        metacell_start_idx=0,
        boundary_columns=["sample_id", "cell_type"],
        boundary_key=("D2", "T_cell"),
        sample_column="sample_id",
        cell_type_column="cell_type",
    )
    assert metacells.n_obs == 5, "Random mode should emit exactly 5 metacells."
    assert membership["metacell_id"].nunique() == 5
    assert next_idx == 5
    return {
        "success": True,
        "message": "Random metacell generation succeeded for D2/T_cell.",
        "input_shape": (group.n_obs, group.n_vars),
        "output_shape": (metacells.n_obs, metacells.n_vars),
    }


def test_process_group_kmeans(adata: ad.AnnData) -> TestResult:
    group = subset_adata(
        adata,
        sample_column="sample_id",
        cell_type_column="cell_type",
        restrict_samples=["D3"],
        restrict_cell_types=["B_cell"],
    )
    params = make_params(
        graph_algorithm="kmeans",
        target_size=4,
        min_size=3,
        max_size=6,
        n_top_genes=30,
        n_pcs=5,
    )
    metacells, membership, next_idx = process_group(
        group=group,
        params=params,
        metacell_start_idx=10,
        boundary_columns=["sample_id", "cell_type"],
        boundary_key=("D3", "B_cell"),
        sample_column="sample_id",
        cell_type_column="cell_type",
    )
    expected = math.ceil(group.n_obs / params.target_size)
    assert metacells.n_obs == expected, "KMeans mode produced unexpected metacell count."
    assert membership["metacell_id"].nunique() == expected
    assert next_idx == 10 + expected
    return {
        "success": True,
        "message": "KMeans metacell aggregation produced expected shapes.",
        "input_shape": (group.n_obs, group.n_vars),
        "output_shape": (metacells.n_obs, metacells.n_vars),
    }


def test_initial_labels_leiden(adata: ad.AnnData) -> TestResult:
    analysis = adata[:12].copy()
    analysis.obsm["X_pca"] = np.random.default_rng(3).normal(size=(analysis.n_obs, 4))
    params = make_params(graph_algorithm="leiden", target_size=4, resolution=None)

    def fake_leiden(adata_obj, resolution, random_state, key_added):
        labels = ["0" if i < 6 else "1" for i in range(adata_obj.n_obs)]
        adata_obj.obs[key_added] = labels

    with mock.patch.object(mc_main.sc.tl, "leiden", side_effect=fake_leiden):
        labels = initial_labels(analysis, params)
    assert len(np.unique(labels)) == 2, "Leiden path should emit two clusters via patched routine."
    return {
        "success": True,
        "message": "Leiden label initialization executed via patched function.",
        "input_shape": (analysis.n_obs, analysis.n_vars),
        "output_shape": (len(labels),),
    }


def test_initial_labels_louvain(adata: ad.AnnData) -> TestResult:
    analysis = adata[:9].copy()
    analysis.obsm["X_pca"] = np.random.default_rng(4).normal(size=(analysis.n_obs, 3))
    params = make_params(graph_algorithm="louvain", target_size=3, resolution=None)

    def fake_louvain(adata_obj, resolution, random_state, key_added):
        labels = ["0", "1", "0", "1", "2", "2", "0", "1", "2"]
        adata_obj.obs[key_added] = labels[: adata_obj.n_obs]

    with mock.patch.object(mc_main.sc.tl, "louvain", side_effect=fake_louvain):
        labels = initial_labels(analysis, params)
    assert len(np.unique(labels)) == 3, "Patched Louvain should produce three labels."
    return {
        "success": True,
        "message": "Louvain label initialization executed via patched function.",
        "input_shape": (analysis.n_obs, analysis.n_vars),
        "output_shape": (len(labels),),
    }


def test_aggregate_expression() -> TestResult:
    matrix = np.array([[1, 2], [3, 4], [5, 6]], dtype=float)
    groups = [np.array([0, 1]), np.array([2])]
    summed = aggregate_expression(matrix, groups, aggregation="sum")
    assert np.allclose(summed, np.array([[4, 6], [5, 6]]))
    return {
        "success": True,
        "message": "Aggregation produced expected summed outputs.",
        "input_shape": (matrix.shape[0], matrix.shape[1]),
        "output_shape": (len(groups), matrix.shape[1]),
    }


def test_build_metacell_obs() -> TestResult:
    groups = [np.array([0, 1, 2]), np.array([3, 4, 5])]
    metacell_ids = ["MC_000001", "MC_000002"]
    cell_obs = pd.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2", "S3", "S3", "S3"],
            "cell_type": ["T", "T", "T", "B", "B", "B"],
        }
    )
    metacell_obs = build_metacell_obs(
        groups=groups,
        metacell_ids=metacell_ids,
        cell_obs=cell_obs,
        sample_column="sample_id",
        cell_type_column="cell_type",
        boundary_info={"boundary": "demo"},
    )
    assert "dominant_sample_id" in metacell_obs.columns
    assert metacell_obs.loc["MC_000001", "dominant_sample_id"] == "S1"
    assert metacell_obs.loc["MC_000002", "dominant_cell_type"] == "B"
    return {
        "success": True,
        "message": "Metacell obs annotation computed dominant labels correctly.",
        "input_shape": (cell_obs.shape[0], 0),
        "output_shape": (metacell_obs.shape[0], metacell_obs.shape[1]),
    }


def test_build_membership_df() -> TestResult:
    groups = [np.array([0, 2]), np.array([1, 3])]
    metacell_ids = ["MC_A", "MC_B"]
    cell_names = np.array(["c0", "c1", "c2", "c3"])
    cell_obs = pd.DataFrame(
        {
            "sample_id": ["S1", "S2", "S1", "S2"],
            "cell_type": ["T", "B", "T", "B"],
        }
    )
    membership = build_membership_df(
        groups=groups,
        metacell_ids=metacell_ids,
        cell_names=cell_names,
        cell_obs=cell_obs,
        sample_column="sample_id",
        cell_type_column="cell_type",
    )
    assert membership.shape[0] == 4
    assert set(membership["metacell_id"]) == {"MC_A", "MC_B"}
    return {
        "success": True,
        "message": "Membership table captures all cell to metacell links.",
        "input_shape": (cell_names.size, 0),
        "output_shape": membership.shape,
    }


def test_assemble_metacells_global(adata: ad.AnnData) -> TestResult:
    params = make_params(
        graph_algorithm="kmeans",
        target_size=12,
        min_size=6,
        max_size=15,
        n_top_genes=40,
        n_pcs=8,
        aggregation="sum",
    )
    metacells, membership = assemble_metacells(
        adata,
        params=params,
        boundary_columns=[],
        sample_column="sample_id",
        cell_type_column="cell_type",
    )
    total_cells = membership.shape[0]
    assert total_cells == adata.n_obs, "Membership table should cover all cells."
    assert metacells.n_obs > 0, "Expected at least one metacell."
    return {
        "success": True,
        "message": "Global assembly succeeded using a kmeans backbone.",
        "input_shape": (adata.n_obs, adata.n_vars),
        "output_shape": (metacells.n_obs, metacells.n_vars),
    }


def test_resolve_boundary_columns() -> TestResult:
    class DummyArgs:
        boundary_columns = None
        mode = "boundary"
        sample_column = "sample_id"
        cell_type_column = "cell_type"

    columns = resolve_boundary_columns(DummyArgs())
    assert columns == ["sample_id", "cell_type"]
    return {
        "success": True,
        "message": "Boundary column resolution fell back to sample/cell_type.",
        "input_shape": None,
        "output_shape": len(columns),
    }


TEST_CASES: List[Tuple[str, TestFunc]] = [
    ("subset_restriction", test_subset_restriction),
    ("group_iterable", test_group_iterable),
    ("random_groups", test_generate_random_groups),
    ("process_group_random", test_process_group_random),
    ("process_group_kmeans", test_process_group_kmeans),
    ("initial_labels_leiden", test_initial_labels_leiden),
    ("initial_labels_louvain", test_initial_labels_louvain),
    ("aggregate_expression", lambda _: test_aggregate_expression()),
    ("build_metacell_obs", lambda _: test_build_metacell_obs()),
    ("build_membership_df", lambda _: test_build_membership_df()),
    ("assemble_metacells_global", test_assemble_metacells_global),
    ("resolve_boundary_columns", lambda _: test_resolve_boundary_columns()),
]


def _configure_logging() -> Path:
    """Prepare log file destination and return its path."""
    component_dir = Path(__file__).resolve().parents[1]
    log_dir = component_dir / "tests" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"metacells_{timestamp}.log"
    return log_path


def run_tests() -> None:
    log_path = _configure_logging()
    adata = build_synthetic_metacell_adata()

    # Ensure deterministic ordering for reproducibility
    adata = adata.copy()

    records: List[Dict[str, object]] = []
    failures: List[str] = []

    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write(
            json.dumps(
                {
                    "generated_at": datetime.now().isoformat(),
                    "input_cells": int(adata.n_obs),
                    "input_genes": int(adata.n_vars),
                }
            )
            + "\n"
        )

        for name, test_func in tqdm(TEST_CASES, desc="Metacell tests", unit="test"):
            start = time.perf_counter()
            try:
                result = test_func(adata)
                success = bool(result.get("success", True))
                message = str(result.get("message", ""))
            except AssertionError as err:
                success = False
                message = f"AssertionError: {err}"
                result = {}
            except Exception as exc:
                success = False
                message = f"Exception: {exc}"
                result = {}
            duration = time.perf_counter() - start

            entry = {
                "test": name,
                "success": success,
                "duration_sec": duration,
                "input_shape": result.get("input_shape"),
                "output_shape": result.get("output_shape"),
                "details": message,
            }
            handle.write(json.dumps(entry) + "\n")
            records.append(entry)
            if not success:
                failures.append(f"{name}: {message}")

    if failures:
        summary = "\n".join(failures)
        raise SystemExit(f"Metacell unit tests failed:\n{summary}\nLog: {log_path}")

    print(f"All metacell tests passed. Detailed log written to {log_path}")


if __name__ == "__main__":
    run_tests()
