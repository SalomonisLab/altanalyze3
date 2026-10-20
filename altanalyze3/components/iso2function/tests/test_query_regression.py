"""Optional read-only parity check against the independently audited HSPC reference.

Set ISO2FUNCTION_QUERY_RESULTS to a standard CLI result directory. Dated artifacts are
validation references only; production queries never import them.
"""
import csv
import json
import os
from pathlib import Path

import pytest

from altanalyze3.components.iso2function import config


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def call_key(row):
    calls = tuple(sorted([(row["clone_A"], row["call_A"] if row["call_A"] in ["True", "False"] else "unknown"),
                          (row["clone_B"], row["call_B"] if row["call_B"] in ["True", "False"] else "unknown")]))
    return row["gene"], calls, row["modality"].lower(), row["target"]


def test_real_hspc_query_matches_independent_audit():
    supplied = os.environ.get("ISO2FUNCTION_QUERY_RESULTS")
    if not supplied:
        pytest.skip("Set ISO2FUNCTION_QUERY_RESULTS for the real-data parity check")
    out = Path(supplied)
    baseline = Path(config.COMPONENT_DIR).parents[1] / "output/hspc_tf_full_ppi_reconciliation_20261020"
    assert baseline.exists(), "Independent baseline artifacts required for this requested regression"
    manifest = json.loads((out / "query_manifest.json").read_text())
    assert manifest["status"] == "complete"
    reference_pairs = read(baseline / "ALL_EXACT_NATIVE_ASSAY_PAIR_COMPARISONS.tsv")
    reference_genes = {r["clone"]: r["ensg"] for r in json.loads((baseline / "RESOLVED_ALL_693_CLONES.json").read_text())}
    actual_genes = {r["clone"]: r["resolved_ensg"] for r in read(out / "clone_gene_audit.tsv")}
    assert actual_genes == reference_genes and len(actual_genes) == 693
    all_pairs = read(out / "all_pair_decisions.tsv")
    supported = [r for r in all_pairs if r["native_supported"] == "True"]
    key = lambda r: tuple(sorted([r["clone_A"], r["clone_B"]]))
    assert {key(r) for r in supported} == {key(r) for r in reference_pairs}
    assert len(supported) == 60 and len(all_pairs) == 803
    ref_reads = {(r["gene"], r["clone_" + s]): float(r["max_single_native_reads_" + s])
                 for r in reference_pairs for s in ["A", "B"]}
    assert all(float(r["max_single_native_reads_" + s]) == ref_reads[r["gene"], r["clone_" + s]]
               for r in supported for s in ["A", "B"])
    supported_keys = {key(r) for r in supported}
    calls = [r for r in read(out / "all_pair_target_calls.tsv") if key(r) in supported_keys]
    expected_calls = read(baseline / "ALL_EXACT_PAIR_TARGET_CALLS.tsv")
    assert {call_key(r) for r in calls} == {call_key(r) for r in expected_calls}
    divergent = read(out / "divergent_interactions.tsv")
    reference_ppi = read(baseline / "partners_gt1_CPM/EXPRESSED_DIFFERENTIAL_PPI_PARTNERS.tsv")
    reference_pdi = [r for r in expected_calls if r["modality"] == "PDI" and r["contrast"] in ["A_only", "B_only"]]
    assert {call_key(r) for r in divergent} == {call_key(r) for r in reference_ppi + reference_pdi}
    assert len(divergent) == 535
    assert sum(r["modality"] == "ppi" for r in divergent) == 342
    assert sum(r["modality"] == "pdi" for r in divergent) == 193
    expression = read(out / "partner_expression_audit.tsv")
    expected_expression = {r["target"]: r for r in read(baseline / "partners_gt1_CPM/ALL_330_PARTNER_EXPRESSION_AUDIT.tsv")}
    assert {r["target"] for r in expression} == set(expected_expression)
    exp = manifest["coverage"]["expression"]
    assert exp["donors"] == 9 and len(exp["selected_observations"]) == 36
    assert exp["cpm_verified_against_counts"]
    reference_cpm = {(r["state"], r["library"], r["ENSG"]): float(r["native_gene_CPM"])
                     for r in read(baseline / "partners_gt1_CPM/GENE_EXPRESSION_DONOR_STATE.tsv")}
    actual_cpm = {(r["state"], r["library"], r["ensg"]): float(r["cpm"])
                  for r in read(out / "gene_expression_donor_state.tsv")}
    assert actual_cpm == reference_cpm
