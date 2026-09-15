"""Confirm scALABLE resolves the lung lipid modality to the lipid-wise bundle.

scALABLE (components/visualization/scalable_viewer) builds its application from
components/cellHarmony/webapp/app.create_app and computes through
components/cellHarmony/flask/pipeline. The lipid modality enters at
pipeline._build_imputed_lipid_adata. This script drives that function with each
lung reference entry, exactly as the job runner does, and reports which bundle
each reference resolves to.

The input matrix is synthetic. This is a wiring test, not a biological test: it
proves which model runs and what shape it returns, and nothing about accuracy.

Run:
  /usr/bin/python3 components/rna2lipid/validation/validate_scalable_wiring.py
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MODULE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = MODULE_DIR.parent.parent.parent
REGISTRY = REPO_ROOT / "altanalyze3" / "components" / "cellHarmony" / "flask" / "reference_config.json"


def _load_flask_pipeline():
    """Load components/cellHarmony/flask/pipeline.py without importing the
    package __init__, which pulls in the Flask web framework. The compute path
    scALABLE uses lives entirely in pipeline.py.
    """
    import importlib.util
    import types

    package_name = "altanalyze3.components.cellHarmony.flask"
    flask_dir = REGISTRY.parent
    if package_name not in sys.modules:
        package = types.ModuleType(package_name)
        package.__path__ = [str(flask_dir)]
        sys.modules[package_name] = package
    spec = importlib.util.spec_from_file_location(
        f"{package_name}.pipeline", flask_dir / "pipeline.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    sys.path.insert(0, str(REPO_ROOT))
    import anndata as ad

    from altanalyze3.components.rna2lipid.api import DEFAULT_BUNDLE_PATH, load_bundle

    flask_pipeline = _load_flask_pipeline()

    default_bundle = load_bundle()
    genes = list(default_bundle.input_genes)

    rng = np.random.default_rng(0)
    n_cells = 40
    # Real model gene symbols plus decoys, so gene matching is exercised.
    var_names = genes[:900] + [f"DECOY{i}" for i in range(100)]
    matrix = rng.lognormal(mean=1.0, sigma=0.5, size=(n_cells, len(var_names))).astype(np.float32)
    query = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"cell_state": ["AT2"] * 20 + ["Cap1"] * 20},
            index=pd.Index([f"cell_{i}" for i in range(n_cells)], dtype=str),
        ),
        var=pd.DataFrame(index=pd.Index(var_names, dtype=str)),
    )
    query.obsm["X_umap"] = rng.normal(size=(n_cells, 2))

    registry = json.loads(REGISTRY.read_text())
    lung_refs = [
        ref["id"]
        for species in registry["species"]
        for ref in species["references"]
        if "lipids" in (ref.get("impute_modalities") or [])
    ]

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "registry": str(REGISTRY),
        "rna2lipid_default_bundle": str(DEFAULT_BUNDLE_PATH),
        "default_bundle_architecture": default_bundle.architecture,
        "input": {
            "cells": n_cells,
            "var_names": len(var_names),
            "model_genes_present": 900,
            "decoy_genes": 100,
            "note": "Synthetic matrix. Wiring test only, no accuracy claim.",
        },
        "references": {},
    }

    for reference_id in lung_refs:
        entry = flask_pipeline._lookup_reference("human", reference_id, REGISTRY)
        lipid_adata, summary = flask_pipeline._build_imputed_lipid_adata(query, entry)
        report["references"][reference_id] = {
            "declared_bundle_path": (entry.get("impute_config") or {}).get("lipids", {}).get("bundle_path"),
            "resolved_bundle_path": summary.get("bundle_path"),
            "resolves_to_default": bool(
                Path(str(summary.get("bundle_path"))).resolve() == DEFAULT_BUNDLE_PATH.resolve()
            ),
            "architecture": summary.get("architecture"),
            "target_scaling_mode": summary.get("target_scaling_mode"),
            "matched_genes": summary.get("matched_genes"),
            "model_gene_count": summary.get("model_gene_count"),
            "output_lipids": int(lipid_adata.n_vars),
            "output_cells": int(lipid_adata.n_obs),
            "expression_scale": lipid_adata.uns.get("expression_scale"),
            "modality": lipid_adata.uns.get("modality"),
            "counts_layer_present": "counts" in lipid_adata.layers,
            "umap_carried_over": "X_umap" in lipid_adata.obsm,
            "pass": bool(
                summary.get("architecture") == "lipidwise"
                and int(lipid_adata.n_vars) == len(default_bundle.output_lipids)
                and int(lipid_adata.n_obs) == n_cells
            ),
        }

    report["all_pass"] = bool(
        lung_refs and all(v["pass"] for v in report["references"].values())
    )
    out_path = MODULE_DIR / "validation" / "scalable_wiring.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"\nWrote {out_path}")
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
