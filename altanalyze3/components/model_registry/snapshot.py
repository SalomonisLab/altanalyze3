"""Explicitly refresh the model catalog from current files, without unpickling."""
import json
from pathlib import Path
from .registry import PACKAGE_ROOT, REGISTRY_PATH, describe_model
from .static_defaults import api_selection

API_MODULES = {
    "rna2lipid": "rna2lipid", "rna2adt": "rna2adt", "rna2grn": "rna2grn",
    "rna2metabolite": "rna2metabolite", "rna2lipid_aml": "rna2lipid/aml",
}

MODALITIES = {"lipids": "rna2lipid", "adt": "rna2adt", "grn": "rna2grn",
              "metabolite": "rna2metabolite", "lipid": "rna2lipid_aml"}


def main():
    config_path = PACKAGE_ROOT / "components/cellHarmony/flask/reference_config.json"
    config = json.loads(config_path.read_text())
    selections = []
    for component, module in API_MODULES.items():
        api_path = "components/" + module + "/api.py"
        def read_source(path):
            source = PACKAGE_ROOT / path
            return source.read_text() if source.is_file() else None
        path = api_selection((PACKAGE_ROOT / api_path).read_text(), api_path, read_source)
        if not path:
            raise ValueError(f"Cannot resolve a static API default in {api_path}")
        selections.append((component, path, "altanalyze3", "api-default"))
    for species in config["species"]:
        for ref in species["references"]:
            for modality, cfg in ref.get("impute_config", {}).items():
                path = (config_path.parent / cfg["bundle_path"]).resolve().relative_to(PACKAGE_ROOT).as_posix()
                selections.append((MODALITIES[modality], path, "scALABLE", ref["id"]))
    # Resource models are versioned using BOTH actual files, just as at runtime.
    for species in ("human", "mouse"):
        folder = PACKAGE_ROOT / f"components/fastComm/training_data/processed/upstream/cellchat_nichenet_{species}"
        artifacts = {"ligand_receptor": folder / "ligand_receptor.tsv"}
        if (folder / "response_signatures.tsv").exists():
            artifacts["response_matrix"] = folder / "response_signatures.tsv"
        if artifacts["ligand_receptor"].exists():
            selections.append(("fastComm", artifacts, "altanalyze3/scALABLE", species))
    models = {}
    defaults = []
    for component, path, application, context in selections:
        paths = path if isinstance(path, dict) else {"bundle": PACKAGE_ROOT / path}
        identity = describe_model(component, paths, catalog_path=Path('/nonexistent-catalog'))
        mid = identity["model_version_id"]
        identity.pop("registry_status", None)
        models[mid] = {**identity, "name": ", ".join(p.name for p in paths.values()),
                       "status": "observed-default", "trained_at": None,
                       "training_provenance_status": "not-certified-by-registry",
                       "artifact_paths": {r: p.relative_to(PACKAGE_ROOT).as_posix() for r, p in paths.items()}}
        defaults.append({"application": application, "context": context,
                         "model_version_id": mid, "observed_on": "2026-10-06",
                         "source": "working-tree configuration", "effective_from": None,
                         "effective_until": None})
    catalog = {"schema_version": "1.0", "models": list(models.values()), "defaults": defaults}
    REGISTRY_PATH.write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")
    destination = PACKAGE_ROOT / "model_repository/catalog.json"
    destination.write_text(REGISTRY_PATH.read_text())
    print(f"Recorded {len(models)} model versions and {len(defaults)} default selections")


if __name__ == "__main__":
    main()
