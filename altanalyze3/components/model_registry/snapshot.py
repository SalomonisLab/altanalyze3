"""Explicitly refresh the model catalog from current files, without unpickling."""
import json
from datetime import datetime, timezone
import argparse
from pathlib import Path
from .registry import PACKAGE_ROOT, REGISTRY_PATH, describe_model, describe_application, _identity
from .static_defaults import api_selection
from .layout import LABEL_FIELDS, write_layout

API_MODULES = {
    "rna2lipid": "rna2lipid", "rna2adt": "rna2adt", "rna2grn": "rna2grn",
    "rna2metabolite": "rna2metabolite", "rna2lipid_aml": "rna2lipid/aml",
}

MODALITIES = {"lipids": "rna2lipid", "adt": "rna2adt", "grn": "rna2grn",
              "metabolite": "rna2metabolite", "lipid": "rna2lipid_aml"}


def build_catalog(observed_on):
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
        defaults.append({"application": application, "context": context, "component": component,
                         "model_version_id": mid, "observed_on": observed_on,
                         "source": "working-tree configuration", "effective_from": None,
                         "effective_until": None})
    catalog = {"schema_version": "1.0", "models": list(models.values()), "defaults": defaults,
               "application_methods": {name: describe_application(name) for name in
                                       ("scALABLE", "scALABLE-discover", "scALABLE-viewer")}}
    return catalog


def merge_catalog(previous, current):
    """Keep every historical artifact and inference identity when defaults move."""
    records = {record["model_version_id"]: dict(record) for record in previous.get("models", [])}
    active = {record["model_version_id"] for record in current["models"]}
    for record in records.values():
        if record["model_version_id"] not in active:
            record["status"] = "historical-default"
    version_keys = ("inference_version_id", "analysis_model_version_id", "inference_code_sha256")
    for record in current["models"]:
        old = records.get(record["model_version_id"], {})
        versions = list(old.get("inference_versions", []))
        if old and not versions:
            versions.append({key: old[key] for key in version_keys})
        new_version = {key: record[key] for key in version_keys}
        if new_version not in versions:
            versions.append(new_version)
        merged = {**old, **record, "inference_versions": versions}
        if "artifacts" in record:
            merged["artifacts"] = {}
            for role, artifact in record["artifacts"].items():
                prior = old.get("artifacts", {}).get(role, {})
                if prior.get("sha256") == artifact["sha256"]:
                    merged["artifacts"][role] = {**prior, **artifact}
                else:
                    merged["artifacts"][role] = dict(artifact)
        records[record["model_version_id"]] = merged
    result = {**current, "models": list(records.values())}
    by_id = {record["model_version_id"]: record for record in result["models"]}
    for default in result.get("defaults", []):
        for key in LABEL_FIELDS:
            if key in by_id[default["model_version_id"]]:
                default[key] = by_id[default["model_version_id"]][key]
    if previous.get("registry_release"):
        result["registry_release"] = previous["registry_release"]
    return result


def archive_catalog(catalog, destination):
    """Archive a complete catalog by its canonical JSON digest without overwriting."""
    path = destination / "snapshots" / ("catalog-" + _identity(catalog) + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if json.loads(path.read_text()) != catalog:
            raise ValueError("Archived catalog content differs from its identity")
    else:
        path.write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")
    return path.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-on", default=datetime.now(timezone.utc).date().isoformat())
    args = parser.parse_args()
    repository = PACKAGE_ROOT / "model_repository"
    previous = json.loads(REGISTRY_PATH.read_text()) if REGISTRY_PATH.exists() else {"models": [], "defaults": []}
    old_snapshot = archive_catalog(previous, repository)
    current = build_catalog(args.observed_on)
    catalog = merge_catalog(previous, current)
    new_snapshot = archive_catalog(catalog, repository)
    catalog["latest_snapshot"] = new_snapshot
    text = json.dumps(catalog, indent=2, sort_keys=True) + "\n"
    REGISTRY_PATH.write_text(text)
    (repository / "catalog.json").write_text(text)
    if all(record.get("registry_path") for record in catalog["models"]):
        write_layout(catalog, repository)
    print(f"Recorded {len(current['models'])} current and {len(catalog['models'])} total model versions; "
          f"{len(current['defaults'])} defaults. Prior snapshot: {old_snapshot}")


if __name__ == "__main__":
    main()
