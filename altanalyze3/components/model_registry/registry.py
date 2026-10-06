from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = Path(__file__).with_name("catalog.json")
INFERENCE_FILES = {
    "rna2lipid": ["rna2lipid/api.py"],
    "rna2adt": ["rna2adt/api.py", "rna2adt/training.py", "rna2adt/rna2lipid_arch.py", "rna2adt/lung/model.py", "rna2adt/mouse/train_mouse.py", "rna2adt/generalizable.py", "rna2adt/centroid.py", "rna2adt/targeted.py"],
    "rna2grn": ["rna2grn/api.py", "rna2grn/model.py"],
    "rna2metabolite": ["rna2metabolite/api.py", "rna2metabolite/_impute.py"],
    "rna2lipid_aml": ["rna2lipid/aml/api.py", "rna2lipid/aml/_impute.py"],
    "fastComm": ["fastComm/api.py", "fastComm/scoring.py", "fastComm/upstream_resources.py"],
}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _runtime_versions():
    versions = {}
    for package in ("altanalyze3", "numpy", "scipy", "scikit-learn", "pandas", "anndata"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _json_default(value):
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize provenance value {type(value).__name__}")


def describe_model(component, artifacts, *, code_paths=None, catalog_path=REGISTRY_PATH):
    """Hash actual bytes, never infer an identity from a name or unpickle a proposal.

    Artifact role names form part of identity; absolute paths do not. Missing optional
    resources must be omitted explicitly by the caller. Required files fail closed.
    """
    resources = {role: {"sha256": sha256_file(path), "filename": Path(path).name}
                 for role, path in sorted(artifacts.items())}
    artifact_hashes = {role: record["sha256"] for role, record in resources.items()}
    if not resources:
        raise ValueError("A model must identify at least one artifact")
    paths = code_paths if code_paths is not None else [
        PACKAGE_ROOT / "components" / rel for rel in INFERENCE_FILES[component]]
    if code_paths is None and component == "rna2lipid":
        release_code = PACKAGE_ROOT / "components/rna2lipid/release.py"
        if release_code.is_file():
            paths.append(release_code)
    code = {(Path(path).relative_to(PACKAGE_ROOT).as_posix().replace("/", ":") if Path(path).is_relative_to(PACKAGE_ROOT) else Path(path).name): sha256_file(path) for path in paths}
    artifact_id = f"{component}:sha256:{_identity(artifact_hashes)}"
    inference_id = _identity(code)
    result = {
        "component": component,
        "model_version_id": artifact_id,
        "inference_version_id": f"sha256:{inference_id}",
        "analysis_model_version_id": f"{component}:sha256:{_identity({'artifacts': artifact_hashes, 'code': code})}",
        "artifacts": resources,
        "inference_code_sha256": code,
        "registry_status": "unregistered",
        "runtime_versions": _runtime_versions(),
    }
    if Path(catalog_path).exists():
        catalog = json.loads(Path(catalog_path).read_text())
        for record in catalog.get("models", []):
            if record["model_version_id"] == artifact_id:
                result["registry_status"] = record["status"]
                result["registry_name"] = record["name"]
                break
    return result


def write_provenance(output_path, provenance):
    """Write a mandatory companion to tabular outputs without changing table shape."""
    path = Path(str(output_path) + ".model_provenance.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"schema_version": "1.0", "models": provenance},
                               indent=2, sort_keys=True, default=_json_default) + "\n")
    return path


def write_run_provenance(output_dir, models, *, application):
    path = Path(output_dir) / "model_provenance.json"
    path.write_text(json.dumps({"schema_version": "1.0", "application": application,
        "recorded_at": datetime.now(timezone.utc).isoformat(), "models": models},
        indent=2, sort_keys=True, default=_json_default) + "\n")
    return path


def read_result_provenance(source):
    """Read stored IDs from a result; never assign today's model to old results."""
    source = Path(source)
    with source.open("rb") as handle:
        is_hdf5 = handle.read(8) == b"\x89HDF\r\n\x1a\n"
    if is_hdf5:
        import h5py
        try:
            from anndata.io import read_elem
        except ImportError:
            from anndata.experimental import read_elem
        with h5py.File(source, "r") as handle:
            for key in ("prediction_summary", "model_provenance", "imputation"):
                if "uns" in handle and key in handle["uns"]:
                    record = read_elem(handle["uns"][key])
                    if record.get("model_version_id"):
                        return json.loads(json.dumps(record, default=_json_default))
        return {}
    companion = Path(str(source) + ".model_provenance.json")
    if companion.exists():
        record = json.loads(companion.read_text()).get("models", {})
        return record if record.get("model_version_id") else {}
    return {}
