"""Readable modality/version directories backed by immutable artifact identities."""
import json
import re
from pathlib import Path

LABEL_FIELDS = ("model_modality", "model_version", "model_variant", "model_display_id", "registry_path")


def apply_labels(catalog, labels):
    used = set()
    for record in catalog["models"]:
        identity = record["model_version_id"]
        label = labels[identity]
        modality, version, variant = (label[key] for key in ("modality", "version", "variant"))
        if not re.fullmatch(r"v[1-9]\d*\.(0|[1-9]\d*)(?:\.(0|[1-9]\d*))?", version):
            raise ValueError(f"Invalid model version: {version}")
        for value in (modality, variant):
            if not re.fullmatch(r"[A-Za-z0-9_-]+", value):
                raise ValueError("Modality and variant must be safe directory names")
        key = (modality, version, variant)
        if key in used:
            raise ValueError(f"Version is assigned to multiple artifacts: {key}")
        used.add(key)
        proposed = {"model_modality": modality, "model_version": version, "model_variant": variant}
        for field, value in proposed.items():
            if record.get(field) is not None and record[field] != value:
                raise ValueError(f"Cannot relabel an existing artifact: {identity}")
        record.update(model_modality=modality, model_version=version, model_variant=variant,
                      model_display_id=f"{modality}/{variant}/{version}",
                      registry_path=f"models/{modality}/{version}/{variant}/model.json")
    by_id = {record["model_version_id"]: record for record in catalog["models"]}
    for default in catalog["defaults"]:
        for key in LABEL_FIELDS:
            default[key] = by_id[default["model_version_id"]][key]
    return catalog


def write_layout(catalog, repository):
    repository = Path(repository)
    grouped = {}
    current = []
    for record in catalog["models"]:
        path = repository / record["registry_path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        # Historical metadata/IDs stay in catalog snapshots; this card includes current annotations.
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        urls = []
        for role, artifact in record["artifacts"].items():
            if artifact.get("url"):
                urls.append(f"- [{role}: {artifact['filename']}]({artifact['url']})")
        text = f"# {record['model_display_id']}\n\nStatus: {record['status']}.\n\n"
        if record.get("declared_method"):
            text += "Recorded method, training scope and authorization are in [model.json](model.json).\n\n"
        if urls:
            text += "Download artifacts:\n\n" + "\n".join(urls) + "\n\n"
        text += "Immutable hashes and inference-code versions: [model.json](model.json).\n"
        (path.parent / "README.md").write_text(text)
        grouped.setdefault(record["model_modality"], []).append(record)
        if record["status"] == "observed-default":
            current.append({key: record[key] for key in LABEL_FIELDS + ("model_version_id",)})
    for modality, models in grouped.items():
        base = repository / "models" / modality
        lines = [f"# {modality}", "", "| Version | Variant | Status | Model |",
                 "| --- | --- | --- | --- |"]
        for record in sorted(models, key=lambda r: (tuple(map(int, r["model_version"][1:].split('.'))), r["model_variant"])):
            relative = f"{record['model_version']}/{record['model_variant']}/README.md"
            lines.append(f"| {record['model_version']} | {record['model_variant']} | {record['status']} | [{record['model_display_id']}]({relative}) |")
        lines += ["", "Versions apply independently to each tissue/species variant within this modality.",
                  "v1.0 identifies the first registered model in that series; later revisions increment v1.1, v1.2, etc.",
                  "These labels do not assign training dates or imply tissue variants are interchangeable.", ""]
        (base / "README.md").write_text("\n".join(lines))
        for version in sorted({record["model_version"] for record in models}):
            version_models = [record for record in models if record["model_version"] == version]
            version_dir = base / version
            links = [f"- [{record['model_variant']}]({record['model_variant']}/README.md): {record['status']}"
                     for record in sorted(version_models, key=lambda r: r["model_variant"])]
            (version_dir / "README.md").write_text(f"# {modality} {version}\n\n" + "\n".join(links)
                    + f"\n\n[Download model artifacts](https://github.com/SalomonisLab/scalable-models/releases/tag/{modality}-{version}).\n")
    (repository / "models" / "current.json").write_text(json.dumps({"models": current}, indent=2, sort_keys=True) + "\n")
