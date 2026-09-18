"""Presentation-only identities for the audited CPTAC AML imputation features."""
from functools import lru_cache
import json
from pathlib import Path


@lru_cache(maxsize=1)
def catalog():
    return json.loads((Path(__file__).with_name("resources") / "aml_feature_annotations.json").read_text())


def for_dataset(meta, modality="metabolite"):
    """Never apply study-local Unknown IDs to an unrelated metabolomics assay."""
    if modality != "metabolite":
        return {}
    entries = (meta.get("modalities") or {}).get("available", [])
    entry = next((e for e in entries if e.get("id") == "metabolite"), {})
    if entry.get("annotation_source") == "PDC000561" or "AML" in str(entry.get("label", "")).upper():
        return catalog()
    return {}


def display_name(feature, annotations):
    name = str(feature)
    if name in annotations:
        return annotations[name]["label"]
    # Marker heatmap row IDs include a cluster prefix; keep it intact.
    prefix, separator, suffix = name.rpartition(":")
    if separator and suffix in annotations:
        return prefix + separator + annotations[suffix]["label"]
    return name


def pdf_payload(payload, meta, modality=None):
    """Copy display fields only, after all feature lookups have completed."""
    modality = modality or (meta.get("differential", {}).get("config", {}) or {}).get("modality", "rna")
    annotations = for_dataset(meta, modality)
    result = dict(payload)
    if "gene" in result:
        result["gene"] = display_name(result["gene"], annotations)
    if "rows" in result:
        result["rows"] = [dict(row, gene=display_name(row.get("gene", ""), annotations)) for row in result["rows"]]
    return result


def label_tsv(content, annotations):
    """Label a visualization export; leave all measurement columns untouched."""
    if not annotations:
        return content
    lines = content.splitlines(keepends=True)
    for i in range(1, len(lines)):
        key, separator, rest = lines[i].partition("\t")
        if separator:
            lines[i] = display_name(key, annotations) + separator + rest
    return "".join(lines)
