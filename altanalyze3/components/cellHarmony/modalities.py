"""Modality identities shared by upload analysis and the precomputed viewer.

GRN edges (`grn`) and predicted factor activity (`grn_tf`) are separate feature spaces.
Selecting GRN imputation produces both; choosing a differential never swaps their inputs.
"""
from typing import Dict

MODALITY_DEFINITIONS: Dict[str, Dict[str, object]] = {
    "rna": {
        "id": "rna",
        "label": "RNA",
        "feature_label": "gene",
        "supports_marker_heatmap": True,
        "supports_marker_network": True,
        "supports_differential_network": True,
        "supports_differential_go": True,
    },
    "lipids": {
        "id": "lipids",
        "label": "Lipids",
        "feature_label": "lipid",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "adt": {
        "id": "adt",
        "label": "ADT (CITE-seq)",
        "feature_label": "ADT",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "metabolite": {
        "id": "metabolite",
        "label": "Metabolite (AML)",
        "feature_label": "metabolite",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "lipid": {
        "id": "lipid",
        "label": "Lipid (AML)",
        "feature_label": "lipid",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "grn": {
        "id": "grn",
        "label": "GRN (edges)",
        "feature_label": "edge",
        "supports_explore": False,
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "grn_tf": {
        "id": "grn_tf", "label": "TF activity (imputed)", "feature_label": "factor",
        "supports_marker_heatmap": True, "supports_marker_network": False,
        "supports_differential_network": False, "supports_differential_go": False,
    },
    "cell_communication": {
        "id": "cell_communication",
        "label": "Cell communication",
        "feature_label": "ligand-receptor interaction",
        "supports_marker_heatmap": False,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
}

def normalize_modality_id(value: object, *, default: str = "rna") -> str:
    raw = str(value or "").strip().lower()
    if not raw or raw in {"none", "null", "false", "off"}:
        return default
    if raw in {"lipids"}:
        return "lipids"
    if raw in {"lipid", "lipid_aml", "aml_lipid", "rna2lipid_aml"}:
        return "lipid"
    if raw in {"metabolite", "metabolites", "rna2metabolite"}:
        return "metabolite"
    if raw in {"grn", "grns", "regulon", "regulons", "gene_regulatory_network", "rna2grn", "grn_edges", "tf_edges"}:
        return "grn"
    if raw in {"grn_tf", "tf_activity", "regulator_activity"}:
        return "grn_tf"
    if raw in {"adt", "adts", "cite", "cite-seq", "citeseq"}:
        return "adt"
    if raw in {"cell_communication", "cell communication", "communication", "fastcomm", "fastcomm_network"}:
        return "cell_communication"
    if raw == "rna":
        return "rna"
    return raw


def modality_artifacts(meta):
    """Read old uploaded jobs without confusing their factor and edge matrices.

    Legacy uploads stored factor enrichment at grn.h5ad and edges at
    grn.differential_h5ad. Their enrichment has no factor pseudobulk differential.
    Precomputed bundles use feature TSVs and are not legacy uploaded jobs.
    """
    entries = {key: dict(value) for key, value in (meta.get("modality_artifacts") or {}).items()}
    grn = entries.get("grn", {})
    if ("grn_tf" not in entries and str(grn.get("h5ad", "")).endswith(".h5ad")
            and grn.get("differential_h5ad")
            and grn.get("h5ad") != grn.get("differential_h5ad")):
        entries["grn_tf"] = {"h5ad": grn["h5ad"], "legacy_enrichment": True}
        entries["grn"]["h5ad"] = grn["differential_h5ad"]
    return entries
