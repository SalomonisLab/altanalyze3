"""Scale and correction-family declarations for the supported lung lipids."""
import numpy as np


def is_native_lipid_log2(adata):
    return adata.uns.get("expression_scale") == "native_relative_log2"


def full_lipid_bh_mask(adata):
    """Return the prespecified full panel, or None for unchanged other modalities."""
    if adata.uns.get("modality") != "lipids" and not is_native_lipid_log2(adata):
        return None
    from altanalyze3.components.rna2lipid.release import release_manifest, exact_identifiers
    expected = release_manifest()["Y_columns"]
    exact_identifiers(expected, adata.var_names, "Lipid differential testing panel")
    return np.ones(adata.n_vars, dtype=bool)


def full_imputed_bh_mask(adata):
    """User-authorized complete-panel BH for all imputed modalities."""
    from .imputed_scale import uses_imputed_scale, prediction_encoding
    lipid_mask = full_lipid_bh_mask(adata)
    if lipid_mask is not None:
        return lipid_mask
    if not uses_imputed_scale(adata):
        return None
    prediction_encoding(adata)
    if not adata.var_names.is_unique:
        raise ValueError("Duplicate imputed feature identities; resolve the panel before testing.")
    return np.ones(adata.n_vars, dtype=bool)
