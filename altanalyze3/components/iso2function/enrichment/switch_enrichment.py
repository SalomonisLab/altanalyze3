"""Isoform-switch enrichment (P4) — structure-keyed, binary-only, high-confidence.

Two products, both keyed on the canonical structure string:

  - ``switch_differential_interactions.tsv`` : for every reference->alternative switch (from
    ``switch_function.tsv``), the PPI partners and PDI DNA targets GAINED and LOST across the switch,
    computed from the authors' binary calls within the CO-TESTED space of both isoforms
    (``network.build.differential_interactome``). This is the mechanistic detail behind the authors'
    ``Data_S8`` PDI/PPI categories — and a built-in cross-check against them.
  - ``isoform_gsea_by_structure.tsv`` : the authors' own alternative-vs-reference isoform GSEA
    (``Data_S10``, MSigDB terms + NES + q-value) re-keyed from clone_ids to structure strings.

Also exposes :func:`enrich_partner_categories` — a hypergeometric test (no external gene sets needed)
asking whether a switch's gained/lost PPI partners are enriched for an interactor class (TF / cofactor
type), the same axis the authors' ``PPI loss: cofactor`` categories use.
"""

import os
import logging

import pandas as pd

from .. import config
from ..network import build
from . import isoform_gsea

logger = logging.getLogger(__name__)


def differential_interactions_by_switch(data_dir=None):
    """Per ref->alt switch: gained/lost PPI partners and PDI targets (binary, co-tested). Returns the
    DataFrame and writes ``switch_differential_interactions.tsv``. Requires ``switch_function.tsv``
    (run associate) and the tidy edge tables (run ingest)."""
    data_dir = data_dir or config.DATA_DIR
    sw_path = os.path.join(data_dir, "switch_function.tsv")
    if not os.path.exists(sw_path):
        raise FileNotFoundError("switch_function.tsv not found — run `iso2function associate` first.")
    sw = pd.read_csv(sw_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    _, edges = build.build_graph(data_dir=data_dir)          # all called edges (Y2H + eY1H)

    rows = []
    for _, r in sw.iterrows():
        ref, alt = r.get("ref_clone"), r.get("alt_clone")
        ppi = build.differential_interactome(alt, ref, edges, edge_type="ppi")   # up=alt, down=ref
        pdi = build.differential_interactome(alt, ref, edges, edge_type="pdi")
        rows.append({
            "gene_symbol": r.get("gene_symbol"),
            "ref_clone": ref, "alt_clone": alt,
            "ref_structure": r.get("ref_structure", ""), "alt_structure": r.get("alt_structure", ""),
            "ref_enst": r.get("ref_enst", ""), "alt_enst": r.get("alt_enst", ""),
            "PDI_category": r.get("PDI_category", ""), "PPI_category": r.get("PPI_category", ""),
            "n_ppi_cotested": len(ppi["co_tested_targets"]),
            "n_ppi_gained": len(ppi["gained"]), "n_ppi_lost": len(ppi["lost"]),
            "ppi_gained": ";".join(ppi["gained"]), "ppi_lost": ";".join(ppi["lost"]),
            "n_pdi_cotested": len(pdi["co_tested_targets"]),
            "n_pdi_gained": len(pdi["gained"]), "n_pdi_lost": len(pdi["lost"]),
            "pdi_gained": ";".join(pdi["gained"]), "pdi_lost": ";".join(pdi["lost"]),
        })
    out = pd.DataFrame(rows)
    path = os.path.join(data_dir, "switch_differential_interactions.tsv")
    out.to_csv(path, sep="\t", index=False)
    logger.info("[iso2function.enrich] switch_differential_interactions.tsv: %d switches "
                "(%d with PPI co-tested data, %d with PDI co-tested data)", len(out),
                int((out["n_ppi_cotested"] > 0).sum()), int((out["n_pdi_cotested"] > 0).sum()))
    return out


def gsea_by_structure(data_dir=None):
    """Re-key the authors' alternative-vs-reference isoform GSEA (Data_S10) from clone_ids to structure
    strings via ``clone_to_structure.tsv``. Returns the DataFrame and writes
    ``isoform_gsea_by_structure.tsv``."""
    data_dir = data_dir or config.DATA_DIR
    g_path = os.path.join(data_dir, "isoform_gsea.tsv")
    c_path = os.path.join(data_dir, "clone_to_structure.tsv")
    if not (os.path.exists(g_path) and os.path.exists(c_path)):
        raise FileNotFoundError("need isoform_gsea.tsv (ingest) + clone_to_structure.tsv (crosswalk).")
    g = pd.read_csv(g_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    cts = pd.read_csv(c_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    s = cts.drop_duplicates("clone_id").set_index("clone_id")["structure"].to_dict()
    g["ref_structure"] = g["ref_iso"].map(s).fillna("")
    g["alt_structure"] = g["alt_iso"].map(s).fillna("")
    g["both_mapped"] = (g["ref_structure"] != "") & (g["alt_structure"] != "")
    for col in ("gsea_nes", "gsea_qval"):
        if col in g.columns:
            g[col] = pd.to_numeric(g[col], errors="coerce")
    path = os.path.join(data_dir, "isoform_gsea_by_structure.tsv")
    g.to_csv(path, sep="\t", index=False)
    logger.info("[iso2function.enrich] isoform_gsea_by_structure.tsv: %d enrichment rows "
                "(%d with both isoforms mapped to a structure)", len(g), int(g["both_mapped"].sum()))
    return g


def _interactor_category_sets(data_dir):
    """Build interactor-class gene sets from the Y2H table itself (no external gene sets): for each
    value of ``db_gene_category`` and ``db_gene_cofactor_type``, the set of partner gene symbols."""
    y2h_path = os.path.join(data_dir, "ppi_y2h.tsv")
    if not os.path.exists(y2h_path):
        return {}, set()
    y = pd.read_csv(y2h_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    universe = set(y["db_gene_symbol"].dropna().astype(str))
    sets = {}
    for col, prefix in (("db_gene_category", "category"), ("db_gene_cofactor_type", "cofactor")):
        if col in y.columns:
            for val, sub in y.groupby(col):
                if str(val).strip():
                    sets[f"{prefix}:{val}"] = set(sub["db_gene_symbol"].dropna().astype(str))
    return sets, universe


def enrich_partner_categories(gained_or_lost, data_dir=None, min_overlap=2):
    """Hypergeometric enrichment of a partner-gene list against interactor-class gene sets (TF /
    cofactor type), within the assayed Y2H partner universe. Demonstrates the isoform-focused
    enrichment layer end-to-end without external gene sets. Returns the enrichment DataFrame."""
    data_dir = data_dir or config.DATA_DIR
    sets, universe = _interactor_category_sets(data_dir)
    if not sets:
        logger.warning("[iso2function.enrich] no Y2H table for category enrichment")
        return pd.DataFrame()
    return isoform_gsea.enrich(gained_or_lost, sets, universe, min_overlap=min_overlap)


def build_all(data_dir=None):
    """Build the structure-keyed enrichment products. Returns (differential, gsea)."""
    return differential_interactions_by_switch(data_dir), gsea_by_structure(data_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    diff, gsea = build_all()
    print(diff[["gene_symbol", "ref_clone", "alt_clone", "PPI_category",
                "n_ppi_gained", "n_ppi_lost"]].head(10).to_string(index=False))
