"""Layer 1 (ingest) — parse the Yang et al. 2016 (Cell, NIHMS754821) isoform PPI dataset.

This is a SECOND, independent, MEASURED isoform PPI resource ("Widespread Expansion of Protein
Interaction Capabilities by Alternative Splicing"). Unlike the TFIso atlas (paper2, TF-only), paper1
includes non-TF genes (e.g. ACTN4, BAG1), so it broadens coverage beyond transcription factors. It also
ships the mechanism behind isoform-specific PPIs: linear motifs (3B) and domain-domain interactions
(3D), which are the domain-contact basis of PPI gain/loss.

Reads the .xlsx workbooks via pandas/openpyxl and writes tidy TSVs (``paper1_*.tsv``) under
``config.DATA_DIR`` with a QC manifest. Isoforms are keyed by ``Isoform_ID`` (``SYMBOL_N``, underscore);
the crosswalk layer maps these to the canonical structure key via the ORF sequences in ``paper1_orfs``.

PPI calls are binary (``Interaction_Found`` = positive/negative); the ``detected`` column is True for
positive, False for negative (assayed-and-not-found), and NA for blank/untested — preserving the
tested-vs-untested distinction (NA rows are dropped from the edge table, kept nowhere as negatives).
"""

import os
import logging

import pandas as pd

from .. import config

logger = logging.getLogger(__name__)


def _read_sheet(modality):
    fname, sheet, _ = config.PAPER1_SHEETS[modality]
    path = config.require(os.path.join(config.PAPER1_DIR, fname), what=f"paper1 workbook {fname}")
    return pd.read_excel(path, sheet_name=sheet, dtype=str)


def _detected_from_found(series):
    def _one(v):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return pd.NA
        s = str(v).strip().lower()
        if s in config.TRUE_TOKENS or s == "positive":
            return True
        if s in config.FALSE_TOKENS or s == "negative":
            return False
        return pd.NA
    return series.map(_one).astype("boolean")


def load_paper1_isoforms():
    """Isoforms tested in Y2H (2A): Isoform_ID, gene/Entrez, GenBank accession, reference/alternative,
    known/novel."""
    return _read_sheet("paper1_isoforms")


def load_paper1_ppi():
    """Isoform PPIs (2B), one row per assayed isoform x interactor pair. Adds a nullable-boolean
    ``detected`` from ``Interaction_Found`` and a ``tested`` flag (True unless the call was blank)."""
    df = _read_sheet("paper1_ppi")
    df["detected"] = _detected_from_found(df["Interaction_Found"])
    df["tested"] = df["detected"].notna()
    return df


def load_paper1_orfs():
    """Reference + alternative ORFs (1B): Isoform_ID, category, NMD_target, novelty, NT_Variants, and
    the full ORF NUCLEOTIDE sequence (``Isoform__Open_Reading_Frame_Sequence``) used for the crosswalk."""
    return _read_sheet("paper1_orfs")


def load_paper1_linear_motifs():
    """Linear motifs in isoforms (3B): the motif (sequence + AA_Start/Stop) that drives an isoform's
    interaction phenotype, with the interacting vs non-interacting isoforms."""
    return _read_sheet("paper1_linear_motifs")


def load_paper1_domain_domain():
    """Domain-domain interactions (3D): the isoform domain (``Domain_In_Isoform``) and partner domain
    (``Domain_In_Partner``) mediating a PPI, with the interacting vs non-interacting isoform and whether
    only the isoform with the intact domain interacts — the domain-contact basis of isoform PPI loss."""
    return _read_sheet("paper1_domain_domain")


_LOADERS = {
    "paper1_isoforms": load_paper1_isoforms,
    "paper1_ppi": load_paper1_ppi,
    "paper1_orfs": load_paper1_orfs,
    "paper1_linear_motifs": load_paper1_linear_motifs,
    "paper1_domain_domain": load_paper1_domain_domain,
}


def parse_all(out_dir=None):
    """Parse every paper1 sheet, QC against the verified expected row count (when known), write tidy
    TSVs + ``paper1_ingest_manifest.tsv``. Returns the manifest DataFrame."""
    out_dir = out_dir or config.DATA_DIR
    os.makedirs(out_dir, exist_ok=True)
    records = []
    for modality, loader in _LOADERS.items():
        _, sheet, expected = config.PAPER1_SHEETS[modality]
        df = loader()
        out_path = os.path.join(out_dir, f"{modality}.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        observed = len(df)
        if expected is None:
            ok, note = (observed > 0), "" if observed > 0 else "EMPTY"
        else:
            ok = (observed == expected)
            note = "" if ok else f"ROW COUNT MISMATCH: expected {expected}, got {observed}"
        (logger.info if ok else logger.error)(
            "[iso2function.ingest.paper1] %-22s %5d rows -> %s  %s", modality, observed, out_path, note)
        records.append({"modality": modality, "sheet": sheet, "expected_rows": expected,
                        "observed_rows": observed, "qc_ok": ok, "note": note, "output": out_path})
    manifest = pd.DataFrame.from_records(records)
    manifest.to_csv(os.path.join(out_dir, "paper1_ingest_manifest.tsv"), sep="\t", index=False)
    # PPI call summary
    ppi = pd.read_csv(os.path.join(out_dir, "paper1_ppi.tsv"), sep="\t", dtype=str)
    pos = (ppi["detected"].astype(str).str.lower() == "true").sum()
    neg = (ppi["detected"].astype(str).str.lower() == "false").sum()
    logger.info("[iso2function.ingest.paper1] PPI calls: %d positive, %d negative (binary, measured)",
                pos, neg)
    return manifest


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    print(parse_all().to_string(index=False))
