"""Layer 1 (ingest) — parse the Mol Cell 2025 TFIso atlas supplementary tables into tidy, validated,
long-form interaction tables.

Reads each modality directly out of its supplementary ``.zip`` (members are streamed, never extracted
into the read-only download folder) and writes tidy TSVs under ``config.DATA_DIR``. Every parse is QC'd
against the verified expected row count (``config.PAPER2_TABLES``) and the result is logged.

Tidy outputs (all keyed on ``clone_id`` = ``SYMBOL-N``, the paper2 namespace; the crosswalk layer maps
these to the canonical structure string):
  - ``ppi_y2h.tsv``        long PPI edges (one row per assayed AD-clone x DB-partner pair)
  - ``ppi_n2h.tsv``        long PPI edges (scored; log2 NLR)
  - ``pdi_ey1h.tsv``       long PDI edges (melted from the clone x DNA-bait matrix)
  - ``pdi_dna_baits.tsv``  bait_id -> DNA sequence
  - ``pdi_validation.tsv`` quantitative PDI validation (replicates, log2 FC)
  - ``activation_m1h.tsv`` per-clone M1H activation (reps + mean)
  - ``condensate.tsv``     per-clone condensate / localization phenotypes
  - ``paralog_divergence.tsv`` pairwise paralog PDI/PPI Jaccard distances + sequence identity
  - ``ingest_manifest.tsv``    per-table QC: expected vs observed row counts, output path

The ``tested`` column on PPI/PDI edges preserves the crucial tested-vs-untested distinction: a Boolean
``False`` is *assayed-and-not-detected*, not "interaction absent in biology". Pairs absent from an
assay matrix are simply absent from the tidy table (never emitted as negatives).
"""

import os
import io
import zipfile
import logging

import pandas as pd

from .. import config

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- low-level readers
def _read_zip_member(zip_name, member, **read_csv_kwargs):
    """Return a DataFrame read from ``member`` inside ``<PAPER2_DIR>/<zip_name>`` (tab-separated by
    default). The zip is opened read-only; nothing is extracted to disk."""
    zip_path = config.require(os.path.join(config.PAPER2_DIR, zip_name), what=f"paper2 archive {zip_name}")
    with zipfile.ZipFile(zip_path) as zf:
        if member not in zf.namelist():
            raise KeyError(f"{member!r} not found in {zip_path} (members: {zf.namelist()})")
        raw = zf.read(member)
    kwargs = dict(sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    kwargs.update(read_csv_kwargs)
    return pd.read_csv(io.BytesIO(raw), **kwargs)


def _parse_bool_series(series):
    """Map a string series of result tokens to nullable boolean: TRUE/FALSE tokens -> True/False,
    anything else (incl. empty/NaN) -> pd.NA (untested)."""
    def _one(v):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return pd.NA
        s = str(v).strip().lower()
        if s in config.TRUE_TOKENS:
            return True
        if s in config.FALSE_TOKENS:
            return False
        return pd.NA
    return series.map(_one).astype("boolean")


# --------------------------------------------------------------------------- per-modality loaders
def load_ppi_y2h():
    """Pairwise Y2H PPIs. AD-clone (the test TF isoform) x DB-partner. One row per assayed pair."""
    zip_name, member, _ = config.PAPER2_TABLES["ppi_y2h"]
    df = _read_zip_member(zip_name, member)
    df = df.rename(columns={"Y2H_result": "y2h_result_raw"})
    df["detected"] = _parse_bool_series(df["y2h_result_raw"])
    df["tested"] = True  # present in the pairwise result table == assayed
    return df


def load_ppi_n2h():
    """N2H PPIs (orthogonal validation, SCORE-ONLY). ``log2 NLR`` is the interaction strength, but the
    paper states NO numeric NLR cutoff in the metadata (only a dashed line in Fig S2I). Under the
    binary-only policy this table is parsed for reference (and to keep the scores available as
    attributes) but is NEVER used as a called interaction edge. There is no per-row Boolean call."""
    zip_name, member, _ = config.PAPER2_TABLES["ppi_n2h"]
    df = _read_zip_member(zip_name, member)
    df = df.rename(columns={"log2 NLR": "log2_nlr"})
    for col in ("score_pair", "score_empty-N1", "score_empty-N2", "log2_nlr"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["tested"] = True
    df["called"] = False  # explicit: no binary call available; excluded from the network
    return df


def load_pdi_ey1h():
    """eY1H PDI matrix (clone x DNA-bait Boolean), melted to long edges. The first two columns are
    ``gene_symbol`` and ``clone_id``; every remaining column is a DNA bait id."""
    zip_name, member, _ = config.PAPER2_TABLES["pdi_ey1h"]
    wide = _read_zip_member(zip_name, member)
    id_cols = ["gene_symbol", "clone_id"]
    missing = [c for c in id_cols if c not in wide.columns]
    if missing:
        raise ValueError(f"eY1H matrix missing expected id columns {missing}; got {list(wide.columns)[:5]}")
    bait_cols = [c for c in wide.columns if c not in id_cols]
    long = wide.melt(id_vars=id_cols, value_vars=bait_cols, var_name="bait_id", value_name="y1h_raw")
    long["detected"] = _parse_bool_series(long["y1h_raw"])
    long["tested"] = True  # every cell in the matrix was assayed
    # keep only assayed cells with a definite call; drop the (rare) blank cells as untested
    long = long[long["detected"].notna()].reset_index(drop=True)
    return long


def load_pdi_dna_baits():
    """bait_id -> DNA bait sequence."""
    zip_name, member, _ = config.PAPER2_TABLES["pdi_dna_baits"]
    return _read_zip_member(zip_name, member)


def load_pdi_validation():
    """Quantitative PDI validation (Y1H replicate activation + log2 FC)."""
    zip_name, member, _ = config.PAPER2_TABLES["pdi_validation"]
    df = _read_zip_member(zip_name, member)
    df = df.rename(columns={"Y1H_result": "y1h_result_raw", "Log2(FC)": "log2_fc"})
    df["detected"] = _parse_bool_series(df["y1h_result_raw"])
    for col in ("Replicate1", "Replicate2", "Replicate3", "Average (empty-pEZY3-VP160)", "log2_fc"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _m1h_call(mean):
    """Apply the paper's verbatim M1H rule to an M1H_mean signal: >= M1H_ACTIVATOR_MIN -> 'activator',
    <= M1H_REPRESSOR_MAX -> 'repressor', else 'neither'. NaN -> '' (untested)."""
    if mean is None or (isinstance(mean, float) and pd.isna(mean)):
        return ""
    if mean >= config.M1H_ACTIVATOR_MIN:
        return "activator"
    if mean <= config.M1H_REPRESSOR_MAX:
        return "repressor"
    return "neither"


def load_activation_m1h():
    """Per-clone M1H transcriptional activation (Data_S3): reps + mean + a paper-defined ``m1h_call``
    (activator/repressor/neither) from the verbatim cutoff M1H signal >= 1 / <= -1 (mmc12.pdf p.e12)."""
    zip_name, member, _ = config.PAPER2_TABLES["activation_m1h"]
    df = _read_zip_member(zip_name, member)
    for col in ("M1H_rep1", "M1H_rep2", "M1H_rep3", "M1H_mean"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["m1h_call"] = df["M1H_mean"].map(_m1h_call) if "M1H_mean" in df.columns else ""
    return df


def load_functional_categories():
    """Authors' OWN per-switch functional consequence categories (Data_S8): one row per
    reference/alternative isoform pair (clone-keyed, SYMBOL-N) with PDI_category, PPI_category,
    M1H_activation_category, localization_category, DBD_pct_lost_in_alt, and the rewirer/negative-
    regulator classification. These are the high-confidence calls used for switch interpretation."""
    zip_name, member, _ = config.PAPER2_TABLES["functional_categories"]
    df = _read_zip_member(zip_name, member)
    if "DBD_pct_lost_in_alt" in df.columns:
        df["DBD_pct_lost_in_alt"] = pd.to_numeric(df["DBD_pct_lost_in_alt"], errors="coerce")
    return df


def load_isoform_gsea():
    """Alternative-vs-reference isoform GSEA (Data_S10): per (alt_iso, ref_iso) clone pair, the MSigDB
    term, normalized enrichment score, and q-value (from Joung et al. TF-ORF overexpression)."""
    zip_name, member, _ = config.PAPER2_TABLES["isoform_gsea"]
    df = _read_zip_member(zip_name, member)
    for col in ("gsea_nes", "gsea_qval"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_paralog_divergence():
    """Pairwise paralog divergence (Data_S6): sequence identity + PDI/PPI Jaccard distances."""
    zip_name, member, _ = config.PAPER2_TABLES["paralog_divergence"]
    df = _read_zip_member(zip_name, member)
    for col in ("aa_seq_pct_identity", "PDI_Jaccard_d", "PPI_Jaccard_d", "activation_abs_fold_change_log2"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_condensate():
    """Per-clone condensate formation + subcellular localization (Data_S7), 2 cell lines x 2 reps."""
    zip_name, member, _ = config.PAPER2_TABLES["condensate"]
    return _read_zip_member(zip_name, member)


# loader registry: modality -> (loader_fn, "matrix" rows are melted so observed != expected wide rows)
_LOADERS = {
    "ppi_y2h": (load_ppi_y2h, "rows"),
    "ppi_n2h": (load_ppi_n2h, "rows"),
    "pdi_ey1h": (load_pdi_ey1h, "melted"),   # long form; QC against melt, not the wide row count
    "pdi_dna_baits": (load_pdi_dna_baits, "rows"),
    "pdi_validation": (load_pdi_validation, "rows"),
    "activation_m1h": (load_activation_m1h, "rows"),
    "paralog_divergence": (load_paralog_divergence, "rows"),
    "condensate": (load_condensate, "rows"),
    "functional_categories": (load_functional_categories, "rows"),
    "isoform_gsea": (load_isoform_gsea, "rows"),
}


# --------------------------------------------------------------------------- driver
def parse_all(out_dir=None):
    """Parse every paper2 modality, QC against expected row counts, write tidy TSVs + a manifest.

    Returns the manifest DataFrame. For ``rows``-type tables the observed data-row count must equal the
    verified expected count; a mismatch is logged as an ERROR (and flagged in the manifest) so a changed
    input never silently corrupts downstream layers. For the melted matrix (``pdi_ey1h``) the wide
    expected count is the number of isoform rows; the long output is larger, so we record both.
    """
    out_dir = out_dir or config.DATA_DIR
    os.makedirs(out_dir, exist_ok=True)
    records = []
    for modality, (loader, qc_kind) in _LOADERS.items():
        zip_name, member, expected = config.PAPER2_TABLES[modality]
        df = loader()
        out_path = os.path.join(out_dir, f"{modality}.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        observed = len(df)
        if qc_kind == "rows":
            ok = (observed == expected)
            note = "" if ok else f"ROW COUNT MISMATCH: expected {expected}, got {observed}"
        else:  # melted matrix: long form is expected to exceed the wide isoform-row count
            ok = (observed >= expected)
            note = f"melted long edges from {expected} isoform rows"
            if not ok:
                note = f"UNEXPECTED: melted rows {observed} < isoform rows {expected}"
        (logger.info if ok else logger.error)(
            "[iso2function.ingest] %-20s %6d rows -> %s  %s", modality, observed, out_path, note
        )
        records.append({
            "modality": modality, "source_zip": zip_name, "source_member": member,
            "expected_rows": expected, "observed_rows": observed, "qc_kind": qc_kind,
            "qc_ok": ok, "note": note, "output": out_path,
        })
    manifest = pd.DataFrame.from_records(records)
    manifest_path = os.path.join(out_dir, "ingest_manifest.tsv")
    manifest.to_csv(manifest_path, sep="\t", index=False)
    n_bad = int((~manifest["qc_ok"]).sum())
    logger.info("[iso2function.ingest] wrote manifest %s (%d/%d tables OK)",
                manifest_path, len(manifest) - n_bad, len(manifest))
    return manifest


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    m = parse_all()
    print(m.to_string(index=False))
