"""Canonical input/output paths and namespace constants for iso2function.

Every path a layer needs is resolved here so the rest of the package never hard-codes a location.
Shared-drive paths are written in their local ``/Volumes/...`` form and transparently mapped to the
cluster ``/data/...`` form (and vice versa) by :func:`resolve_shared` based on what actually exists.

Nothing here reads or writes; callers decide. Paths are validated lazily via :func:`require`.
"""

import os

# --------------------------------------------------------------------------- package layout
_HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT_DIR = _HERE
DATA_DIR = os.path.join(_HERE, "data")          # gitignored: extracted raw + parsed tidy tables
EXTERNAL_DIR = os.path.join(_HERE, "external")  # gitignored: downloaded reference snapshots
ARTIFACTS_DIR = os.path.join(_HERE, "artifacts")  # gitignored: scratch / figures / per-run outputs
RESOURCES_DIR = os.path.join(_HERE, "resources")  # tracked: small bundled lookups

# altanalyze3 repo root (…/altanalyze3/altanalyze3/components/iso2function -> repo root)
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))


# --------------------------------------------------------------------------- shared-drive resolver
def resolve_shared(path):
    """Return ``path`` if it exists; else try swapping the ``/Volumes/salomonis*`` (local mount) and
    ``/data/salomonis*`` (cluster) prefixes and return whichever exists. If neither exists, return the
    original ``path`` unchanged (callers surface a clear error via :func:`require`)."""
    if os.path.exists(path):
        return path
    swaps = [("/Volumes/", "/data/"), ("/data/", "/Volumes/")]
    for src, dst in swaps:
        if path.startswith(src):
            alt = dst + path[len(src):]
            if os.path.exists(alt):
                return alt
    return path


def require(path, what="input"):
    """Return ``path`` after asserting it exists (post shared-drive resolution); else raise with a
    clear, actionable message. Use at the point of consumption, never at import time."""
    resolved = resolve_shared(path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(
            f"iso2function: required {what} not found: {path}"
            + (f" (also tried {resolve_shared(path)})" if resolve_shared(path) != path else "")
        )
    return resolved


# --------------------------------------------------------------------------- v1 inputs: paper2 atlas
# Mol Cell 2025 TFIso atlas supplementary archives. Each modality lives as a named member inside a
# zip; ingest streams members out of the zip (read-only) into DATA_DIR rather than extracting in place.
TFISO_DOWNLOADS = "/Users/saljh8/Downloads/TFIso"
PAPER2_DIR = os.path.join(TFISO_DOWNLOADS, "paper2")
PAPER1_DIR = os.path.join(TFISO_DOWNLOADS, "paper1")

# paper1 = Yang et al. 2016 (Cell, NIHMS754821) "Widespread Expansion of Protein Interaction
# Capabilities by Alternative Splicing" — a SECOND, MEASURED isoform PPI dataset that INCLUDES non-TF
# genes (e.g. ACTN4, BAG1). xlsx workbooks (read via pandas/openpyxl). modality -> (file, sheet,
# expected_rows). Isoforms keyed by Isoform_ID = SYMBOL_N (underscore); full ORF nucleotide sequences
# in 1B enable a sequence-based crosswalk to the structure key. Headers verified 2026-06-18.
PAPER1_SHEETS = {
    "paper1_isoforms":      ("NIHMS754821-supplement-11.xlsx", "2A-Isoforms tested in Y2H", 1035),
    "paper1_ppi":           ("NIHMS754821-supplement-11.xlsx", "2B-Isoform PPIs", 2026),  # binary
    "paper1_orfs":          ("NIHMS754821-supplement-10.xlsx", "1B-Ref+Alt ORFs", None),  # ORF seqs
    "paper1_linear_motifs": ("NIHMS754821-supplement-12.xlsx", "3B-Linear Motifs in Isoforms", 151),
    "paper1_domain_domain": ("NIHMS754821-supplement-12.xlsx", "3D-Domain-Domain Interactions", 60),
}

# modality -> (zip filename, member filename, expected_data_rows). Headers + counts VERIFIED by
# streaming each member 2026-06-17 (data rows = nrows - 1 header). NOTE: Data_S3/S6/S7 contents were
# mislabeled in an early inventory; the assignments below are confirmed from the actual headers.
PAPER2_TABLES = {
    "pdi_ey1h":        ("mmc3.zip", "SuppTable_eY1HResults.txt", 171),   # 171 iso x 186 baits (Bool)
    "pdi_dna_baits":   ("mmc3.zip", "SuppTable_DNABaits.txt", 186),      # bait_id, seq
    "pdi_validation":  ("mmc3.zip", "SuppTable_PDI_validation.txt", 113),
    "ppi_y2h":         ("mmc5.zip", "SuppTable_PairwiseY2HResults.txt", 9562),
    "ppi_n2h":         ("mmc5.zip", "SuppTable_N2HResults.txt", 765),
    "activation_m1h":  ("mmc4.zip", "Data_S3.tsv", 622),                 # clone_id, M1H_rep1..3, mean
    "paralog_divergence": ("mmc7.zip", "Data_S6.tsv", 2139),            # *_a/*_b, Jaccard distances
    "condensate":      ("mmc8.zip", "Data_S7.tsv", 189),                # condensates/localization
    # authors' OWN high-confidence per-switch functional categories (ref vs alt isoform), clone-keyed:
    "functional_categories": ("mmc9.zip", "Data_S8.tsv", 447),          # PDI/PPI/M1H/localization category
    "isoform_gsea":    ("mmc11.zip", "Data_S10.tsv", 1156),             # alt vs ref isoform MSigDB GSEA
    # NOTE: mmc6 PBM tables (CREB1/TBX5 8-mer DNA specificity, 32896 rows each) are available but niche;
    # not ingested in v1. N2H (ppi_n2h) is score-only with NO numeric cutoff in the metadata, so it is
    # parsed for reference but NEVER used as a called interaction edge (binary-only policy).
}

# --------------------------------------------------------------------------- bridge / crosswalk inputs
# Prior TFIso annotation: GFF clone -> structure -> observed isoform (+ protein consequences).
TFISO_ASSIGNMENT = (
    "/Users/saljh8/Claude_Project/TFiso_AltAnalyze3_annotation/"
    "tfiso_to_final_assignment.protein.tsv"
)
# TFIso ORFeome clones aligned to the genome (source of GFF clone IDs; exon coords per clone).
TFISO_GFF = "/Users/saljh8/Dropbox/Revio/ENCODE/TFiso/c_6k_unique_acc_aligns.gff"
# GFF-annotated Ens91 structure for EVERY 2025-catalogue clone (the authoritative atlas structure
# source; transcript_id = the GFF "final name" SYMBOL|i/n|well, structure = E/I tokens, or a genomic
# coordinate for the 76 UNK-locus clones). The 2025 interaction IDs are harmonized to these names.
TFISO_GFF_STRUCTURES = ("/Users/saljh8/Dropbox/Revio/ENCODE/TFiso/TFiso_AltAnalyze3_annotation/"
                        "gff-output/transcript_associations.txt")
# Authoritative 2025 TFIso1.0 clone library (github CCSB-DFCI/TF_isoforms_paper/supp): clone_id ->
# ensembl_transcript_ids + cds_seq + aa_seq. MANE Select summary (NCBI): ENSG/symbol -> the gene's
# standard reference transcript, used as the fallback ENST for clones with no specific transcript.
SUPP_CLONELIST = os.path.join(COMPONENT_DIR, "external", "SuppTable_CloneList.txt")
MANE_SUMMARY = os.path.join(COMPONENT_DIR, "external", "MANE.summary.txt.gz")
# WikiPathways GPML release (Homo sapiens) for pathway-constrained network contrasts.
WIKIPATHWAYS_ZIP = "/Users/saljh8/Downloads/wikipathways-20260610-gpml-Homo_sapiens.zip"
INTERACTIONS_TXT = os.path.join(DATA_DIR, "isoform_interactions.txt")

# Higher-accuracy isoform->ENST->Ens91-structure mappings produced in the TFiso annotation project
# (Yang 2016 ORFs and UniProt proteins mapped to Ensembl116 transcripts, GFF-derived, re-annotated to
# Ensembl91 structures). See that project's README.md. Relational ID->ENST tables preserve the source
# (Vidal lab 2016 / UniProt) isoform identity. Preferred over the component's own crosswalk.
TFISO_ANNOT_DIR = os.path.dirname(TFISO_ASSIGNMENT)
YANG_ID_TO_ENST = os.path.join(TFISO_ANNOT_DIR, "yang_isoformID_to_ENST.txt")
YANG_STRUCTURES = os.path.join(TFISO_ANNOT_DIR, "yang_isoform_structures.txt")
UNIPROT_ID_TO_ENST = os.path.join(TFISO_ANNOT_DIR, "uniprot_to_ENST.txt")
UNIPROT_STRUCTURES_MAP = os.path.join(TFISO_ANNOT_DIR, "uniprot_isoform_structures.txt")
# Ensembl116 peptide FASTA the session matched against (header: 'gene:ENSG... transcript:ENST...');
# resolves ENST->ENSG for the Ensembl116-only transcripts absent from Ens91/GENCODE.
ENSEMBL_PEP_FASTA = os.path.join(TFISO_ANNOT_DIR, "ensembl_ref", "Homo_sapiens.GRCh38.pep.all.fa.gz")

# --------------------------------------------------------------------------- reference inputs
ENSEMBL91_EXON = (
    "/Volumes/salomonis-archive/LabFiles/Nathan/Revio/Hs/Ensembl91/Hs_Ensembl_exon.txt"
)
ENSEMBL91_ANNOT = (
    "/Volumes/salomonis-archive/LabFiles/Nathan/Revio/Hs/Ensembl91/Hs_Ensembl-annotations.txt"
)
GENOME_FASTA = "/Volumes/salomonis-archive/LabFiles/Nathan/Revio/Hs/genome.fa"

# MDS-AML-KINNEX-4 collapse outputs (observed-isoform catalog + protein summary).
MDS_AML_GFF_OUTPUT = (
    "/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-4/gff-output"
)
# Reference ENST -> Ensembl91 structure (244k transcripts) and structure -> genomic coordinates store
# (packed exon coords, decoded with bam.coord_store.unpack_exons) + all-isoform protein FASTA. Produced
# by the MDS-AML collapse; reused to map any matched ENST to its Ensembl91 structure + genomic splice
# sites (the portable, build-independent annotation).
ENST_REFERENCE_STRUCTURES = os.path.join(MDS_AML_GFF_OUTPUT, "ENST_reference_structures.tsv")
STRUCTURE_COORDS_SQLITE = os.path.join(MDS_AML_GFF_OUTPUT, "structure_coords.sqlite")
MDS_AML_PROTEIN_FASTA = os.path.join(MDS_AML_GFF_OUTPUT, "protein_sequences.fasta")        # AA, by final_id
MDS_AML_TRANSCRIPT_FASTA = os.path.join(MDS_AML_GFF_OUTPUT, "transcript_sequences.fasta")  # nt, by final_id
# gene <-> final_isoform_id <-> Ens91 structure for ALL cohort full-length isoforms (known + novel)
TRANSCRIPT_ASSOCIATIONS = os.path.join(MDS_AML_GFF_OUTPUT, "transcript_associations.txt")

# --------------------------------------------------------------------------- reused reference resources
# Processed protein/feature/interaction resources already assembled by the Daedalus component (no
# rebuild/download). Used to crosswalk external isoforms (paper1, UniProt) to ENST by protein-sequence
# match and to annotate domain/interface features for domain-contact PPI inference.
_DAEDALUS = os.path.join(REPO_ROOT, "altanalyze3", "components", "Daedalus", "data_ingest", "data")
GENCODE_PROTEIN_FASTA = os.path.join(_DAEDALUS, "raw", "gencode", "gencode.v49.pc_translations.fa.gz")
UNIPROT_ISOFORM_PROTEINS = os.path.join(_DAEDALUS, "interim", "uniprot_isoform_proteins.tsv")
UNIPROT_ISOFORM_FEATURES = os.path.join(_DAEDALUS, "interim", "uniprot_isoform_features.tsv")
UNIPROT_FEATURES = os.path.join(_DAEDALUS, "interim", "uniprot_features.tsv")
UNIPROT_ENTRIES = os.path.join(_DAEDALUS, "interim", "uniprot_entries.tsv")
BIOGRID_ZIP = os.path.join(_DAEDALUS, "raw", "BIOGRID-ALL-LATEST.tab3.zip")

# --------------------------------------------------------------------------- namespace constants
# Boolean-y string tokens seen in the atlas (Y2H/Y1H results). Anything not in TRUE/FALSE is "untested".
TRUE_TOKENS = {"true", "1", "yes", "positive", "pos", "y"}
FALSE_TOKENS = {"false", "0", "no", "negative", "neg", "n"}

# M1H activation cutoff — the ONLY numeric interaction/function threshold stated verbatim in the paper2
# STAR Methods (Lambourne et al., Mol Cell 2025, mmc12.pdf p.e12): "M1H signal >= 1 (activation) or
# <= -1 (repression)" (== >2-fold vs baseline). Y2H (PPI) and eY1H (PDI) calls are the authors' own
# Boolean columns (Y2H growth score >=2 vs auto-activator; eY1H manual 3-researcher consensus) and need
# no numeric cutoff here. N2H log2 NLR and luciferase Log2(FC) cutoffs are NOT stated in the metadata
# (only dashed lines in Figs S2I/S2K) -> per the binary-only policy they are not thresholded.
M1H_ACTIVATOR_MIN = 1.0     # M1H_mean >= this  -> activator
M1H_REPRESSOR_MAX = -1.0    # M1H_mean <= this  -> repressor


def ensure_dirs():
    """Create the gitignored output dirs if missing (safe; never touches input folders)."""
    for d in (DATA_DIR, ARTIFACTS_DIR, EXTERNAL_DIR):
        os.makedirs(d, exist_ok=True)
