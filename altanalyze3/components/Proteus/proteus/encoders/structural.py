"""
StructuralFeatureEncoder: Compiles protein structural and functional features
from Daedalus Phase A UniProt annotation tables and sequence-based predictors.

Produces a 96-dimensional feature vector per protein encoding:
  - UniProt domain/PTM/TM features (from long-format features + ENSP bridge)
  - Sequence-predicted TM topology and signal peptide (novel isoform capable)
  - BioGRID PPI partner counts and UniProt PPI interface residues
  - Phospho-PPI regulation sites and DNA binding domain annotations
  - HPA surfaceome/extracellular annotation

The first 48 dimensions are backward-compatible with the v1 layout.
Dimensions 48-95 are new and populated from additional Daedalus tables.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch


class StructuralFeatureEncoder:
    """
    Compiles 96-dimensional structural/functional feature vectors.

    Feature layout [96]:
    --- UniProt domain/modification block [0-18] (unchanged from v1) ---
    [0]  n_domains (log1p, from uniprot_features long-format)
    [1]  n_active_sites (log1p)
    [2]  n_binding_sites (log1p)
    [3]  n_ptm_sites (log1p)
    [4]  n_tm_helices (log1p, from uniprot_features Transmembrane entries)
    [5]  has_signal_peptide (binary)
    [6]  is_kinase (binary, from pair TSV family_class)
    [7]  is_transcription_factor (binary)
    [8]  is_membrane_protein (binary, proxy: has TM helices)
    [9]  n_topological_domains (log1p)
    [10] has_coiled_coil (binary)
    [11] has_zinc_finger (binary)
    [12] has_dna_binding (binary)
    [13] protein_length_log1p  (from gencode_protein_sequence_features)
    [14] n_disulfide_bonds (log1p)
    [15] n_glycosylation_sites (log1p)
    [16] n_phosphorylation_sites (log1p)
    [17] n_ubiquitination_sites (log1p)
    [18] n_interaction_partners (log1p, from uniprot_region_priors)

    --- Sequence-predicted topology block [19-27] ---
    [19] seq_tm_count (log1p, Kyte-Doolittle predictor — works on novel isoforms)
    [20] seq_signal_candidate (binary)
    [21] seq_signal_score (normalised 0-1)
    [22] seq_tm_span_fraction (TM span / protein length)
    [23] seq_max_hydropathy_19 (normalised 0-1)
    [24] seq_n_terminal_hydropathy (normalised 0-1)
    [25] seq_glyco_motif_count (log1p)
    [26] seq_cysteine_fraction
    [27] seq_topology_type (0-5 / 5)

    --- BioGRID PPI block [28-37] ---
    [28] biogrid_partner_count (log1p)
    [29] biogrid_physical_partner_count (log1p)
    [30] biogrid_interaction_count (log1p)
    [31] biogrid_physical_interaction_count (log1p)
    [32] uniprot_ppi_binding_feature_count (log1p)
    [33] uniprot_phospho_interaction_feature_count (log1p)
    [34] uniprot_extracellular_topology_count (log1p)
    [35] uniprot_extracellular_topology_aa_fraction
    [36] uniprot_tm_total_span (log1p, from region_priors — most accurate)
    [37] uniprot_signal_count (log1p)

    --- Phospho / DNA-binding / surfaceome block [38-47] ---
    [38] uniprot_phospho_feature_count (log1p)
    [39] uniprot_dna_binding_feature_count (log1p)
    [40] uniprot_zinc_finger_feature_count (log1p)
    [41] uniprot_dna_region_feature_count (log1p)
    [42] uniprot_coiled_coil_count (log1p, from region_priors)
    [43] uniprot_cytoplasmic_topology_aa_fraction
    [44] uniprot_intramembrane_count (log1p)
    [45] has_hpa_extracellular (binary)
    [46] is_hpa_membrane_or_surface (binary)
    [47] uniprot_nuclear_topology_count (log1p)

    --- Extended block [48-95] (additional detail) ---
    [48] uniprot_motif_count (log1p)
    [49] uniprot_glycosylation_count (log1p, from region_priors)
    [50] uniprot_disulfide_count (log1p, from region_priors)
    [51] uniprot_active_site_count (log1p)
    [52] uniprot_domain_count (log1p, from region_priors)
    [53] uniprot_modified_residue_count (log1p)
    [54] uniprot_cytoplasmic_topology_count (log1p)
    [55] uniprot_signal_max_end (normalised 0-1 over 50 aa)
    [56] seq_topology_type_onehot_0 (soluble)
    [57] seq_topology_type_onehot_1 (single-pass type I)
    [58] seq_topology_type_onehot_2 (single-pass type II)
    [59] seq_topology_type_onehot_3 (single-pass type III)
    [60] seq_topology_type_onehot_4 (multi-pass)
    [61] seq_topology_type_onehot_5 (GPI proxy)
    [62-95] reserved (zeros)
    """

    FEATURE_DIM: int = 96

    def __init__(
        self,
        uniprot_features_path: Optional[Path] = None,
        uniprot_priors_path: Optional[Path] = None,
    ) -> None:
        """
        Parameters
        ----------
        uniprot_features_path : Path | None
            Path to uniprot_features.tsv from Daedalus Phase A interim.
        uniprot_priors_path : Path | None
            Path to uniprot_functional_priors.tsv from Daedalus Phase A interim.
        """
        self._features_df = None
        self._priors_df = None

        if uniprot_features_path is not None:
            self._load_features(Path(uniprot_features_path))
        if uniprot_priors_path is not None:
            self._load_priors(Path(uniprot_priors_path))

    def _load_features(self, path: Path) -> None:
        """Load uniprot_features.tsv."""
        try:
            import pandas as pd
            self._features_df = pd.read_csv(path, sep="\t", low_memory=False)
            print(f"[StructuralFeatureEncoder] Loaded features: {path}, {len(self._features_df)} rows.")
        except Exception as exc:
            warnings.warn(f"Failed to load uniprot features from {path}: {exc}", stacklevel=2)

    def _load_priors(self, path: Path) -> None:
        """Load uniprot_functional_priors.tsv."""
        try:
            import pandas as pd
            self._priors_df = pd.read_csv(path, sep="\t", low_memory=False)
            print(f"[StructuralFeatureEncoder] Loaded priors: {path}, {len(self._priors_df)} rows.")
        except Exception as exc:
            warnings.warn(f"Failed to load uniprot priors from {path}: {exc}", stacklevel=2)

    def _get_protein_row(self, protein_id: str):
        """Retrieve row for protein_id from features table."""
        if self._features_df is None:
            return None
        df = self._features_df
        # Try multiple ID columns
        for col in ("protein_id", "uniprot_id", "entry", "accession", "Entry"):
            if col in df.columns:
                mask = df[col] == protein_id
                if mask.any():
                    return df[mask].iloc[0]
        return None

    def _get_prior_row(self, protein_id: str, gene_id: str):
        """Retrieve row for protein/gene from priors table."""
        if self._priors_df is None:
            return None
        df = self._priors_df
        # Try by protein_id first, then gene_id
        for col, val in [("protein_id", protein_id), ("uniprot_id", protein_id),
                          ("gene_id", gene_id), ("gene", gene_id)]:
            if col in df.columns:
                mask = df[col] == val
                if mask.any():
                    return df[mask].iloc[0]
        return None

    @staticmethod
    def _safe_float(row, key: str, default: float = 0.0) -> float:
        """Safely extract a numeric value from a row dict/Series."""
        if row is None:
            return default
        try:
            val = row[key] if hasattr(row, "__getitem__") else getattr(row, key, default)
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return default
            return float(val)
        except (KeyError, AttributeError, TypeError, ValueError):
            return default

    @staticmethod
    def _safe_bool(row, key: str, true_values=("yes", "1", "true", "y")) -> float:
        """Safely extract a binary flag from a row."""
        if row is None:
            return 0.0
        try:
            val = row[key] if hasattr(row, "__getitem__") else getattr(row, key, 0)
            if isinstance(val, (bool, np.bool_)):
                return 1.0 if val else 0.0
            if isinstance(val, (int, float)):
                return 1.0 if val != 0 else 0.0
            return 1.0 if str(val).lower().strip() in true_values else 0.0
        except (KeyError, AttributeError, TypeError):
            return 0.0

    def compile_features(self, protein_id: str, gene_id: str = "") -> np.ndarray:
        """
        Compile the 48-dimensional structural feature vector for a protein.

        Parameters
        ----------
        protein_id : str
            UniProt accession or Ensembl protein ID.
        gene_id : str
            Ensembl gene ID (used as fallback lookup key).

        Returns
        -------
        np.ndarray
            Shape [48], dtype float32. Uses 0.0 for unavailable features.
        """
        feat_row = self._get_protein_row(protein_id)
        prior_row = self._get_prior_row(protein_id, gene_id)

        features = np.zeros(self.FEATURE_DIM, dtype=np.float32)

        # [0] n_domains
        features[0] = np.log1p(self._safe_float(feat_row, "n_domains"))
        # [1] n_active_sites
        features[1] = np.log1p(self._safe_float(feat_row, "n_active_sites"))
        # [2] n_binding_sites
        features[2] = np.log1p(self._safe_float(feat_row, "n_binding_sites"))
        # [3] n_ptm_sites
        n_ptm = (
            self._safe_float(feat_row, "n_ptm_sites")
            or self._safe_float(feat_row, "n_modified_residues")
        )
        features[3] = np.log1p(n_ptm)
        # [4] n_tm_helices
        features[4] = np.log1p(
            self._safe_float(feat_row, "n_tm_helices")
            or self._safe_float(feat_row, "tm_helix_count")
            or self._safe_float(feat_row, "n_transmembrane")
        )
        # [5] has_signal_peptide
        features[5] = self._safe_bool(feat_row, "has_signal_peptide") or \
                      self._safe_bool(feat_row, "signal_peptide")
        # [6] is_kinase
        features[6] = self._safe_bool(prior_row, "is_kinase") or \
                      self._safe_bool(feat_row, "is_kinase")
        # [7] is_transcription_factor
        features[7] = self._safe_bool(prior_row, "is_transcription_factor") or \
                      self._safe_bool(feat_row, "is_transcription_factor") or \
                      self._safe_bool(prior_row, "is_tf")
        # [8] is_membrane_protein
        features[8] = self._safe_bool(feat_row, "is_membrane_protein") or \
                      (features[4] > 0).astype(float)  # proxy: has TM helices
        # [9] n_topological_domains
        features[9] = np.log1p(self._safe_float(feat_row, "n_topological_domains"))
        # [10] has_coiled_coil
        features[10] = self._safe_bool(feat_row, "has_coiled_coil") or \
                       self._safe_bool(feat_row, "coiled_coil")
        # [11] has_zinc_finger
        features[11] = self._safe_bool(feat_row, "has_zinc_finger") or \
                       self._safe_bool(feat_row, "zinc_finger")
        # [12] has_dna_binding
        features[12] = self._safe_bool(feat_row, "has_dna_binding") or \
                       self._safe_bool(feat_row, "dna_binding") or \
                       float(features[7])  # TFs have DNA binding
        # [13] protein_length_log1p
        features[13] = np.log1p(self._safe_float(feat_row, "protein_length") or
                                self._safe_float(feat_row, "length") or
                                self._safe_float(feat_row, "sequence_length"))
        # [14] n_disulfide_bonds
        features[14] = np.log1p(self._safe_float(feat_row, "n_disulfide_bonds") or
                                self._safe_float(feat_row, "disulfide_bonds"))
        # [15] n_glycosylation_sites
        features[15] = np.log1p(self._safe_float(feat_row, "n_glycosylation_sites") or
                                self._safe_float(feat_row, "glycosylation"))
        # [16] n_phosphorylation_sites
        features[16] = np.log1p(self._safe_float(feat_row, "n_phosphorylation_sites") or
                                self._safe_float(feat_row, "phosphorylation"))
        # [17] n_ubiquitination_sites
        features[17] = np.log1p(self._safe_float(feat_row, "n_ubiquitination_sites") or
                                self._safe_float(feat_row, "ubiquitination"))
        # [18] n_interaction_partners
        features[18] = np.log1p(self._safe_float(prior_row, "n_interaction_partners") or
                                self._safe_float(feat_row, "n_interaction_partners") or
                                self._safe_float(prior_row, "ppi_degree"))

        # [19-28] reserved DeepTMHMM (zeros)
        # [29-38] reserved SignalP (zeros)
        # [39-47] reserved future (zeros)
        # Already initialized to zero

        return features

    @staticmethod
    def _build_ensp_lookup(
        daedalus_interim_dir: Path,
    ) -> Dict[str, Dict]:
        """
        Build a {ensp_id: wide_feature_dict} lookup by:
        1. Aggregating long-format uniprot_features.tsv per accession.
        2. Bridging accessions → ENSP/ENST IDs via uniprot_gencode_map.tsv.
        Also loads per-pair flags (is_kinase, is_tf) from isoform_pair_candidates if available.
        """
        import pandas as pd

        ensp_lookup: Dict[str, Dict] = {}

        feat_path = daedalus_interim_dir / "uniprot_features.tsv"
        map_path  = daedalus_interim_dir / "uniprot_gencode_map.tsv"

        if not feat_path.exists():
            return ensp_lookup

        feat_df = pd.read_csv(feat_path, sep="\t", low_memory=False)

        # Aggregate per-accession feature counts from long-format table
        agg: Dict[str, Dict] = {}
        for _, row in feat_df.iterrows():
            acc = str(row.get("primary_accession", ""))
            if not acc:
                continue
            if acc not in agg:
                agg[acc] = {
                    "n_domains": 0.0, "n_active_sites": 0.0, "n_binding_sites": 0.0,
                    "n_ptm_sites": 0.0, "n_tm_helices": 0.0, "has_signal_peptide": 0.0,
                    "n_topological_domains": 0.0, "has_coiled_coil": 0.0,
                    "has_zinc_finger": 0.0, "has_dna_binding": 0.0,
                    "n_disulfide_bonds": 0.0, "n_glycosylation_sites": 0.0,
                    "n_phosphorylation_sites": 0.0, "n_ubiquitination_sites": 0.0,
                    "is_kinase": 0.0, "is_transcription_factor": 0.0,
                    "n_interaction_partners": 0.0,
                }
            ftype = str(row.get("feature_type", "")).lower().strip()
            desc  = str(row.get("description", "")).lower()
            d = agg[acc]
            if ftype == "domain":
                d["n_domains"] += 1
            elif ftype == "active site":
                d["n_active_sites"] += 1
            elif ftype == "binding site":
                d["n_binding_sites"] += 1
            elif ftype == "modified residue":
                d["n_ptm_sites"] += 1
                if "phospho" in desc:
                    d["n_phosphorylation_sites"] += 1
                elif "ubiquitin" in desc:
                    d["n_ubiquitination_sites"] += 1
            elif ftype == "transmembrane":
                d["n_tm_helices"] += 1
            elif ftype == "signal peptide":
                d["has_signal_peptide"] = 1.0
            elif ftype == "topological domain":
                d["n_topological_domains"] += 1
            elif ftype == "coiled coil":
                d["has_coiled_coil"] = 1.0
            elif ftype == "zinc finger":
                d["has_zinc_finger"] = 1.0
            elif ftype == "dna binding":
                d["has_dna_binding"] = 1.0
            elif ftype == "disulfide bond":
                d["n_disulfide_bonds"] += 1
            elif ftype == "glycosylation":
                d["n_glycosylation_sites"] += 1

        if not map_path.exists():
            return ensp_lookup

        map_df = pd.read_csv(map_path, sep="\t", low_memory=False)
        for _, mrow in map_df.iterrows():
            acc = str(mrow.get("primary_accession", ""))
            if acc not in agg:
                continue
            feat_dict = dict(agg[acc])
            for id_col in ("gencode_protein_id", "gencode_transcript_id"):
                gid = str(mrow.get(id_col, ""))
                if gid and gid != "nan":
                    ensp_lookup[gid] = feat_dict

        print(f"[StructuralFeatureEncoder] Built ENSP lookup: {len(ensp_lookup)} entries "
              f"from {len(agg)} UniProt accessions.")
        return ensp_lookup

    def extract_and_cache(
        self,
        protein_ids: List[str],
        gene_ids: List[str],
        daedalus_interim_dir: Path,
        output_dir: Path,
    ) -> None:
        """
        Batch compile structural features and save to disk.

        Parameters
        ----------
        protein_ids : list[str]
            Protein identifiers (ENSP IDs or UniProt accessions).
        gene_ids : list[str]
            Corresponding gene identifiers.
        daedalus_interim_dir : Path
            Path to Daedalus Phase A interim directory.
        output_dir : Path
            Directory to save {protein_id}.pt files.
        """
        import pandas as pd
        from tqdm import tqdm

        daedalus_interim_dir = Path(daedalus_interim_dir)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Build ENSP → feature_dict lookup
        ensp_lookup = self._build_ensp_lookup(daedalus_interim_dir)

        # Also load per-pair annotation flags (kinase, TF) from pair TSV if available
        pair_flags: Dict[str, Dict] = {}
        pair_path = daedalus_interim_dir.parent / "processed" / "proteus_train.tsv"
        for split in ("train", "val", "test"):
            p = daedalus_interim_dir.parent / "processed" / f"proteus_{split}.tsv"
            if not p.exists():
                p = daedalus_interim_dir.parent.parent / "Proteus" / "data" / "processed" / f"proteus_{split}.tsv"
            if not p.exists():
                continue
            try:
                df = pd.read_csv(p, sep="\t", low_memory=False, usecols=lambda c: c in (
                    "reference_protein_id", "alternative_protein_id",
                    "family_class", "reference_is_kinase", "alternative_is_kinase",
                    "reference_is_transcription_factor", "alternative_is_transcription_factor",
                ))
                for _, row in df.iterrows():
                    is_kinase = str(row.get("family_class", "")).lower() == "kinase"
                    is_tf     = str(row.get("family_class", "")).lower() == "transcription_factor"
                    for col in ("reference_protein_id", "alternative_protein_id"):
                        pid_val = str(row.get(col, ""))
                        if pid_val and pid_val not in ("nan", "None"):
                            pair_flags[pid_val] = {
                                "is_kinase": float(is_kinase or bool(row.get(f"{col.split('_')[0]}_is_kinase", False))),
                                "is_transcription_factor": float(is_tf or bool(row.get(f"{col.split('_')[0]}_is_transcription_factor", False))),
                            }
            except Exception:
                pass

        # Load supplementary tables
        # 1. gencode_protein_sequence_features.tsv → ENSP-keyed TM/signal predictions
        seq_feat_lookup: Dict[str, Dict] = {}
        seq_feat_path = daedalus_interim_dir / "gencode_protein_sequence_features.tsv"
        if seq_feat_path.exists():
            sf = pd.read_csv(seq_feat_path, sep="\t", low_memory=False)
            for _, row in sf.iterrows():
                pid_val = str(row.get("protein_id", ""))
                if pid_val and pid_val != "nan":
                    seq_feat_lookup[pid_val] = row.to_dict()
            print(f"[StructuralFeatureEncoder] Loaded {len(seq_feat_lookup)} sequence topology entries.")

        # 2. uniprot_region_priors.tsv → ENSP-keyed region/PPI/topology features
        region_priors_lookup: Dict[str, Dict] = {}
        rp_path  = daedalus_interim_dir / "uniprot_region_priors.tsv"
        map_path2 = daedalus_interim_dir / "uniprot_gencode_map.tsv"
        if rp_path.exists() and map_path2.exists():
            rp = pd.read_csv(rp_path, sep="\t", low_memory=False)
            gmap2 = pd.read_csv(map_path2, sep="\t", low_memory=False,
                                usecols=lambda c: c in ("primary_accession",
                                                         "gencode_protein_id",
                                                         "gencode_transcript_id"))
            acc_rp: Dict[str, Dict] = {str(r["primary_accession"]): r.to_dict() for _, r in rp.iterrows()}
            for _, mrow in gmap2.iterrows():
                acc = str(mrow.get("primary_accession", ""))
                if acc not in acc_rp:
                    continue
                for id_col in ("gencode_protein_id", "gencode_transcript_id"):
                    gid2 = str(mrow.get(id_col, ""))
                    if gid2 and gid2 != "nan":
                        region_priors_lookup[gid2] = acc_rp[acc]
            print(f"[StructuralFeatureEncoder] Loaded {len(region_priors_lookup)} region_priors entries.")

        # 3. BioGRID (gene_name → counts) → bridge via gene_supervision_catalog
        biogrid_lookup: Dict[str, Dict] = {}
        bg_path  = daedalus_interim_dir / "biogrid_gene_interactions.tsv"
        sup_path = daedalus_interim_dir / "gene_supervision_catalog.tsv"
        gene_map3_path = daedalus_interim_dir / "uniprot_gencode_map.tsv"
        if bg_path.exists():
            bg_df = pd.read_csv(bg_path, sep="\t", low_memory=False)
            bg_by_name: Dict[str, Dict] = {str(r["gene_name"]): r.to_dict() for _, r in bg_df.iterrows()}
            # Build ENSP → gene_name via uniprot_gencode_map
            ensp_to_gene: Dict[str, str] = {}
            if gene_map3_path.exists():
                gmap3 = pd.read_csv(gene_map3_path, sep="\t", low_memory=False,
                                    usecols=lambda c: c in ("gene_name", "gencode_protein_id"))
                for _, row in gmap3.iterrows():
                    pid_v = str(row.get("gencode_protein_id", ""))
                    gn    = str(row.get("gene_name", ""))
                    if pid_v and gn and pid_v != "nan":
                        ensp_to_gene[pid_v] = gn
            for ensp, gn in ensp_to_gene.items():
                if gn in bg_by_name:
                    biogrid_lookup[ensp] = bg_by_name[gn]
            print(f"[StructuralFeatureEncoder] Loaded {len(biogrid_lookup)} BioGRID entries.")

        # 4. HPA surfaceome
        hpa_lookup: Dict[str, Dict] = {}
        hpa_path  = daedalus_interim_dir / "hpa_localization.tsv"
        gp_path   = daedalus_interim_dir / "gencode_protein_reference.tsv"
        if hpa_path.exists() and gp_path.exists():
            hpa_df = pd.read_csv(hpa_path, sep="\t", low_memory=False)
            gp_df  = pd.read_csv(gp_path, sep="\t", low_memory=False,
                                 usecols=lambda c: c in ("gene_id", "protein_id"))
            hpa_by_gid: Dict[str, Dict] = {str(r["ensembl_gene_id"]): r.to_dict() for _, r in hpa_df.iterrows()}
            for _, row in gp_df.iterrows():
                pid_v = str(row.get("protein_id", ""))
                gid_v = str(row.get("gene_id", ""))
                gid_base = gid_v.split(".")[0]  # HPA uses unversioned IDs
                hpa_row = hpa_by_gid.get(gid_v) or hpa_by_gid.get(gid_base)
                if pid_v and hpa_row:
                    hpa_lookup[pid_v] = hpa_row
            print(f"[StructuralFeatureEncoder] Loaded {len(hpa_lookup)} HPA entries.")

        _lookup    = ensp_lookup
        _seq       = seq_feat_lookup
        _rp        = region_priors_lookup
        _biogrid   = biogrid_lookup
        _hpa       = hpa_lookup

        to_process = []
        cached_count = 0
        for pid, gid in zip(protein_ids, gene_ids):
            out_path = output_dir / f"{pid}.pt"
            if out_path.exists():
                cached_count += 1
            else:
                to_process.append((pid, gid))

        print(
            f"[StructuralFeatureEncoder] {cached_count} already cached, "
            f"{len(to_process)} to compile."
        )

        from proteus.encoders.topology import (
            SequenceTopologyEncoder, predict_tm_segments,
            predict_signal_peptide, predict_topology_type,
        )
        _topo_enc = SequenceTopologyEncoder()

        def _f(d: dict, k: str) -> float:
            try:
                v = d.get(k, 0)
                if v is None:
                    return 0.0
                fv = float(v)
                return 0.0 if np.isnan(fv) else fv
            except (TypeError, ValueError):
                return 0.0

        for pid, gid in tqdm(to_process, desc="Compiling structural features"):
            out_path = output_dir / f"{pid}.pt"
            try:
                feat_dict = _lookup.get(pid, {})
                flags = pair_flags.get(pid, {})
                feat_dict = {**feat_dict, **flags}
                sf_row = _seq.get(pid, {})
                rp_row = _rp.get(pid, {})
                bg_row = _biogrid.get(pid, {})
                hpa_row = _hpa.get(pid, {})

                features = np.zeros(self.FEATURE_DIM, dtype=np.float32)

                # --- Block 0-18: UniProt domain/PTM ---
                if feat_dict:
                    features[0]  = np.log1p(feat_dict.get("n_domains", 0.0))
                    features[1]  = np.log1p(feat_dict.get("n_active_sites", 0.0))
                    features[2]  = np.log1p(feat_dict.get("n_binding_sites", 0.0))
                    features[3]  = np.log1p(feat_dict.get("n_ptm_sites", 0.0))
                    features[4]  = np.log1p(feat_dict.get("n_tm_helices", 0.0))
                    features[5]  = feat_dict.get("has_signal_peptide", 0.0)
                    features[6]  = feat_dict.get("is_kinase", 0.0)
                    features[7]  = feat_dict.get("is_transcription_factor", 0.0)
                    features[8]  = float(features[4] > 0)
                    features[9]  = np.log1p(feat_dict.get("n_topological_domains", 0.0))
                    features[10] = feat_dict.get("has_coiled_coil", 0.0)
                    features[11] = feat_dict.get("has_zinc_finger", 0.0)
                    features[12] = max(feat_dict.get("has_dna_binding", 0.0), features[7])
                    features[14] = np.log1p(feat_dict.get("n_disulfide_bonds", 0.0))
                    features[15] = np.log1p(feat_dict.get("n_glycosylation_sites", 0.0))
                    features[16] = np.log1p(feat_dict.get("n_phosphorylation_sites", 0.0))
                    features[17] = np.log1p(feat_dict.get("n_ubiquitination_sites", 0.0))

                # [13] protein length from gencode sequence features
                prot_len = _f(sf_row, "protein_length")
                features[13] = float(np.log1p(prot_len))

                # [18] interaction_partner_count from region_priors (most complete)
                features[18] = float(np.log1p(_f(rp_row, "interaction_partner_count")))

                # --- Block 19-27: Sequence-predicted topology ---
                features[19] = float(np.log1p(_f(sf_row, "predicted_tm_count")))
                features[20] = float(_f(sf_row, "predicted_signal_candidate"))
                raw_sig_score = _f(sf_row, "predicted_signal_score")
                features[21] = float(min(raw_sig_score / 3.0, 1.0))
                tm_span = _f(sf_row, "predicted_tm_total_span")
                features[22] = float(tm_span / prot_len) if prot_len > 0 else 0.0
                max_h = _f(sf_row, "predicted_tm_max_hydropathy_19")
                features[23] = float(np.clip((max_h + 4.5) / 9.0, 0, 1))
                nterm_h = _f(sf_row, "n_terminal_hydropathy_max8")
                features[24] = float(np.clip((nterm_h + 4.5) / 9.0, 0, 1))
                features[25] = float(np.log1p(_f(sf_row, "glyco_motif_count")))
                features[26] = float(_f(sf_row, "cysteine_fraction"))
                # topology_type encoded as normalised int (derived from TM count + position)
                tm_cnt = int(_f(sf_row, "predicted_tm_count"))
                # approximate type from count (without position info at this stage)
                topo_type = 0 if tm_cnt == 0 else (4 if tm_cnt >= 2 else 1)
                features[27] = float(topo_type / 5.0)

                # --- Block 28-37: BioGRID + region_priors topology ---
                features[28] = float(np.log1p(_f(bg_row, "biogrid_partner_count")))
                features[29] = float(np.log1p(_f(bg_row, "biogrid_physical_partner_count")))
                features[30] = float(np.log1p(_f(bg_row, "biogrid_interaction_count")))
                features[31] = float(np.log1p(_f(bg_row, "biogrid_physical_interaction_count")))
                features[32] = float(np.log1p(_f(rp_row, "ppi_binding_feature_count")))
                features[33] = float(np.log1p(_f(rp_row, "phospho_interaction_feature_count")))
                features[34] = float(np.log1p(_f(rp_row, "extracellular_topology_count")))
                ec_aa = _f(rp_row, "extracellular_topology_aa")
                cy_aa = _f(rp_row, "cytoplasmic_topology_aa")
                features[35] = float(min(ec_aa / (ec_aa + cy_aa + 1e-8), 1.0))
                features[36] = float(np.log1p(_f(rp_row, "tm_total_span")))
                features[37] = float(np.log1p(_f(rp_row, "signal_count")))

                # --- Block 38-47: Phospho / DNA / surfaceome ---
                features[38] = float(np.log1p(_f(rp_row, "phospho_feature_count")))
                features[39] = float(np.log1p(_f(rp_row, "dna_binding_feature_count")))
                features[40] = float(np.log1p(_f(rp_row, "zinc_finger_feature_count")))
                features[41] = float(np.log1p(_f(rp_row, "dna_region_feature_count")))
                features[42] = float(np.log1p(_f(rp_row, "coiled_coil_count")))
                features[43] = float(min(cy_aa / (ec_aa + cy_aa + 1e-8), 1.0))
                features[44] = float(np.log1p(_f(rp_row, "intramembrane_count")))
                features[45] = float(bool(hpa_row.get("has_extracellular_annotation", False)))
                features[46] = float(bool(hpa_row.get("is_membrane_or_surface_annotated", False)))
                features[47] = float(np.log1p(_f(rp_row, "nuclear_topology_count")))

                # --- Extended block 48-61 ---
                features[48] = float(np.log1p(_f(rp_row, "motif_count")))
                features[49] = float(np.log1p(_f(rp_row, "glycosylation_count")))
                features[50] = float(np.log1p(_f(rp_row, "disulfide_count")))
                features[51] = float(np.log1p(_f(rp_row, "active_site_count")))
                features[52] = float(np.log1p(_f(rp_row, "domain_count")))
                features[53] = float(np.log1p(_f(rp_row, "modified_residue_count")))
                features[54] = float(np.log1p(_f(rp_row, "cytoplasmic_topology_count")))
                features[55] = float(min(_f(rp_row, "signal_max_end") / 50.0, 1.0))
                # One-hot for topology type [56-61]
                features[56 + topo_type] = 1.0

                torch.save(torch.from_numpy(features), out_path)
            except Exception as exc:
                warnings.warn(f"Failed to compile features for {pid!r}: {exc}", stacklevel=2)
                torch.save(torch.zeros(self.FEATURE_DIM), out_path)

    @staticmethod
    def load_cached(protein_id: str, cache_dir: Path) -> Optional[np.ndarray]:
        """
        Load cached structural features.

        Parameters
        ----------
        protein_id : str
            Protein identifier.
        cache_dir : Path
            Directory containing {protein_id}.pt files.

        Returns
        -------
        np.ndarray | None
            Shape [96], or None if not found.
        """
        path = Path(cache_dir) / f"{protein_id}.pt"
        if not path.exists():
            return None
        tensor = torch.load(path, weights_only=True)
        return tensor.numpy().astype(np.float32)
