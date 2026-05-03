"""
ppi_interface.py: BioGRID protein-protein interaction and domain interface features.

Builds per-protein feature vectors encoding:
- BioGRID interaction partner counts (physical vs. genetic)
- UniProt PPI binding region and phospho-PPI residue counts
- DNA binding domain annotation counts
- Region_priors topology features (extracellular, cytoplasmic, intramembrane)
- HPA surfaceome annotation

All features are derived from pre-built Daedalus Phase A interim tables and
the raw BioGRID TAB3 file.  CPU-only.
"""

from __future__ import annotations

import warnings
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np


class PPIInterfaceEncoder:
    """
    Compiles a 29-dimensional PPI/interface feature vector per protein.

    Sources:
    - biogrid_gene_interactions.tsv  (gene_name → partner/interaction counts)
    - uniprot_region_priors.tsv      (primary_accession → topology + phospho + DNA features)
    - hpa_localization.tsv           (ensembl_gene_id → surfaceome annotation)
    - uniprot_gencode_map.tsv        (primary_accession ↔ ENSP/ENST)
    - gencode_protein_reference.tsv  (ENSP → gene_name for BioGRID bridge)

    Feature layout [29]:
    [0]  biogrid_partner_count (log1p)
    [1]  biogrid_physical_partner_count (log1p)
    [2]  biogrid_interaction_count (log1p)
    [3]  biogrid_physical_interaction_count (log1p)
    [4]  uniprot_interaction_partner_count (log1p)          ← region_priors
    [5]  uniprot_ppi_binding_feature_count (log1p)          ← residues at PPI interfaces
    [6]  uniprot_phospho_interaction_feature_count (log1p)  ← phospho sites regulating PPIs
    [7]  uniprot_phospho_feature_count (log1p)              ← all phospho sites
    [8]  uniprot_extracellular_topology_count (log1p)
    [9]  uniprot_extracellular_topology_aa_fraction          ← extracellular AA / total
    [10] uniprot_cytoplasmic_topology_count (log1p)
    [11] uniprot_cytoplasmic_topology_aa_fraction
    [12] uniprot_intramembrane_count (log1p)
    [13] uniprot_tm_total_span (log1p)
    [14] uniprot_signal_count (log1p)
    [15] uniprot_signal_max_end (normalised 0-1 over 50 aa)
    [16] uniprot_dna_binding_feature_count (log1p)
    [17] uniprot_zinc_finger_feature_count (log1p)
    [18] uniprot_dna_region_feature_count (log1p)
    [19] uniprot_coiled_coil_count (log1p)
    [20] uniprot_glycosylation_count (log1p)
    [21] uniprot_disulfide_count (log1p)
    [22] uniprot_active_site_count (log1p)
    [23] uniprot_motif_count (log1p)
    [24] uniprot_nuclear_topology_count (log1p)
    [25] has_hpa_extracellular (binary)
    [26] is_hpa_membrane_or_surface (binary)
    [27] uniprot_modified_residue_count (log1p)
    [28] uniprot_domain_count (log1p)                       ← from region_priors
    """

    FEATURE_DIM: int = 29

    def __init__(self) -> None:
        # Lazy-loaded lookup dicts: ENSP → feature_dict
        self._biogrid: Dict[str, Dict] = {}
        self._region_priors: Dict[str, Dict] = {}
        self._hpa: Dict[str, Dict] = {}
        self._loaded = False

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(self, daedalus_interim_dir: Path) -> None:
        """
        Load all four source tables and build ENSP-keyed lookup dicts.

        Parameters
        ----------
        daedalus_interim_dir : Path
            Path to Daedalus Phase A data/interim/ directory.
        """
        import pandas as pd

        d = Path(daedalus_interim_dir)

        # --- BioGRID (gene_name → counts) → bridge via gencode protein reference ---
        biogrid_path = d / "biogrid_gene_interactions.tsv"
        gencode_prot_path = d / "gencode_protein_reference.tsv"
        if biogrid_path.exists() and gencode_prot_path.exists():
            bg = pd.read_csv(biogrid_path, sep="\t", low_memory=False)
            gp = pd.read_csv(gencode_prot_path, sep="\t", low_memory=False,
                             usecols=lambda c: c in ("species", "gene_id", "transcript_id",
                                                     "protein_id", "protein_length"))
            # gene_id → gene_name bridge via uniprot_gencode_map
            gene_map_path = d / "uniprot_gencode_map.tsv"
            gene_name_lookup: Dict[str, str] = {}
            if gene_map_path.exists():
                gmap = pd.read_csv(gene_map_path, sep="\t", low_memory=False,
                                   usecols=lambda c: c in ("gene_name", "gencode_gene_id",
                                                            "gencode_protein_id"))
                for _, row in gmap.iterrows():
                    pid = str(row.get("gencode_protein_id", ""))
                    gn  = str(row.get("gene_name", ""))
                    if pid and gn and pid != "nan":
                        gene_name_lookup[pid] = gn

            # Also build from gencode_protein_reference via gene_supervision_catalog if available
            sup_path = d / "gene_supervision_catalog.tsv"
            if sup_path.exists():
                try:
                    sup = pd.read_csv(sup_path, sep="\t", low_memory=False,
                                      usecols=lambda c: c in ("gene_id", "gene_name"))
                    gene_id_to_name = dict(zip(sup["gene_id"].astype(str),
                                               sup["gene_name"].astype(str)))
                except Exception:
                    gene_id_to_name = {}
            else:
                gene_id_to_name = {}

            bg_lookup: Dict[str, Dict] = {}
            for _, row in bg.iterrows():
                bg_lookup[str(row.get("gene_name", ""))] = row.to_dict()

            # Map ENSP → BioGRID row
            for _, row in gp.iterrows():
                pid = str(row.get("protein_id", ""))
                if not pid or pid == "nan":
                    continue
                gn = gene_name_lookup.get(pid, "")
                if not gn:
                    gid = str(row.get("gene_id", ""))
                    gn = gene_id_to_name.get(gid, "")
                if gn and gn in bg_lookup:
                    self._biogrid[pid] = bg_lookup[gn]

            print(f"[PPIInterfaceEncoder] BioGRID: {len(self._biogrid):,} ENSP entries mapped.")
        else:
            warnings.warn(f"BioGRID or gencode_protein_reference not found in {d}")

        # --- UniProt region priors (primary_accession → wide features) ---
        region_priors_path = d / "uniprot_region_priors.tsv"
        map_path = d / "uniprot_gencode_map.tsv"
        if region_priors_path.exists() and map_path.exists():
            rp = pd.read_csv(region_priors_path, sep="\t", low_memory=False)
            gmap = pd.read_csv(map_path, sep="\t", low_memory=False,
                               usecols=lambda c: c in ("primary_accession", "gencode_protein_id",
                                                       "gencode_transcript_id"))
            acc_lookup: Dict[str, Dict] = {
                str(row["primary_accession"]): row.to_dict()
                for _, row in rp.iterrows()
            }
            for _, mrow in gmap.iterrows():
                acc = str(mrow.get("primary_accession", ""))
                if acc not in acc_lookup:
                    continue
                feat = acc_lookup[acc]
                for id_col in ("gencode_protein_id", "gencode_transcript_id"):
                    gid = str(mrow.get(id_col, ""))
                    if gid and gid != "nan":
                        self._region_priors[gid] = feat

            print(f"[PPIInterfaceEncoder] Region priors: {len(self._region_priors):,} ENSP entries.")
        else:
            warnings.warn(f"uniprot_region_priors.tsv or gencode_map not found in {d}")

        # --- HPA surfaceome (ensembl_gene_id → extracellular annotation) ---
        hpa_path = d / "hpa_localization.tsv"
        if hpa_path.exists():
            hpa = pd.read_csv(hpa_path, sep="\t", low_memory=False)
            hpa_by_gene: Dict[str, Dict] = {}
            for _, row in hpa.iterrows():
                gid = str(row.get("ensembl_gene_id", ""))
                if gid:
                    hpa_by_gene[gid] = row.to_dict()

            # Bridge: ENSP → gene_id → HPA row via gencode_protein_reference
            if gencode_prot_path.exists():
                gp2 = pd.read_csv(gencode_prot_path, sep="\t", low_memory=False,
                                  usecols=lambda c: c in ("gene_id", "protein_id"))
                for _, row in gp2.iterrows():
                    pid = str(row.get("protein_id", ""))
                    gid_versioned = str(row.get("gene_id", ""))
                    gid_base = gid_versioned.split(".")[0]  # HPA uses unversioned IDs
                    hpa_row = hpa_by_gene.get(gid_versioned) or hpa_by_gene.get(gid_base)
                    if pid and hpa_row:
                        self._hpa[pid] = hpa_row

            print(f"[PPIInterfaceEncoder] HPA: {len(self._hpa):,} ENSP entries.")
        else:
            warnings.warn(f"hpa_localization.tsv not found in {d}")

        self._loaded = True

    # ------------------------------------------------------------------
    # Feature compilation
    # ------------------------------------------------------------------

    def compile_features(self, protein_id: str) -> np.ndarray:
        """
        Compile 29-dimensional PPI/interface feature vector for a protein.

        Parameters
        ----------
        protein_id : str
            Ensembl protein ID (versioned, e.g. ENSP00000360644.5).

        Returns
        -------
        np.ndarray
            Shape [29], dtype float32.
        """
        features = np.zeros(self.FEATURE_DIM, dtype=np.float32)

        bg  = self._biogrid.get(protein_id, {})
        rp  = self._region_priors.get(protein_id, {})
        hpa = self._hpa.get(protein_id, {})

        def f(d: dict, k: str) -> float:
            try:
                v = d.get(k, 0)
                if v is None:
                    return 0.0
                fv = float(v)
                return 0.0 if np.isnan(fv) else fv
            except (TypeError, ValueError):
                return 0.0

        # [0-3] BioGRID counts
        features[0]  = float(np.log1p(f(bg, "biogrid_partner_count")))
        features[1]  = float(np.log1p(f(bg, "biogrid_physical_partner_count")))
        features[2]  = float(np.log1p(f(bg, "biogrid_interaction_count")))
        features[3]  = float(np.log1p(f(bg, "biogrid_physical_interaction_count")))

        # [4-7] UniProt PPI / phospho
        features[4]  = float(np.log1p(f(rp, "interaction_partner_count")))
        features[5]  = float(np.log1p(f(rp, "ppi_binding_feature_count")))
        features[6]  = float(np.log1p(f(rp, "phospho_interaction_feature_count")))
        features[7]  = float(np.log1p(f(rp, "phospho_feature_count")))

        # [8-15] Topology (extracellular, cytoplasmic, intramembrane, signal, TM)
        total_topo = f(rp, "extracellular_topology_aa") + f(rp, "cytoplasmic_topology_aa") + 1e-8
        features[8]  = float(np.log1p(f(rp, "extracellular_topology_count")))
        features[9]  = float(
            min(f(rp, "extracellular_topology_aa") / (total_topo + f(rp, "cytoplasmic_topology_aa") + 1e-8), 1.0)
        )
        features[10] = float(np.log1p(f(rp, "cytoplasmic_topology_count")))
        features[11] = float(
            min(f(rp, "cytoplasmic_topology_aa") / (total_topo + 1e-8), 1.0)
        )
        features[12] = float(np.log1p(f(rp, "intramembrane_count")))
        features[13] = float(np.log1p(f(rp, "tm_total_span")))
        features[14] = float(np.log1p(f(rp, "signal_count")))
        features[15] = float(min(f(rp, "signal_max_end") / 50.0, 1.0))

        # [16-19] DNA binding features
        features[16] = float(np.log1p(f(rp, "dna_binding_feature_count")))
        features[17] = float(np.log1p(f(rp, "zinc_finger_feature_count")))
        features[18] = float(np.log1p(f(rp, "dna_region_feature_count")))
        features[19] = float(np.log1p(f(rp, "coiled_coil_count")))

        # [20-24] Additional modification / structural features
        features[20] = float(np.log1p(f(rp, "glycosylation_count")))
        features[21] = float(np.log1p(f(rp, "disulfide_count")))
        features[22] = float(np.log1p(f(rp, "active_site_count")))
        features[23] = float(np.log1p(f(rp, "motif_count")))
        features[24] = float(np.log1p(f(rp, "nuclear_topology_count")))

        # [25-26] HPA surfaceome
        features[25] = float(bool(hpa.get("has_extracellular_annotation", False)))
        features[26] = float(bool(hpa.get("is_membrane_or_surface_annotated", False)))

        # [27-28]
        features[27] = float(np.log1p(f(rp, "modified_residue_count")))
        features[28] = float(np.log1p(f(rp, "domain_count")))

        return features

    # ------------------------------------------------------------------
    # Batch extraction
    # ------------------------------------------------------------------

    def extract_and_cache(
        self,
        protein_ids: List[str],
        daedalus_interim_dir: Path,
        output_dir: Path,
    ) -> None:
        """
        Batch compile PPI/interface features and save to {protein_id}.pt files.

        Parameters
        ----------
        protein_ids : list[str]
            Protein identifiers (ENSP versioned IDs).
        daedalus_interim_dir : Path
            Path to Daedalus Phase A interim directory.
        output_dir : Path
            Directory to save .pt files.
        """
        from tqdm import tqdm

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self._loaded:
            self.load(Path(daedalus_interim_dir))

        to_process = []
        cached = 0
        for pid in protein_ids:
            out = output_dir / f"{pid}.pt"
            if out.exists():
                cached += 1
            else:
                to_process.append(pid)

        print(f"[PPIInterfaceEncoder] {cached} cached, {len(to_process)} to compile.")

        import torch
        for pid in tqdm(to_process, desc="Compiling PPI interface features"):
            out = output_dir / f"{pid}.pt"
            try:
                feats = self.compile_features(pid)
                torch.save(torch.from_numpy(feats), out)
            except Exception as exc:
                warnings.warn(f"Failed for {pid!r}: {exc}")
                torch.save(torch.zeros(self.FEATURE_DIM), out)
