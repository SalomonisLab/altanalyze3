from __future__ import annotations

from dataclasses import dataclass
import json
import shutil
from pathlib import Path
import re
from typing import Any, Dict, Iterable, Mapping, Optional
from urllib.request import urlopen

import pandas as pd

from .data_sources import ensure_artifact_dirs
from .training import prune_response_matrix


PACKAGE_DIR = Path(__file__).resolve().parent

CELLCHAT_URLS = {
    "human": "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.human.rda",
    "mouse": "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.mouse.rda",
}

NICHENET_LR_URLS = {
    "human": "https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds?download=1",
    "mouse": "https://zenodo.org/records/7074291/files/lr_network_mouse_21122021.rds?download=1",
}

NICHENET_LIGAND_TARGET_URLS = {
    "human": "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds?download=1",
    "mouse": "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds?download=1",
}


@dataclass(frozen=True)
class UpstreamResourceBuildParams:
    species: str
    output_dir: Optional[Path] = None
    raw_dir: Optional[Path] = None
    top_targets_per_signature: int = 250
    min_target_weight: float = 0.001
    include_nichenet_lr_network: bool = True
    overwrite: bool = False


def _import_rdata():
    try:
        import rdata  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only outside dev env
        raise ImportError(
            "Building fastComm upstream resources requires the optional 'rdata' package. "
            "Install it with `pip install rdata`."
        ) from exc
    return rdata


def _download_file(url: str, destination: Path, *, overwrite: bool = False) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        return destination
    with urlopen(url) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return destination


def _load_r_object(path: Path) -> Any:
    rdata = _import_rdata()
    parsed = rdata.parser.parse_file(str(path))
    return rdata.conversion.convert(parsed)


def _first_value(mapping: Mapping[Any, Any]) -> Any:
    return next(iter(mapping.values()))


def _clean_text(value: object) -> str:
    text = str(value or "").strip()
    return "" if text == "<NA>" else text


def _normalize_symbol_text(value: object) -> str:
    text = _clean_text(value)
    if not text:
        return ""
    parts = [part.strip() for part in re.split(r"[,+;|]", text) if part.strip()]
    return "+".join(parts)


def _complex_lookup(complex_table: pd.DataFrame) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for complex_name, row in complex_table.iterrows():
        parts = [_clean_text(value) for value in row.tolist()]
        members = [part for part in parts if part]
        if members:
            mapping[str(complex_name).strip()] = "+".join(members)
    return mapping


def _resolve_part(
    alias: object,
    symbol_text: object,
    complex_map: Mapping[str, str],
) -> str:
    normalized_symbol = _normalize_symbol_text(symbol_text)
    if normalized_symbol:
        return normalized_symbol
    alias_text = _clean_text(alias)
    if alias_text in complex_map:
        return complex_map[alias_text]
    return _normalize_symbol_text(alias_text)


def _as_dataframe(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        out = value.copy()
    else:
        out = pd.DataFrame(value)
    out.index = out.index.map(str)
    out.columns = out.columns.map(str)
    return out


def build_cellchat_lr_table(cellchat_db: Mapping[str, Any]) -> pd.DataFrame:
    interactions = _as_dataframe(cellchat_db["interaction"])
    complex_table = _as_dataframe(cellchat_db["complex"])
    complex_map = _complex_lookup(complex_table)

    rows = []
    for _, row in interactions.iterrows():
        ligand = _resolve_part(row.get("ligand", ""), row.get("ligand.symbol", ""), complex_map)
        receptor = _resolve_part(row.get("receptor", ""), row.get("receptor.symbol", ""), complex_map)
        if not ligand or not receptor:
            continue
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "pathway": _clean_text(row.get("pathway_name", "")),
                "evidence_weight": 1.0,
                "interaction_class": _clean_text(row.get("annotation", "")) or "cellchat_curated",
                "source": "CellChatDB",
                "source_detail": _clean_text(row.get("interaction_name", "")),
                "evidence": _clean_text(row.get("evidence", "")),
                "agonist": _clean_text(row.get("agonist", "")),
                "antagonist": _clean_text(row.get("antagonist", "")),
                "co_A_receptor": _clean_text(row.get("co_A_receptor", "")),
                "co_I_receptor": _clean_text(row.get("co_I_receptor", "")),
                "version": _clean_text(row.get("version", "")),
            }
        )
    out = pd.DataFrame(rows)
    out = out.drop_duplicates(subset=["ligand", "receptor", "pathway"], keep="first").reset_index(drop=True)
    return out


def build_nichenet_lr_table(lr_network: pd.DataFrame) -> pd.DataFrame:
    network = _as_dataframe(lr_network)
    rows = []
    for _, row in network.iterrows():
        ligand = _normalize_symbol_text(row.get("from", ""))
        receptor = _normalize_symbol_text(row.get("to", ""))
        if not ligand or not receptor:
            continue
        database = _clean_text(row.get("database", ""))
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "pathway": "",
                "evidence_weight": 0.9,
                "interaction_class": "nichenet_prior",
                "source": "NicheNet",
                "source_detail": database,
                "evidence": _clean_text(row.get("source", "")),
            }
        )
    return pd.DataFrame(rows).drop_duplicates(subset=["ligand", "receptor"], keep="first").reset_index(drop=True)


def merge_ligand_receptor_tables(
    cellchat_lr: pd.DataFrame,
    nichenet_lr: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    out = cellchat_lr.copy()
    if nichenet_lr is None or nichenet_lr.empty:
        return out.reset_index(drop=True)

    dedup = out.groupby(["ligand", "receptor"], sort=False)["source"].agg(lambda values: ",".join(sorted(set(values))))
    existing_pairs = set(dedup.index.tolist())
    novel = nichenet_lr.loc[
        ~nichenet_lr[["ligand", "receptor"]].apply(tuple, axis=1).isin(existing_pairs)
    ].copy()
    if novel.empty:
        return out.reset_index(drop=True)
    merged = pd.concat([out, novel], ignore_index=True, sort=False)
    return merged.reset_index(drop=True)


def build_nichenet_response_matrix(
    lr_table: pd.DataFrame,
    ligand_target_matrix: pd.DataFrame,
    *,
    top_targets_per_signature: int = 250,
    min_target_weight: float = 0.001,
) -> pd.DataFrame:
    matrix = _as_dataframe(ligand_target_matrix).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix.index = matrix.index.map(str)
    matrix.columns = matrix.columns.map(str)

    ligands = pd.Index(sorted(set(lr_table["ligand"].astype(str))))
    matched_ligands = [ligand for ligand in ligands if ligand in matrix.columns]
    if not matched_ligands:
        raise ValueError("No fastComm ligands overlapped the NicheNet ligand-target matrix")

    ligand_rows: Dict[str, pd.Series] = {}
    for ligand in matched_ligands:
        series = matrix[ligand].astype(float)
        series = series.loc[series > float(min_target_weight)]
        if series.empty:
            continue
        if top_targets_per_signature > 0:
            series = series.sort_values(ascending=False).head(top_targets_per_signature)
        ligand_rows[ligand] = series

    if not ligand_rows:
        raise ValueError("NicheNet ligand-target matrix did not yield any signatures after filtering")

    pathway_rows: Dict[str, pd.Series] = {}
    pathway_to_ligands = (
        lr_table.loc[lr_table["pathway"].astype(str).str.strip().ne(""), ["pathway", "ligand"]]
        .drop_duplicates()
        .groupby("pathway", sort=False)["ligand"]
        .agg(list)
    )
    for pathway, pathway_ligands in pathway_to_ligands.items():
        members = [ligand_rows[ligand] for ligand in pathway_ligands if ligand in ligand_rows]
        if not members:
            continue
        aggregated = pd.concat(members, axis=1).fillna(0.0).mean(axis=1)
        aggregated = aggregated.loc[aggregated > float(min_target_weight)]
        if aggregated.empty:
            continue
        if top_targets_per_signature > 0:
            aggregated = aggregated.sort_values(ascending=False).head(top_targets_per_signature)
        pathway_rows[str(pathway)] = aggregated

    all_rows = {**ligand_rows, **pathway_rows}
    all_genes = sorted({gene for series in all_rows.values() for gene in series.index})
    response = pd.DataFrame(0.0, index=list(all_rows), columns=all_genes)
    for signature, series in all_rows.items():
        response.loc[signature, series.index] = series.to_numpy(dtype=float)
    response = prune_response_matrix(
        response,
        top_genes=top_targets_per_signature,
        min_abs_score=min_target_weight,
        l2_normalize=True,
    )
    return response


def _bundle_paths(params: UpstreamResourceBuildParams) -> Dict[str, Path]:
    created = ensure_artifact_dirs(keys=["raw_downloads", "processed_tables"])
    raw_root = params.raw_dir or created["raw_downloads"] / "upstream"
    processed_root = params.output_dir or created["processed_tables"] / "upstream" / f"cellchat_nichenet_{params.species}"
    processed_root.mkdir(parents=True, exist_ok=True)
    raw_root.mkdir(parents=True, exist_ok=True)
    return {
        "raw_root": raw_root,
        "processed_root": processed_root,
        "lr_table": processed_root / "ligand_receptor.tsv",
        "response_matrix": processed_root / "response_signatures.tsv",
        "manifest": processed_root / "manifest.json",
    }


def _download_inputs(params: UpstreamResourceBuildParams, raw_root: Path) -> Dict[str, Path]:
    species = params.species
    paths = {
        "cellchat": raw_root / f"CellChatDB.{species}.rda",
        "nichenet_lr": raw_root / f"lr_network_{species}_21122021.rds",
        "nichenet_targets": raw_root / (
            "ligand_target_matrix_nsga2r_final.rds"
            if species == "human"
            else "ligand_target_matrix_nsga2r_final_mouse.rds"
        ),
    }
    _download_file(CELLCHAT_URLS[species], paths["cellchat"], overwrite=params.overwrite)
    if params.include_nichenet_lr_network:
        _download_file(NICHENET_LR_URLS[species], paths["nichenet_lr"], overwrite=params.overwrite)
    _download_file(NICHENET_LIGAND_TARGET_URLS[species], paths["nichenet_targets"], overwrite=params.overwrite)
    return paths


def build_upstream_resources(params: UpstreamResourceBuildParams) -> Dict[str, object]:
    species = params.species.strip().lower()
    if species not in {"human", "mouse"}:
        raise ValueError("species must be 'human' or 'mouse'")
    params = UpstreamResourceBuildParams(
        species=species,
        output_dir=params.output_dir,
        raw_dir=params.raw_dir,
        top_targets_per_signature=params.top_targets_per_signature,
        min_target_weight=params.min_target_weight,
        include_nichenet_lr_network=params.include_nichenet_lr_network,
        overwrite=params.overwrite,
    )

    paths = _bundle_paths(params)
    raw_paths = _download_inputs(params, paths["raw_root"])

    cellchat_obj = _load_r_object(raw_paths["cellchat"])
    cellchat_db = _first_value(cellchat_obj)
    cellchat_lr = build_cellchat_lr_table(cellchat_db)

    nichenet_lr = None
    if params.include_nichenet_lr_network and raw_paths["nichenet_lr"].exists():
        nichenet_lr = build_nichenet_lr_table(_load_r_object(raw_paths["nichenet_lr"]))
    merged_lr = merge_ligand_receptor_tables(cellchat_lr, nichenet_lr)

    ligand_target_matrix = _load_r_object(raw_paths["nichenet_targets"])
    response_matrix = build_nichenet_response_matrix(
        merged_lr,
        ligand_target_matrix,
        top_targets_per_signature=params.top_targets_per_signature,
        min_target_weight=params.min_target_weight,
    )

    merged_lr.to_csv(paths["lr_table"], sep="\t", index=False)
    response_matrix.to_csv(paths["response_matrix"], sep="\t")

    ligand_signature_names = set(merged_lr["ligand"].astype(str))
    manifest = {
        "species": species,
        "lr_table": str(paths["lr_table"]),
        "response_matrix": str(paths["response_matrix"]),
        "manifest": str(paths["manifest"]),
        "cellchat_url": CELLCHAT_URLS[species],
        "nichenet_lr_url": NICHENET_LR_URLS[species],
        "nichenet_ligand_target_url": NICHENET_LIGAND_TARGET_URLS[species],
        "n_cellchat_interactions": int(cellchat_lr.shape[0]),
        "n_merged_interactions": int(merged_lr.shape[0]),
        "n_ligand_signatures": int(sum(idx in ligand_signature_names for idx in response_matrix.index)),
        "n_total_signatures": int(response_matrix.shape[0]),
        "n_response_genes": int(response_matrix.shape[1]),
        "top_targets_per_signature": int(params.top_targets_per_signature),
        "min_target_weight": float(params.min_target_weight),
        "used_nichenet_lr_network": bool(params.include_nichenet_lr_network),
        "raw_files": {key: str(value) for key, value in raw_paths.items()},
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
