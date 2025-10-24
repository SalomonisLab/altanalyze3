"""
Resource management utilities for GO-Elite.

These helpers download, cache, and reload GO ontology DAGs together with
species-specific GOA gene associations so runs can reuse preprocessed data.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import shutil
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd

from .parser import ParsedGO, build_tree_with_annotations
from .structures import GOTermNode, GOTree

DEFAULT_CACHE_ENV = "ALTA_GOELITE_DATA"
DEFAULT_CACHE_ROOT = Path.home() / ".altanalyze3" / "goelite"

GO_ONTOLOGY_URL = "https://current.geneontology.org/ontology/go-basic.obo.gz"

SPECIES_CONFIG: Dict[str, Dict[str, str]] = {
    "human": {
        "gaf_url": "https://current.geneontology.org/annotations/goa_human.gaf.gz",
        "gaf_filename": "goa_human.gaf.gz",
    },
    "mouse": {
        "gaf_url": "https://current.geneontology.org/annotations/goa_mouse.gaf.gz",
        "gaf_filename": "goa_mouse.gaf.gz",
    },
}

AVAILABLE_SPECIES = tuple(sorted(SPECIES_CONFIG.keys()))


def resolve_cache_dir(cache_dir: Optional[str] = None) -> Path:
    """
    Determine the root cache directory, creating it if necessary.
    """

    if cache_dir:
        base = Path(cache_dir).expanduser()
    elif os.getenv(DEFAULT_CACHE_ENV):
        base = Path(os.environ[DEFAULT_CACHE_ENV]).expanduser()
    else:
        base = DEFAULT_CACHE_ROOT

    base.mkdir(parents=True, exist_ok=True)
    return base


def _download(url: str, dest: Path, chunk_size: int = 1 << 18) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, open(dest, "wb") as handle:
        shutil.copyfileobj(response, handle, length=chunk_size)


def _copy_or_download(source_path: Optional[str], url: Optional[str], dest: Path) -> Path:
    if source_path:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(Path(source_path), dest)
        return dest
    if not url:
        raise ValueError("Either a source path or URL must be provided for download.")
    _download(url, dest)
    return dest


def _open_text(path: Path) -> io.TextIOWrapper:
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")


def _sanitize_version(value: str) -> str:
    clean = value.strip().replace(" ", "_").replace("/", "_")
    return clean or datetime.utcnow().strftime("%Y%m%d")


def _extract_obo_version(path: Path) -> str:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("data-version:"):
                return _sanitize_version(line.split("data-version:")[1])
    return datetime.utcnow().strftime("%Y%m%d")


def _extract_gaf_version(path: Path) -> str:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("!date:"):
                return _sanitize_version(line.split("!date:")[1])
            if line.startswith("!Generated:"):
                return _sanitize_version(line.split("!Generated:")[1])
    return datetime.utcnow().strftime("%Y%m%d")


def _write_tree_json(tree: GOTree, path: Path) -> None:
    payload = {
        "roots": list(tree.roots),
        "namespace_index": {ns: list(ids) for ns, ids in tree.namespace_index.items()},
        "nodes": [
            {
                "term_id": node.term_id,
                "name": node.name,
                "namespace": node.namespace,
                "parents": list(node.parents),
                "children": list(node.children),
                "depth": node.depth,
            }
            for node in tree.nodes.values()
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _read_tree_json(path: Path) -> GOTree:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    nodes = {
        entry["term_id"]: GOTermNode(
            term_id=entry["term_id"],
            name=entry["name"],
            namespace=entry["namespace"],
            parents=tuple(entry.get("parents", [])),
            children=tuple(entry.get("children", [])),
            depth=int(entry.get("depth", 0)),
        )
        for entry in payload["nodes"]
    }
    return GOTree(
        nodes=nodes,
        roots=tuple(payload.get("roots", [])),
        namespace_index={ns: tuple(ids) for ns, ids in payload.get("namespace_index", {}).items()},
    )


def _write_term_gene_table(mapping: Dict[str, Iterable[str]], base_path: Path) -> Tuple[Path, str]:
    records = [
        {"term_id": term_id, "gene_symbol": str(gene)}
        for term_id, genes in mapping.items()
        for gene in genes
    ]
    df = pd.DataFrame.from_records(records).drop_duplicates()
    parquet_path = base_path.with_suffix(".parquet")
    try:
        df.to_parquet(parquet_path, index=False)
        return parquet_path, "parquet"
    except Exception:
        csv_path = base_path.with_suffix(".tsv.gz")
        df.to_csv(csv_path, sep="\t", index=False, compression="gzip")
        return csv_path, "tsv.gz"


def _read_term_gene_table(base_path: Path) -> Dict[str, set]:
    parquet_path = base_path.with_suffix(".parquet")
    if parquet_path.exists():
        try:
            df = pd.read_parquet(parquet_path)
        except Exception as exc:  # pragma: no cover
            raise RuntimeError(f"Unable to read cached Parquet table at {parquet_path}: {exc}") from exc
    else:
        csv_path = base_path.with_suffix(".tsv.gz")
        if not csv_path.exists():
            raise FileNotFoundError(f"Gene table not found for cache base {base_path}")
        df = pd.read_csv(csv_path, sep="\t")

    mapping: Dict[str, set] = {}
    for term_id, group in df.groupby("term_id"):
        mapping[str(term_id)] = set(group["gene_symbol"].astype(str))
    return mapping


def prepare_species_resources(
    species: str,
    *,
    cache_dir: Optional[str] = None,
    version: Optional[str] = None,
    obo_path: Optional[str] = None,
    gaf_path: Optional[str] = None,
    force: bool = False,
    acceptable_evidence: Optional[set[str]] = None,
) -> ParsedGO:
    cache_root = resolve_cache_dir(cache_dir)
    species_key = species.lower()
    if not force:
        try:
            return load_cached_resources(species_key, cache_dir=cache_root, version=version)
        except FileNotFoundError:
            pass

    return _build_species_resources(
        species_key,
        cache_root=cache_root,
        user_version=version,
        obo_path=obo_path,
        gaf_path=gaf_path,
        acceptable_evidence=acceptable_evidence,
    )


def load_cached_resources(
    species: str,
    *,
    cache_dir: Optional[str] = None,
    version: Optional[str] = None,
) -> ParsedGO:
    cache_root = resolve_cache_dir(cache_dir)
    species_key = species.lower()
    species_dir = cache_root / species_key

    if not species_dir.exists():
        raise FileNotFoundError(f"No GO-Elite cache found for species '{species_key}' at {species_dir}")

    if version is None:
        latest_path = species_dir / "latest.json"
        if not latest_path.exists():
            raise FileNotFoundError(f"No cached version metadata found for species '{species_key}'.")
        with open(latest_path, "r", encoding="utf-8") as handle:
            version = json.load(handle)["version"]

    version_dir = species_dir / version
    if not version_dir.exists():
        raise FileNotFoundError(f"Cached version '{version}' not found for species '{species_key}'.")

    data_dir = version_dir / "data"
    tree_path = data_dir / "go_tree.json"
    term_gene_base = data_dir / "term_gene"

    tree = _read_tree_json(tree_path)
    term_to_genes = _read_term_gene_table(term_gene_base)

    for term_id, genes in term_to_genes.items():
        if term_id in tree.nodes:
            tree.nodes[term_id] = tree.nodes[term_id].with_genes(genes)

    return ParsedGO(tree=tree, term_to_genes=term_to_genes)


def _build_species_resources(
    species: str,
    *,
    cache_root: Path,
    user_version: Optional[str],
    obo_path: Optional[str],
    gaf_path: Optional[str],
    acceptable_evidence: Optional[set[str]],
) -> ParsedGO:
    if species not in SPECIES_CONFIG:
        raise ValueError(f"Unsupported species '{species}'. Available: {sorted(SPECIES_CONFIG)}")

    species_dir = cache_root / species
    species_dir.mkdir(parents=True, exist_ok=True)

    working_tag = _sanitize_version(user_version or datetime.utcnow().strftime("%Y%m%d%H%M%S"))
    working_dir = species_dir / f"{working_tag}_build"
    if working_dir.exists():
        shutil.rmtree(working_dir)
    temp_dir = working_dir / "downloads"
    temp_dir.mkdir(parents=True, exist_ok=True)

    obo_target_name = Path(obo_path).name if obo_path else "go-basic.obo.gz"
    gaf_target_name = Path(gaf_path).name if gaf_path else SPECIES_CONFIG[species]["gaf_filename"]

    obo_target = temp_dir / obo_target_name
    gaf_target = temp_dir / gaf_target_name

    obo_source_url = None if obo_path else GO_ONTOLOGY_URL
    gaf_source_url = None if gaf_path else SPECIES_CONFIG[species]["gaf_url"]

    _copy_or_download(obo_path, obo_source_url, obo_target)
    _copy_or_download(gaf_path, gaf_source_url, gaf_target)

    parsed = build_tree_with_annotations(str(obo_target), str(gaf_target), acceptable_evidence=acceptable_evidence)

    obo_version = _extract_obo_version(obo_target)
    gaf_version = _extract_gaf_version(gaf_target)

    version = _sanitize_version(user_version or gaf_version or obo_version)
    version_dir = species_dir / version
    if version_dir.exists():
        shutil.rmtree(version_dir)
    shutil.move(str(working_dir), str(version_dir))

    temp_dir = version_dir / "downloads"
    obo_target = temp_dir / obo_target_name
    gaf_target = temp_dir / gaf_target_name

    data_dir = version_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    tree_path = data_dir / "go_tree.json"
    term_gene_base = data_dir / "term_gene"

    _write_tree_json(parsed.tree, tree_path)
    term_table_path, storage_format = _write_term_gene_table(parsed.term_to_genes, term_gene_base)

    metadata = {
        "species": species,
        "version": version,
        "built_at": datetime.utcnow().isoformat(),
        "source": {
            "obo": {
                "path": str(obo_target),
                "url": obo_source_url or "local",
                "release": obo_version,
            },
            "gaf": {
                "path": str(gaf_target),
                "url": gaf_source_url or "local",
                "release": gaf_version,
            },
        },
        "storage": {
            "tree": str(tree_path.relative_to(version_dir)),
            "term_gene": str(term_table_path.relative_to(version_dir)),
            "term_gene_format": storage_format,
        },
        "counts": {
            "terms": len(parsed.tree.nodes),
            "annotations": sum(len(genes) for genes in parsed.term_to_genes.values()),
        },
    }

    metadata_path = version_dir / "metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    latest_path = species_dir / "latest.json"
    with open(latest_path, "w", encoding="utf-8") as handle:
        json.dump({"version": version}, handle, indent=2)

    versions_path = species_dir / "versions.json"
    versions = []
    if versions_path.exists():
        with open(versions_path, "r", encoding="utf-8") as handle:
            versions = json.load(handle)
    versions = [entry for entry in versions if entry.get("version") != version]
    versions.append({"version": version, "built_at": metadata["built_at"]})
    with open(versions_path, "w", encoding="utf-8") as handle:
        json.dump(sorted(versions, key=lambda x: x["version"]), handle, indent=2)

    return parsed
