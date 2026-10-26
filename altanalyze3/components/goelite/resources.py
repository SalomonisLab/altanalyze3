"""
Resource management utilities for GO-Elite.

These helpers download, cache, and reload GO ontology DAGs together with
species-specific GOA gene associations so runs can reuse preprocessed data.
"""

from __future__ import annotations

import gzip
import io
import json
import logging
import os
import re
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
BUNDLED_CACHE_ROOT = Path(__file__).resolve().parent.parent / "cellHarmony" / "goelite_cache"

GO_ONTOLOGY_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

SPECIES_CACHE_ALIAS: Dict[str, str] = {
    "human": "Hs",
    "mouse": "Mm",
}

# Default build protocol (Nathan Salomonis, 2026-10-24): the annotation of a species is the UNION of the
# UniProt-GOA GAF and NCBI gene2go (GeneIDs mapped to symbols through NCBI gene_info), PLUS the same union
# for the ortholog species carried over through NCBI gene_orthologs (ToppFun does this: IKZF1, IGHM and LYST
# reach "regulation of leukocyte activation" only through their mouse orthologs), then propagated to parent
# terms over is_a, part_of and regulates. When a source cannot be obtained the build proceeds without it and
# records the problem in metadata.json. Nathan's benchmark for the human build is ToppFun (2026-10-24).
NCBI_GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
NCBI_GENE_ORTHOLOGS_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
SPECIES_CONFIG: Dict[str, Dict[str, str]] = {
    "human": {
        "gaf_url": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
        "gaf_filename": "goa_human.gaf.gz",
        "taxon_id": "9606",
        "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        "gene_info_filename": "Homo_sapiens.gene_info.gz",
        "ortholog_species": "mouse",
    },
    "mouse": {
        "gaf_url": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz",
        "gaf_filename": "goa_mouse.gaf.gz",
        "taxon_id": "10090",
        "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
        "gene_info_filename": "Mus_musculus.gene_info.gz",
        "ortholog_species": "human",
    },
}

AVAILABLE_SPECIES = tuple(sorted(SPECIES_CONFIG.keys()))


def _species_cache_name(species: str) -> str:
    return SPECIES_CACHE_ALIAS.get(species.lower(), species.lower())


def candidate_cache_dirs(cache_dir: Optional[str] = None) -> Tuple[Path, ...]:
    """
    Return cache roots to search for prebuilt GO-Elite assets.

    The bundled cellHarmony cache is preferred by default so Docker images can
    use repository-local data without attempting a download. If an explicit
    cache_dir is supplied, that path is treated as the authoritative cache.
    """

    roots = []
    if cache_dir:
        roots.append(Path(cache_dir).expanduser())
    else:
        roots.append(BUNDLED_CACHE_ROOT)
        if os.getenv(DEFAULT_CACHE_ENV):
            roots.append(Path(os.environ[DEFAULT_CACHE_ENV]).expanduser())
        else:
            roots.append(DEFAULT_CACHE_ROOT)

    seen = set()
    ordered = []
    for root in roots:
        key = str(root.expanduser())
        if key in seen:
            continue
        seen.add(key)
        ordered.append(root)
    return tuple(ordered)


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


def resolve_species_cache_dir(species: str, cache_dir: Optional[str] = None) -> Path:
    cache_root = resolve_cache_dir(cache_dir)
    cache_name = _species_cache_name(species)
    return cache_root / cache_name


def _resolve_existing_species_dir(species: str, cache_root: Path) -> Optional[Path]:
    species_key = species.lower()
    species_dir = cache_root / _species_cache_name(species_key)
    if species_dir.exists():
        return species_dir
    legacy_dir = cache_root / species_key
    if legacy_dir.exists():
        return legacy_dir
    return None


def _download(url: str, dest: Path, chunk_size: int = 1 << 18) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    # purl.obolibrary.org and the NCBI FTP mirror refuse Python's default User-Agent with HTTP 403.
    request = urllib.request.Request(url, headers={"User-Agent": "AltAnalyze3-GO-Elite/1.0 (+https://altanalyze.org)"})
    with urllib.request.urlopen(request) as response, open(dest, "wb") as handle:
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
            low = line.lower()
            for key in ("!date-generated:", "!date:", "!generated:"):
                if low.startswith(key):
                    # keep the date part only (e.g. "2026-10-11 13:43" -> 20261011)
                    return _sanitize_version(re.sub(r"[^0-9A-Za-z]", "", line.split(":", 1)[1].strip().split()[0]))
            if not line.startswith("!"):
                break
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
                "regulates": list(node.regulates),
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
            regulates=tuple(entry.get("regulates", [])),
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
    gene2go_path: Optional[str] = None,
    gene_info_path: Optional[str] = None,
    use_gene2go: bool = True,
    use_orthologs: bool = True,
    use_regulates: bool = True,
    gene_orthologs_path: Optional[str] = None,
    ortholog_gaf_path: Optional[str] = None,
    ortholog_gene_info_path: Optional[str] = None,
) -> ParsedGO:
    species_key = species.lower()
    if not force:
        try:
            return load_cached_resources(species_key, cache_dir=cache_dir, version=version)
        except FileNotFoundError:
            pass

    return _build_species_resources(
        species_key,
        cache_root=resolve_cache_dir(cache_dir),
        user_version=version,
        obo_path=obo_path,
        gaf_path=gaf_path,
        acceptable_evidence=acceptable_evidence,
        gene2go_path=gene2go_path,
        gene_info_path=gene_info_path,
        use_gene2go=use_gene2go,
        use_orthologs=use_orthologs,
        use_regulates=use_regulates,
        gene_orthologs_path=gene_orthologs_path,
        ortholog_gaf_path=ortholog_gaf_path,
        ortholog_gene_info_path=ortholog_gene_info_path,
    )


def load_cached_resources(
    species: str,
    *,
    cache_dir: Optional[str] = None,
    version: Optional[str] = None,
) -> ParsedGO:
    species_key = species.lower()
    last_error: Optional[FileNotFoundError] = None

    # Collect every usable cache; prefer the richest protocol: ortholog transfer over none, the union
    # (GOA + gene2go) over a single source, then keep the root order (bundled first, then the user cache).
    found: list = []
    for cache_root in candidate_cache_dirs(cache_dir):
        species_dir = _resolve_existing_species_dir(species_key, cache_root)
        if species_dir is None:
            continue
        # Candidate versions of this root: the requested one, else latest.json first and then every
        # version listed in versions.json (a richer protocol built earlier still wins over a later plain one).
        candidates: list = []
        if version is not None:
            candidates.append(version)
        else:
            latest_path = species_dir / "latest.json"
            if not latest_path.exists():
                last_error = FileNotFoundError(
                    f"No cached version metadata found for species '{species_key}' in {species_dir}."
                )
                continue
            with open(latest_path, "r", encoding="utf-8") as handle:
                candidates.append(json.load(handle)["version"])
            versions_path = species_dir / "versions.json"
            if versions_path.exists():
                try:
                    with open(versions_path, "r", encoding="utf-8") as handle:
                        listed = [entry.get("version") for entry in json.load(handle)]
                except Exception:
                    listed = []
                candidates.extend(v for v in listed if v and v not in candidates)
        root_index = len({f[3].parent for f in found}) if found else 0
        for rank_in_root, selected_version in enumerate(candidates):
            version_dir = species_dir / selected_version
            if not version_dir.exists() or not (version_dir / "data").exists():
                if rank_in_root == 0:
                    last_error = FileNotFoundError(
                        f"Cached version '{selected_version}' not found for species '{species_key}' in {species_dir}."
                    )
                continue
            meta_path = version_dir / "metadata.json"
            protocol, propagation = "", "is_a,part_of"
            if meta_path.exists():
                try:
                    meta = json.load(open(meta_path))
                    protocol = meta.get("annotation_protocol", "")
                    propagation = meta.get("propagation", propagation)
                except Exception:
                    protocol = ""
            found.append((0 if "+ortholog(" in protocol else 1, 0 if protocol.startswith("union") else 1,
                          (root_index, rank_in_root), version_dir, protocol, propagation))
    for _, _, _, version_dir, protocol, propagation in sorted(found):
        logging.getLogger("go_elite").info("GO-Elite cache: %s (%s)", version_dir, protocol or "GAF only, pre-union build")
        data_dir = version_dir / "data"
        tree_path = data_dir / "go_tree.json"
        term_gene_base = data_dir / "term_gene"

        tree = _read_tree_json(tree_path)
        term_to_genes = _read_term_gene_table(term_gene_base)
        # Caches built before 2026-10-24 hold direct annotations only; propagation is idempotent, and the
        # relations used are the ones recorded at build time (old caches: is_a/part_of).
        from .parser import propagate_annotations
        term_to_genes = propagate_annotations(tree, term_to_genes, use_regulates="regulates" in propagation)

        for term_id, genes in term_to_genes.items():
            if term_id in tree.nodes:
                tree.nodes[term_id] = tree.nodes[term_id].with_genes(genes)

        parsed = ParsedGO(tree=tree, term_to_genes=term_to_genes)
        parsed.version_dir = str(version_dir)  # type: ignore[attr-defined]
        return parsed

    if last_error is not None:
        raise last_error
    searched = ", ".join(str(path) for path in candidate_cache_dirs(cache_dir))
    raise FileNotFoundError(
        f"No GO-Elite cache found for species '{species_key}'. Checked: {searched}"
    )


def _set_aside(path: Path) -> None:
    """Never remove a cache directory: move it under delete/ with a timestamp."""
    if not path.exists():
        return
    trash = path.parent / "delete"
    trash.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(trash / f"{path.name}_{datetime.utcnow().strftime('%Y%m%d%H%M%S')}"))


def _build_species_resources(
    species: str,
    *,
    cache_root: Path,
    user_version: Optional[str],
    obo_path: Optional[str],
    gaf_path: Optional[str],
    acceptable_evidence: Optional[set[str]],
    gene2go_path: Optional[str] = None,
    gene_info_path: Optional[str] = None,
    use_gene2go: bool = True,
    use_orthologs: bool = True,
    use_regulates: bool = True,
    gene_orthologs_path: Optional[str] = None,
    ortholog_gaf_path: Optional[str] = None,
    ortholog_gene_info_path: Optional[str] = None,
) -> ParsedGO:
    if species not in SPECIES_CONFIG:
        raise ValueError(f"Unsupported species '{species}'. Available: {sorted(SPECIES_CONFIG)}")
    from .parser import (parse_gaf_synonyms, parse_gene_orthologs, transfer_ortholog_annotations,
                         build_tree_with_annotation_frames)
    logger = logging.getLogger("go_elite")
    cfg = SPECIES_CONFIG[species]

    species_dir = cache_root / _species_cache_name(species)
    species_dir.mkdir(parents=True, exist_ok=True)

    working_tag = _sanitize_version(user_version or datetime.utcnow().strftime("%Y%m%d%H%M%S"))
    working_dir = species_dir / f"{working_tag}_build"
    _set_aside(working_dir)
    temp_dir = working_dir / "downloads"
    temp_dir.mkdir(parents=True, exist_ok=True)

    obo_target_name = Path(obo_path).name if obo_path else Path(GO_ONTOLOGY_URL).name or "go-basic.obo"
    gaf_target_name = Path(gaf_path).name if gaf_path else cfg["gaf_filename"]
    obo_target = temp_dir / obo_target_name
    gaf_target = temp_dir / gaf_target_name
    obo_source_url = None if obo_path else GO_ONTOLOGY_URL
    gaf_source_url = None if gaf_path else cfg["gaf_url"]
    _copy_or_download(obo_path, obo_source_url, obo_target)

    sources: Dict[str, Dict[str, str]] = {"obo": {"path": str(obo_target), "url": obo_source_url or "local"}}
    problems: list = []
    own = _collect_species_annotations(species, cfg, temp_dir=temp_dir, cache_root=cache_root, gaf_path=gaf_path,
                                       gaf_target=gaf_target, gene_info_path=gene_info_path, gene2go_path=gene2go_path,
                                       acceptable_evidence=acceptable_evidence, use_gene2go=use_gene2go, logger=logger)
    frames = dict(own["frames"])
    sources.update(own["sources"])
    problems.extend(own["problems"])
    gene_info_synonyms: Dict[str, str] = own["gene_info_synonyms"]
    if not frames:
        raise RuntimeError("GO-Elite build failed: no annotation source could be obtained: " + "; ".join(problems))
    protocol = own["protocol"]

    # Source 3: the ortholog species' annotations (same union protocol) carried over through NCBI gene_orthologs.
    ortholog_meta: Dict[str, object] = {"enabled": bool(use_orthologs)}
    other_species = cfg.get("ortholog_species")
    if use_orthologs and other_species and other_species in SPECIES_CONFIG:
        other_cfg = SPECIES_CONFIG[other_species]
        try:
            if not own["gene_id_to_symbol"]:
                raise RuntimeError(f"{species} gene_info unavailable, GeneIDs cannot be mapped to symbols")
            other_dir = temp_dir / f"ortholog_{other_species}"
            other_dir.mkdir(parents=True, exist_ok=True)
            other = _collect_species_annotations(
                other_species, other_cfg, temp_dir=other_dir, cache_root=cache_root, gaf_path=ortholog_gaf_path,
                gaf_target=other_dir / (Path(ortholog_gaf_path).name if ortholog_gaf_path else other_cfg["gaf_filename"]),
                gene_info_path=ortholog_gene_info_path, gene2go_path=gene2go_path,
                acceptable_evidence=acceptable_evidence, use_gene2go=use_gene2go, logger=logger)
            problems.extend(f"ortholog {other_species}: {p}" for p in other["problems"])
            if not other["frames"]:
                raise RuntimeError(f"no {other_species} annotation source could be obtained")
            if not other["gene_id_to_symbol"]:
                raise RuntimeError(f"{other_species} gene_info unavailable, GeneIDs cannot be mapped to symbols")
            # gene_orthologs covers every taxon (about 130 MB); it is kept once at the cache root and reused.
            if gene_orthologs_path:
                orth_target = temp_dir / Path(gene_orthologs_path).name
                _copy_or_download(gene_orthologs_path, None, orth_target)
            else:
                orth_target = cache_root / "downloads" / "gene_orthologs.gz"
                if not orth_target.exists():
                    _download(NCBI_GENE_ORTHOLOGS_URL, orth_target)
            orthologs = parse_gene_orthologs(str(orth_target), int(cfg["taxon_id"]), int(other_cfg["taxon_id"]),
                                             own["gene_id_to_symbol"], other["gene_id_to_symbol"])
            unmapped_pairs = orthologs.pop("__unmapped__", {"0"})
            other_union = pd.concat([df[["gene_symbol", "go_id"]] for df in other["frames"].values()],
                                    ignore_index=True).drop_duplicates()
            transferred = transfer_ortholog_annotations(other_union, orthologs)
            frames[f"ortholog_{other_species}"] = transferred
            for key, val in other["sources"].items():
                sources[f"ortholog_{other_species}_{key}"] = val
            sources["gene_orthologs"] = {"path": str(orth_target), "url": "local" if gene_orthologs_path else NCBI_GENE_ORTHOLOGS_URL,
                                         "unmapped_pairs": next(iter(unmapped_pairs))}
            ortholog_meta.update({"species": other_species, "protocol": other["protocol"],
                                  "ortholog_pairs": int(sum(len(v) for v in orthologs.values())),
                                  "other_species_rows": int(len(other_union)),
                                  "transferred_rows": int(len(transferred)),
                                  "transferred_genes": int(transferred.attrs.get("transferred_genes", 0))})
            protocol = f"{protocol}+ortholog({other_species}:{other['protocol']})"
        except Exception as exc:
            problems.append(f"ortholog transfer from {other_species} unavailable: {exc}")
            ortholog_meta["error"] = str(exc)
            logger.warning("Ortholog transfer from %s unavailable for %s (%s); building without it", other_species, species, exc)
    elif use_orthologs:
        problems.append("no ortholog_species configured for this species; no ortholog transfer")

    parsed, counts = build_tree_with_annotation_frames(str(obo_target), frames, use_regulates=use_regulates)
    propagation = "is_a,part_of,regulates" if use_regulates else "is_a,part_of"
    logger.info("GO-Elite %s annotation protocol: %s; propagation %s; rows %s", species, protocol, propagation, counts)

    obo_version = _extract_obo_version(obo_target)
    gaf_version = sources.get("gaf", {}).get("release")
    base_version = _sanitize_version(user_version or gaf_version or obo_version)
    if user_version:
        version = base_version
    else:
        version = base_version + ("_union" if len(own["frames"]) == 2 else "_" + next(iter(own["frames"])))
        version += ("_orth" if "+ortholog(" in protocol else "") + ("_reg" if use_regulates else "")
    version_dir = species_dir / version
    _set_aside(version_dir)
    shutil.move(str(working_dir), str(version_dir))
    for key in sources:
        sources[key]["path"] = sources[key]["path"].replace(str(working_dir), str(version_dir))

    # synonyms (GAF synonyms and gene_info synonyms) for mapping outdated query symbols
    synonyms: Dict[str, str] = {}
    if "gaf" in sources:
        synonyms.update(parse_gaf_synonyms(sources["gaf"]["path"]))
    for k, v in gene_info_synonyms.items():
        synonyms.setdefault(k, v)

    data_dir = version_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    tree_path = data_dir / "go_tree.json"
    term_gene_base = data_dir / "term_gene"

    _write_tree_json(parsed.tree, tree_path)
    term_table_path, storage_format = _write_term_gene_table(parsed.term_to_genes, term_gene_base)
    pd.DataFrame({"synonym": list(synonyms), "symbol": list(synonyms.values())}).to_csv(
        data_dir / "synonyms.tsv.gz", sep="\t", index=False)
    sources["obo"]["release"] = obo_version

    metadata = {
        "species": species,
        "version": version,
        "built_at": datetime.utcnow().isoformat(),
        "annotation_protocol": protocol,
        "propagation": propagation,
        "ortholog_transfer": ortholog_meta,
        "annotation_rows": counts,
        "build_problems": problems,
        "source": sources,
        "storage": {
            "tree": str(tree_path.relative_to(version_dir)),
            "term_gene": str(term_table_path.relative_to(version_dir)),
            "term_gene_format": storage_format,
        },
        "counts": {
            "terms": len(parsed.tree.nodes),
            "annotations": sum(len(genes) for genes in parsed.term_to_genes.values()),
            "propagated": True,
            "synonyms": len(synonyms),
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

    parsed.version_dir = str(version_dir)  # type: ignore[attr-defined]
    return parsed


def _collect_species_annotations(species: str, cfg: Dict[str, str], *, temp_dir: Path, cache_root: Path,
                                 gaf_path: Optional[str], gaf_target: Path, gene_info_path: Optional[str],
                                 gene2go_path: Optional[str], acceptable_evidence: Optional[set[str]],
                                 use_gene2go: bool, logger: logging.Logger) -> Dict[str, object]:
    """The union protocol for one species: UniProt-GOA GAF + NCBI gene2go (through gene_info). Each
    source that cannot be obtained is recorded in ``problems`` and skipped; the caller decides whether
    an empty result is fatal."""
    from .parser import parse_gaf, parse_gene_info, parse_gene2go
    frames: Dict[str, pd.DataFrame] = {}
    sources: Dict[str, Dict[str, str]] = {}
    problems: list = []
    gene_info_synonyms: Dict[str, str] = {}
    gene_id_to_symbol: Dict[str, str] = {}
    gaf_source_url = None if gaf_path else cfg["gaf_url"]
    # Source 1: UniProt-GOA GAF
    try:
        _copy_or_download(gaf_path, gaf_source_url, gaf_target)
        frames["goa_gaf"] = parse_gaf(str(gaf_target), evidence=acceptable_evidence)
        sources["gaf"] = {"path": str(gaf_target), "url": gaf_source_url or "local", "release": _extract_gaf_version(gaf_target)}
    except Exception as exc:
        problems.append(f"GOA GAF unavailable: {exc}")
        logger.warning("GOA GAF unavailable for %s (%s); building from gene2go only", species, exc)
    # Source 2: NCBI gene2go + gene_info
    if use_gene2go and cfg.get("taxon_id") and cfg.get("gene_info_url"):
        try:
            gi_target = temp_dir / (Path(gene_info_path).name if gene_info_path else cfg["gene_info_filename"])
            _copy_or_download(gene_info_path, None if gene_info_path else cfg["gene_info_url"], gi_target)
            # gene2go covers every taxon (about 1.4 GB); it is kept once at the cache root and reused.
            if gene2go_path:
                g2g_target = Path(gene2go_path)
            else:
                g2g_target = cache_root / "downloads" / "gene2go.gz"
                if not g2g_target.exists():
                    _download(NCBI_GENE2GO_URL, g2g_target)
            gene_id_to_symbol, gene_info_synonyms = parse_gene_info(str(gi_target), int(cfg["taxon_id"]))
            frames["ncbi_gene2go"] = parse_gene2go(str(g2g_target), gene_id_to_symbol, int(cfg["taxon_id"]))
            sources["gene2go"] = {"path": str(g2g_target), "url": "local" if gene2go_path else NCBI_GENE2GO_URL,
                                  "taxon_id": cfg["taxon_id"],
                                  "unmapped_gene_ids": str(frames["ncbi_gene2go"].attrs.get("unmapped_gene_ids", 0))}
            sources["gene_info"] = {"path": str(gi_target), "url": "local" if gene_info_path else cfg["gene_info_url"]}
        except Exception as exc:
            problems.append(f"NCBI gene2go/gene_info unavailable: {exc}")
            logger.warning("NCBI gene2go unavailable for %s (%s); building from the GOA GAF only", species, exc)
    elif use_gene2go:
        problems.append("no taxon_id/gene_info_url configured for this species; GAF only")
    protocol = "union(goa_gaf,ncbi_gene2go)" if len(frames) == 2 else (f"single({next(iter(frames))})" if frames else "none")
    return {"frames": frames, "sources": sources, "problems": problems, "gene_info_synonyms": gene_info_synonyms,
            "gene_id_to_symbol": gene_id_to_symbol, "protocol": protocol}


def load_symbol_synonyms(version_dir: str) -> Dict[str, str]:
    """Synonym -> current symbol map for a cached species version, built once from the GAF in
    ``<version_dir>/downloads`` (or the species-level ``downloads``) and cached as data/synonyms.tsv.gz."""
    vdir = Path(version_dir)
    species_name = vdir.parent.name
    # The GAF may sit beside this version, beside the species, or in another cache root's species folder
    # (the repository bundles a cache without its downloads; the user cache keeps them).
    search = [vdir / "downloads", vdir.parent / "downloads"]
    for root in candidate_cache_dirs(None):
        search += [root / species_name / vdir.name / "downloads", root / species_name / "downloads"]
    for d in search:
        cache_path = d.parent / vdir.name / "data" / "synonyms.tsv.gz" if d.name == "downloads" and d.parent != vdir else vdir / "data" / "synonyms.tsv.gz"
        if cache_path.exists():
            df = pd.read_csv(cache_path, sep="\t", dtype=str)
            return dict(zip(df["synonym"], df["symbol"]))
    candidates = [g for d in search if d.exists() for g in sorted(d.glob("*.gaf*"))]
    if not candidates:
        return {}
    from .parser import parse_gaf_synonyms
    mapping = parse_gaf_synonyms(str(candidates[0]))
    out_dir = candidates[0].parent.parent / vdir.name / "data"
    if not out_dir.exists():
        out_dir = candidates[0].parent.parent / "data"
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"synonym": list(mapping), "symbol": list(mapping.values())}).to_csv(out_dir / "synonyms.tsv.gz", sep="\t", index=False)
    except OSError:
        pass
    return mapping
