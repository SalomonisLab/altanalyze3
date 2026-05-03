from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, Iterable, List


PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE_REGISTRY = PACKAGE_DIR / "configs" / "source_registry.json"


@dataclass(frozen=True)
class FastCommSource:
    source_id: str
    kind: str
    access: str
    url: str
    citation: str
    priority: str
    notes: str = ""


def load_source_registry(path: Path | str = DEFAULT_SOURCE_REGISTRY) -> Dict[str, object]:
    registry_path = Path(path)
    with registry_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object in {registry_path}")
    return data


def list_sources(path: Path | str = DEFAULT_SOURCE_REGISTRY) -> List[FastCommSource]:
    data = load_source_registry(path)
    sources = data.get("sources", [])
    if not isinstance(sources, list):
        raise ValueError("Source registry 'sources' must be a list")
    parsed: List[FastCommSource] = []
    for source in sources:
        parsed.append(
            FastCommSource(
                source_id=str(source["source_id"]),
                kind=str(source["kind"]),
                access=str(source["access"]),
                url=str(source["url"]),
                citation=str(source.get("citation", "")),
                priority=str(source.get("priority", "")),
                notes=str(source.get("notes", "")),
            )
        )
    return parsed


def artifact_paths(path: Path | str = DEFAULT_SOURCE_REGISTRY) -> Dict[str, Path]:
    data = load_source_registry(path)
    policy = data.get("artifact_policy", {})
    if not isinstance(policy, dict):
        raise ValueError("Source registry 'artifact_policy' must be an object")
    root = PACKAGE_DIR.parents[1]
    resolved: Dict[str, Path] = {}
    for key, value in policy.items():
        if key == "gitignore_required":
            continue
        target = Path(str(value))
        if not target.is_absolute():
            target = root / target
        resolved[key] = target
    return resolved


def ensure_artifact_dirs(keys: Iterable[str] | None = None) -> Dict[str, Path]:
    paths = artifact_paths()
    selected = set(keys) if keys is not None else set(paths)
    created: Dict[str, Path] = {}
    for key, path in paths.items():
        if key in selected:
            path.mkdir(parents=True, exist_ok=True)
            created[key] = path
    return created
