#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path


def load_manifest(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, destination.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Daedalus Phase A resources.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "resources" / "manifest.json",
        help="Path to the resource manifest JSON."
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "raw",
        help="Destination directory for downloaded files."
    )
    parser.add_argument(
        "--include-large",
        action="store_true",
        help="Also download resources disabled by default."
    )
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    resources = manifest.get("resources", [])
    if not resources:
        print("[error] No resources defined in manifest.", file=sys.stderr)
        return 1

    for resource in resources:
        if not resource.get("enabled_by_default", False) and not args.include_large:
            print(f"[skip] {resource['name']} disabled by default")
            continue
        url = resource["url"]
        target = args.dest / resource["filename"]
        if target.exists():
            print(f"[keep] {resource['name']} already present at {target}")
            continue
        print(f"[download] {resource['name']} -> {target}")
        download(url, target)
    print("[ok] Resource download pass complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
