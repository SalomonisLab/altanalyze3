#!/usr/bin/env python3
from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path


RESOURCES = [
    {
        "name": "appris_human_principal",
        "url": "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt",
        "filename": "appris_human_GRCh38_principal.txt"
    },
    {
        "name": "appris_human_appris_scores",
        "url": "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.appris.txt",
        "filename": "appris_human_GRCh38_appris.txt"
    },
    {
        "name": "appris_mouse_principal",
        "url": "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/mus_musculus/GRCm39/appris_data.principal.txt",
        "filename": "appris_mouse_GRCm39_principal.txt"
    },
    {
        "name": "appris_mouse_appris_scores",
        "url": "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/mus_musculus/GRCm39/appris_data.appris.txt",
        "filename": "appris_mouse_GRCm39_appris.txt"
    }
]


def _download(url: str, destination: Path, retries: int = 4) -> None:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url) as response, destination.open("wb") as out:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            return
        except Exception as exc:
            last_error = exc
            if destination.exists():
                destination.unlink()
            print(f"[warn] download attempt {attempt}/{retries} failed for {url}: {exc}")
            if attempt < retries:
                time.sleep(2.0)
    assert last_error is not None
    raise last_error


def main() -> int:
    parser = argparse.ArgumentParser(description="Download APPRIS current-release principal and score files.")
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "raw" / "appris"
    )
    args = parser.parse_args()
    args.dest.mkdir(parents=True, exist_ok=True)
    for resource in RESOURCES:
        target = args.dest / resource["filename"]
        if target.exists():
            print(f"[keep] {resource['name']} already present at {target}")
            continue
        print(f"[download] {resource['name']} -> {target}")
        _download(resource["url"], target)
    print("[ok] APPRIS download pass complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
