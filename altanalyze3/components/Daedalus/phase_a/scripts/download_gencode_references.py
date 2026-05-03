#!/usr/bin/env python3
from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path


RESOURCES = [
    {
        "name": "gencode_human_gtf",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
        "filename": "gencode.v49.annotation.gtf.gz"
    },
    {
        "name": "gencode_human_pc_transcripts",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_transcripts.fa.gz",
        "filename": "gencode.v49.pc_transcripts.fa.gz"
    },
    {
        "name": "gencode_human_pc_translations",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_translations.fa.gz",
        "filename": "gencode.v49.pc_translations.fa.gz"
    },
    {
        "name": "gencode_mouse_gtf",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz",
        "filename": "gencode.vM38.annotation.gtf.gz"
    },
    {
        "name": "gencode_mouse_pc_transcripts",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.pc_transcripts.fa.gz",
        "filename": "gencode.vM38.pc_transcripts.fa.gz"
    },
    {
        "name": "gencode_mouse_pc_translations",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.pc_translations.fa.gz",
        "filename": "gencode.vM38.pc_translations.fa.gz"
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
    parser = argparse.ArgumentParser(description="Download GENCODE Phase A references.")
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "raw" / "gencode"
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
    print("[ok] GENCODE download pass complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
