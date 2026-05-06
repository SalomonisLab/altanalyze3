#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import http.client
import json
import re
import time
import urllib.parse
import urllib.error
import urllib.request
from pathlib import Path


DEFAULT_QUERY = "(reviewed:true) AND (organism_id:9606 OR organism_id:10090)"
DEFAULT_FIELDS = ",".join(
    [
        "accession",
        "gene_names",
        "organism_name",
        "xref_ensembl"
    ]
)


def _extract_next_link(link_header: str | None) -> str | None:
    if not link_header:
        return None
    match = re.search(r"<([^>]+)>;\s*rel=\"next\"", link_header)
    if match:
        return match.group(1)
    return None


def _make_url(query: str, size: int) -> str:
    params = {
        "query": query,
        "format": "json",
        "size": str(size)
    }
    return "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)


def _fetch_json(url: str, retries: int = 5, sleep_seconds: float = 2.0) -> tuple[dict, str | None]:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(request) as response:
                payload = json.loads(response.read().decode("utf-8"))
                next_link = _extract_next_link(response.headers.get("Link"))
            return payload, next_link
        except (urllib.error.URLError, http.client.IncompleteRead, json.JSONDecodeError) as exc:
            last_error = exc
            print(f"[warn] fetch attempt {attempt}/{retries} failed for {url}: {exc}")
            if attempt < retries:
                time.sleep(sleep_seconds)
    assert last_error is not None
    raise last_error


def main() -> int:
    parser = argparse.ArgumentParser(description="Download reviewed UniProtKB entries for Phase A.")
    parser.add_argument("--query", default=DEFAULT_QUERY)
    parser.add_argument("--size", type=int, default=200)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
    )
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    url = _make_url(args.query, args.size)
    page = 0
    total = 0
    with gzip.open(args.out, "wt", encoding="utf-8") as out_handle:
        while url:
            page += 1
            payload, url = _fetch_json(url)
            for record in payload.get("results", []):
                out_handle.write(json.dumps(record, ensure_ascii=True) + "\n")
                total += 1
            print(f"[page {page}] cumulative_entries={total}")
    print(f"[ok] Wrote {args.out} with {total} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
