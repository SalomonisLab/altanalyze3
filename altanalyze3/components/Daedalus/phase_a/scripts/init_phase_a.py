#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    for relative in [
        "data/raw",
        "data/interim",
        "data/processed",
        "benchmarks",
        "logs"
    ]:
        (root / relative).mkdir(parents=True, exist_ok=True)
    print(f"[ok] Initialized Phase A workspace under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
