from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable, Optional

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parent))
    from config import load_config  # type: ignore
else:  # pragma: no cover
    from .config import load_config


TERMINAL_JOB_STATUSES = {"completed", "failed", "cancelled", "canceled"}
ACTIVE_DIFFERENTIAL_STATUSES = {"queued", "processing"}


@dataclass(frozen=True)
class JobRecord:
    job_id: str
    job_dir: Path
    metadata_path: Path
    metadata: dict
    status: str
    updated_at: datetime
    differential_status: str


def _parse_args() -> argparse.Namespace:
    config = load_config()
    parser = argparse.ArgumentParser(
        description="Purge old terminal cellHarmony web jobs from the filesystem.",
    )
    parser.add_argument(
        "--job-root",
        default=config["JOB_STORAGE"],
        help="Directory containing per-job folders. Defaults to CELLHARMONY_JOB_STORAGE.",
    )
    parser.add_argument(
        "--retain-days",
        type=int,
        default=7,
        help="Keep jobs updated within the last N days. Default: 7.",
    )
    parser.add_argument(
        "--keep-latest",
        type=int,
        default=0,
        help="Always retain the newest N purge-eligible jobs, regardless of age. Default: 0.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would be deleted without removing any files.",
    )
    return parser.parse_args()


def _parse_timestamp(raw: object) -> Optional[datetime]:
    if not raw:
        return None
    text = str(raw).strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _iter_job_records(job_root: Path) -> Iterable[JobRecord]:
    for job_dir in sorted(job_root.iterdir()):
        if not job_dir.is_dir():
            continue
        metadata_path = job_dir / "job.json"
        if not metadata_path.exists():
            print(f"skip {job_dir.name}: missing job.json")
            continue
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            print(f"skip {job_dir.name}: invalid job.json ({exc})")
            continue

        updated_at = _parse_timestamp(metadata.get("updated_at")) or _parse_timestamp(metadata.get("created_at"))
        if updated_at is None:
            print(f"skip {job_dir.name}: missing created_at/updated_at")
            continue

        status = str(metadata.get("status") or "").strip().lower()
        differential = metadata.get("differential") or {}
        differential_status = str(differential.get("status") or "").strip().lower()
        yield JobRecord(
            job_id=str(metadata.get("job_id") or job_dir.name),
            job_dir=job_dir,
            metadata_path=metadata_path,
            metadata=metadata,
            status=status,
            updated_at=updated_at,
            differential_status=differential_status,
        )


def _retention_reason(record: JobRecord, cutoff: datetime) -> Optional[str]:
    if record.status not in TERMINAL_JOB_STATUSES:
        return f"top-level status '{record.status or 'unknown'}' is not terminal"
    if record.differential_status in ACTIVE_DIFFERENTIAL_STATUSES:
        return f"differential status '{record.differential_status}' is still active"
    if record.updated_at >= cutoff:
        age_days = (datetime.now(timezone.utc) - record.updated_at).days
        return f"retained ({age_days}d old, cutoff {cutoff.date().isoformat()})"
    return None


def _delete_job(record: JobRecord, dry_run: bool) -> None:
    action = "would delete" if dry_run else "delete"
    print(
        f"{action} {record.job_id}: status={record.status}, "
        f"differential={record.differential_status or 'none'}, "
        f"updated_at={record.updated_at.isoformat()}, path={record.job_dir}"
    )
    if not dry_run:
        shutil.rmtree(record.job_dir)


def main() -> int:
    args = _parse_args()
    job_root = Path(args.job_root).expanduser().resolve()
    if not job_root.exists():
        print(f"job root does not exist: {job_root}", file=sys.stderr)
        return 1
    if not job_root.is_dir():
        print(f"job root is not a directory: {job_root}", file=sys.stderr)
        return 1
    if args.retain_days < 0:
        print("--retain-days must be >= 0", file=sys.stderr)
        return 1
    if args.keep_latest < 0:
        print("--keep-latest must be >= 0", file=sys.stderr)
        return 1

    now = datetime.now(timezone.utc)
    cutoff = now - timedelta(days=args.retain_days)
    records = list(_iter_job_records(job_root))
    print(
        f"scanned {len(records)} jobs under {job_root} "
        f"(retain_days={args.retain_days}, keep_latest={args.keep_latest}, dry_run={args.dry_run})"
    )

    purgeable: list[JobRecord] = []
    for record in records:
        reason = _retention_reason(record, cutoff)
        if reason is not None:
            print(f"keep {record.job_id}: {reason}")
            continue
        purgeable.append(record)

    purgeable.sort(key=lambda record: record.updated_at, reverse=True)
    kept_by_latest = purgeable[: args.keep_latest]
    to_delete = purgeable[args.keep_latest :]

    for record in kept_by_latest:
        print(f"keep {record.job_id}: protected by --keep-latest ({args.keep_latest})")

    for record in to_delete:
        _delete_job(record, dry_run=args.dry_run)

    print(
        f"summary: scanned={len(records)} purgeable={len(purgeable)} "
        f"kept_latest={len(kept_by_latest)} deleted={len(to_delete) if not args.dry_run else 0} "
        f"would_delete={len(to_delete) if args.dry_run else 0}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
