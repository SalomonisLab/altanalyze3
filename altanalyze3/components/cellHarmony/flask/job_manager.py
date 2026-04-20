from __future__ import annotations

import json
import shutil
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, List, Optional
from uuid import uuid4


TERMINAL_JOB_STATUSES = {"completed", "failed", "cancelled", "canceled"}
ACTIVE_DIFFERENTIAL_STATUSES = {"queued", "processing"}


class JobStore:
    """
    Lightweight filesystem-backed job registry used by the Flask application.

    Each job lives in ``<root>/<job_id>`` and contains ``uploads/``, ``outputs/``
    and ``logs/`` directories alongside a ``job.json`` metadata file.
    """

    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    # ------------------------------------------------------------------ helpers
    def _job_dir(self, job_id: str) -> Path:
        return self.root / job_id

    def _metadata_path(self, job_id: str) -> Path:
        return self._job_dir(job_id) / "job.json"

    def _read_metadata(self, job_id: str) -> Dict:
        with self._metadata_path(job_id).open("r", encoding="utf-8") as handle:
            return json.load(handle)

    def _write_metadata(self, job_id: str, data: Dict) -> None:
        tmp_path = self._metadata_path(job_id)
        with tmp_path.open("w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, sort_keys=True)

    def _parse_timestamp(self, raw: object) -> Optional[datetime]:
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

    def _is_purge_eligible(self, metadata: Dict, cutoff: datetime) -> bool:
        status = str(metadata.get("status") or "").strip().lower()
        if status not in TERMINAL_JOB_STATUSES:
            return False
        differential = metadata.get("differential") or {}
        differential_status = str(differential.get("status") or "").strip().lower()
        if differential_status in ACTIVE_DIFFERENTIAL_STATUSES:
            return False
        updated_at = self._parse_timestamp(metadata.get("updated_at")) or self._parse_timestamp(metadata.get("created_at"))
        if updated_at is None:
            return False
        return updated_at < cutoff

    def _purge_old_jobs_unlocked(self, *, max_age_hours: float) -> int:
        cutoff = datetime.now(timezone.utc) - timedelta(hours=max_age_hours)
        deleted = 0
        for job_dir in sorted(self.root.iterdir()):
            if not job_dir.is_dir():
                continue
            meta_path = job_dir / "job.json"
            if not meta_path.exists():
                continue
            try:
                metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            if not self._is_purge_eligible(metadata, cutoff):
                continue
            shutil.rmtree(job_dir, ignore_errors=False)
            deleted += 1
        return deleted

    # ------------------------------------------------------------------- public
    def purge_old_jobs(self, *, max_age_hours: float = 8.0) -> int:
        with self._lock:
            return self._purge_old_jobs_unlocked(max_age_hours=max_age_hours)

    def create_job(
        self,
        species: str,
        reference: str,
        soupx_option: Optional[str],
        files: List[Dict[str, str]],
    ) -> Dict:
        with self._lock:
            self._purge_old_jobs_unlocked(max_age_hours=8.0)
            job_id = uuid4().hex
            job_dir = self._job_dir(job_id)
            job_dir.mkdir(parents=True, exist_ok=True)
            for sub in ("uploads", "outputs", "logs"):
                (job_dir / sub).mkdir(parents=True, exist_ok=True)

            metadata = {
                "job_id": job_id,
                "species": species,
                "reference": reference,
                "soupx_option": soupx_option,
                "files": files,
                "status": "uploaded",
                "message": "Files uploaded, awaiting QC.",
                "progress": 0,
                "created_at": datetime.utcnow().isoformat() + "Z",
                "updated_at": datetime.utcnow().isoformat() + "Z",
                "qc": {},
                "artifacts": {},
            }
            self._write_metadata(job_id, metadata)
            return metadata

    def get_job(self, job_id: str) -> Dict:
        return self._read_metadata(job_id)

    def update_job(self, job_id: str, **changes) -> Dict:
        with self._lock:
            meta = self._read_metadata(job_id)
            meta.update(changes)
            meta["updated_at"] = datetime.utcnow().isoformat() + "Z"
            self._write_metadata(job_id, meta)
            return meta

    def add_artifact(self, job_id: str, key: str, path: Path) -> None:
        with self._lock:
            meta = self._read_metadata(job_id)
            meta.setdefault("artifacts", {})[key] = str(path)
            meta["updated_at"] = datetime.utcnow().isoformat() + "Z"
            self._write_metadata(job_id, meta)

    def append_log(self, job_id: str, message: str) -> None:
        log_path = self._job_dir(job_id) / "logs" / "pipeline.log"
        timestamp = datetime.utcnow().isoformat()
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(f"[{timestamp}] {message}\n")

    def job_exists(self, job_id: str) -> bool:
        return self._metadata_path(job_id).exists()

    def uploads_dir(self, job_id: str) -> Path:
        return self._job_dir(job_id) / "uploads"

    def outputs_dir(self, job_id: str) -> Path:
        return self._job_dir(job_id) / "outputs"

    def logs_dir(self, job_id: str) -> Path:
        return self._job_dir(job_id) / "logs"

    def list_jobs(self) -> List[Dict]:
        jobs: List[Dict] = []
        for job_dir in sorted(self.root.iterdir()):
            meta_path = job_dir / "job.json"
            if not meta_path.exists():
                continue
            jobs.append(self._read_metadata(job_dir.name))
        return jobs
