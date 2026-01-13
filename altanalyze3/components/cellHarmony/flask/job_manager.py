from __future__ import annotations

import json
import threading
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
from uuid import uuid4


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

    # ------------------------------------------------------------------- public
    def create_job(
        self,
        species: str,
        reference: str,
        soupx_option: Optional[str],
        files: List[Dict[str, str]],
    ) -> Dict:
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
            "progress": 5,
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
