"""Validate community metadata without downloading or executing model payloads."""
import argparse
import json
import re
from datetime import date
from pathlib import Path
from .registry import _identity


def validate_proposal(record):
    required = ("name", "component", "submitter", "license", "species", "tissue",
                "trained_at", "method", "training_sources", "input_roster", "output_roster",
                "sample_roster", "evaluation", "artifacts", "baseline_comparison")
    for key in required:
        if key not in record or record[key] in (None, "", [], {}):
            raise ValueError(f"Missing required field: {key}")
    for key in ("artifacts", "method", "evaluation", "baseline_comparison"):
        if not isinstance(record[key], dict):
            raise ValueError(f"{key} must be an object")
    if not isinstance(record["training_sources"], list):
        raise ValueError("training_sources must list the original sources")
    try:
        date.fromisoformat(record["trained_at"])
    except (TypeError, ValueError) as exc:
        raise ValueError("trained_at must be an ISO date") from exc
    if record.get("status") != "proposed":
        raise ValueError("Community submissions must have status=proposed")
    if record.get("default_history") or record.get("approval"):
        raise ValueError("A proposal cannot approve itself or declare default history")
    for role, artifact in record["artifacts"].items():
        if not re.fullmatch(r"[a-f0-9]{64}", artifact.get("sha256", "")):
            raise ValueError(f"Invalid SHA-256 for artifact {role}")
        if not str(artifact.get("url", "")).startswith("https://"):
            raise ValueError(f"Artifact {role} requires an HTTPS download URL")
    expected = record["component"] + ":sha256:" + _identity({
        role: artifact["sha256"] for role, artifact in record["artifacts"].items()})
    if record.get("model_version_id") != expected:
        raise ValueError(f"model_version_id must be {expected}")
    for roster in ("input_roster", "output_roster", "sample_roster"):
        value = record[roster]
        if not isinstance(value, dict) or not re.fullmatch(r"[a-f0-9]{64}", value.get("sha256", "")):
            raise ValueError(f"{roster} requires the hash of a complete ordered identifier file")
        if not isinstance(value.get("count"), int) or isinstance(value["count"], bool) or value["count"] < 1 or not value.get("url"):
            raise ValueError(f"{roster} requires a positive count and source URL")
    method = record["method"]
    for key in ("estimator", "preprocessing", "feature_selection", "hyperparameters",
                "target_transform", "missing_input_handling", "code_url", "code_commit"):
        if key not in method:
            raise ValueError(f"method.{key} is required (use explicit 'none' when appropriate)")
    comparison = record["baseline_comparison"]
    for key in ("model_version_id", "coverage_report_url", "method_differences", "evaluation_report_url"):
        if key not in comparison:
            raise ValueError(f"baseline_comparison.{key} is required")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.paths:
        validate_proposal(json.loads(path.read_text()))
        print(f"Valid proposal: {path}")


if __name__ == "__main__":
    main()
