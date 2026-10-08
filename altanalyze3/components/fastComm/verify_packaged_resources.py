"""Validate the immutable LR deployment resources without loading model dependencies."""
import csv
import hashlib
import json
from pathlib import Path


def verify_resources(root=None):
    root = Path(root) if root else Path(__file__).with_name("resources")
    result = {}
    for species in ("human", "mouse"):
        folder = root / f"cellchat_nichenet_{species}"
        manifest = json.loads((folder / "manifest.json").read_text())
        path = folder / "ligand_receptor.tsv"
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"ligand", "receptor", "source"}.issubset(reader.fieldnames):
                raise RuntimeError(f"Invalid fastComm LR columns: {path}")
            count = sum(1 for _ in reader)
        if digest != manifest["sha256"] or count != manifest["n_merged_interactions"]:
            raise RuntimeError(f"fastComm resource altered or incomplete: {path}")
        result[species] = {"sha256": digest, "interactions": count}
    return result


if __name__ == "__main__":
    print(json.dumps(verify_resources(), sort_keys=True))
