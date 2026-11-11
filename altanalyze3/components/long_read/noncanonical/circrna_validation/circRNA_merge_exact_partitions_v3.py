"""Admit complete pinned exact-model partitions without changing read alignment science.

The output is a new, explicitly partitioned event directory. Original serial and
partition files remain immutable; per-model payloads are linked, not rewritten.
"""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def query_id(library, query):
    return "Q" + hashlib.sha256(
        json.dumps([library, query], separators=(",", ":")).encode()
    ).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def main(args):
    root = Path(args.partitions)
    model_dir = Path(args.model_dir)
    source = Path(args.source_shard)
    admission = Path(args.admission_manifest)
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(output)
    stage_path = model_dir / "EVENT_STAGE_RESULT.json"
    bundle_path = model_dir / "PREPARED_MODEL_BUNDLE.json"
    minimap_path = model_dir / "alignments/ALIGNMENT_EXECUTION.json"
    stage, bundle, admitted = map(read_json, (stage_path, bundle_path, admission))
    if admitted["status"] != "INDEPENDENT_EVENT_MODEL_ARTIFACT_ADMISSION":
        raise ValueError("Missing independent model admission")
    for relative, digest in admitted["model_artifacts"].items():
        if sha(model_dir / relative) != digest:
            raise ValueError("Model admission artifact changed: " + relative)
    source_hash = sha(source)
    if source_hash != admitted["source_shard_sha256"] or source_hash != stage["input_sha256"]:
        raise ValueError("Source shard changed")
    models = bundle["templates"]
    groups = collections.OrderedDict()
    for model in models:
        seqhash = hashlib.sha256(model["sequence"].encode()).hexdigest()
        if seqhash != model["sequence_sha256"]:
            raise ValueError("Prepared model sequence changed")
        groups.setdefault((seqhash, tuple(model.get("boundaries", []))), []).append(model)
    if not 0 <= args.serial_completed_through <= len(groups):
        raise ValueError("Invalid claimed completed serial group bound")
    expected_queries = []
    pair_rows = 0
    with source.open() as handle:
        for line in handle:
            row = json.loads(line)
            qid = query_id(row["library"], row["query"])
            if qid != row["alignment_query_id"]:
                raise ValueError("Source query identity changed")
            expected_queries.append((qid, row["library"], row["query"]))
            pair_rows += len(row["baseline_pairs"])
    ids = [row[0] for row in expected_queries]
    ids_set = set(ids)
    if len(ids) != len(set(ids)) or len(ids) != stage["source_contexts"]:
        raise ValueError("Duplicate or missing source queries")
    if pair_rows != stage["raw_pair_rows"]:
        raise ValueError("Raw pair count changed")
    pins = dict(
        independent_admission=sha(admission), source_shard=source_hash,
        event_stage=sha(stage_path), prepared_model_bundle=sha(bundle_path),
        minimap_execution=sha(minimap_path),
    )
    all_parts = sorted(root.glob("*/GROUP_PART_COMPLETE.json"))
    if not all_parts:
        raise ValueError("No completed partitions")
    found = {}
    part_provenance = []
    for manifest_path in all_parts:
        manifest = read_json(manifest_path)
        if manifest["status"] != "PARTITION_ALL_SELECTED_MODELS_ALL_SOURCE_QUERIES_EXACT":
            raise ValueError("Incomplete partition: " + str(manifest_path))
        if manifest["event"] != stage["event"] or manifest["source_contexts"] != len(ids):
            raise ValueError("Partition event/query mismatch")
        if manifest["raw_pair_rows"] != pair_rows or manifest["all_model_count"] != len(models):
            raise ValueError("Partition pair/model count mismatch")
        if manifest["distinct_sequence_boundary_groups"] != len(groups):
            raise ValueError("Partition group universe mismatch")
        for key, value in pins.items():
            if manifest["source_pins"][key] != value:
                raise ValueError("Partition pin mismatch: " + key)
        if manifest["backend_module_sha256"] != sha(Path(args.backend_module)):
            raise ValueError("Alignment backend changed")
        if manifest["serial_dispatch_module_sha256"] != sha(Path(args.serial_module)):
            raise ValueError("Serial model evaluator changed")
        expected_range = range(manifest["requested_group_start"], manifest["requested_group_end"] + 1)
        listed = [r["group_index"] for r in manifest["models"]]
        if listed != list(expected_range):
            raise ValueError("Partition range/model rows mismatch")
        for row in manifest["models"]:
            index = row["group_index"]
            if index in found:
                raise ValueError("Overlapping partition group " + str(index))
            aliases = list(groups.values())[index - 1]
            if row["canonical_model_id"] != aliases[0]["model_id"]:
                raise ValueError("Canonical model mismatch")
            if row["all_model_aliases"] != [m["model_id"] for m in aliases]:
                raise ValueError("Model alias mismatch")
            if row["sequence_sha256"] != aliases[0]["sequence_sha256"]:
                raise ValueError("Model sequence mismatch")
            if row["boundaries"] != aliases[0].get("boundaries", []):
                raise ValueError("Model boundaries mismatch")
            if row["query_outcomes"] != len(ids):
                raise ValueError("Per-model query count mismatch")
            path = manifest_path.parent / row["path"]
            if sha(path) != row["sha256"]:
                raise ValueError("Per-model payload SHA mismatch")
            seen = set()
            statuses = collections.Counter()
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    outcome = json.loads(line)
                    qid = outcome["alignment_query_id"]
                    if qid not in ids_set or qid in seen:
                        raise ValueError("Unknown or duplicated query in " + str(path))
                    if outcome["canonical_model_id"] != row["canonical_model_id"]:
                        raise ValueError("Per-query model mismatch")
                    seen.add(qid)
                    statuses[outcome["evidence"]["status"]] += 1
            if seen != ids_set or dict(statuses) != row["status_counts"]:
                raise ValueError("Per-model query/status conservation failed")
            found[index] = (row, path)
        part_provenance.append(dict(path=str(manifest_path), sha256=sha(manifest_path),
                                    start=manifest["requested_group_start"],
                                    end=manifest["requested_group_end"]))
    serial_dir = Path(args.serial_dir)
    serial_reused = []
    for index in range(1, args.serial_completed_through + 1):
        if index in found:
            continue
        aliases = list(groups.values())[index - 1]
        model = aliases[0]
        path = serial_dir / ("model_%05d.outcomes.jsonl.gz" % index)
        if not path.is_file():
            raise ValueError("Missing claimed completed serial group " + str(index))
        seen = set()
        statuses = collections.Counter()
        with gzip.open(path, "rt") as handle:
            for line in handle:
                outcome = json.loads(line)
                qid = outcome["alignment_query_id"]
                if qid not in ids_set or qid in seen:
                    raise ValueError("Unknown or duplicated serial query " + str(index))
                if outcome["canonical_model_id"] != model["model_id"]:
                    raise ValueError("Serial model identity mismatch " + str(index))
                seen.add(qid)
                statuses[outcome["evidence"]["status"]] += 1
        if seen != ids_set:
            raise ValueError("Incomplete claimed serial group " + str(index))
        row = dict(group_index=index, canonical_model_id=model["model_id"],
                   all_model_aliases=[m["model_id"] for m in aliases],
                   sequence_sha256=model["sequence_sha256"],
                   boundaries=model.get("boundaries", []),
                   query_outcomes=len(ids), status_counts=dict(statuses),
                   path=path.name, sha256=sha(path),
                   execution_origin="VERIFIED_PRESERVED_SERIAL_PARTIAL")
        found[index] = (row, path)
        serial_reused.append(dict(index=index, path=str(path), sha256=row["sha256"],
                                  queries=len(ids)))
    if set(found) != set(range(1, len(groups) + 1)):
        missing = sorted(set(range(1, len(groups) + 1)) - set(found))
        raise ValueError("Not all model groups complete; first missing=" + str(missing[:10]))
    serial_parity = []
    for index in map(int, args.serial_parity_groups.split(",")):
        if index not in found:
            raise ValueError("Parity group absent")
        original = serial_dir / ("model_%05d.outcomes.jsonl.gz" % index)
        if not original.is_file():
            raise ValueError("Original serial parity group absent")
        with gzip.open(original, "rt") as a, gzip.open(found[index][1], "rt") as b:
            old = [json.loads(line) for line in a]
            new = [json.loads(line) for line in b]
        if old != new:
            raise ValueError("Serial/partition decoded parity mismatch")
        serial_parity.append(dict(index=index, queries=len(old), serial_sha256=sha(original),
                                  partition_sha256=sha(found[index][1])))
    output.mkdir()
    source_keys = [dict(id=qid, library=lib, query=query) for qid, lib, query in expected_queries]
    query_path = output / "QUERY_IDENTITIES.json"
    query_path.write_text(json.dumps(source_keys, indent=2) + "\n")
    records = []
    for index in sorted(found):
        row, path = found[index]
        name = "model_%05d.outcomes.jsonl.gz" % index
        (output / name).symlink_to(path.resolve())
        records.append({key: value for key, value in row.items() if key != "group_index"})
    report = dict(
        status="ALL_EVENT_QUERIES_ALL_PREPARED_MODELS_EXACT_SEMIGLOBAL_EVIDENCE_RETAINED",
        execution_mode="PINNED_PARTITIONED_EXACT_WITH_INDEPENDENT_MERGE_AND_SERIAL_PILOT_PARITY",
        event=stage["event"], source_contexts=len(ids), raw_pair_rows=pair_rows,
        all_model_count=len(models), distinct_sequence_boundary_groups=len(groups),
        models=records, source_pins={**pins, "query_identities": sha(query_path)},
        code_sha256=sha(Path(__file__)),
        backend_module_sha256=sha(Path(args.backend_module)),
        serial_dispatch_module_sha256=sha(Path(args.serial_module)),
        partition_manifests=part_provenance, serial_reused=serial_reused,
        serial_completed_through=args.serial_completed_through,
        serial_parity=serial_parity,
        limits=["This is sequence fit, not proof of physical circularity.",
                "All source queries and all prepared model groups are retained.",
                "Model results were generated in separate partitions and linked without modification.",
                "Preserved serial groups are reused only through the explicitly claimed completed index after full gzip, query-roster and model-identity validation.",
                "Serial parity was checked for listed groups; remaining groups conserve exact source/model pins and query rosters."],
    )
    (output / "SEMIGLOBAL_EVENT_COMPLETE.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(status=report["status"], execution_mode=report["execution_mode"],
                          event=report["event"], queries=len(ids), groups=len(groups),
                          partitions=len(part_provenance), serial_parity=serial_parity)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("partitions", "model_dir", "source_shard", "admission_manifest",
                "output", "backend_module", "serial_module", "serial_dir"):
        parser.add_argument("--" + key.replace("_", "-"), required=True)
    parser.add_argument("--serial-parity-groups", default="1")
    parser.add_argument("--serial-completed-through", type=int, default=0)
    options = parser.parse_args()
    main(options)
