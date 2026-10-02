"""Calculate donor-weighted relative folds from native log2 lipid predictions.

Input: sample-by-lipid CSV and metadata CSV with sample, donor, condition.
Use only after validating that targets are log2 abundance, with all target
standardization inverted. Z-scores and log2(1+c*A) require separate inversion.
No p-values are inferred from model predictions.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict


def log2_mean_exp(values):
    """Stable log2 of the arithmetic mean of positive relative abundances."""
    peak = max(values)
    return peak + math.log2(math.fsum(2 ** (v - peak) for v in values) / len(values))


def contrast(case_by_donor, control_by_donor):
    # Each list contains native log2 predictions for one donor and cell state.
    case_log = [statistics.mean(v) for v in case_by_donor.values() if v]
    ctrl_log = [statistics.mean(v) for v in control_by_donor.values() if v]
    if not case_log or not ctrl_log:
        return None
    geometric = statistics.mean(case_log) - statistics.mean(ctrl_log)
    case_linear_log = [log2_mean_exp(v) for v in case_by_donor.values() if v]
    ctrl_linear_log = [log2_mean_exp(v) for v in control_by_donor.values() if v]
    arithmetic = log2_mean_exp(case_linear_log) - log2_mean_exp(ctrl_linear_log)
    return geometric, arithmetic, len(case_log), len(ctrl_log)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--control", required=True)
    parser.add_argument("--state", help="Restrict metadata cell_type to this value.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.case == args.control:
        parser.error("Case and control must differ.")
    with open(args.metadata, newline="") as handle:
        records = list(csv.DictReader(handle))
    metadata = {}
    for row in records:
        if not all(row.get(k) for k in ["sample", "donor", "condition"]):
            raise ValueError("Metadata requires sample, donor and condition.")
        if row["sample"] in metadata:
            raise ValueError("Duplicate sample in metadata.")
        metadata[row["sample"]] = row
    with open(args.predictions, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or reader.fieldnames[0] != "sample":
            raise ValueError("First predictions column must be sample.")
        lipids = reader.fieldnames[1:]
        if len(set(lipids)) != len(lipids):
            raise ValueError("Duplicate lipid columns.")
        groups = {arm: {lipid: defaultdict(list) for lipid in lipids}
                  for arm in [args.case, args.control]}
        seen, states = set(), set()
        for row in reader:
            sample = row["sample"]
            if sample in seen:
                raise ValueError("Duplicate prediction sample.")
            seen.add(sample)
            if sample not in metadata:
                raise ValueError("Prediction sample missing from metadata: " + sample)
            info = metadata[sample]
            arm = info["condition"]
            if arm not in groups or (args.state is not None and info.get("cell_type") != args.state):
                continue
            states.add(info.get("cell_type", "unspecified"))
            for lipid in lipids:
                if row[lipid] in ["", "NA", "NaN", "nan"]:
                    continue
                value = float(row[lipid])
                if not math.isfinite(value):
                    raise ValueError("Nonfinite prediction.")
                groups[arm][lipid][info["donor"]].append(value)
    if len(states) > 1:
        raise ValueError("Multiple cell types selected; use --state for a within-state comparison.")
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["lipid", "case", "control", "log2fc_geometric", "log2fc_arithmetic",
                         "n_case_donors", "n_control_donors", "shared_donors"])
        for lipid in lipids:
            result = contrast(groups[args.case][lipid], groups[args.control][lipid])
            if result is not None:
                shared = len(set(groups[args.case][lipid]) & set(groups[args.control][lipid]))
                writer.writerow([lipid, args.case, args.control, *result, shared])


if __name__ == "__main__":
    main()
