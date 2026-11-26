#!/usr/bin/env python3
"""Combine allelic_tests outputs written per chromosome or per region chunk, and recompute Benjamini-Hochberg
over the combined family.

allelic_tests computes BH inside one run (one chromosome or one region). A family spans every chunk, so
this step concatenates the chunk tables, refuses duplicated tests (overlapping regions), and adds
``BH_q_combined``: within ``strata`` for the four-family tables, over the whole table otherwise.

Run: python -m altanalyze3.components.snp.combine_tests --stats STATS_DIR --out STATS_DIR/combined
Outputs: <out>/<table>.tsv.gz for every table found, <out>/logs/combine_tests.parameters.json
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import sys

import pandas as pd

try:
    from altanalyze3.components.snp.allelic_tests import bh
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from allelic_tests import bh

# table -> (p-value column, key columns that identify one test, grouping column for BH or None)
TABLES = {
    "isoform_tests": ("CMH_p", ["site_id", "gene", "final", "strata"], "strata"),
    "isoform_state_tests": ("CMH_p", ["site_id", "gene", "final", "cell_state"], None),
    "isoform_covariate_tests": ("CMH_p", ["site_id", "gene", "final", "covariate"], None),
    "isoform_state_heterogeneity": ("Q_p", ["key"], None),
    "isoform_covariate_heterogeneity": ("Q_p", ["key"], None),
    "junction_tests": ("CMH_p", ["site_id", "gene", "junction", "strata"], "strata"),
    "junction_state_tests": ("CMH_p", ["site_id", "gene", "junction", "cell_state"], None),
    "junction_covariate_tests": ("CMH_p", ["site_id", "gene", "junction", "covariate"], None),
    "ase_site_donor": ("binom_p", ["site_id", "donor"], None),
    "ase_site_donor_state": ("binom_p", ["site_id", "donor", "cell_state"], None),
}
for _t in list(TABLES):
    if not _t.startswith(("ase", "junction")):
        TABLES[_t + "_haplotype"] = TABLES[_t]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def combine(stats_dir, out_dir):
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
    chunks = sorted(d for d in glob.glob(os.path.join(stats_dir, "chr*")) if os.path.isdir(d))
    log = dict(module="altanalyze3.components.snp.combine_tests", stats=stats_dir, chunks=[os.path.basename(c) for c in chunks],
               tables={})
    for table, (pcol, key, group) in TABLES.items():
        parts = []
        for c in chunks:
            p = os.path.join(c, table + ".tsv.gz")
            if os.path.exists(p):
                df = pd.read_csv(p, sep="\t", dtype={"site_id": str})
                df.insert(0, "chunk", os.path.basename(c))
                parts.append(df)
        if not parts:
            continue
        df = pd.concat(parts, ignore_index=True)
        dup = df.duplicated(key, keep=False)
        if dup.any():
            raise SystemExit("{}: {} rows test the same {} in more than one chunk (overlapping regions?)".format(
                table, int(dup.sum()), key))
        df["BH_q_combined"] = float("nan")
        if group:
            for _, g in df.groupby(group):
                df.loc[g.index, "BH_q_combined"] = bh(g[pcol].values)
        else:
            df["BH_q_combined"] = bh(df[pcol].values)
        path = os.path.join(out_dir, table + ".tsv.gz")
        df.to_csv(path + ".tmp", sep="\t", index=False, compression="gzip")
        os.replace(path + ".tmp", path)
        log["tables"][table] = dict(chunks=len(parts), rows=int(len(df)), with_p=int(df[pcol].notna().sum()),
                                    q_below_0_05=int((df.BH_q_combined < 0.05).sum()))
    log.update(finished=datetime.datetime.now().isoformat(timespec="seconds"),
               module_sha256=_sha256(os.path.abspath(__file__)))
    with open(os.path.join(out_dir, "logs", "combine_tests.parameters.json"), "w") as fh:
        json.dump(log, fh, indent=2)
    return log


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m altanalyze3.components.snp.combine_tests",
                                 description="Concatenate allelic_tests chunks and recompute BH per family.")
    ap.add_argument("--stats", required=True, help="allelic_tests --out folder holding chr* chunk folders")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    log = combine(a.stats, a.out)
    print("[combine_tests] {} chunks: {}".format(len(log["chunks"]),
                                                {k: v["rows"] for k, v in log["tables"].items()}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
