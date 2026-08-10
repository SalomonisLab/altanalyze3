"""Differential splicing for bulk RNA-Seq, on top of the long-read statistics engine.

The long-read path in ``components/long_read/comparisons.py`` already computes two-group
differentials per cell state and writes ``<cond1>-<cond2>-<cluster>_stats.txt``. It selects the
samples of a state by an index prefix, ``state.sample``, which bulk data does not carry: a bulk
sample is its own unit.

This module supplies that missing adapter and nothing else. The statistics come from
``run_metadataAnalysis``, so a bulk run and a long-read run share one engine, one column schema and
one file name, which is what makes the two comparable.
"""

import logging
import os

import pandas as pd

from ..long_read.comparisons import run_metadataAnalysis

### Name of the pseudo-state a bulk cohort forms. The long-read files carry a real state here.
BULK_CLUSTER = "bulk"


def load_psi_matrix(psi_file):
    """Read a PSI file written by ``altanalyze3 psi`` into events by samples.

    The writer emits a header of sample names with no label over the event column, so the event
    identifier arrives as the frame's first column.
    """
    frame = pd.read_csv(psi_file, sep="\t", index_col=0)
    ### The long-read stats files name this column 'Feature'. Matching it keeps a bulk stats file
    ### and a long-read stats file directly comparable, which is the point of running both.
    frame.index.name = "Feature"
    return frame


def load_groups(groups_file, sample_column=None, group_column=None):
    """Read a two-column groups file into the frame ``run_metadataAnalysis`` expects.

    Accepts a headerless ``sample<TAB>group`` file, which is the common form, or a file whose
    columns are named. The result is indexed by sample and carries one ``grp`` column.
    """
    frame = pd.read_csv(groups_file, sep="\t", header=None, dtype=str)
    if str(frame.iloc[0, 0]).lower() in ("sample", "sample_id", "uid", "library"):
        frame = pd.read_csv(groups_file, sep="\t", dtype=str)
        sample_column = sample_column or frame.columns[0]
        group_column = group_column or frame.columns[1]
    else:
        frame.columns = ["sample", "grp"] + list(frame.columns[2:])
        sample_column, group_column = "sample", "grp"
    frame = frame[[sample_column, group_column]].copy()
    frame.columns = ["sample", "grp"]
    frame = frame.set_index("sample")
    frame["grp"] = frame["grp"].astype(str)
    return frame


def run_bulk_differential(psi_file, groups_file, condition1, condition2, outdir,
                          method="limma", cluster=BULK_CLUSTER):
    """Compare two groups of bulk samples and write the long-read stats file.

    ``method`` is 'limma' for the empirical-Bayes moderated t-test, or 'mwu' for the Mann-Whitney
    rank test. Both write the same columns, so downstream annotation does not change.

    Returns the path written, or None when the comparison could not run.
    """
    psi_matrix = load_psi_matrix(psi_file) if isinstance(psi_file, str) else psi_file
    groups = load_groups(groups_file) if isinstance(groups_file, str) else groups_file

    selected = groups[groups["grp"].isin([condition1, condition2])]
    present = selected.index.intersection(psi_matrix.columns)
    logging.info(
        f"{len(groups)} samples in the groups file; {len(selected)} carry '{condition1}' or "
        f"'{condition2}'; {len(present)} of those appear in the PSI matrix"
    )
    ### Say which samples were dropped rather than letting them disappear into a smaller N.
    missing = sorted(set(selected.index) - set(psi_matrix.columns))
    if missing:
        logging.warning(f"{len(missing)} grouped samples are absent from the PSI matrix: {missing}")
    unassigned = sorted(set(psi_matrix.columns) - set(groups.index))
    if unassigned:
        logging.warning(f"{len(unassigned)} PSI columns have no group: {unassigned}")

    counts = selected.loc[present, "grp"].value_counts().to_dict()
    logging.info(f"group sizes: {counts}")
    if len(counts) < 2:
        logging.error(f"Need both groups; found {counts}")
        return None
    if min(counts.values()) < 2:
        logging.error(f"Each group needs at least 2 samples; found {counts}")
        return None

    os.makedirs(outdir, exist_ok=True)
    run_metadataAnalysis(
        cluster, psi_matrix, condition1, condition2, selected.loc[present],
        outdir, method=method,
    )
    stats_file = os.path.join(outdir, f"{condition1}-{condition2}-{cluster}_stats.txt")
    if not os.path.exists(stats_file):
        logging.error(f"Expected stats file was not written: {stats_file}")
        return None
    rows = sum(1 for _ in open(stats_file)) - 1
    logging.info(f"Wrote {rows} tested events to {stats_file}")
    return stats_file


def run_bulk_differential_cli(args):
    """CLI entry point for ``python -m altanalyze3 diff-splice``."""
    stats_file = run_bulk_differential(
        str(args.psi), str(args.groups), args.condition1, args.condition2,
        str(args.output), method=args.method,
    )
    if stats_file is None:
        raise RuntimeError(
            "No differential stats were produced. Check that both conditions name groups in the "
            "groups file and that each has at least 2 samples present in the PSI matrix."
        )
