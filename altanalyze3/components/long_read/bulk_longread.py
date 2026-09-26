"""Bulk long-read isoform workflow: BAM -> collapsed isoform catalog -> sample-level counts.

This module runs the VALIDATED single-cell long-read collapse on BULK long-read BAMs (PacBio CCS,
ONT). It adds no algorithm. Every structural step stays the one the single-cell workflow uses:

  bam/isoform_structure_extract.parallel_extract_isoform_structures   BAM -> read GFF + molecule h5ad
  gff_process.consolidateLongReadGFFs(mode='collapse')                read GFF -> exon structures
  isoform_collapse.pipeline.run_pipeline                              cross-sample scored collapse
  isoform_collapse.pipeline.rekey_one_sample                          molecules -> final isoforms
  isoform_matrix.pseudo_cluster_counts_optimized                      counts / CPM / ratio matrices

ONE thing differs: the cell dimension. A bulk BAM carries no CB/CR/BC/BX tag, so
``resolve_barcode`` returned None for every read and the extractor discarded the whole file. Bulk
mode gives the extractor a ``default_barcode`` equal to the library name, which makes each BAM ONE
pseudo-cell. Every library then carries the single cluster label ``bulk``, so the existing pseudobulk
step sums all libraries of a ``uid`` into one sample-level column. A ``uid`` with three BAM rows
therefore collapses to one sample, which is the behaviour the single-cell path already had for
multi-BAM samples.

There is no UMI. PacBio CCS and ONT emit one record per source molecule, so the molecule id is the
read name and the count is a molecule count. Nothing is deduplicated, and nothing needs to be.

Junction quantification and PSI are NOT run by this module. Bulk callers ask for isoform structures
and isoform counts; the junction path stays in the ``sclr`` commands.
"""

import os
import glob
import logging
from pathlib import Path

from . import gff_process
from . import io_utils as _io
from . import isoform_automate as isoa
from . import isoform_matrix as iso

#: Single cluster label every bulk library carries. The downstream pseudobulk groups by
#: ``obs['cluster']``, so one shared label makes every library of a uid sum into one column named
#: ``bulk.<uid>``.
BULK_CLUSTER = 'bulk'

#: Name of the synthesized barcode->cluster file bulk mode writes (the same format
#: ``isoform_matrix.import_barcode_clusters`` reads for cellHarmony annotations).
BULK_ANNOTATION_NAME = 'bulk_barcode_clusters.txt'


def bulk_barcode(library):
    """The pseudo-barcode for one bulk library.

    ``import_barcode_clusters`` splits its first column ONCE on '.' into (barcode, sample_name), so a
    library name containing '.' would produce three fields and corrupt the mapping. Reject it here
    with a clear message rather than letting the split fail downstream.
    """
    library = str(library)
    if '.' in library:
        raise ValueError(
            f"Bulk library name '{library}' contains '.', which collides with the "
            f"'<barcode>.<sample>' convention of the barcode->cluster file. Rename the library in "
            f"the metadata (for example replace '.' with '_')."
        )
    return library


def write_bulk_cell_annotation(sample_dict, out_dir, cluster_label=BULK_CLUSTER):
    """Write the one-row-per-library barcode->cluster file bulk mode needs.

    Format matches what ``isoform_matrix.import_barcode_clusters`` / ``return_cluster_order`` read
    from cellHarmony: ``<barcode>.<sample_name><TAB><cluster>``. For bulk the barcode and the sample
    name are both the library, and every library shares one cluster, so the existing cluster-aware
    pseudobulk, filtering and cluster-order code runs unmodified.

    Returns the absolute path written.
    """
    out_dir = str(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, BULK_ANNOTATION_NAME)
    rows = []
    for uid, libs in sample_dict.items():
        for s in libs:
            bc = bulk_barcode(s['library'])
            rows.append(f"{bc}.{bc}\t{cluster_label}")
    # ATOMIC. `blr --sample <uid>` rewrites this whole file from the FULL metadata, so a cluster
    # fan-out runs one writer per sample against one path. Writing in place would let a reader
    # see a truncated file, and two writers interleave. Every writer produces identical bytes, so
    # a temp file plus os.replace (atomic on POSIX within a directory) makes the race harmless.
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, 'w') as handle:
        handle.write("\n".join(rows) + "\n")
    os.replace(tmp, path)
    logging.info("bulk: wrote %d library cluster label(s) -> %s", len(rows), path)
    return path


def annotate_sample_structures(sample_dict, exon_annot, force=False):
    """Per library: read GFF -> ``gff-output/transcript_associations.txt`` (gene, strand, structure,
    molecule, source), the table the collapse consumes.

    In the single-cell workflow this file is a side effect of ``exportJunctionMatrix``. Bulk skips
    junction quantification, so the SAME ``gff_process.consolidateLongReadGFFs`` call is made here
    directly, on ONE GFF at a time, exactly as the single-cell path calls it.

    Returns {library: transcript_associations_path}.
    """
    exon_annot = str(exon_annot)
    produced = {}
    for uid, libs in sample_dict.items():
        for s in libs:
            gff = str(s['gff'])
            ta = os.path.join(os.path.dirname(gff), 'gff-output', 'transcript_associations.txt')
            # Phase 1 gzips this table (io_utils.compress), so the file on disk may be
            # transcript_associations.txt OR transcript_associations.txt.gz. Resolve either form;
            # testing only the uncompressed name re-annotated every library on every re-run.
            existing = _io.resolve(ta, missing_ok=True)
            if existing and os.path.getsize(existing) > 0 and not force:
                logging.info("bulk: structures exist for %s, reusing %s", s['library'], existing)
                produced[s['library']] = existing
                continue
            logging.info("bulk: annotating exon structures for %s from %s", s['library'], gff)
            written = gff_process.consolidateLongReadGFFs(gff, exon_annot, mode='collapse')
            produced[s['library']] = str(written)
    return produced


def bulk_annotation_path(metadata_file, cell_annot=None):
    """Locate the bulk barcode->cluster file for a run: an explicit --cell_annot, else the one
    written beside the metadata by the extract phase."""
    if cell_annot:
        return str(cell_annot)
    base = os.path.dirname(os.path.abspath(str(metadata_file))) or os.getcwd()
    candidate = os.path.join(base, BULK_ANNOTATION_NAME)
    if os.path.exists(candidate):
        return candidate
    found = sorted(glob.glob(os.path.join(base, '*', BULK_ANNOTATION_NAME)))
    if found:
        return found[0]
    raise FileNotFoundError(
        f"No bulk cluster annotation found at {candidate}. Run `altanalyze3 blr` first, or pass "
        f"--cell_annot."
    )


def assert_no_gff_output_collision(sample_dict, work_dir=None):
    """Fail loudly when a sample's gff-output/ is the SAME directory as the run's gff-output/.

    Two different things are written to ``<dir>/gff-output/transcript_associations.txt``:

      * the PER-SAMPLE molecule structures, written next to each sample's read GFF by
        ``gff_process.consolidateLongReadGFFs(<one gff>, ...)``;
      * the REFERENCE transcript structures, written by
        ``isoform_collapse.reference.annotate_reference``, which calls
        ``consolidateLongReadGFFs([ref_gff], ...)``. Given a LIST, consolidateLongReadGFFs writes
        into ``os.getcwd()/gff-output`` (gff_process.py:731), not next to the input.

    When a sample's read GFF sits directly in the working directory, those two paths are identical
    and the reference annotation silently REPLACES the sample's molecule table. The collapse then
    reads reference rows where it expects millions of molecule rows.

    Observed: a 1,417,117-row PC3 molecule table replaced by GENCODE v45 reference rows.

    The fix is a directory layout, not a code change: keep each sample's BAM and read GFF in a
    SUBDIRECTORY of the working directory, which is what the single-cell metadata layout already
    does.
    """
    work_dir = os.path.abspath(str(work_dir or os.getcwd()))
    run_gff_output = os.path.join(work_dir, 'gff-output')
    clashing = []
    for uid, libs in sample_dict.items():
        for s in libs:
            sample_gff_output = os.path.abspath(
                os.path.join(os.path.dirname(os.path.abspath(str(s['gff']))), 'gff-output'))
            if sample_gff_output == run_gff_output:
                clashing.append((s['library'], str(s['gff'])))
    if clashing:
        listing = "\n  ".join(f"{lib}: {gff}" for lib, gff in clashing)
        raise ValueError(
            "Sample gff-output/ collides with the run's gff-output/ at\n"
            f"  {run_gff_output}\n"
            "so the reference annotation would overwrite the per-sample molecule structures.\n"
            f"Colliding sample(s):\n  {listing}\n"
            "Move each BAM into a SUBDIRECTORY of the working directory (for example "
            "<run>/bam/<library>.bam) and update the metadata 'bam' column, then re-run."
        )
    return True


def assert_extract_complete(sample_dict):
    """Every library must have its read GFF, molecule h5ad and structure table before the collapse."""
    missing = []
    for uid, libs in sample_dict.items():
        for s in libs:
            gff = str(s['gff'])
            ta = os.path.join(os.path.dirname(gff), 'gff-output', 'transcript_associations.txt')
            for path in (gff, str(s['matrix']), ta):
                # io_utils.exists resolves <path> OR <path>.gz. Phase 1 compresses
                # transcript_associations.txt, so os.path.exists alone reported every completed
                # sample as missing and blocked the collapse.
                if not _io.exists(path):
                    missing.append(path)
    if missing:
        raise FileNotFoundError(
            "Bulk extract outputs missing for {} path(s); run `altanalyze3 blr` first:\n  {}"
            .format(len(missing), "\n  ".join(missing))
        )


def report_extract_stats(stats, library, log=logging.info):
    """Report the extractor's own read counters with denominators (RULE 4: a run that exits zero is
    not a run that succeeded). Every rate is quoted against BOTH the total records and the primary
    mapped records, because the two differ whenever a BAM carries unmapped or supplementary rows."""
    total = stats.get('total_reads', 0)
    primary = stats.get('mapped_primary_reads', 0)
    lines = [f"[bulk:{library}] read accounting"]
    lines.append(f"  total records                 : {total:,}")
    lines.append(f"  mapped primary                : {primary:,}"
                 + (f" ({100.0 * primary / total:.1f}% of total)" if total else ""))
    for key, label in (('spliced_reads', 'spliced (>=1 N in CIGAR)'),
                       ('barcode_reads', 'barcode assigned'),
                       ('gene_assigned_spliced_reads', 'gene assigned'),
                       ('known_splice_reads', '>=1 known Ensembl splice site'),
                       ('kept_reads', 'kept (written to GFF/h5ad)')):
        value = stats.get(key, 0)
        pct_primary = f" ({100.0 * value / primary:.1f}% of primary)" if primary else ""
        lines.append(f"  {label:<30}: {value:,}{pct_primary}")
    kept = stats.get('kept_reads', 0)
    if primary and kept / primary < 0.5:
        lines.append(f"  *** RETENTION {100.0 * kept / primary:.1f}% of primary reads. Below 50%. "
                     f"Check the exon reference build and the BAM chromosome naming.")
    log("\n".join(lines))
    return dict(stats)
