#!/usr/bin/env python3
"""Build a sample's isoform-structure -> genomic-coordinate store from its GFF.

``gff_process`` writes ``<sample>/gff-output/structure_coords.sqlite`` while it annotates a sample,
but a run that only STAGES a sample's reads copies ``transcript_associations.txt`` and the gene
indexes without it. ISV then has to rebuild each read's blocks from its exon-ID string, which
cannot be done: ``exon_str`` drops the first and last token's coordinate and ``uniqueList`` drops a
repeated token, so two different block layouts share one string.

Everything the store needs is already on disk, so it is rebuilt by a join instead of a re-run:

  * ``transcript_associations.txt`` gives ``gene, strand, structure, transcript_id, file`` -- the
    structure string is exactly the key ``gff_process`` stores under;
  * the sample GFF gives ``transcript_id -> exon blocks`` in source order.

The result is the same table ``gff_process`` would have written, deduplicated by (gene, structure).

    python -m altanalyze3.components.visualization.isv_web.precompute_struct_coords \\
        --gff <sample>.gff.gz --associations <sample>/gff-output/transcript_associations.txt
"""
import argparse
import gzip
import os
import re
import sys

from ...bam import coord_store

_Q = re.compile(r'(\w+)\s+"([^"]*)"')
_E = re.compile(r'(\w+)=([^;]*)')


def _read_associations(path):
    """transcript_id -> (gene, structure). The file is gene, strand, structure, transcript, source."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            gene, structure, tid = f[0], f[2], f[3]
            if gene and structure and tid:
                out[tid] = (gene, structure)
    return out


def _gff_attrs(field):
    a = dict(_Q.findall(field))
    if "transcript_id" not in a:
        a.update(dict(_E.findall(field)))
    return a


def build(gff_path, associations_path, db_path, verbose=True):
    """Write db_path. Returns (structures written, transcripts read, transcripts with no structure)."""
    assoc = _read_associations(associations_path)
    if verbose:
        print(f"[struct_coords] associations: {len(assoc):,} transcripts", flush=True)

    tmp = db_path + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    writer = coord_store.StructCoordWriter(tmp)

    opener = gzip.open if gff_path.endswith(".gz") else open
    cur_tid = None
    cur = None
    n_tx = n_missing = n_written = 0
    seen = set()

    def flush(tid, rec):
        nonlocal n_written, n_missing
        if tid is None or not rec or not rec["exons"]:
            return
        hit = assoc.get(tid)
        if hit is None:
            n_missing += 1
            return
        gene, structure = hit
        if (gene, structure) in seen:
            return
        seen.add((gene, structure))
        writer.add(gene, structure, rec["chrom"], rec["strand"], rec["source"], rec["exons"])
        n_written += 1

    with opener(gff_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = _gff_attrs(f[8])
            tid = a.get("transcript_id")
            if not tid:
                continue
            if tid != cur_tid:
                flush(cur_tid, cur)
                cur_tid = tid
                cur = {"chrom": f[0], "strand": f[6], "source": f[1], "exons": []}
                n_tx += 1
            cur["exons"].append((int(f[3]), int(f[4])))
    flush(cur_tid, cur)
    writer.close()
    os.replace(tmp, db_path)
    if verbose:
        print(f"[struct_coords] transcripts in GFF: {n_tx:,}", flush=True)
        print(f"[struct_coords] structures written: {n_written:,}", flush=True)
        print(f"[struct_coords] transcripts with no association row: {n_missing:,}", flush=True)
        print(f"[struct_coords] wrote {db_path}", flush=True)
    return n_written, n_tx, n_missing


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="the sample GFF the reads came from")
    ap.add_argument("--associations", default=None,
                    help="transcript_associations.txt (default: <gff dir>/gff-output/transcript_associations.txt)")
    ap.add_argument("--out", default=None,
                    help="output db (default: beside the associations file, structure_coords.sqlite)")
    ap.add_argument("--force", action="store_true", help="rebuild even when the db already exists")
    a = ap.parse_args(argv)

    assoc = a.associations or os.path.join(os.path.dirname(a.gff), "gff-output",
                                           "transcript_associations.txt")
    out = a.out or os.path.join(os.path.dirname(assoc), "structure_coords.sqlite")
    if os.path.exists(out) and not a.force:
        print(f"[struct_coords] already present, keeping {out}", flush=True)
        return 0
    if not os.path.exists(assoc):
        print(f"[struct_coords] no associations file: {assoc}", file=sys.stderr, flush=True)
        return 2
    build(a.gff, assoc, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
