#!/usr/bin/env python3
"""Index the REGISTERED final-isoform structures out of ``gff-output/combined.gff.gz``.

Why this exists
---------------
The ISV structure track used to draw, for each final isoform, the merged exons of one
REPRESENTATIVE READ picked by ``render()``. That read is arbitrary with respect to
completeness -- when every read carries count 1 the choice degenerates to lexicographic
molecule-id order -- so a 5'-truncated read could stand in for a full-length transcript.
Measured on FLNB (ENSG00000136068): ENST00000295956 was drawn with 5 merged exon blocks
spanning 58,159,623-58,170,783 (11,161 bp) while the registered transcript has 46 exons
spanning 58,008,422-58,172,251 (163,830 bp) -- 151,201 bp of 5' structure missing.

``combined.gff.gz`` is the authoritative registration of every final isoform, but it has
no tabix index, so random access needs this one-pass extract.

Output
------
``<root>/_isv_web_cache/final_structures.db``

    structures(gene, tid, tid_base, chrom, strand, source, exons)
    exons = "start-end,start-end,..." ascending, 1-based inclusive, as written in the GFF.
    INDEX on (gene, tid_base) and on (tid_base).

Key forms -- the GFF carries two attribute dialects and two id conventions:

    bam     gene_id "ENSG00000223972";transcript_id "21863493";
            -> novel finals are '<molecule>.<library>' in the viewer; the GFF holds the
               BARE '<molecule>'. tid_base = '21863493'.
    HAVANA  ID=ENST00000295956;gene_id=ENSG00000136068.16;transcript_id=ENST00000295956.9
            -> known finals are unversioned in the viewer. tid_base = 'ENST00000295956'.

Both gene_id and transcript_id are stored version-stripped in the *_base columns so the
viewer's unversioned ids match directly.

Run:
  python -m altanalyze3.components.visualization.isv_web.precompute_gff_structures \\
      --root /path/to/dataset
"""
import os
import re
import gzip
import sqlite3
import argparse
from collections import defaultdict

CACHE_DIRNAME = "_isv_web_cache"
DB_NAME = "final_structures.db"

_ATTR_QUOTED = re.compile(r'(\w+)\s+"([^"]*)"')
_ATTR_EQUALS = re.compile(r'(\w+)=([^;]*)')


def parse_attrs(field):
    """Parse both GFF dialects: `key "value";` (bam) and `key=value;` (HAVANA/GENCODE)."""
    out = dict(_ATTR_QUOTED.findall(field))
    if not out or "transcript_id" not in out:
        out.update(dict(_ATTR_EQUALS.findall(field)))
    return out


def strip_version(x):
    """ENSG00000136068.16 -> ENSG00000136068 ; ENST00000295956.9 -> ENST00000295956.
    A bare numeric bam id ('21863493') has no version and is returned unchanged."""
    x = str(x or "")
    if "." in x and x.split(".", 1)[0].startswith("ENS"):
        return x.split(".", 1)[0]
    return x


def gff_path(root, explicit=None):
    if explicit:
        return explicit
    for name in ("combined.gff.gz", "combined.gff"):
        p = os.path.join(root, "gff-output", name)
        if os.path.exists(p):
            return p
    return os.path.join(root, "gff-output", "combined.gff.gz")


def db_path(root):
    return os.path.join(root, CACHE_DIRNAME, DB_NAME)


def build(root, gff=None, log=print, batch=200_000):
    """One streaming pass over the GFF -> SQLite. Memory stays flat.

    An earlier version accumulated every transcript in a Python dict and was OOM-killed on the
    login node at 764k transcripts / 6M exon rows. Exon rows are now staged into a table as they
    are read, then re-read through ``ORDER BY gene, tid, s`` -- which makes each transcript's rows
    contiguous AND already sorted -- and folded one transcript at a time into ``structures``.
    Two invariants fail loudly: staged row count == parsed row count, and the exon count implied
    by ``structures`` == the staged row count.
    """
    src = gff_path(root, gff)
    if not os.path.exists(src):
        raise FileNotFoundError(f"GFF not found: {src}")
    out = db_path(root)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)

    log(f"[gff-structures] reading {src}", flush=True)
    con = sqlite3.connect(tmp)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA cache_size=-262144")
    con.execute("""CREATE TABLE ex (gene TEXT, tid TEXT, chrom TEXT, strand TEXT,
                   source TEXT, s INTEGER, e INTEGER)""")

    opener = gzip.open if src.endswith(".gz") else open
    n_lines = n_exon = 0
    buf = []
    with opener(src, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            n_lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            a = parse_attrs(f[8])
            tid = a.get("transcript_id")
            gid = a.get("gene_id")
            if not tid or not gid:
                continue
            n_exon += 1
            buf.append((strip_version(gid), tid, f[0], f[6], f[1], int(f[3]), int(f[4])))
            if len(buf) >= batch:
                con.executemany("INSERT INTO ex VALUES (?,?,?,?,?,?,?)", buf)
                buf = []
                if n_exon % 2_000_000 < batch:
                    log(f"[gff-structures]   staged {n_exon:,} exon rows", flush=True)
    if buf:
        con.executemany("INSERT INTO ex VALUES (?,?,?,?,?,?,?)", buf)
    con.commit()
    staged = con.execute("SELECT COUNT(*) FROM ex").fetchone()[0]
    log(f"[gff-structures] scanned {n_lines:,} feature lines; staged {staged:,} exon rows", flush=True)
    if staged != n_exon:
        raise RuntimeError(f"staging lost rows: parsed {n_exon}, stored {staged}")
    if not staged:
        raise RuntimeError(f"no exon features parsed from {src}; refusing to write an empty index")

    con.execute("""CREATE TABLE structures (gene TEXT, tid TEXT, tid_base TEXT,
                   chrom TEXT, strand TEXT, source TEXT, exons TEXT)""")
    log("[gff-structures] folding exons per transcript (ORDER BY gene, tid, s) ...", flush=True)
    cur = con.execute("SELECT gene, tid, chrom, strand, source, s, e FROM ex ORDER BY gene, tid, s")
    ins = con.cursor()
    state = {"rows": [], "n_tx": 0, "key": None, "meta": None, "segs": []}

    def flush_tx():
        if state["key"] is None:
            return
        g, t = state["key"]
        c, sd, so = state["meta"]
        state["rows"].append((g, t, strip_version(t), c, sd, so, ",".join(state["segs"])))
        state["n_tx"] += 1
        if len(state["rows"]) >= 100_000:
            ins.executemany("INSERT INTO structures VALUES (?,?,?,?,?,?,?)", state["rows"])
            state["rows"] = []

    for gene, tid, chrom, strand, source, st, en in cur:
        k = (gene, tid)
        if k != state["key"]:
            flush_tx()
            state["key"], state["meta"], state["segs"] = k, (chrom, strand, source), []
        state["segs"].append(f"{st}-{en}")
    flush_tx()
    if state["rows"]:
        ins.executemany("INSERT INTO structures VALUES (?,?,?,?,?,?,?)", state["rows"])
    con.commit()

    n = con.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
    if n != state["n_tx"]:
        raise RuntimeError(f"row-count mismatch: folded {state['n_tx']} transcripts, wrote {n}")
    exploded = con.execute(
        "SELECT SUM(LENGTH(exons) - LENGTH(REPLACE(exons, ',', '')) + 1) FROM structures").fetchone()[0]
    if exploded != staged:
        raise RuntimeError(f"exon conservation failed: staged {staged}, in structures {exploded}")
    log(f"[gff-structures] {n:,} transcripts, {exploded:,} exons (conserved)", flush=True)

    con.execute("DROP TABLE ex")
    con.execute("CREATE INDEX ix_gene_base ON structures(gene, tid_base)")
    con.execute("CREATE INDEX ix_base ON structures(tid_base)")
    con.commit()
    con.execute("VACUUM")
    con.close()
    os.replace(tmp, out)
    log(f"[gff-structures] wrote {n:,} transcripts -> {out} ({os.path.getsize(out):,} bytes)", flush=True)
    return out


# ---------------------------------------------------------------- read side

def final_to_tid_base(final_isoform_id):
    """Viewer final-isoform id -> the transcript_id key used in combined.gff.gz.

      'ENST00000295956'            -> 'ENST00000295956'        (HAVANA, version stripped in the index)
      'ENST00000295956.9'          -> 'ENST00000295956'
      '24303983.Y2982'             -> '24303983'               (bam records hold the BARE molecule id)
      '11897914.BF21-CD34'         -> '11897914'
    """
    f = str(final_isoform_id or "")
    if not f:
        return None
    if f.startswith("ENS"):
        return f.split(".", 1)[0]
    return f.split(".", 1)[0] if "." in f else f


class StructureIndex:
    """Read-only lookup of registered exon coordinates by (gene, final_isoform_id)."""

    def __init__(self, path):
        self.path = path
        self._con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, check_same_thread=False)
        self._cache = {}

    @classmethod
    def open(cls, root, log=print):
        p = db_path(root)
        if not os.path.exists(p):
            log(f"[gff-structures] index absent ({p}); structure track falls back to a representative read")
            return None
        return cls(p)

    def exons(self, gene, final_isoform_id):
        """-> [(start, end), ...] ascending, or None when the isoform is not registered."""
        base = final_to_tid_base(final_isoform_id)
        if not base:
            return None
        key = (gene, base)
        if key in self._cache:
            return self._cache[key]
        row = self._con.execute(
            "SELECT exons FROM structures WHERE gene=? AND tid_base=?", (gene, base)).fetchone()
        if row is None:
            row = self._con.execute(
                "SELECT exons FROM structures WHERE tid_base=?", (base,)).fetchone()
        val = None
        if row and row[0]:
            val = [tuple(int(v) for v in seg.split("-")) for seg in row[0].split(",") if seg]
        self._cache[key] = val
        return val


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m altanalyze3.components.visualization.isv_web.precompute_gff_structures",
        description="Index registered final-isoform structures from gff-output/combined.gff.gz.")
    ap.add_argument("--root", required=True, help="dataset dir holding gff-output/ and _isv_web_cache/")
    ap.add_argument("--gff", default=None, help="explicit path to combined.gff.gz (default: <root>/gff-output/)")
    a = ap.parse_args(argv)
    build(a.root, a.gff)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
