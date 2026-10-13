#!/usr/bin/env python3
"""Build the ISV web MOLECULE -> FINAL-ISOFORM index for the "group reads by collapsed isoform" molecule
view. The isoform-collapse protocol assigns every read to a final isoform; this precomputes that lookup so
the molecule view can group/colour reads (and key their FASTA link-outs) by the final isoform they collapsed
into, instead of clustering structures live.

Read -> final chain (verified): per-sample ``tier1/<sample>.mol2struct.tsv`` (molecule_id -> structure, in
the collapse's terminal-normalised token form) JOINed to ``gff-output/FINAL_structure_to_exemplar.tsv``
(gene, structure -> final_isoform_id). The raw per-sample ``transcripts.db`` structures do NOT match FINAL
(they keep terminal coords); ``mol2struct.tsv`` is the authoritative key (~97% of reads resolve; the rest
were dropped in collapse and fall back to their own id in the viewer).

Output (alongside the reads index, ``<root>/_isv_web_cache/mol_index/``):
    <library>.mol2final.db   table m2f(mol TEXT PRIMARY KEY, final TEXT)

Idempotent: a sample db newer than both its mol2struct.tsv and FINAL tsv is skipped; the shared
``_s2f.db`` (FINAL loaded once, indexed by (gene,structure)) is built once and reused across samples.

Usage:
    python3 -m altanalyze3.components.visualization.isv_web.precompute_final_isoform_index \
        --root /path/to/dataset_root [--gff_output DIR]
or programmatically: build_mol2final_index(lib2mol2struct, final_tsv, out_dir)
"""
import os
import glob
import sqlite3
import time


def mol2final_db_path(out_dir, library):
    return os.path.join(out_dir, f"{library}.mol2final.db")


def reads_final_db_path(out_dir, library):
    return os.path.join(out_dir, f"{library}.reads_final.db")


def build_reads_final(reads_db, m2f_db, out_db, workdir=None, log=print):
    """<library>.reads_final.db: one row per (gene, molecule, cell type) with its read count and FINAL isoform, in a
    WITHOUT ROWID table keyed (gene, mol, ct), so one gene's reads and finals come back in one contiguous range read.
    final is NULL where <library>.mol2final.db has no entry (the query then keeps the read under its own id).
    workdir: a job-private scratch directory (LSF $TMPDIR); inputs are copied there because the join probes m2f once
    per read, which is slow on network storage. Nothing is deleted here; the scheduler clears its scratch."""
    import shutil
    if (os.path.exists(out_db) and os.path.getmtime(out_db) >= os.path.getmtime(reads_db)
            and os.path.getmtime(out_db) >= os.path.getmtime(m2f_db)):
        log(f"[reads_final] {os.path.basename(out_db)} up to date; skip"); return -1
    t0 = time.time()
    src_r, src_f, tmp = reads_db, m2f_db, out_db + ".tmp"
    if workdir:
        src_r = os.path.join(workdir, os.path.basename(reads_db)); shutil.copyfile(reads_db, src_r)
        src_f = os.path.join(workdir, os.path.basename(m2f_db)); shutil.copyfile(m2f_db, src_f)
        tmp = os.path.join(workdir, os.path.basename(out_db) + ".tmp")
        log(f"[reads_final] copied inputs to {workdir} ({time.time() - t0:.0f}s)")
    for p_ in {tmp, out_db + ".tmp"}:
        if os.path.exists(p_):
            raise FileExistsError(f"[reads_final] leftover temporary file {p_}; inspect it before rebuilding")
    con = sqlite3.connect(tmp); cur = con.cursor()
    for pragma in ("journal_mode=OFF", "synchronous=OFF", "temp_store=MEMORY", "cache_size=-4000000"):
        cur.execute(f"PRAGMA {pragma}")
    cur.execute("ATTACH DATABASE ? AS r", (src_r,)); cur.execute("ATTACH DATABASE ? AS f", (src_f,))
    cur.execute("CREATE TABLE rf (gene TEXT NOT NULL, mol TEXT NOT NULL, ct TEXT NOT NULL, cnt REAL, final TEXT, "
                "PRIMARY KEY (gene, mol, ct)) WITHOUT ROWID")
    cur.execute("INSERT INTO rf SELECT x.gene, x.mol, x.ct, SUM(x.cnt), y.final FROM r.mol x "
                "LEFT JOIN f.m2f y ON y.mol = x.mol GROUP BY x.gene, x.mol, x.ct ORDER BY x.gene, x.mol, x.ct")
    n_src, s_src = cur.execute("SELECT COUNT(*), SUM(cnt) FROM (SELECT gene, mol, ct, SUM(cnt) AS cnt FROM r.mol "
                               "GROUP BY gene, mol, ct)").fetchone()
    n_out, s_out, n_fin = cur.execute("SELECT COUNT(*), SUM(cnt), COUNT(final) FROM rf").fetchone()
    con.commit(); cur.execute("DETACH DATABASE r"); cur.execute("DETACH DATABASE f"); con.close()
    if n_out != n_src or abs((s_out or 0) - (s_src or 0)) > 1e-6 * max(1.0, s_src or 0):
        raise AssertionError(f"[reads_final] {os.path.basename(out_db)}: rows {n_out} vs {n_src}, reads {s_out} vs {s_src}")
    if workdir:
        shutil.copyfile(tmp, out_db + ".tmp")
    os.replace(out_db + ".tmp", out_db)
    log(f"[reads_final] {os.path.basename(out_db)}: {n_out:,} rows, {s_out:,.0f} reads, final resolved {n_fin:,} "
        f"({100 * n_fin / max(1, n_out):.1f}%) ({time.time() - t0:.0f}s)")
    return n_out


def _stream_tsv(path, ncols):
    with open(path) as f:
        next(f, None)                                            # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= ncols:
                yield p


def build_s2f(final_tsv, s2f_db, log=print):
    """Load FINAL_structure_to_exemplar.tsv into SQLite once: s2f(gene, structure, final), indexed."""
    if os.path.exists(s2f_db) and os.path.getmtime(s2f_db) >= os.path.getmtime(final_tsv):
        log(f"[mol2final] {os.path.basename(s2f_db)} up to date; reuse"); return s2f_db
    t0 = time.time(); tmp = s2f_db + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp); cur = con.cursor()
    cur.execute("PRAGMA journal_mode=OFF"); cur.execute("PRAGMA synchronous=OFF")
    cur.execute("CREATE TABLE s2f (gene TEXT, structure TEXT, final TEXT)")
    batch, n = [], 0
    for p in _stream_tsv(final_tsv, 3):
        batch.append((p[0], p[1], p[2])); n += 1
        if len(batch) >= 100000:
            cur.executemany("INSERT INTO s2f VALUES (?,?,?)", batch); batch = []
    if batch:
        cur.executemany("INSERT INTO s2f VALUES (?,?,?)", batch)
    cur.execute("CREATE INDEX ix_s2f ON s2f(gene, structure)")
    con.commit(); con.close(); os.replace(tmp, s2f_db)
    log(f"[mol2final] s2f built: {n:,} structures ({time.time() - t0:.1f}s)")
    return s2f_db


def build_one(mol2struct_tsv, s2f_db, out_db, log=print):
    """molecule_id -> final_isoform_id for one sample, via JOIN on (gene, structure)."""
    if (os.path.exists(out_db) and os.path.getmtime(out_db) >= os.path.getmtime(mol2struct_tsv)
            and os.path.getmtime(out_db) >= os.path.getmtime(s2f_db)):
        log(f"[mol2final] {os.path.basename(out_db)} up to date; skip"); return -1
    t0 = time.time(); tmp = out_db + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp); cur = con.cursor()
    cur.execute("PRAGMA journal_mode=OFF"); cur.execute("PRAGMA synchronous=OFF")
    cur.execute("CREATE TABLE m (mol TEXT, gene TEXT, structure TEXT)")
    batch = []
    for p in _stream_tsv(mol2struct_tsv, 3):
        batch.append((p[0], p[1], p[2]))
        if len(batch) >= 100000:
            cur.executemany("INSERT INTO m VALUES (?,?,?)", batch); batch = []
    if batch:
        cur.executemany("INSERT INTO m VALUES (?,?,?)", batch)
    n_mol = cur.execute("SELECT COUNT(*) FROM m").fetchone()[0]
    cur.execute(f"ATTACH DATABASE '{s2f_db}' AS s")
    cur.execute("CREATE TABLE m2f (mol TEXT PRIMARY KEY, final TEXT)")
    cur.execute("INSERT OR IGNORE INTO m2f SELECT m.mol, s.final FROM m "
                "JOIN s.s2f s ON m.gene = s.gene AND m.structure = s.structure")
    n_res = cur.execute("SELECT COUNT(*) FROM m2f").fetchone()[0]
    cur.execute("DROP TABLE m"); con.commit(); cur.execute("DETACH DATABASE s"); con.close()
    os.replace(tmp, out_db)
    log(f"[mol2final] {os.path.basename(out_db)}: {n_res:,}/{n_mol:,} reads resolved "
        f"({100 * n_res / max(1, n_mol):.1f}%) ({time.time() - t0:.1f}s)")
    return n_res


def build_mol2final_index(lib2mol2struct, final_tsv, out_dir, log=print):
    """lib2mol2struct: {library: mol2struct.tsv path}. Builds <library>.mol2final.db in out_dir."""
    if not os.path.exists(final_tsv):
        log(f"[mol2final] FINAL_structure_to_exemplar.tsv not found ({final_tsv}); skip"); return {}
    os.makedirs(out_dir, exist_ok=True)
    s2f_db = build_s2f(final_tsv, os.path.join(out_dir, "_s2f.db"), log=log)
    built = {}
    for lib, m2s in lib2mol2struct.items():
        if not (m2s and os.path.exists(m2s)):
            log(f"[mol2final] {lib}: no mol2struct.tsv; skip"); continue
        built[lib] = build_one(m2s, s2f_db, mol2final_db_path(out_dir, lib), log=log)
    return built


def discover(root, gff_output=None):
    """Find FINAL tsv + {library: mol2struct.tsv} under a dataset root."""
    gff = gff_output or os.path.join(root, "gff-output")
    final_tsv = os.path.join(gff, "FINAL_structure_to_exemplar.tsv")
    lib2m2s = {}
    for f in glob.glob(os.path.join(root, "**", "tier1", "*.mol2struct.tsv"), recursive=True):
        lib2m2s.setdefault(os.path.basename(f)[:-len(".mol2struct.tsv")], f)
    return final_tsv, lib2m2s


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Build isv_web molecule->final-isoform index.")
    ap.add_argument("--root", required=True, help="dataset root (holds gff-output/ and */tier1/)")
    ap.add_argument("--gff_output", default=None)
    a = ap.parse_args(argv)
    final_tsv, lib2m2s = discover(a.root, a.gff_output)
    print(f"[mol2final] FINAL={final_tsv}\n[mol2final] libraries: {sorted(lib2m2s)}", flush=True)
    out_dir = os.path.join(a.root, "_isv_web_cache", "mol_index")
    build_mol2final_index(lib2m2s, final_tsv, out_dir)


if __name__ == "__main__":
    main()
