"""Read-only replay of LSP1/RUNX1 annotation and collapse evidence loss.

Uses the actual source function bodies (AST extraction avoids optional scientific
dependencies); never changes source data, fits models, or writes collapse outputs.
Every supplied row is retained. Coordinates in renderer exports are descriptive,
not raw alignments; raw exon replay uses the upstream coordinate store separately.
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import gzip
from pathlib import Path
import re
import sqlite3
import struct
import sys


GENES = {"lsp1": "ENSG00000130592", "runx1": "ENSG00000159216"}
TARGETS = {"lsp1": "52795516.MDS5801-1229IR_CD34",
           "runx1": "28775276.EP_SRSF2_Ctrl"}


def load_functions(path, namespace, names=None, nested=False):
    tree = ast.parse(path.read_text())
    nodes = ast.walk(tree) if nested else tree.body
    funcs = [n for n in nodes if isinstance(n, ast.FunctionDef)
             and (names is None or n.name in names)]
    if names is not None and set(names) != {n.name for n in funcs}:
        raise RuntimeError(f"Required source functions unresolved: {path}: {names}")
    module = ast.Module(body=[ast.ImportFrom(module="__future__",
                        names=[ast.alias(name="annotations")], level=0)] + funcs,
                        type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), str(path), "exec"), namespace)


def load_runtime(lr):
    ns = dict(collections=collections, sys=sys, struct=struct,
              novelGene=0, novelGeneLoci=[1, 1], additional_junctions={})
    score_tree = ast.parse((lr / "isoform_collapse/scored_collapse.py").read_text())
    constants = [n for n in score_tree.body if isinstance(n, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id in
                         {"DEFAULT_MIN_READS", "LONG_BIN_FRACTION"} for t in n.targets)]
    exec(compile(ast.Module(body=constants, type_ignores=[]), str(lr), "exec"), ns)
    load_functions(lr / "gff_process.py", ns,
                   ["_split_token", "_strip_terminal_coords", "_filter_exon_intron_tokens",
                    "exon_str", "exonAnnotate", "findNovelSpliceSite", "rangeCheck",
                    "importEnsemblGenes", "uniqueList"], nested=True)
    ns.update(split_token=ns["_split_token"], strip_terminal_coords=ns["_strip_terminal_coords"],
              filter_exon_intron_tokens=ns["_filter_exon_intron_tokens"])
    load_functions(lr / "isoform_collapse_utils.py", ns)
    load_functions(lr / "isoform_collapse/scored_collapse.py", ns)
    load_functions(lr.parent / "bam/coord_store.py", ns, ["_decode_varint", "unpack_exons"])
    return ns


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--exon-annot", type=Path, required=True)
    ap.add_argument("--runx1-coord-db", type=Path, required=True)
    ap.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    args = ap.parse_args()
    ns = load_runtime(args.source)
    args.out.mkdir(parents=True, exist_ok=True)
    sources = {str(p): digest(p) for p in [args.source / "gff_process.py",
               args.source / "isoform_collapse_utils.py",
               args.source / "isoform_collapse/scored_collapse.py", args.exon_annot]}
    # Authoritative final transcript text, rather than a potentially stale viewer index.
    parents = {tag: {} for tag in GENES}
    with (args.root / "gff-output/transcript_associations.txt").open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                for tag, gene in GENES.items():
                    if p[0] == gene:
                        parents[tag][p[3]] = p[2]
    for tag in GENES:
        if TARGETS[tag] not in parents[tag]:
            raise RuntimeError(f"Required exemplar unresolved: {TARGETS[tag]}")
    summary = {"source_sha256": sources, "scope": "all rows in both supplied bin exports",
               "production_changed": False, "genes": {}}
    for tag, gene in GENES.items():
        p = args.root / f"analysis/hspc_regulator_splicing_switches/{tag}_isv/novel_isoform/{tag}_bin_reads_ALL.tsv"
        sources[str(p)] = digest(p)
        counts = collections.Counter()
        examples = {}
        parent = parents[tag][TARGETS[tag]]
        parent_core = ns["junction_core_tokens"](parent, gene)
        with p.open() as f, (args.out / f"{tag}_every_read.tsv").open("w") as out:
            reader = csv.DictReader(f, delimiter="\t")
            extra = ["legacy_structure", "E_only_trim_structure", "legacy_core",
                     "E_only_trim_core", "legacy_fits_assigned_core", "E_only_trim_fits_assigned_core"]
            writer = csv.DictWriter(out, fieldnames=reader.fieldnames + extra, delimiter="\t")
            writer.writeheader()
            for row in reader:
                raw = row["tokens"].split("|")
                legacy = ns["exon_str"](list(raw))[0]
                repaired = "|".join(ns["_strip_terminal_coords"](raw, gene))
                lc = ns["junction_core_tokens"](legacy, gene)
                rc = ns["junction_core_tokens"](repaired, gene)
                assigned = parents[tag].get(row["final_isoform"])
                if assigned is None:
                    raise RuntimeError(f"Required assigned exemplar unresolved: {row['final_isoform']}")
                pc = ns["junction_core_tokens"](assigned, gene)
                fits = ns["is_contiguous_subsequence"]
                row.update(legacy_structure=legacy, E_only_trim_structure=repaired,
                           legacy_core="|".join(lc), E_only_trim_core="|".join(rc),
                           legacy_fits_assigned_core=int(fits(lc, pc)),
                           E_only_trim_fits_assigned_core=int(fits(rc, pc)))
                writer.writerow(row)
                counts["rows"] += 1
                counts["legacy_core_fits_assigned"] += int(fits(lc, pc))
                counts["E_only_trim_core_fits_assigned"] += int(fits(rc, pc))
                if row["final_isoform"] == TARGETS[tag]:
                    counts["target_bin_rows"] += 1
                    counts["target_legacy_core_fits"] += int(fits(lc, pc))
                    counts["target_E_only_trim_core_fits"] += int(fits(rc, pc))
                    intronic_end = raw[0].startswith("I34.1_") if tag == "lsp1" else raw[-1].startswith("I42.1_")
                    counts["target_intronic_end_rows"] += int(intronic_end)
                    if intronic_end and "intronic_end" not in examples:
                        examples["intronic_end"] = row
        summary["genes"][tag] = dict(counts=dict(counts), target_structure=parent,
                                     target_core=parent_core, examples=examples)
    # Point lookup of a genuine stored RUNX1 structure: raw exon tuples are upstream
    # records, unlike the rendering segments in the bin export.
    key = "E42.4|E42.7"
    with sqlite3.connect(f"file:{args.runx1_coord_db}?mode=ro", uri=True) as con:
        record = con.execute("SELECT chrom,strand,blob FROM structures WHERE gene=? AND structure=?",
                             (GENES["runx1"], key)).fetchone()
    if record is None:
        raise RuntimeError(f"Required upstream coordinate source unresolved: {key}")
    chrom, strand, blob = record
    exons = ns["unpack_exons"](blob)
    ns["importEnsemblGenes"](str(args.exon_annot))
    raw = ns["exonAnnotate"](chrom, list(exons), strand, "diagnostic_stored_structure")[1]
    ref = {e: (a, b) for a, b, e in ns["geneData"][GENES["runx1"]]}
    summary["runx1_upstream_replay"] = dict(source=str(args.runx1_coord_db), key=key,
            exons=exons, raw_tokens=raw, legacy_structure=ns["exon_str"](list(raw))[0],
            E_only_trim_structure="|".join(ns["_strip_terminal_coords"](raw, GENES["runx1"])),
            reference_regions={k: v for k, v in ref.items() if k in ["E42.4", "E42.7", "I42.1", "E43.1"]})
    # Trace the actual novel2 branch for this stored read, then contrast with a
    # clearly labelled constructed control at the first subregion of E42.
    gid = GENES["runx1"]
    start_index = ns["exonCoordinates"][chrom, exons[0][1], strand, 1][-1]
    intronic_base = ns["findNovelSpliceSite"](gid, exons[0][0], strand).split("_")[0]
    proxy_exon = intronic_base.replace("I", "E")
    proxy_start, proxy_stop = ns["exonData"][gid, proxy_exon]
    end_index = ns["exonCoordinates"][chrom, proxy_stop, strand, 2][-1]
    control = ns["exonAnnotate"](chrom, [(exons[0][0], proxy_start)], strand,
                                  "constructed_first_subregion_control")[1]
    summary["runx1_novel2_branch"] = dict(start_index=start_index, end_index=end_index,
        intronic_base=intronic_base, proxy_exon=proxy_exon,
        same_index_guard_triggered=start_index == end_index,
        in_between_regions=ns["geneData"][gid][start_index:end_index],
        constructed_control_raw=control,
        constructed_control_legacy=ns["exon_str"](list(control))[0])
    # Interior paired intronic endpoints are also erased by the default core;
    # explicitly test the supplied interior LSP1 exon, not just terminal examples.
    interior = "E34.2|I34.1_1888012|I34.1_1888084|E35.1|E35.2"
    summary["interior_lsp1_replay"] = dict(structure=interior,
        core=ns["junction_core_tokens"](interior, GENES["lsp1"]))
    # Validate the exemplars directly against the final GFF, bypassing viewer caches.
    wanted = {"52795516": GENES["lsp1"], "28634443": GENES["lsp1"],
              "28775276": GENES["runx1"]}
    final_exons = {tid: [] for tid in wanted}
    with gzip.open(args.root / "gff-output/combined.gff.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] != "exon":
                continue
            if not any(gene in p[8] for gene in GENES.values()):
                continue
            match = re.search(r'transcript_id[ =]+["\']?([^;"\']+)', p[8])
            if match and match[1] in wanted and wanted[match[1]] in p[8]:
                final_exons[match[1]].append((int(p[3]), int(p[4])))
    if any(not ex for ex in final_exons.values()):
        raise RuntimeError(f"Required final GFF records unresolved: {final_exons}")
    summary["final_gff_exons"] = final_exons
    # True exon gaps distinguish actual splices from block changes inside a
    # continuously covered exon (the core currently treats both as junctions).
    summary["lsp1_terminal_junctions"] = {tid: [(a[1], b[0]) for a, b in
        zip(sorted(final_exons[tid]), sorted(final_exons[tid])[1:])][-2:]
        for tid in ["52795516", "28634443"]}
    # Reproduce fold edges with the actual collapse function on explicit two-structure
    # diagnostic fixtures. This is not a pooled rerun or evidence of a final winner.
    replay = {}
    for tag in GENES:
        gene = GENES[tag]
        e = summary["genes"][tag]["examples"]["intronic_end"]
        parent = summary["genes"][tag]["target_structure"]
        replay[tag] = {}
        for variant in ["legacy_structure", "E_only_trim_structure"]:
            child = e[variant]
            sr = {parent: 10, child: 20}
            tok = {s: tuple(ns["structure_tokens_for_containment"](s, gene)) for s in sr}
            core = {s: tuple(ns["junction_core_tokens"](s, gene)) for s in sr}
            result = ns["collapse_gene_jc"](sr, tok, core, {s: len(tok[s]) for s in sr})
            replay[tag][variant] = dict(child=child, child_core=core[child],
                                         assigned_to_target=result["assignment"].get(child) == parent)
    summary["two_structure_replay"] = replay
    (args.out / "diagnosis.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"counts": {k: v["counts"] for k, v in summary["genes"].items()},
                      "upstream_replay": summary["runx1_upstream_replay"],
                      "novel2_branch": summary["runx1_novel2_branch"],
                      "two_structure_replay": replay}, indent=2))


if __name__ == "__main__":
    main()
