#!/usr/bin/env python3
import argparse
import csv
import os
import sys
import heapq
import shutil
from pathlib import Path
from itertools import groupby


def _iter_input_files(inputs, input_dir):
    files = []
    if input_dir:
        for path in Path(input_dir).rglob("transcript_associations.txt"):
            files.append(str(path))
    if inputs:
        files.extend(inputs)
    return sorted(set(files))


def _infer_sample_from_path(path):
    base = os.path.basename(path)
    if base == "transcript_associations.txt":
        parent = os.path.basename(os.path.dirname(path))
        return parent or "sample"
    return os.path.splitext(base)[0]


def _write_sorted_runs(raw_path, runs_dir, max_records):
    runs_dir.mkdir(parents=True, exist_ok=True)
    runs = []
    buffer_records = []
    with open(raw_path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            buffer_records.append(row)
            if len(buffer_records) >= max_records:
                run_path = runs_dir / f"run_{len(runs):05d}.tsv"
                buffer_records.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
                with open(run_path, "w", newline="") as out_handle:
                    writer = csv.writer(out_handle, delimiter="\t")
                    writer.writerows(buffer_records)
                runs.append(run_path)
                buffer_records.clear()
    if buffer_records:
        run_path = runs_dir / f"run_{len(runs):05d}.tsv"
        buffer_records.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
        with open(run_path, "w", newline="") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t")
            writer.writerows(buffer_records)
        runs.append(run_path)
    return runs


def _merge_sorted_runs(run_paths, sorted_path):
    handles = [open(path, "r", newline="") for path in run_paths]
    heap = []
    for idx, handle in enumerate(handles):
        line = handle.readline()
        if line:
            row = line.rstrip("\n").split("\t")
            heapq.heappush(heap, (row[0], row[1], row[2], row[3], row[4], idx))
    with open(sorted_path, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["gene", "strand", "structure", "sample", "molecule_id"])
        while heap:
            gene, strand, structure, sample, molecule_id, idx = heapq.heappop(heap)
            writer.writerow([gene, strand, structure, sample, molecule_id])
            line = handles[idx].readline()
            if line:
                row = line.rstrip("\n").split("\t")
                heapq.heappush(heap, (row[0], row[1], row[2], row[3], row[4], idx))
    for handle in handles:
        handle.close()


def _longest_common_substring_length(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0
    if len(seq_a) > len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    prev = [0] * (len(seq_a) + 1)
    best = 0
    for token in seq_b:
        current = [0]
        for idx, a_token in enumerate(seq_a, start=1):
            if token == a_token:
                val = prev[idx - 1] + 1
            else:
                val = 0
            current.append(val)
            if val > best:
                best = val
        prev = current
    return best


def _substring_similarity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    short = seq_a if len(seq_a) <= len(seq_b) else seq_b
    long = seq_b if short is seq_a else seq_a
    lcs = _longest_common_substring_length(short, long)
    return lcs / float(len(short))


def _jaccard_similarity(tokens_a, tokens_b):
    set_a = set(tokens_a)
    set_b = set(tokens_b)
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def _tokenize(structure):
    return [tok for tok in structure.split("|") if tok]


def parse_and_write_raw(inputs, output_path):
    sample_names = set()
    with open(output_path, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["gene", "strand", "structure", "sample", "molecule_id"])
        for path in inputs:
            fallback_sample = _infer_sample_from_path(path)
            with open(path, "r") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    gene, strand, structure, molecule_id = parts[:4]
                    sample = parts[4] if len(parts) > 4 and parts[4] else fallback_sample
                    if not gene or not structure:
                        continue
                    writer.writerow([gene, strand, structure, sample, molecule_id])
                    sample_names.add(sample)
    return sorted(sample_names)


def aggregate_structure_sample_counts(sorted_raw, out_counts):
    with open(sorted_raw, "r", newline="") as handle, open(out_counts, "w", newline="") as out_handle:
        reader = csv.reader(handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t")
        header = next(reader, None)
        writer.writerow(["gene", "strand", "structure", "sample", "count"])
        prev_key = None
        count = 0
        for row in reader:
            if not row:
                continue
            key = (row[0], row[1], row[2], row[3])
            if prev_key is None:
                prev_key = key
                count = 1
                continue
            if key == prev_key:
                count += 1
                continue
            writer.writerow([prev_key[0], prev_key[1], prev_key[2], prev_key[3], count])
            prev_key = key
            count = 1
        if prev_key is not None:
            writer.writerow([prev_key[0], prev_key[1], prev_key[2], prev_key[3], count])


def aggregate_structure_totals(counts_path, totals_path, gene_counts_path):
    with open(counts_path, "r", newline="") as handle, \
            open(totals_path, "w", newline="") as totals_handle, \
            open(gene_counts_path, "w", newline="") as gene_handle:
        reader = csv.reader(handle, delimiter="\t")
        totals_writer = csv.writer(totals_handle, delimiter="\t")
        gene_writer = csv.writer(gene_handle, delimiter="\t")
        header = next(reader, None)
        totals_writer.writerow(["gene", "strand", "structure", "total_count", "sample_count"])
        gene_writer.writerow(["gene", "sample_count"])

        current_gene = None
        gene_samples = set()
        current_key = None
        total_count = 0
        sample_count = 0
        for row in reader:
            if not row:
                continue
            gene, strand, structure, sample, count = row[0], row[1], row[2], row[3], int(row[4])
            key = (gene, strand, structure)
            if current_gene is None:
                current_gene = gene
            if gene != current_gene:
                gene_writer.writerow([current_gene, len(gene_samples)])
                gene_samples = set()
                current_gene = gene
            gene_samples.add(sample)
            if current_key is None:
                current_key = key
                total_count = count
                sample_count = 1
                prev_sample = sample
                continue
            if key == current_key:
                total_count += count
                if sample != prev_sample:
                    sample_count += 1
                    prev_sample = sample
                continue
            totals_writer.writerow([current_key[0], current_key[1], current_key[2], total_count, sample_count])
            current_key = key
            total_count = count
            sample_count = 1
            prev_sample = sample
        if current_key is not None:
            totals_writer.writerow([current_key[0], current_key[1], current_key[2], total_count, sample_count])
        if current_gene is not None:
            gene_writer.writerow([current_gene, len(gene_samples)])


def _iter_gene_maps(path):
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        for gene, rows in groupby(reader, key=lambda r: r[0]):
            mapping = {}
            for row in rows:
                mapping[(row[1], row[2])] = row[3]
            yield gene, mapping


def collapse_structures(totals_path, gene_counts_path, collapse_map_path,
                        collapsed_totals_path, low_per_sample, min_similarity, strategy):
    gene_sample_counts = {}
    with open(gene_counts_path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        _ = next(reader, None)
        for row in reader:
            if not row:
                continue
            gene_sample_counts[row[0]] = int(row[1])

    def similarity(a_tokens, b_tokens):
        if strategy == "jaccard":
            return _jaccard_similarity(a_tokens, b_tokens)
        return _substring_similarity(a_tokens, b_tokens)

    with open(totals_path, "r", newline="") as handle, \
            open(collapse_map_path, "w", newline="") as map_handle, \
            open(collapsed_totals_path, "w", newline="") as totals_handle:
        reader = csv.reader(handle, delimiter="\t")
        map_writer = csv.writer(map_handle, delimiter="\t")
        totals_writer = csv.writer(totals_handle, delimiter="\t")
        _ = next(reader, None)
        map_writer.writerow(["gene", "strand", "structure", "collapsed_structure", "similarity", "total_count"])
        totals_writer.writerow(["gene", "strand", "structure", "total_count"])

        for gene, rows in groupby(reader, key=lambda r: r[0]):
            records = []
            for row in rows:
                strand, structure, total_count, sample_count = row[1], row[2], int(row[3]), int(row[4])
                records.append((strand, structure, total_count, sample_count))

            if not records:
                continue

            sample_total = gene_sample_counts.get(gene, 1)
            low_threshold = max(1, int(round(sample_total * low_per_sample)))
            high_conf = [r for r in records if r[2] >= low_threshold]
            if not high_conf:
                high_conf = [max(records, key=lambda r: r[2])]

            high_tokens = {r[1]: _tokenize(r[1]) for r in high_conf}
            collapsed_totals = {}

            for strand, structure, total_count, sample_count in records:
                if any(structure == h[1] for h in high_conf):
                    collapsed = structure
                    score = 1.0
                else:
                    tokens = _tokenize(structure)
                    best_score = -1.0
                    best_struct = structure
                    for h_struct, h_tokens in high_tokens.items():
                        score = similarity(tokens, h_tokens)
                        if score > best_score:
                            best_score = score
                            best_struct = h_struct
                    if best_score >= min_similarity:
                        collapsed = best_struct
                        score = best_score
                    else:
                        collapsed = structure
                        score = best_score
                map_writer.writerow([gene, strand, structure, collapsed, f"{score:.3f}", total_count])
                collapsed_totals[(strand, collapsed)] = collapsed_totals.get((strand, collapsed), 0) + total_count

            for (strand, collapsed), total_count in collapsed_totals.items():
                totals_writer.writerow([gene, strand, collapsed, total_count])


def collapse_sample_counts(counts_path, collapse_map_path, out_path):
    map_iter = _iter_gene_maps(collapse_map_path)
    current_map_gene = None
    current_map = {}
    try:
        current_map_gene, current_map = next(map_iter)
    except StopIteration:
        current_map_gene, current_map = None, {}

    with open(counts_path, "r", newline="") as handle, open(out_path, "w", newline="") as out_handle:
        reader = csv.reader(handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t")
        _ = next(reader, None)
        writer.writerow(["gene", "strand", "structure", "sample", "count"])

        prev_key = None
        acc = 0
        for row in reader:
            if not row:
                continue
            gene, strand, structure, sample, count = row[0], row[1], row[2], row[3], int(row[4])
            while current_map_gene is not None and current_map_gene < gene:
                try:
                    current_map_gene, current_map = next(map_iter)
                except StopIteration:
                    current_map_gene, current_map = None, {}
            if current_map_gene == gene:
                collapsed = current_map.get((strand, structure), structure)
            else:
                collapsed = structure
            key = (gene, strand, collapsed, sample)
            if prev_key is None:
                prev_key = key
                acc = count
                continue
            if key == prev_key:
                acc += count
                continue
            writer.writerow([prev_key[0], prev_key[1], prev_key[2], prev_key[3], acc])
            prev_key = key
            acc = count
        if prev_key is not None:
            writer.writerow([prev_key[0], prev_key[1], prev_key[2], prev_key[3], acc])


def collapse_associations(sorted_raw, collapse_map_path, out_path):
    map_iter = _iter_gene_maps(collapse_map_path)
    current_map_gene = None
    current_map = {}
    try:
        current_map_gene, current_map = next(map_iter)
    except StopIteration:
        current_map_gene, current_map = None, {}

    with open(sorted_raw, "r", newline="") as handle, open(out_path, "w", newline="") as out_handle:
        reader = csv.reader(handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t")
        _ = next(reader, None)
        writer.writerow(["gene", "strand", "structure", "molecule_id", "sample"])
        for row in reader:
            if not row:
                continue
            gene, strand, structure, sample, molecule_id = row
            while current_map_gene is not None and current_map_gene < gene:
                try:
                    current_map_gene, current_map = next(map_iter)
                except StopIteration:
                    current_map_gene, current_map = None, {}
            if current_map_gene == gene:
                collapsed = current_map.get((strand, structure), structure)
            else:
                collapsed = structure
            writer.writerow([gene, strand, collapsed, molecule_id, sample])


def main():
    parser = argparse.ArgumentParser(
        description="Collapse transcript_associations.txt files into frequency-aware isoform models."
    )
    parser.add_argument("--input", dest="inputs", nargs="*", default=None,
                        help="Transcript association files (transcript_associations.txt).")
    parser.add_argument("--input-dir", dest="input_dir", default=None,
                        help="Directory to scan for transcript_associations.txt files.")
    parser.add_argument("--outdir", default="isoform_collapse_out",
                        help="Output directory.")
    parser.add_argument("--tmpdir", default=None,
                        help="Temporary directory for sort runs.")
    parser.add_argument("--run-size", type=int, default=500000,
                        help="Records per sorted run (Python sort/merge).")
    parser.add_argument("--low-evidence-per-sample", type=float, default=2.0,
                        help="Low-evidence threshold as reads per sample (default 2.0).")
    parser.add_argument("--min-similarity", type=float, default=0.85,
                        help="Minimum similarity to collapse low-evidence isoforms (default 0.85).")
    parser.add_argument("--cluster-strategy", choices=["substring", "jaccard"], default="substring",
                        help="Similarity strategy for low-evidence collapsing.")
    parser.add_argument("--keep-raw-sorted", action="store_true",
                        help="Keep sorted raw records.")
    parser.add_argument("--keep-runs", action="store_true",
                        help="Keep sorted run files.")
    args = parser.parse_args()

    inputs = _iter_input_files(args.inputs, args.input_dir)
    if not inputs:
        print("No transcript_associations.txt files found.")
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(args.tmpdir) if args.tmpdir else outdir / "tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)

    raw_path = outdir / "raw_transcript_associations.tsv"
    sorted_raw = outdir / "raw_transcript_associations.sorted.tsv"
    runs_dir = tmpdir / "runs"

    print(f"Parsing {len(inputs)} input files...")
    parse_and_write_raw(inputs, raw_path)

    print("Sorting raw records...")
    run_paths = _write_sorted_runs(raw_path, runs_dir, args.run_size)
    _merge_sorted_runs(run_paths, sorted_raw)

    counts_path = outdir / "structure_sample_counts.tsv"
    totals_path = outdir / "structure_totals.tsv"
    gene_counts_path = outdir / "gene_sample_counts.tsv"
    collapse_map_path = outdir / "structure_collapse_map.tsv"
    collapsed_totals_path = outdir / "collapsed_structure_totals.tsv"
    collapsed_sample_counts_path = outdir / "collapsed_structure_sample_counts.tsv"
    collapsed_associations_path = outdir / "collapsed_transcript_associations.tsv"

    print("Aggregating counts per structure/sample...")
    aggregate_structure_sample_counts(sorted_raw, counts_path)
    print("Aggregating totals per structure...")
    aggregate_structure_totals(counts_path, totals_path, gene_counts_path)
    print("Collapsing low-evidence structures...")
    collapse_structures(
        totals_path,
        gene_counts_path,
        collapse_map_path,
        collapsed_totals_path,
        args.low_evidence_per_sample,
        args.min_similarity,
        args.cluster_strategy,
    )
    print("Aggregating collapsed sample counts...")
    collapse_sample_counts(counts_path, collapse_map_path, collapsed_sample_counts_path)
    print("Writing collapsed molecule associations...")
    collapse_associations(sorted_raw, collapse_map_path, collapsed_associations_path)

    if not args.keep_raw_sorted:
        try:
            os.remove(sorted_raw)
        except OSError:
            pass
    if not args.keep_runs:
        shutil.rmtree(runs_dir, ignore_errors=True)
    if not args.keep_raw_sorted:
        try:
            os.remove(raw_path)
        except OSError:
            pass

    print("Done.")


if __name__ == "__main__":
    main()
