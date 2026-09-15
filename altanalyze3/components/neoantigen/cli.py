"""Dependency-light command-line interface for portable SNAF workflow stages."""
import argparse
import csv
import math
from pathlib import Path
from .io import manifest


def filter_counts(source, output, min_reads=20):
    if min_reads < 0 or not math.isfinite(min_reads):
        raise ValueError('min_reads must be finite and nonnegative')
    if Path(source).resolve() == Path(output).resolve() or (Path(output).exists() and Path(source).samefile(output)):
        raise ValueError('Input and output count paths must differ')
    kept = total = 0
    with Path(source).open(newline='') as inp, Path(output).open('w', newline='') as out:
        reader, writer = csv.reader(inp, delimiter='\t'), csv.writer(out, delimiter='\t')
        header = next(reader)
        if len(header) < 2 or len(set(header)) != len(header):
            raise ValueError('Count matrix needs distinct sample columns')
        writer.writerow(header)
        for row in reader:
            total += 1
            if len(row) != len(header):
                raise ValueError(f'Ragged count row {total}')
            values = [float(v) for v in row[1:]]
            if any(not math.isfinite(v) or v < 0 for v in values):
                raise ValueError(f'Invalid count in row {total}')
            if max(values) >= min_reads:
                kept += 1; writer.writerow(row)
    manifest(str(output)+'.json', min_reads=min_reads, input_rows=total, retained_rows=kept)
    return kept


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    export = commands.add_parser('export')
    export.add_argument('--candidates', nargs='+', required=True)
    export.add_argument('--outdir', required=True)
    export.add_argument('--canonical-fasta')
    export.add_argument('--predictor', default='')
    merge = commands.add_parser('merge')
    merge.add_argument('--candidates', required=True)
    merge.add_argument('--evidence', required=True)
    merge.add_argument('--output', required=True)
    merge.add_argument('--q-threshold', type=float, default=.01)
    run = commands.add_parser('proteomics')
    run.add_argument('--bundle', required=True)
    run.add_argument('--psm-table', required=True)
    run.add_argument('--sample', required=True)
    run.add_argument('--outdir', required=True)
    run.add_argument('--executable', default='pyneoquant')
    run.add_argument('--search-format', default='generic')
    run.add_argument('--q-threshold', type=float, default=.01)
    hla = commands.add_parser('hla')
    hla.add_argument('--sample', required=True)
    hla.add_argument('--bam')
    hla.add_argument('--supplied')
    hla.add_argument('--output', required=True)
    hla.add_argument('--qc', required=True)
    hla.add_argument('--build', default='auto', choices=['auto', 'hg19', 'hg38'])
    hla.add_argument('--min-depth', type=int, default=8)
    hla.add_argument('--require-all', action='store_true')
    combine = commands.add_parser('combine-hla')
    combine.add_argument('--inputs', nargs='+', required=True)
    combine.add_argument('--output', required=True)
    binding = commands.add_parser('bind')
    binding.add_argument('--evidence', required=True)
    binding.add_argument('--hla', required=True)
    binding.add_argument('--output', required=True)
    binding.add_argument('--method', choices=['MHCflurry', 'netMHCpan'], default='MHCflurry')
    binding.add_argument('--software-path')
    filt = commands.add_parser('filter')
    filt.add_argument('--counts', required=True)
    filt.add_argument('--output', required=True)
    filt.add_argument('--min-reads', type=float, default=20)
    unpack = commands.add_parser('unpack')
    unpack.add_argument('--archive', required=True)
    unpack.add_argument('--outdir', required=True)
    args = parser.parse_args(argv)
    if args.command == 'unpack':
        from .archive import unpack_reference
        unpack_reference(args.archive, args.outdir)
    elif args.command == 'export':
        from .proteomics import export_candidates
        export_candidates(args.candidates, args.outdir, args.canonical_fasta, args.predictor)
    elif args.command == 'merge':
        from .proteomics import merge_evidence
        merge_evidence(args.candidates, args.evidence, args.output, args.q_threshold)
    elif args.command == 'proteomics':
        from .proteomics import run_pyneoquant
        run_pyneoquant(args.bundle, args.psm_table, args.sample, args.outdir, args.executable, args.search_format, args.q_threshold)
    elif args.command == 'hla':
        from .hla import prepare_hla
        prepare_hla(args.sample, args.output, args.qc, args.bam, args.supplied, args.build, args.min_depth, args.require_all)
    elif args.command == 'combine-hla':
        from .hla import combine_hla
        combine_hla(args.inputs, args.output)
    elif args.command == 'bind':
        from .hla import predict_binding
        predict_binding(args.evidence, args.hla, args.output, args.method, args.software_path)
    elif args.command == 'filter':
        filter_counts(args.counts, args.output, args.min_reads)

if __name__ == '__main__':
    main()
