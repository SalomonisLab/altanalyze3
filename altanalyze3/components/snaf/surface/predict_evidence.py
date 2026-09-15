"""CLI for short-read-supported SNAF-B synthetic isoform hypotheses.

Run with python -m altanalyze3.components.snaf.surface.predict_evidence --help.
The evidence cache contains (sample_names, {gene: {coordinates: {counts, uid, chrom}}}).
Prepare it with components/snaf/dev/prepare_isoform_evidence.py. Model pickles are
local artifacts from benchmark_evidence_isoforms.py; no long-read input is needed here.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import pickle

from . import evidence_isoform as E
from . import predict_isoform as P


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--junctions', required=True,
                    help='TSV: junction_id, gene, chrom, junction_start, junction_end; optional first_exon (true/false)')
    ap.add_argument('--models', required=True,
                    help='Reference-only pickle: {gene: {transcript: (strand, exon_chain)}}')
    ap.add_argument('--evidence', required=True)
    ap.add_argument('--exon-annotation', required=True, help='AltAnalyze Hs_Ensembl_exon.txt')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--method', choices=['learned', 'evidence'], default='learned',
                    help='learned (default): bundled ranker; evidence: explicit two-edit rule')
    ap.add_argument('--ranker', help='Optional JSON/pickle model overriding the bundled learned ranker')
    ap.add_argument('--sample', help='Restrict novel junction combinations to this sample column')
    ap.add_argument('--min-reads', type=float, default=3)
    ap.add_argument('--max-edits', type=int, default=None)
    ap.add_argument('--search-mode', choices=['union', 'depth'], default=None,
                    help='union preserves shorter searches; depth reproduces a single benchmark search depth')
    ap.add_argument('--max-extension-nt', type=int, default=500)
    ap.add_argument('--top-k', type=int, default=5)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    args.max_edits = args.max_edits if args.max_edits is not None else (4 if args.method == 'learned' else 2)
    args.search_mode = args.search_mode or ('union' if args.method == 'learned' else 'depth')
    if args.ranker and args.method != 'learned':
        ap.error('--ranker requires --method learned')
    if args.top_k < 1 or args.max_edits < 1:
        ap.error('--top-k and --max-edits must be positive')
    import pysam
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.junctions) as fh:
        junctions = list(csv.DictReader(fh, delimiter='\t'))
    models = pickle.load(open(args.models, 'rb'))
    sample_names, counts = pickle.load(open(args.evidence, 'rb'))
    if args.sample and args.sample not in sample_names:
        ap.error('Sample absent from evidence cache: ' + args.sample)
    sample_index = sample_names.index(args.sample) if args.sample else None
    annotation = E.load_exon_annotation(args.exon_annotation, {r['gene'] for r in junctions})
    model = E.load_ranker(args.ranker) if args.method == 'learned' else None
    fa = pysam.FastaFile(args.genome)
    records, gtfs, mrnas, proteins = [], [], [], []
    for r in junctions:
        gene, uid = r['gene'], r['junction_id']
        first_exon = r.get('first_exon', '').strip().lower()
        if first_exon not in ('', 'true', 'false', '1', '0'):
            raise ValueError('first_exon must be true/false or 1/0: ' + uid)
        target = tuple(sorted((int(r['junction_start']), int(r['junction_end']))))
        chrom = r['chrom'] if r['chrom'] in fa.references else r['chrom'].removeprefix('chr')
        if chrom not in fa.references:
            raise ValueError('Chromosome absent from genome: ' + r['chrom'])
        gc = counts.get(gene, {})
        ev = E.Evidence({j: v['counts'] for j, v in gc.items()}, args.min_reads,
                        sample_index=sample_index if gc else None)
        fetch = lambda s, e: fa.fetch(chrom, s - 1, e)
        builder = E.generate_hypotheses if args.search_mode == 'union' else E.generate_candidates
        cs = builder(models.get(gene, {}), target, fetch, ev,
             first_exon=first_exon in ('true', '1'),
             max_edits=args.max_edits, max_extension_nt=args.max_extension_nt,
             exon_annotation=annotation.get(gene), junction_label=uid,
             labels={j: v['uid'] for j, v in gc.items()}) if (gc or not args.sample) else []
        ranked = E.rank_candidates(cs, model)
        if not ranked:
            records.append(dict(junction_id=uid, gene=gene, rank=0, score='',
                status='unresolved', joint_sample_columns=0, supporting_samples='',
                construction='', candidate_count=0, artifact_id='',
                **{k: '' for k in P.PREDICTION_COLUMNS}))
            continue
        for rank, (score, c) in enumerate(ranked[:args.top_k], 1):
            p = c.prediction
            # Multiple ORFs can have identical parent/junction labels. Artifact IDs
            # include strand, chain and CDS bounds to prevent overwriting alternatives.
            signature = (gene, p.strand, p.chain, p.cds_start, p.cds_end)
            pid = 'SNAF_SYN_' + hashlib.sha256(repr(signature).encode()).hexdigest()[:20]
            records.append(dict(junction_id=uid, gene=gene, rank=rank, score=score,
                status='hypothesis' if c.joint_samples else 'reference_only_hypothesis',
                joint_sample_columns=len(c.joint_samples),
                supporting_samples=','.join(sample_names[i] for i in c.joint_samples),
                construction=c.construction, candidate_count=len(cs), artifact_id=pid,
                **P.prediction_row(p, uid)))
            if rank == 1:
                for s, e in p.chain:
                    gtfs.append('\t'.join([r['chrom'], 'SNAF_evidence', 'exon', str(s), str(e),
                        '.', p.strand, '.', 'gene_id "%s"; transcript_id "%s"; first_exon_boundary_source "%s";'
                        % (gene, pid, p.first_exon_boundary_source)]))
                mrnas.append('>' + pid + '\n' + p.mrna)
                proteins.append('>' + pid + '\n' + p.protein)
    with open(out / 'predictions.tsv', 'w') as fh:
        if records:
            writer = csv.DictWriter(fh, fieldnames=list(records[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
    (out / 'top1.gtf').write_text('\n'.join(dict.fromkeys(gtfs)) + '\n')
    (out / 'top1.mrna.fa').write_text('\n'.join(dict.fromkeys(mrnas)) + '\n')
    (out / 'top1.protein.fa').write_text('\n'.join(dict.fromkeys(proteins)) + '\n')
    manifest = vars(args).copy()
    manifest['input_files'] = {p: {'size': Path(p).stat().st_size,
                                 'mtime_ns': Path(p).stat().st_mtime_ns}
                               for p in [args.junctions, args.models, args.evidence,
                                         args.exon_annotation, args.genome] +
                               ([args.ranker] if args.ranker else [])}
    manifest['feature_names'] = E.FEATURE_NAMES
    manifest['co_detection_is_phasing'] = False
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    fa.close()
    print('Wrote', len(records), 'ranked hypotheses to', out)


if __name__ == '__main__':
    main()
