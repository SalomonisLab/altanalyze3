"""Rebuild and audit SLC7A5, MPZL1, and SLC24A4 using actual cohort counts.

Outputs are separate from the existing viewer/reference bundles. A top hypothesis
is not a validated isoform; preserve evidence, alternatives and NMD annotations.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import pickle
import sys

os.environ.setdefault('SNAF_OFFLINE', '1')
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np
import pysam
from altanalyze3.components.snaf.surface import evidence_isoform as E
from altanalyze3.components.snaf.surface import predict_isoform as P
from benchmark_evidence_isoforms import EXAMPLES, write_tsv


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', required=True)
    ap.add_argument('--evidence', required=True)
    ap.add_argument('--ogawa-evidence')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--ranker')
    ap.add_argument('--exon-annotation')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    models = pickle.load(open(Path(args.project) / 'results/index/ens91_models.pkl', 'rb'))
    annotation = E.load_exon_annotation(args.exon_annotation, set(models)) if args.exon_annotation else {}
    samples, counts = pickle.load(open(args.evidence, 'rb'))
    ranker = pickle.load(open(args.ranker, 'rb')) if args.ranker else None
    fa = pysam.FastaFile(args.genome)
    predictions, reports, gtf, mrna, protein = [], [], [], [], []
    for uid, (symbol, chrom, target) in EXAMPLES.items():
        gene = uid.split(':')[0]
        sample_names, gene_counts = samples, counts[gene]
        cohort = 'pedAML'
        if symbol == 'SLC24A4' and args.ogawa_evidence:
            sample_names, ogawa = pickle.load(open(args.ogawa_evidence, 'rb'))
            gene_counts = ogawa[gene]
            cohort = 'MDS_Ogawa'
        gm = models[gene]
        fetch = lambda s, e: fa.fetch(chrom.removeprefix('chr'), s - 1, e)
        ref_introns = {j for _, ch in gm.values() for j in E.introns(ch)}
        hosts = annotation.get(gene, {}).get('introns', ref_introns)
        ev = E.Evidence({j: v['counts'] for j, v in gene_counts.items()})
        partners = []
        for j in gene_counts:
            if j == target:
                continue
            if E.paired_exon(target, j, hosts, ev):
                left, right = sorted([target, j])
                partners.append({'uid': gene_counts[j]['uid'], 'coordinates': j,
                                 'exon_nt': right[0] - left[1] + 1,
                                 'both_gt0': E.Evidence(ev.counts, min_reads=1e-12).recurrence([target, j]),
                                 'both_ge3': ev.recurrence([target, j]),
                                 'both_ge10': E.Evidence(ev.counts, min_reads=10).recurrence([target, j])})
        # Legacy uses nearby candidate junctions from the old published frequency list;
        # recover the ACTUAL selected historical row instead of fabricating a baseline.
        bundle = 'MDS_Ogawa' if symbol == 'SLC24A4' else 'pedAML'
        old_path = Path(args.project) / ('results/%s_synthetic/%s_predicted_isoforms.tsv' % (bundle, bundle))
        csv.field_size_limit(10**9)
        historical = []
        if old_path.exists():
            with open(old_path) as fh:
                historical = [r for r in csv.DictReader(fh, delimiter='\t') if r['junction_id'] == uid]
        # These exemplar genes have no supervised labels in the 201-event benchmark;
        # fail loudly if a future benchmark starts including one of them.
        model = None
        if ranker:
            if gene in ranker['training_genes']:
                raise ValueError('Exemplar gene was in training: ' + gene)
            model = ranker['model']
        alternatives = []
        for edits in [1, 2, 4]:
            candidates = E.generate_candidates(gm, target, fetch, ev, max_edits=edits,
                junction_label=uid, labels={j: v['uid'] for j, v in gene_counts.items()},
                exon_annotation=annotation.get(gene))
            ranked = E.rank_candidates(candidates, model)
            alternatives.append({'max_edits': edits, 'candidates': len(candidates),
                                 'top_backbone': ranked[0][1].prediction.backbone if ranked else None,
                                 'top_chain': ranked[0][1].prediction.chain if ranked else [],
                                 'top_protein_nt': 3 * len(ranked[0][1].prediction.protein) if ranked else 0})
        candidates = E.generate_hypotheses(gm, target, fetch, ev,
            junction_label=uid, labels={j: v['uid'] for j, v in gene_counts.items()},
            exon_annotation=annotation.get(gene))
        ranked = E.rank_candidates(candidates, model)
        for rank, (score, c) in enumerate(ranked[:10], 1):
            p = c.prediction
            # Learned NMD is retained as a hypothesis annotation, never renamed non-NMD.
            row = dict(symbol=symbol, cohort=cohort, rank=rank, score=score,
                       construction=c.construction, joint_sample_columns=len(c.joint_samples),
                       sample_names=','.join(sample_names[i] for i in c.joint_samples),
                       **P.prediction_row(p, uid))
            predictions.append(row)
            if rank == 1:
                pid = symbol + '_evidence_top1'
                for s, e in p.chain:
                    gtf.append('\t'.join([chrom, 'SNAF_evidence', 'exon', str(s), str(e), '.',
                                          p.strand, '.', 'gene_id "%s"; transcript_id "%s";' % (gene, pid)]))
                mrna.append('>' + pid + '\n' + p.mrna)
                protein.append('>' + pid + '\n' + p.protein)
        report = dict(symbol=symbol, junction=uid, cohort=cohort, sample_columns=len(sample_names),
                      target_gt0=E.Evidence(ev.counts, min_reads=1e-12).recurrence([target]),
                      target_ge3=ev.recurrence([target]), paired_exon_options=partners,
                      iterations=alternatives,
                      historical=[{k: r.get(k) for k in ['predicted_isoform_id', 'exon_chain',
                                                       'protein_length', 'nmd']} for r in historical])
        reports.append(report)
        print(symbol, 'target samples', report['target_ge3'], 'partner options', len(partners),
              'candidates', len(ranked), flush=True)
    write_tsv(out / 'exemplar_alternatives.tsv', predictions)
    (out / 'exemplars.json').write_text(json.dumps(reports, indent=2))
    (out / 'exemplar_top1.gtf').write_text('\n'.join(gtf) + '\n')
    (out / 'exemplar_top1.mrna.fa').write_text('\n'.join(mrna) + '\n')
    (out / 'exemplar_top1.protein.fa').write_text('\n'.join(protein) + '\n')


if __name__ == '__main__':
    main()
