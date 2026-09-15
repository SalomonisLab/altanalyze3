"""Known-isoform control: annotated intron chains also observed in long reads.

Select one well-detected annotated junction per gene, without inspecting predictions.
Accept any annotated, long-read-observed chain containing that junction. Evaluate
intron structure only: the long-read chain index does not retain transcript ends or
proteins, so a full-length protein accuracy claim would be unsupported here.
"""
import argparse
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
from benchmark_evidence_isoforms import fold, rows, structure_scores, write_tsv


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', required=True)
    ap.add_argument('--evidence', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--rankers', required=True)
    ap.add_argument('--exon-annotation', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    root = Path(args.project) / 'results'
    models = pickle.load(open(root / 'index/ens91_models.pkl', 'rb'))
    lr = pickle.load(open(root / 'index/longread_chains.pkl', 'rb'))
    _, counts = pickle.load(open(args.evidence, 'rb'))
    genes = sorted({r['gene'] for r in rows(root / '03_truth.tsv')})
    annotations = E.load_exon_annotation(args.exon_annotation, set(genes))
    rankers = {f: pickle.load(open(Path(args.rankers) / ('fold%d_ranker.pkl' % f), 'rb'))
               for f in range(5)}
    fa = pysam.FastaFile(args.genome)
    records = []
    for gene in genes:
        gm = models.get(gene, {})
        known = {E.introns(ch) for _, ch in gm.values() if E.introns(ch) in lr.get(gene, {})}
        ev = E.Evidence({j: v['counts'] for j, v in counts.get(gene, {}).items()})
        eligible = {j for ch in known for j in ch if ev.recurrence([j]) >= 3}
        if not eligible:
            continue
        # The junction is chosen by read evidence alone, then all matching truths accepted.
        target = max(sorted(eligible), key=lambda j: ev.recurrence([j]))
        truth = {ch for ch in known if target in ch}
        chrom = counts[gene][target]['chrom'].removeprefix('chr')
        fetch = lambda s, e: fa.fetch(chrom, s - 1, e)
        baseline = P.predict_isoform(gm, *target, fetch)
        picks = {'legacy_single': baseline}
        for edits in (1, 2, 4):
            cs = E.generate_candidates(gm, target, fetch, ev, max_edits=edits,
                                       exon_annotation=annotations.get(gene))
            ranked = E.rank_candidates(cs)
            picks['local%d_evidence' % edits] = ranked[0][1].prediction if ranked else None
        cs = E.generate_hypotheses(gm, target, fetch, ev, exon_annotation=annotations.get(gene))
        bundle = rankers[fold(gene)]
        assert gene not in bundle['training_genes'], 'Known-control gene leaked into training'
        ranked = E.rank_candidates(cs, bundle['model'])
        picks['nested_learned'] = ranked[0][1].prediction if ranked else None
        for method, p in picks.items():
            pi = E.introns(p.chain) if p else ()
            f1 = max((structure_scores(pi, t)[0] for t in truth), default=0) if p else 0
            records.append(dict(gene=gene, junction='%d-%d' % target, method=method,
                                eligible_truth_chains=len(truth), intron_f1=f1,
                                exact_known_chain=int(p is not None and pi in truth),
                                any_lr_exact=int(p is not None and pi in lr.get(gene, {})),
                                prediction=int(p is not None)))
    write_tsv(out / 'known_per_gene.tsv', records)
    summaries = []
    for method in sorted({r['method'] for r in records}):
        rr = [r for r in records if r['method'] == method]
        summaries.append(dict(method=method, n=len(rr), **{
            k: float(np.mean([r[k] for r in rr])) for k in
            ['intron_f1', 'exact_known_chain', 'any_lr_exact', 'prediction']}))
    write_tsv(out / 'known_summary.tsv', summaries)
    (out / 'known_manifest.json').write_text(json.dumps(vars(args), indent=2))
    print(json.dumps(summaries, indent=2), flush=True)


if __name__ == '__main__':
    main()
