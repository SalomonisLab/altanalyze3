"""Render the completed isoform experiment as a report and standalone figures."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def table(records, columns):
    lines = ['| ' + ' | '.join(label for _, label in columns) + ' |',
             '| ' + ' | '.join('---' for _ in columns) + ' |']
    for r in records:
        values = []
        for k, _ in columns:
            v = r[k]
            if k in ['protein_identity', 'intron_f1']:
                v = '%.3f' % float(v)
            elif k in ['exact_chain', 'epitope', 'any_lr_exact', 'exact_known_chain']:
                v = '%d/%s (%.1f%%)' % (round(float(v) * int(r['n'])), r['n'], 100 * float(v))
            values.append(str(v))
        lines.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', required=True)
    args = ap.parse_args()
    root = Path(args.root)
    ev = root / 'evaluation'
    summaries = read(ev / 'summary.tsv')
    wanted = ['legacy_pair', 'local2_evidence', 'learned_pairwise', 'nested_learned', 'oracle_joint']
    names = {'legacy_pair': 'Legacy: one target + one candidate-list partner',
             'local2_evidence': 'Measured pairs + explicit evidence rule',
             'learned_pairwise': 'Measured combinations + linear ranker (fixed joint objective)',
             'nested_learned': 'Measured combinations + nested model selection',
             'oracle_joint': 'Oracle over generated candidates (not deployable)'}
    cols = [('method', 'Method'), ('n', 'N'), ('protein_identity', 'Protein identity'),
            ('intron_f1', 'Intron F1'), ('exact_chain', 'Exact intron chain'),
            ('epitope', 'Reported peptide recovered')]
    def selected(subset):
        return [dict(next(r for r in summaries if r['method'] == m and r['subset'] == subset),
                     method=names[m]) for m in wanted]
    known = read(root / 'known_evaluation/known_summary.tsv')
    examples = json.loads((root / 'exemplars_final/exemplars.json').read_text())
    example_rows = read(root / 'exemplars_final/exemplar_alternatives.tsv')
    top = {r['symbol']: r for r in example_rows if r['rank'] == '1'}
    audit = read(ev / 'truth_audit.tsv')
    metrics = read(ev / 'metrics_per_junction.tsv')
    verified = [r for r in metrics if r['truth_verified'] == '1' and r['truth_source'] == 'longread']
    base = {r['uid']: r for r in verified if r['method'] == 'legacy_pair'}
    nested = {r['uid']: r for r in verified if r['method'] == 'nested_learned'}
    genes = sorted({r['gene'] for r in base.values()})
    intervals = {}
    rng = np.random.default_rng(193)
    for metric in ['protein_identity', 'intron_f1', 'exact_chain', 'epitope']:
        by_gene = {g: [float(nested[u][metric]) - float(r[metric])
                       for u, r in base.items() if r['gene'] == g] for g in genes}
        boot = [np.mean([v for g in rng.choice(genes, len(genes), replace=True) for v in by_gene[g]])
                for _ in range(2000)]
        intervals[metric] = {'delta': float(np.mean([v for vv in by_gene.values() for v in vv])),
                             'ci95': np.quantile(boot, [.025, .975]).tolist()}
    (ev / 'uncertainty_verified_longread.json').write_text(json.dumps(intervals, indent=2))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))
    methods = ['legacy_pair', 'local2_evidence', 'learned_pairwise', 'nested_learned']
    for ax, (metric, title) in zip(axes, [('protein_identity', 'Protein identity'),
                                         ('intron_f1', 'Intron-chain F1'),
                                         ('exact_chain', 'Exact intron chain'),
                                         ('epitope', 'Peptide recovery')]):
        vals = [float(next(r for r in summaries if r['method'] == m and r['subset'] == 'verified_longread')[metric])
                for m in methods]
        ax.bar(range(4), vals, color=['#98a0a8', '#55a5a5', '#3c6b99', '#253951'])
        ax.set_xticks(range(4), ['Legacy', 'Evidence', 'Linear', 'Nested'], rotation=35, ha='right')
        ax.set_ylim(0, 1)
        ax.set_title(title)
        for i, v in enumerate(vals):
            ax.text(i, v + .02, '%.3f' % v, ha='center', fontsize=9)
    fig.suptitle('191 verified long-read cases · held-out genes · unresolved predictions count as zero')
    fig.tight_layout()
    fig.savefig(root / 'benchmark.png', dpi=180)
    fig.savefig(root / 'benchmark.svg')
    plt.close(fig)
    regions = {'SLC7A5': (87840000, 87853000, 87843643, 87843819),
               'MPZL1': (167721500, 167766000, 167744903, 167744987),
               'SLC24A4': (92322500, 92326500, 92324060, 92324200)}
    fig, axes = plt.subplots(3, 1, figsize=(12, 7.3))
    for ax, d in zip(axes, examples):
        sym = d['symbol']
        lo, hi, ns, ne = regions[sym]
        old = d['historical'][0]['exon_chain'] if d['historical'] else ''
        new = top[sym]['exon_chain']
        for y, text, color in [(1, old, '#d69256'), (0, new, '#3c6b99')]:
            exons = [tuple(map(int, x.split('-'))) for x in text.split(',') if x]
            ax.plot([(max(lo, exons[0][0]) - lo) / 1000, (min(hi, exons[-1][1]) - lo) / 1000],
                    [y, y], color=color, linewidth=1)
            for s, e in exons:
                a, b = max(lo, s), min(hi, e)
                if a <= b:
                    ax.add_patch(Rectangle(((a - lo) / 1000, y - .14), max((b - a + 1) / 1000, .015),
                                           .28, color='#bf4852' if (s, e) == (ns, ne) and y == 0 else color))
        ax.annotate('%d-base exon' % (ne - ns + 1), xy=((ns - lo) / 1000, .1),
                    xytext=((ns - lo) / 1000, .5), ha='center',
                    arrowprops={'arrowstyle': '->', 'color': '#bf4852'}, color='#a13741')
        ax.set_xlim(0, (hi - lo) / 1000)
        ax.set_ylim(-.35, 1.45)
        ax.set_yticks([0, 1], ['New hypothesis', 'Historical prediction'])
        ax.set_xlabel('Genomic distance from %s (kb; local window)' % format(lo, ','))
        ax.set_title('%s · %s · %s jointly supporting sample columns at ≥3 reads/arm' %
                     (sym, d['cohort'], top[sym]['joint_sample_columns']), loc='left', fontweight='bold')
    fig.tight_layout()
    fig.savefig(root / 'exemplars.png', dpi=180)
    fig.savefig(root / 'exemplars.svg')
    plt.close(fig)
    text = ['# SNAF-B synthetic isoform inference — September 15, 2026',
            '\n## Outcome\n',
            'Implemented an inference method that actually uses measured sample-level junction combinations, '
            'plus a small learned candidate ranker. The final construction checks repair both extension-based '
            'and inserted-block exon fusion artifacts. The fitted model reads no long-read catalog at inference. '
            'The comparison includes known reference and long-read-observed isoforms, rather than treating every '
            'historical row as novel long-read truth.',
            '\n## Verified long-read comparison\n', table(selected('verified_longread'), cols),
            '\n![Benchmark](benchmark.png)\n',
            'Protein identity is identical aligned residues divided by the longer full-protein length. '
            'Intron F1 compares the complete predicted and associated splice-gap sets. Exact chain means '
            'the entire intron chain, not exact transcript ends. Peptide recovery checks the reported insert '
            'sequence inside the predicted protein; it does not establish surface presentation. All missing '
            'predictions receive zero. NMD-positive predictions remain labeled as such.',
            '\n### Uncertainty\n',
            'Paired gene-bootstrap differences for nested selection versus legacy (2,000 resamples):',
            '\n| Metric | Difference | 95% interval |\n|---|---:|---:|']
    for k, d in intervals.items():
        text.append('| %s | %+.3f | [%+.3f, %+.3f] |' % (k, d['delta'], *d['ci95']))
    text += ['\nThese are exploratory internal results. Gene-held-out training prevents direct same-gene '
             'memorization, but construction rules and model alternatives were iterated on this collection. '
             'An independent cohort remains necessary to measure generalization after method selection.',
             '\n## Historical 201-row denominator\n', table(selected('all'), cols),
             '\n## Seven known-reference associations\n', table(selected('reference_annotation'), cols),
             '\nThese seven rows name ENST reference transcripts, not KINNEX transcripts. Their chains are '
             'recovered from the reference annotation. The original harness stored their chains as empty, '
             'which incorrectly made them appear to be novel isoforms and distorted structural scores.',
             '\n## Separate known-isoform control\n',
             table(known, [('method', 'Method'), ('n', 'Genes'), ('intron_f1', 'Intron F1'),
                           ('exact_known_chain', 'Exact known chain')]),
             '\nThe 77 control genes have reference intron chains independently present in the long-read catalog. '
             'One well-detected annotated junction is selected per gene before examining predictions; any '
             'annotated, long-read-observed chain containing it is accepted. The trained model excludes that gene. '
             'This control evaluates structure only because the chain index lacks full transcript ends/proteins. '
             'The explicit evidence rule is stronger here than a model trained mainly on novel events.',
             '\n## Three exemplars\n', '\n![Exemplar exon structures](exemplars.png)\n']
    for d in examples:
        sym = d['symbol']; r = top[sym]; lo, hi, ns, ne = regions[sym]
        text += ['\n### %s\n' % sym,
                 '`%s` (%s).' % (d['junction'], d['cohort']),
                 '\nTop model: **%d-base exon at %s–%s**, supported by **%s sample columns** with '
                 'at least three reads at every applied junction; %s amino acids, predicted NMD=%s. '
                 'Backbone: `%s`.' % (ne - ns + 1, format(ns, ','), format(ne, ','),
                     r['joint_sample_columns'], r['protein_length'], r['nmd'], r['parent_isoform_id']),
                 '\nApplied junctions: `%s`.' % r['junctions_applied']]
    text += ['\nSLC7A5 retains the supported 177-base cassette exon; the previously reported 75 co-detected '
             'columns use >0 reads, while 68 meet the new ≥3-read threshold. MPZL1 uses a measured '
             'E3.7-to-I6.1 partner to make an 85-base exon: the intervening sequence is spliced out instead '
             'of becoming one huge exon. SLC24A4 uses the 141-base cassette exon in 31 Ogawa columns '
             '(76 at >0 reads). Its top reference backbone has an alternative first exon; this is an '
             'annotated transcript-start hypothesis, not evidence that short reads resolved the full 5′ end. '
             'An alternative backbone preserving E2/E3 is retained in the ranked output.',
             '\nAll three exemplar genes are excluded from ranker training. Sample co-detection is not '
             'molecular phasing, and sample columns are not necessarily independent patients.',
             '\n## Method and experiments\n',
             'The generator edits every available backbone, searches compatible junction combinations observed '
             'jointly in samples, and preserves candidates from smaller searches when the larger beam is capped. '
             'It distinguishes exon/intron coverage markers from splice gaps, checks every junction after editing, '
             'and requires complementary novel-exon boundaries in the same local intron. Broader AltAnalyze exon '
             'annotation is essential: MPZL1 E4/E5/E6 have historical mRNA accessions absent from the Ensembl '
             'transcript models. New blocks cannot silently absorb these introns.',
             '\nDefaults: ≥3 reads; ≤4 edits; 16 neighbors; beam 24; 160 chains per search depth; '
             'three distinct longest ORFs plus a non-NMD alternative. Unbounded extensions are limited to '
             '500 bases. Novel exons >500 bases need ≥3 jointly supporting sample columns. These bounds '
             'are explicit priors and search limits; unresolved events are not declared biologically invalid.',
             '\nThe linear model learns within-junction differences in 20 candidate features. Ridge, '
             'extra trees, boosting, and linear ranking were compared with five gene-held-out outer folds and '
             'an inner gene split. The objectives were structural F1 or mean structural F1/protein identity. '
             'The fixed linear model is the deployable model family selected most often internally. '
             'Nested-selection estimates show the cost of choosing a model without its outer-test labels.',
             '\nEarlier permissive iterations scored higher on protein similarity but admitted invalid MPZL1 '
             'models. Those results are retained in the numbered/experimental folders for inspection, but '
             '**evaluation/** is the final corrected benchmark. The explicit evidence rule leads on exact '
             'chain recovery in several comparisons; a learned score is not a universal improvement.',
             '\n## Truth audit\n',
             '\n| Junction | Audit result |\n|---|---|']
    for r in audit:
        if not r['status'].startswith('verified_longread'):
            text.append('| `%s` | %s |' % (r['junction_id'], r['status']))
    text += ['\nFusion cases remain zero in the historical denominator and are outside this single-locus '
             'generator. The chromosome-mismatched association is excluded from supervised training and the '
             'verified comparison. Reference cases are scored against recovered annotation, with provenance '
             'separate from observed long-read truth.',
             '\n## Deliverables and validation\n',
             '- [Per-event results](evaluation/metrics_per_junction.tsv), [summary](evaluation/summary.tsv), '
             '[truth audit](evaluation/truth_audit.tsv), and [held-out predictions](evaluation/heldout_predictions.tsv).',
             '- [Fitted ranker](evaluation/ranker.pkl), [fold definitions](evaluation/folds.json), '
             '[manifest](evaluation/manifest.json).',
             '- [Exemplar alternatives](exemplars_final/exemplar_alternatives.tsv), '
             '[exemplar GTF](exemplars_final/exemplar_top1.gtf), '
             '[mRNA FASTA](exemplars_final/exemplar_top1.mrna.fa), '
             '[protein FASTA](exemplars_final/exemplar_top1.protein.fa).',
             '- [Known controls](known_evaluation/known_summary.tsv).',
             '- [Test log](tests.log): full SNAF suite, 65 passing tests; four existing dependency deprecation warnings.',
             '\nThe repository exposes `predict_isoform(..., evidence=..., ranking_model=...)` and a '
             '`surface.predict_evidence` CLI. Use the new evidence argument or CLI explicitly; existing '
             'reference-only calls retain their legacy behavior. See the repository guide '
             '`components/snaf/docs/SYNTHETIC_ISOFORM_INFERENCE.md` and [reproduction commands](reproduce.sh).',
             '\n## Remaining limits\n',
             'Whole-chain recovery remains modest. Transcript ends, distant-event phasing, translation starts, '
             'genuine intron retention, and unseen exon boundaries remain ambiguous. The oracle measures only '
             'what this bounded candidate set can represent; it is not a predictor. Learning from a larger '
             'set of matched long/short-read samples and independent patients is the next meaningful test. '
             'The three rebuilt models are hypotheses with measured splice support, not experimentally '
             'confirmed full-length isoforms.',
             '\n## Background sources\n',
             'Splice-graph transcript reconstruction: [StringTie](https://www.nature.com/articles/nbt.3122). '
             'NMD predictions are annotations rather than measured fate: '
             '[experimental NMD boundary exceptions](https://pmc.ncbi.nlm.nih.gov/articles/PMC1084009/). '
             'All numerical results above come from the local experiments, not these papers.']
    (root / 'REPORT.md').write_text('\n'.join(text) + '\n')
    print(root / 'REPORT.md')


if __name__ == '__main__':
    main()
