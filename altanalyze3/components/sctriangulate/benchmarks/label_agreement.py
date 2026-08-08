"""Compare the label sets produced by the end_to_end runs, with ARI and NMI.

WHAT THERE IS TO VALIDATE AGAINST
---------------------------------
The upstream repository ships no expected results. `test/mini_test.py` contains
no assert statement; it runs the pipeline and checks only that nothing raises.
The demo h5ad carries no ground-truth cell type column, only the three leiden
annotations that scTriangulate is asked to reconcile. So there is no published
answer key, and "correct" can only mean "what upstream's own code produces".

That makes upstream's run-to-run variability the yardstick, which is why the
end_to_end harness runs each implementation twice.

THREE MEASURES, BECAUSE THEY ANSWER DIFFERENT QUESTIONS
-------------------------------------------------------
exact  : fraction of cells given the identical label string. Strictest. Sensitive
         to a cluster being renamed even when the grouping is unchanged.
ARI    : adjusted Rand index. Compares the PARTITIONS and ignores label names,
         adjusted so that a chance agreement scores 0 and an identical partition
         scores 1.
NMI    : normalised mutual information. Also label-invariant, and less punishing
         than ARI when one run splits a cluster the other keeps whole.

Usage:
  python3.11 label_agreement.py
"""

import os
import json
import itertools

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

HERE = os.path.dirname(os.path.abspath(__file__))
COLUMNS = ['final_annotation', 'raw', 'pruned']
RUNS = [
    ('upstream  seq  A', 'out_e2e_upstream_sequential_1000'),
    ('upstream  seq  B', 'out_e2e_upstream_sequential_2000'),
    ('upstream  par   ', 'out_e2e_upstream_parallel_1000'),
    ('optimized seq  A', 'out_e2e_optimized_sequential_1000'),
    ('optimized seq  B', 'out_e2e_optimized_sequential_2000'),
    ('optimized par   ', 'out_e2e_optimized_parallel_1000'),
]

PAIRS = [
    ('upstream reproduces itself', 'upstream  seq  A', 'upstream  seq  B'),
    ('optimized reproduces itself', 'optimized seq  A', 'optimized seq  B'),
    ('optimized seq vs parallel', 'optimized seq  A', 'optimized par   '),
    ('optimized vs upstream', 'optimized seq  A', 'upstream  seq  A'),
    ('optimized vs upstream (other seeds)', 'optimized seq  B', 'upstream  seq  B'),
]


def load():
    out = {}
    for label, d in RUNS:
        path = os.path.join(HERE, d, 'end_to_end_result.json')
        if not os.path.exists(path):
            print('missing: {}'.format(path))
            continue
        out[label] = json.load(open(path))
    return out


def main():
    runs = load()
    if not runs:
        raise SystemExit('no end_to_end results found; run end_to_end.py first')

    print('{:<38s} {:<18s} {:>12s} {:>8s} {:>8s}'.format(
        'COMPARISON', 'COLUMN', 'EXACT', 'ARI', 'NMI'))
    print('-' * 92)
    for title, a_label, b_label in PAIRS:
        if a_label not in runs or b_label not in runs:
            continue
        a, b = runs[a_label], runs[b_label]
        for col in COLUMNS:
            x, y = a[col], b[col]
            n_match = sum(1 for p, q in zip(x, y) if p == q)
            ari = adjusted_rand_score(x, y)
            nmi = normalized_mutual_info_score(x, y)
            print('{:<38s} {:<18s} {:>6d}/{:<5d} {:>8.4f} {:>8.4f}'.format(
                title if col == COLUMNS[0] else '', col, n_match, len(x), ari, nmi))
        print()

    print('Cluster counts per run:')
    for label, _ in RUNS:
        if label not in runs:
            continue
        r = runs[label]
        print('  {}  {}'.format(label, {c: len(set(r[c])) for c in COLUMNS}))

    rows = []
    for title, a_label, b_label in PAIRS:
        if a_label not in runs or b_label not in runs:
            continue
        a, b = runs[a_label], runs[b_label]
        rows.append({'comparison': title, 'a': a_label.strip(), 'b': b_label.strip(),
                     **{col: {'exact': sum(1 for p, q in zip(a[col], b[col]) if p == q),
                              'n': len(a[col]),
                              'ari': adjusted_rand_score(a[col], b[col]),
                              'nmi': normalized_mutual_info_score(a[col], b[col])}
                        for col in COLUMNS}})
    out = os.path.join(HERE, 'label_agreement.json')
    with open(out, 'w') as f:
        json.dump(rows, f, indent=1)
    print('\nwrote {}'.format(out))


if __name__ == '__main__':
    main()
