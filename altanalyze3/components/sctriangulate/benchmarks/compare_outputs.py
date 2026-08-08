"""Compare the optimized run against the pristine upstream run, field by field.

Reads baseline_reference.json (upstream 0.13.0) and optimized_reference.json
(AltAnalyze3 integration), both produced from the same input h5ad by
baseline_profile.py / optimized_profile.py.

Reports, for every field, how many entries match out of how many compared.
Exits non-zero if any field is not a complete match.

Usage:
  python3.11 compare_outputs.py out_baseline/baseline_reference.json \
                                out_optimized/optimized_reference.json
"""

import sys
import json

import numpy as np

TOL = 1e-9


def close(a, b, tol=TOL):
    """Equality that treats inf==inf and nan==nan as a match.

    gseapy returns NES = inf when the permutation null has zero spread, and
    abs(inf - inf) is nan, which would otherwise read as a mismatch.
    """
    a, b = float(a), float(b)
    if np.isnan(a) or np.isnan(b):
        return np.isnan(a) and np.isnan(b)
    if np.isinf(a) or np.isinf(b):
        return a == b
    return abs(a - b) <= tol


def finite_absdiff(a, b):
    """abs(a-b), or 0.0 when both sides are the same non-finite value."""
    a, b = float(a), float(b)
    if np.isfinite(a) and np.isfinite(b):
        return abs(a - b)
    return 0.0 if close(a, b) else float('inf')


class Report:
    def __init__(self):
        self.rows = []
        self.failed = 0

    def add(self, name, n_match, n_total, detail=''):
        ok = (n_match == n_total)
        if not ok:
            self.failed += 1
        self.rows.append((name, n_match, n_total, ok, detail))

    def show(self):
        print('{:<52s} {:>18s}  {:<6s} {}'.format('FIELD', 'MATCH/TOTAL', 'OK', 'NOTE'))
        print('-' * 110)
        for name, m, t, ok, detail in self.rows:
            print('{:<52s} {:>8d}/{:<9d}  {:<6s} {}'.format(
                name, m, t, 'yes' if ok else 'NO', detail))
        print('-' * 110)
        print('{} of {} fields matched completely'.format(
            len(self.rows) - self.failed, len(self.rows)))


def cmp_float_dict(a, b):
    """Return (n_match, n_total, max_abs_diff) over the union of keys."""
    keys = sorted(set(a) | set(b))
    n_match, worst = 0, 0.0
    for k in keys:
        if k not in a or k not in b:
            continue
        d = finite_absdiff(a[k], b[k])
        worst = max(worst, d)
        if close(a[k], b[k]):
            n_match += 1
    return n_match, len(keys), worst


def cmp_list_dict(a, b):
    keys = sorted(set(a) | set(b))
    n_match = sum(1 for k in keys if k in a and k in b and list(a[k]) == list(b[k]))
    return n_match, len(keys)


def main():
    base = json.load(open(sys.argv[1]))
    opt = json.load(open(sys.argv[2]))
    rep = Report()

    assert base['_keys'] == opt['_keys'], 'different query keys'
    keys = base['_keys']

    # the score matrix fed to shapley must be identical before shapley runs
    b_data = np.asarray(base['_data_matrix'])
    o_data = np.asarray(opt['_data_matrix'])
    same = int(np.isclose(b_data, o_data, rtol=0, atol=TOL).sum())
    rep.add('metric matrix fed to shapley', same, b_data.size,
            'max abs diff {:.3e}'.format(float(np.abs(b_data - o_data).max())))

    for key in keys:
        b, o = base[key], opt[key]
        for metric in ['reassign', 'tfidf10', 'tfidf5', 'SCCAF', 'doublet']:
            m, t, worst = cmp_float_dict(b[metric], o[metric])
            rep.add('{}: {}'.format(key, metric), m, t, 'max abs diff {:.3e}'.format(worst))

        m, t = cmp_list_dict(b['marker_genes'], o['marker_genes'])
        rep.add('{}: marker gene lists (per cluster)'.format(key), m, t)
        m, t = cmp_list_dict(b['purify'], o['purify'])
        rep.add('{}: purified gene lists'.format(key), m, t)
        m, t = cmp_list_dict(b['exclusive_genes_head'], o['exclusive_genes_head'])
        rep.add('{}: exclusive genes (top 50)'.format(key), m, t)

        # enrichr: -log10 adjusted p per artifact class per cluster
        n_match = n_total = 0
        worst = 0.0
        for c in b['enrichr']:
            for cls in b['enrichr'][c]:
                n_total += 1
                bv, ov = b['enrichr'][c][cls], o['enrichr'][c].get(cls, np.nan)
                worst = max(worst, finite_absdiff(bv, ov))
                if close(bv, ov):
                    n_match += 1
        rep.add('{}: enrichr -log10(adj p)'.format(key), n_match, n_total,
                'max abs diff {:.3e}'.format(worst))

        # gsea: (nes, matched_size) per artifact class per cluster
        n_match = n_total = 0
        worst = 0.0
        for c in b['gsea']:
            for cls in b['gsea'][c]:
                n_total += 1
                bv, ov = b['gsea'][c][cls], o['gsea'][c].get(cls)
                worst = max(worst, finite_absdiff(bv[0], ov[0]))
                if close(bv[0], ov[0]) and int(bv[1]) == int(ov[1]):
                    n_match += 1
        rep.add('{}: gsea (NES, n_hits)'.format(key), n_match, n_total,
                'max abs NES diff {:.3e}'.format(worst))

        for cm in ['confusion_reassign', 'confusion_sccaf']:
            bm, om = np.asarray(b[cm]), np.asarray(o[cm])
            if bm.shape != om.shape:
                rep.add('{}: {}'.format(key, cm), 0, bm.size, 'shape mismatch')
                continue
            rep.add('{}: {}'.format(key, cm), int((bm == om).sum()), bm.size)

    # shapley: scores and the resulting winner, per mode
    for mode in ['shapley_all_or_none', 'shapley', 'rank', 'rank_all_or_none']:
        bv = np.asarray(base[mode]['values'])
        ov = np.asarray(opt[mode]['values'])
        agree = np.isclose(bv, ov, rtol=0, atol=1e-9, equal_nan=True) | (bv == ov)
        finite = np.isfinite(bv) & np.isfinite(ov)
        worst = float(np.abs(bv[finite] - ov[finite]).max()) if finite.any() else 0.0
        rep.add('shapley values [{}]'.format(mode), int(agree.sum()), bv.size,
                'max abs diff {:.3e}'.format(worst))
        bf, of = base[mode]['final'], opt[mode]['final']
        n_match = sum(1 for x, y in zip(bf, of) if x == y)
        rep.add('winning annotation [{}]'.format(mode), n_match, len(bf))

    n_match = sum(1 for x, y in zip(base['raw'], opt['raw']) if x == y)
    rep.add('raw label (annotation@cluster)', n_match, len(base['raw']))

    rep.show()

    bt, ot = sum(base['_timings'].values()), sum(opt['_timings'].values())
    print('\nbaseline  total {:8.2f} s'.format(bt))
    print('optimized total {:8.2f} s'.format(ot))
    print('speedup         {:8.2f} x'.format(bt / ot))

    return 1 if rep.failed else 0


if __name__ == '__main__':
    sys.exit(main())
