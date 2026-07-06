#!/usr/bin/env python3
"""
Ready-to-run local validation on 7 diverse MDS-Ogawa-2018 samples whose OptiType
2-field calls are embedded below (captured from aggregated_hla_genotypes.txt).
Uses direct BAM paths (no directory glob) + a per-sample timeout, so it tolerates
the flaky archive mount.  Run whenever the mount is connected:

  python validate_subset.py --bam-dir <hg38 bams dir> [--build hg38] [--timeout 120]
"""
import os, sys, argparse, signal, time
from collections import defaultdict
import bam2hla as B

# OptiType 2-field truth (from aggregated_hla_genotypes.txt)
TRUTH = {
    '293T_Ctrl_None': 'A*02:01,A*02:01,B*07:02,B*07:02,C*07:02,C*07:02',
    'PV689_BMMNC':    'A*26:08,A*01:01,B*08:01,B*35:03,C*04:01,C*07:01',
    'PV838_BMMNC':    'A*02:01,A*02:01,B*15:01,B*44:02,C*03:04,C*05:01',
    'PV870_CD34':     'A*32:01,A*03:01,B*07:02,B*15:18,C*07:04,C*07:02',
    'PV663_BMMNC':    'A*03:01,A*11:01,B*35:01,B*35:01,C*02:02,C*04:01',
    'PV359_BMMNC':    'A*24:02,A*01:01,B*18:01,B*37:01,C*07:01,C*06:02',
    'PV965_CD34':     'A*02:01,A*02:01,B*35:01,B*27:05,C*01:02,C*04:01',
}
GENES = ['HLA-A', 'HLA-B', 'HLA-C']


class _TO(Exception):
    pass


def _alarm(s, f):
    raise _TO()


def truth_by_gene(s):
    g = defaultdict(list)
    for a in s.split(','):
        g[f'HLA-{a.split("*")[0]}'].append(a)
    return g


def match(pred, gold):
    g = list(gold); m = 0
    for a in pred:
        if a in g:
            g.remove(a); m += 1
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam-dir', required=True)
    ap.add_argument('--build', default='hg38')
    ap.add_argument('--timeout', type=int, default=120)
    a = ap.parse_args()
    signal.signal(signal.SIGALRM, _alarm)
    per_gene = {g: [0, 0] for g in GENES}
    n_ok = n_all = 0
    for sample, tstr in TRUTH.items():
        bam = os.path.join(a.bam_dir, f'{sample}.bam')
        if not os.path.exists(bam):
            print(f'{sample:18s} MISSING {bam}'); continue
        tr = truth_by_gene(tstr)
        signal.alarm(a.timeout); t0 = time.time()
        try:
            out = B.type_bam(bam, build=a.build, verbose=False)
        except _TO:
            print(f'{sample:18s} TIMEOUT'); continue
        except Exception as e:
            print(f'{sample:18s} ERR {e}'); continue
        finally:
            signal.alarm(0)
        parts = []
        for g in GENES:
            call = out['genes'][g].get('call') or ['NA', 'NA']
            pred = sorted(call); gold = sorted(tr[g])
            m = match(pred, gold); per_gene[g][0] += m; per_gene[g][1] += 2
            n_ok += m; n_all += 2
            pp = '/'.join(x.split('*')[1] for x in pred)
            gg = '/'.join(x.split('*')[1] for x in gold)
            parts.append(f'{g[-1]} {pp} vs {gg} [{m}/2]')
        print(f'{sample:18s} ' + '  '.join(parts) + f'   ({time.time()-t0:.0f}s)', flush=True)
    print('\nper-gene 2-field concordance vs OptiType:')
    for g in GENES:
        m, t = per_gene[g]
        if t:
            print(f'  {g}: {m}/{t} = {100*m/t:.0f}%')
    if n_all:
        print(f'OVERALL: {n_ok}/{n_all} = {100*n_ok/n_all:.1f}%')


if __name__ == '__main__':
    main()
