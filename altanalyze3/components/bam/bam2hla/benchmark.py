#!/usr/bin/env python3
"""
Benchmark bam2hla against OptiType calls.

OptiType truth: <sample>.bed \t HLA-A*..,HLA-A*..,HLA-B*..,..,HLA-C*..,HLA-C*..

Two phases (so the slow, network-bound part runs once):
  cache : one pileup pass per BAM -> results/cache/<sample>.<build>.json.gz
          (per-sample SIGALRM timeout; resumable; skips cached)
  eval  : type every cached sample and score 2-field concordance vs OptiType

  python benchmark.py cache --build hg38 --bam-dir <dir> [--n 40] [--timeout 60]
  python benchmark.py eval  --build hg38
"""
import os, sys, json, gzip, glob, signal, argparse, time
from collections import defaultdict

import bam2hla as B

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = '/Volumes/salomonis2/LabFiles/Frank-Li/neoantigen/revision/blood/MDS-Ogawa-2018/aggregated_hla_genotypes.txt'
CACHE = os.path.join(HERE, 'results', 'cache')


def norm(a):
    return a.replace('HLA-', '').strip()


def parse_truth(path=TRUTH):
    truth = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            s, hla = line.rstrip('\n').split('\t')
            sample = s.rsplit('.', 1)[0]
            alleles = [norm(x) for x in hla.split(',')]
            g = defaultdict(list)
            for a in alleles:
                g[a.split('*')[0]].append(a)
            truth[sample] = {f'HLA-{k}': v for k, v in g.items()}
    return truth


class _TO(Exception):
    pass


def _alarm(signum, frame):
    raise _TO()


def cache_sample(bam_path, build, timeout=60):
    """One network pileup pass -> local JSON of per-position base counts."""
    sample = os.path.basename(bam_path).rsplit('.bam', 1)[0]
    out = os.path.join(CACHE, f'{sample}.{build}.json.gz')
    if os.path.exists(out):
        return 'cached'
    sig = B.load_signatures(build)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout)
    try:
        bam = B.pysam.AlignmentFile(bam_path, 'rb')
        chrom, _ = B._chrom_for(bam, B.CHR6_LEN[build])
        data = {'sample': sample, 'build': build, 'genes': {}}
        for g, gs in sig['genes'].items():
            counts = B.pileup_positions(bam, chrom, gs['positions'])
            data['genes'][g] = {str(p): [c['A'], c['C'], c['G'], c['T']]
                                for p, c in counts.items()}
        bam.close()
    except _TO:
        return 'timeout'
    except Exception as e:
        return f'err:{type(e).__name__}:{e}'
    finally:
        signal.alarm(0)
    os.makedirs(CACHE, exist_ok=True)
    with gzip.open(out, 'wt') as fh:
        json.dump(data, fh)
    return 'ok'


def type_from_cache(path):
    with gzip.open(path, 'rt') as fh:
        data = json.load(fh)
    sig = B.load_signatures(data['build'])
    out = {'genes': {}}
    for g, gs in sig['genes'].items():
        raw = data['genes'].get(g, {})
        counts = {int(p): dict(zip(B._BASES, v)) for p, v in raw.items()}
        out['genes'][g] = B.type_gene(gs, counts)
    return data['sample'], data['build'], out


def cmd_cache(args):
    truth = parse_truth()
    os.makedirs(CACHE, exist_ok=True)
    bams = sorted(glob.glob(os.path.join(args.bam_dir, '*.bam')))
    # restrict to samples present in truth
    keep = []
    for b in bams:
        s = os.path.basename(b).rsplit('.bam', 1)[0]
        if s in truth:
            keep.append(b)
    if args.n:
        keep = keep[:args.n]
    print(f'[cache] {len(keep)} BAMs (build={args.build}, timeout={args.timeout}s)')
    tally = defaultdict(int)
    for i, b in enumerate(keep, 1):
        t0 = time.time()
        r = cache_sample(b, args.build, args.timeout)
        tally[r.split(':')[0]] += 1
        print(f'  [{i}/{len(keep)}] {os.path.basename(b):35s} {r:12s} {time.time()-t0:5.1f}s',
              flush=True)
    print(f'[cache] done: {dict(tally)}')


def cmd_eval(args):
    truth = parse_truth()
    files = sorted(glob.glob(os.path.join(CACHE, f'*.{args.build}.json.gz')))
    print(f'[eval] {len(files)} cached samples (build={args.build})')
    genes = ['HLA-A', 'HLA-B', 'HLA-C']
    n_all = n_ok = 0
    per_gene = {g: [0, 0] for g in genes}          # [matched_alleles, total_alleles]
    rows = []
    for f in files:
        sample, build, out = type_from_cache(f)
        tr = truth.get(sample)
        if not tr:
            continue
        line = [sample]
        for g in genes:
            call = out['genes'][g].get('call')
            pred = sorted(call) if call else ['NA', 'NA']
            gold = sorted(tr.get(g, []))
            while len(gold) < 2:
                gold.append('NA')
            # multiset match
            m = _multiset_match(pred, gold)
            per_gene[g][0] += m; per_gene[g][1] += 2
            n_all += 2; n_ok += m
            line.append(f'{g[-1]}:{"/".join(x.split("*")[1] if "*" in x else x for x in pred)}'
                        f' vs {"/".join(x.split("*")[1] if "*" in x else x for x in gold)}'
                        f' [{m}/2]')
        rows.append(line)
    # report
    print('\n'.join('  ' + '  '.join(r) for r in rows))
    print('\n[eval] per-gene 2-field allele concordance:')
    for g in genes:
        m, t = per_gene[g]
        print(f'  {g}: {m}/{t} = {100*m/t:.1f}%' if t else f'  {g}: n/a')
    print(f'[eval] OVERALL: {n_ok}/{n_all} = {100*n_ok/n_all:.1f}%' if n_all else 'no data')


def cmd_run(args):
    """Type BAMs directly (no cache) and score vs OptiType. For fast local
    disk (e.g. the cluster).  Writes a per-sample TSV of calls + concordance."""
    truth = parse_truth(args.truth) if args.truth else parse_truth()
    bams = sorted(glob.glob(os.path.join(args.bam_dir, '*.bam')))
    genes = ['HLA-A', 'HLA-B', 'HLA-C']
    n_all = n_ok = 0
    per_gene = {g: [0, 0] for g in genes}
    outtsv = args.out or os.path.join(HERE, 'results', f'benchmark_{args.build}.tsv')
    os.makedirs(os.path.dirname(outtsv), exist_ok=True)
    fh = open(outtsv, 'w')
    fh.write('sample\t' + '\t'.join(f'{g}_pred\t{g}_optitype\t{g}_match' for g in genes) + '\n')
    for i, b in enumerate(bams, 1):
        sample = os.path.basename(b).rsplit('.bam', 1)[0]
        tr = truth.get(sample)
        if tr is None and not args.all:
            continue
        try:
            out = B.type_bam(b, build=args.build, verbose=False)
        except Exception as e:
            print(f'  [{i}] {sample}: ERROR {e}', flush=True); continue
        row = [sample]
        for g in genes:
            call = out['genes'][g].get('call')
            pred = sorted(call) if call else ['NA', 'NA']
            gold = sorted(tr.get(g, [])) if tr else []
            while len(gold) < 2:
                gold.append('NA')
            m = _multiset_match(pred, gold) if tr else 0
            if tr:
                per_gene[g][0] += m; per_gene[g][1] += 2; n_all += 2; n_ok += m
            row += ['/'.join(pred), '/'.join(gold), str(m)]
        fh.write('\t'.join(row) + '\n'); fh.flush()
        print(f'  [{i}/{len(bams)}] {sample}: ' +
              ' '.join(f'{g[-1]}={"/".join(a.split("*")[-1] for a in sorted(out["genes"][g].get("call") or ["NA","NA"]))}'
                       for g in genes), flush=True)
    fh.close()
    if n_all:
        print('\n[run] per-gene 2-field concordance vs OptiType:')
        for g in genes:
            m, t = per_gene[g]
            print(f'  {g}: {m}/{t} = {100*m/t:.1f}%')
        print(f'[run] OVERALL: {n_ok}/{n_all} = {100*n_ok/n_all:.1f}%')
    print(f'[run] wrote {outtsv}')


def _multiset_match(pred, gold):
    g = list(gold); m = 0
    for a in pred:
        if a in g:
            g.remove(a); m += 1
    return m


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('cache'); c.add_argument('--build', required=True)
    c.add_argument('--bam-dir', required=True); c.add_argument('--n', type=int, default=0)
    c.add_argument('--timeout', type=int, default=60); c.set_defaults(func=cmd_cache)
    e = sub.add_parser('eval'); e.add_argument('--build', default='hg38')
    e.set_defaults(func=cmd_eval)
    r = sub.add_parser('run'); r.add_argument('--build', default='auto')
    r.add_argument('--bam-dir', required=True); r.add_argument('--truth', default=None)
    r.add_argument('--out', default=None); r.add_argument('--all', action='store_true')
    r.set_defaults(func=cmd_run)
    args = ap.parse_args()
    args.func(args)
