#!/usr/bin/env python3
"""Build the bundled BayesTS score tables SNAF ships with.

BayesTS is expensive: the GTEx run took hours of pyro SVI over ~2.5M junctions, and the
Tabula Sapiens run covers ~2.85M. Neither ever needs to run again -- the per-junction
answer is fixed by the reference cohort, not by the tumor cohort being analysed. This
extracts JUST the final scores into two small tables that SNAF loads directly.

Sources (read-only; neither is modified)
  GTEx            /Users/saljh8/Dropbox/SNAF/MDS/snaf_test_2026-07-04/db/data/controls/
                  GTEx_junction_counts.snaf_stats.tsv.gz
                  built 2026-07-06 by `altanalyze3 snaf-precompute-control`, 2,476,734
                  junctions, 2,629 samples over 51 tissues.
  Tabula Sapiens  /Users/saljh8/Dropbox/Revio/Healthy-ExNeoEpitopes/BayesTS/
                  Run01b_full_results.zip
                  2,860,898 rows / 2,852,713 unique junctions, 762 pseudobulk samples over
                  112 cell types. Confirmed to be Tabula Sapiens by an exact feature match
                  against counts.TabulaSapiens.h5ad (100.0% UID overlap) and by X_mean
                  taking the modal value 1/112.

Output, one row per junction, in the gzipped-TSV shape control_stats.load_control_stats
already reads:

    uid <tab> bayests_sigma <tab> bayests_percentile

sigma is in [0,1]; LOWER means MORE tumor-specific (less expressed across normal tissue).
percentile is sigma's rank in [0,1] within the reference. UIDs are stored WITHOUT the
'=chr...' coordinate suffix, matching the key SNAF sifts on.
"""
import os
import gzip
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
GTEX_SRC = ('/Users/saljh8/Dropbox/SNAF/MDS/snaf_test_2026-07-04/db/data/controls/'
            'GTEx_junction_counts.snaf_stats.tsv.gz')
TS_SRC = ('/Users/saljh8/Dropbox/Revio/Healthy-ExNeoEpitopes/BayesTS/'
          'Run01b_full_results.zip')
TS_MEMBER = 'Run01b_full_results.txt'

GTEX_OUT = os.path.join(HERE, 'GTEx_BayesTS.tsv.gz')
TS_OUT = os.path.join(HERE, 'TabulaSapiens_BayesTS.tsv.gz')

HEADER = 'uid\tbayests_sigma\tbayests_percentile\n'


def _fmt(v):
    """6 significant digits: sigma and percentile are both in [0,1] and are only ever
    compared against a threshold, so this is far finer than any decision they drive, and
    it roughly halves the stored size versus full repr."""
    try:
        return '{:.6g}'.format(float(v))
    except (TypeError, ValueError):
        return ''


def build_gtex():
    n_in = n_out = n_skip = 0
    with gzip.open(GTEX_SRC, 'rt') as fin, gzip.open(GTEX_OUT, 'wt') as fout:
        head = fin.readline().rstrip('\n').split('\t')
        try:
            i_sig = head.index('bayests_sigma')
            i_pct = head.index('bayests_percentile')
        except ValueError:
            raise SystemExit('{} has no BayesTS columns: {}'.format(GTEX_SRC, head))
        fout.write(HEADER)
        for line in fin:
            n_in += 1
            t = line.rstrip('\n').split('\t')
            if len(t) <= max(i_sig, i_pct):
                n_skip += 1
                continue
            sig, pct = _fmt(t[i_sig]), _fmt(t[i_pct])
            if sig == '' or sig.lower() == 'nan':
                n_skip += 1
                continue
            fout.write('{}\t{}\t{}\n'.format(t[0].split('=')[0], sig, pct))
            n_out += 1
    return n_in, n_out, n_skip


def build_tabula_sapiens():
    n_in = n_out = n_skip = 0
    seen = set()
    with zipfile.ZipFile(TS_SRC) as z, gzip.open(TS_OUT, 'wt') as fout:
        with z.open(TS_MEMBER) as fin:
            head = fin.readline().decode('utf-8').rstrip('\n').split('\t')
            try:
                i_sig = head.index('mean_sigma')
                i_pct = head.index('percentile')
            except ValueError:
                raise SystemExit('{} has unexpected columns: {}'.format(TS_SRC, head))
            fout.write(HEADER)
            for raw in fin:
                n_in += 1
                t = raw.decode('utf-8', 'replace').rstrip('\n').split('\t')
                if len(t) <= max(i_sig, i_pct):
                    n_skip += 1
                    continue
                uid = t[0].split('=')[0]
                if uid in seen:      # the source has 2,860,898 rows / 2,852,713 unique UIDs
                    n_skip += 1
                    continue
                sig, pct = _fmt(t[i_sig]), _fmt(t[i_pct])
                if sig == '' or sig.lower() == 'nan':
                    n_skip += 1
                    continue
                seen.add(uid)
                fout.write('{}\t{}\t{}\n'.format(uid, sig, pct))
                n_out += 1
    return n_in, n_out, n_skip


if __name__ == '__main__':
    for label, fn, out in (('GTEx', build_gtex, GTEX_OUT),
                           ('TabulaSapiens', build_tabula_sapiens, TS_OUT)):
        n_in, n_out, n_skip = fn()
        size = os.path.getsize(out) / 1e6
        print('{:<15s} rows in={:>9d}  written={:>9d}  skipped={:>7d}  -> {} ({:.1f} MB)'
              .format(label, n_in, n_out, n_skip, out, size))
