"""Extract a small, reproducible gene subset from a junction count matrix.

No long-read features are extracted. Duplicate coordinate rows are merged by maximum,
not sum, to avoid counting alternative labels of the same measurement twice.
"""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np


def _open_text(path):
    """Open a counts matrix whether it is plain text or gzip.

    A network share can unmount mid-run, so every matrix is kept on local disk as .gz.
    Reading the compressed copy must give the same bytes as the uncompressed original.
    """
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def extract(path, genes):
    evidence = {}
    with _open_text(path) as fh:
        samples = fh.readline().rstrip('\n').split('\t')[1:]
        for line in fh:
            key, _, rest = line.partition('\t')
            if key.split(':', 1)[0] not in genes:
                continue
            m = re.search(r'=(chr[^:]+):(\d+)-(\d+)', key)
            if not m:
                continue
            uid = key.split('=')[0]
            a, b = sorted(map(int, m.group(2, 3)))
            counts = np.asarray(rest.rstrip('\n').split('\t'), dtype=float)
            if len(counts) != len(samples) or not np.isfinite(counts).all() or (counts < 0).any():
                raise ValueError('Invalid counts: ' + uid)
            gene = uid.split(':')[0]
            record = evidence.setdefault(gene, {}).setdefault((a, b),
                {'uid': uid, 'chrom': m[1], 'counts': counts})
            record['counts'] = np.maximum(record['counts'], counts)
    return samples, evidence


def extract_h5ad(path, genes):
    import anndata
    a = anndata.read_h5ad(path, backed='r')
    try:
        samples = a.obs_names.tolist()
        indices = [i for i, name in enumerate(a.var_names) if name.split(':')[0] in genes]
        sub = a[:, indices].X
        if hasattr(sub, 'toarray'):
            sub = sub.toarray()
        evidence = {}
        for col, idx in enumerate(indices):
            key = a.var_names[idx]
            m = re.search(r'=(chr[^:]+):(\d+)-(\d+)', key)
            if not m:
                continue
            uid = key.split('=')[0]
            gene = uid.split(':')[0]
            j = tuple(sorted(map(int, m.group(2, 3))))
            v = np.asarray(sub[:, col], dtype=float)
            if not np.isfinite(v).all() or (v < 0).any():
                raise ValueError('Invalid counts: ' + uid)
            record = evidence.setdefault(gene, {}).setdefault(j,
                {'uid': uid, 'chrom': m[1], 'counts': v})
            record['counts'] = np.maximum(record['counts'], v)
        return samples, evidence
    finally:
        a.file.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--genes-table', '--truth', dest='truth', required=True,
                    help='TSV with a gene column; only gene IDs are used, no long-read data is needed')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    genes = {r['gene'] for r in csv.DictReader(open(args.truth), delimiter='\t')}
    genes.update(['ENSG00000103257', 'ENSG00000197965', 'ENSG00000140090'])
    extractor = extract_h5ad if args.counts.endswith('.h5ad') else extract
    samples, evidence = extractor(args.counts, genes)
    import pickle
    with open(args.out, 'wb') as fh:
        pickle.dump((samples, evidence), fh)
    p = Path(args.counts)
    manifest = {'counts': str(p.resolve()), 'bytes': p.stat().st_size,
                'mtime_ns': p.stat().st_mtime_ns, 'genes': sorted(genes),
                'samples': len(samples), 'junctions': sum(map(len, evidence.values()))}
    Path(args.out + '.json').write_text(json.dumps(manifest, indent=2))
    print(manifest['samples'], 'samples;', manifest['junctions'], 'junctions', flush=True)


if __name__ == '__main__':
    main()
