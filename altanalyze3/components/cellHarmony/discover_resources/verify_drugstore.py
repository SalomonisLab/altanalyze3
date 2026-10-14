"""Deployment check for the complete native Discover signature resource.

Reads metadata and memory-maps arrays; does not rebuild, subset or score sources.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from . import drug_matching as dm


def verify(path=None):
    base = Path(path or dm._store_dir())
    manifest_bytes = (base / 'manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get('complete') is not True:
        raise ValueError('The source-signature resource manifest must declare complete=true.')
    if manifest.get('method_contract') != dm.METHOD_VERSION:
        raise ValueError('The DrugMatch resource method contract differs from the native engine.')
    resource = dm.Store(str(base))
    if resource.n_signatures != manifest['n_signatures']:
        raise ValueError('The declared complete signature roster does not match its partitions.')
    partitions = []
    for spec in manifest['partitions']:
        part = resource.partitions[spec['partition_id']]
        if len(part.genes) != spec['n_genes'] or len(set(part.keys)) != part.n_signatures:
            raise ValueError(f'{part.id}: gene-axis size or unique signature roster mismatch.')
        records = Path(part._records_path)
        if (len(part.offsets) != part.n_signatures + 1 or part.offsets[0] != 0
                or part.offsets[-1] != records.stat().st_size or np.any(np.diff(part.offsets) <= 0)):
            raise ValueError(f'{part.id}: signature metadata offsets do not cover the full records file.')
        if part.pvalue_available and (part.neglog10p is None or part.neglog10p.shape != part.values.shape):
            raise ValueError(f'{part.id}: declared source p-value array is missing or has the wrong shape.')
        if spec.get('landmark_mask') and (part.is_landmark is None or part.is_landmark.shape != (len(part.genes),)
                                         or int(np.count_nonzero(part.is_landmark)) != spec['n_landmark']):
            raise ValueError(f'{part.id}: declared measured-gene mask is missing or incomplete.')
        partitions.append({'id': part.id, 'genes': len(part.genes), 'signatures': part.n_signatures,
                           'signature_type': part.signature_type, 'shape': list(part.values.shape)})
    return {'status': 'ready', 'path': str(base), 'store_version': resource.version,
            'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
            'n_signatures': resource.n_signatures, 'partitions': partitions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--path', help='Complete existing drugstore_v1 directory; defaults to CELLREF2_DRUGSTORE.')
    args = parser.parse_args()
    print(json.dumps(verify(args.path), indent=2))


if __name__ == '__main__':
    main()
