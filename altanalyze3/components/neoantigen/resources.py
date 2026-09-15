"""Record and validate external reference/model bundles without downloading data."""
import argparse
import json
from pathlib import Path
from .io import checksum, manifest


def record_resources(roots, output):
    entries = []
    for name, root in roots.items():
        root = Path(root).resolve()
        if not root.is_dir():
            raise ValueError(f'Resource directory missing: {root}')
        for path in sorted(root.rglob('*')):
            if path.is_file():
                entries.append({'bundle': name, 'path': str(path.relative_to(root)),
                                'bytes': path.stat().st_size, 'sha256': checksum(path)})
    manifest(output, resources=entries)


def verify_resources(roots, source):
    data = json.loads(Path(source).read_text())
    for row in data['resources']:
        root = Path(roots[row['bundle']]).resolve()
        path = (root / row['path']).resolve()
        if not path.is_relative_to(root) or not path.is_file() or checksum(path) != row['sha256']:
            raise ValueError(f'Resource checksum mismatch: {row["bundle"]}/{row["path"]}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['record','verify'])
    p.add_argument('--reference', required=True)
    p.add_argument('--mhcflurry', required=True)
    p.add_argument('--manifest', required=True)
    a=p.parse_args()
    roots={'reference':a.reference, 'mhcflurry':a.mhcflurry}
    (record_resources if a.action=='record' else verify_resources)(roots,a.manifest)

if __name__=='__main__':main()
