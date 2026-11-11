"""Pinned inputs, source identities and atomic completion receipts."""
from __future__ import annotations
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

MODALITIES = ('standard', 'circrna', 'oncofusion', 'trans-splicing')

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def code_digest():
    """Fingerprint installed Python and resource bytes, independent of checkout path."""
    root = Path(__file__).resolve().parents[3]
    h = hashlib.sha256()
    for path in sorted(root.rglob('*')):
        if path.is_file() and path.suffix in ('.py', '.json', '.txt', '.gff3') and '__pycache__' not in path.parts:
            h.update(str(path.relative_to(root)).encode())
            h.update(bytes.fromhex(digest(path)))
    return h.hexdigest()


def verify_reference_pins(config):
    pins = config.get('reference_pins')
    if not pins:
        raise ValueError('Independently frozen reference_pins are required before execution')
    for path, expected in pins.items():
        if digest(path) != expected:
            raise ValueError('Reference bytes changed: ' + path)
    return pins


def verify_barcode_roster(path, source_ids):
    """Reject ambiguous annotations before the native reader can deduplicate them."""
    seen=set();covered=set()
    with open(path) as handle:
        for number,row in enumerate(csv.reader(handle,delimiter='\t'),1):
            if number==1 and row and row[0].lower() in ('cellbarcode','barcode','barcode_cluster'):
                continue
            if len(row)<2 or not row[0] or not row[1] or row[0].count('.')!=1:
                raise ValueError(f'Unresolved barcode annotation at {path}:{number}')
            if row[0] in seen:raise ValueError('Duplicate barcode annotation identity: '+row[0])
            seen.add(row[0]);covered.add(row[0].split('.')[1])
    if not set(source_ids)<=covered:
        raise ValueError('Sources missing from declared barcode roster: '+str(sorted(set(source_ids)-covered)))
    return len(seen)


def verify_executables(config):
    required={'samtools'} if 'oncofusion' in config['modalities'] else set()
    pins=config.get('runtime_executables',{})
    if not required<=pins.keys():raise ValueError('Pinned native executable required: '+str(sorted(required-pins.keys())))
    for command,pin in pins.items():
        executable=shutil.which(command)
        if not executable or digest(executable)!=pin['sha256']:
            raise ValueError('Pinned native executable missing or changed: '+command)
        version=subprocess.check_output([executable,'--version'],text=True,errors='replace').splitlines()[0]
        if version!=pin['version']:raise ValueError('Native executable version changed: '+command)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with tmp.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')
    os.replace(tmp, path)

def read_tsv(path):
    import gzip
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path, fields, rows):
    import gzip
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)

def load_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    config['_path'] = str(path)
    config['_sha256'] = digest(path)
    config['output'] = str(Path(config['output']).resolve())
    modalities = config.get('modalities', list(MODALITIES))
    if not modalities or len(set(modalities)) != len(modalities) or set(modalities) - set(MODALITIES):
        raise ValueError('Unknown, duplicate or empty modality selection')
    config['modalities'] = modalities
    sources = config['sources']
    if not sources or len({r['source_id'] for r in sources}) != len(sources):
        raise ValueError('An independent, nonempty unique source roster is required')
    for source in sources:
        sid = source['source_id']
        if not sid or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in sid):
            raise ValueError('source_id must be a safe stable identifier: ' + sid)
        for key in ('bam', 'uid', 'library', 'group', 'reverse'):
            if key not in source or source[key] == '':
                raise ValueError(f'{sid}: missing {key}')
        if not Path(source['bam']).is_file():
            raise FileNotFoundError(source['bam'])
        if not source.get('index_sha256') or not source.get('index') or not Path(source['index']).is_file():
            raise ValueError(f'{sid}: pinned BAM index required')
        if not source.get('bam_sha256'):
            raise ValueError(f'{sid}: independent BAM hash required')
    required = {'exon_model', 'barcode_clusters'}
    if 'standard' in modalities:
        required |= {'genome','reference_gff','gene_symbols'}
    if set(modalities) - {'standard'}:
        required.add('genome')
    if 'oncofusion' in modalities:
        required |= {'reference', 'reference_gff'}
    if 'trans-splicing' in modalities:
        required |= {'reference', 'reference_gff', 'gene_symbols', 'hgnc', 'barcode_clusters', 'patients', 'expression_dir'}
    for key in required:
        if key not in config or not Path(config[key]).exists():
            raise FileNotFoundError(f'Required {key} unavailable: {config.get(key)}')
    config['_code_sha256'] = code_digest()
    if config.get('code_sha256') != config['_code_sha256']:
        raise ValueError('Package code differs from the pinned implementation')
    expected_pins={str(Path(config[k])) for k in required if Path(config[k]).is_file()}
    if 'reference' in required:
        expected_pins.update(str(Path(config['reference'])/name) for name in ('genes.tsv.gz','fusion_db.tsv','blacklist_members.tsv','genome.splice.mmi','homology_null.json','transcripts.tsv.gz'))
    if 'expression_dir' in required:
        expected_pins.update(str(p) for p in Path(config['expression_dir']).glob('*.h5ad'))
    if not expected_pins<=set(config.get('reference_pins',{})):
        raise ValueError('Required reference resource hashes are missing')
    verify_reference_pins(config)
    verify_barcode_roster(config['barcode_clusters'],[s['source_id'] for s in sources])
    verify_executables(config)
    if config.get('runtime_versions'):
        import importlib.metadata
        for package,version in config['runtime_versions'].items():
            if importlib.metadata.version(package)!=version:
                raise ValueError('Pinned runtime version differs: '+package)
    # Sources are inspected afresh at each execution boundary; receipts bind the original config.
    return config

def verify_source(source):
    if digest(source['index']) != source['index_sha256']:
        raise ValueError('BAM index content changed: ' + source['source_id'])
    if digest(source['bam']) != source['bam_sha256']:
        raise ValueError('BAM content changed: ' + source['source_id'])

def receipt(path, config, inputs, outputs, **extra):
    atomic_json(path, dict(config_sha256=config['_sha256'], code_sha256=(config['_code_sha256'] if '_code_sha256' in config else code_digest()), inputs=inputs,
                           outputs={str(Path(p).resolve()): digest(p) for p in outputs}, **extra))

def verify_receipt(path, config):
    record = json.loads(Path(path).read_text())
    if record['config_sha256'] != config['_sha256']:
        raise ValueError('Existing output belongs to another configuration: ' + str(path))
    if record.get('code_sha256') != (config['_code_sha256'] if '_code_sha256' in config else code_digest()):
        raise ValueError('Completed output belongs to changed package code: ' + str(path))
    for name, hashed in record['outputs'].items():
        if digest(name) != hashed:
            raise ValueError('Completed output changed: ' + name)
    return record
