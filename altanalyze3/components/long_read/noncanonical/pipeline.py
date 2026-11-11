"""Run the validated detectors with a complete, source-qualified manifest.

Each subprocess uses this installed package, not a project-specific fork. Expensive
steps have independent completion receipts; a failed branch cannot masquerade as an
empty successful result. Ordinary analysis entry points remain the native functions.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from .manifest import (MODALITIES, atomic_json, digest, load_config, receipt,
                       verify_receipt, verify_source, write_tsv)

HERE = Path(__file__).resolve().parent

def run_command(command, log, config, cwd=None):
    env = dict(os.environ, ALTANALYZE3_NONCANONICAL_CONFIG=config['_path'],
               PYTHONHASHSEED='0', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    # Include this checkout for child module invocations; preserve installed dependencies.
    package_root = str(HERE.parents[3])
    env['PYTHONPATH'] = package_root + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with open(log, 'w') as handle:
        handle.write(json.dumps([str(x) for x in command]) + '\n')
        handle.flush()
        subprocess.run([str(x) for x in command], check=True, stdout=handle,
                       stderr=subprocess.STDOUT, env=env, cwd=cwd)

def module_command(module, *args):
    return [sys.executable, '-m', module, *map(str, args)]

def _complete_or_run(done, config, inputs, call, paths):
    if done.exists():
        return verify_receipt(done, config)
    call()
    return receipt(done, config, inputs, [p for p in paths() if p.is_file()], status='COMPLETE')

def extract_source(config, source, standard=True):
    root = Path(config['output'])
    sid = source['source_id']
    verify_source(source)
    inputs = {'bam': source['bam_sha256']}
    branches = [x for x in config['modalities'] if standard or x != 'standard']
    for branch in branches:
        folder = {'standard':'compatibility', 'circrna':'circrna',
                  'oncofusion':'oncofusion', 'trans-splicing':'trans_splicing'}[branch]
        out = root / folder / 'samples' / sid
        out.mkdir(parents=True, exist_ok=True)
        done = out / 'COMPLETE.json'
        if done.exists():
            verify_receipt(done, config)
            continue
        prefix = out / sid
        if branch == 'standard':
            cmd = module_command('altanalyze3.components.bam.isoform_structure_extract', source['bam'],
                                 '--gene-model', config['exon_model'], '--output-prefix', prefix, '--max-processes', '1')
        elif branch == 'circrna':
            cmd = module_command('altanalyze3.components.bam.noncolinear_rna', '--bam', source['bam'],
                                 '--gene-model', config['exon_model'], '--genome-fasta', config['genome'],
                                 '--output-prefix', prefix, '--library', sid)
        elif branch == 'trans-splicing':
            ref = Path(config['reference'])
            cmd = module_command('altanalyze3.components.bam.isoform_structure_extract', source['bam'],
                                 '--gene-model', config['exon_model'], '--output-prefix', prefix,
                                 '--cross-chromosome-only', '--library', sid,
                                 '--molecule-id-prefix', sid + '__',
                                 '--blacklist-members', ref/'blacklist_members.tsv',
                                 '--genome-index', ref/'genome.splice.mmi', '--genome-fasta', config['genome'])
        else:
            cmd = module_command('altanalyze3.components.long_read.noncanonical.oncofusion.scan_bam',
                                 sid, source['bam'], out, Path(config['reference'])/'genes.tsv.gz',
                                 config.get('fusion_panel', HERE/'oncofusion/aml_mds_fusion_genes.txt'))
        run_command(cmd, out/'extract.log', config)
        if branch == 'trans-splicing':
            stats = json.loads(Path(str(prefix)+'.cross_chromosome_stats.json').read_text())
            atomic_json(out/'STAGE1_COMPLETE.json', {'library':sid, 'extract_stats':stats})
        # A successful empty source still has its summary and raw table with a header.
        files = sorted(p for p in out.iterdir() if p.is_file() and p.name != 'COMPLETE.json')
        receipt(done, config, inputs, files, status='COMPLETE', source_id=sid, modality=branch)
    verify_source(source)

def prepare_metadata(config):
    root = Path(config['output'])
    rows = []
    for i, s in enumerate(config['sources'], 1):
        rows.append(dict(index=i, library=s['source_id'], cohort=s.get('cohort',s['group']),
                         group=s['group'], original_uid=s['uid'], original_library=s['library'],
                         reverse=s['reverse'], bam=s['bam']))
    fields=['index','library','cohort','group','original_uid','original_library','reverse','bam']
    write_tsv(root/'trans_splicing/config/samples.tsv', fields, rows)
    # The scan's uid is source-qualified. Its legacy per-source report never merges
    # source reads by repeated library labels; biological pooling is done in ledger.py.
    write_tsv(root/'oncofusion/metadata.tsv', ['uid','library','groups'],
              ({'uid':s['source_id'],'library':'','groups':s['group']} for s in config['sources']))
    atomic_json(root/'provenance/SOURCES.json', config['sources'])

def require_sources(config, branches=None):
    branches = branches or config['modalities']
    folder = {'standard':'compatibility','circrna':'circrna','oncofusion':'oncofusion','trans-splicing':'trans_splicing'}
    for source in config['sources']:verify_source(source)
    for branch in branches:
        for source in config['sources']:
            done=Path(config['output'])/folder[branch]/'samples'/source['source_id']/'COMPLETE.json'
            if not done.exists():
                raise ValueError(f'Incomplete requested modality: {branch} / {source["source_id"]}')
            verify_receipt(done,config)

def aggregate(config):
    root=Path(config['output']);require_sources(config)
    prepare_metadata(config)
    if 'circrna' in config['modalities']:
        prepare_circrna_gate(config)
        out=root/'circrna/evidence';done=root/'circrna/AGGREGATE_COMPLETE.json'
        _complete_or_run(done,config,{},lambda:run_command(module_command(
            'altanalyze3.components.bam.noncolinear_evidence','--root',root/'circrna','--out',out),
            root/'circrna/aggregate.log',config),lambda:list(out.iterdir()))
    if 'oncofusion' in config['modalities']:
        from .onco_profile import prepare
        out=prepare(config);done=out/'AGGREGATE_COMPLETE.json'
        _complete_or_run(done,config,{},lambda:run_command(module_command(
            'altanalyze3.components.long_read.noncanonical.oncofusion.combine',out,config['reference'],
            config['genome'],out/'metadata.tsv',config['reference_gff']),out/'aggregate.log',config),
            lambda:list(out.glob('*.tsv'))+[out/'SUMMARY.json'])
    if 'trans-splicing' in config['modalities']:
        out=root/'trans_splicing'
        run_command(module_command('altanalyze3.components.long_read.noncanonical.xchr_catalog','--config',config['_path']),out/'catalog.log',config)
        done=out/'REPORT_COMPLETE.json'
        _complete_or_run(done,config,{},lambda:run_command(module_command(
            'altanalyze3.components.long_read.noncanonical.xchr.report_stream',out),out/'report.log',config),
            lambda:list((out/'results').iterdir()))
        done=out/'CLASSIFICATION_COMPLETE.json'
        _complete_or_run(done,config,{},lambda:run_command(module_command(
            'altanalyze3.components.long_read.noncanonical.xchr.classify_joints',out,
            config.get('workers',1),'results','classification',config['patients']),out/'classify.log',config),
            lambda:list((out/'classification').iterdir()))
    from .ledger import build_ledger
    build_ledger(config)


def prepare_circrna_gate(config):
    """Bind the existing aggregator's complete-roster gate to the independent manifest."""
    import gzip
    root=Path(config['output'])/'circrna'
    expected=[s['source_id'] for s in config['sources']]
    run=root/'RUN_MANIFEST.json'
    atomic_json(run,dict(expected_ordered_libraries=expected,configuration_sha256=config['_sha256'],
                         source_roster=config['sources'],scope='declared manifest sources'))
    libraries=[]
    for sid in expected:
        folder=root/'samples'/sid
        proof=verify_receipt(folder/'COMPLETE.json',config)
        raw=folder/(sid+'.noncolinear_joints.tsv.gz')
        with gzip.open(raw,'rt') as handle:count=sum(1 for _ in handle)-1
        if count<0:raise ValueError('Missing raw table header: '+sid)
        libraries.append(dict(library=sid,pair_rows=count,files={'raw':{'sha256':digest(raw)}}))
    atomic_json(root/'VALIDATED_EVIDENCE.json',dict(expected_ordered_libraries=expected,
        manifest_sha256=digest(run),libraries=libraries,validation_scope='source completion, hashes and pair retention; not physical topology'))


def finalize(config):
    root=Path(config['output'])
    require_sources(config)
    proofs=[root/'reconciliation/COMPLETE.json',root/'integrated/events/COMPLETE.json']
    if 'circrna' in config['modalities']:proofs.append(root/'circrna/AGGREGATE_COMPLETE.json')
    if 'oncofusion' in config['modalities']:proofs.append(root/'oncofusion/native_profile/AGGREGATE_COMPLETE.json')
    if 'trans-splicing' in config['modalities']:proofs.extend([root/'trans_splicing/CATALOG_COMPLETE.json',root/'trans_splicing/REPORT_COMPLETE.json',root/'trans_splicing/CLASSIFICATION_COMPLETE.json'])
    if 'standard' in config['modalities']:
        proofs += [root/'compatibility/STANDARD_COMPLETE.json',root/'integrated/standard/STANDARD_COMPLETE.json']
    for path in proofs:verify_receipt(path,config)
    receipt(root/'COMPLETE.json',config,{},proofs,status='COMPLETE_REQUESTED_DETECTION_AND_COUNTING',
        modalities=config['modalities'],
        validation_scope='discovery, native classification, read-accounting and standard views; sequence-assessment outputs are separately scoped',
        statistical_scope='native standard analyses from standard_parameters; noncanonical descriptive counts preserve every event',
        physical_topology='candidate labels do not establish RNA circularity or genomic rearrangement')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True)
    parser.add_argument('--stage',choices=['extract','aggregate','standard','all'],default='all')
    parser.add_argument('--source',help='One exact manifest source ID (extract stage only)')
    args=parser.parse_args(argv);config=load_config(args.config)
    if args.source and args.stage!='extract':
        parser.error('--source is valid only with --stage extract')
    sources=config['sources']
    if args.source:
        sources=[s for s in sources if s['source_id']==args.source]
        if len(sources)!=1:raise ValueError('Unknown source ID: '+args.source)
    if args.stage in ('all','extract'):
        for s in sources:extract_source(config,s)
    if args.stage in ('all','aggregate'):aggregate(config)
    if args.stage in ('all','standard') and 'standard' in config['modalities']:
        from .standard import run_standard
        run_standard(config,'compatibility')
        run_standard(config,'integrated')
    if args.stage in ('all','standard'):
        from .event_counts import run as event_counts
        event_counts(config)
        finalize(config)

if __name__=='__main__':main()
