"""Opt-in SNAF exports for the supplied OneClick iPepGen Galaxy workflow.

Only file contracts are implemented here; no Galaxy jobs or prediction services
are invoked. Junction coordinates are never treated as peptide coordinates.
"""
import csv
import json
import re
import tempfile
from pathlib import Path
from .hla import read_hla_table
from .io import checksum, read_table, write_table, manifest
from .proteomics import export_candidates, stable_id

PROFILE = Path(__file__).with_name('data') / 'ipepgen_contract.json'
ANNOTATION_FIELDS = ['Peptide', 'Chromosome', 'Start', 'End', 'Strand', 'Annotation',
                     'IGV_Genome_Coordinate', 'UCSC_Genome_Browser']


def extract_contract(path):
    """Extract executable input requirements from an embedded Galaxy .ga file."""
    workflow = json.loads(Path(path).read_text())
    if workflow.get('a_galaxy_workflow') != 'true' and workflow.get('a_galaxy_workflow') is not True:
        raise ValueError('Not a Galaxy workflow export')
    found = {}
    def visit(w, prefix=''):
        for key, step in w.get('steps', {}).items():
            where = prefix + str(key)
            if step.get('subworkflow'):
                visit(step['subworkflow'], where + '/')
            tool = step.get('tool_id') or ''
            for name, token in [('fragpipe', '/fragpipe/fragpipe/'), ('pepquery', '/pepquery2/pepquery2/'), ('iedb', '/iedb_api/iedb_api/')]:
                if token in tool:
                    if name in found:
                        raise ValueError(f'Ambiguous workflow: multiple {name} tools')
                    state = step.get('tool_state') or {}
                    found[name] = (where, tool, json.loads(state) if isinstance(state, str) else state)
    visit(workflow)
    if set(found) != {'fragpipe', 'pepquery', 'iedb'}:
        raise ValueError('Workflow must contain FragPipe, PepQuery2, and IEDB tools')
    fp, pq, hla = (found[n][2] for n in ('fragpipe', 'pepquery', 'iedb'))
    try:
        digestion = fp['wf']['msfragger']['digestion']
        if pq['req_inputs']['input_type']['input_type_selector'] != 'peptide':
            raise ValueError('PepQuery must use peptide input')
        if hla['prediction']['tool'] != 'mhci' or hla['sequence']['seqsrc'] != 'fasta' or hla['prediction']['alleles']['allelesrc'] != 'history':
            raise ValueError('Expected IEDB class-I FASTA and history allele inputs')
        lengths = sorted({int(n) for n in hla['prediction']['lengths']})
        bounds = {'fragpipe': [int(digestion['digest_min_length']), int(digestion['digest_max_length'])],
                  'pepquery': [int(pq['ms_params']['search']['min_length']), int(pq['ms_params']['search']['max_length'])]}
        if not lengths or min(lengths) < 1 or any(a < 1 or b < a for a,b in bounds.values()):
            raise ValueError('Invalid workflow peptide length settings')
    except (KeyError, TypeError) as exc:
        raise ValueError('Unsupported iPepGen workflow tool-state schema') from exc
    return dict(schema_version='1.0', name=workflow['name'], source_sha256=checksum(path),
                source_filename=Path(path).name, bounds=bounds, iedb_lengths=lengths,
                tools={n: dict(step=v[0], tool_id=v[1]) for n,v in found.items()},
                fasta_header='generic|SNAF_<stable source hash>|description',
                pepquery_format='one unmodified peptide per line; no header',
                iedb_alleles_format='one HLA-A/B/C*NN:NN allele per line; no header',
                annotation_columns=ANNOTATION_FIELDS,
                annotation_requirement='SNAF accession join required; existing StringTie/variant positional parser does not support SNAF')


def candidate_reports(outdir):
    """The combined report supersedes its per-sample inputs (it adds in_db)."""
    directory = Path(outdir) / 'T_candidates'
    combined = directory / 'T_antigen_candidates_all.txt'
    paths = [combined] if combined.exists() else sorted(directory.glob('T_antigen_candidates_*.txt'))
    if not paths:
        raise ValueError(f'No SNAF candidate reports in {directory}')
    return paths


def _bed_rows(path, sources):
    """Validate supplied peptide BED12 mappings; names are SNAF accessions.

    Coding blocks must account for the complete peptide. This is structural
    validation of supplied mappings, not independent biological confirmation.
    """
    result = []
    if not path:
        return result
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            row = line.rstrip('\r\n').split('\t')
            if len(row) != 12 or row[3] not in sources:
                raise ValueError(f'BED line {number}: require 12 columns and a known SNAF accession')
            chrom, start, end, accession, score, strand, thick_start, thick_end, rgb, count, sizes, starts = row
            try:
                start,end,score,thick_start,thick_end,count=map(int,[start,end,score,thick_start,thick_end,count])
                sizes=[int(x) for x in sizes.rstrip(',').split(',')]
                starts=[int(x) for x in starts.rstrip(',').split(',')]
            except ValueError as exc:
                raise ValueError(f'BED line {number}: invalid numeric field') from exc
            if (not re.fullmatch(r'[A-Za-z0-9_.-]+',chrom) or strand not in ('+','-') or start<0 or end<=start
                or not 0<=score<=1000 or thick_start!=start or thick_end!=end or count<1
                or len(sizes)!=count or len(starts)!=count or starts[0]!=0
                or any(x<=0 for x in sizes) or any(x<0 for x in starts)
                or any(starts[i]+sizes[i]>starts[i+1] for i in range(count-1))
                or starts[-1]+sizes[-1]!=end-start or sum(sizes)!=3*len(sources[accession]['peptide'])):
                raise ValueError(f'BED line {number}: invalid complete-peptide coding blocks')
            if not (rgb == '0' or re.fullmatch(r'\d{1,3},\d{1,3},\d{1,3}',rgb) and all(int(x)<=255 for x in rgb.split(','))):
                raise ValueError(f'BED line {number}: invalid RGB')
            result.append(row)
    return sorted({tuple(row) for row in result}, key=lambda r:(r[0],int(r[1]),int(r[2]),r[3]))


def _fasta(path, records):
    with Path(path).open('w') as handle:
        for title, sequence in records:
            handle.write(f'>{title}\n{sequence}\n')


def export_galaxy(paths, hla, outdir, workflow=None, predictor='', sample=None, peptide_bed=None, assembly='hg38', sample_ids=None):
    """Write cohort search files and separate sample files for raw-MS/HLA routes."""
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', assembly):
        raise ValueError('Invalid genome assembly identifier')
    contract = extract_contract(workflow) if workflow else json.loads(PROFILE.read_text())
    alleles = read_hla_table(hla)
    if sample_ids is not None:
        missing = set(sample_ids) - set(alleles)
        if missing:
            raise ValueError(f'Missing HLA samples: {sorted(missing)}')
        alleles = {name: alleles[name] for name in sample_ids}
    if sample is not None:
        if sample not in alleles:
            raise ValueError(f'Sample {sample!r} is absent from HLA table')
        alleles = {sample: alleles[sample]}
    if not alleles or any(not v for v in alleles.values()):
        raise ValueError('Galaxy integration requires nonempty HLA calls for each exported sample')
    paths = list(paths)
    # Validate before writing outputs, and reuse the independent exchange identity contract.
    with tempfile.TemporaryDirectory(prefix='snaf-galaxy-export-') as temp:
        bundle = export_candidates(paths, temp, predictor=predictor)
        candidates, fields = read_table(bundle/'candidates.tsv')
    if sample is not None:
        candidates = [r for r in candidates if r['sample_id']==sample]
    missing = {r['sample_id'] for r in candidates} - set(alleles)
    if missing:
        raise ValueError(f'Missing HLA samples: {sorted(missing)}')
    for row in candidates:
        if row['hla_allele'] and row['hla_allele'] not in alleles[row['sample_id']]:
            raise ValueError(f"Candidate allele {row['hla_allele']} is absent from supplied calls for {row['sample_id']}")
    sources = {}
    for row in candidates:
        row['accession'] = 'SNAF_' + row['source_id'].removeprefix('snaf_')
        row['fasta_id'] = f"generic|{row['accession']}|SNAF_splice_candidate"
        sources[row['accession']] = row
    bed = _bed_rows(peptide_bed, sources)
    outdir = Path(outdir)
    # A fresh destination prevents stale sample files from being discovered by Galaxy.
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f'Galaxy export directory must be empty: {outdir}')
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'by_sample').mkdir()
    fp_min,fp_max = contract['bounds']['fragpipe']; pq_min,pq_max = contract['bounds']['pepquery']
    def export_files(rows, base):
        selected = {r['accession']:r for r in rows}
        _fasta(str(base)+'.database.fasta',[(f'generic|{a}|SNAF_splice_candidate',selected[a]['peptide'])
               for a in sorted(selected) if fp_min<=len(selected[a]['peptide'])<=fp_max])
        peptides = sorted({r['peptide'] for r in rows})
        Path(str(base)+'.pepquery.txt').write_text(''.join(p+'\n' for p in peptides if pq_min<=len(p)<=pq_max))
        _fasta(str(base)+'.iedb.fasta',[(stable_id('peptide_',p),p) for p in peptides if len(p) in contract['iedb_lengths']])
    export_files(candidates,outdir/'snaf')
    sample_rows=[]
    for name, calls in sorted(alleles.items()):
        token = re.sub(r'[^A-Za-z0-9_-]','_',name)[:60] + '_' + stable_id('',name)[:10]
        selected=[r for r in candidates if r['sample_id']==name]
        base=outdir/'by_sample'/token
        export_files(selected,base)
        Path(str(base)+'.alleles.txt').write_text(''.join(a+'\n' for a in sorted(calls)))
        sample_rows.append(dict(sample_id=name,collection_id=token,candidate_rows=len(selected),
            database=f'by_sample/{token}.database.fasta',pepquery=f'by_sample/{token}.pepquery.txt',
            iedb_fasta=f'by_sample/{token}.iedb.fasta',iedb_alleles=f'by_sample/{token}.alleles.txt'))
    write_table(outdir/'samples.tsv',sample_rows,['sample_id','collection_id','candidate_rows','database','pepquery','iedb_fasta','iedb_alleles'])
    write_table(outdir/'candidate_map.tsv',candidates,['accession','fasta_id']+[f for f in fields if f not in ('accession','fasta_id')])
    length_rows=[]
    for row in candidates:
        n=len(row['peptide'])
        length_rows.append(dict(candidate_id=row['candidate_id'],sample_id=row['sample_id'],peptide=row['peptide'],length=n,
            fragpipe_eligible=fp_min<=n<=fp_max,pepquery_eligible=pq_min<=n<=pq_max,iedb_eligible=n in contract['iedb_lengths']))
    write_table(outdir/'length_eligibility.tsv',length_rows,['candidate_id','sample_id','peptide','length','fragpipe_eligible','pepquery_eligible','iedb_eligible'])
    with (outdir/'peptides.bed').open('w') as handle:
        for row in bed:
            handle.write('\t'.join(row)+'\n')
    annotations=[]; mapped={r[3] for r in bed}
    for row in bed:
        source=sources[row[3]]; chrom,start,end=row[0],int(row[1]),int(row[2])
        annotations.append(dict(zip(ANNOTATION_FIELDS,[source['peptide'],chrom,start,end,row[5],
            f"SNAF;accession={row[3]};event={source['event_id']};supplied_peptide_BED12",
            f'{chrom}:{start+1}-{end}',f'https://genome.ucsc.edu/cgi-bin/hgTracks?db={assembly}&position={chrom}%3A{start+1}-{end}'])))
    for accession in sorted(set(sources)-mapped):
        row=sources[accession]
        annotations.append(dict(Peptide=row['peptide'],Annotation=f"SNAF;accession={accession};event={row['event_id']};peptide_coordinates_unavailable"))
    write_table(outdir/'annotations.tsv',annotations,ANNOTATION_FIELDS)
    write_table(outdir/'hla_alleles.tsv',[dict(sample_id=s,allele=a) for s in sorted(alleles) for a in sorted(alleles[s])],['sample_id','allele'])
    manifest(outdir/'integration.json',kind='snaf_galaxy_ipepgen',contract=contract,assembly=assembly,
        inputs=[dict(path=str(p),sha256=checksum(p)) for p in paths],hla_sha256=checksum(hla),
        peptide_bed_sha256=checksum(peptide_bed) if peptide_bed else None,
        candidates=len(candidates),sources=len(sources),samples=len(sample_rows),mapped_sources=len(mapped),
        semantics='Reported candidate peptides; supplied BED12 mappings only; no inferred peptide coordinates',
        routes={'fragpipe':'Add snaf.database.fasta as an extra input to FASTA merge at step 11/4; retain reference/nonreference/fusion inputs and FragPipe decoy generation',
                'pepquery':'Connect matching by_sample/*.pepquery.txt directly to step 12 input 0; no header removal',
                'iedb':'Pair matching by_sample/*.iedb.fasta and *.alleles.txt at step 14 inputs 1 and 0',
                'annotation':'Join FragPipe accessions to candidate_map.tsv; SNAF requires a separate branch from step 13 StringTie/variant parsing'},
        outputs={str(p.relative_to(outdir)):checksum(p) for p in sorted(outdir.rglob('*')) if p.is_file()})
    return outdir
