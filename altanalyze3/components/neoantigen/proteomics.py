"""SNAF TSV/FASTA exchange with an independently installed pyNeoQuant library."""
import hashlib
import importlib.metadata
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

from .io import read_table, write_table, checksum, manifest

CANDIDATE_FIELDS = ['candidate_id', 'source_id', 'sample_id', 'peptide', 'event_id',
                    'gene', 'isoform_id', 'hla_allele', 'hla_binding_rank',
                    'hla_predictor', 'hla_score_type', 'immunogenicity']
SOURCE_FIELDS = ['source_id', 'source_type', 'fasta_path', 'fasta_record', 'gene',
                 'transcript_id', 'protein_id', 'isoform_id', 'event_id']


def stable_id(prefix, *values):
    return prefix + hashlib.sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()[:24]


def export_candidates(paths, outdir, canonical_fasta=None, predictor=''):
    """Export reported SNAF peptides; never infer a full protein or junction geometry.

    Each source is the event/isoform/peptide association; each candidate additionally
    identifies its sample and allele. Candidate metadata retains all original fields.
    """
    from .hla import normalize_class_i
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    sources, candidates, originals = {}, {}, []
    for path in paths:
        rows, fields = read_table(path)
        originals.extend(f for f in fields if f not in originals)
        for row in rows:
            sample = row.get('sample_id') or row.get('sample')
            peptide = (row.get('peptide') or '').upper()
            event = row.get('event_id') or row.get('uid')
            isoform = row.get('isoform_id', '')
            if not sample or not event or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', peptide):
                raise ValueError(f'{path}: every candidate needs sample, event_id/uid and an unmodified peptide')
            allele = normalize_class_i(row.get('hla_allele') or row.get('hla') or '')
            source_id = stable_id('snaf_', event, isoform, peptide)
            candidate_id = stable_id('candidate_', source_id, sample, allele)
            gene = row.get('gene') or event.split(':')[0]
            source = dict(source_id=source_id, source_type='junction', fasta_path='candidates.fasta',
                          fasta_record=source_id, gene=gene, transcript_id=row.get('transcript_id', ''),
                          protein_id=row.get('protein_id', ''), isoform_id=isoform, event_id=event)
            if source_id in sources and sources[source_id][0] != source:
                raise ValueError(f'Conflicting source metadata: {source_id}')
            sources[source_id] = (source, peptide)
            record = {**row, 'candidate_id': candidate_id, 'source_id': source_id,
                      'sample_id': sample, 'peptide': peptide, 'event_id': event,
                      'gene': gene, 'isoform_id': isoform, 'hla_allele': allele,
                      'hla_binding_rank': row.get('hla_binding_rank') or row.get('binding_affinity', ''),
                      'hla_predictor': row.get('hla_predictor') or predictor,
                      'hla_score_type': row.get('hla_score_type') or ('presentation_percentile' if predictor == 'MHCflurry' else 'el_percentile' if predictor == 'netMHCpan' else 'unspecified'),
                      'immunogenicity': row.get('immunogenicity', '')}
            if candidate_id in candidates and candidates[candidate_id] != record:
                raise ValueError(f'Conflicting duplicate candidate: {candidate_id}')
            candidates[candidate_id] = record
    with (outdir / 'candidates.fasta').open('w') as f:
        for sid in sorted(sources):
            f.write(f'>{sid}\n{sources[sid][1]}\n')
    source_rows = [sources[sid][0] for sid in sorted(sources)]
    if canonical_fasta:
        shutil.copyfile(canonical_fasta, outdir / 'canonical.fasta')
        source_rows.append(dict(source_id='{fasta_title}', source_type='canonical', fasta_path='canonical.fasta'))
    write_table(outdir / 'source_manifest.tsv', source_rows, SOURCE_FIELDS)
    fields = CANDIDATE_FIELDS + [f for f in originals if f not in CANDIDATE_FIELDS]
    write_table(outdir / 'candidates.tsv', [candidates[k] for k in sorted(candidates)], fields)
    manifest(outdir / 'export_manifest.json', kind='snaf_peptide_export',
             source_semantics='reported_candidate_peptide; geometry not inferred',
             inputs=[{'path': str(p), 'sha256': checksum(p)} for p in paths],
             canonical_sha256=checksum(canonical_fasta) if canonical_fasta else None,
             sources=len(sources), candidates=len(candidates))
    return outdir


def merge_evidence(candidates_path, evidence_path, output, q_threshold=0.01):
    """Keep all candidates and summarize only matching sample/source evidence.

    MS supports a peptide/source association, not an HLA assignment. Evidence with
    no allele applies to each corresponding candidate; explicit alleles must match.
    Decoys/entrapments and missing or nonfinite q-values never count as passing.
    """
    from .hla import normalize_class_i
    if not 0 <= q_threshold <= 1:
        raise ValueError('q_threshold must be in [0, 1]')
    candidates, fields = read_table(candidates_path)
    evidence, _ = read_table(evidence_path)
    idx = {}
    for row in evidence:
        key = (row.get('sample_id'), row.get('peptide'), row.get('source_id'))
        if not all(key):
            continue
        idx.setdefault(key, []).append(row)
    output_rows = []
    for candidate in candidates:
        matches = idx.get((candidate['sample_id'], candidate['peptide'], candidate['source_id']), [])
        allele = normalize_class_i(candidate.get('hla_allele', ''))
        matches = [r for r in matches if not r.get('hla_allele') or normalize_class_i(r['hla_allele']) == allele]
        passing, qvalues, spectra, ambiguity = [], [], set(), set()
        for row in matches:
            if row.get('source_type', '').lower() in ('decoy', 'entrapment'):
                continue
            try:
                q = float(row.get('q_value', ''))
            except (ValueError, TypeError):
                continue
            if not math.isfinite(q) or not 0 <= q <= 1:
                continue
            qvalues.append(q)
            if q <= q_threshold:
                passing.append(row)
                spectrum = row.get('spectrum_id') or row.get('scan')
                if spectrum:
                    spectra.add(spectrum)
                ambiguity.update(v for v in row.get('competing_sources', '').split(';') if v)
        output_rows.append({**candidate, 'ms_evidence_rows': len(matches),
                            'ms_passing_rows': len(passing), 'ms_distinct_spectra': len(spectra),
                            'ms_min_q_value': min(qvalues) if qvalues else '',
                            'ms_status': 'peptide_source_supported' if passing else 'no_passing_evidence',
                            'ms_competing_sources': ';'.join(sorted(ambiguity))})
    write_table(output, output_rows, fields + ['ms_evidence_rows', 'ms_passing_rows', 'ms_distinct_spectra',
                'ms_min_q_value', 'ms_status', 'ms_competing_sources'])
    return output


def run_pyneoquant(bundle, psm_table, sample_id, outdir, executable='pyneoquant', search_format='generic', q_threshold=0.01):
    """Use the external CLI and versioned file contract; no pyNeoQuant import."""
    executable_path = shutil.which(executable) or (str(Path(executable).resolve()) if Path(executable).is_file() else None)
    if not executable_path:
        raise RuntimeError('Optional pyNeoQuant is not installed. Install its separate library and put pyneoquant on PATH, or pass --executable.')
    capabilities = subprocess.run([executable_path, 'schema-version'], capture_output=True, text=True, check=True)
    if 'fasta_record' not in json.loads(capabilities.stdout).get('features', []):
        raise RuntimeError('pyNeoQuant lacks the fasta_record contract; install version 0.1.0a1 or newer')
    bundle, outdir = Path(bundle).resolve(), Path(outdir).resolve()
    candidates, _ = read_table(bundle / 'candidates.tsv')
    if candidates and sample_id not in {r['sample_id'] for r in candidates}:
        raise ValueError(f'MS sample {sample_id} is absent from the SNAF export')
    # Require the updated generic manifest selector before executing a search.
    outdir.mkdir(parents=True, exist_ok=True)
    config = dict(run_name='snaf_evidence', output_dir=str(outdir / 'analysis'), sample_id=sample_id,
                  assay_type='immunopeptidomics', source_manifest=str(bundle / 'source_manifest.tsv'),
                  psm_table=str(Path(psm_table).resolve()), search_format=search_format,
                  min_length=7, max_length=25, q_threshold=q_threshold)
    config_path = outdir / 'pyneoquant_config.json'
    config_path.write_text(json.dumps(config, indent=2)+'\n')
    command = [executable_path, 'run-workflow', '--config', str(config_path)]
    with (outdir / 'pyneoquant.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    evidence = outdir / 'analysis/pyneoquant_outputs/SNAF_proteome_peptide_evidence.tsv'
    if not evidence.is_file():
        raise RuntimeError(f'pyNeoQuant did not produce the required evidence contract: {evidence}')
    output = outdir / 'candidates_with_ms.tsv'
    merge_evidence(bundle / 'candidates.tsv', evidence, output, q_threshold)
    manifest(outdir / 'integration_manifest.json', command=command, sample_id=sample_id,
             input_sha256=checksum(psm_table), evidence_sha256=checksum(evidence),
             candidate_sha256=checksum(bundle / 'candidates.tsv'),
             pyneoquant_capabilities=json.loads(capabilities.stdout))
    return output
