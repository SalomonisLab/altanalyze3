"""Class-I HLA input adapters. No pyNeoQuant dependency."""
import importlib.metadata
import math
import re
from pathlib import Path
from .io import read_table, write_table, manifest


def normalize_class_i(value):
    value = str(value).strip().upper().removeprefix('HLA-')
    if not value:
        return ''
    m = re.fullmatch(r'([ABC])\*?(\d{2,3}):(\d{2,3})(?::\d{2,3})*(?:[NLSCAQ])?', value)
    if not m:
        m = re.fullmatch(r'([ABC])\*?(\d{2})(\d{2})(?:[NLSCAQ])?', value)
    if not m:
        raise ValueError(f'Unsupported class-I allele {value!r}; expected HLA-A/B/C*NN:NN')
    if value[-1:].isalpha():
        raise ValueError(f'Expression-status HLA suffix requires explicit review: {value}')
    return f'HLA-{m[1]}*{m[2]}:{m[3]}'


def read_hla_table(path):
    rows, fields = read_table(path)
    if len(fields) < 2:
        raise ValueError('HLA table needs sample and allele columns')
    result = {}
    for row in rows:
        sample = row[fields[0]]
        if not sample or sample in result:
            raise ValueError(f'Empty or duplicate HLA sample: {sample!r}')
        tokens = re.split('[,;\\s]+', ' '.join(row[f] for f in fields[1:]))
        result[sample] = list(dict.fromkeys(normalize_class_i(t) for t in tokens if t and t.upper() not in ('NA', 'NAN', 'NONE')))
    return result


def prepare_hla(sample, output, qc, bam=None, supplied=None, build='auto', min_depth=8, require_all=False):
    """Supplied alleles override inference; empty/no-call supplied rows can infer."""
    values = read_hla_table(supplied).get(sample, []) if supplied else []
    quality = []
    if values:
        for gene in ('A', 'B', 'C'):
            alleles = [a for a in values if a.startswith(f'HLA-{gene}*')]
            quality.append(dict(sample_id=sample, gene=f'HLA-{gene}', origin='supplied',
                                status='called' if alleles else 'no_call', alleles=','.join(alleles)))
    else:
        if not bam:
            raise ValueError(f'{sample}: no supplied HLA alleles and no BAM for inference')
        from altanalyze3.components.bam.bam2hla.bam2hla import type_bam
        typed = type_bam(str(bam), build=build, min_depth=min_depth, verbose=False)
        for gene, result in typed['genes'].items():
            alleles = [normalize_class_i(a) for a in result.get('call') or []]
            values.extend(alleles)
            quality.append({**result, 'sample_id': sample, 'gene': gene, 'origin': 'bam2hla',
                            'build': typed['build'], 'status': 'called' if alleles else 'no_call',
                            'alleles': ','.join(alleles)})
    values = list(dict.fromkeys(values))
    qc_fields = ['sample_id', 'gene', 'origin', 'build', 'status', 'alleles', 'reason',
                 'n_positions', 'mean_depth', 'explained_frac', 'minor_support']
    write_table(qc, quality, qc_fields)
    if not values or (require_all and any(r['status'] == 'no_call' for r in quality)):
        raise ValueError(f'{sample}: insufficient HLA calls; inspect {qc}')
    write_table(output, [{'sample': sample, 'hla': ','.join(values)}], ['sample', 'hla'])
    return values


def combine_hla(paths, output):
    combined = {}
    for path in paths:
        for sample, alleles in read_hla_table(path).items():
            if sample in combined:
                raise ValueError(f'Duplicate HLA sample: {sample}')
            if not alleles:
                raise ValueError(f'No HLA calls for {sample}')
            combined[sample] = alleles
    write_table(output, [{'sample': s, 'hla': ','.join(combined[s])} for s in sorted(combined)], ['sample', 'hla'])


def predict_binding(evidence, hla, output, method='MHCflurry', software_path=None):
    """Emit all sample/peptide/allele pairs for the pyNeoQuant prediction contract."""
    from altanalyze3.components.snaf.binding import run_MHCflurry, run_netMHCpan
    if method not in ('MHCflurry', 'netMHCpan'):
        raise ValueError(f'Unsupported binding method: {method}')
    rows, _ = read_table(evidence)
    sample_hla = read_hla_table(hla)
    groups = {}
    for row in rows:
        sample = row.get('sample_id') or row.get('sample')
        peptide = row.get('peptide', '')
        if sample not in sample_hla or not sample_hla[sample]:
            raise ValueError(f'Missing HLA alleles for sample {sample}')
        if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', peptide):
            raise ValueError(f'Invalid unmodified peptide {peptide!r}')
        groups.setdefault(sample, set()).add(peptide)
    output_rows = []
    try:
        version = importlib.metadata.version('mhcflurry') if method == 'MHCflurry' else 'user-installed'
    except importlib.metadata.PackageNotFoundError:
        version = 'unknown'
    for sample, peptides in sorted(groups.items()):
        if method == 'MHCflurry':
            tables = [run_MHCflurry(sorted(peptides), sample_hla[sample])]
        else:
            if not software_path:
                raise ValueError('netMHCpan requires --software-path')
            tables = [run_netMHCpan(software_path, sorted(p for p in peptides if len(p) == n), sample_hla[sample], n)
                      for n in sorted({len(p) for p in peptides})]
        for table in tables:
            for record in table.to_dict('records'):
                rank = float(record['score'])
                if not math.isfinite(rank) or not 0 <= rank <= 100:
                    raise ValueError(f'Invalid binding percentile: {rank}')
                output_rows.append(dict(sample_id=sample, peptide=record['peptide'], hla_allele=record['hla'],
                    hla_class='I', hla_binding_rank=rank, hla_predictor=method,
                    hla_predictor_version=version,
                    hla_score_type='presentation_percentile' if method == 'MHCflurry' else 'el_percentile'))
    expected = {(s, p, h) for s, peptides in groups.items() for p in peptides for h in sample_hla[s]}
    observed = {(r['sample_id'], r['peptide'], normalize_class_i(r['hla_allele'])) for r in output_rows}
    if expected != observed:
        raise RuntimeError(f'Incomplete binding predictions: {len(expected-observed)} missing pairs, {len(observed-expected)} unexpected pairs')
    write_table(output, output_rows, ['sample_id', 'peptide', 'hla_allele', 'hla_class', 'hla_binding_rank',
                                    'hla_predictor', 'hla_predictor_version', 'hla_score_type'])
