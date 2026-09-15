import csv
from pathlib import Path
import pytest
from altanalyze3.components.neoantigen.io import read_table, write_table
from altanalyze3.components.neoantigen.proteomics import export_candidates, merge_evidence, run_pyneoquant
from altanalyze3.components.neoantigen.hla import prepare_hla, read_hla_table
from altanalyze3.components.neoantigen.cli import filter_counts


def test_roundtrip_preserves_event_sample_and_allele(tmp_path):
    table = tmp_path / 'c.tsv'
    table.write_text('sample\tpeptide\tuid\thla\tbinding_affinity\n'
                     'S1\tPEPTIDE\tENSG1:E1-E3\tA0201\t0.1\n'
                     'S1\tPEPTIDE\tENSG2:E1-E3\tA0201\t0.2\n'
                     'S2\tPEPTIDE\tENSG1:E1-E3\tA0201\t0.1\n')
    out = export_candidates([table], tmp_path / 'bundle', predictor='MHCflurry')
    candidates, _ = read_table(out / 'candidates.tsv')
    source = next(c['source_id'] for c in candidates if c['event_id'].startswith('ENSG1:'))
    evidence = tmp_path / 'e.tsv'
    evidence.write_text('sample_id\tpeptide\tsource_id\thla_allele\tq_value\tspectrum_id\n'
                         f'S1\tPEPTIDE\t{source}\tHLA-A*02:01\t0.001\tsp1\n'
                         f'S1\tPEPTIDE\t{source}\tHLA-B*07:02\t0.002\tsp2\n')
    merge_evidence(out / 'candidates.tsv', evidence, tmp_path / 'merged.tsv')
    rows, _ = read_table(tmp_path / 'merged.tsv')
    assert sum(r['ms_status'] == 'peptide_source_supported' for r in rows) == 1
    assert next(r for r in rows if r['ms_status'] == 'peptide_source_supported')['ms_distinct_spectra'] == '1'
    reversed_table = tmp_path / 'reverse.tsv'
    header, *data = table.read_text().splitlines()
    reversed_table.write_text('\n'.join([header]+list(reversed(data)))+'\n')
    other = export_candidates([reversed_table], tmp_path / 'other', predictor='MHCflurry')
    assert (other / 'candidates.tsv').read_bytes() == (out / 'candidates.tsv').read_bytes()


def test_supplied_hla_overrides_bam_and_keeps_no_calls(tmp_path):
    supplied = tmp_path / 'supplied.tsv'
    supplied.write_text('sample\thla\nS1\tA0201,A*02:01,B0702\n')
    output, qc = tmp_path / 'hla.tsv', tmp_path / 'qc.tsv'
    prepare_hla('S1', output, qc, bam='does-not-exist.bam', supplied=supplied)
    assert read_hla_table(output)['S1'] == ['HLA-A*02:01', 'HLA-B*07:02']
    assert read_table(qc)[0][-1]['status'] == 'no_call'
    with pytest.raises(ValueError, match='insufficient'):
        prepare_hla('S1', output, qc, supplied=supplied, require_all=True)


def test_filter_boundary_and_invalid_counts(tmp_path):
    p, o = tmp_path / 'counts.tsv', tmp_path / 'filtered.tsv'
    p.write_text('uid\tS1\tS2\na\t19\t0\nb\t0\t20\nc\t21\t0\n')
    with pytest.raises(ValueError, match='must differ'):
        filter_counts(p, p)
    assert filter_counts(p, o) == 2
    assert [r['uid'] for r in read_table(o)[0]] == ['b', 'c']
    p.write_text('uid\tS1\na\tnan\n')
    with pytest.raises(ValueError, match='Invalid count'):
        filter_counts(p, o)


def test_pyneoquant_missing_is_optional(tmp_path):
    with pytest.raises(RuntimeError, match='Optional pyNeoQuant'):
        run_pyneoquant(tmp_path, 'absent', 'S1', tmp_path, executable='missing-pyneoquant-executable')


def test_optional_external_library_roundtrip(tmp_path):
    import os
    executable = os.environ.get('PYNEOQUANT_EXECUTABLE')
    if not executable:
        pytest.skip('Set PYNEOQUANT_EXECUTABLE to test the separately installed library')
    table = tmp_path / 'c.tsv'
    table.write_text('sample\tpeptide\tuid\thla\nS1\tPEPTIDEK\tENSG1:E1-E3\tA0201\n')
    canonical = tmp_path / 'canonical.fasta'
    canonical.write_text('>canonical\nXXPEPTIDEKXX\n')
    bundle = export_candidates([table], tmp_path / 'bundle', canonical)
    psm = tmp_path / 'psms.tsv'
    psm.write_text('sample_id\tpeptide\tq_value\tspectrum_id\nS1\tPEPTIDEK\t0.001\tsp1\n')
    output = run_pyneoquant(bundle, psm, 'S1', tmp_path / 'ms', executable)
    rows, _ = read_table(output)
    assert rows[0]['ms_status'] == 'peptide_source_supported'
    assert 'canonical' in rows[0]['ms_competing_sources']


def test_reference_archive_rejects_traversal(tmp_path):
    import tarfile, io
    from altanalyze3.components.neoantigen.archive import unpack_reference
    archive = tmp_path / 'bad.tar'
    with tarfile.open(archive, 'w') as tf:
        info=tarfile.TarInfo('../outside');info.size=3;tf.addfile(info,io.BytesIO(b'bad'))
    with pytest.raises(ValueError,match='Unsafe'):
        unpack_reference(archive,tmp_path/'ref')
    assert not (tmp_path/'outside').exists()


def test_standard_gtf_exon_ids_and_custom_gene(tmp_path):
    import pandas as pd
    from altanalyze3.components.gene_model.main import build_gene_model
    for strand in ['+', '-']:
        data = pd.DataFrame([
            ['chr1','gene',100,400,strand,'custom_gene','gene_id "custom_gene";'],
            ['chr1','transcript',100,400,strand,'custom_gene','transcript_id "t";'],
            ['chr1','exon',100,200,strand,'custom_gene','exon_id "exon-one";'],
            ['chr1','exon',300,400,strand,'custom_gene',''],
        ],columns=['chr','type','start','end','strand','gene','attrs'])
        result=build_gene_model(data)
        assert 'custom_gene\tI1.1\tchr1\t'+strand in result
        assert '\t201\t299\t' in result
        assert 'exon-one' in result


def test_hla_genotype_adapter_preserves_inferred_no_call(tmp_path, monkeypatch):
    from altanalyze3.components.bam.bam2hla import bam2hla
    def infer(*args, **kwargs):
        return {'build':'hg38','genes':{
            'HLA-A':{'gene':'HLA-A','call':['A*02:01','A*02:01'],'n_positions':50,'mean_depth':20},
            'HLA-B':{'gene':'HLA-B','call':None,'reason':'insufficient coverage','n_positions':0},
            'HLA-C':{'gene':'HLA-C','call':['C*07:02','C*07:02'],'n_positions':50,'mean_depth':20}}}
    monkeypatch.setattr(bam2hla,'type_bam',infer)
    output,qc=tmp_path/'hla.tsv',tmp_path/'qc.tsv'
    prepare_hla('S1',output,qc,bam='fixture.bam')
    rows,_=read_table(qc)
    assert rows[1]['status']=='no_call' and rows[1]['reason']=='insufficient coverage'
    assert read_hla_table(output)['S1']==['HLA-A*02:01','HLA-C*07:02']


def test_binding_adapter_keeps_sample_allele_pairs_and_rejects_missing(tmp_path, monkeypatch):
    import pandas as pd
    from altanalyze3.components.snaf import binding
    from altanalyze3.components.neoantigen.hla import predict_binding
    evidence, hla, output = tmp_path/'peptides.tsv', tmp_path/'hla.tsv', tmp_path/'binding.tsv'
    evidence.write_text('sample_id\tpeptide\nS1\tPEPTIDEK\nS2\tPEPTIDEK\n')
    hla.write_text('sample\thla\nS1\tA0201,B0702\nS2\tC0702\n')
    def predict(peptides, alleles):
        return pd.DataFrame([dict(peptide=p, hla=h, score=0.5) for p in peptides for h in alleles])
    monkeypatch.setattr(binding, 'run_MHCflurry', predict)
    predict_binding(evidence, hla, output)
    rows, _ = read_table(output)
    assert {(r['sample_id'], r['hla_allele']) for r in rows} == {
        ('S1', 'HLA-A*02:01'), ('S1', 'HLA-B*07:02'), ('S2', 'HLA-C*07:02')}
    assert all(r['hla_score_type'] == 'presentation_percentile' for r in rows)
    monkeypatch.setattr(binding, 'run_MHCflurry', lambda p,h: predict(p,h).iloc[:0])
    with pytest.raises(RuntimeError, match='Incomplete binding'):
        predict_binding(evidence, hla, output)
    monkeypatch.setattr(binding, 'run_MHCflurry', lambda p,h: predict(p,h).assign(score=float('nan')))
    with pytest.raises(ValueError, match='Invalid binding percentile'):
        predict_binding(evidence, hla, output)
