import json
import shutil
from pathlib import Path
from types import SimpleNamespace
import pytest
from altanalyze3.components.neoantigen.galaxy import (
    ANNOTATION_FIELDS, PROFILE, export_galaxy, extract_contract, candidate_reports)
from altanalyze3.components.neoantigen.io import read_table

FIXTURES=Path(__file__).resolve().parents[3]/'deployment/galaxy/test-data'


def export(tmp_path, **kwargs):
    return export_galaxy([FIXTURES/'ipepgen_candidates.tsv'], FIXTURES/'ipepgen_hla.tsv',
                         tmp_path/'out', predictor='MHCflurry', **kwargs)


def test_workflow_extraction_uses_nested_tool_settings(tmp_path):
    profile=json.loads(PROFILE.read_text())
    states={
        'fragpipe': {'wf':{'msfragger':{'digestion':{'digest_min_length':'7','digest_max_length':'25'}}}},
        'pepquery': {'req_inputs':{'input_type':{'input_type_selector':'peptide'}},'ms_params':{'search':{'min_length':'8','max_length':'10'}}},
        'iedb': {'prediction':{'tool':'mhci','alleles':{'allelesrc':'history'},'lengths':['8','10']},'sequence':{'seqsrc':'fasta'}}}
    steps={str(i):{'tool_id':profile['tools'][n]['tool_id'],'tool_state':json.dumps(s)} for i,(n,s) in enumerate(states.items())}
    workflow=tmp_path/'workflow.ga'
    workflow.write_text(json.dumps({'a_galaxy_workflow':'true','name':'contract test',
        'steps':{'1':{'subworkflow':{'steps':steps}}}}))
    actual=extract_contract(workflow)
    assert actual['bounds']['pepquery']==[8,10] and actual['iedb_lengths']==[8,10]
    out=export(tmp_path,workflow=workflow)
    assert (out/'snaf.pepquery.txt').read_text()=='PEPTIDEK\nSYFPEITHI\n'
    states['iedb']['prediction']['tool']='mhcii'
    steps['2']['tool_state']=json.dumps(states['iedb'])
    workflow.write_text(json.dumps({'a_galaxy_workflow':'true','name':'invalid','steps':steps}))
    with pytest.raises(ValueError,match='class-I'):
        extract_contract(workflow)


def test_sample_collections_keep_alleles_and_length_boundaries(tmp_path):
    out=export(tmp_path)
    samples,_=read_table(out/'samples.tsv')
    by_sample={r['sample_id']:r for r in samples}
    assert len(samples)==2
    assert (out/by_sample['S1']['iedb_alleles']).read_text()=='HLA-A*02:01\nHLA-B*07:02\n'
    assert (out/by_sample['S2']['iedb_alleles']).read_text()=='HLA-B*07:02\n'
    assert (out/by_sample['S1']['pepquery']).read_text()=='SYFPEITHI\n'
    assert 'PEPTIDEK' in (out/by_sample['S1']['iedb_fasta']).read_text()
    assert 'PEPTIDEK' not in (out/by_sample['S2']['iedb_fasta']).read_text()
    mapping,_=read_table(out/'candidate_map.tsv')
    assert len({r['accession'] for r in mapping if r['peptide']=='SYFPEITHI'})==2
    titles=[l[1:].split('|')[1] for l in (out/'snaf.database.fasta').read_text().splitlines() if l.startswith('>')]
    assert set(titles)=={r['accession'] for r in mapping}
    assert (out/'peptides.bed').read_text()==''
    annotations,fields=read_table(out/'annotations.tsv')
    assert fields==ANNOTATION_FIELDS
    assert all(r['Start']=='' and 'coordinates_unavailable' in r['Annotation'] for r in annotations)
    assert all(r['coord'] for r in mapping)  # event context remains separate
    provenance=json.loads((out/'integration.json').read_text())
    assert provenance['contract']['bounds']['pepquery']==[9,11]
    assert provenance['mapped_sources']==0


def test_explicit_peptide_bed_has_split_blocks_and_correct_browser_coordinates(tmp_path):
    initial=export(tmp_path/'initial')
    rows,_=read_table(initial/'candidate_map.tsv')
    source=next(r for r in rows if r['peptide']=='SYFPEITHI')
    bed=tmp_path/'peptides.bed'
    bed.write_text(f"chr1\t99\t215\t{source['accession']}\t0\t-\t99\t215\t0\t2\t12,15\t0,101\n")
    out=export(tmp_path/'mapped',peptide_bed=bed)
    rows,_=read_table(out/'annotations.tsv')
    row=next(r for r in rows if r['Start'])
    assert row['Start']=='99' and row['End']=='215' and row['Strand']=='-'
    assert row['IGV_Genome_Coordinate']=='chr1:100-215'
    assert '100-215' in row['UCSC_Genome_Browser']
    assert (out/'peptides.bed').read_text()==bed.read_text()
    bed.write_text(bed.read_text().replace('12,15','12,14'))
    with pytest.raises(ValueError,match='coding blocks'):
        export(tmp_path/'invalid',peptide_bed=bed)
    assert not (tmp_path/'invalid/out').exists()


def test_missing_sample_calls_rejected_and_empty_candidates_supported(tmp_path):
    hla=tmp_path/'hla.tsv';hla.write_text('sample\thla\nS1\tA0201\n')
    with pytest.raises(ValueError,match='Missing HLA samples'):
        export_galaxy([FIXTURES/'ipepgen_candidates.tsv'],hla,tmp_path/'bad')
    candidates=tmp_path/'empty.tsv';candidates.write_text('sample\tpeptide\tuid\thla\n')
    out=export_galaxy([candidates],hla,tmp_path/'empty')
    assert (out/'snaf.database.fasta').read_text()==''
    assert len(read_table(out/'samples.tsv')[0])==1
    assert read_table(out/'candidate_map.tsv')[0]==[]
    with pytest.raises(ValueError,match='must be empty'):
        export_galaxy([candidates],hla,out)


def test_report_selection_and_snaf_exports_are_opt_in(tmp_path):
    from altanalyze3.components.snaf.cli import _export_snaf_outputs
    _export_snaf_outputs(SimpleNamespace(),tmp_path)
    assert list(tmp_path.iterdir())==[]
    directory=tmp_path/'T_candidates';directory.mkdir()
    combined=directory/'T_antigen_candidates_all.txt'
    shutil.copyfile(FIXTURES/'ipepgen_candidates.tsv',combined)
    (directory/'T_antigen_candidates_S1.txt').write_text('intentionally invalid duplicate\n')
    assert candidate_reports(tmp_path)==[combined]
    _export_snaf_outputs(SimpleNamespace(galaxy_integration=True,hla=FIXTURES/'ipepgen_hla.tsv',binding_method='MHCflurry'),tmp_path)
    assert (tmp_path/'galaxy_export/snaf.database.fasta').exists()
    assert not (tmp_path/'proteomics_export').exists()


def test_stable_ids_with_sample_names_unsafe_for_paths(tmp_path):
    table=tmp_path/'c.tsv';hla=tmp_path/'h.tsv'
    table.write_text('sample\tpeptide\tuid\thla\n../A B\tSYFPEITHI\tENSG1:E1-E3\tA0201\n')
    hla.write_text('sample\thla\n../A B\tA0201\n')
    out=export_galaxy([table],hla,tmp_path/'out')
    row=read_table(out/'samples.tsv')[0][0]
    assert row['sample_id']=='../A B'
    assert '/' not in row['collection_id'] and ' ' not in row['collection_id']
    assert (out/row['iedb_alleles']).is_file()


def test_three_digit_hla_fields_and_sample_subset(tmp_path):
    from altanalyze3.components.neoantigen.hla import normalize_class_i
    assert normalize_class_i('B*15:153') == 'HLA-B*15:153'
    assert normalize_class_i('HLA-A*02:101:01') == 'HLA-A*02:101'
    assert normalize_class_i('A0201') == 'HLA-A*02:01'
    with pytest.raises(ValueError):
        normalize_class_i('A02101')
    table=tmp_path/'empty.tsv';hla=tmp_path/'h.tsv'
    table.write_text('sample\tpeptide\tuid\thla\n')
    hla.write_text('sample\thla\nS1\tB*15:153\nS2\tA0201\n')
    out=export_galaxy([table],hla,tmp_path/'out',sample_ids=['S1'])
    samples,_=read_table(out/'samples.tsv')
    assert len(samples)==1 and samples[0]['sample_id']=='S1'
    assert (out/samples[0]['iedb_alleles']).read_text()=='HLA-B*15:153\n'
