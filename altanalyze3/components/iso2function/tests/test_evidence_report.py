import json
from types import SimpleNamespace
import pytest
from altanalyze3.components.iso2function.protein_impact import write
from altanalyze3.components.iso2function.evidence_report import run_report


def test_report_retains_mismatched_pair_without_interaction_differences(tmp_path):
    b=tmp_path/'assay';b.mkdir()
    (b/'assay_match_manifest.json').write_text(json.dumps(dict(status='complete')))
    write(b/'all_native_pair_evidence.tsv',[dict(gene='G',both_native_proteins_exact='False',native_A='A',native_B='B')])
    write(b/'native_assayed_protein_comparison.tsv',[dict(gene='G')])
    write(b/'amino_acid_differences.tsv',[dict(gene='G',native='A',native_sequence='<C>',assayed_sequence='Y')])
    args=SimpleNamespace(impact_results=None,query_results=None,assay_results=str(b),artifacts_dir=str(tmp_path/'report'))
    assert run_report(args)==0
    text=(tmp_path/'report/iso2function_evidence.html').read_text()
    assert '<h2>G</h2>' in text and '&lt;C&gt;' in text
    assert '<C>' not in text
    assert json.loads((tmp_path/'report/report_manifest.json').read_text())['genes']==1
    args.artifacts_dir=str(tmp_path/'other')
    (b/'assay_match_manifest.json').write_text(json.dumps(dict(status='incomplete')))
    with pytest.raises(ValueError,match='Incomplete'):
        run_report(args)


def test_new_query_only_gene_shows_the_specific_substitution(tmp_path):
    b=tmp_path/'query';b.mkdir()
    (b/'query_manifest.json').write_text(json.dumps(dict(status='complete',parameters=dict(min_isoform_reads=30))))
    write(b/'all_pair_decisions.tsv',[dict(gene='NEW',native_supported='True',clone_A='NEW-1',clone_B='NEW-2',
        native_aliases_A=json.dumps([dict(isoform='A',reads=40,match_type='exact',substitutions=[])]),
        native_aliases_B=json.dumps([dict(isoform='B',reads=50,match_type='AA_substitutions',substitutions=[dict(position=266,native='W',assayed='C')])]),
        native_evidence_A='exact',native_evidence_B='AA_substitutions',native_pair_evidence='putative_transfer_with_AA_substitutions')])
    write(b/'divergent_interactions.tsv',[dict(gene='NEW',modality='pdi',target='DNA_BAIT')])
    args=SimpleNamespace(impact_results=None,query_results=str(b),assay_results=None,artifacts_dir=str(tmp_path/'report'))
    assert run_report(args)==0
    text=(tmp_path/'report/iso2function_evidence.html').read_text()
    assert '266: W to C' in text and 'putative_transfer_with_AA_substitutions' in text
    assert 'DNA_BAIT' in text
