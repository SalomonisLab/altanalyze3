"""Optional independent checks of saved real-data CLI figure exports."""
import json
import os
from pathlib import Path

import pytest

from altanalyze3.components.iso2function.protein_impact import sha,table


def test_verified_figures_reproduce_all_five_claims_and_export_bytes():
    root=os.environ.get('ISO2FUNCTION_PROTEIN_FIGURES')
    repeat=os.environ.get('ISO2FUNCTION_PROTEIN_FIGURES_REPEAT')
    if not root or not repeat:
        pytest.skip('Set ISO2FUNCTION_PROTEIN_FIGURES and ISO2FUNCTION_PROTEIN_FIGURES_REPEAT')
    root,repeat=Path(root),Path(repeat)
    assert {p.name:sha(p) for p in root.iterdir()}=={p.name:sha(p) for p in repeat.iterdir()}
    m=json.loads((root/'figure_manifest.json').read_text())
    assert m['status']=='complete' and m['pairs']==m['pages']==5
    assert m['genes']==['TOX','CITED2','HOXB5','NR4A2','LYL1']
    assert not m['generative_model_used'] and not m['analysis_changed']
    assert m['projection_disagreements']==0
    assert m['selection']['all_other_annotations_preserved_in_display_audit']
    for name,digest in m['outputs'].items():assert sha(root/name)==digest
    rows=table(root/'plotted_feature_coordinates.tsv')
    def feature(gene,start,end):
        found=[r for r in rows if r['gene_symbol']==gene and (r['A_start'],r['A_end'])==(str(start),str(end))]
        assert len(found)==1
        return found[0]
    for gene,start,end,n in [('TOX',261,329,69),('HOXB5',176,181,6),('CITED2',216,269,54)]:
        r=feature(gene,start,end)
        assert r['status']=='absent' and int(r['absent_residues'])==n
    assert feature('CITED2',216,269)['provider']=='literature_construct_snapshot'
    for gene,start,end,b_start,b_end in [('HOXB5',194,253,6,65),('LYL1',150,202,76,128),('NR4A2',260,335,260,335)]:
        r=feature(gene,start,end)
        assert r['status']=='retained' and int(r['identical_residues'])==end-start+1
        assert (int(r['B_native_start']),int(r['B_native_end']))==(b_start,b_end)
    r=feature('NR4A2',360,595)
    assert r['status']=='partial' and (int(r['identical_residues']),int(r['feature_length']))==(95,236)
    assert (r['B_native_start'],r['B_native_end'])==('360','454')
    assert len(list(root.glob('*.pdf')))==6
    assert len(list(root.glob('*.svg')))==len(list(root.glob('*.png')))==5
    layouts=json.loads((root/'label_layout_audit.json').read_text())
    assert len(layouts)==5 and all(page['labels'] for page in layouts)
    if m.get('plot_style')=='compact':
        assert m['presentation']==dict(amino_acid_substitutions_highlighted=False,
                                       narrative_footer=False,domain_colors='annotation identity')
        native=table(root/'native_architecture_coordinates.tsv')
        assert len(native)==44
        def native_feature(gene,description):
            found=[r for r in native if r['gene_symbol']==gene and r['side']=='B' and r['description']==description]
            assert len(found)==1
            return found[0]
        for gene,desc,start,end in [('HOXB5','Homeobox',6,65),('LYL1','bHLH',76,128),('NR4A2','NR LBD',360,454)]:
            row=native_feature(gene,desc)
            assert (row['native_start'],row['native_end'])==(str(start),str(end))
        assert native_feature('NR4A2','NR LBD')['placement']=='partial'
        baseline=root.parent/'figures_v2'
        if baseline.exists():
            for name in ['plot_data.json','annotation_display_audit.tsv','alignment_events.tsv','plotted_feature_coordinates.tsv']:
                assert sha(root/name)==sha(baseline/name),name
        if m.get('coding_exon_map'):
            drawn=table(root/'coding_exon_coordinates.tsv')
            assert {r['gene_symbol'] for r in drawn}==set(m['genes'])
            assert len({r['isoform_id'] for r in drawn})==10
            from collections import defaultdict
            lengths=defaultdict(float)
            for r in drawn:lengths[r['isoform_id']]+=float(r['aa_end0'])-float(r['aa_start0'])
            expected={'ENST00000361421':526,'6896443.WF40_ND21-251_HSC_3k__source000':148,
                      'ENST00000367651':270,'8515942.BM27_CD34':135,
                      'ENST00000239151':269,'1995169.WF40_ND21-251_HSC_3k__source000':81,
                      'ENST00000339562':598,'41458418.AM72-N-0206_CD34':454,
                      'ENST00000264824':280,'82405157.AML-7_1_CITE_GEX__source020':206}
            assert lengths==pytest.approx(expected)
            if 'exon_label' in drawn[0]:
                observed=defaultdict(list)
                for row in drawn:
                    assert row['exon_label_status']=='assigned'
                    observed[row['isoform_id']].append(row['exon_label'])
                assert observed['ENST00000361421']==['E'+str(n) for n in range(1,10)]
                assert observed['6896443.WF40_ND21-251_HSC_3k__source000']==['E1','E2','E3/I3']
                assert observed['ENST00000367651']==['E4']
                assert observed['8515942.BM27_CD34']==['E4','E4']
                assert observed['ENST00000339562']==['E5','E6','E7','E8','E9','E10']
                assert observed['41458418.AM72-N-0206_CD34']==['E5','E6','E7','E8/I8']
                assert observed['ENST00000239151']==['E1','E2']
                assert observed['1995169.WF40_ND21-251_HSC_3k__source000']==['E2']
                for tid in ['ENST00000264824','82405157.AML-7_1_CITE_GEX__source020']:
                    assert observed[tid]==['E2','E3','E4']
                for page in layouts:
                    assert any(r['role'].startswith('exon_label_') for r in page['labels'])
