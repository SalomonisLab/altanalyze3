"""Integrity tests for the default publication adapters (no scientific inference changes)."""
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from altanalyze3.components.visualization.clone_alluvial import render
from altanalyze3.components.visualization.clone_umap_mapping import prepare
from altanalyze3.components.visualization.clone_figure_workflow import run_manifest

class CloneFigureIntegrityTests(unittest.TestCase):
    def fixture(self, root):
        cells=pd.DataFrame(dict(cell_id=['c1','c2','c3'],patient=['p']*3,sample=['s1','s2','xeno'],stage=['Diagnosis','Post','Graft'],stage_order=[0,1,0],is_xenograft=[False,False,True],cell_state=['new state','ERP-1','other'],node=[0,1,1]))
        cells.to_csv(root/'cell_assignments.tsv.gz',sep='\t',index=False)
        pd.DataFrame(dict(node=[0,1],clone=['WT','SC1'],color=['#cccccc','#ff0000'],mutation_combination=['No direct or imputed mutant call','MUT1'])).to_csv(root/'clone_key.tsv',sep='\t',index=False)
        pd.DataFrame(dict(sample=['s1','s2','xeno'],clone=['WT','SC1','SC1'],cells=[1,1,1])).to_csv(root/'clone_frequencies.tsv',sep='\t',index=False)
        return cells

    def test_alluvial_cannot_silently_exclude_clinical_specimen(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);self.fixture(p);(p/'config.json').write_text(json.dumps({'samples':['s1']}))
            with self.assertRaisesRegex(ValueError,'clinical specimens'):
                render(p,p/'out',p/'config.json')

    def test_unknown_state_retained_and_graft_audited_separately(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);self.fixture(p)
            with patch('matplotlib.figure.Figure.savefig'):
                result=render(p,p/'out')
            self.assertEqual(result['clinical_cells'],2)
            self.assertEqual(result['MISTRG_separate_from_clinical_alluvial'],1)
            mapped=pd.read_csv(p/'out/cell_state_aggregation.tsv',sep='\t')
            self.assertIn('new state',set(mapped.cell_state))
            self.assertEqual(pd.read_csv(p/'out/clinical_cell_assignments.tsv.gz',sep='\t').cell_id.tolist(),['c1','c2'])

    def test_corrupt_manifest_blocks_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);(p/'input').write_text('changed');(p/'manifest.json').write_text(json.dumps({'patients':{},'sha256':{'input':'incorrect'}}))
            with self.assertRaises(Exception):run_manifest(p/'manifest.json',p/'out')

    def test_approved_preset_keeps_mpp_mep_separate_and_merges_eryth(self):
        p=Path(__file__).parents[1]/'presets/kinnex_5801_broad_states.json';c=json.loads(p.read_text())
        self.assertEqual(c['state_map']['MPP-MEP'],'MPP-MEP')
        for s in ['MEP-Eryth-1','MEP-Eryth-2','Erythroblast-3']:
            self.assertEqual(c['state_map'][s],'ERP')
        self.assertEqual(c['state_colors']['MEP'],'#ff0000')

if __name__=='__main__':unittest.main()
