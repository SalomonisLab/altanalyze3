import unittest, tempfile, json
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from altanalyze3.components.visualization.clone_alluvial import render, add_arguments
import argparse
class UnassignedTests(unittest.TestCase):
 def fixture(self,p):
  a=pd.DataFrame(dict(cell_id=['a','b','u','c','v'],patient=['p']*5,sample=['s1']*3+['s2']*2,stage=['Pre']*3+['Post']*2,stage_order=[0]*3+[1]*2,is_xenograft=[False]*5,cell_state=['HSC']*5,node=[0,1,-1,1,-1]))
  a.to_csv(p/'cell_assignments.tsv.gz',sep='\t',index=False)
  pd.DataFrame(dict(node=[0,1,-1],clone=['WT','SC1','U'],color=['#cccccc','#ff0000','#999999'],mutation_combination=['No direct or imputed mutant call','MUT1','Unresolved'])).to_csv(p/'clone_key.tsv',sep='\t',index=False)
  a['clone']=a.node.map({0:'WT',1:'SC1',-1:'U'})
  a.groupby(['sample','clone']).size().rename('cells').reset_index().to_csv(p/'clone_frequencies.tsv',sep='\t',index=False)
 def test_both_modes_and_accounting(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);self.fixture(p)
   for include in [False,True]:
    o=p/str(include)
    with patch('matplotlib.figure.Figure.savefig'), patch('altanalyze3.components.bam.clonal_inference.detection_plotting.require_arial'):
     v=render(p,o,include_unassigned=include)
    d=pd.read_csv(o/'clinical_cell_assignments.tsv.gz',sep='\t')
    self.assertEqual(len(d),5 if include else 3)
    self.assertEqual(set(d.clone),{'WT','SC1','U'} if include else {'WT','SC1'})
    self.assertEqual(v['excluded_unassigned_cells'],0 if include else 2)
    f=pd.read_csv(o/'cell_state_clone_frequencies.tsv',sep='\t');col='percent_all_cells' if include else 'percent_assigned_cells'
    self.assertTrue((f.groupby('sample')[col].sum().round(8)==100).all())
    if not include:self.assertEqual(f.loc[(f['sample']=='s1')&(f.clone=='WT'),col].iloc[0],50)
 def test_empty_sample_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);self.fixture(p)
   a=pd.read_csv(p/'cell_assignments.tsv.gz',sep='\t');a.loc[a['sample']=='s2','node']=-1;a.to_csv(p/'cell_assignments.tsv.gz',sep='\t',index=False)
   a['clone']=a.node.map({0:'WT',1:'SC1',-1:'U'});a.groupby(['sample','clone']).size().rename('cells').reset_index().to_csv(p/'clone_frequencies.tsv',sep='\t',index=False)
   with patch('altanalyze3.components.bam.clonal_inference.detection_plotting.require_arial'),self.assertRaisesRegex(ValueError,'No assigned cells'):render(p,p/'out')
 def test_default_cli(self):
  p=argparse.ArgumentParser();add_arguments(p)
  self.assertFalse(p.parse_args(['--figures','x','--output','y']).include_unassigned)
  self.assertTrue(p.parse_args(['--figures','x','--output','y','--include-unassigned']).include_unassigned)
if __name__=='__main__':unittest.main()
