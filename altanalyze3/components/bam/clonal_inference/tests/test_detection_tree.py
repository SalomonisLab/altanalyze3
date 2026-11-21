import unittest,tempfile,json
from pathlib import Path
import numpy as np,pandas as pd
from altanalyze3.components.bam.clonal_inference.__main__ import main
from altanalyze3.components.bam.clonal_inference.detection_workflow import specification
from altanalyze3.components.bam.clonal_inference.integrity import sha256
from altanalyze3.components.bam.clonal_inference.detection_tree import evidence,genotypes,TreeSearch

class DetectionTests(unittest.TestCase):
 def fixture(self,p):
  X=np.repeat([[0,0,0],[1,1,0],[1,1,1]],[8,30,8],axis=0).astype(np.int8);ids=np.array([str(i) for i in range(len(X))]);features=np.array(['a','b','c']);np.savez(p/'calls.npz',X=X,A=X*5,T=np.full_like(X,5),cell_ids=ids,features=features)
  np.savez(p/'detector.npz',Z=X,available=np.ones_like(X,bool),cell_ids=ids,features=features)
  pd.DataFrame(dict(cell_id=ids,sample='s',patient='p')).to_csv(p/'cells.tsv',sep='\t',index=False);pd.DataFrame(dict(feature=features,event_label=features,chrom='chr1')).to_csv(p/'events.tsv',sep='\t',index=False);(p/'independent.txt').write_text('fixture independently defined complete roster a/b/c and 46 cells')
  con=dict(schema_version=1,method='detection-tree',algorithm=specification(sex='F'),source_precedence='fixture',authorization='test',expected_cell_ids=ids.tolist(),expected_features=features.tolist(),expected_samples=['s']*len(X),input_sha256={k:sha256(p/f) for k,f in [('matrix','calls.npz'),('cells','cells.tsv'),('events','events.tsv')]},detector_sha256=sha256(p/'detector.npz'),provenance_sources=[dict(path='independent.txt',sha256=sha256(p/'independent.txt'))]);(p/'contract.json').write_text(json.dumps(con))
  return ['--matrix',str(p/'calls.npz'),'--cells',str(p/'cells.tsv'),'--events',str(p/'events.tsv'),'--contract',str(p/'contract.json'),'--detector',str(p/'detector.npz'),'--sex','F','--output',str(p/'out')]
 def test_default_and_complete_roster(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);argv=self.fixture(p);main(argv);s=json.loads((p/'out/summary.json').read_text());self.assertEqual(s['seed'],222);self.assertEqual(s['cells'],46);f=pd.read_csv(p/'out/clone_frequencies.tsv',sep='\t');self.assertAlmostEqual(f.posterior_expected_percent_all.sum(),100);self.assertIn('WT',f.composition.tolist())
   with self.assertRaisesRegex(ValueError,'overwrite'):main(argv)
 def test_detector_tamper_blocks_before_output(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);argv=self.fixture(p);(p/'detector.npz').write_bytes(b'bad')
   with self.assertRaisesRegex(ValueError,'Detector hash'):main(argv)
   self.assertFalse((p/'out').exists())
 def test_algorithm_change_requires_contract(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);argv=self.fixture(p)
   with self.assertRaisesRegex(ValueError,'algorithm differs'):main(argv+['--pooling-strength','10'])
   self.assertFalse((p/'out').exists())
 def test_no_inheritance_loss(self):
  g=genotypes([-1,0,1,2]);self.assertTrue(np.all(g[3]>=g[2]));self.assertEqual(g[3].sum(),3)
 def test_low_depth_reference_is_uncertain(self):
  a,b=evidence(np.array([[0],[-1]]),np.zeros((2,1)),np.array([[1],[0]]),[.5]);self.assertAlmostEqual(np.exp(b[0,0]),.5);self.assertEqual(a[1,0],b[1,0])

class DetectionFigureTests(unittest.TestCase):
 def test_fish_preserves_knots_and_total(self):
  from altanalyze3.components.bam.clonal_inference.detection_plotting import interpolate_exclusive
  f=np.array([[80,20,0],[15,70,95],[5,10,5.]])
  x,v=interpolate_exclusive(f,points=11)
  np.testing.assert_allclose(v[:,[0,10,20]],f)
  np.testing.assert_allclose(v.sum(0),100)
  self.assertTrue((v>=0).all())
  with self.assertRaisesRegex(ValueError,'sum to 100'):interpolate_exclusive(f[:2])
 def test_complete_vector_figures_and_source_gate(self):
  from altanalyze3.components.bam.clonal_inference.detection_plotting import render
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);argv=DetectionTests().fixture(p)
   c=pd.read_csv(p/'cells.tsv',sep='\t',dtype=str);c['stage']='Diagnosis';c['cell_state']='State1';c['stage_order']=0;c['is_xenograft']=False;c.to_csv(p/'cells.tsv',sep='\t',index=False)
   con=json.loads((p/'contract.json').read_text());con['input_sha256']['cells']=sha256(p/'cells.tsv');(p/'contract.json').write_text(json.dumps(con));main(argv)
   z=np.load(p/'calls.npz');np.savez(p/'source_evidence.npz',matrix=(z['X']+1).T,cell_ids=z['cell_ids'],features=z['features'])
   pd.DataFrame(dict(event=['a','b','c'],label=['a | direct','b | direct','c | direct'],channel=['Direct reads']*3)).to_csv(p/'evidence_rows.tsv',sep='\t',index=False)
   pc=dict(sha256={n:sha256(p/n) for n in ['source_evidence.npz','evidence_rows.tsv']},verified_source_provenance='independent synthetic fixture')
   (p/'plot_source_contract.json').write_text(json.dumps(pc));checks=render(p/'out',p,p/'figures');self.assertEqual(checks['cells'],46);self.assertEqual(len(list((p/'figures').glob('*.svg'))),7)
   (p/'evidence_rows.tsv').write_text('tampered')
   with self.assertRaisesRegex(ValueError,'Plot source hash'):render(p/'out',p,p/'should_not_exist')
   self.assertFalse((p/'should_not_exist').exists())
