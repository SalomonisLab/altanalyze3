import unittest,tempfile,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np,pandas as pd
from altanalyze3.components.bam.clonal_inference.shared import fit,expectation
from altanalyze3.components.bam.clonal_inference.shared_workflow import merge_profiles,mask_calls,run_shared

class SharedTests(unittest.TestCase):
 def test_missing_and_exact_likelihood(self):
  X=np.array([[1,0,-1],[-1,-1,-1]])
  p=np.array([[.9,.1,.7],[.1,.9,.3]]);w=np.array([[.3,.7]])
  r,n=expectation(X,np.zeros(2,int),p,w);q=np.array([.3*.9*.9,.7*.1*.1])
  np.testing.assert_allclose(r[0],q/q.sum());np.testing.assert_allclose(r[1],w[0]);self.assertAlmostEqual(n[1],0)
 def test_duplicate_and_unobserved_marginalization(self):
  f=dict(G=np.array([[1,0,0],[1,0,1],[0,0,1]]),posterior=np.array([[.2,.3,.5],[.4,.3,.3]]),w=np.array([[.2,.3,.5]]))
  g,q,w=merge_profiles(f,np.array([True,True,False]));self.assertEqual(g.tolist(),[[0,0,-1],[1,0,-1]])
  np.testing.assert_allclose(q,[[.5,.5],[.3,.7]]);np.testing.assert_allclose(w,[[.5,.5]])
  perm=[2,0,1];f2=dict(G=f['G'][perm],posterior=f['posterior'][:,perm],w=f['w'][:,perm]);g2,q2,w2=merge_profiles(f2,np.array([True,True,False]))
  np.testing.assert_array_equal(g,g2);np.testing.assert_allclose(q,q2)
 def test_known_clones_and_error_caps(self):
  G=np.array([[0,0,0],[1,1,0],[1,1,1]])
  X=np.repeat(G,[100,300,100],axis=0);s=np.zeros(len(X),int)
  fits=[fit(X,s,3,'learned',seed,fp_cap=.1,fn_cap=.3,maxiter=10000) for seed in range(5)]
  f=max(fits,key=lambda z:z['objective']);self.assertEqual({tuple(g) for g in f['G']},{tuple(g) for g in G})
  self.assertTrue((f['errors'][0]<=.1).all() and (f['errors'][1]>=.7).all());self.assertTrue(f['converged'])
 def test_grouped_mask_no_cross_event_leakage(self):
  X=np.ones((100,3),int);events=pd.DataFrame(dict(feature=['RUNX1','LOH','other'],validation_group=['RUNX1','RUNX1','other']))
  m=mask_calls(X,events,7);np.testing.assert_array_equal(m[:,0],m[:,1]);self.assertTrue((m[:,0]!=m[:,2]).any())
 def test_full_roundtrip_preserves_unknown_event_wt_and_no_evidence(self):
  G=np.array([[0,0,-1],[1,0,-1],[1,1,-1]])
  X=np.repeat(G,[80,100,80],axis=0);X=np.vstack([X,[-1,-1,-1]])
  cells=pd.DataFrame(dict(cell_id=[f'c{i}' for i in range(len(X))],sample=['s']*len(X),patient=['p']*len(X)))
  events=pd.DataFrame(dict(feature=['a','b','never_covered']))
  args=SimpleNamespace(seed=7,starts=2,k_grid=[2,3,4],validation_seed=17,max_iterations=10000,fp_cap=.1,fn_cap=.3)
  with tempfile.TemporaryDirectory() as folder:
   root=Path(folder);selection=run_shared(X,cells,events,root,args,{})
   out=pd.read_csv(root/'cell_assignments.tsv.gz',sep='\t');self.assertEqual(out.cell_id.tolist(),cells.cell_id.tolist());self.assertEqual(out.clone.iloc[-1],'NE')
   geno=pd.read_csv(root/'clone_genotypes.tsv',sep='\t');self.assertTrue((geno.loc[geno.feature=='never_covered','state']==-1).all())
   freq=pd.read_csv(root/'clone_frequencies.tsv',sep='\t');self.assertAlmostEqual(freq.percent_all.sum(),100);self.assertEqual(freq.cells.sum(),len(X))
   self.assertEqual(selection['globally_unobserved_features'],['never_covered'])
   self.assertTrue(out.genotype.str.endswith('?').all())

if __name__=='__main__':unittest.main()
