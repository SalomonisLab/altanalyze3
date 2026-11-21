import unittest
import numpy as np
from altanalyze3.components.bam.clonal_inference.fishplot import smooth_composition

class FishTests(unittest.TestCase):
 def test_exact_knots_and_all_cell_mass(self):
  values=np.array([[0.,40,10],[75,10,60],[25,50,30]])
  x,y=smooth_composition(values)
  self.assertTrue((y>=0).all());np.testing.assert_allclose(y.sum(0),100)
  for j in range(values.shape[1]):np.testing.assert_allclose(y[:,np.argmin(abs(x-j))],values[:,j])
 def test_reject_single_stage(self):
  with self.assertRaises(ValueError):smooth_composition(np.array([[100.]]))

 def test_clone_leadin_preserves_measured_knots_and_wt(self):
  from altanalyze3.components.bam.clonal_inference.clone_fish import fish_coordinates
  nodes=[0,1,3,-1];parents=np.array([0,0,1,2])
  values=np.array([[48.,80.,5.],[40.,5.,20.],[10.,10.,70.],[2.,5.,5.]])
  x,y=fish_coordinates(values,nodes,parents)
  np.testing.assert_allclose(y[:,0],[100,0,0,0])
  np.testing.assert_allclose(y.sum(0),100)
  self.assertTrue((y>=0).all())
  for j in range(values.shape[1]):np.testing.assert_allclose(y[:,np.argmin(abs(x-j))],values[:,j])
  # A descendant appears only after the displayed founder in the synthetic lead-in.
  self.assertLess(x[np.flatnonzero(y[1]>0)[0]],x[np.flatnonzero(y[2]>0)[0]])
  # The outer expansion has one inflection, not repeated synthetic plateaus.
  burden=100-y[0,x<=0]
  self.assertTrue((np.diff(burden)>0).all())
  curvature=np.diff(burden,n=2)
  signs=np.sign(curvature[np.abs(curvature)>1e-7])
  self.assertEqual(np.count_nonzero(signs[1:]!=signs[:-1]),1)

 def test_clone_renderer_keeps_background_vector(self):
  from altanalyze3.components.bam.clonal_inference.clone_fish import draw_clone_fish
  import matplotlib.pyplot as plt
  import io
  fig,ax=plt.subplots()
  draw_clone_fish(ax,np.array([[50.,75.],[45.,20.],[5.,5.]]),[0,1,-1],
                  np.array([0,0]),{0:'#edcf91',1:'#e41a1c',-1:'#bbbbbb'},['Diagnosis','Post'])
  s=io.StringIO();fig.savefig(s,format='svg');plt.close(fig)
  self.assertNotIn('<image',s.getvalue())
