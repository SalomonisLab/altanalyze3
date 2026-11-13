import copy
import unittest
import numpy as np
import pandas as pd
from altanalyze3.components.bam.circRNA_native_statistics import native_differential, NORMALIZATION, DECISIONS
from altanalyze3.components.bam.circRNA_cohort_pipeline import barcode_state, source_group, build_contrast_roster
from altanalyze3.components.oncosplice.metadataAnalysis import limmaCompute


def fixture():
    units=['h1','h2','d1','d2'];features=['signal','null','zero']
    counts=pd.DataFrame([[4,6,40,60],[5,5,5,5],[0,0,0,0]],index=features,columns=units)
    depth=pd.Series([1000,1000,1000,1000],index=units)
    groups=pd.DataFrame({'grp':['control','control','case','case'],'biological_donor':units},index=units)
    decision={'pending_decisions':[],'normalization':NORMALIZATION,'correction_family':'all_testable_features_per_state_per_contrast','authorizations':{k:{'status':'user_approved','user_verbatim':'SYNTHETIC TEST FIXTURE ONLY: this is not real analysis authorization.'} for k in DECISIONS}}
    return counts,depth,groups,features,units,'case','control',decision


class Tests(unittest.TestCase):
    def test_native_parity_and_full_feature_retention(self):
        args=fixture();r,x=native_differential(*args)
        native=limmaCompute(np.log2(args[0].divide(args[1],axis=1)*1e6+1),args[2][['grp']],min_group_size=2,shrink_factor=.2)
        np.testing.assert_allclose(r['ebayesPval'],native['mwuPval']);np.testing.assert_allclose(r['ebayesAdjPval'],native['mwuAdjPval'])
        self.assertEqual(list(r.index),args[3]);self.assertGreater(r.loc['signal','case_minus_control_mean_log2_CPM1'],0)
    def test_unresolved_donor_identity_does_not_run_test(self):
        a=fixture();r,x=native_differential(*a,inference_unavailable_reason='UNAVAILABLE_UNRESOLVED_DONOR_IDENTITY')
        self.assertTrue(r['ebayesPval'].isna().all());self.assertTrue(r['ebayesAdjPval'].isna().all());self.assertEqual(len(r),3);self.assertTrue((r['analysis_status']=='UNAVAILABLE_UNRESOLVED_DONOR_IDENTITY').all())
    def test_pending_rejected(self):
        a=list(fixture());a[-1]['pending_decisions']=['normalization']
        with self.assertRaises(ValueError):native_differential(*a)
    def test_identity_swap_rejected(self):
        a=list(fixture());a[3]=['null','signal','zero']
        with self.assertRaises(ValueError):native_differential(*a)
    def test_duplicate_donor_rejected(self):
        a=list(fixture());a[2].loc['d2','biological_donor']='d1'
        with self.assertRaises(ValueError):native_differential(*a)
    def test_positive_without_depth_rejected(self):
        a=list(fixture());a[1].loc['h1']=0
        with self.assertRaises(ValueError):native_differential(*a)
    def test_missing_depth_retained_not_zero(self):
        a=list(fixture());a[0].loc[:,'h1']=0;a[1].loc['h1']=np.nan
        r,x=native_differential(*a);self.assertTrue(x['h1'].isna().all());self.assertTrue(r['ebayesPval'].isna().all());self.assertEqual(len(r),3)
    def test_allzero_and_flags_no_loss(self):
        a=fixture();r,x=native_differential(*a);self.assertEqual(x.loc['zero'].sum(),0);self.assertIn('zero',r.index)
    def test_barcode_dual_conflict_unavailable(self):
        self.assertEqual(barcode_state('AAA-1',{'AAA-1':'HSC','TTT-1':'B'}),'__AMBIGUOUS_BARCODE__')
        self.assertEqual(barcode_state('AAA-1',{'TTT-1':'HSC'}),'HSC')
    def test_models_not_pooled_as_independent_human(self):
        self.assertIsNone(source_group({'group':'MDS-pre-SRSF2-iPSC'}));self.assertIsNone(source_group({'group':'MDS-pre-SRSF2-MISTRG'}));self.assertEqual(source_group({'group':'MDS-post-SRSF2'}),'MDS-SRSF2')
    def test_no_mutation_propagation(self):
        self.assertIsNone(source_group({'group':'AML'}));self.assertEqual(source_group({'group':'AML-SRSF2'}),'AML-SRSF2')
if __name__=='__main__':unittest.main()
