"""Focused display-filter checks, using the observed CFAP45 splice coordinates."""
import unittest
from copy import deepcopy
from altanalyze3.components.visualization.isv_web.read_filters import bleeding_intervals,filter_bleeding_reads,cap_reads_per_isoform

class BleedingReadFilterTests(unittest.TestCase):
    # Relevant blocks from the saved TAGLN2/CFAP45 representative and read 10800002.
    parent={'1':[(159890480,159890622),(159893180,159893305),(159918763,159918941)]}
    retained={'1':[(159892584,159893305),(159918763,159918941)]}
    def test_actual_retained_interval_is_flagged(self):
        hits=bleeding_intervals(self.retained,self.parent)
        self.assertEqual(hits[0]['intronic_extension'],[159892584,159893179])
    def test_actual_spliced_read_is_kept(self):
        self.assertEqual(bleeding_intervals(self.parent,self.parent),[])
    def test_shorter_read_does_not_create_missing_exons(self):
        self.assertEqual(bleeding_intervals({'1':[(159893200,159893305)]},self.parent),[])
    def test_free_outer_end_extension_is_kept(self):
        self.assertEqual(bleeding_intervals({'1':[(159890000,159890622)]},self.parent),[])
    def test_same_coordinates_on_another_chromosome_are_not_hidden(self):
        self.assertEqual(bleeding_intervals({'2':self.retained['1']},self.parent),[])
    def test_wholly_intronic_exon_is_not_called_edge_bleeding(self):
        self.assertEqual(bleeding_intervals({'1':[(159892584,159892900)]},self.parent),[])
    def test_retention_in_the_assigned_model_is_allowed(self):
        self.assertEqual(bleeding_intervals(self.retained,self.retained),[])
    def test_exact_tolerance(self):
        self.assertTrue(bleeding_intervals({'1':[(159893179,159893305)]},self.parent,0))
        self.assertFalse(bleeding_intervals({'1':[(159893179,159893305)]},self.parent,1))
    def test_drawing_copy_keeps_source_counts_and_hidden_audit(self):
        seg=lambda blocks:[dict(gene='G',chrom='chr1',start=a,end=b,type='E') for a,b in blocks]
        rows=[dict(isoform_id='10800002',final_isoform_id='model',count=1,exon_segments=seg(self.retained['1']),registered_exon_segments=seg(self.parent['1'])),dict(isoform_id='8684604',final_isoform_id='model',count=1,exon_segments=seg(self.parent['1']),registered_exon_segments=seg(self.parent['1']))]
        source=dict(gene='G',panels=[dict(condition='pooled',molecules=rows,n_molecules=2,structure_reads={'model':2})]);before=deepcopy(source)
        self.assertIs(filter_bleeding_reads(source,enabled=False),source)
        result=filter_bleeding_reads(source);self.assertEqual(source,before)
        self.assertEqual(result['panels'][0]['structure_reads'],{'model':2})
        self.assertEqual(result['hidden_bleeding_reads'][0]['isoform_id'],'10800002')
        self.assertEqual(result['panels'][0]['molecules'][0]['isoform_id'],'8684604')
        capped=cap_reads_per_isoform(result,200);self.assertEqual(capped['panels'][0]['n_molecules'],1)
    def test_missing_model_is_not_silently_excluded(self):
        with self.assertRaisesRegex(ValueError,'registered model'):
            filter_bleeding_reads(dict(gene='G',panels=[dict(molecules=[dict(isoform_id='r',final_isoform_id='unknown')])]))

if __name__=='__main__':unittest.main()
