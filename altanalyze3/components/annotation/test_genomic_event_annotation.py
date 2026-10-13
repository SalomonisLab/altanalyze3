import gzip,json,tempfile,unittest
from pathlib import Path
from . import junction_isoform as J,splice_event as S
from .trans_splicing import JunctionGeometry


class AnnotationIntegrationChecks(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        header='gene\texon-id\tchromosome\tstrand\texon-region-start(s)\texon-region-stop(s)\n'
        data=header+'ENSGA\tE1.1\tchr1\t+\t100\t200\nENSGA\tE2.1\tchr1\t+\t400\t500\nENSGOVER\tE1.1\tchr1\t+\t100\t200\nENSGOVER\tE2.1\tchr1\t+\t400\t500\nENSGDIST\tE1.1\tchr1\t+\t1000\t1100\nENSGDIST\tE2.1\tchr1\t+\t1400\t1500\n'
        self.exons=self.write_gz('Hs_Ensembl_exon.txt.gz',data);self.geometry=JunctionGeometry(self.exons)
    def tearDown(self):self.temp.cleanup()
    def write_gz(self,name,data):
        p=self.root/name
        with gzip.open(p,'wt') as f:f.write(data)
        return p
    def classify(self,event,c1,c2):return S.classify_splicing_event(event,'+',['E1.1','E2.1'],['E1.1','E2.1'],c1,c2,geometry=self.geometry)
    def test_overlapping_LAT_like_event_is_not_trans(self):
        self.assertEqual(self.classify('ENSGA:E1.1-ENSGOVER:E2.1|ENSGA:E1.1-E2.1','chr1:200-400','chr1:200-400')[0],'Overlapping-Locus Annotation')
    def test_distinct_physical_loci_event_is_trans_candidate(self):
        self.assertEqual(self.classify('ENSGA:E1.1-ENSGDIST:E2.1|ENSGA:E1.1-E2.1','chr1:200-1400','chr1:200-400'),('Trans-Splicing','distal'))
    def test_missing_geometry_never_uses_gene_prefix_count(self):
        r=S.classify_splicing_event('ENSGA:E1.1-ENSGDIST:E2.1|ENSGA:E1.1-E2.1','+','E1.1|E2.1','E1.1|E2.1','chr1:200-1400','chr1:200-400')
        self.assertEqual(r[0],'Unresolved Multi-Gene Annotation')
    def test_same_gene_qualified_twice_is_not_trans(self):
        self.assertNotEqual(self.classify('ENSGA:E1.1-ENSGA:E2.1|ENSGA:E1.1-E2.1','chr1:200-400','chr1:200-400')[0],'Trans-Splicing')
    def test_gzip_logical_and_explicit_readers_and_exact_protein_ids(self):
        ta=self.write_gz('transcript_associations.txt.gz','ENSGA\t+\tE1.1|E2.1\t11815285.WF40_ND21-251_HSC_3k__source000\tsource\n')
        prot=self.write_gz('protein_summary.txt.gz','Gene ID\tTranscript ID\tProtein Length\tNMD Status\tIntron Retention\tLongest Isoform Length\nENSGA\t11815285.WF40_ND21-251_HSC_3k__source000\t57\tNot-NMD\tTrue\t742\nENSGA\tENSGA_ENST1\t100\tNot-NMD\tFalse\t742\n')
        a=J.load_transcript_associations(ta.with_suffix(''),dataType='isoform');b=J.load_transcript_associations(ta,dataType='isoform');self.assertEqual(a,b)
        p=J.load_protein_summary(prot.with_suffix(''));self.assertEqual(J._pd_lookup('ENSGA','11815285.WF40_ND21-251_HSC_3k__source000',p)[0],57);self.assertEqual(p['ENSGA','ENST1'][0],100)
    def test_compressed_statistics_preserve_every_value_and_record_geometry(self):
        stats=self.write_gz('contrast_stats.txt.gz','Feature\tebayesSign\tDiffMeans\tpFDR\nENSGA:E1.1-ENSGOVER:E2.1|ENSGA:E1.1-E2.1\t1\t0.2\t0.01\n')
        J.event_isoform_db={}
        td={'ENSGA':{'E1.1-ENSGOVER:E2.1':[('ENST1',['E1.1','ENSGOVER:E2.1'])],'E1.1-E2.1':[('ENST2',['E1.1','E2.1'])]}}
        pd={('ENSGA','ENST1'):(57,'Not-NMD',100),('ENSGA','ENST2'):(100,'Not-NMD',100)}
        cd={'ENSGA:E1.1-ENSGOVER:E2.1':'chr1:200-400','ENSGA:E1.1-E2.1':'chr1:200-400'}
        J.annotate_junction_stats_file(stats,td,{'ENSGA':'LAT'},cd,pd,self.root,geometry=self.geometry)
        import pandas as pandas
        result=pandas.read_csv(self.root/'contrast_stats-annotated.txt',sep='\t');self.assertEqual(result.loc[0,'pFDR'],0.01);self.assertEqual(result.loc[0,'DiffMeans'],0.2);self.assertEqual(result.loc[0,'Event_Type'],'Overlapping-Locus Annotation');self.assertEqual(json.loads(result.loc[0,'Trans_Splicing_Evidence'])[0]['candidate'],False)


if __name__=='__main__':unittest.main()
