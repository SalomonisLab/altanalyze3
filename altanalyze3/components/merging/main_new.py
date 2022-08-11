from genericpath import exists
from operator import truediv
import pandas as pd
import os
from functools import cache

class   JunctionAnnotation:
    def __init__(self):
        self.start_annotation = {},
        self.stop_annotation = {},
        self.junction_coordinate_BAM_dict = {}

    def each_junction_annotation(self,junction_dir, gene_model_all):
        '''
        Input args: path to all BAM derived sample txt files, gene_model_all txt file(hg38)
        Output: map of annotations for a given junction
        '''
        junction_files = os.listdir(junction_dir)
        junction_dict = self.generate_junction_dict(gene_model_all)
        junction_coordinate_BAM_listof_dicts = []
        for junction_file in junction_files:
            print("Reading junction file" + '' + junction_file)
            with open(junction_dir + junction_file) as f:
                
                for junction in f:
                    chr = junction.split('\t')[0]
                    junction_start = int(junction.split('\t')[1])
                    junction_stop = int(junction.split('\t')[2])
                    start_tar_tup = (chr, junction_start)
                    stop_tar_tup = (chr, junction_stop)
                    #print(start_tar_tup)
                    if junction_dict.get(start_tar_tup) != None:
                        
                        self.start_annotation = {'exon_region_id':junction_dict[start_tar_tup]['exon_region_id'],
                        'gene_id':junction_dict[start_tar_tup]['gene_id'],'start':junction_start}
                        #print(self.start_annotation)
                    if junction_dict.get(stop_tar_tup) != None:
                        self.stop_annotation = {'exon_region_id':junction_dict[stop_tar_tup]['exon_region_id'],
                        'gene_id':junction_dict[stop_tar_tup]['gene_id'],'stop':junction_stop}
                    #print(self.start_annotation,self.stop_annotation)
                    self.junction_coordinate_BAM_dict[(chr,junction_start,junction_stop)] = self.splice_annotations(self.start_annotation,self.stop_annotation)
        #print(self.junction_coordinate_BAM_dict)
                   
        return self.junction_coordinate_BAM_dict

   
    @cache
    def generate_junction_dict(self,gene_model_all):
        junction_dict = {}
        gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
                "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"]) 
        print("Generating junction dictionary from gene model(hg38).....")
        
        for idx,row in gene_model_df.iterrows():
            
            junction_dict[(row.chr,row.start)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id}
            junction_dict[(row.chr,row.stop)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id}
        
        print("Finished generating junction map from gene model")
        return junction_dict
    
    def splice_annotations(self,start_ann,stop_ann):
        annotations = {}
        #print("splice annotation")
        return annotations
    

gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
junction_annot = JunctionAnnotation()
junction_annot.each_junction_annotation(junction_dir=junction_dir,gene_model_all=gene_model_all)
#junction_annot.generate_junction_dict('/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt')
