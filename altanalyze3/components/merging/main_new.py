from genericpath import exists
from operator import truediv
import pandas as pd
import os
from functools import cache
import logging

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
        gene_model_dict = self.generate_gene_model_dict(gene_model_all)
        junction_coordinate_BAM_listof_dicts = []
        for junction_file in junction_files:
            logging.info("Reading junction file" + '' + junction_file)
            with open(junction_dir + junction_file) as f:
                
                for junction in f:
                    chr = junction.split('\t')[0]
                    junction_start = int(junction.split('\t')[1])
                    junction_stop = int(junction.split('\t')[2])
                    start_tar_tup = (chr, junction_start)
                    stop_tar_tup = (chr, junction_stop)
                    
                    if gene_model_dict.get(start_tar_tup) != None:
                        self.start_annotation = {'exon_region_id':gene_model_dict[start_tar_tup]['exon_region_id'],
                        'gene_id':gene_model_dict[start_tar_tup]['gene_id'],'junction_start':junction_start, 'ref_stop':gene_model_dict[start_tar_tup]['ref_stop']}
                        
                    if gene_model_dict.get(stop_tar_tup) != None:
                        self.stop_annotation = {'exon_region_id':gene_model_dict[stop_tar_tup]['exon_region_id'],
                        'gene_id':gene_model_dict[stop_tar_tup]['gene_id'],'junction_stop':junction_stop, 'ref_start':junction_start}
                    
                    self.junction_coordinate_BAM_dict[(chr,junction_start,junction_stop)] = Splice_Site_Annotation(self.start_annotation,self.stop_annotation)
                   
        return self.junction_coordinate_BAM_dict

   
    @cache
    def generate_gene_model_dict(self,gene_model_all):
        '''
        This helper function genr
        '''
        gene_model_dict = {}
        gene_model_exon_dict = {}
        ea = []
        gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
                "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"])
        print(gene_model_df)
        logging.info("Generating junction dictionary from gene model(hg38).....")
        
        for idx,row in gene_model_df.iterrows():           
            gene_model_dict[(row.chr,row.start)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id, 'start':row.start}
            gene_model_dict[(row.chr,row.stop)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id,'stop':row.stop}
            ea = [row.chr,row.exon_id,row.start,row.stop]
            if(row.gene_id in gene_model_exon_dict):
                gene_model_exon_dict[row.gene_id].append(ea)
            else:
                gene_model_exon_dict[row.gene_id] = [ea]
        print("hello")
        logging.info("Finished generating junction map from gene model")
        return gene_model_dict
    
   


class Splice_Site_Annotation:

    def __init__(self,chr,start,stop,startannotation,stopannotation):
        self.chr = chr
        self.start = start
        self.stop = stop
        self.startannotation = startannotation
        self.stopannotation = stopannotation
    
    def start_annotation(self):
        if self.stopannotation == None and self.stop_annotation() == None:
            self.start_annotation = 'NA'
        return self.start_annotation
    
    def stop_annotation(self):
        if self.stop_annotation == None and self.start_annotation() == None:
            self.stop_annotation = 'NA'
        elif self.stop_annotation == None:
            candidategene = self.find_gene(self.start_annotation())
            self.stop_annot = self.splice_annotations(self.chr,self.stop,candidategene)
        return self.stop_annot
    
    def find_gene(start_annotation):
        return start_annotation['gene_id']
    
    



    def splice_annotations(self,start_ann,stop_ann):
        annotations = {}
        junction_start = start_ann['junction_start']
        junction_stop = stop_ann['junction_stop']
        genemodel_stop = start_ann['ref_stop']
        genemodel_start = stop_ann['ref_start']
        inRange = self.coordinate_in_range(junction_start, junction_stop, genemodel_start,genemodel_stop)
        if(inRange):
            annotation = start_ann['exon_region_id'] + '_' + start_ann['junction_start']   
        return annotation
    
    def coordinate_in_range(self,start,stop,ref_start,ref_stop):
        if(start == ref_start and stop == ref_stop):
            return True

gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
junction_annot = JunctionAnnotation()
#junction_annot.each_junction_annotation(junction_dir=junction_dir,gene_model_all=gene_model_all)
junction_annot.generate_gene_model_dict('/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt')
