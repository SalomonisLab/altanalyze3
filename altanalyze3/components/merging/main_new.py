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
        self.gene_model_dict = {}
        self.gene_model_exon_dict = {}

    def each_junction_annotation(self,junction_dir, gene_model_all):
        '''
        Input args: path to all BAM derived sample txt files, gene_model_all txt file(hg38)
        Output: map of annotations for a given junction
        '''
        junction_files = os.listdir(junction_dir)
        gene_model_dict, gene_model_exon_dict = self.generate_gene_model_dict(gene_model_all)
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
                        'junction_start':junction_start,'ref_stop':gene_model_dict[start_tar_tup]['ref_stop']}
                        
                    if gene_model_dict.get(stop_tar_tup) != None:
                        self.stop_annotation = {'exon_region_id':gene_model_dict[stop_tar_tup]['exon_region_id'],
                        'junction_stop':junction_stop,'ref_start':gene_model_dict[stop_tar_tup]['ref_start']}

                    if self.stop_annotation == None and self.Start_Annotation(self) == None:
                        self.stop_annotation = 'NA'

                    elif self.stop_annotation == None:
                        candidate_gene = gene_model_dict[start_tar_tup]['gene_id']
                        self.stop_annotation = self.annotate_splice_site(chr = self.chr,junction_coordinate = self.stop,candidate_gene = candidate_gene)
                    
                    #self.junction_coordinate_BAM_dict[(chr,junction_start,junction_stop)] = Splice_Site_Annotation(self.start_annotation,self.stop_annotation)
                   
        return self.junction_coordinate_BAM_dict
    
    def annotate_splice_site(self,chr,junction_coordinate,candidate_gene):
        annotations = {}
        ref_gene_stop = self.gene_model_exon_dict[candidate_gene][0]
        ref_gene_start = self.gene_model_exon_dict[candidate_gene][-1]
        
        inRange = self.coordinate_in_range(junction_coordinate, ref_gene_stop,ref_gene_start)
        
        ea = self.gene_model_exon_dict[candidate_gene]
        
        
        # if(inRange):
        #     annotation = start_ann['exon_region_id'] + '_' + start_ann['junction_start']   
        return annotations
   
    def Start_Annotation(self):
        if self.stop_annotation == None and self.Start_Annotation() == None:
            self.start_annotation = 'NA'
        return self.stop_annotation

    def generate_gene_model_dict(self,gene_model_all):
        '''
        This helper function genr
        '''
        ea = []
        gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
                "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"])
        #print(gene_model_df)
        logging.info("Generating junction dictionary from gene model(hg38).....")
        
        for idx,row in gene_model_df.iterrows():           
            self.gene_model_dict[(row.chr,row.start)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id, 'start':row.start}
            self.gene_model_dict[(row.chr,row.stop)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id,'stop':row.stop}
            #ea = [row.chr,row.exon_region_id,row.start,row.stop]
            if(row.gene_id in self.gene_model_exon_dict):
                self.gene_model_exon_dict[row.gene_id.split(' ')[0]].append(int(row.start))
                self.gene_model_exon_dict[row.gene_id.split(' ')[0]].append(int(row.stop))
            else:
                self.gene_model_exon_dict[row.gene_id.split(' ')[0]] = [int(row.start),int(row.stop)]
        
        print(self.gene_model_exon_dict['ENSG00000223972'][0], self.gene_model_exon_dict['ENSG00000223972'][-1])
        logging.info("Finished generating junction map from gene model")
        return self.gene_model_dict,self.gene_model_exon_dict


gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
gene_model_ENSG00000223972 = '/Users/sin9gp/altanalyze3/tests/data/gene_model_ENSG00000223972.txt'
junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
junction_annot = JunctionAnnotation()
#junction_annot.each_junction_annotation(junction_dir=junction_dir,gene_model_all=gene_model_all)
junction_annot.generate_gene_model_dict(gene_model_ENSG00000223972)
