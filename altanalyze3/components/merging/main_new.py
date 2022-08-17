from genericpath import exists
from operator import truediv
from tracemalloc import stop
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
        self.generate_gene_model_dict(gene_model_all)
        junction_coordinate_BAM_listof_dicts = []
        print(junction_files)
        annotation_keys = []
        annotations = []
        for junction_file in junction_files:
            if not junction_file.startswith('.'):
                print("Reading junction file " + ' ' + junction_file)
                with open(junction_dir + junction_file) as f:
                    junction_df = pd.read_csv(f,sep='\t',header=None, 
                    names=["chr", "start", "stop", "annotation", "splice_count"])
                    
                    for idx,row in junction_df.iterrows():
                        
                        start_annotation = {}
                        stop_annotation = {}
                        chr = row.chr
                        junction_start = row.start
                        junction_stop = row.stop + 1
                        start_tar_tup = (chr, junction_start)
                        stop_tar_tup = (chr, junction_stop)
                        
                        if self.gene_model_dict.get(start_tar_tup) != None:
                            start_annotation[start_tar_tup] = { 'exon_region_id':self.gene_model_dict[start_tar_tup]['exon_region_id'],
                            'junction_start':junction_start, 'candidate_gene':self.gene_model_dict[start_tar_tup]['gene_id']}
                         
                        elif self.gene_model_dict.get(stop_tar_tup) != None:
                            stop_annotation[stop_tar_tup] = { 'exon_region_id':self.gene_model_dict[stop_tar_tup]['exon_region_id'],
                            'junction_stop':junction_stop, 'candidate_gene':self.gene_model_dict[stop_tar_tup]['gene_id']}
                            
                        else:
                            print("junction cannot be associated with a gene - do not annotate")
                        

                        if not bool(stop_annotation) and bool(start_annotation):
                            print("stop annotation is empty and start annotation exists")
                            print(start_annotation)
                            print(stop_annotation)
                            print(bool(start_annotation))
                            print(bool(stop_annotation))
                           
                            annotation = self.annotate_splice_site(chr = row.chr, junction_coordinate = junction_start, candidate_gene = start_annotation[start_tar_tup]['candidate_gene'])
                            annotation_key = 'chr' + str(chr) + ':' + str(junction_stop) + '-' + str(junction_start)
                            annotation_keys.append(annotation_key)
                            annotations.append(annotation)

                        elif not bool(start_annotation) and bool(stop_annotation):
                            annotation = self.annotate_splice_site(chr = row.chr, junction_coordinate = junction_stop, candidate_gene = stop_annotation[stop_tar_tup]['candidate_gene'])
                            annotation_key = 'chr' + str(chr) + ':' + str(junction_stop) + '-' + str(junction_start)
                            annotation_keys.append(annotation_key)
                            annotations.append(annotation)
                        
                        else:
                            print(start_annotation)
                            print(stop_annotation)
                        
        #print(annotation_keys)
        data = { 'annotation_key':annotation_keys,'annotations':annotations}
        df1 = pd.DataFrame(data)
        df1.to_csv('annotations.txt',sep='\t')

        return self.junction_coordinate_BAM_dict
    
    def annotate_splice_site(self,chr,junction_coordinate,candidate_gene):
       
        ref_gene_start = self.gene_model_exon_dict[candidate_gene][0][0]
        ref_gene_stop = self.gene_model_exon_dict[candidate_gene][-1][-1]
        status = self.coordinate_in_range(junction_coordinate, ref_gene_stop,ref_gene_start, buffer = 0)
        
        #? - this might run indefinitely until status is True ?? - should we just return the annotation?
        if status == False:
            status = self.coordinate_in_range(junction_coordinate, ref_gene_start,ref_gene_stop, buffer = 2000)
            if status == True:
                #TO_DO if strand is - or + then
                if(ref_gene_start > ref_gene_stop): #if upstream
                    annotation = 'U1' + '.1_' + str(junction_coordinate)
                    
                else:
                    annotation = 'U100' + '.1_' + str(junction_coordinate)
        else: #yes candidate gene found
            candidate_found = False #preset when initially looking (always false to start with)
            
            for ea in self.gene_model_exon_dict[candidate_gene]:
                exon_start = ea[0]
                exon_stop = ea[-1]
                exon_id = ea[2]

                if 'I' in exon_id:
                    intron_status = True
                    buffer = 0
                else:
                    intron_status = False
                    buffer = 50
                status = self.coordinate_in_range(junction_coordinate, exon_start, exon_stop, buffer)
                if(status):
                    annotation = candidate_gene + ':'+ str(exon_id) + '_' + str(junction_coordinate)
                    
                    if intron_status == False:
                        return annotation
                    else:
                        candidate_found = True
                else:
                    if candidate_found:
                        return annotation
        
        print(annotation)
        return annotation
    
    def coordinate_in_range(self,junction_coordinate, ref_gene_start,ref_gene_stop, buffer):
        #this will take care if the strand is negative as well.
        if ref_gene_stop < ref_gene_start:
            ref_gene_start, ref_gene_stop = ref_gene_stop, ref_gene_start
        if junction_coordinate in range(ref_gene_start - buffer ,ref_gene_stop + buffer):
            return True
        else:
            return False

    def generate_gene_model_dict(self,gene_model_all):
        '''
        This helper function genr
        '''
        ea = []
        gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
                "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"])
        #print(gene_model_df)
        print("Generating reference junction dictionary from gene model(hg38).....")
        
        for idx,row in gene_model_df.iterrows():           
            self.gene_model_dict[(row.chr,row.start)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id, 'start':row.start}
            self.gene_model_dict[(row.chr,row.stop)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id,'stop':row.stop}
            ea = [row.start,row.chr,row.exon_region_id,row.stop]
            #ea store as object?
            if(row.gene_id in self.gene_model_exon_dict):
                self.gene_model_exon_dict[row.gene_id].append(ea)
            else:
                self.gene_model_exon_dict[row.gene_id] = [ea]
        
        #print(self.gene_model_exon_dict['ENSG00000117335'][0], self.gene_model_exon_dict['ENSG00000117335'][-1])
        print("Finished generating junction map from gene model")
        return self.gene_model_dict,self.gene_model_exon_dict


gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
gene_model_ENSG00000223972 = '/Users/sin9gp/altanalyze3/tests/data/gene_model_ENSG00000223972.txt'
junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
junction_annot = JunctionAnnotation()
junction_annot.each_junction_annotation(junction_dir=junction_dir,gene_model_all=gene_model_all)
#junction_annot.generate_gene_model_dict(gene_model_all)
