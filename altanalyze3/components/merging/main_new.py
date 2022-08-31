from genericpath import exists
from json.encoder import INFINITY
from operator import truediv
from tracemalloc import start, stop
from unittest import skip
import pandas as pd
import os
import re
from functools import cache
import logging

class   JunctionAnnotation:
    def __init__(self):
        
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
                        if row.chr == 'X':
                            chr = 23
                        elif row.chr == 'Y':
                            chr = 24
                        
                        start_tar_tup = (int(chr), row.start)
                        stop_tar_tup = (int(chr), row.stop + 1)
                        
                        
                        if self.gene_model_dict.get(start_tar_tup) != None:
                            start_annotation[start_tar_tup] = { 'exon_region_id':self.gene_model_dict[start_tar_tup]['exon_region_id'],
                            'junction_start':row.start, 'candidate_gene':self.gene_model_dict[start_tar_tup]['gene_id']}
                         
                        
                        if self.gene_model_dict.get(stop_tar_tup) != None:
                            stop_annotation[stop_tar_tup] = { 'exon_region_id':self.gene_model_dict[stop_tar_tup]['exon_region_id'],
                            'junction_stop':row.stop + 1, 'candidate_gene':self.gene_model_dict[stop_tar_tup]['gene_id']}
                        

                        elif self.gene_model_dict.get(start_tar_tup) == None and self.gene_model_dict.get(stop_tar_tup) == None:
                            continue

                        
                        if start_annotation and stop_annotation:
                            gene_id = ''
                            if start_annotation[(start_tar_tup)]['candidate_gene'] == stop_annotation[(stop_tar_tup)]['candidate_gene']:
                                gene_id = self.gene_model_dict[(start_tar_tup)]['gene_id']

                                if(self.gene_model_dict[start_tar_tup]['strand'] == '-'):
                                    annotation = gene_id + ':' + stop_annotation[(stop_tar_tup)]['exon_region_id'] + '-' + start_annotation[(start_tar_tup)]['exon_region_id']
                                else:
                                    annotation = gene_id + ':' + start_annotation[(start_tar_tup)]['exon_region_id'] + '-' + stop_annotation[(stop_tar_tup)]['exon_region_id']
                            else:
                                start_gene_id = self.gene_model_dict[(start_tar_tup)]['gene_id']
                                stop_gene_id = self.gene_model_dict[(stop_tar_tup)]['gene_id']
                                annotation = start_gene_id + ':' + self.gene_model_dict[(start_tar_tup)]['exon_region_id'] + '-' + stop_gene_id + ":" + self.gene_model_dict[(stop_tar_tup)]['exon_region_id']


                            annotations.append(annotation)
                            annotation_key = 'chr' + str(chr) + ':' + str(row.start) + '-' + str(row.stop + 1)
                            annotation_keys.append(annotation_key)     

                        
                        if start_annotation and not stop_annotation:
                            strand = self.gene_model_dict[start_tar_tup]['strand']
                            splice_site = self.annotate_splice_site(chr = row.chr, junction_coordinate = row.stop + 1,
                            candidate_gene = start_annotation[start_tar_tup]['candidate_gene'], exon_region_id = self.gene_model_dict[start_tar_tup]['exon_region_id'], strand = strand)                         
                            
                           
                            annotations.append(start_annotation[start_tar_tup]['candidate_gene'] + ':' + start_annotation[start_tar_tup]['exon_region_id'] + '-'+  splice_site)
                            annotation_key = 'chr' + str(chr) + ':' + str(row.start) + '-' + str(row.stop + 1)
                            annotation_keys.append(annotation_key)

                        if  stop_annotation and not start_annotation:
                            splice_site = self.annotate_splice_site(chr = row.chr, junction_coordinate = row.start,
                            candidate_gene = stop_annotation[stop_tar_tup]['candidate_gene'], exon_region_id = self.gene_model_dict[stop_tar_tup]['exon_region_id'], strand = self.gene_model_dict[stop_tar_tup]['strand'])

                            annotations.append(stop_annotation[stop_tar_tup]['candidate_gene'] + ':' + stop_annotation[stop_tar_tup]['exon_region_id'] + '-'+  splice_site)
                            annotation_key = 'chr' + str(chr) + ':' + str(row.start) + '-' + str(row.stop + 1)
                            annotation_keys.append(annotation_key)
                            
                        elif not stop_annotation and not start_annotation:
                            continue

        data = { 'annotation_key':annotation_keys,'annotations':annotations}
        df1 = pd.DataFrame(data)
        df1.to_csv('annotations_all.txt',sep='\t')
    
    def annotate_splice_site(self,chr,junction_coordinate,candidate_gene, exon_region_id, strand):
        annotation = ''
        buffer = 0
        ref_gene_start = self.gene_model_exon_dict[candidate_gene][0][0]
        ref_gene_stop = self.gene_model_exon_dict[candidate_gene][-1][-1]
        status = self.coordinate_in_range(junction_coordinate, ref_gene_start,ref_gene_stop, buffer = 0, strand=strand)
        
        if status == True:
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
                status = self.coordinate_in_range(junction_coordinate, exon_start, exon_stop, buffer, strand)
                if(status):
                    annotation = str(exon_id) + "_" + str(junction_coordinate)
                    
                    if intron_status == False:
                        return annotation
                    else:
                        candidate_found = True
                else:
                    if candidate_found:
                        return annotation
        else:
            status = self.coordinate_in_range(junction_coordinate, ref_gene_start,ref_gene_stop, buffer = 20000, strand = strand)
            region_numbers = []
            for ea in self.gene_model_exon_dict[candidate_gene]:
                region_numbers.append(int(str.split(ea[2][1:],'.')[0]))
            
            if status == True:
                if(strand == '+'): #applicable to 3'UTR
                    annotation = 'U'+str(region_numbers[-1])+'.1_'+str(junction_coordinate)
                else: #5' UTR
                    annotation = 'U0.1' + '_' + str(junction_coordinate)
            else:
                #iterate over all the genes and fine which junction in gene model is nearest to this coordinate
                # should be in the same gene               
                for (chrom, junction), value in self.gene_model_dict.items():
                    if chr == chrom:
                        if junction_coordinate >=self.gene_model_dict[(chrom,junction)]['start'] and junction_coordinate <=self.gene_model_dict[(chrom,junction)]['stop']:
                            print("this worked")
                            print(self.gene_model_dict[(chrom,junction)]['gene_id'] + ':' + self.gene_model_dict[(chrom,junction)]['exon_region_id'] + '_' + str(junction_coordinate))
                            return self.gene_model_dict[(chrom,junction)]['gene_id'] + ':' + self.gene_model_dict[(chrom,junction)]['exon_region_id'] + '_' + str(junction_coordinate)
                        else:
                            
                            continue
                         
                    else:
                        continue  
                            
                print("annotated the weird junction", annotation)

        return annotation

                
    
    def coordinate_in_range(self,junction_coordinate, ref_gene_start,ref_gene_stop, buffer, strand):
        if ref_gene_start <= ref_gene_stop or strand == '+':
            start = int(ref_gene_start) - buffer
            stop = int(ref_gene_stop) + buffer
            if junction_coordinate in range(start, stop):
                return True
            else:
                return False
        elif ref_gene_start > ref_gene_stop or strand == '-':
            start = int(ref_gene_start) + buffer
            stop = int(ref_gene_stop) - buffer
            if junction_coordinate in range(stop, start):
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
            chr = row.chr
            if row.chr == 'X':
                chr = 23
            elif row.chr == 'Y':
                chr = 24
            elif re.search("^[a-zA-Z]", row.chr) is not None:
                continue
            self.gene_model_dict[(int(chr),row.start)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id, 'start':row.start, 'stop':row.stop,'strand':row.strand}
            self.gene_model_dict[(int(chr),row.stop)] = {'gene_id':row.gene_id, 'exon_region_id':row.exon_region_id,'stop':row.stop, 'start':row.start, 'strand':row.strand}
            ea = [row.start,chr,row.exon_region_id,row.strand,row.stop]
            #ea store as object?
            if(row.gene_id in self.gene_model_exon_dict):
                self.gene_model_exon_dict[row.gene_id].append(ea)
            else:
                self.gene_model_exon_dict[row.gene_id] = [ea]
        
        #print(self.gene_model_exon_dict['ENSG00000117335'][0], self.gene_model_exon_dict['ENSG00000117335'][-1])
        print("Finished generating junction map from gene model")
        # print("gene dict for (1, 17619445)", self.gene_model_dict[(1, 17619445)])
        # print("gene dict for (1, 17621863)",self.gene_model_dict[(1, 17621863)])
        return self.gene_model_dict,self.gene_model_exon_dict


gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
gene_model_ENSG00000223972 = '/Users/sin9gp/altanalyze3/tests/data/gene_model_ENSG00000223972.txt'
gene_model_1 = '/Users/sin9gp/altanalyze3/tests/data/gene_model_chr17_U.txt'
junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
subset_dir = '/Users/sin9gp/altanalyze3/tests/data/subset/newsubset/new_newsubset/'
junction_annot = JunctionAnnotation()
junction_annot.each_junction_annotation(junction_dir=subset_dir,gene_model_all=gene_model_all)
#junction_annot.generate_gene_model_dict(gene_model_all)
