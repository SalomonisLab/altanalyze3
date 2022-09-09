from xml.dom.minidom import TypeInfo
import pandas as pd
import os
import timeit
import logging
import multiprocessing as mp

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
                logging.info("Reading junction file " + ' ' + junction_file)
                with open(junction_dir + junction_file) as f:
                    
                    for line in f:
                        (chr, start, stop, annotation, splice_count) = line.split('\t')
                    
                    
                        start_annotation = {}
                        stop_annotation = {}
                        stop = int(stop)
                        start_tar_tup = (str(chr), int(start))
                        stop_tar_tup = (str(chr), int(stop) + 1)
                        
                        
                        annotation_key = ''.join(['chr', str(chr), ':', str(start),'-',str(stop + 1)])
                        if self.gene_model_dict.get(start_tar_tup) != None:
                            #print("i found start")
                            start_annotation[start_tar_tup] = { 'exon_region_id':self.gene_model_dict[start_tar_tup]['exon_region_id'],
                            'junction_start':start, 'candidate_gene':self.gene_model_dict[start_tar_tup]['gene_id']}
                         
                        
                        if self.gene_model_dict.get(stop_tar_tup) != None:
                            #print("i found stop")
                            stop_annotation[stop_tar_tup] = { 'exon_region_id':self.gene_model_dict[stop_tar_tup]['exon_region_id'],
                            'junction_stop':stop + 1, 'candidate_gene':self.gene_model_dict[stop_tar_tup]['gene_id']}
                        

                        elif self.gene_model_dict.get(start_tar_tup) == None and self.gene_model_dict.get(stop_tar_tup) == None:
                            
                            continue

                        
                        if start_annotation and stop_annotation:
                            gene_id = ''
                            if start_annotation[(start_tar_tup)]['candidate_gene'] == stop_annotation[(stop_tar_tup)]['candidate_gene']:
                                gene_id = self.gene_model_dict[(start_tar_tup)]['gene_id']

                                if(self.gene_model_dict[start_tar_tup]['strand'] == '-'):
                                    annotation = ''.join([gene_id, ':', stop_annotation[(stop_tar_tup)]['exon_region_id'], '-', start_annotation[(start_tar_tup)]['exon_region_id']])
                                else:
                                    annotation = ''.join([gene_id,':',start_annotation[(start_tar_tup)]['exon_region_id'],'-', stop_annotation[(stop_tar_tup)]['exon_region_id']])
                            else:
                                start_gene_id = self.gene_model_dict[(start_tar_tup)]['gene_id']
                                stop_gene_id = self.gene_model_dict[(stop_tar_tup)]['gene_id']
                                annotation = ''.join([start_gene_id, ':', self.gene_model_dict[(start_tar_tup)]['exon_region_id'], '-', stop_gene_id, ":",self.gene_model_dict[(stop_tar_tup)]['exon_region_id']])
                            annotations.append(annotation)
                            annotation_keys.append(annotation_key)
                            print("both start and stop were found")
                            print(annotation_key)     
                            continue

                        
                        if start_annotation and not stop_annotation:
                            annotated_splice_site = self.annotate_splice_site(chr = chr, junction_coordinate = int(stop) + 1,
                            candidate_gene = start_annotation[start_tar_tup]['candidate_gene'], exon_region_id = self.gene_model_dict[start_tar_tup]['exon_region_id'], strand = self.gene_model_dict[start_tar_tup]['strand'])                         
                            annotations.append(''.join([start_annotation[start_tar_tup]['candidate_gene'],':', start_annotation[start_tar_tup]['exon_region_id'],'-',  annotated_splice_site]))
                            annotation_keys.append(annotation_key)
                            print("start was found not stop")
                            print(annotation_key)  
                            continue

                        if  stop_annotation and not start_annotation:
                            annotated_splice_site = self.annotate_splice_site(chr = chr, junction_coordinate = int(start),
                            candidate_gene = stop_annotation[stop_tar_tup]['candidate_gene'], exon_region_id = self.gene_model_dict[stop_tar_tup]['exon_region_id'], strand = self.gene_model_dict[stop_tar_tup]['strand'])
                            annotations.append(''.join([stop_annotation[stop_tar_tup]['candidate_gene'],':', stop_annotation[stop_tar_tup]['exon_region_id'], '-',  annotated_splice_site]))
                            annotation_keys.append(annotation_key)
                            print(annotation_key)  
                            print("stop was found but not start")
                            continue
                            
                        elif not stop_annotation and not start_annotation:
                            continue

        self.write_to_file(annotation_keys, annotations)
    
    def write_to_file(self, annotation_keys, annotations):
        data = { 'annotation_key':annotation_keys,'annotations':annotations}
        df1 = pd.DataFrame(data)
        df1.to_csv('annotations_all.txt',sep='\t')


    def annotate_splice_site(self,chr,junction_coordinate,candidate_gene, exon_region_id, strand):
        annotation = ''
        buffer = 0
        ref_gene_start = self.gene_model_exon_dict[candidate_gene][0]['start']
        ref_gene_stop = self.gene_model_exon_dict[candidate_gene][-1]['stop']
        status = self.coordinate_in_range(junction_coordinate, ref_gene_start,ref_gene_stop, buffer = 0, strand=strand)
        
        if status == True:
            # preset when initially looking (always false to start with)
            candidate_found = False 
            
            for ea in self.gene_model_exon_dict[candidate_gene]:
                exon_start = ea['start']
                exon_stop = ea['stop']
                exon_id = ea['exon_region_id']
                if 'I' in exon_id:
                    intron_status = True
                    buffer = 0
                else:
                    intron_status = False
                    buffer = 50
                status = self.coordinate_in_range(junction_coordinate, exon_start, exon_stop, buffer, strand)
                if(status):
                    annotation = ''.join([str(exon_id), "_", str(junction_coordinate)])
                    
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
                region_numbers.append(int(str.split(ea['exon_region_id'][1:],'.')[0]))
            
            if status == True:
                # applicable to 3'UTR
                if(strand == '+'): 
                    annotation =  ''.join(['U', str(region_numbers[-1]),'.1_', str(junction_coordinate)])
                # applicable to 5'UTR
                else: 
                    annotation =  ''.join(['U0.1', '_', str(junction_coordinate)])
            else:
                # iterate over all the genes on same chromosome and find which junction in gene model 
                # is nearest to this coordinate         
                for (chrom, junction), value in self.gene_model_dict.items():
                    if chr == chrom:
                        if int(junction_coordinate) >= int(self.gene_model_dict[(chrom,junction)]['start']) and int(junction_coordinate) <= int(self.gene_model_dict[(chrom,junction)]['stop']):
                            annotation =  ''.join([self.gene_model_dict[(chrom,junction)]['gene_id'], ':', self.gene_model_dict[(chrom,junction)]['exon_region_id'], '_', str(junction_coordinate)])
                        else:
                            continue
                    else:
                        continue  

        return annotation
    
    def coordinate_in_range(self,junction_coordinate, ref_gene_start,ref_gene_stop, buffer, strand):
        #needs refactoring
        if ref_gene_start <= ref_gene_stop or strand == '+':
            start = int(ref_gene_start) - buffer
            stop = int(ref_gene_stop) + buffer
            
        else:
            start = int(ref_gene_start) + buffer
            stop = int(ref_gene_stop) - buffer

        if int(start) <= int(junction_coordinate) <= int(stop):
            return True
        else:
            return False

        


    def generate_gene_model_dict(self,gene_model_all):
        '''
        This helper function takes gene reference file (generated from Frank's gene_model)
        and creates a map of each junction coordinate which will be later used to find splice 
        site annotations for junctions in sample files. 

        Args: gene_model.txt file
        Output: 
        1. gene_model dictionary where key is (chr, coordinate) and value has gene_id, exon_region_id, strand
        2. gene_model_exon dictionary where key is (gene_id)
        '''
        
        logging.info("Generating reference junction dictionary from gene model(hg38).....")
        starttime = timeit.default_timer()
        print("The start time is :",starttime)
        
        with open(gene_model_all) as f:
            for line in f:
                (gene_id,chr,strand,exon_region_id, start,stop,annotation) = line.split('\t')
                
                self.gene_model_dict[(chr,int(start))] =  {'gene_id':gene_id, 'exon_region_id':exon_region_id, 'start':start, 'stop':stop,'strand':strand}
                self.gene_model_dict[(chr,int(stop))] =  {'gene_id':gene_id, 'exon_region_id':exon_region_id, 'start':start, 'stop':stop,'strand':strand}
                ea = {'exon_region_id':exon_region_id, 'start':start, 'stop':stop,'strand':strand, 'chr':chr}
            #ea store as map - more efficient than list.
                if(gene_id in self.gene_model_exon_dict):
                    self.gene_model_exon_dict[gene_id].append(ea)
                else:
                    self.gene_model_exon_dict[gene_id] = [ea]
        print(self.gene_model_dict[('X', 103707210)])
        logging.info("Gene Model dictionary generated")
        print("The time difference is :", timeit.default_timer() - starttime)

if __name__ == '__main__':
    print("Number of processors: ", mp.cpu_count())
    pool = mp.Pool(mp.cpu_count())
    logging.basicConfig(level=logging.INFO)
    starttime = timeit.default_timer()
    print("The start time is :",starttime)
    gene_model_all = '/Users/sin9gp/altanalyze3/tests/data/gene_model_all.txt'
    junction_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'
    # pysam.sort("-o", str(args.output.with_suffix(".bam")), str(tmp_bam))
    subset_dir = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/subset/'
    junction_annot = JunctionAnnotation()
    #junction_annot.generate_gene_model_dict(gene_model_all=gene_model_all)
    junction_annot.each_junction_annotation(junction_dir=subset_dir,gene_model_all=gene_model_all)

    print("The time difference is :", timeit.default_timer() - starttime)

