from operator import truediv
import pandas as pd
import os

#path for gene_model and junction files

junction_dict = {}
exon_dict = {}
gene_dict = {}

def sampleBEDFileAnnotation(gene_model_all, junction_dir, intron_dir):
    exon_dict, gene_dict,junction_dict = importJunctionInfo(gene_model_all)
    importJunctionFromBED(junction_dir, junction_dict)


def importJunctionInfo(gene_model_all):
    '''
    This method will take gene_model as input and create three dictionaries
    1. junction_dict
    2. exon_dict
    3. gene_dict
    '''
    

    gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
            "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"])    
    #junction_dict = gene_model_df.set_index('chr','start').to_dict()['gene_id']
    for idx,row in gene_model_df.iterrows():
        junction_dict[(row.chr,row.start)] = row.gene_id + ':' + row.exon_region_id
        junction_dict[(row.chr,row.stop)] = row.gene_id + ':' + row.exon_region_id
        exon_annotation = get_exon_annotation(row.chr,row.start,row.stop,row.exon_region_id)
        if row.gene_id in exon_dict.keys():
            exon_dict[row.gene_id].append(exon_annotation)
            gene_dict[row.gene_id].append(row.start)
            gene_dict[row.gene_id].append(row.stop)
        else:
            exon_dict[row.gene_id] = [exon_annotation]
            gene_dict[row.gene_id] = [row.start,row.stop]
    # for gene in gene_dict:
    #     gene_dict[gene].sort()
    #print(junction_dict)
    return exon_dict, gene_dict, junction_dict

   
def get_exon_annotation(chr, start,stop,exonid):
    exon_start = chr + ':' + str(start) 
    exon_stop = chr + ':' + str(stop)
    exon_id =  exonid + '.' + str(chr)
    return {'exon_start':exon_start, 'exon_stop': exon_stop,'exon_id':exon_id}


def importJunctionFromBED(junction_dir, junction_dict):
    '''
    This function imports junction data from BAM derived exon junction text files
    '''
    junction_coordinate_BAM_dict = {}

    for junctionfilename in os.listdir(junction_dir):
        junctionfile = os.path.join(junction_dir, junctionfilename)
        junction_coordinate_BAM_dict = initialBAMJunctionImport(junctionfile)

    for (chr,start, stop) in junction_coordinate_BAM_dict.keys():
        unprefixedchr = chr.removeprefix('chr')
        startAnnotation = None
        stopAnnotation = None
        #print(unprefixedchr,start)
        #print(junction_dict)
        tar_start_tup = (unprefixedchr,start)
        tar_stop_tup = (unprefixedchr,stop)
        if junction_dict.get(tar_start_tup) != None:
            print("exists")
            startAnnotation = junction_dict[chr,start]
        if junction_dict.get(tar_stop_tup) != None:
            stopAnnotation = junction_dict[chr,stop]
        junctionAnnotation = JunctionAnnotation(startAnnotation,stopAnnotation)
        junction_coordinate_BAM_dict[chr,start,stop] = junctionAnnotation
        print(junctionAnnotation)
    return junction_coordinate_BAM_dict


def initialBAMJunctionImport(junctionfile):
    interimJunctionCoordinates={}
    junction_df = pd.read_csv(junctionfile,sep='\t',header=None, names=[
            "chr", "start", "stop","annotation", "counts","strand"])
    for idx,row in junction_df.iterrows():
        # check if it's a valid chr
        if(len(row.chr) <= 5):
            interimJunctionCoordinates[(row.chr,row.start,row.stop)] = []
    return interimJunctionCoordinates
    
class JunctionAnnotation:
    def __init__(self,chr,start,stop,start_annotation,stop_annotation):
        self.chr = chr,
        self.start = start,
        self.stop = stop,
        self.start_annotation = start_annotation
        self.stop_annotation = stop_annotation
    
    def stopAnnotation(self):
        if self.stop_annotation == None and self.startAnnotation() == None:
            self.stop_annotation = 'NA'
        elif self.stop_annotation == None:
            candidateGene = self.findGene(self.startAnnotation())
            self.stop_annotation = self.annotateSpliceSite(self.chr,self.stop,candidateGene)
        return self.stop_annotation
    
    def startAnnotation(self):
        if self.stop_annotation == None and self.start_annotation == None:
            self.start_annotation = 'NA'
        print('start annotation', self.start_annotation)
        return self.start_annotation

    def findGene(self,exonAnnotation):
        return str.split(exonAnnotation,':')[0]
    
    def annotateSpliceSite(self,chr,coordinate,candidateGene):
        gene_start = gene_dict[candidateGene][0]
        gene_stop = gene_dict[candidateGene][-1]
        print("we are gene start and stop inside annotat splice site", gene_start, gene_stop)

        inRange = self.coordinateInRange(coordinate,gene_start,gene_stop, buffer = 0)

        if inRange == False:
            self.coordinateInRange(coordinate,gene_start,gene_stop,buffer = 2000)

        if inRange: ### Yes the candidate gene is the correct gene
            candidateFound = False ### preset when initially looking (always false to start)
            for ea in exon_dict[candidateGene]:
                exon_start = ea['exon_start']
                exon_stop = ea['exon_stop']
                exon_id = ea['exon_id']
                if exon_id.startswith('I'):
                    intron_status = True
                    buffer = 0
                else:
                    intron_status = False
                    buffer = 50
                inRange = self.coordinateInRange(coordinate,exon_start,exon_stop, buffer)
                if inRange:
                    annotation = exon_id + '_' + coordinate
                    if intron_status == False:
                        return annotation
                        break
                    else:
                        ### This means that the junction is within intron 
                        # but we need to check if it is in the next exon with 50nt spacer
                        candidateFound = True
                        continue
                else:
                    if candidateFound:
                        ### annotated as a novel intron splice in the last loop
                        return annotation 
                    pass

    def coordinateInRange(coordinate, start,stop, buffer):
        if coordinate >= start - buffer and coordinate <= stop + buffer:
            return True
        return False



if __name__ == '__main__':
    genemodelpath = "/Users/sin9gp/altanalyze3/tests/data/" 
    file = 'gene_model_ENSG00000223972.txt'
    gene_model_file = f"{genemodelpath}/{file}"
    junction_dir = "/Users/sin9gp/altanalyze3/tests/data/junction_dir/"
    intron_dir = "/Users/sin9gp/altanalyze3/tests/data/intron_dir/"
    sampleBEDFileAnnotation(gene_model_file, junction_dir, intron_dir)
    # importJunctionInfo(gene_model_all=gene_model_file)
    # importJunctionFromBED(junctiondir)