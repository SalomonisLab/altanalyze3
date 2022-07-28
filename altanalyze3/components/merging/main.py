import pandas as pd
import os

#path for gene_model and junction files


def sampleBEDFileAnnotation(gene_model_all, junction_dir, intron_dir):
    junction_dict,exon_dict,gene_dict = importJunctionInfo(gene_model_all)
    importJunctionFromBED(junction_dir, junction_dict)


def importJunctionInfo(gene_model_all):
    '''
    This method will take gene_model as input and create three dictionaries
    1. junction_dict
    2. exon_dict
    3. gene_dict
    '''
    junction_dict = {}
    exon_dict = {}
    gene_dict = {}


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
    for gene in gene_dict:
        gene_dict[gene].sort()
    return exon_dict, gene_dict, junction_dict

   
def get_exon_annotation(chr, start,stop,exonid):
    if(exonid.startswith('E')):
        return exonid + '.' + str(chr) + '_' + str(start) + '_'+  str(stop)


def importJunctionFromBED(junction_dir, junction_dict):
    '''
    This function imports junction data from BAM derived exon junction text files
    '''
    junction_coordinate_BAM_dict = {}
    for junctionfilename in os.listdir(junction_dir):
        junctionfile = os.path.join(junction_dir, junctionfilename)
        initialBAMJunctionImport(junctionfile)


def initialBAMJunctionImport(junctionfile):
    interimJunctionCoordinates={}
    junction_df = pd.read_csv(junctionfile,sep='\t',header=None, names=[
            "chr", "start", "stop","annotation", "counts","strand"])
    for idx,row in junction_df.iterrows():
        #check if it's chrfehrkh
        if(len(row.chr) <= 5):
            interimJunctionCoordinates[(row.chr,row.start,row.stop)] = []
    
    print(interimJunctionCoordinates)
    return interimJunctionCoordinates
    


if __name__ == '__main__':
    genemodelpath = "/Users/sin9gp/altanalyze3/tests/data/" 
    file = 'gene_model_ENSG00000223972.txt'
    gene_model_file = f"{genemodelpath}/{file}"
    junction_dir = "/Users/sin9gp/altanalyze3/tests/data/junction_dir/"
    intron_dir = "/Users/sin9gp/altanalyze3/tests/data/intron_dir/"
    sampleBEDFileAnnotation(gene_model_file, junction_dir, intron_dir)
    # importJunctionInfo(gene_model_all=gene_model_file)
    # importJunctionFromBED(junctiondir)