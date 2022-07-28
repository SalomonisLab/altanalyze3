import pandas as pd
import os

#path for gene_model and junction files
genemodelpath = "/Users/sin9gp/altanalyze3/tests/data/"

def importJunctionInfo(gene_model_all):
    '''
    This method will take gene_model as input and create three dictionaries
    1. junction_dict
    2. exon_dict
    3. gene_dict
    '''
    junction_dict = {}

    gene_model_df = pd.read_csv(gene_model_all,sep='\t',header=None, names=[
            "gene_id", "chr", "strand","exon_region_id", "start", "stop", "exon_annotations"])    
    #junction_dict = gene_model_df.set_index('chr','start').to_dict()['gene_id']
    for idx,row in gene_model_df.iterrows():
        junction_dict[(row.chr,row.start)] = row.gene_id + ':' + row.exon_region_id
    print(junction_dict)
   


if __name__ == '__main__':

    file = 'gene_model_ENSG00000223972.txt'
    gene_model_file = f"{genemodelpath}/{file}"
    importJunctionInfo(gene_model_all=gene_model_file)