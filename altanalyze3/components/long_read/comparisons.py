import os, sys, csv, copy
import pandas as pd
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_automate as isoa
from ..annotation import junction_isoform as ji
from ..oncosplice import metadataAnalysis as ma
import numpy as np

def compare_one_cluster_to_many(condition,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction'):
    # Restrict by a specific condition (e.g., aged, young)
    print(f'Performing differential {dataType} analysis for every cluster within a single condition...')
    transcript_assoc_file = 'gff-output/transcript_associations.txt'
    protein_summary_file = 'protein_summary.txt'
    if dataType == 'junction':
        stats_folder = 'dPSI-cluster-events/'
    if dataType == 'isoform':
        stats_folder = 'diff-cluster-isoform/'
    if dataType == 'isoform-ratio':
        stats_folder = 'diff-cluster-ratio/'

    for cluster in cluster_order:
        # Filter samples for the current cluster
        current_cluster_samples = groups_file[
            (groups_file.index.str.startswith(cluster)) &  # Match the specific cluster
            (groups_file['grp'].str.endswith(condition))  # Ensure it belongs to the queried condition
        ]

        # Filter other samples within the same condition but in different clusters
        other_samples = groups_file[
            (~groups_file.index.str.startswith(cluster)) &  # Exclude the current cluster
            (groups_file['grp'].str.endswith(condition))  # Ensure it belongs to the queried condition
        ]

        # Label the groups appropriately
        current_cluster_samples = current_cluster_samples.copy()  # Optional: Ensure it is a copy
        current_cluster_samples.loc[:, 'grp'] = cluster  # Label as the current cluster (e.g., 'HSC-1')
        other_samples = other_samples.copy()  # Optional: Make a copy to avoid warnings
        other_samples.loc[:, 'grp'] = 'Others'  # Label as 'Others'

        # Combine the two groups
        filtered_groups_file = pd.concat([current_cluster_samples, other_samples])
        # Proceed with analysis if both groups have samples
        if not current_cluster_samples.empty and not other_samples.empty:
            print(f"Performing analysis for cluster: {cluster}")
            run_metadataAnalysis(cluster,psi_matrix,cluster,'Others',filtered_groups_file,stats_folder,filter_condition=condition)

    print ('Annotating the test results by gene, isoform and protein...')
    ji.annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType=dataType)

def compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction'):
    print (f'Performing differential {dataType} analysis between two user supplied conditions...')
    transcript_assoc_file = 'gff-output/transcript_associations.txt'
    protein_summary_file = 'protein_summary.txt'
    if dataType == 'junction':
        stats_folder = 'dPSI-covariate-events/'
    if dataType == 'isoform':
        stats_folder = 'diff-covariate-isoform/'
    if dataType == 'isoform-ratio':
        stats_folder = 'diff-covariate-ratio/'
        
    for cluster in cluster_order:
        filtered_groups_file = groups_file[
            (groups_file.index.str.startswith(cluster)) &  # Filter by cell type prefix
            (groups_file['grp'].str.endswith((condition1, condition2)))  # Keep only 'aged' and 'young' groups
        ]
        try: run_metadataAnalysis(cluster,psi_matrix,condition1,condition2,filtered_groups_file,stats_folder)
        except:
            print (f'Failed for {condition1} vs. {condition2} for {cluster}')
    print ('Annotating the test results by gene, isoform and protein...')
    ji.annotate(gene_symbol_file,transcript_assoc_file,junction_coords_file,protein_summary_file,stats_folder,dataType=dataType)

def run_metadataAnalysis(cluster,psi_matrix,condition1,condition2,filtered_groups_file,stats_folder,filter_condition=False):
    matching_uids = filtered_groups_file.index.intersection(psi_matrix.columns)
    # Subset the PSI matrix to only include the matching sample columns
    subset_psi_matrix = psi_matrix[matching_uids]
    print (f'Exporting {condition1} vs. {condition2} for {cluster}')
    diff_stats = ma.mwuCompute(subset_psi_matrix, filtered_groups_file, grpvar='grp', min_group_size=2)

    if filter_condition == False:
        stats_file = f'{stats_folder}/{condition1}-{condition2}-{cluster}_stats.txt'
    else:
        stats_file = f'{stats_folder}/{filter_condition}-{condition1}-{condition2}-{cluster}_stats.txt'
    os.makedirs(os.path.dirname(stats_file), exist_ok=True)    
    diff_stats.to_csv(stats_file, sep='\t', index=True)

def compute_differentials(sample_dict,conditions,cluster_order,gene_symbol_file,analyses=['junction','isoform','isoform-ratio']):
    analyzed_intial_conditions = []
    for pair in conditions:
        (condition1,condition2) = pair
        print (f'Analyzing {condition2} vs. {condition1} ')

        if 'junction' in analyses:
            psi_dir = 'psi_combined_pseudo_cluster_counts.txt'
            junction_coords_file = 'junction_combined_pseudo_cluster_counts-filtered.txt'

            psi_matrix = pd.read_csv(psi_dir, sep='\t', index_col=0)
            print ('PSI matrix imported')

            # Extract 'groups' for each UID from the sample_dict
            uid_to_group = isoa.get_sample_to_group(sample_dict,'junction')

            # Create a DataFrame for the group mapping using the PSI matrix columns
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(col.split('.')[1], 'unknown') for col in psi_matrix.columns],  # Extract UID and map to group
                'uid': psi_matrix.columns
            }).set_index('uid')

            # Compute Differential Splicing between groups and cell types
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction')
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,psi_matrix,junction_coords_file,gene_symbol_file,dataType='junction')

        if 'isoform' in analyses:
            # Compute Differential Isoform abundance between groups and cell types
            isoform_dir = 'isoform_combined_pseudo_cluster_tpm-filtered.txt'
            isoform_matrix = pd.read_csv(isoform_dir, sep='\t', index_col=0)
            log2_iso_matrix = np.log2(isoform_matrix + 1)
            uid_to_group = isoa.get_sample_to_group(sample_dict,'isoform')
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(col.split('.')[1], 'unknown') for col in isoform_matrix.columns],  # Extract UID and map to group
                'uid': isoform_matrix.columns
            }).set_index('uid')
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,log2_iso_matrix,None,gene_symbol_file,dataType='isoform')
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,log2_iso_matrix,None,gene_symbol_file,dataType='isoform')

        if 'isoform-ratio' in analyses:
            isoform_dir = 'isoform_combined_pseudo_cluster_ratio-filtered.txt'
            isoform_matrix = pd.read_csv(isoform_dir, sep='\t', index_col=0)
            uid_to_group = isoa.get_sample_to_group(sample_dict,'isoform')
            groups_file = pd.DataFrame({
                'grp': [uid_to_group.get(col.split('.')[1], 'unknown') for col in isoform_matrix.columns],  # Extract UID and map to group
                'uid': isoform_matrix.columns
            }).set_index('uid')
            if condition1 not in analyzed_intial_conditions:
                compare_one_cluster_to_many(condition1,cluster_order,groups_file,isoform_matrix,None,gene_symbol_file,dataType='isoform-ratio')
            compare_two_groups_per_cluster(condition1,condition2,cluster_order,groups_file,isoform_matrix,None,gene_symbol_file,dataType='isoform-ratio')
            analyzed_intial_conditions.append(condition1)
