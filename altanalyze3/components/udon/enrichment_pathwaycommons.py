import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import math
from statsmodels.stats.multitest import multipletests


def call_pathway_database(species, database_file_path=None):

    # read gtf file
    if database_file_path is not None:
        # read gtf file
        pathways_df = read_database_file(database_file_path)
    else:
        # Filter based on species
        if species == "Hs":
            pathways_df = pd.read_csv('PathwayCommons-Hs.txt', sep="\t", names=["GeneID", "Pathway"], header=0)
            # Filter uppercase entries or those starting with "ENSG"
            pathways_df = pathways_df[(pathways_df['GeneID'].str.isupper()) | (pathways_df['GeneID'].str.startswith('ENSG'))]
        elif species == "Mm":
            # Filter uppercase entries or those starting with "ENSG"
            pathways_df = pd.read_csv('PathwayCommons-Mm.txt', sep="\t", names=["GeneID", "Pathway"], header=0)
            pathways_df = pathways_df[~(pathways_df['GeneID'].str.isupper()) | (pathways_df['GeneID'].str.startswith('ENSMU'))]

    return pathways_df


def read_database_file(pathways_database_file_path):

    columns = ["GeneID", "Pathway"]
    df = pd.read_csv(pathways_database_file_path, sep='\t', comment='#', names=columns, dtype={"seqname": str}, header=0)

    return df


# Define functions for Z-score and Fisher's exact test
def calculate_z_score_aa(associated_in_group, in_static_interval, total, in_flexible_interval):
    r = float(associated_in_group)  # genes regulated in pathway
    _n = float(in_static_interval)  # genes regulated
    N = float(total)  # measured genes in the genome - total_count (#all genes evaluated on pathways)
    R = float(in_flexible_interval)  #genes in pathway
    if (R - N) == 0:
        return 0
    elif r == 0 and _n == 0:
        return 0
    else:
        try:
            # try:
            z = (r - _n * (R / N)) / math.sqrt(_n * (R / N) * (1 - (R / N)) * (1 - ((_n - 1) / (N - 1))))
            return z
            # except ZeroDivisionError: print 'r,_n,R,N: ', r,_n,R,N;kill
        except ValueError:
            print(r - _n * (R / N)), _n * (R / N) * (1 - (R / N)) * (1 - ((_n - 1) / (N - 1))), r, _n, N, R;kill
    # if expected == 0:
        # return 0
    # return (observed - expected) / np.sqrt(expected * (1 - (expected / total_genes)))


def calculate_z_score(observed, expected, total_genes):
    if expected == 0:
        return 0
    return (observed - expected) / np.sqrt(expected * (1 - (expected / total_genes)))


def fisher_test(contingency_table):

    _, p_value = fisher_exact(contingency_table, alternative='greater')
    return p_value


def pathway_enrichment(adata, pathways_df):

    # Load gene sets DataFrame
    gene_sets_df = adata.uns['udon_marker_genes_top_n']
    gene_sets_df = gene_sets_df[['marker', 'top_cluster']]

    # Convert DataFrame to dictionary
    query_gene_sets = gene_sets_df.groupby('top_cluster')['marker'].apply(list).to_dict()

    # Calculate pathway statistics for each query set
    results = []

    total_genes = pathways_df['GeneID'].nunique()
    total_pathways = pathways_df['Pathway'].nunique()

    for set_name, gene_list in query_gene_sets.items():
        query_gene_df = pathways_df[pathways_df['GeneID'].isin(gene_list)]
        pathway_counts = query_gene_df['Pathway'].value_counts()

        for pathway in pathways_df['Pathway'].unique():
            observed = pathway_counts.get(pathway, 0)
            total_in_pathway = pathways_df[pathways_df['Pathway'] == pathway]['GeneID'].nunique()
            expected = (total_in_pathway / total_genes) * len(gene_list)

            z_score = calculate_z_score(observed, expected, total_genes)

            # Contingency table for Fisher's exact test
            not_in_pathway = total_genes - total_in_pathway
            in_query_set_not_in_pathway = len(gene_list) - observed
            contingency_table = [[observed, total_in_pathway - observed],
                                 [in_query_set_not_in_pathway, not_in_pathway - in_query_set_not_in_pathway]]

            p_value = fisher_test(contingency_table)
            results.append({
                'query_set': set_name,
                'pathway': pathway,
                'observed': observed,
                'expected': expected,
                'z_score': z_score,
                'p_value': p_value
            })


    # Save the results
    results_df = pd.DataFrame(results)

    return results_df


def enrichment_wrapper(adata, species, database_file_path=None, p_val=0.05):

    # Load the database
    pathways_df = call_pathway_database(species, database_file_path)

    # Perform pathway enrichment analysis
    results_df = pathway_enrichment(adata, pathways_df)

    # Apply Benjamini-Hochberg FDR correction within each query gene set
    def apply_fdr_correction(group):
        p_values = group['p_value'].values
        corrected_p_values = multipletests(p_values, method='fdr_bh')[1]
        group['adj_p_value'] = corrected_p_values
        return group

    # Group by 'query_gene_set' and apply the correction
    results_df = results_df.groupby('query_set').apply(apply_fdr_correction).reset_index(drop=True)
    results_df_clean = results_df[results_df['p_value'] < p_val]
    results_df_clean = results_df_clean[results_df_clean['observed'] > 2]

    # Function to get top 5 entries with p_val < 0.05
    def top_n_entries(group, n=5):
        # Filter entries with p_val < 0.05 and get top n by p_val
        return group.nlargest(n, 'z_score')

    # Apply the function to each group
    top_entries = results_df_clean.groupby('query_set').apply(top_n_entries).reset_index(drop=True)

    top_entries = top_entries[top_entries['observed'] > 2]

    adata.uns['significant_pathway_enrichment'] = results_df
    adata.uns['significant_pathway_enrichment_top_n'] = top_entries

    return adata
