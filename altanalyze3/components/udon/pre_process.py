import pandas as pd
import numpy as np
from tqdm import tqdm


def calculate_pseudobulk_folds(uadata, groups='cell_types', donors='donor_id', disease_type_col='disease_type', min_cells=5, baseline='control', query='disease', control_donor_specificity='collective'):
    """
        Calculate pseudobulk fold values for each gene in the adata object

        :param uadata: input ann data which contains cell by gene gene expression in adata.X
        :param groups: the column name in adata.obs that contains the cell type cluster
        :param donors: the column name in adata.obs that contains the donor id of the cell
        :param disease_type_col: the column name in adata.obs that contains which cells are “Disease” and which cells are “Healthy” for fold calculation
        :param min_cells: minimum number of cells present for calculation of any pseudobulk vector
        :param baseline:in the disease_type_col, the character string that refers to the baseline (denominator) for pseudobulk fold calculation. For udon, the baseline are cells labeled as “Control” in the disease_type_col
        :param query: : in the disease_type_col, the character string that refers to the query (numerator) for pseudobulk fold calculation. For udon, the baseline are cells labeled as “Disease” in the disease_type_col
        :param control_donor_specificity: how the baseline for pseudobulk-fold should be calculated? If “collective” (default), then the baseline (typically the control samples for UDON) will be calculated agnostic of donors.
        :return: adata with slot adata.varm['pseudobulk_folds'] containing the pseudobulk fold values

        """

    # separate the cells from control and disease samples
    baseline_adata = uadata[uadata.obs[disease_type_col] == baseline, :]
    disease_adata = uadata[uadata.obs[disease_type_col] == query, :]

    # finalize which clusters and samples to include for pseudobulk calculation
    baseline_qc = pre_pseudobulk_qc(baseline_adata, groups=groups, donors=donors, min_cells=min_cells)
    disease_qc = pre_pseudobulk_qc(disease_adata, groups=groups, donors=donors, min_cells=min_cells)

    # if there are differences in cell clusters between the control and disease, inform the user, but carry on with the assumption one group has 0 expression
    # Find the intersection
    common_elements = set(baseline_qc[groups]) & set(disease_qc[groups])

    # Filter the DataFrames to only include common elements
    baseline_qc = baseline_qc[baseline_qc[groups].isin(common_elements)]
    disease_qc = disease_qc[disease_qc[groups].isin(common_elements)]

    # calculate the pseudobulk for each donor and each cell type in each group (baseline and query)
    baseline_pseudobulks_mat = calculate_pseudobulks(adata=baseline_adata, pseudobulk_groups=baseline_qc, clusters=groups, donors=donors, donor_specificity=control_donor_specificity)
    disease_pseudobulks_mat = calculate_pseudobulks(adata=disease_adata, pseudobulk_groups=disease_qc, clusters=groups, donors=donors, donor_specificity='donor-specific')

    # calculate psueodbulk folds between baseline and query
    pseudobulk_folds = calculate_folds(baseline_pseudobulks_df=baseline_pseudobulks_mat, disease_pseudobulks_df=disease_pseudobulks_mat)

    # finally, ensure the entries in the pseudobulk folds are positive for ICGS
    non_neg_folds = remove_negatives(pseudobulk_folds=pseudobulk_folds)

    # remove columns with missing values in non_neg_folds
    non_neg_folds = non_neg_folds.dropna(axis=1)

    # store the pseudobulk folds in the adata object in the varp slot
    uadata.varm['pseudobulk_folds'] = non_neg_folds

    return uadata


def remove_negatives(pseudobulk_folds):

    # Find the minimum values across each row
    min_vals = pseudobulk_folds.min(axis=1)

    # Subtract min_vals from each column
    for i in range(pseudobulk_folds.shape[1]):
        pseudobulk_folds.iloc[:, i] = pseudobulk_folds.iloc[:, i] - min_vals

    return pseudobulk_folds


def pre_pseudobulk_qc(adata, groups='cell_types', donors='donor_id', min_cells=5):

    baseline_qc = adata.obs.groupby([donors, groups], observed=True).size().reset_index(
        name='cell_count')

    # Identify and store dropped rows
    dropped_rows = baseline_qc[baseline_qc['cell_count'] < min_cells]

    # Write out which rows were dropped
    if len(dropped_rows) > 0:
        print(f"Rows dropped due to cell count < {min_cells}:")
        print(dropped_rows)

    # Filter out rows where 'cell_count' is less than 5
    baseline_qc = baseline_qc[baseline_qc['cell_count'] >= min_cells]

    return baseline_qc


def calculate_pseudobulks(adata, pseudobulk_groups, clusters='cell_type', donors='donor_id', donor_specificity='donor-specific'):
    """
    Calculate the pseudobulks for the given cell types and donors
    :param adata: input AnnData object containing the gene expression matrix (adata.X) for which to compute the pseudobulks
    :param pseudobulk_groups: a DataFrame containing the cells and the associated cell types for which to compute pseudobulk values for.
    this is a good place to select the exact cells to include you want to compute pseudobulks for
    :param clusters: the column name in adata.obs AND in pseudobulk_groups that contains the cell type cluster
    :param donors: the column name in adata.obs that contains the donor id of the cell
    :param donor_specificity: how the pseudobulks should be calculated. If 'donor-specific', then the pseudobulks will be calculated for each cell type and donor combination.
    If 'collective', then the pseudobulks will be calculated for each cell type, regardless of donor.
    :return: a DataFrame containing the pseudobulk values for the given cell types and donors
    """

    # Filter for the cell types and donors for which to compute the pseudobulks
    adata = adata[adata.obs[clusters].isin(pseudobulk_groups[clusters])]

    # Extract the gene expression matrix directly from adata
    gex_mat = adata.to_df()

    # Set the cell types and donors for whom to create pseudobulks
    selected_cell_types = pd.unique(pseudobulk_groups[clusters])
    selected_donors = pd.unique(pseudobulk_groups[donors])

    # Initialize an empty dictionary to store results
    pseudobulk_dict = {}

    if donor_specificity == 'donor-specific':
        # Iterate over selected cell types and donors
        for cell_type in tqdm(selected_cell_types, desc="Calculating donor-specific pseudobulks"):
            for donor in selected_donors:
                # Filter cells based on cell type and donor
                selected_cells = adata.obs_names[(adata.obs[clusters] == cell_type) & (adata.obs[donors] == donor)]

                # Calculate average gene expression for selected cells
                avg_expression = gex_mat.loc[selected_cells].mean(axis=0)

                # Add average expression values to the result dictionary
                column_name = f"{cell_type}__{donor}"
                pseudobulk_dict[column_name] = avg_expression

    elif donor_specificity == 'collective':
        # Iterate over selected cell types
        for cell_type in tqdm(selected_cell_types, desc="Calculating collective pseudobulks"):
            # Filter cells based on cell type
            selected_cells = adata.obs_names[adata.obs[clusters] == cell_type]

            # Calculate average gene expression for selected cells
            avg_expression = gex_mat.loc[selected_cells].mean(axis=0)

            # Add average expression values to the result dictionary
            column_name = f"{cell_type}__collective"
            pseudobulk_dict[column_name] = avg_expression

    # Convert the dictionary to a DataFrame
    pseudobulk_matrix = pd.DataFrame(pseudobulk_dict)

    return pseudobulk_matrix

def calculate_folds(baseline_pseudobulks_df, disease_pseudobulks_df):

    # Initialize an empty DataFrame to store results
    folds_matrix = pd.DataFrame(index=disease_pseudobulks_df.index)  # Rows for genes
    baseline_cols = pd.Series(baseline_pseudobulks_df.columns)

    for col in disease_pseudobulks_df.columns:
        celltype_match = col.split('__')[0]
        matching_bcol = baseline_cols.str.contains(celltype_match, na=False)

        if matching_bcol.sum() == 0:
            print(f"Cell type {celltype_match} not found in baseline pseudobulks")
            continue
        baseline_mask = baseline_pseudobulks_df.loc[:, matching_bcol.to_list()]
        # join the difference vector as a column in the folds_matrix
        folds_matrix = pd.concat([folds_matrix, disease_pseudobulks_df[col] - baseline_mask.iloc[:, 0]], axis=1)
        folds_matrix = folds_matrix.set_axis([*folds_matrix.columns[:-1], col], axis=1)
        # folds_matrix[col] = disease_pseudobulks_df[col] - baseline_mask.iloc[:, 0]

    return folds_matrix
