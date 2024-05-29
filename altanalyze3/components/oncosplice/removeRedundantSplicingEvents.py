
import pandas as pd
import numpy as np


def format_psi_matrix(psi_file_path, n=11):

    psi_file = pd.read_csv(psi_file_path, sep="\t", header=0, index_col="UID")
    metadata = psi_file.iloc[:, 0:n]  # this needs to change maybe because not all PSI files have exact 11 columns of metadata
    formatted_psi_file = psi_file.iloc[:, n:]

    return formatted_psi_file, metadata


def remove_redundant_splicing_events(formatted_psi_file, corr_threshold, metadata):
    # print("corr.threshold for corr level being tested: " + str(corr_threshold))

    # create an empty list which will be appended with every gene run
    filtered_rows = []

    # get the list of genes
    rows_withonlygene = metadata.iloc[:, 0]
    rows_withonlygene = rows_withonlygene + ":"
    rows_withonlygene_np = rows_withonlygene.to_numpy()
    gene_list = pd.unique(rows_withonlygene)

    mat = formatted_psi_file

    for i in gene_list:
        # print(i)
        logical_vec = np.isin(rows_withonlygene_np, i)  # determine which are the rows in the PSI matrix associated with that gene
        gene_i_matrix = mat.iloc[logical_vec, ]  # get the events associated with the gene from PSI matrix
        gene_i_matrix = gene_i_matrix.transpose()  # transpose the filtered PSI matrix for the cor function
        cor_i_matrix = gene_i_matrix.corr(method="pearson")  # compute correlation which takes care of NA values
        cor_i_matrix_np = cor_i_matrix.to_numpy()  # to work with downstream steps, we need np array of cor matrix
        cor_i_matrix_np[np.triu_indices(cor_i_matrix_np.shape[0], 0)] = np.nan  # set upper triangular cor matrix to be NA because its symmetrical matrix
        idx_cor_i_binary = np.absolute(cor_i_matrix_np) > corr_threshold  # find indices/logical np array for values above threshold for cor values
        idx_cor_i_binary = 1 * idx_cor_i_binary  # binarize T/F matrix for easy column/row sum calculations
        sum_idx_cor_i_binary = np.nansum(idx_cor_i_binary, axis=0)  # determine the column sum of the binarized matrix
        # greater_than_zero_sums = sum_idx_cor_i_binary[sum_idx_cor_i_binary > 0]

        if any(sum_idx_cor_i_binary > 0):  # check if there are any correlated events,i.e., if any col sum is > 0
            if any(sum_idx_cor_i_binary == 1):  # check if there are exactly two events correlated with each other
                for col_idx in np.where(sum_idx_cor_i_binary == 1)[0]:  # from col idx we get ONE event out of the 2 correlated events
                    row_idx = np.where(idx_cor_i_binary[:, col_idx] == 1)[0][0]  # get the second event correlated with col idx event

                    num_of_nas_col_idx = gene_i_matrix.iloc[:, col_idx].isna().sum()  # determine # of NAs for the col idx event
                    num_of_nas_row_idx = gene_i_matrix.iloc[:, row_idx].isna().sum()  # determine # of NAs for the row idx event

                    number_of_nas = np.array([num_of_nas_col_idx, num_of_nas_row_idx])
                    event_to_keep = gene_i_matrix.columns.values[col_idx] if np.argmin(number_of_nas) == 0 else gene_i_matrix.columns.values[row_idx]  # keep the event with less # of NAs

                    # filtered_rows_i = event_to_keep #keep that event for this run
                    # filtered_rows = [filtered_rows, filtered_rows_i] # add that event to the initialized empty list
                    filtered_rows.append(event_to_keep)
                sum_idx_cor_i_binary = sum_idx_cor_i_binary[:len(sum_idx_cor_i_binary) - 1]  # last col is always full of NAs but the binarized matrix shows NAs to be 0 so we exclude the last col
                events_to_keep = gene_i_matrix.columns.values[np.where(sum_idx_cor_i_binary == 0)[0]]  # two events with 0 value means they are mutually exclusive so we retain that/those event(s)
                for j in events_to_keep:  # DO FOR LOOP FOR EACH EVENT IN EVENTS TO KEEP AND DO APPEND FOR THAT I
                    filtered_rows.append(j)

            else:
                cols_to_remove = np.where(np.nansum(idx_cor_i_binary, axis=0) > 0)[0]  # set of highly correlated events
                rows_to_remove = np.where(np.nansum(idx_cor_i_binary, axis=1) > 0)[0]  # prev step would not include the last col so we do the same for rows

                events_to_compare = np.union1d(cols_to_remove, rows_to_remove)  # identify the collected set of highly correlated events

                number_of_nas = gene_i_matrix.iloc[:, events_to_compare].isna().sum()  # determine the # of NAs in those events
                idx_to_keep = np.argmin(number_of_nas)  # find the event with least # of NAs
                event_to_keep = gene_i_matrix.columns.values[events_to_compare[idx_to_keep]]  # keep that event as the representative of that set of cor events
                filtered_rows.append(event_to_keep)
                cols_to_keep = np.where(np.nansum(idx_cor_i_binary, axis=0) == 0)[0]  # determine mutually exclusive events from cols
                rows_to_keep = np.where(np.nansum(idx_cor_i_binary, axis=1) == 0)[0]  # determine mutually exclusive events from rows to include for last col
                idxs_to_keep = np.intersect1d(cols_to_keep, rows_to_keep)  # here we have intersect as opposed to union because otherwise it will include repetitive events from correlated events
                # only if an event has a sum of 0 across both row and columns will it be an exclusive event
                events_to_keep = gene_i_matrix.columns.values[idxs_to_keep]
                for k in events_to_keep:  # DO FOR LOOP FOR EACH EVENT IN EVENTS TO KEEP AND DO APPEND FOR THAT I
                    filtered_rows.append(k)
        else:  # if the column sum of binarized matrix is 0 for all events then all events are mutually exclusive so we save them
            filtered_rows_i = gene_i_matrix.columns.values.tolist()
            for m in filtered_rows_i:  # DO FOR LOOP FOR EACH EVENT IN EVENTS TO KEEP AND DO APPEND FOR THAT I
                filtered_rows.append(m)
    # print("Number of filtered rows: " + str(len(filtered_rows)-1))

    return filtered_rows
