
import numpy as np
import nimfa

# write the pseudocode here ##
'''
Notes:
1. the input for the nmf analysis cannot contain NA values so we use imputed PSI matrix 
2. Doing NMF on a matrix (in this case the imputed PSI matrix) decomposes the input matrix into W (basis matrix of size # of samples by rank/components) and H (coefficient matrix) matrices 
3. we will consider W matrix for all downstream analyses


Main steps in the run_nmf function: 
1. Run NMF and get the W matrix with "rank" parameter indicating the number of clusters NMF should determine. The rank parameter is determined before this analysis. 
2. Binarize W matrix in the following way so that it indicates the cluster membership of a sample (1 indicates the sample belongs to that NMF Cluster): 
    
3. Remove redundant components/clusters in the following way:
4. Finally, refine the clusters in either a "stringent" or a "conservative" way:
    a. if "stringent" 
    b. if "conservative":
'''


def run_nmf(imputed_psi_file, rank, strategy="stringent", sum_threshold=0.10):

    imputed_psi_file = imputed_psi_file.to_numpy()
    imputed_psi_file_t = imputed_psi_file.transpose()
    # can values be greater than 1 after imp? no
    # sd cannot be 0 across all rows/cols? ##what does nimfa snmf do when some columns or even rows are all 0s? ## see source for their checks
    # enter try statement here
    w = None
    try:
        nmf = nimfa.Snmf(imputed_psi_file_t, seed="nndsvd", rank=int(rank), max_iter=20, n_run=5, track_factor=True)  # n_run was originally 10 to keep the results stable
        nmf_fit = nmf()
        w = nmf_fit.basis()
        # w_a = np.array(w)  # maybe update the variable name here???
        # h = nmf_fit.coef()
        # H = np.array(h)
    except ValueError: 
        w = imputed_psi_file # this needs to change -- how are you taking care of errors?

    finally:
        nmf_matrix = w

    if int(rank) == 2:
        binarized_output = binarize_nmf(w, sum_threshold, par=1, lower_bound=3, upper_bound=40)  # adjusted way of determining outliers with par = 1
    else:
        binarized_output = binarize_nmf(w, sum_threshold, par=2, lower_bound=3, upper_bound=40)  # standard way of determining the outliers

    binarized_output_for_refining = np.array(binarized_output)
    
    # what if binarized output is all 0's? can reduce the the sum_threshold in the except statement within the main function call or if it is user provided, print the info and let the user know about it
    refined_binarized_output = remove_redundant_components(binarized_output_for_refining)

    # refined_binarized_output = finalize_binarized_output(refined_binarized_output, strategy=strategy) # this line not needed here -- MV does nothing in this function in original code
    # ...contd. from MV's dissertation, I believe this is a part of the ExpandSampleClusters/LinearSVC stuff so ignore it here for sure.

    return nmf_matrix, binarized_output, refined_binarized_output


def binarize_nmf(w, sum_threshold, par=1, lower_bound=3, upper_bound=40):

    w = w.transpose()
    # components = w.shape[0]  # this is number of rows/components
    samples = w.shape[1]  # number of columns/samples
    
    z_sum_threshold = np.where(w > sum_threshold, 1, 0)  # binarize W matrix based on sum_threshold (greater than sum_threshold (0.1) means 1)
    num_vec = np.sum(z_sum_threshold, axis=1)  # for every component/row, number of samples that have a value greater than sum_threshold (0.1)

    compstd_vec_criteria1 = num_vec < lower_bound
    compstd_vec_criteria2 = num_vec > upper_bound
    compstd_vec = np.where(np.logical_or(compstd_vec_criteria1, compstd_vec_criteria2) == True) # 40 can be a fraction of the num value #compstd is True if value in num_vec is < 3 or > 40

    # compstd_vec = np.where((num_vec < 3 or num_vec > 40), True, False)
    subset_compstd_w = w[compstd_vec, :].reshape(np.shape(compstd_vec)[1], samples)  # reshape it into a matrix of # of components (where compstd_vec == TRUE) by # of samples
    mean_vec = np.mean(subset_compstd_w, axis=1).transpose()
    st_vec = np.std(subset_compstd_w, axis=1).transpose()

    for i in range(samples):
        #print(i)
        z_sum_threshold[compstd_vec, i] = np.where((w[compstd_vec, i] >= mean_vec + par*st_vec), 1, 0)

    binarized_output = z_sum_threshold

    return binarized_output


def remove_redundant_components(binarized_output):

    components = binarized_output.shape[0]  # this is number of rows/components
    refined_binarized_output = []

    for i in range(components):
        component1 = binarized_output[i, :]
        sum1 = sum(component1) # number of samples belonging to the component/cluster
        flag = False
        component1_membership = np.where(component1 == 1) # samples that belong to the i component
        for j in range(components):
            if i != j:
                component2 = binarized_output[j, :]
                
                sum_of_overlapping_labels = sum(component2[component1_membership])
                sum2 = sum(component2)
                # print(float(sum_of_overlapping_labels) / sum1)
                if float(sum_of_overlapping_labels)/sum1 > 0.5:
                    if sum2 > sum1:
                        flag = True

        if not flag:
            refined_binarized_output.append(component1)

    refined_binarized_output = np.array(refined_binarized_output)

    return refined_binarized_output


def finalize_binarized_output(refined_binarized_output, strategy="stringent"):

    if strategy == "stringent":
        print("stringent")
    else:
        print("conservative")

    return refined_binarized_output
