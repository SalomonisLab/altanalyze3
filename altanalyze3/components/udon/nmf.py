import pandas as pd
import numpy as np
import nimfa
from sklearn.preprocessing import scale
from numpy import linalg as LA


# determine the rank of the NMF decomposition
def determine_nmf_ranks(df):
    # "To estimate the rank of the matrix (i.e. clusters) for SNMF, the ICGS Guide3 matrix is z-score normalized and its eigenvalues are calculated."

    # Convert X to a numpy array
    X = np.array(df)
    g = float(X.shape[0])  # Number of rows/genes
    c = float(X.shape[1])  # Number of columns/pseudobulks
    print(g)
    print(c)

    # Scale the data (zero mean and unit variance)
    X = scale(X)
    Xt = np.transpose(X)  # Transpose of X/df

    # Calculate muTW and sigmaTW for Tracy-Widom distribution
    muTW = (np.sqrt(g - 1) + np.sqrt(c)) ** 2.0
    sigmaTW = (np.sqrt(g - 1) + np.sqrt(c)) * (1.0 / np.sqrt(g - 1) + 1.0 / np.sqrt(c)) ** (1.0 / 3.0)

    # Compute the covariance matrix
    sigmaHat = np.dot(Xt, X)

    # Calculate the threshold boundary
    boundary = 3.273 * sigmaTW + muTW
    print(boundary)

    # Compute the eigenvalues of the covariance matrix
    w, v = LA.eig(sigmaHat)
    w = w.tolist()

    k = 0
    for i in range(len(w)):
        try:
            # Count the number of eigenvalues greater than the boundary
            if w[i] > boundary:
                k += 1
        except Exception:
            if w[i].real > boundary:
                k += 1

    est_k = 2*k   # Estimate the rank of the matrix

    return est_k


# write the pseudocode here ##
'''
Notes:
1. the input for the nmf analysis cannot contain NA values 
2. Doing NMF on a matrix (in this case the imputed PSI matrix) decomposes the input matrix into W (basis matrix of size # of samples by rank/components) and H (coefficient matrix) matrices 
3. we will consider W matrix for all downstream analyses

Main steps in the run_nmf function: 
1. Run NMF and get the W matrix with "rank" parameter indicating the number of clusters NMF should determine. The rank parameter is determined before this analysis. 

'''


def run_nmf(df, rank):
    mat = df.to_numpy()
    mat_t = mat.transpose()

    # sd cannot be 0 across all rows/cols? ##what does nimfa snmf do when some columns or even rows are all 0s? ## see source for their checks
    # enter try statement here
    w = None
    try:
        nmf = nimfa.Snmf(mat_t, seed="nndsvd", rank=int(rank), max_iter=20, n_run=5,
                         track_factor=True)  # n_run was originally 10 to keep the results stable
        nmf_fit = nmf()
        w = nmf_fit.basis()
    except ValueError:
        w = mat  # this needs to change -- how are errors handled?

    finally:
        nmf_matrix = pd.DataFrame(w.transpose(), columns=df.columns)

    nmf_clusters = binarize_nmf(w)
    nmf_clusters.index = df.columns

    return nmf_matrix, nmf_clusters


def binarize_nmf(w):
    w = w.transpose()
    # components = w.shape[0]  # this is number of rows/components
    samples = w.shape[1]  # number of columns/samples

    max_decision_value = np.zeros((samples, 1))

    for i in range(samples):
        # print(i)
        # i-th entry in max_decision_value is the row number with the maximum value in the i-th column of w
        max_decision_value[i] = np.argmax(w[:, i])

    nmf_clusters = pd.DataFrame(data=max_decision_value, columns=['cluster'], dtype='int64')

    return nmf_clusters
