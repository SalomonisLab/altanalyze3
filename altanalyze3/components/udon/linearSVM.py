import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC
from tqdm import tqdm


def generate_train_data(fold_matrix_with_mf_genes, groups):
    # for each of the nmf clusters provided, compute the centroid for each cluster for the differential splicing events using the output_file arg (provided in the output file?)

    # Set the cell types and donors for whom to create pseudobulks
    selected_cell_types = pd.unique(groups['cluster'])

    # Initialize an empty dictionary to store results
    pseudobulk_dict = {}

    # Iterate over selected cell types
    for cell_type in tqdm(selected_cell_types, desc="generating train data using centroids"):
        # Filter cells based on cell type
        selected_cells = groups.index[groups['cluster'] == cell_type]

        # Calculate average gene expression for selected cells
        avg_expression = fold_matrix_with_mf_genes.loc[:, selected_cells].mean(axis=1)

        # Add average expression values to the result dictionary
        column_name = f"U{cell_type}"
        pseudobulk_dict[column_name] = avg_expression

    # Convert the dictionary to a DataFrame
    centroids = pd.DataFrame(pseudobulk_dict)

    return centroids


def classify(train, fold_matrix_with_mf_genes, groups):

    print("using random state 42...")
    # regression fit on the train dataset
    regression = LinearSVC(random_state=42)  # alternative is OneVsRestClassifier function in sklearn -- why not this one?
    regression.fit(train.transpose(), y=train.columns)  # train should be # of ranks by # of mf genes ("samples" as per python documentation are the ranks/clusters) and "features" are genes and groups is # of ranks/clusters by # of pseudobulks/samples

    prob = regression.decision_function(fold_matrix_with_mf_genes.transpose())  # it should be a # of samples/patients by # of features format so we are transposing it

    # Determine the cluster with the maximum decision value for each sample.
    # Binary case (2 clusters): decision_function is 1-D -> sign selects the class.
    # Use regression.classes_ for the canonical decision-function column order.
    if np.ndim(prob) == 1:
        pred = np.where(prob > 0, regression.classes_[1], regression.classes_[0])
        win_score = np.abs(prob)                                    # binary margin (>=0)
    else:
        pred = np.asarray(regression.classes_)[np.argmax(prob, axis=1)]
        win_score = np.max(prob, axis=1)                            # winning class decision value
    final_clusters = pd.DataFrame(pred, index=fold_matrix_with_mf_genes.columns, columns=['cluster'])
    # exclude low-confidence assignments: winning SVM decision score < 0
    keep = win_score >= 0
    n_drop = int((~keep).sum())
    if n_drop:
        print(f"linear_svm: excluding {n_drop} pseudobulks with decision score < 0 (low confidence)")
    final_clusters = final_clusters[keep].sort_values('cluster')

    return final_clusters

