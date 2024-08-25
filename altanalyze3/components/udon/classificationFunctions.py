import math

import scanpy as sc
import numpy as np
from sklearn.svm import OneClassSVM
import scanpy.external as sce
from sklearn.ensemble import IsolationForest

from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity


def gene_cosine_similarity(train_adata, test_adata):

    # for plotting purposes
    n_cells_train = np.shape(train_adata)[0]
    sc.tl.pca(train_adata, n_comps=min(30, n_cells_train - 1))
    sc.tl.pca(test_adata, n_comps=min(30, n_cells_train - 1))

    # compute similarity matrix
    similarity_matrix = cosine_similarity(test_adata.X, train_adata.X)
    print("computed similarity matrix", flush=True)

    # average "similarity score" for each cell to all the controls
    row_mean_similarity = np.mean(similarity_matrix, axis=1)

    predicted_labels = np.zeros(np.shape(similarity_matrix)[0])
    predicted_labels[row_mean_similarity < np.mean(row_mean_similarity)] = 1

    return train_adata, test_adata, predicted_labels, row_mean_similarity


def one_class_svm(train_adata, test_adata, nu=0.1, matrix="gex"):

    model = OneClassSVM(gamma='scale', nu=nu)

    if matrix == "gex":

        model.fit(train_adata.X)
        y_pred = model.predict(test_adata.X)
        decision_values = model.decision_function(test_adata.X)

    elif matrix == "pca":
        n_cells_train = np.shape(train_adata)[0]
        sc.tl.pca(train_adata, n_comps=min(30, n_cells_train - 1))
        sc.tl.pca(test_adata, n_comps=min(30, n_cells_train - 1))

        pca_train = train_adata.obsm['X_pca']
        pca_test = test_adata.obsm['X_pca']

        model.fit(pca_train)
        y_pred = model.predict(pca_test)
        decision_values = model.decision_function(pca_test)

    elif matrix == "harmony_pca":
        n_cells_train = np.shape(train_adata)[0]
        sc.tl.pca(train_adata, n_comps=min(30, n_cells_train - 1))
        nclust_harmony = np.min([np.round(n_cells_train/ 30.0), 100]).astype(int)
        if nclust_harmony <= 1:
            nclust = 2
        else:
            nclust = None
        sce.pp.harmony_integrate(train_adata, 'patient_id', nclust=nclust)
        harmony_pca_train = train_adata.obsm['X_pca_harmony']

        sc.tl.pca(test_adata, n_comps=np.shape(harmony_pca_train)[1])
        pca_test = test_adata.obsm['X_pca']
        # harmony_pca_test = test_adata.obsm['X_pca_harmony']

        model.fit(harmony_pca_train)
        y_pred = model.predict(pca_test)
        decision_values = model.decision_function(pca_test)

    # convert the outliers (-1) to disease (1)
    y_pred[y_pred == 1] = 0
    y_pred[y_pred == -1] = 1
    predicted_labels = y_pred

    return train_adata, test_adata, predicted_labels, decision_values


def isolation_forest(train_adata, test_adata, nu=0.1, matrix="gex"):

    model = IsolationForest()

    if matrix == "gex":

        model.fit(train_adata.X)
        y_pred = model.predict(test_adata.X)
        decision_values = model.decision_function(test_adata.X)

    elif matrix == "pca":
        n_cells_train = np.shape(train_adata)[0]
        sc.tl.pca(train_adata, n_comps=min(30, n_cells_train - 1))
        sc.tl.pca(test_adata, n_comps=min(30, n_cells_train - 1))

        pca_train = train_adata.obsm['X_pca']
        pca_test = test_adata.obsm['X_pca']

        model.fit(pca_train)
        y_pred = model.predict(pca_test)
        decision_values = model.decision_function(pca_test)

    elif matrix == "harmony_pca":
        n_cells_train = np.shape(train_adata)[0]
        sc.tl.pca(train_adata, n_comps=min(30, n_cells_train - 1))
        nclust_harmony = np.min([np.round(n_cells_train/ 30.0), 100]).astype(int)
        if nclust_harmony <= 1:
            nclust = 2
        else:
            nclust = None
        sce.pp.harmony_integrate(train_adata, 'patient_id', nclust=nclust)
        harmony_pca_train = train_adata.obsm['X_pca_harmony']

        sc.tl.pca(test_adata, n_comps=np.shape(harmony_pca_train)[1])
        pca_test = test_adata.obsm['X_pca']
        # harmony_pca_test = test_adata.obsm['X_pca_harmony']

        model.fit(harmony_pca_train)
        y_pred = model.predict(pca_test)

    # convert the outliers (-1) to disease (1)
    y_pred[y_pred == 1] = 0
    y_pred[y_pred == -1] = 1
    predicted_labels = y_pred

    return predicted_labels, decision_values


def clf_metrics_all(true_labels, predicted_labels):

    # Compute the ROC curve
    fpr, tpr, _ = metrics.roc_curve(true_labels, predicted_labels)

    # Compute the AUC (Area Under the Curve) value
    auc = metrics.roc_auc_score(true_labels, predicted_labels)

    # Compute the confusion matrix
    conf_matrix = metrics.confusion_matrix(true_labels, predicted_labels)
    specificity_healthy = conf_matrix[0, 0]/(conf_matrix[0, 0] + conf_matrix[0, 1])
    sensitivity_disease = conf_matrix[1, 1] / (conf_matrix[1, 0] + conf_matrix[1, 1])

    # Generate the classification report
    classification_report = metrics.classification_report(true_labels, predicted_labels, target_names=['Healthy', 'Disease'])

    return fpr, tpr, auc, specificity_healthy, sensitivity_disease, conf_matrix, classification_report


def find_top_bottom_cells(metadata, decision_obs_col, top_perc=0.15, bottom_perc=0.15):

    metadata = metadata.sort_values(by=[decision_obs_col])
    metadata["new_predicted_labels"] = 2
    n_cells = np.shape(metadata)[0]
    top_n_cells = math.ceil(top_perc*n_cells)
    bottom_n_cells = math.ceil(bottom_perc*n_cells)

    metadata.iloc[:top_n_cells, 13] = 1   # disease
    metadata.iloc[-bottom_n_cells:, 13] = 0  # healthy

    return metadata
