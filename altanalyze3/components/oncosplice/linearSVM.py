import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC


def generate_train_data(psi_file_with_diff_events, nmf_binarized_output):

    # for each of the nmf clusters provided, compute the centroid for each cluster for the differential splicing events using the output_file arg (provided in the output file?)

    centroids = pd.DataFrame()
    for col in range(nmf_binarized_output.shape[1]):

        col_all_samples = nmf_binarized_output[:, col]
        # col_all_samples = col_all_samples.to_numpy()
        col_cluster_samples = np.where(col_all_samples == 1)  # this is a tuple which is ridiculously annoying. it looks like: (array([ 63,  65,  71, 234, 349, 364]),) so need I am extracting "first" element in the next line to compute centroid to get only the array which is indices.

        centroid_vec = psi_file_with_diff_events.iloc[:, col_cluster_samples[0]].mean(axis=1)  # row means (rows are events) for each group
        centroids[col] = centroid_vec
                    
    return centroids


def classify(train, imputed_psi_file_with_diff_events, groups, conservation="stringent"):

    # Y is median imputed psi file and X is "grplst"/groups in the original function by MV
    # regression fit on the train dataset
    regression = LinearSVC()  # alternative is OneVsRestClassifier function in sklearn -- why not this one?
    regression.fit(train.transpose(), groups)  # train should be # of ranks by # of differential splicing events ("samples" as per python documentation are the ranks/clusters) and "features" are splicing events and groups is # of ranks/clusters by # of patients/samples
    # in the previous regression fit, why do we choose just the first column (it refers to the first patient/sample only, so what does the reproducibility looks like here?
    # i checked this by selecting a column or a patient/sample that has exactly 1 cluster assignment when k =2 and the decision function values are not identical but they are similar so yes, there is some issues regarding reproducibility here. here is an example from 3 different columns for the groups[:,"mysterious column"]

    # prob[1:10]
    # Out[287]:
    # array([-0.21518692, 0.1017927, 0.00177943, 0.89216488, 0.2264813,
    #        1.24986209, 0.97950145, 1.44557026, 1.15520506])
    # prob4[1:10]
    # Out[288]:
    # array([-2.15928144e-01, 1.01157955e-01, 1.12518304e-03, 8.91793519e-01,
    #        2.25907636e-01, 1.24954443e+00, 9.79122635e-01, 1.44533371e+00,
    #        1.15487401e+00])
    # prob10[1:10]
    # Out[289]:
    # array([-0.21470678, 0.10214063, 0.00216086, 0.89218438, 0.22676602,
    #        1.24977117, 0.979506, 1.44538875, 1.15514289])

    prob = regression.decision_function(imputed_psi_file_with_diff_events.transpose())  # in the PSI file here, it should be a # of samples/patients by # of features format so we are transposing it

    if len(groups) > 2:
        final_clusters = prob > 0
        #print(prob)
    else:  # this is when NMF rank = 2;
        # if number of groups/rank=2 then choose the cluster for which the decision function score is greater than 0.5 or less than -0.5 ("stringent") (> 0.25 and < -0.25 in case of "conservative"/"expanded")
        if conservation == "stringent":
            decision_score_threshold = 0.5
        else:
            decision_score_threshold = 0.25
        final_cluster_1 = prob > decision_score_threshold
        final_cluster_2 = prob < -decision_score_threshold
        final_clusters = np.stack([final_cluster_1, final_cluster_2])
        final_clusters = final_clusters.transpose()

    final_clusters = final_clusters.astype(int)
    final_clusters = pd.DataFrame(final_clusters, index=imputed_psi_file_with_diff_events.columns)

    return final_clusters
