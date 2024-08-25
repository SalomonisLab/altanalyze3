import pandas as pd
import os
from determineNMFRank import *
from runNMF import *
import time
from metadataAnalysis import *
from linearSVM import *
from correlationDepletion import *
from correlationDepletion_vectorized import *
from visualizations import *


def round_wrapper(filename, full_psi_file, full_imputed_psi_file, highly_variable_events, metadata, rank=None, min_group_size=5,
                  dPSI=0.1, dPSI_p_val=0.05, min_differential_events=100, top_n_differential_events=150,
                  conservation="stringent", strictness="easy", depletion_corr_threshold=0.3, write_files=True,
                  speed_corr="og"):
    # a wrapper function that does for ONE iteration 1) nmf analysis 2) finds unique differential splicing events per NMF component and 3) linear SVC and 4) depletes the events for next round
    # full_psi_file: splicing events by samples dataframe of all splicing events (events after 1) 75% presence of samples  and 2) removing redundant events using inter-correlation analysis)
    # full_imputed_file: splicing events by samples dataframe of all splicing events (events with 75% presence of samples) but missing values imputed with median values for each splicing events
    # highly_variable_events: list of splicing events that are a subset of the events in the full psi file/full imputed file (usually events after feature selection step)
    # rank: optional integer value for the number of clusters for NMF analysis. alternatively, if user wants to detect automatically then, ignore.
    # min_group_size: integer value indicating the minimum number of samples that should belong to each cluster for valid differential expression analysis for the events
    # dPSI: fold change threshold. if events have a dPSI above this value and are significant based on dPSI_p_val value, then they are considered as differentially expressed events for a given cluster
    # dPSI_p_val: adjusted p value to use as reference to determine differential events
    # min_differential_events: integer value. if a cluster has less than this number of differential events
    # conservation: "stringent" or "expanded" to determine decision function threshold in linearSVM classification when NMF rank=2 or # of valid clusters = 2.
    # depletion_corr_threshold: correlation threshold to use to determine which events to deplete at the end of the round
    # write_files: writes out all the relevant files in txt form if set to True. otherwise, saves time by not writing/overwriting files.

    if write_files is True:
        os.mkdir(filename)
        cd = os.getcwd()

    # filter the full psi files (one with na and one with median imputed) with the highly variable events (currently comes from the file PCAbasedFeatureSelection.py ran before this function)
    formatted_psi_file_hve = full_psi_file.loc[highly_variable_events, :]
    formatted_psi_file_imp_hve = full_imputed_psi_file.loc[highly_variable_events, :]

    # when rank is not 2 or not provided by the user (that they trust more), then automatically calculate number of clusters for NMF to determine (aka rank)
    # requires the determineNMFRank.py script
    if rank is None:
        k_results = determine_nmf_ranks(PSIdf=formatted_psi_file_imp_hve.transpose(), MADA=True, CORRINT=True, MOM=True,
                                        PCAFO=True, TLE=True, FS=True)
        est_k = round(pd.Series.median(k_results, skipna=True))
        rank = est_k

    # run and time the NMF analysis using the runNMF.py script
    start_time = time.time()
    print(("Running NMF analysis with rank set to " + str(rank) + "..."))
    basis_matrix, binarized_output_test, refined_binarized_output = run_nmf(imputed_psi_file=formatted_psi_file_imp_hve,
                                                                            rank=rank)
    elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
    print("--- %.3f minutes ---" % elapsed_time)

    # for metadata analysis, which determines the differential splicing events, format the binary groups file (an output from NMF analysis) to match the expectations of the mwuCompute function
    groups_file = pd.DataFrame(
        refined_binarized_output.transpose())  # dimensions should be number of patients/samples (rows) by number of clusters (cols)
    groups_file.index = formatted_psi_file_hve.columns
    groups_file.columns = groups_file.columns.map(str)
    groups_file_cols = groups_file.columns

    if write_files is True:
        initial_clusters_path = os.path.join(cd, filename, "initial_round_clusters.txt")
        groups_file.to_csv(initial_clusters_path, sep="\t")

    remove_cols_from_nmf_clusters = []
    # check whether the data being provided is in good format and matches the sanity requirements for determining differential splicing events
    for col in range(len(groups_file_cols)):
        # print(col)
        mwu_check = mwuChecks(full_psi_file, groups_file, grpvar=groups_file_cols[col], min_group_size=min_group_size)
        # print(mwu_check)
        if mwu_check == 0:
            # print(col)
            # print(groups_file_cols[col])
            print(("differential expression not compatible for cluster " + groups_file_cols[col]))
            groups_file = groups_file.drop(groups_file_cols[col], axis='columns')
            # print(groups_file.columns)
            # refined_binarized_output = np.delete(refined_binarized_output, col, axis=0)
            remove_cols_from_nmf_clusters.append(col)

    refined_binarized_output = np.delete(refined_binarized_output, remove_cols_from_nmf_clusters, axis=0)

    start_time = time.time()
    print(("Finding top " + str(top_n_differential_events) + " differential splicing events per NMF cluster..."))
    diff_events, remove_clusters, results_all = find_top_differential_events(exp_file=full_psi_file, groups_file=groups_file,
                                                                             metadata=metadata,
                                                                             output_loc=os.path.join(cd, filename),
                                                                             min_group_size=min_group_size, dPSI=dPSI,
                                                                             dPSI_p_val=dPSI_p_val,
                                                                             min_differential_events=min_differential_events,
                                                                             top_n=top_n_differential_events)
    elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
    print("--- %.3f minutes ---" % elapsed_time)

    start_time = time.time()
    print(("Determined " + str(len(diff_events)) + " unique differential events across all clusters"))
    psi_file_with_diff_events = full_psi_file.loc[diff_events,]
    psi_file_with_diff_events_imp = full_imputed_psi_file.loc[diff_events,]
    elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
    print("--- %.3f minutes ---" % elapsed_time)

    refined_binarized_output_f = np.delete(refined_binarized_output, remove_clusters,
                                           axis=0)  # remove clusters (rows) with number of differential splicing events less than threshold (100)

    # Linear SVM/SVC training dataset generation by creating centroids for the NMF clusters

    centroids = generate_train_data(psi_file_with_diff_events=psi_file_with_diff_events,
                                    nmf_binarized_output=refined_binarized_output_f.transpose())  # splicing events by clusters matrix
    centroids = centroids.dropna(axis=0)
    psi_file_with_diff_events_imp = psi_file_with_diff_events_imp.loc[centroids.index, :]

    start_time = time.time()
    print(("Running Linear SVM on final " + str(np.shape(refined_binarized_output_f)[0]) + " clusters with..."))
    final_clusters = classify(train=centroids, imputed_psi_file_with_diff_events=psi_file_with_diff_events_imp,
                              groups=np.arange(np.shape(refined_binarized_output_f)[0]),
                              conservation=conservation)  # need to figure out a way to give more meaningful cluster name instead of an integer after this step (such as C1, etc etc).
    elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
    print("--- %.3f minutes ---" % elapsed_time)

    print(("Depleting events using a correlation threshold of " + str(depletion_corr_threshold) + "..."))
    if speed_corr == "og":
        start_time = time.time()
        depleted_psi_file_after_round = deplete_events(nmf_basis_matrix=basis_matrix, full_psi_file=full_psi_file,
                                                       full_imputed_psi_file=full_imputed_psi_file, corr_threshold=depletion_corr_threshold, strictness=strictness, speed="og")
        elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
        print("--- %.3f minutes ---" % elapsed_time)
    elif speed_corr == "efficient":
        start_time = time.time()
        depleted_psi_file_after_round = deplete_events(nmf_basis_matrix=basis_matrix, full_psi_file=full_psi_file,
                                                       full_imputed_psi_file=full_imputed_psi_file, corr_threshold=depletion_corr_threshold, strictness=strictness, speed="vectorized")
        elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
        print("--- %.3f minutes ---" % elapsed_time)

    depleted_psi_file_after_round_imp = full_imputed_psi_file.filter(depleted_psi_file_after_round.index, axis="rows")
    depleted_psi_file_after_round = full_psi_file.filter(depleted_psi_file_after_round.index, axis="rows")

    if write_files is True:
        diff_events_path = os.path.join(cd, filename, "input_linear_svm_features.txt")
        final_clusters_path = os.path.join(cd, filename, "round_clusters.txt")
        diff_events.to_csv(diff_events_path, sep="\t")
        final_clusters.to_csv(final_clusters_path, sep="\t")

    return final_clusters, depleted_psi_file_after_round_imp, depleted_psi_file_after_round, results_all
