import pandas as pd

from RoundWrapper import *
from preprocess import *
from removeRedundantSplicingEvents import *
from PCAbasedFeatureSelection import *
from medianImpute import *
from visualizations import *
import os
import argparse
import time
import matplotlib.pyplot as plt
import warnings
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

def oncosplice(args):
    psi_file_path = args.psi_file
    n_metadata_cols = args.n_metadata_cols
    corr_threshold = args.remove_redundant_events_threshold
    wd = args.wd
    dir = args.dir
    pca_corr_threshold = args.pca_threshold
    npcs = args.npcs
    rank_str = args.rank
    if rank_str == "k30":
        rank = 30
    else:
        rank = None
    print("NMF Rank used:")
    print(rank)
    rank_round1 = args.force_broad
    if rank_round1 == "on":
        rank1 = 2
    else:
        rank1 = rank
    print("Force broad used:")
    print(rank1)
    min_group_size = args.min_group_size
    dPSI = args.dpsi
    dPSI_p_val = args.dpsi_pval
    min_differential_events = args.min_differential_events
    top_n_differential_events = args.top_n_events
    conservation = args.conservation
    depletion_corr_threshold = args.depletion_corr
    speed = args.speed
    n_rounds = args.n_rounds
    fold_threshold = args.fold_threshold  # 0.3
    samples_differing = args.samples_differing  # 4
    corr_threshold_intercorr = args.corr_feature_selection  # 0.2
    corr_n_events = args.corr_n_events  # 10 or 5
    print("number of rounds to iterate through:")
    print(n_rounds)

    # read the psi file in and separate out the metadata columns from numeric PSI data
    formatted_psi_file, metadata = format_psi_matrix(psi_file_path, n=n_metadata_cols)
    metadata.index = formatted_psi_file.index
    formatted_psi_file = formatted_psi_file.dropna(how='all')
    metadata = metadata.loc[formatted_psi_file.index, :]

    st = time.time()
    metadata = variance_based_feature_selection(formatted_psi_file, metadata, fold_threshold=fold_threshold, samples_differing=samples_differing)
    et = time.time()
    el_t = (et-st) / 60
    print(f"Time taken for variance based feature selection: {el_t} minutes")
    metadata = metadata.loc[metadata['high_variance'] == True, :]
    formatted_psi_file = formatted_psi_file.loc[metadata.index, :]

    # inter-feature correlation module
    st = time.time()
    formatted_psi_file = intercorrelation_based_feature_selection(formatted_psi_file, corr_threshold=corr_threshold_intercorr, corr_n_events=corr_n_events)
    et = time.time()
    el_t = (et - st) / 60
    print(f"Time taken for intercorrelated based feature selection: {el_t} minutes")
    metadata = metadata.loc[formatted_psi_file.index, :]

    # impute the psi matrix for certain downstream steps
    # formatted_psi_file_imp, metadata_imp = format_psi_matrix(imputed_psi_file_path)
    formatted_psi_file_imp = median_impute(formatted_psi_file)
    formatted_psi_file_imp.index = formatted_psi_file.index

    # Find the non-redundant splicing events
    st = time.time()
    list_of_events = remove_redundant_splicing_events(formatted_psi_file, corr_threshold, metadata)
    et = time.time()
    el_t = (et - st) / 60
    print(f"Time taken for redundant splicing event removal: {el_t} minutes")
    list_of_events = list(set(list_of_events))


    # Subset the PSI matrix to only the non-redundant events
    formatted_psi_file_0 = formatted_psi_file.loc[list_of_events, :]
    formatted_psi_file_imp_0 = formatted_psi_file_imp.loc[list_of_events, :]

    # Set working directory to save evaluations
    path = os.path.join(wd, dir)
    os.mkdir(path)
    os.chdir(path)

    # write the essential files out for future assessments
    list_of_events_df = pd.DataFrame(list_of_events)
    list_of_events_df.to_csv("list_of_non_redundant_events.txt", sep="\t")

    metadata.to_csv("final_metadata.txt", sep="\t")

    # Initialize an empty dictionaries to store DEG dataframes,
    final_clusters_dict = {}
    depleted_psi_file_after_round_imp_dict = {}
    depleted_psi_file_after_round_dict = {}
    deg_results_all_dict = {}

    print("Number of PCs used for feature selection")
    print(npcs)

    # feature selection prior to Round 1
    st = time.time()
    pca_events_round = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, list_of_events, corr_threshold=pca_corr_threshold, n_components=npcs)
    et = time.time()
    el_t = (et - st) / 60
    print(f"Time taken for PCA based feature selection round 1: {el_t} minutes")

    pca_events_round_df = pd.DataFrame(pca_events_round)
    pca_events_round_df.to_csv("pca_events_round_1.txt", sep="\t")

    print("STARING ROUND 1...")
    print("Number of events prior to entering ROUND 1: ")
    print(len(formatted_psi_file_0.index))

    # Round 1 OncoSplice
    st = time.time()
    final_clusters_i, depleted_psi_file_after_round_imp_i, depleted_psi_file_after_round_i, deg_results_all_i = round_wrapper(
        filename="Round1", full_psi_file=formatted_psi_file_0, full_imputed_psi_file=formatted_psi_file_imp_0, metadata=metadata,
        highly_variable_events=pca_events_round, rank=rank1, min_group_size=min_group_size, dPSI=dPSI, dPSI_p_val=dPSI_p_val,
        min_differential_events=min_differential_events, top_n_differential_events=top_n_differential_events, conservation=conservation, strictness="tough",
        depletion_corr_threshold=depletion_corr_threshold, write_files=True, speed_corr=speed)

    et = time.time()
    el_t = (et - st) / 60
    print(f"Time taken for round 1 clustering: {el_t} minutes")

    # Round i OncoSplice (tough)
    depleted_events_round_i = depleted_psi_file_after_round_i.index.to_list()

    print("Number of events removed after ROUND 1: ")
    print(len(formatted_psi_file_0.index) - len(depleted_psi_file_after_round_i.index))

    deg_results_all_i['cluster'] = 'R1_' + deg_results_all_i['cluster']
    # Add a meaninful prefix to all column names
    prefix = "R1_C"
    final_clusters_i = final_clusters_i.rename(columns=lambda x: prefix + str(x))

    # store the mandatory round 1 outputs in respective dictionaries
    key = "Round1"
    final_clusters_dict[key] = final_clusters_i
    # depleted_psi_file_after_round_imp_dict[key] = depleted_psi_file_after_round_imp_i
    # depleted_psi_file_after_round_dict[key] = depleted_psi_file_after_round_i
    deg_results_all_dict[key] = deg_results_all_i

    for round_j in range(2, n_rounds+1):

        # feature selection prior to Round J
        st = time.time()
        pca_events_round_j = get_events_from_pca(formatted_psi_file_imp, formatted_psi_file, depleted_events_round_i,
                                                 corr_threshold=pca_corr_threshold, n_components=npcs)

        et = time.time()
        el_t = (et - st) / 60
        print(f"Time taken for PCA based feature selection round {round_j}: {el_t} minutes")

        pca_events_round_df_j = pd.DataFrame(pca_events_round_j)
        pca_events_round_df_j.to_csv(f"pca_events_round_{round_j}.txt", sep="\t")

        print(f"STARING ROUND {round_j}...")
        print(f"Number of events prior to entering ROUND {round_j}: ")
        print(len(depleted_psi_file_after_round_i.index))

        # Round j Oncosplice
        st = time.time()
        final_clusters_j, depleted_psi_file_after_round_imp_j, depleted_psi_file_after_round_j, deg_results_all_j = round_wrapper(
            filename=f"Round{round_j}", full_psi_file=depleted_psi_file_after_round_i,
            full_imputed_psi_file=depleted_psi_file_after_round_imp_i, highly_variable_events=pca_events_round_j, metadata=metadata, rank=rank,
            min_group_size=min_group_size, dPSI=dPSI, dPSI_p_val=dPSI_p_val, min_differential_events=min_differential_events,
            top_n_differential_events=top_n_differential_events, conservation=conservation, strictness="tough",
            depletion_corr_threshold=depletion_corr_threshold, write_files=True, speed_corr=speed)

        et = time.time()
        el_t = (et - st) / 60
        print(f"Time taken for round {round_j}: {el_t} minutes")

        depleted_events_round_j = depleted_psi_file_after_round_j.index.to_list()
        print(f"Number of events removed after ROUND {round_j}: ")
        print(len(depleted_psi_file_after_round_i.index) - len(depleted_psi_file_after_round_j.index))

        deg_results_all_j['cluster'] = f'R{round_j}_' + deg_results_all_j['cluster']
        # Add a meaninful prefix to all column names
        prefix = f"R{round_j}_C"
        final_clusters_j = final_clusters_j.rename(columns=lambda x: prefix + str(x))

        # store the outputs in respective dictionaries
        key = f"Round{round_j}"
        final_clusters_dict[key] = final_clusters_j
        # depleted_psi_file_after_round_imp_dict[key] = depleted_psi_file_after_round_imp_j
        # depleted_psi_file_after_round_dict[key] = depleted_psi_file_after_round_j
        deg_results_all_dict[key] = deg_results_all_j

        # reset vars for next round
        depleted_events_round_i = depleted_events_round_j
        depleted_psi_file_after_round_i = depleted_psi_file_after_round_j
        depleted_psi_file_after_round_imp_i = depleted_psi_file_after_round_imp_j

        if len(depleted_psi_file_after_round_i) < top_n_differential_events:
            warnings.warn(f"Less than 50 splicing events pending after round {round_j}. Halting further iterative clustering process...")
            break

    path = os.path.join(path, "FinalResults")
    os.mkdir(path)
    os.chdir(path)

    final_clusters_all_rounds = pd.concat(final_clusters_dict.values(), axis=1)
    final_clusters_all_rounds.to_csv("MergedResults.txt", sep="\t")

    deg_results_all_rounds = pd.concat(deg_results_all_dict.values(), axis=0, ignore_index=True)

    # Grouping by cluster assignment and type of feature, and counting occurrences
    grouped_deg_df = deg_results_all_rounds.groupby(['cluster', 'EventAnnotation']).size().unstack(fill_value=0)
    grouped_deg_df.to_csv("event_annotations_numbers.txt", sep='\t')

    # Calculating percentage of each type of feature in each cluster
    grouped_annotation_percentage = grouped_deg_df.div(grouped_deg_df.sum(axis=1), axis=0) * 100
    grouped_annotation_percentage.to_csv("event_annotations_percentage.txt", sep='\t')

    # Define color palette for the types of features
    colors = sns.color_palette("Set1")

    # Plotting horizontal stacked bar graph
    plt.figure(figsize=(10, 6))
    grouped_annotation_percentage.plot(kind='barh', stacked=True, color=colors)
    plt.title('Percentage of Features in Each Type by Cluster Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Percentage')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(path, "event_annotations_percentage.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')

    plt.figure(figsize=(10, 6))
    grouped_deg_df.plot(kind='barh', stacked=True, color=colors)
    plt.title('Number of Features in Each Type by Cluster Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Number of Features')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(path, "event_annotations_numbers.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')

    # Grouping by cluster assignment and type of feature, and counting occurrences
    grouped_deg_df = deg_results_all_rounds.groupby(['cluster', 'event_direction'])['ClusterID'].nunique().unstack(fill_value=0)
    grouped_deg_df.to_csv("event_direction_numbers.txt", sep='\t')

    # Calculating percentage of each type of feature in each cluster
    grouped_annotation_percentage = grouped_deg_df.div(grouped_deg_df.sum(axis=1), axis=0) * 100
    grouped_annotation_percentage.to_csv("event_direction_percentage.txt", sep='\t')

    # Plotting horizontal stacked bar graph
    plt.figure(figsize=(10, 6))
    grouped_annotation_percentage.plot(kind='barh', stacked=True, color=colors)
    plt.title('Percentage of Cluster IDs in Each Type by Cluster/Group Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Percentage')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(path, "event_direction_percentage.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')

    plt.figure(figsize=(10, 6))
    grouped_deg_df.plot(kind='barh', stacked=True, color=colors)
    plt.title('Number of ClusterID in Each Type by Cluster Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Number of Features')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(path, "event_direction_numbers.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Oncosplice on the PSI file')
    parser.add_argument('--wd', type=str, default=None, help='set working directory')  # dont change wd; instead have this as the path to write files to
    parser.add_argument('--dir', type=str, default=None, help='directory name to create for writing the results')  
    parser.add_argument('--psi_file', type=str, default=None, help='the path to the PSI file (with missing values)')
    parser.add_argument('--n_metadata_cols', type=int, default=11, help='expected number of metadata columns appended in front of PSI file')
    parser.add_argument('--remove_redundant_events_threshold', type=float, default=0.8, help='correlation threshold for removing redundant splicing events')
    parser.add_argument('--pca_threshold', type=float, default=0.4, help='correlation threshold for doing PCA based feature selection')
    parser.add_argument('--npcs', type=int, default=30, help='number of PCs to account for in PCA based feature selection')
    parser.add_argument('--rank', type=str, default="k30", help='rank for NMF based clustering. Set to "k30" for "overclustering" (original MV algorithm) or ignore it for automatic clustering')
    parser.add_argument('--force_broad', type=str, default="off", help='Whether to find broad 2 clusters in Round 1 of Oncosplice')
    parser.add_argument('--min_group_size', type=int, default=5, help='minimum group size required for calculating differential events')
    parser.add_argument('--dpsi', type=float, default=0.1, help='delta psi required to be considered as differential event')
    parser.add_argument('--dpsi_pval', type=float, default=0.05, help='p value threshold for differential events')
    parser.add_argument('--min_differential_events', type=int, default=100, help='minimum number of differential events required for a cluster to not be removed')
    parser.add_argument('--top_n_events', type=int, default=150, help='top n differential events; max number of marker events per clusters for LinearSVM')  # add more details
    parser.add_argument('--conservation', type=str, default="stringent", help='strategy for cluster removal. see code comments for more details.')
    parser.add_argument('--depletion_corr', type=float, default=0.4, help='depletion correlation for keeping only the unique events after a round')
    parser.add_argument('--speed', type=str, default="og", help='which algorithm to use for correlation depeletion')
    parser.add_argument('--n_rounds', type=int, default=3, help='number of rounds of iterative clustering. default: 3 sets')
    parser.add_argument('--fold_threshold', type=float, default=0.2,
                        help='fold threshold for variance based filtering. default: 0.3')
    parser.add_argument('--samples_differing', type=int, default=3,
                        help='n corresponding to n_highest - n_lowest in fold calculation in variance based filtering. default: 4')
    parser.add_argument('--corr_feature_selection', type=float, default=0.2,
                        help='intercorrelation feature selection correlation. default: 0.2')
    parser.add_argument('--corr_n_events', type=int, default=5,
                        help='minimum number of events a splicing event should be correlated with. default: 5')
    arg = parser.parse_args()
    start_time = time.time()  # Record the start time
    oncosplice(arg)
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time  # Calculate the elapsed time in seconds
    elapsed_minutes = elapsed_time / 60  # Convert elapsed time to minutes
    print(f"Code execution time: {elapsed_minutes:.2f} minutes")
