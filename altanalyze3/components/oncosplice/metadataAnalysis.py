import pandas as pd
import numpy as np
import scipy.stats
import statsmodels.stats.multitest as ssm
import os


# ### Define function to perform some data checks before mwu computation

def mwuChecks(dataf, grpdf, grpvar='grp', min_group_size=6):
    #
    #   This function performs some checks before calling mwuCompute().
    #
    #   dataf is a pandas data frame with samples as columns and events (features) as rows.
    #   Event names must be the row names (index), not in a column.
    #   Missing values (np.nan='NaN') are allowed
    #
    #   grpdf is a pandas data frame with a group membership column ('grp' is the default ).
    #   Sample names must be the row names (index), not in a column.
    #   Aside from ordering, must be 1-1 match between dataf.columns and grpdf.index
    #
    #   min_group_size is minimum group size for kw test to be computed.
    #   min_group_size = 6 is based on statistical rule-of-thumb for Mann-Whitney 2-group test.
    #
    #  Returns integer to indicate course of action:
    #     0 = Problem, do not run mwuCompute
    #     1 = Run mwuCompute

    #    print(grpvar)
    #    print(type(grpvar))

    # --- Check for 1-1 match between sample names provided in grpdf and the data matrix

    if set(dataf.columns) != set(grpdf.index):
        print('Group information is inconsistent with names in data matrix')
        return 0

    grpCounts = grpdf[grpvar].value_counts().to_frame()
    # returns a pandas df of counts, indexed by grp, decreasing by count
    grpCounts.columns = ['grpRawN']
    nGrps = grpCounts.shape[0]
    minGrpN = min(grpCounts['grpRawN'])
    # print('nGrps, minGrpN',nGrps,minGrpN)

    # -- Handle groups <> 2 --

    if nGrps < 2:
        print('Number of sample groups is < 2; Mann-Whitney test not conducted')
        return 0

    if nGrps > 2:
        print('Number of sample groups is > 2; Mann-Whitney test not conducted')
        return 2

    # -- Don't proceed if already know a group size is < minimum --

    if minGrpN < min_group_size:
        print('Mann-Whitney test not conducted: A group has fewer samples than minimum group size: ', minGrpN, ' < ',
              min_group_size)
        return 0

    return 1


# ### Define function to compute the mwu test and collect summary stats (Ns, medians)

def mwuCompute(dataf, grpdf, metadata_df, grpvar='grp', min_group_size=6):
    #
    #   This function performs the Mann-Whitney NP (rank) two independent groups test, aka Wilcoxon 2 grp indep test.
    #
    #   dataf is a pandas data frame with samples as columns and events (features) as rows.
    #   Event names must be the row names (index), not in a column.
    #   Missing values (np.nan='NaN') are allowed
    #
    #   grpdf is a pandas data frame with a group membership column.
    #   Sample names must be the row names (index), not in a column.
    #   Aside from ordering, must be 1-1 match between dataf.columns and grpdf.index
    #
    #   min_group_size is minimum group size for mwu test to be computed.
    #   6 is based on statistical rule-of-thumb for Mann-Whitney 2-group test.
    #
    #   Direction of comparison is by sorted grp levels, e.g. A vs. B, 1 vs 2.
    #   Sign of result is determined by comparing test statistic value to the null mean (not comparing group medians or means)
    #
    #   Returns pandas dataframe with columns containing group Ns (xNaN), group medians and means, value of mwu test statistic, sign and p-value.
    #   Returned df has row for each event, even if test was not computed (in this case mwu values will be NaN)

    groups = dataf.T.groupby(grpdf[grpvar])
    #    print(groups.groups)

    grpCounts = grpdf[grpvar].value_counts().to_frame()
    #    display(grpCounts)

    grpLevels = grpCounts.index
    grpLevels = grpLevels.tolist()
    #    print(grpLevels)
    SgrpLevels = sorted(grpLevels)
    #    print("SgrpLevels",SgrpLevels)

    g1Level = SgrpLevels[0]
    g2Level = SgrpLevels[1]
    #    print('g1Level',g1Level)
    #    print('g2Level',g2Level)

    g1 = groups.get_group(g1Level)
    g2 = groups.get_group(g2Level)
    #    print('g1 dat',g1)
    #    print('g2 dat',g2)
    #   count() method counts number of non-missing entries
    eventsToRun = ((g1.count() >= min_group_size) & (g2.count() >= min_group_size)).values
    #    print('eventsToRun, sum',eventsToRun,sum(eventsToRun))
    #    print(g1.iloc[:,eventsToRun])
    #    print(g2.iloc[:,eventsToRun])

    #   Set up df for results
    # tempdf.columns=[('N_' + str(gID)),('Median_' + str(gID))]

    result = pd.DataFrame(
        {
            ('N_' + str(g1Level)): g1.count(),
            ('N_' + str(g2Level)): g2.count(),
            ('Median_' + str(g1Level)): g1.median(),
            ('Median_' + str(g2Level)): g2.median(),
            ('Mean_' + str(g1Level)): g1.mean(),
            ('Mean_' + str(g2Level)): g2.mean(),
            'mwuStat': np.nan,
            'mwuPval': np.nan,
            'mwuAdjPval': np.nan,
            'mwuSign': np.nan,
            'Abs_Mean_Diff': np.nan,
        },
        index=dataf.index)

    #    display(result)

    #   Handle situation of no events meeting min_group_size criterion

    if sum(eventsToRun) == 0:
        print('Mann-Whitney test not conducted: No events meet minimum group size criterion of ', min_group_size, '.')
        return result

    #   Compute the mwu test statistic & p-value for eventsToRun

    mwu = scipy.stats.mannwhitneyu(
        g1.iloc[:, eventsToRun].values,
        g2.iloc[:, eventsToRun].values,
        alternative="two-sided",
        method="asymptotic",
        use_continuity=True,
        nan_policy='omit'
    )
    #    print(mwu)

    #   Populate stat & p-value columns
    result.loc[eventsToRun, "mwuStat"] = mwu.statistic
    result.loc[eventsToRun, "mwuPval"] = mwu.pvalue
    fdr_bh_test = ssm.multipletests(result.loc[eventsToRun, "mwuPval"].to_numpy(), method="fdr_bh")
    result.loc[eventsToRun, "mwuAdjPval"] = fdr_bh_test[1]

    #    display(result)

    #   Determine sign (direction) of test result based on expected mean of test stat under null
    nullMean = (result[('N_' + str(g1Level))] * result[('N_' + str(g2Level))]) / 2.0 + 0.5
    #    print(nullMean)
    result["mwuSign"] = np.select(
        [
            (result["mwuStat"] < nullMean),
            (result["mwuStat"] == nullMean),
            (result["mwuStat"] > nullMean),
        ],
        [-1, 0, 1],
        default=np.nan)

    #    display(result)

    result['cluster'] = 'c' + grpvar
    metadata_df = metadata_df.loc[result.index, ]
    result = pd.concat([metadata_df, result], axis=1)

    return result


def find_top_differential_events(exp_file, groups_file, min_group_size, metadata, output_loc, dPSI=0.1, dPSI_p_val=0.05, min_differential_events=100, top_n=150):

    all_reg_events = []
    remove_clusters = []

    groups_file_cols = groups_file.columns
    deg_cols = pd.Index(['N_cluster',
                         'N_others',
                         'Median_cluster',
                         'Median_others',
                         'Mean_cluster',
                         'Mean_others',
                         'mwuStat',
                         'mwuPval',
                         'mwuAdjPval',
                         'mwuSign',
                         'Abs_Mean_Diff',
                         'cluster'])
    metadata_cols = metadata.columns
    results_cols = metadata_cols.union(deg_cols, sort=False)

    results_all = pd.DataFrame(columns=results_cols)

    for cluster in range(len(groups_file_cols)):
        results_c0 = mwuCompute(exp_file, groups_file, grpvar=groups_file_cols[cluster], metadata_df=metadata, min_group_size=min_group_size)
        # results_c0 = results_c0[results_c0['mwuPval'] < pval_threshold]  # for using raw p value to determining significant events
        results_c0.loc[:, "Abs_Mean_Diff"] = abs(results_c0['Mean_0'] - results_c0['Mean_1'])
        results_c0.columns = results_cols
        results_c0 = calculate_coordinate_lengths(results_c0)
        results_c0 = determine_event_direction(results_c0)
        # filter based on deg requirements
        results_c0 = results_c0[results_c0['mwuAdjPval'] < dPSI_p_val]
        results_c0 = results_c0[results_c0["Abs_Mean_Diff"] > dPSI]
        if results_c0.shape[0] < min_differential_events:
            print(("Less than " + str(min_differential_events) + " significant diff events for Cluster " + groups_file_cols[cluster]))
            remove_clusters.append(cluster)
        else:
            # results_c0 = results_c0.sort_values(by=['Abs_Mean_Diff'], ascending=False)  # for using top n events sorted by fold change instead of p val
            results_c0 = results_c0.sort_values(by=['mwuAdjPval'])  # does not matter much whether you sort by raw or adjusted p value here -- results will be identical most likely
            results_all = pd.concat([results_all, results_c0], axis=0)  # equivalent to rbind
            file_path = os.path.join(output_loc, 'c_' + groups_file_cols[cluster] + "_vs_others_splicing_events.txt")
            results_c0.to_csv(file_path, sep="\t")
            results_c0 = results_c0.iloc[:top_n, :]
            all_reg_events = all_reg_events + results_c0.index.to_list()
            # top150_reg_events_file = "/Users/tha8tf/Documents/testsOncosplice/differential_events_per_cluster/all_reg_events_force_broad_nmf_cluster_adjp_0.05_" + cluster + ".txt"
            # results_c0.to_csv(top150_reg_events_file, sep='\t')

    all_reg_events_unique = pd.Series(all_reg_events).drop_duplicates()
    # print(("Number of unique differential events found for this round is " + str(len(all_reg_events_unique.tolist()))))

    remove_clusters = list(map(int, remove_clusters))  # converting string to integer for the clusters to remove so I can use these integers for filtering the NMF binary groups file

    return all_reg_events_unique, remove_clusters, results_all


def calculate_coordinate_lengths(deg_result):

    # Function to split the string by "|" or ":"
    def split_string(entry, sep):
        return entry.split(sep)

    # Apply the function to split the coordinates into numerator and denominator junctions (| is the separator)
    split_data = deg_result['Coordinates'].apply(split_string, sep="|")

    # Assign the split values to new columns
    deg_result['num_coordinates'] = split_data.apply(lambda x: x[0])
    deg_result['den_coordinates'] = split_data.apply(lambda x: x[1])

    num_coordinates = pd.DataFrame(deg_result['num_coordinates'].apply(split_string, sep=":"))
    num_coordinates['num_coordinates'] = num_coordinates['num_coordinates'].apply(lambda x: x[1])
    num_coordinates['num_coordinates'] = num_coordinates['num_coordinates'].apply(split_string, sep="-")
    num_coordinates['Coordinate1'] = pd.to_numeric(num_coordinates['num_coordinates'].apply(lambda x: x[0]))
    num_coordinates['Coordinate2'] = pd.to_numeric(num_coordinates['num_coordinates'].apply(lambda x: x[1]))
    num_coordinates['num_coordinate_diff'] = abs(num_coordinates['Coordinate1'] - num_coordinates['Coordinate2'])

    den_coordinates = pd.DataFrame(deg_result['den_coordinates'].apply(split_string, sep=":"))
    den_coordinates['den_coordinates'] = den_coordinates['den_coordinates'].apply(lambda x: x[1])
    den_coordinates['den_coordinates'] = den_coordinates['den_coordinates'].apply(split_string, sep="-")
    den_coordinates['Coordinate1'] = pd.to_numeric(den_coordinates['den_coordinates'].apply(lambda x: x[0]))
    den_coordinates['Coordinate2'] = pd.to_numeric(den_coordinates['den_coordinates'].apply(lambda x: x[1]))
    den_coordinates['den_coordinate_diff'] = abs(den_coordinates['Coordinate1'] - den_coordinates['Coordinate2'])

    deg_result['num_coordinate_diff'] = num_coordinates['num_coordinate_diff']
    deg_result['den_coordinate_diff'] = den_coordinates['den_coordinate_diff']

    return deg_result


def determine_event_direction(deg_result):

    # Junction1 | Junction2
    # Example: chr3: 129575563 - 129575562 | chr3:129575766 - 129575562
    # A = abs(129575563 - 129575562) = 1
    # B = abs(129575766 - 129575562) = 204
    # A = num_coordinate_diff
    # B = den_coordinate_diff
    # If A < B and dPSI > 0: type = inclusion
    # If A > B and dPSI < 0: type = inclusion
    # If A < B and dPSI < 0: type = exclusion
    # If A > B and dPSI > 0: type = exclusion

    deg_result['event_direction'] = 'inclusion'

    mask_junction = deg_result['num_coordinate_diff'] < deg_result['den_coordinate_diff']
    mask_regulation = deg_result['mwuSign'] < 0

    deg_result.loc[mask_junction & mask_regulation, 'event_direction'] = 'exclusion'
    deg_result.loc[~mask_junction & ~mask_regulation, 'event_direction'] = 'exclusion'

    deg_result.drop(columns=['num_coordinates', 'den_coordinates'], inplace=True)

    return deg_result

