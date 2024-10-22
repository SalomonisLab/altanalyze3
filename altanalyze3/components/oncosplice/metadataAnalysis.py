#!/usr/bin/env python
# coding: utf-8

# ##### mwu6.py
# 
# May 3, 2023
# 
# Mann-Whitney U non-parametric (rank-based) test (aka Wilcoxon for 2 independent groups)
# 


# ### Define function to perform some data checks before mwu computation

def mwuChecks( dataf, grpdf, *args, grpvar='grp',min_group_size=6 ): 
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


    import pandas as pd
    import numpy as np

    
#    print(grpvar)
#    print(type(grpvar))
    
    # --- Check for 1-1 match between sample names provided in grpdf and the data matrix

    if set(dataf.columns)!= set(grpdf.index) :
        print('Group information is inconsistent with names in data matrix')
        return( 0 )
  
        
    grpCounts = grpdf[grpvar].value_counts().to_frame()  
# returns a pandas df of counts, indexed by grp, decreasing by count
    grpCounts.columns=['grpRawN']
    nGrps = grpCounts.shape[0]
    minGrpN = min( grpCounts['grpRawN'] )
    # print('nGrps, minGrpN',nGrps,minGrpN)
    
    
    # -- Handle groups <> 2 --

    if nGrps < 2 :
        print('Number of sample groups is < 2; Mann-Whitney test not conducted')
        return( 0 )
        
    if nGrps > 2 :
        print('Number of sample groups is > 2; Mann-Whitney test not conducted')
        return( 2 )
    
    # -- Don't proceed if already know a group size is < minimum --   

    if minGrpN < min_group_size:
        print('Mann-Whitney test not conducted: A group has fewer samples than minimum group size: ',minGrpN,' < ',min_group_size)
        return( 0 )
    
    return( 1 )



# ### Define function to compute the mwu test and collect summary stats (Ns, medians)

def mwuCompute(dataf, grpdf, *args, grpvar='grp' ,min_group_size=6 ): 
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
    #   Returns pandas dataframe with columns containing group Ns (xNaN), group medians and means, difference of group means ,value of mwu test statistic, sign based on test stat, p-value and FDR-adjusted p-value.
    #   Returned df has row for each event, even if test was not computed (in this case mwu values will be NaN)
    
    import pandas as pd
    import numpy as np

    import scipy.stats

    import statsmodels.stats.multitest as smm
    
   
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

    diffGrpMeans = g1.mean() - g2.mean()

#   Set up df for results

    result = pd.DataFrame(
    {
        ('N_' + str(g1Level)): g1.count(),
        ('N_' + str(g2Level)): g2.count(),
        ('Median_' + str(g1Level)): g1.median(),
        ('Median_' + str(g2Level)): g2.median(),
        ('Mean_' + str(g1Level)): g1.mean(),
        ('Mean_' + str(g2Level)): g2.mean(),
        ('DiffMeans_' + str(g1Level) + '_m_' + str(g2Level)): diffGrpMeans,
        'mwuStat': np.nan,
        'mwuSign': np.nan,
        'mwuPval': np.nan,
        'mwuAdjPval':np.nan
    },
        index=dataf.index )
    
#    display(result)

#   Handle situation of no events meeting min_group_size criterion

    if sum(eventsToRun) == 0:
        print('Mann-Whitney test not conducted: No events meet minimum group size criterion of ',min_group_size,'.')
        return( result)

#   Compute the mwu test statistic & p-value for eventsToRun
    
    mwu = scipy.stats.mannwhitneyu(
    g1.iloc[:,eventsToRun].values, 
    g2.iloc[:,eventsToRun].values,
    alternative="two-sided",
    method="asymptotic",
    use_continuity=True,
    nan_policy='omit'
    )  
#    print(mwu)

#   Populate stat & p-value columns    
    result.loc[eventsToRun, "mwuStat"] = mwu.statistic
    result.loc[eventsToRun, "mwuPval"] = mwu.pvalue
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
    default=np.nan )

#   Calculate FDR of mwuPval
#   Re-create mask in case any of the calculated p-values were NaN
    eventsForFDR = (eventsToRun & (~ result['mwuPval'].isnull().values ) )

#   Do not provide alpha cutoff value to function since want adjusted p-values, not significance calls
    FDRres = smm.fdrcorrection(result.loc[eventsForFDR, "mwuPval"].to_numpy(), 
                               method="indep",is_sorted=False)

    result.loc[eventsForFDR,'mwuAdjPval'] = FDRres[1]

#    display(result)

    return( result)

#  End of file