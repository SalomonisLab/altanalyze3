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


def limmaCompute(dataf, grpdf, *args, grpvar='grp', min_group_size=2, shrink_factor=0.2):
    """limma-like moderated two-group test -- a drop-in alternative to mwuCompute.

    Same I/O contract as mwuCompute (so downstream annotation/filtering is unchanged):
      dataf : features (index) x samples (columns) matrix; NaN allowed.
      grpdf : sample (index) -> group (column ``grpvar``); exactly two groups.
    Returns the SAME column schema as mwuCompute -- group Ns, medians, means,
    DiffMeans, and ``mwuStat`` / ``mwuSign`` / ``mwuPval`` / ``mwuAdjPval`` -- but the statistic and
    p-value come from the empirical-Bayes MODERATED t-test used in cellHarmony-differential for
    pseudobulks (``cellHarmony_differential._moderated_t_test``): pooled variance shrunk toward its
    median (eBayes), moderated t, t-distribution p-values, BH-FDR over tested events. ``mwuStat`` then
    holds the moderated t-statistic; ``mwuSign`` = sign of (group1 mean - group2 mean).

    This makes the parametric pseudobulk test selectable wherever the MWU rank test was used, with no
    change to the consumers of the stats file.
    """
    import pandas as pd
    import numpy as np
    import scipy.stats as scs
    import statsmodels.stats.multitest as smm

    groups = dataf.T.groupby(grpdf[grpvar])
    grpLevels = sorted(grpdf[grpvar].value_counts().index.tolist())
    if len(grpLevels) < 2:
        raise ValueError(f"limmaCompute needs exactly two groups; got {grpLevels}")
    g1Level, g2Level = grpLevels[0], grpLevels[1]
    g1 = groups.get_group(g1Level)   # samples x features
    g2 = groups.get_group(g2Level)

    diffGrpMeans = g1.mean() - g2.mean()
    result = pd.DataFrame(
        {
            ('N_' + str(g1Level)): g1.count(),
            ('N_' + str(g2Level)): g2.count(),
            ('Median_' + str(g1Level)): g1.median(),
            ('Median_' + str(g2Level)): g2.median(),
            ('Mean_' + str(g1Level)): g1.mean(),
            ('Mean_' + str(g2Level)): g2.mean(),
            ('DiffMeans_' + str(g1Level) + '_m_' + str(g2Level)): diffGrpMeans,
            'mwuStat': np.nan, 'mwuSign': np.nan, 'mwuPval': np.nan, 'mwuAdjPval': np.nan,
        },
        index=dataf.index)

    # events with enough non-missing samples in BOTH groups
    eventsToRun = ((g1.count() >= min_group_size) & (g2.count() >= min_group_size)).values
    if eventsToRun.sum() == 0:
        print(f"limma moderated t-test not conducted: no events meet min group size {min_group_size}.")
        return result

    X1 = g1.to_numpy(dtype=float)   # samples1 x features
    X2 = g2.to_numpy(dtype=float)   # samples2 x features
    n1 = np.sum(~np.isnan(X1), axis=0).astype(float)
    n2 = np.sum(~np.isnan(X2), axis=0).astype(float)
    mean1 = np.nanmean(X1, axis=0)
    mean2 = np.nanmean(X2, axis=0)
    # ddof=1 variance, NaN-aware
    var1 = np.nanvar(X1, axis=0, ddof=1)
    var2 = np.nanvar(X2, axis=0, ddof=1)
    s2 = (var1 + var2) / 2.0

    # empirical Bayes shrinkage of the variance toward its median (limma-like)
    s2_prior = np.nanmedian(s2[eventsToRun])
    s2_shrunk = (1 - shrink_factor) * s2 + shrink_factor * s2_prior
    with np.errstate(invalid='ignore', divide='ignore'):
        se = np.sqrt(s2_shrunk / np.maximum(n1, 1) + s2_shrunk / np.maximum(n2, 1))
    se = np.maximum(se, 1e-12)
    diff = mean1 - mean2
    with np.errstate(invalid='ignore', divide='ignore'):
        tvals = diff / se
    df = np.maximum(n1 + n2 - 2, 1)
    pvals = 2 * scs.t.sf(np.abs(tvals), df=df)

    run = eventsToRun & np.isfinite(pvals)
    result.loc[run, 'mwuStat'] = tvals[run]            # moderated t-statistic
    result.loc[run, 'mwuPval'] = pvals[run]
    result['mwuSign'] = np.select([diff < 0, diff == 0, diff > 0], [-1, 0, 1], default=np.nan)

    # BH-FDR over the tested events only
    eventsForFDR = run & (~np.isnan(result['mwuPval'].to_numpy()))
    if eventsForFDR.sum() > 0:
        fdr = smm.fdrcorrection(result.loc[eventsForFDR, 'mwuPval'].to_numpy(),
                                method='indep', is_sorted=False)[1]
        result.loc[eventsForFDR, 'mwuAdjPval'] = fdr
    return result

#  End of file