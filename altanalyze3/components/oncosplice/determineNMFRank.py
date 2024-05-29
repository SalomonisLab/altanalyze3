#  IDKfunction2.txt
#
#  IDK function
#
#  August 16, 2022
#
#  Can identify Intrinsic Dimension K of numeric matrix using 1 or more of the following methods:
#   MADA, CORRINT, MOM, PCAF0, TLE

import skdim
import numpy as np
import pandas as pd


# **** IDK function *******
#
#
def determine_nmf_ranks(PSIdf, *args, MADA=False, CORRINT=False, MOM=False, PCAFO=False, TLE=False, FS=False):
    #
    #   PSIdf is a pandas data frame with events (features) as columns, no missing values.
    #   event names are not needed; can be row names but not a column.
    #   Kmethods is a pandas series with True/False for each available method.
    #   Available methods are: MADA, CORRINT, MOM, PCAF0, TLE, FS.
    #
    # MADA is Manifold-Adaptive Dimension Estimation
    # CORRINT is Correlation Dimension
    # MOM is Method of Moments
    # PCAFO is the Fukunaga-Olsen version of the linear PCA method
    # TLE is Tight Local intrinsic dimensionality Estimator
    # FS is the FisherS method (PCA w/ l-descrim)

    #   Set up pandas series for estimates of K.  NaN for not done or error.
    Kests = {'MADA': np.nan, 'CORRINT': np.nan, 'MOM': np.nan, 'PCAFO': np.nan, 'TLE': np.nan, 'FS': np.nan}
    #
    Kests = pd.Series(Kests)
    #
    if MADA:
        mada = skdim.id.MADA().fit(PSIdf)
        Kests['MADA'] = mada.dimension_

    if CORRINT:
        corrint = skdim.id.CorrInt().fit(PSIdf)
        Kests['CORRINT'] = corrint.dimension_
    #
    if MOM:
        Mom = skdim.id.MOM().fit(PSIdf)
        Kests['MOM'] = Mom.dimension_
    #
    if PCAFO:
        pcafo = skdim.id.lPCA(ver='FO').fit(PSIdf)
        Kests['PCAFO'] = pcafo.dimension_
    #
    if TLE:
        tle = skdim.id.TLE(epsilon=0.0001).fit(PSIdf)
        Kests['TLE'] = tle.dimension_

    if FS:
        Fs = skdim.id.FisherS().fit(PSIdf)
        Kests['FS'] = Fs.dimension_
    #
    return Kests