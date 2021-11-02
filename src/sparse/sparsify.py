#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import anndata as ad
import pandas as pd
import numpy as np
import sys,os
from scipy.sparse import csr_matrix
import h5py

# inspection
dense = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/migration/sparse/event_annotation.txt',sep='\t')
cleaned = pd.DataFrame(data=dense.iloc[:,11:].values,index=dense['UID'].values,columns=dense.columns[11:].values)

frac_zero = 1- np.count_nonzero(cleaned.values) / cleaned.values.size  # 0.17
frac_nan = np.count_nonzero(np.isnan(cleaned.values)) / cleaned.values.size  # 0.16

# each component
'''psi value'''
tmp = dense.iloc[:,11:].values
psi = csr_matrix(tmp)

'''psi column'''
psi_column = dense.columns[11:].values
psi_column_max_length = max([len(item) for item in psi_column])

'''psi index'''
psi_index = dense['UID'].values
psi_index_max_length = max([len(item) for item in psi_index])

# transformation
with h5py.File('./test.h5alt','w') as f:   
    g_psi = f.create_group('psi')
    g_psi.create_dataset('psi_data',data=psi.data.astype(np.float16))
    g_psi.create_dataset('psi_indices',data=psi.indices.astype(np.int16))
    g_psi.create_dataset('psi_indptr',data=psi.indptr.astype(np.int32))
    g_psi.create_dataset('psi_column',data=psi_column.astype(np.dtype('|S{}'.format(psi_column_max_length+1))))
    g_psi.create_dataset('psi_index',data=psi_index.astype(np.dtype('|S{}'.format(psi_index_max_length+1))))

# benchmark with anndata
adata = ad.AnnData(X=psi,var=pd.DataFrame(index=psi_column),obs=pd.DataFrame(index=psi_index))
adata.write('test.h5ad')


