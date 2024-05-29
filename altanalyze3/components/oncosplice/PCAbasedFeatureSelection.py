
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from itertools import chain


def get_events_from_pca(formatted_imp_psi_file, formatted_psi_file, filtered_events, corr_threshold=0.4, n_components=30):

    # Do PCA on the list of events provided by user
    pc_psi_file = formatted_imp_psi_file.filter(items=filtered_events, axis=0)

    pca_obj = PCA(n_components, svd_solver='full')

    pca_obj.fit_transform(pc_psi_file.transpose())
    loadings = pd.DataFrame(pca_obj.components_.T, index=pc_psi_file.index)
    correlated_events_all_pc = []

    for pc in range(1, n_components + 1):
        #print(pc)
        pc_i = pd.Series(loadings.loc[:, pc - 1])
        pc_i_sorted = pc_i.sort_values(ascending=False)

        pc_top200 = pc_i_sorted.head(200).index.tolist()  # top 200 events
        pc_bottom200 = pc_i_sorted.tail(200).index.tolist()  # bottom 200 events

        pc_400events = pc_top200 + pc_bottom200

        psi_400_events_pc = formatted_psi_file.loc[pc_400events, :]

        cor_psi_400events_matrix = psi_400_events_pc.T.corr(method="pearson")
        cor_psi_400events_matrix = np.tril(cor_psi_400events_matrix, k=-1)
        idx_cor_psi_400events_binary = abs(cor_psi_400events_matrix) > corr_threshold
        idx_cor_psi_400events_binary = idx_cor_psi_400events_binary.astype(int)
        idx_cor_psi_400events_binary = pd.DataFrame(idx_cor_psi_400events_binary)
        idx_cor_psi_400events_binary.index = psi_400_events_pc.index.tolist()
        idx_cor_psi_400events_binary.columns = psi_400_events_pc.index.tolist()

        z = idx_cor_psi_400events_binary.sum(axis=0)
        z = z[z > 0].index.tolist()
        z2 = idx_cor_psi_400events_binary.sum(axis=1)
        z2 = z2[z2 > 0].index.tolist()

        correlated_events_pc = list(set(z) & set(z2))
        correlated_events_all_pc.append(correlated_events_pc)

    correlated_events_all_pc_unique = list(chain.from_iterable(correlated_events_all_pc))
    correlated_events_all_pc_unique = list(set(correlated_events_all_pc_unique))

    return correlated_events_all_pc_unique



# correvents_allpc_30pc = []
    #
    # for pc in range(1, ncomps):
    #     pci = principal_df.iloc[:, pc]
    #     pcsorted = pci.sort_values(ascending=False)
    #
    #     pc_top200 = pcsorted.head(200)
    #     pc_bottom200 = pcsorted.tail(200)
    #
    #     # pc_400events = pd.concat(pc_top200,pc_bottom200,axis=0)
    #     pc_400events = pd.DataFrame(pc_top200.append(pc_bottom200, ignore_index=False))
    #     psi_400_events_pc = formatted_psi_file.filter(items=pc_400events.index, axis=0)
    #     cor_psi_400events_matrix = psi_400_events_pc.transpose().corr(method="pearson")
    #     cor_psi_400events_matrix = cor_psi_400events_matrix.to_numpy()
    #     cor_psi_400events_matrix[np.triu_indices(cor_psi_400events_matrix.shape[0], 0)] = np.nan
    #     idx_cor_psi_400events_binary = abs(cor_psi_400events_matrix) > corr_threshold
    #     idx_cor_psi_400events_binary = 1 * idx_cor_psi_400events_binary
    #
    #     z = np.where(np.nansum(idx_cor_psi_400events_binary, axis=0) > 0)[0]
    #     z2 = np.where(np.nansum(idx_cor_psi_400events_binary, axis=1) > 0)[0]
    #
    #     correvents_pc_idx = np.intersect1d(z, z2)
    #     correvents_pc = psi_400_events_pc.index[correvents_pc_idx]
    #     correvents_allpc_30pc.append(correvents_pc)

    # correvents_allpc_30pc_unique = np.unique(correvents_allpc_30pc)

