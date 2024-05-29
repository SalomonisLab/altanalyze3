import pandas as pd
import numpy as np

def calculate_sample_weights_vect(full_psi_file):
    percentage_non_missing = (1 - full_psi_file.isna().mean()) * 100
    weights = percentage_non_missing.values
    return weights

def calculate_weighted_correlation_vect(nmf_basis_matrix, full_psi_file, weights):
    nmf_basis_matrix = pd.DataFrame(nmf_basis_matrix).transpose()
    nmf_basis_matrix.columns = full_psi_file.columns

    weights_matrix = np.outer(weights, weights)

    centered_nmf = nmf_basis_matrix - np.average(nmf_basis_matrix, axis=0, weights=weights)
    centered_full_psi = full_psi_file - np.average(full_psi_file, axis=0, weights=weights)

    cov_matrix = np.dot(centered_nmf * weights[:, np.newaxis], centered_full_psi.T * weights[np.newaxis, :])

    std_nmf = np.sqrt(np.dot(centered_nmf ** 2, weights))
    std_full_psi = np.sqrt(np.dot(centered_full_psi ** 2, weights))

    corr_matrix = cov_matrix / np.outer(std_nmf, std_full_psi)

    return corr_matrix

def deplete_events_vectorized(nmf_basis_matrix, full_psi_file, corr_threshold=0.3, strictness="tough"):
    weights = calculate_sample_weights_vect(full_psi_file)
    corr_matrix = calculate_weighted_correlation_vect(nmf_basis_matrix, full_psi_file, weights)

    mask = np.abs(corr_matrix) < corr_threshold

    if strictness == "tough":
        has_low_corr = np.all(mask, axis=0)
    else:
        has_low_corr = np.any(mask, axis=0)

    depleted_psi_file = full_psi_file.loc[has_low_corr]

    return depleted_psi_file
