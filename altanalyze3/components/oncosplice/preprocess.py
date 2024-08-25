import pandas as pd
import numpy as np

# Variance based feature selection
def variance_based_feature_selection(formatted_psi_file, metadata, fold_threshold=0.2, samples_differing=3):

    metadata['high_variance'] = True

    for feature in formatted_psi_file.index:
        feature_vector = formatted_psi_file.loc[feature, :]
        feature_vector = feature_vector.dropna()
        feature_vector = feature_vector.sort_values(ascending=True)

        n_lowest = feature_vector.head(samples_differing).iloc[samples_differing-1]
        n_highest = feature_vector.tail(samples_differing).iloc[0]
        fold = n_highest - n_lowest
        metadata.loc[feature, 'high_variance'] = fold > fold_threshold

    return metadata


# Inter correlation step do this step if only less than 8000 to start with -- cut off for this 0.4 -- should bring it down to 1500
def intercorrelation_based_feature_selection(df, corr_threshold=0.5, corr_n_events=10):

    # Calculate correlation matrix
    df = df.transpose()

    corr_matrix = df.corr(method="pearson")

    # Get columns with correlation > 0.4
    correlated_features = (corr_matrix > corr_threshold).sum() > corr_n_events

    # Filter dataframe based on correlated features
    filtered_df = df[correlated_features.index[correlated_features]]

    filtered_df = filtered_df.transpose()

    return filtered_df


def fast_intercorrelation_based_feature_selection(df, corr_threshold=0.5, corr_n_events=10):
    # Calculate correlation matrix
    df = df.transpose()

    corr_matrix = fast_corr_matrix(data_matrix=df.to_numpy())

    corr_matrix = pd.DataFrame(corr_matrix, index=df.columns, columns=df.columns)

    # Get columns with correlation > 0.4
    correlated_features = (corr_matrix > corr_threshold).sum() > corr_n_events

    # Filter dataframe based on correlated features
    filtered_df = df[correlated_features.index[correlated_features]]

    filtered_df = filtered_df.transpose()

    return filtered_df


def test_func(df1, df2):
    norm1 = df1 - df1.mean(axis=0)
    norm2 = df2 - df2.mean(axis=0)
    sqsum1 = (norm1**2).sum(axis=0)
    sqsum2 = (norm2**2).sum(axis=0)
    return((norm1.T @ norm2) / np.sqrt(sqsum1.apply(lambda x: x*sqsum2)))


import numpy as np


def fast_corr_matrix(data_matrix):

    # Mask out NaNs by replacing them with the column mean, if necessary
    valid_data = np.nan_to_num(data_matrix, nan=np.nanmean(data_matrix, axis=0))

    # Subtract the mean from each column (centering the data)
    centered_data = valid_data - np.nanmean(data_matrix, axis=0)

    # Calculate the dot product of the centered data matrix with itself
    dot_product = np.dot(centered_data.T, centered_data)

    # Calculate the sum of squared deviations (like variance)
    sum_squared = np.sqrt(np.sum(centered_data ** 2, axis=0))

    # Outer product of sum_squared to get all combinations of feature variances
    denom = np.outer(sum_squared, sum_squared)

    # Calculate the correlation matrix
    corr_matrix = dot_product / denom

    # Set the diagonal to 1 (perfect correlation with self)
    np.fill_diagonal(corr_matrix, 1)

    return corr_matrix

