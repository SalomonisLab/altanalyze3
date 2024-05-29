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
