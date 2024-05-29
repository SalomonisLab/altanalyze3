
def median_impute(df):
    """
    Impute missing values in a DataFrame with the median value of the respective row using a vectorized approach.

    Parameters:
    df (pandas.DataFrame): The input DataFrame with missing values.

    Returns:
    pandas.DataFrame: DataFrame with missing values imputed using row-wise medians.
    """

    # Define a custom function to impute missing values in each row with the row's median
    def impute_row(row):
        return row.fillna(row.median())

    # Apply the custom function to each row in the DataFrame
    imputed_df = df.apply(impute_row, axis=1)

    return imputed_df

# Example usage:
# Assuming you have a DataFrame called 'your_dataframe' with missing values
# imputed_df = medianImpute(your_dataframe)
