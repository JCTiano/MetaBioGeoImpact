import numpy as np
import pandas as pd


def to_number(df, columns):
    df = df[pd.to_numeric(df[columns], errors='coerce').notnull()]
    df[columns] = df[columns].astype('float64')
    return df


def _weighted_sd(df, data_cols='Response value', weight_col='replicates'):
    if df[data_cols].count() == 0:
        return np.nan
    numer = np.sum(df[weight_col] * (df[data_cols] - df[data_cols].mean()) ** 2)
    denom = ((df[data_cols].count() - 1) / df[data_cols].count()) * np.sum(df[weight_col])
    if denom == 0:
        return np.nan
    else:
        return np.sqrt(numer / denom)


def _weighted_avg(df, data_cols, weight_col='replicates'):
    if isinstance(data_cols, list) and len(data_cols) > 1:
        return pd.Series(np.average(df[data_cols], weights=df[weight_col], axis=0), data_cols)
    else:
        return np.average(df[data_cols], weights=df[weight_col], axis=0)


def weighted_avg(df, data_cols, by_cols):
    if 'Sampling date' in data_cols:
        df['Sampling date'] = pd.to_datetime(df['Sampling date'], errors='coerce').values.astype(np.int64)
    for data_col in data_cols:
        if df[data_col].dtype != 'float64' or df[data_col].dtype != 'int64':
            df[data_col] = pd.to_numeric(df.loc[:, data_col], errors='coerce')
    g = df.groupby(by_cols)
    # Create a new dataframe to store the weighted average, SD, and the new weight
    result = g.apply(_weighted_avg, data_cols)
    result.reset_index(inplace=True)
    if 'Sampling date' in result.columns:
        result['Sampling date'] = pd.to_datetime(result['Sampling date'])
    return result


def weighted_values(df, data_col, weight_col, by_cols):
    result = pd.DataFrame()
    g = df.groupby(by_cols)
    if not g.size().empty:
        # Create a new dataframe to store the weighted average, SD, and the new weight
        result[f'{weight_col}'] = g[weight_col].sum()
        result[f'{data_col}'] = g.apply(_weighted_avg, data_col)
        result[f'Error estimate'] = g.apply(_weighted_sd)
        result[f'Error estimate'].fillna(g['Error estimate'].first(), inplace=True)
        result['Error estimate type'] = np.where(result['Error estimate'] == g['Error estimate'].first(),
                                                 g['Error estimate type'].first(),
                                                 'SD')
        result.reset_index(inplace=True)
    return result


def group_studies(df, variables, numerical_covariates, categorical_covariates, data_col, weight_col, by_cols):
    """
    Calculate weighed mean and standard deviation based on grouping_id, iterating through each variable.
    Iterating through each variable allows to eliminate a variable that does not provide enough information.
    :param df: DataFrame where the data is stored
    :param variables: Reponse variables to iterate through
    :param numerical_covariates: Additional variables (metadata) that needs to be averaged while grouping, such as
    latitude and longitude, water depth, and the numerical co-variates that will be included in the meta-analysis
    :param categorical_covariates: Categorical covariates that need to be appended to the output.
    This is done separate from the grouping because if they have empty
    :param data_col:
    :param weight_col:
    :param by_cols:
    :return:
    """
    df_out = pd.DataFrame()
    variables_excluded = []
    for variable in variables:
        print(f'Grouping studies based on {variable}')
        df_temp = df[df['Response variable'] == variable]
        df_temp_weighted_values = weighted_values(df_temp, data_col=data_col, weight_col=weight_col,
                                                  by_cols=by_cols)
        if not df_temp_weighted_values.empty:
            df_temp_weighted_values['Response variable'] = variable
            df_out = df_out.append(df_temp_weighted_values, ignore_index=True)
        else:
            variables_excluded.append(variable)
    print(f'Could not find enough data for variables {variables_excluded}')
    # Calculate the standard deviation based on other error types
    df_out['SD'] = df_out.apply(sd_from_error_estimates, axis=1)
    # Add missing columns (metadata) to the dataframe
    # For numerical metadata, calculate their weighed average
    df_weighted_values = weighted_avg(df, numerical_covariates, by_cols)
    df_out.set_index(by_cols, inplace=True)
    # For categorical metadata, check that they are all the same (should be), and take the first value
    # Check if the metadata is the same in each group
    df_categorical_group = df.groupby(by_cols)
    df_unequal_categories = df_categorical_group[categorical_covariates].nunique().eq(1)
    for categorical_covariate in categorical_covariates:
        unequal_metadata = df_unequal_categories.index[df_unequal_categories[categorical_covariate] == False].tolist()
        if len(unequal_metadata) > 0:
            print(categorical_covariate)
            # print(unequal_metadata)
    df_categories = df_categorical_group.first()
    df_out = pd.merge(left=df_out,
                      right=df_categories[categorical_covariates],
                      how='inner', left_index=True, right_index=True)
    df_out = pd.merge(left=df_out,
                      right=df_weighted_values,
                      on=by_cols,
                      how='left')
    return df_out


def empty_columns(df):
    return df.columns[df.columns.str.startswith('Unnamed')]


def drop_nans(df):
    df.dropna(axis=0, how='all', inplace=True)
    df.dropna(axis=1, how='all', inplace=True)


def studies_variables(df, variable_column, variable_dict, count_column, df_response_variables):
    variables = list(variable_dict.keys())
    df_out = pd.DataFrame(columns=variables, index=[False, True])
    df = df.copy()
    for variable, variable_list in variable_dict.items():
        df[variable] = df_response_variables[variable_column].isin(variable_list)
        grouped = df.groupby(variable)[count_column].nunique()
        df_out.at[True, variable] = grouped.loc[True]
        df_out.at[False, variable] = df[count_column].nunique() - grouped.loc[True]
        df_out.at['Total', variable] = df[count_column].nunique()
    return df_out


def sd_from_error_estimates(df, error_type_column='Error estimate type', error_column='Error estimate'):
    if df[error_type_column] == 'SD':
        return df[error_column]
    elif df[error_type_column] == 'SE':
        return df[error_column] * np.sqrt(df['replicates'])
    elif df[error_type_column] == 'range':
        # Estimate the standard deviation using the range rule (95% of the data can be found within 2 standard
        # deviations)
        max = df['Response value'] + df[error_column]
        min = df['Response value'] - df[error_column]
        return (max - min) / 4
    elif df[error_type_column] == '95 % confidence intervals':
        CI_upper = df['Response value'] * (1 + df[error_column])
        CI_lower = df['Response value'] * df[error_column]
        sd = ((CI_upper - CI_lower) * np.sqrt(df['replicates'])) / (2 * 1.96)
        return sd
    else:
        return np.nan


def trawled_covariates(df, covariates, grouping_types, grouping_id_map, cat_covariates):
    grouping_id_impact = []
    for grouping_type in grouping_types:
        if pd.notna(grouping_id_map[grouping_type]):
            grouping_id_impact.append(grouping_id_map[grouping_type])
    df_trawled = df[df['Grouping_id'].isin(grouping_id_impact)]
    df_trawled_weighted = weighted_avg(df=df_trawled.copy(), data_cols=covariates, by_cols='Grouping_id')
    if not df_trawled_weighted.empty:
        df_trawled_weighted = df_trawled_weighted[covariates].mean()
    else:
        df_trawled_weighted[covariates] = np.nan
    if not set(cat_covariates).issubset(df_trawled.columns):
        print('stop')
    if not df_trawled.empty:
        df_trawled_cat_covariates = df_trawled[cat_covariates].iloc[0]
    else:
        df_trawled_cat_covariates = pd.DataFrame(columns=cat_covariates)
        df_trawled_cat_covariates[cat_covariates] = np.nan
    grouping_id_control =  grouping_id_map['After_control']
    if pd.notna(grouping_id_control):
        df_trawled_cat_covariates['Control_historically_trawled'] = df[df['Grouping_id'] == grouping_id_control]['Historically fished'].iloc[0]
    else:
        df_trawled_cat_covariates['Control_historically_trawled'] = np.nan
    df_out = pd.concat([df_trawled_weighted, df_trawled_cat_covariates])
    if df_out.empty:
        df_out = df_out.reindex(index = range(1))
    else:
        return df_out


def assert_data_comparability(df, row, columns, grouping_types):
    print(f'Checking data comparability of {row["article_id"]}')
    data_dict = {grouping_type: df[df['Grouping_id'] == row[grouping_type]][columns].reset_index(drop=True)
                 for grouping_type in grouping_types}
    if row['Study type'] == 'Control Impact':
        assert data_dict['Before_control'].empty and data_dict['Before_control'].empty, 'Before variables are not empty'
    elif row['Study type'] == 'Before After':
        assert data_dict['Before_control'].empty and data_dict['After_control'].empty, 'Control variables are not empty'
    elif row['Study type'] == 'Before After Control Impact':
        assert not data_dict['Before_control'].empty and not data_dict['Before_impact'].empty and \
               not data_dict['Before_impact'].empty and not data_dict['After_control'].empty


def map_studies(df_sample_mapping, df_all, columns, grouping_types, covariates, cat_covariates):
    df_out = pd.DataFrame()
    for i, row in df_sample_mapping.iterrows():
        # Assert that all the information is comparable for all grouped variables
        assert_data_comparability(df_all, row, columns, grouping_types)
        if row['Study type'] in ['Control Impact', 'Before After', 'Before After Control Impact']:
            df_comparison = pd.DataFrame()
            for grouping_type in grouping_types:
                if pd.notna(row[grouping_type]):
                    data = df_all[df_all['Grouping_id'] == row[grouping_type]]
                    if not data.empty:
                        df_temp = pd.DataFrame()
                        df_temp['Response variable'] = data['Response variable']
                        df_temp['Sample core depth slice'] = data['Sample core depth slice']
                        df_temp[f'Mean_{grouping_type}'] = data['Response value']
                        df_temp[f'Term_var_{grouping_type}'] = data['Error estimate type']
                        df_temp[f'Value_var_{grouping_type}'] = data['Error estimate']
                        df_temp[f'SD_{grouping_type}'] = data['SD']
                        df_temp[f'n_{grouping_type}'] = data['replicates']
                        df_temp[f'SD_{grouping_type}'].fillna(df_temp[f'Value_var_{grouping_type}'], inplace=True)
                        df_temp[columns] = data[columns]
                        df_temp[covariates + cat_covariates + ['Control_historically_trawled']] = trawled_covariates(
                            df=df_all, covariates=covariates,
                            grouping_types=['After_impact'],
                            grouping_id_map=row,
                            cat_covariates=cat_covariates)
                        if not df_comparison.empty:
                            variable_meta_analysis = [f'Mean_{grouping_type}', f'Term_var_{grouping_type}',
                                                      f'Value_var_{grouping_type}', f'SD_{grouping_type}',
                                                      f'n_{grouping_type}', 'Response variable',
                                                      'Sample core depth slice']
                            df_comparison = pd.merge(df_comparison.reset_index(drop=True),
                                                     df_temp[variable_meta_analysis].reset_index(drop=True),
                                                     on=['Response variable', 'Sample core depth slice'],
                                                     how='inner')
                        else:
                            df_comparison = df_temp
                    else:
                        print(f'No data for Grouping_id {row[grouping_type]} and article_id {row["article_id"]}')
            df_out = df_out.append(df_comparison)
            df_out.reset_index(drop=True)
    print(df_out.shape)
    df_out.drop_duplicates(inplace=True)
    print(df_out.shape)
    return df_out
