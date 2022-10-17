import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import datetime

date = str(datetime.date.today().day) + '_' + str(datetime.date.today().month) + '_' + str(datetime.date.today().year)

# Clean up data functions
from data_processing_functions import to_number, group_studies, empty_columns, drop_nans, studies_variables, map_studies

dir = '../Data'
db_file_name = 'Bottom_Trawling_DB_2022_09_15.xlsx'

# Opening data from the database
df_articles = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='Meta_analysis_papers')
drop_nans(df_articles)
df_response_variables = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='MA_response_variables')
drop_nans(df_response_variables)
df_study_type = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='MA_study_type')
drop_nans(df_study_type)
df_sampling_location = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='MA_sampling_location_covariates')
drop_nans(df_sampling_location)
df_sample_mapping = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='MA_sampling_mapping')
drop_nans(df_sample_mapping)
df_study_region = pd.read_excel(os.path.join(dir, db_file_name), sheet_name='MA_region')
drop_nans(df_study_region)

### Grouping tables
# Group study type and article metadata tables
df_papers_study_type = pd.merge(left=df_articles[['article_id', 'First author', 'Year', 'Journal', 'DOI', 'Reviewer',
                                                  'Quality checked', 'Included in meta-analysis', 'Notes_paper']],
                                right=df_study_type[['article_id', 'article_type_id', 'Experimental/observational',
                                                     'Harmonized_study type', 'Study type', 'Notes_study']],
                                how='inner',
                                on=['article_id'])
# Take only studies that are included in the meta-analysis (second screening)
df_papers_study_type = df_papers_study_type[df_papers_study_type['Included in meta-analysis'] == True]

# Group with study region table
df_papers_study_type = pd.merge(left=df_papers_study_type, right=df_study_region,
                                how='inner',
                                on=['article_id', 'article_type_id'])
df_papers_study_type.rename(columns={'Latitude': 'Latitude_study', 'Longitude': 'Longitude_study'}, inplace=True)

# Group sampling_location and covariates
df_response_variables_sampling = pd.merge(left=df_sampling_location, right=df_response_variables,
                                          how='inner',
                                          on=['article_id', 'article_type_id', 'sampling_id'])
# Convert response value to number (in case it didn't read the column properly)
df_response_variables_sampling['Response value'] = pd.to_numeric(df_response_variables_sampling['Response value'],
                                                                 errors='coerce')
# Fill empty replicates row to 1 (if no data is provided, it is because they only had 1 replicate)
df_response_variables_sampling['replicates'].fillna(1, inplace=True)

# Merge with the previous tables
df_response_variables_sampling_study_type = pd.merge(left=df_response_variables_sampling,
                                                     right=df_papers_study_type[['article_id',
                                                                                 'article_type_id',
                                                                                 'Experimental/observational',
                                                                                 'Study type', 'Harmonized_study type',
                                                                                 'Location', 'Region',
                                                                                 'Longhurst province',
                                                                                 'Latitude_study', 'Longitude_study']],
                                                     how='inner',
                                                     on=['article_id', 'article_type_id'])

# Group samples and calculate weighed mean and standard deviation based on the grouping strategy
# ("Grouping_id" column in "MA_sampling_location_covariates" table)
all_variables = df_response_variables_sampling_study_type['Response variable'].unique()
num_covariates = ['Water depth (m)', 'Trawling effort_numerical_harmonized',
                  'Trawling_effort_GFW',
                  'Model_Chl_bottom (mg/m3)', 'Model_Current_velocity (m/s)',
                  'Model_Dissolved_Oxygen (mol/m3)', ' Model_NO3_bottom (mol/m3)',
                  'Model_PO4_bottom (mol/m3)', 'Model_Sal_bottom',
                  'Model_Silicate_bottom (mol/m3)', 'Model_Temp_bottom',
                  'Model_NPP_surface (g/m3/day)', 'dist_shore (m)', 'Time since first disturbance (years)']
cat_covariates = ['Habitat type_harmonized', 'Seasonality_harmonized', 'Trawling gear type_harmonized',
                  'Historically fished', 'Trawling effort_categorical', 'Trawling effort_units_harmonized']
metadata_columns = ['article_id', 'article_type_id', 'Trawled/untrawled', 'Before/After',
                    'Experimental/observational', 'Study type',
                    'Harmonized_study type', 'Location', 'Region', 'Longhurst province', 'Latitude_study',
                    'Longitude_study']
grouping_cols = ['Grouping_id', 'Sampling_type', 'Sample core depth slice']
df_group_studies = group_studies(df=df_response_variables_sampling_study_type, variables=all_variables,
                                 numerical_covariates=['Latitude',
                                                       'Longitude',
                                                       'Sampling date'] + num_covariates,
                                 categorical_covariates = metadata_columns + cat_covariates,
                                 data_col='Response value',
                                 weight_col='replicates',
                                 by_cols=grouping_cols)

# Map out the studies into individual rows for the meta-analysis
grouping_types = ['Before_control', 'Before_impact', 'After_control', 'After_impact']
df_out = map_studies(df_sample_mapping=df_sample_mapping, df_all=df_group_studies, columns=metadata_columns,
                     grouping_types=grouping_types, cat_covariates=cat_covariates, covariates=num_covariates)

df_out.to_excel(os.path.join(dir, 'trawling_meta-analysis_grouped' + date + '.xlsx'), index=False)
