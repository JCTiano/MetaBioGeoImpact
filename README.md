# Metastudy
 Workflow from dataset to figures and results shown in the article.

# Folder structures

*/dataset* contains the dataset as a flatfile and associated metadata

*/database_grouping* contains python scripts to group observations together (calculate weighed averages if needed) into individual observations (rows) for the meta-analysis

*/markdown_scripts* contains r-markdown scripts in which calculations are perfomed/output is generated step-by-step.

*/r_objects* contains .rda files etc needed to support running the code.

*/source_code* contains .R files containing functions used in the markdown scripts.


# Order of operations

## 1_add_sd

- Dataset is read.
- SD are made positive.
- Add columns to add results to later.
- Add missing standard deviations based on known means.

## 2_EF_calculation

- Assign response, sd, and number of replicates for log response calculation.
- Lines with SD = NA are removed, they cause errors later anyway and have no means.
  - This is the data from study 169.
- Rows where no response ratio is calculated are removed (130).
- Response ratio is calculated.
- Standard deviations and weights are calculated
  - VarLnRR is NaN are set to 0.
  - Minimal weight is set to 0.1 (35 values).
  - Inf weights where VarLnRR is 0 (1172 values).
  - Some weights are extremely high, I don't know what to do with this. 
    - For now nothing.
- Cohen's d calculation.
  - Using the "add value" option in escalc to get weights.
- This is converted to Hedges g.

## 3_EF_Gradient_studies

X Later

## 4_Modelfitting

- Make column names easier to work with.
- Refactoring some of the depth categories to
  - 0-1, 1-2, 2-5, 5-10, 10+, though this leaves uncertainties.
- Subset dataframe in chunks to work with.
- Start modelling.

## 5_Metaplots

- 














































