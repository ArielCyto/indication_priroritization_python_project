import pandas as pd
import numpy as np
from scipy.stats import bootstrap

# load input data
expression_matrix = pd.read_csv("algorithm_data/small_expression_matrix.csv",index_col=[0])
annotation_table = pd.read_csv("algorithm_data/small_annotation_table.csv",index_col=[0])
expression_matrix.head()
annotation_table.head()

# define effector genes by the user
effector_genes = ['ABCA1',	'ABCA2', 'ABCA3', 'NAT2']

# check if the genes are available:
submit = False
for gene in effector_genes:
    gene_status = (gene in expression_matrix.columns)
    if gene_status:
        submit = True
    else:
       print("gene", gene, "is an invalid input, please remove it and try again")
       submit = False
       # check for alternative gene names
       all_genes = pd.Series(expression_matrix.columns.values)
       all_genes.index = expression_matrix.columns.values
       result = all_genes.str.contains(pat=gene.upper())
       if result.any():
           alternative_gene_names = result[np.where(result)[0]].index
           print("optional alternatives: ")
           print(alternative_gene_names.to_list())
       break


# run analysis
if (submit == True):
    # per sample (row), calculate the selected effector_genes average expression
    mini_expression_matrix = expression_matrix[expression_matrix.columns.intersection(effector_genes)]
    mini_expression_matrix['effector_genes_average'] = mini_expression_matrix.mean(axis=1)

    # work with the annotation_table object with all the necessary columns
    annotation_table['effector_genes_average'] = mini_expression_matrix['effector_genes_average']

    # calculate the group median
    groups = annotation_table['group'].unique()
    groups_df = pd.DataFrame(index=[groups], columns=['median', 'CI_low', 'CI_high'])
    for group in groups:
       tmp = annotation_table.loc[annotation_table['group'] == group]
       tmp_median = np.median(tmp['effector_genes_average'])
       groups_df.loc[group, 'median'] = tmp_median

    # calculate the group CI_low, CI_high
    for group in groups:
        group_samples = annotation_table.loc[annotation_table['group'] == group]['sample_id']
        group_average_vector = annotation_table.loc[annotation_table['sample_id'].isin(group_samples)]['effector_genes_average']
        group_average_vector = (group_average_vector.to_numpy(),) # samples must be in a sequence
        res = bootstrap(group_average_vector, np.median, confidence_level=0.9)
        CI_low, CI_high = res.confidence_interval
        groups_df.loc[group, 'CI_low'] = CI_low
        groups_df.loc[group, 'CI_high'] = CI_high

    # Generate final table for the customer
    final_table = ((annotation_table.drop(['sample_id', 'effector_genes_average'], axis=1)).drop_duplicates())
    final_table.index = [groups]
    final_table['median'] = groups_df['median']
    final_table['CI_low'] = groups_df['CI_low']
    final_table['CI_high'] = groups_df['CI_high']
    final_table = final_table.sort_values(by='CI_low',ascending=False)
    final_table['rank'] = np.arange(1, len(final_table)+1, 1)
    final_table.pop("n_per_indication")
    final_table.insert(0, "rank", final_table.pop("rank"))
    print(final_table.to_string())









