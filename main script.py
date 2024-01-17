import pandas as pd
import numpy as np
from scipy.stats import bootstrap

# load input data
expression_matrix = pd.read_csv("algorithm_data/small_expression_matrix.csv",index_col=[0])
annotation_table = pd.read_csv("algorithm_data/small_annotation_table.csv",index_col=[0])

# define effector genes by the user
effector_genes = ['ABCA1', 'ABCA2', 'ABCA3', 'NAT2']

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
    # then put the calculated mean values in a new column named "effector_genes_average" in the annotation table
    effector_genes_expression_matrix = expression_matrix[expression_matrix.columns.intersection(effector_genes)]
    annotation_table['effector_genes_average'] = effector_genes_expression_matrix.mean(axis=1)

    # prepper a groups_df object, which will hold the calculated 'median', 'CI_low', 'CI_high' per group
    groups = annotation_table['group'].unique()
    groups_df = pd.DataFrame(index=[groups], columns=['median', 'CI_low', 'CI_high'])

    # calculate the median average expression per group
    for group in groups:
       group_median = np.median((annotation_table.loc[annotation_table['group'] == group])['effector_genes_average'])
       groups_df.loc[group, 'median'] = group_median

    # calculate CI_low, CI_high per group
    for group in groups:
        # collect all the group's sample_id into group_samples object
        group_samples = annotation_table.loc[annotation_table['group'] == group]['sample_id']
        # collect the average expression of these samples into group_average_vector
        group_average_vector = annotation_table.loc[annotation_table['sample_id'].isin(group_samples)]['effector_genes_average']
        group_average_vector = (group_average_vector.to_numpy(),) # samples must be in a sequence
        # run bootstrapping on group_average_vector and extract the CI_low, CI_high results
        res = bootstrap(group_average_vector, np.median, confidence_level=0.9)
        CI_low, CI_high = res.confidence_interval
        groups_df.loc[group, 'CI_low'] = CI_low
        groups_df.loc[group, 'CI_high'] = CI_high

    # Generate final table to display for the customer
    final_table = ((annotation_table.drop(['sample_id', 'effector_genes_average'], axis=1)).drop_duplicates())
    final_table.index = [groups]
    final_table['median'] = groups_df['median']
    final_table['CI_low'] = groups_df['CI_low']
    final_table['CI_high'] = groups_df['CI_high']
    final_table = final_table.sort_values(by='CI_low',ascending=False)
    final_table.pop("n_per_indication")
    # generate the ranking column
    final_table['rank'] = np.arange(1, len(final_table)+1, 1)
    final_table.insert(0, "rank", final_table.pop("rank"))
    # print results
    print(final_table.to_string())










